"""
VCF file operations.

Provides functions for reading and processing VCF files to extract SNPs.
Uses cyvcf2 for fast indexed access to compressed VCF files.

Note: Always extracts first alt allele (alt[0]) for multi-allelic sites.
"""

from __future__ import annotations
from pathlib import Path
from typing import List, Tuple, Optional
import logging
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

from cyvcf2 import VCF

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNP

logger = logging.getLogger(__name__)


def _count_snps_for_chrom(vcf_path: Path, chrom: str) -> int:
    """Count biallelic SNPs in one chromosome."""
    count = 0
    try:
        vcf = VCF(str(vcf_path))
        for variant in vcf(chrom):
            # Check if biallelic SNP
            if len(variant.REF) == 1 and variant.ALT and len(variant.ALT) > 0:
                first_alt = variant.ALT[0]
                if len(first_alt) == 1 and first_alt in "ACGT":
                    count += 1
        vcf.close()
    except Exception:
        pass
    return count


def count_snps_in_vcf(vcf_path: Path, threads: int = 1) -> Result[int, str]:
    """
    Fast count of total biallelic SNPs in VCF file.
    
    Only counts, does not extract SNP data. Much faster than full extraction.
    
    Args:
        vcf_path: Path to VCF file (.vcf.gz with .tbi index)
        threads: Number of parallel processes
        
    Returns:
        Result containing total SNP count
    """
    try:
        vcf = VCF(str(vcf_path))
        chromosomes = vcf.seqnames
        vcf.close()
    except Exception as e:
        return Err(f"Failed to read VCF: {e}")
    
    total = 0
    if threads == 1:
        for chrom in chromosomes:
            total += _count_snps_for_chrom(vcf_path, chrom)
    else:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(_count_snps_for_chrom, vcf_path, chrom) for chrom in chromosomes]
            for future in as_completed(futures):
                total += future.result()
    
    return Ok(total)


def invert_bed_regions(bed_path: Path, fai_path: Path) -> List[Tuple[str, int, int]]:
    """
    Generate inverted BED regions (complement of input BED).
    
    Uses reference .fai file to get chromosome lengths, then returns
    all regions NOT covered by the input BED.
    
    Args:
        bed_path: Input BED file
        fai_path: Reference genome .fai index file
        
    Returns:
        List of (chrom, start, end) tuples for inverted regions
    """
    # Read chromosome lengths from .fai
    chrom_lengths = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_lengths[parts[0]] = int(parts[1])
    
    # Read and sort BED regions by chromosome and start
    bed_regions = {}
    with open(bed_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                if chrom not in bed_regions:
                    bed_regions[chrom] = []
                bed_regions[chrom].append((start, end))
    
    # Sort regions within each chromosome
    for chrom in bed_regions:
        bed_regions[chrom].sort()
    
    # Generate inverted regions
    inverted = []
    for chrom, length in chrom_lengths.items():
        if chrom not in bed_regions:
            # No BED regions on this chromosome, keep entire chromosome
            inverted.append((chrom, 0, length))
            continue
        
        regions = bed_regions[chrom]
        current_pos = 0
        
        for start, end in regions:
            # Add gap before this region
            if current_pos < start:
                inverted.append((chrom, current_pos, start))
            current_pos = max(current_pos, end)
        
        # Add remaining region after last BED entry
        if current_pos < length:
            inverted.append((chrom, current_pos, length))
    
    return inverted


def _extract_snps_for_chromosome(
    vcf_path: Path,
    chrom: str,
    bed_regions: Optional[List[Tuple[str, int, int]]],
) -> Tuple[str, List[SNP], int, int]:
    """
    Extract SNPs from a single chromosome using cyvcf2 (for multiprocessing).
    
    Always extracts first alt allele (alt[0]) regardless of ploidy.
    
    Args:
        vcf_path: Path to VCF file
        chrom: Chromosome name
        bed_regions: Optional list of (chrom, start, end) tuples to extract from
        
    Returns:
        Tuple of (chromosome, snp_list, multi_allelic_count, indel_count)
    """
    logger.info(f"→ Processing chromosome {chrom}...")
    
    snps = []
    multi_allelic = 0
    indels = 0
    
    # Batch append optimization
    batch = []
    batch_size = 10000
    
    try:
        vcf = VCF(str(vcf_path))
        
        # Fetch variants
        if bed_regions:
            # Filter by BED regions using indexed VCF access
            # Suppress cyvcf2 warnings for empty intervals
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="no intervals found")
                for region_chrom, region_start, region_end in bed_regions:
                    if region_chrom != chrom:
                        continue
                    for variant in vcf(f"{region_chrom}:{region_start}-{region_end}"):
                        if variant.CHROM != chrom:
                            continue
                        
                        # Check if SNP (single base ref and alt)
                        if len(variant.REF) != 1:
                            indels += 1
                            continue
                        
                        if not variant.ALT or len(variant.ALT) == 0:
                            continue
                        
                        first_alt = variant.ALT[0]
                        if len(first_alt) != 1 or first_alt not in "ACGT":
                            indels += 1
                            continue
                        
                        # Count multi-allelic
                        if len(variant.ALT) > 1:
                            multi_allelic += 1
                        
                        # Create SNP (always use first alt)
                        snp = SNP.from_vcf_record(
                            chrom=variant.CHROM,
                            pos=variant.POS,
                            ref=variant.REF,
                            alt=first_alt,
                        )
                        batch.append(snp)
                        
                        # Batch extend
                        if len(batch) >= batch_size:
                            snps.extend(batch)
                            batch.clear()
        else:
            # Fetch entire chromosome
            for variant in vcf(chrom):
                # Check if SNP (single base ref and alt)
                if len(variant.REF) != 1:
                    indels += 1
                    continue
                
                if not variant.ALT or len(variant.ALT) == 0:
                    continue
                
                first_alt = variant.ALT[0]
                if len(first_alt) != 1 or first_alt not in "ACGT":
                    indels += 1
                    continue
                
                # Count multi-allelic
                if len(variant.ALT) > 1:
                    multi_allelic += 1
                
                # Create SNP (always use first alt)
                snp = SNP.from_vcf_record(
                    chrom=variant.CHROM,
                    pos=variant.POS,
                    ref=variant.REF,
                    alt=first_alt,
                )
                batch.append(snp)
                
                # Batch extend
                if len(batch) >= batch_size:
                    snps.extend(batch)
                    batch.clear()
        
        # Flush remaining batch
        if batch:
            snps.extend(batch)
            batch.clear()
        
        logger.info(f"✓ Chromosome {chrom}: {len(snps)} SNPs extracted")
        if multi_allelic > 0:
            logger.info(f"  └─ {multi_allelic} multi-allelic sites (extracted first alt)")
        if indels > 0:
            logger.debug(f"  └─ {indels} indels skipped")
        
        return (chrom, snps, multi_allelic, indels)
        
    except Exception as e:
        logger.error(f"Error processing {chrom}: {e}")
        return (chrom, [], 0, 0)


def extract_snps_from_vcf(
    vcf_path: Path,
    bed_path: Optional[Path] = None,
    bed_mode: str = "keep",
    fai_path: Optional[Path] = None,
    threads: int = 1,
) -> Result[list[SNP], str]:
    """
    Extract SNPs from VCF file using cyvcf2 (with multiprocessing).
    
    Always extracts first alt allele (alt[0]) for multi-allelic sites.
    Skips indels (variants where REF or ALT length != 1).
    
    For remove mode: generates inverted BED regions from .fai file,
    then uses indexed VCF access (same speed as keep mode).
    
    Args:
        vcf_path: Path to VCF file (.vcf.gz with .tbi index)
        bed_path: Optional BED file to filter regions
        bed_mode: 'keep' to keep BED regions, 'remove' to exclude them
        fai_path: Reference .fai file (required for remove mode)
        threads: Number of parallel processes (default: 1)
        
    Returns:
        Result containing dict with keys: 'snps', 'bed_mode', 'total_before_bed', 'total_after_bed'
    """
    # Get chromosome list from VCF header
    try:
        vcf = VCF(str(vcf_path))
        chromosomes = vcf.seqnames
        vcf.close()
        
        logger.info(f"Found {len(chromosomes)} chromosomes in VCF")
        logger.info(f"Using {threads} threads for parallel extraction")
        
    except Exception as e:
        return Err(f"Failed to read VCF header: {e}")
    
    # Process BED regions based on mode
    bed_regions = None
    bed_applied = bed_path is not None
    bed_mode_used = bed_mode if bed_applied else None
    
    if bed_path:
        if bed_mode == "remove":
            # Generate inverted BED for remove mode
            if fai_path is None:
                return Err("remove mode requires fai_path parameter")
            if not fai_path.exists():
                return Err(f"Reference index not found: {fai_path}")
            
            logger.info(f"Generating inverted BED regions from {bed_path}")
            bed_regions = invert_bed_regions(bed_path, fai_path)
            logger.info(f"Generated {len(bed_regions)} inverted regions")
        else:
            # Read BED regions for keep mode
            bed_regions = []
            with open(bed_path) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        bed_regions.append((parts[0], int(parts[1]), int(parts[2])))
            logger.info(f"Loaded {len(bed_regions)} BED regions (keep mode)")
    
    all_snps = []
    total_multi_allelic = 0
    total_indels = 0
    
    if threads == 1:
        # Serial processing
        for chrom in chromosomes:
            chrom_name, snps, multi_allelic, indels = _extract_snps_for_chromosome(
                vcf_path, chrom, bed_regions
            )
            all_snps.extend(snps)
            total_multi_allelic += multi_allelic
            total_indels += indels
    else:
        # Parallel processing by chromosome
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {}
            for chrom in chromosomes:
                future = executor.submit(
                    _extract_snps_for_chromosome,
                    vcf_path, chrom, bed_regions
                )
                futures[future] = chrom
            
            for future in as_completed(futures):
                chrom_name = futures[future]
                try:
                    chrom, snps, multi_allelic, indels = future.result()
                    all_snps.extend(snps)
                    total_multi_allelic += multi_allelic
                    total_indels += indels
                except Exception as e:
                    logger.error(f"Failed to process {chrom_name}: {e}")
    
    # Summary statistics
    if total_multi_allelic > 0:
        logger.info(f"⚠ Found {total_multi_allelic} multi-allelic sites (extracted first alt)")
    if total_indels > 0:
        logger.debug(f"Skipped {total_indels} indels")
    if total_indels > 0:
        logger.info(f"Skipped {total_indels} indels (only SNPs extracted)")
    
    result = {
        'snps': all_snps,
        'bed_applied': bed_applied,
        'bed_mode': bed_mode_used,
    }
    
    return Ok(result)
