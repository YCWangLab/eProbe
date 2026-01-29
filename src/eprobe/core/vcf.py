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


def _extract_snps_for_chromosome(
    vcf_path: Path,
    chrom: str,
    bed_path: Optional[Path],
) -> Tuple[str, List[SNP], int, int]:
    """
    Extract SNPs from a single chromosome using cyvcf2 (for multiprocessing).
    
    Always extracts first alt allele (alt[0]) regardless of ploidy.
    
    Args:
        vcf_path: Path to VCF file
        chrom: Chromosome name
        bed_path: Optional BED file to filter regions
        
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
        
        # Read BED regions if provided
        regions = []
        if bed_path:
            with open(bed_path) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        regions.append((parts[0], int(parts[1]), int(parts[2])))
        
        # Fetch variants
        if regions:
            # Filter by BED regions
            # Suppress cyvcf2 warnings for empty intervals
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="no intervals found")
                for region_chrom, region_start, region_end in regions:
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
    threads: int = 1,
) -> Result[list[SNP], str]:
    """
    Extract SNPs from VCF file using cyvcf2 (with multiprocessing).
    
    Always extracts first alt allele (alt[0]) for multi-allelic sites.
    Skips indels (variants where REF or ALT length != 1).
    
    Args:
        vcf_path: Path to VCF file (.vcf.gz with .tbi index)
        bed_path: Optional BED file to filter regions
        threads: Number of parallel processes (default: 1)
        
    Returns:
        Result containing list of SNPs
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
    
    all_snps = []
    total_multi_allelic = 0
    total_indels = 0
    
    if threads == 1:
        # Serial processing
        for chrom in chromosomes:
            chrom_name, snps, multi_allelic, indels = _extract_snps_for_chromosome(
                vcf_path, chrom, bed_path
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
                    vcf_path, chrom, bed_path
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
        logger.warning(f"⚠ Found {total_multi_allelic} multi-allelic sites (extracted first alt allele)")
    if total_indels > 0:
        logger.info(f"Skipped {total_indels} indels (only SNPs extracted)")
    
    return Ok(all_snps)
