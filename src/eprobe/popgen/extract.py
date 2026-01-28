"""
SNP extraction from VCF files.

Extracts biallelic SNPs from VCF files with flanking sequences from
reference genome. Optionally applies cluster filtering to remove
SNPs in dense regions.

Optimizations:
- Multiprocessing by chromosome
- On-demand reference loading with pysam
- O(n log n) clustering algorithm
- Batch DataFrame write
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import numpy as np
import pysam

from eprobe.core.result import Result, Ok, Err, try_except
from eprobe.core.vcf import extract_snps_from_vcf, VCFReader
from eprobe.core.models import SNP

logger = logging.getLogger(__name__)


@dataclass
class ClusterConfig:
    """Configuration for cluster filtering."""
    flank: int = 60
    max_snp: int = 3
    enabled: bool = True


def detect_clusters(
    snps: List[SNP],
    flank: int,
    max_snp: int,
) -> List[int]:
    """
    Detect SNP clusters using optimized O(n log n) algorithm.
    
    Uses sliding window with two pointers instead of nested loops.
    
    Args:
        snps: List of SNPs (unsorted)
        flank: Flanking distance for cluster detection
        max_snp: Maximum allowed SNPs in a 2*flank window
        
    Returns:
        List of indices of SNPs to keep (non-clustered)
    """
    if not snps:
        return []
    
    # Group by chromosome
    chrom_snps: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for idx, snp in enumerate(snps):
        chrom_snps[snp.chrom].append((snp.pos, idx))
    
    keep_indices: List[int] = []
    
    for chrom, pos_idx_list in chrom_snps.items():
        # Sort by position: O(n log n)
        pos_idx_list.sort()
        
        positions = [pos for pos, _ in pos_idx_list]
        indices = [idx for _, idx in pos_idx_list]
        
        # Sliding window: O(n)
        left = 0
        for right in range(len(positions)):
            pos = positions[right]
            
            # Move left pointer forward
            while left < right and positions[left] < pos - flank:
                left += 1
            
            # Count SNPs in window [pos-flank, pos+flank]
            window_count = 0
            for i in range(left, len(positions)):
                if positions[i] > pos + flank:
                    break
                window_count += 1
            
            if window_count <= max_snp:
                keep_indices.append(indices[right])
            else:
                logger.debug(f"Filtered cluster SNP: {chrom}:{pos} (cluster size: {window_count})")
    
    return sorted(keep_indices)


def apply_cluster_filter(
    snps: List[SNP],
    config: ClusterConfig,
) -> Result[List[SNP], str]:
    """
    Apply cluster filtering to remove SNPs in dense regions.
    
    Args:
        snps: Input list of SNPs
        config: Cluster filtering configuration
        
    Returns:
        Result containing filtered SNP list
    """
    if not config.enabled:
        return Ok(snps)
    
    logger.info(f"Applying cluster filter (flank={config.flank}, max_snp={config.max_snp})")
    
    keep_indices = detect_clusters(snps, config.flank, config.max_snp)
    filtered_snps = [snps[i] for i in keep_indices]
    
    removed = len(snps) - len(filtered_snps)
    logger.info(f"Cluster filter: removed {removed} SNPs, {len(filtered_snps)} remaining")
    
    return Ok(filtered_snps)


def extract_snps_with_flanks_for_chrom(
    chrom: str,
    chrom_snps: List[SNP],
    reference_path: Path,
    flank: int,
) -> List[SNP]:
    """
    Extract flanking sequences for SNPs on one chromosome.
    
    Uses pysam.FastaFile for on-demand loading (memory efficient).
    Designed for multiprocessing - processes one chromosome at a time.
    
    Args:
        chrom: Chromosome name
        chrom_snps: SNPs on this chromosome
        reference_path: Path to indexed reference FASTA
        flank: Flanking region length
        
    Returns:
        List of SNPs with flanking sequences
    """
    snps_with_flanks = []
    
    try:
        ref_fasta = pysam.FastaFile(str(reference_path))
        
        # Get chromosome length
        if chrom not in ref_fasta.references:
            logger.warning(f"Chromosome {chrom} not in reference")
            return []
        
        chrom_len = ref_fasta.get_reference_length(chrom)
        
        for snp in chrom_snps:
            # Calculate flanking region (0-based for pysam)
            start = max(0, snp.pos - 1 - flank)
            end = min(chrom_len, snp.pos + flank)
            
            # Check if we have enough flanking sequence
            actual_left_flank = snp.pos - 1 - start
            actual_right_flank = end - snp.pos
            
            if actual_left_flank < flank or actual_right_flank < flank:
                continue
            
            # Fetch sequences on-demand (very memory efficient)
            left_flank_seq = ref_fasta.fetch(chrom, start, snp.pos - 1)
            right_flank_seq = ref_fasta.fetch(chrom, snp.pos, end)
            
            # Validate reference allele
            ref_at_pos = ref_fasta.fetch(chrom, snp.pos - 1, snp.pos)
            if ref_at_pos.upper() != snp.ref.upper():
                logger.warning(
                    f"Reference mismatch at {chrom}:{snp.pos}: "
                    f"VCF={snp.ref}, Ref={ref_at_pos}"
                )
                continue
            
            # Create updated SNP with flanks
            updated_snp = SNP(
                chrom=snp.chrom,
                pos=snp.pos,
                ref=snp.ref,
                alt=snp.alt,
                snp_id=snp.snp_id,
                left_flank=left_flank_seq.upper(),
                right_flank=right_flank_seq.upper(),
            )
            snps_with_flanks.append(updated_snp)
        
        ref_fasta.close()
        
    except Exception as e:
        logger.error(f"Error processing chromosome {chrom}: {e}")
        return []
    
    return snps_with_flanks


def extract_snps_with_flanks_parallel(
    snps: List[SNP],
    reference_path: Path,
    flank: int,
    threads: int = 1,
) -> Result[List[SNP], str]:
    """
    Add flanking sequences using multiprocessing by chromosome.
    
    Args:
        snps: List of SNPs
        reference_path: Path to indexed reference FASTA (.fai required)
        flank: Flanking region length
        threads: Number of parallel processes
        
    Returns:
        Result containing SNPs with flanking sequences
    """
    # Group SNPs by chromosome
    chrom_snps_dict: Dict[str, List[SNP]] = defaultdict(list)
    for snp in snps:
        chrom_snps_dict[snp.chrom].append(snp)
    
    logger.info(f"Processing {len(chrom_snps_dict)} chromosomes with {threads} threads")
    
    all_snps_with_flanks = []
    
    if threads == 1:
        # Serial processing
        for chrom, chrom_snps in chrom_snps_dict.items():
            result = extract_snps_with_flanks_for_chrom(chrom, chrom_snps, reference_path, flank)
            all_snps_with_flanks.extend(result)
    else:
        # Parallel processing by chromosome
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {}
            for chrom, chrom_snps in chrom_snps_dict.items():
                future = executor.submit(
                    extract_snps_with_flanks_for_chrom,
                    chrom, chrom_snps, reference_path, flank
                )
                futures[future] = chrom
            
            for future in as_completed(futures):
                chrom = futures[future]
                try:
                    result = future.result()
                    all_snps_with_flanks.extend(result)
                    logger.info(f"Completed chromosome {chrom}: {len(result)} SNPs")
                except Exception as e:
                    logger.error(f"Failed processing chromosome {chrom}: {e}")
    
    return Ok(all_snps_with_flanks)


def run_extract(
    vcf_path: Path,
    reference_path: Path,
    output_prefix: Path,
    flank: int = 60,
    cluster_mode: bool = True,
    cluster_flank: int = 60,
    max_cluster_snp: int = 3,
    bed_path: Optional[Path] = None,
    threads: int = 1,
    force_biallelic: bool = False,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Extract SNPs from VCF file with flanking sequences (optimized).
    
    Optimizations:
    - Multiprocessing by chromosome (if threads > 1)
    - On-demand reference loading with pysam (low memory)
    - O(n log n) clustering algorithm
    - Batch DataFrame write
    
    Args:
        vcf_path: Input VCF file path
        reference_path: Reference genome FASTA path (requires .fai index)
        output_prefix: Output file prefix
        flank: Flanking region length
        cluster_mode: Enable cluster filtering
        cluster_flank: Flanking distance for cluster detection
        max_cluster_snp: Maximum SNPs allowed in cluster window
        bed_path: Optional BED file to restrict regions
        threads: Number of parallel processes (recommended: num chromosomes)
        force_biallelic: Force biallelic by taking first alt if multiple
        verbose: Enable verbose logging
        
    Returns:
        Result containing extraction statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Starting SNP extraction from {vcf_path}")
    logger.info(f"Threads: {threads}, Force biallelic: {force_biallelic}")
    
    # Check reference index exists
    fai_path = Path(str(reference_path) + ".fai")
    if not fai_path.exists():
        return Err(f"Reference index not found: {fai_path}. Run 'samtools faidx {reference_path}'")
    
    # Step 1: Extract raw SNPs from VCF
    logger.info("Extracting SNPs from VCF...")
    vcf_result = extract_snps_from_vcf(vcf_path, bed_path, force_biallelic=force_biallelic)
    if vcf_result.is_err():
        return Err(f"VCF extraction failed: {vcf_result.unwrap_err()}")
    
    raw_snps = vcf_result.unwrap()
    logger.info(f"Extracted {len(raw_snps)} SNPs from VCF")
    
    if not raw_snps:
        return Err("No SNPs extracted from VCF file")
    
    # Step 2: Add flanking sequences (using multiprocessing)
    logger.info(f"Extracting flanking sequences (flank={flank}, threads={threads})...")
    flank_result = extract_snps_with_flanks_parallel(
        raw_snps, reference_path, flank, threads
    )
    if flank_result.is_err():
        return Err(f"Flanking extraction failed: {flank_result.unwrap_err()}")
    
    snps_with_flanks = flank_result.unwrap()
    logger.info(f"{len(snps_with_flanks)} SNPs have complete flanking sequences")
    
    # Step 3: Apply cluster filter
    cluster_config = ClusterConfig(
        flank=cluster_flank,
        max_snp=max_cluster_snp,
        enabled=cluster_mode,
    )
    
    filter_result = apply_cluster_filter(snps_with_flanks, cluster_config)
    if filter_result.is_err():
        return Err(f"Cluster filtering failed: {filter_result.unwrap_err()}")
    
    final_snps = filter_result.unwrap()
    
    # Step 4: Save output (batch write)
    output_path = Path(str(output_prefix) + ".snps.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Direct DataFrame construction (more efficient)
    df = pd.DataFrame([
        {
            'chrom': snp.chrom,
            'pos': snp.pos,
            'ref': snp.ref,
            'alt': snp.alt,
            'snp_id': snp.snp_id or f"{snp.chrom}:{snp.pos}",
            'left_flank': snp.left_flank,
            'right_flank': snp.right_flank,
        }
        for snp in final_snps
    ])
    
    try:
        df.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Saved {len(final_snps)} SNPs to {output_path}")
    except Exception as e:
        return Err(f"Failed to save output: {e}")
    
    # Return statistics
    stats = {
        "raw_snp_count": len(raw_snps),
        "flanked_snp_count": len(snps_with_flanks),
        "final_snp_count": len(final_snps),
        "cluster_removed": len(snps_with_flanks) - len(final_snps) if cluster_mode else 0,
        "output_file": str(output_path),
    }
    
    return Ok(stats)
