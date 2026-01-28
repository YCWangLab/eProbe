"""
SNP extraction from VCF files.

Extracts biallelic SNPs from VCF files with flanking sequences from
reference genome. Optionally applies cluster filtering to remove
SNPs in dense regions.

This module corresponds to the original SNP_preprocessor.py functionality.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass
from collections import defaultdict

import pandas as pd
import numpy as np

from eprobe.core.result import Result, Ok, Err, try_except
from eprobe.core.fasta import read_fasta, extract_sequence
from eprobe.core.vcf import extract_snps_from_vcf, VCFReader
from eprobe.core.models import SNP, SNPDataFrame

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
    Detect SNP clusters and return indices of non-clustered SNPs.
    
    A cluster is defined as a region where more than max_snp SNPs
    fall within a 2*flank window. The cluster boundary is determined
    by the flank distance.
    
    Args:
        snps: List of SNPs sorted by position
        flank: Flanking distance for cluster detection
        max_snp: Maximum allowed SNPs in a 2*flank window
        
    Returns:
        List of indices of SNPs to keep (non-clustered)
    """
    if not snps:
        return []
    
    # Group SNPs by chromosome
    chrom_snps: Dict[str, List[Tuple[int, SNP]]] = defaultdict(list)
    for idx, snp in enumerate(snps):
        chrom_snps[snp.chrom].append((idx, snp))
    
    keep_indices: List[int] = []
    window_size = 2 * flank
    
    for chrom, indexed_snps in chrom_snps.items():
        # Sort by position within chromosome
        indexed_snps.sort(key=lambda x: x[1].pos)
        
        # Sliding window to count SNPs
        positions = [snp.pos for _, snp in indexed_snps]
        n = len(positions)
        
        # For each SNP, check if it's in a cluster
        for i, (idx, snp) in enumerate(indexed_snps):
            # Count SNPs within window centered on this SNP
            left_bound = snp.pos - flank
            right_bound = snp.pos + flank
            
            count = sum(1 for pos in positions if left_bound <= pos <= right_bound)
            
            if count <= max_snp:
                keep_indices.append(idx)
            else:
                logger.debug(f"Filtered cluster SNP: {chrom}:{snp.pos} (cluster size: {count})")
    
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


def extract_snps_with_flanks(
    snps: List[SNP],
    reference_sequences: Dict[str, str],
    flank: int,
) -> Result[List[SNP], str]:
    """
    Add flanking sequences to SNPs from reference genome.
    
    Args:
        snps: List of SNPs (may have partial or no flanks)
        reference_sequences: Reference genome sequences by chromosome
        flank: Flanking region length
        
    Returns:
        Result containing SNPs with flanking sequences added
    """
    snps_with_flanks = []
    
    for snp in snps:
        chrom_seq = reference_sequences.get(snp.chrom)
        if chrom_seq is None:
            logger.warning(f"Chromosome {snp.chrom} not in reference, skipping SNP at {snp.pos}")
            continue
        
        chrom_len = len(chrom_seq)
        
        # Calculate flanking region (0-based)
        start = max(0, snp.pos - 1 - flank)
        end = min(chrom_len, snp.pos + flank)
        
        # Check if we have enough flanking sequence
        actual_left_flank = snp.pos - 1 - start
        actual_right_flank = end - snp.pos
        
        if actual_left_flank < flank or actual_right_flank < flank:
            logger.debug(f"Insufficient flanking for {snp.chrom}:{snp.pos}, skipping")
            continue
        
        # Extract sequences
        left_flank_seq = chrom_seq[start:snp.pos - 1]
        right_flank_seq = chrom_seq[snp.pos:end]
        
        # Validate reference allele
        ref_at_pos = chrom_seq[snp.pos - 1]
        if ref_at_pos.upper() != snp.ref.upper():
            logger.warning(
                f"Reference mismatch at {snp.chrom}:{snp.pos}: "
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
    
    return Ok(snps_with_flanks)


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
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Extract SNPs from VCF file with flanking sequences.
    
    Main entry point for the extract command. Performs:
    1. VCF parsing and biallelic SNP extraction
    2. BED region filtering (optional)
    3. Reference sequence retrieval for flanking regions
    4. Cluster filtering (optional)
    5. Output TSV generation
    
    Args:
        vcf_path: Input VCF file path
        reference_path: Reference genome FASTA path
        output_prefix: Output file prefix
        flank: Flanking region length
        cluster_mode: Enable cluster filtering
        cluster_flank: Flanking distance for cluster detection
        max_cluster_snp: Maximum SNPs allowed in cluster window
        bed_path: Optional BED file to restrict regions
        threads: Number of threads (for future parallelization)
        verbose: Enable verbose logging
        
    Returns:
        Result containing extraction statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Starting SNP extraction from {vcf_path}")
    
    # Step 1: Extract raw SNPs from VCF
    logger.info("Extracting SNPs from VCF...")
    vcf_result = extract_snps_from_vcf(vcf_path, bed_path)
    if vcf_result.is_err():
        return Err(f"VCF extraction failed: {vcf_result.unwrap_err()}")
    
    raw_snps = vcf_result.unwrap()
    logger.info(f"Extracted {len(raw_snps)} biallelic SNPs from VCF")
    
    if not raw_snps:
        return Err("No SNPs extracted from VCF file")
    
    # Step 2: Read reference genome
    logger.info(f"Loading reference genome from {reference_path}...")
    ref_result = read_fasta(reference_path)
    if ref_result.is_err():
        return Err(f"Failed to read reference: {ref_result.unwrap_err()}")
    
    reference_sequences = ref_result.unwrap()
    logger.info(f"Loaded {len(reference_sequences)} chromosomes from reference")
    
    # Step 3: Add flanking sequences
    logger.info(f"Extracting flanking sequences (flank={flank})...")
    flank_result = extract_snps_with_flanks(raw_snps, reference_sequences, flank)
    if flank_result.is_err():
        return Err(f"Flanking extraction failed: {flank_result.unwrap_err()}")
    
    snps_with_flanks = flank_result.unwrap()
    logger.info(f"{len(snps_with_flanks)} SNPs have complete flanking sequences")
    
    # Step 4: Apply cluster filter
    cluster_config = ClusterConfig(
        flank=cluster_flank,
        max_snp=max_cluster_snp,
        enabled=cluster_mode,
    )
    
    filter_result = apply_cluster_filter(snps_with_flanks, cluster_config)
    if filter_result.is_err():
        return Err(f"Cluster filtering failed: {filter_result.unwrap_err()}")
    
    final_snps = filter_result.unwrap()
    
    # Step 5: Save output
    output_path = Path(str(output_prefix) + ".snps.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Create SNPDataFrame and save
    snp_df = SNPDataFrame.from_snps(final_snps)
    save_result = snp_df.to_tsv(output_path)
    if save_result.is_err():
        return Err(f"Failed to save output: {save_result.unwrap_err()}")
    
    logger.info(f"Saved {len(final_snps)} SNPs to {output_path}")
    
    # Return statistics
    stats = {
        "raw_snp_count": len(raw_snps),
        "flanked_snp_count": len(snps_with_flanks),
        "final_snp_count": len(final_snps),
        "cluster_removed": len(snps_with_flanks) - len(final_snps) if cluster_mode else 0,
        "output_file": str(output_path),
    }
    
    return Ok(stats)
