"""
SNP extraction from VCF files - Optimized Version 2.

Architecture:
- Per-chromosome parallel processing (one process per chromosome)
- Each process: extract → BED filter → cluster filter → return results
- NumPy vectorized cluster detection
- No flanking sequence extraction (extracted on-demand in later steps)
- BED statistics collected during extraction (no separate count pass)

Performance targets:
- 15M SNPs: < 15 minutes (vs ~1 hour in v1)
- Memory: ~2GB for 15M SNPs
"""

import logging
import warnings
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from cyvcf2 import VCF

from eprobe.core.result import Result, Ok, Err

logger = logging.getLogger(__name__)


@dataclass
class ChromosomeResult:
    """Result from processing one chromosome."""
    chrom: str
    snps: List[Dict[str, Any]]  # List of SNP dicts
    total_in_vcf: int  # Total SNPs in this chromosome (before BED)
    after_bed: int  # After BED filter
    after_cluster: int  # After cluster filter
    multi_allelic: int
    indels_skipped: int


def detect_clusters_numpy(
    positions: np.ndarray,
    flank: int,
    max_snp: int,
) -> np.ndarray:
    """
    Detect SNP clusters using NumPy vectorization.
    
    Algorithm: For each SNP, count neighbors within ±flank using searchsorted.
    Time complexity: O(n log n) for sorting + O(n) for counting.
    
    Args:
        positions: Sorted 1D numpy array of SNP positions
        flank: Flanking distance for cluster window
        max_snp: Maximum SNPs allowed in window (inclusive)
        
    Returns:
        Boolean mask of SNPs to keep (True = keep, False = filter)
    """
    if len(positions) == 0:
        return np.array([], dtype=bool)
    
    n = len(positions)
    
    # Use searchsorted to find window boundaries efficiently
    # For each position p, find count of SNPs in [p-flank, p+flank]
    left_bounds = np.searchsorted(positions, positions - flank, side='left')
    right_bounds = np.searchsorted(positions, positions + flank, side='right')
    
    # Count SNPs in each window
    window_counts = right_bounds - left_bounds
    
    # Keep SNPs where window count <= max_snp
    keep_mask = window_counts <= max_snp
    
    return keep_mask


def _process_chromosome(
    vcf_path: Path,
    chrom: str,
    bed_regions: Optional[List[Tuple[int, int]]],  # Only regions for this chrom
    cluster_flank: int,
    max_cluster_snp: int,
    cluster_enabled: bool,
) -> ChromosomeResult:
    """
    Process a single chromosome: extract SNPs, apply BED filter, apply cluster filter.
    
    This function runs in a separate process.
    
    Args:
        vcf_path: Path to VCF file
        chrom: Chromosome name
        bed_regions: List of (start, end) tuples for this chromosome only
                    None = no BED filter
                    Empty list = no regions to extract (all filtered)
        cluster_flank: Flanking distance for cluster detection
        max_cluster_snp: Maximum SNPs in cluster window
        cluster_enabled: Whether to apply cluster filtering
        
    Returns:
        ChromosomeResult with SNPs and statistics
    """
    snps = []
    total_in_vcf = 0
    after_bed = 0
    multi_allelic = 0
    indels_skipped = 0
    
    try:
        vcf = VCF(str(vcf_path))
        
        # First, count total variants in the entire chromosome for statistics
        for variant in vcf(chrom):
            total_in_vcf += 1
        
        # Reset VCF for extraction
        vcf.close()
        vcf = VCF(str(vcf_path))
        
        # Suppress cyvcf2 warnings for empty intervals
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="no intervals found")
            
            if bed_regions is not None:
                # BED filter mode: iterate through specified regions
                for start, end in bed_regions:
                    region_str = f"{chrom}:{start}-{end}"
                    for variant in vcf(region_str):
                        if variant.CHROM != chrom:
                            continue
                        
                        # Check if biallelic SNP
                        if len(variant.REF) != 1:
                            indels_skipped += 1
                            continue
                        
                        if not variant.ALT or len(variant.ALT) == 0:
                            continue
                        
                        first_alt = variant.ALT[0]
                        if len(first_alt) != 1 or first_alt not in "ACGT":
                            indels_skipped += 1
                            continue
                        
                        if len(variant.ALT) > 1:
                            multi_allelic += 1
                        
                        # Determine mutation type
                        transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
                        mut_type = 'ts' if (variant.REF.upper(), first_alt.upper()) in transitions else 'tv'
                        
                        snps.append({
                            'chr': variant.CHROM,
                            'pos': variant.POS,
                            'ref': variant.REF,
                            'alt': first_alt,
                            'type': mut_type,
                        })
            else:
                # No BED filter: extract entire chromosome
                for variant in vcf(chrom):
                    total_in_vcf += 1
                    
                    if len(variant.REF) != 1:
                        indels_skipped += 1
                        continue
                    
                    if not variant.ALT or len(variant.ALT) == 0:
                        continue
                    
                    first_alt = variant.ALT[0]
                    if len(first_alt) != 1 or first_alt not in "ACGT":
                        indels_skipped += 1
                        continue
                    
                    if len(variant.ALT) > 1:
                        multi_allelic += 1
                    
                    transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
                    mut_type = 'ts' if (variant.REF.upper(), first_alt.upper()) in transitions else 'tv'
                    
                    snps.append({
                        'chr': variant.CHROM,
                        'pos': variant.POS,
                        'ref': variant.REF,
                        'alt': first_alt,
                        'type': mut_type,
                    })
        
        vcf.close()
        after_bed = len(snps)
        
        # Apply cluster filter using NumPy
        if cluster_enabled and len(snps) > 0:
            # Extract positions as numpy array
            positions = np.array([s['pos'] for s in snps], dtype=np.int64)
            
            # Sort by position (required for cluster detection)
            sort_idx = np.argsort(positions)
            sorted_positions = positions[sort_idx]
            
            # Detect clusters
            keep_mask = detect_clusters_numpy(sorted_positions, cluster_flank, max_cluster_snp)
            
            # Map back to original order and filter
            keep_original_idx = sort_idx[keep_mask]
            snps = [snps[i] for i in sorted(keep_original_idx)]
        
        after_cluster = len(snps)
        
        return ChromosomeResult(
            chrom=chrom,
            snps=snps,
            total_in_vcf=total_in_vcf,
            after_bed=after_bed,
            after_cluster=after_cluster,
            multi_allelic=multi_allelic,
            indels_skipped=indels_skipped,
        )
        
    except Exception as e:
        logger.error(f"Error processing {chrom}: {e}")
        return ChromosomeResult(
            chrom=chrom,
            snps=[],
            total_in_vcf=0,
            after_bed=0,
            after_cluster=0,
            multi_allelic=0,
            indels_skipped=0,
        )


def invert_bed_regions(bed_path: Path, fai_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """
    Generate inverted BED regions grouped by chromosome.
    
    Returns a dict where keys are chromosome names and values are
    lists of (start, end) tuples representing regions OUTSIDE the input BED.
    
    Args:
        bed_path: Input BED file
        fai_path: Reference .fai index file
        
    Returns:
        Dict mapping chromosome -> list of (start, end) regions
    """
    # Read chromosome lengths from .fai
    chrom_lengths = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_lengths[parts[0]] = int(parts[1])
    
    # Read and group BED regions by chromosome
    bed_regions: Dict[str, List[Tuple[int, int]]] = {}
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
    
    # Generate inverted regions by chromosome
    inverted: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, length in chrom_lengths.items():
        inverted[chrom] = []
        
        if chrom not in bed_regions:
            # No BED regions on this chromosome, keep entire chromosome
            inverted[chrom].append((0, length))
            continue
        
        regions = bed_regions[chrom]
        current_pos = 0
        
        for start, end in regions:
            if current_pos < start:
                inverted[chrom].append((current_pos, start))
            current_pos = max(current_pos, end)
        
        if current_pos < length:
            inverted[chrom].append((current_pos, length))
    
    return inverted


def read_bed_regions(bed_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """
    Read BED regions grouped by chromosome.
    
    Args:
        bed_path: Input BED file
        
    Returns:
        Dict mapping chromosome -> list of (start, end) regions
    """
    bed_regions: Dict[str, List[Tuple[int, int]]] = {}
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
    
    return bed_regions


def run_extract(
    vcf_path: Path,
    reference_path: Path,
    output_prefix: Path,
    cluster_flank: int = 60,
    max_cluster_snp: int = 3,
    cluster_mode: bool = True,
    bed_path: Optional[Path] = None,
    bed_mode: str = "keep",
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Extract SNPs from VCF file - Optimized Version 2.
    
    Single-pass extraction with per-chromosome parallelization.
    Each chromosome process handles: extract → BED filter → cluster filter.
    
    Args:
        vcf_path: Input VCF file (.vcf.gz with .tbi index)
        reference_path: Reference genome FASTA (for .fai file)
        output_prefix: Output file prefix
        cluster_flank: Flanking distance for cluster detection
        max_cluster_snp: Maximum SNPs allowed in cluster window
        cluster_mode: Enable cluster filtering
        bed_path: Optional BED file
        bed_mode: 'keep' to keep BED regions, 'remove' to exclude them
        threads: Number of parallel processes
        verbose: Enable verbose logging
        
    Returns:
        Result containing extraction statistics
    """
    if verbose:
        logging.basicConfig(level=logging.DEBUG)
    
    logger.info(f"Starting SNP extraction from {vcf_path}")
    logger.info(f"Using {threads} threads")
    
    # Check indices exist
    fai_path = Path(str(reference_path) + ".fai")
    if not fai_path.exists():
        return Err(f"Reference index not found: {fai_path}")
    
    vcf_index = Path(str(vcf_path) + ".tbi")
    if not vcf_index.exists():
        return Err(f"VCF index not found: {vcf_index}")
    
    # Get chromosome list from VCF
    try:
        vcf = VCF(str(vcf_path))
        chromosomes = vcf.seqnames
        vcf.close()
        logger.info(f"Found {len(chromosomes)} chromosomes in VCF")
    except Exception as e:
        return Err(f"Failed to read VCF: {e}")
    
    # Prepare BED regions by chromosome
    bed_regions_by_chrom: Optional[Dict[str, List[Tuple[int, int]]]] = None
    bed_applied = bed_path is not None
    
    if bed_path:
        if bed_mode == "remove":
            logger.info(f"Generating inverted BED regions (remove mode)...")
            bed_regions_by_chrom = invert_bed_regions(bed_path, fai_path)
            total_regions = sum(len(v) for v in bed_regions_by_chrom.values())
            logger.info(f"Generated {total_regions} inverted regions across {len(bed_regions_by_chrom)} chromosomes")
        else:
            logger.info(f"Loading BED regions (keep mode)...")
            bed_regions_by_chrom = read_bed_regions(bed_path)
            total_regions = sum(len(v) for v in bed_regions_by_chrom.values())
            logger.info(f"Loaded {total_regions} regions across {len(bed_regions_by_chrom)} chromosomes")
    
    # Process chromosomes in parallel
    all_results: List[ChromosomeResult] = []
    
    if threads == 1:
        for chrom in chromosomes:
            chrom_bed = bed_regions_by_chrom.get(chrom, []) if bed_regions_by_chrom else None
            result = _process_chromosome(
                vcf_path, chrom, chrom_bed,
                cluster_flank, max_cluster_snp, cluster_mode
            )
            all_results.append(result)
            if result.after_cluster > 0:
                logger.info(f"✓ {chrom}: {result.after_cluster} SNPs (BED: {result.after_bed}, cluster removed: {result.after_bed - result.after_cluster})")
    else:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {}
            for chrom in chromosomes:
                chrom_bed = bed_regions_by_chrom.get(chrom, []) if bed_regions_by_chrom else None
                future = executor.submit(
                    _process_chromosome,
                    vcf_path, chrom, chrom_bed,
                    cluster_flank, max_cluster_snp, cluster_mode
                )
                futures[future] = chrom
            
            for future in as_completed(futures):
                chrom = futures[future]
                try:
                    result = future.result()
                    all_results.append(result)
                    if result.after_cluster > 0:
                        logger.info(f"✓ {chrom}: {result.after_cluster} SNPs")
                except Exception as e:
                    logger.error(f"Failed to process {chrom}: {e}")
    
    # Aggregate results
    all_snps = []
    total_in_vcf = 0
    total_after_bed = 0
    total_after_cluster = 0
    total_multi_allelic = 0
    total_indels = 0
    
    for result in all_results:
        all_snps.extend(result.snps)
        total_in_vcf += result.total_in_vcf
        total_after_bed += result.after_bed
        total_after_cluster += result.after_cluster
        total_multi_allelic += result.multi_allelic
        total_indels += result.indels_skipped
    
    # Sort by chromosome and position
    all_snps.sort(key=lambda x: (x['chr'], x['pos']))
    
    # Save output
    output_path = Path(str(output_prefix) + ".snps.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    df = pd.DataFrame(all_snps)
    try:
        df.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Saved {len(all_snps)} SNPs to {output_path}")
    except Exception as e:
        return Err(f"Failed to save output: {e}")
    
    # Build statistics
    stats = {
        "total_in_vcf": total_in_vcf,
        "after_bed": total_after_bed,
        "after_cluster": total_after_cluster,
        "final_snp_count": len(all_snps),
        "multi_allelic": total_multi_allelic,
        "indels_skipped": total_indels,
        "output_file": str(output_path),
        # BED filter statistics
        "bed_applied": bed_applied,
        "bed_mode": bed_mode if bed_applied else None,
        "bed_kept": total_after_bed if bed_applied else None,
        "bed_removed": (total_in_vcf - total_after_bed) if bed_applied else None,
        # Cluster filter statistics
        "cluster_enabled": cluster_mode,
        "cluster_removed": total_after_bed - total_after_cluster if cluster_mode else 0,
    }
    
    return Ok(stats)
