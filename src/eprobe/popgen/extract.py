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


def filter_by_bed_regions_numpy(
    positions: np.ndarray,
    bed_starts: np.ndarray,
    bed_ends: np.ndarray,
    mode: str = "keep"
) -> np.ndarray:
    """
    Vectorized BED region filtering using NumPy.
    
    Algorithm: Use searchsorted to find which BED interval each SNP might overlap,
    then check if the SNP position is actually within that interval.
    
    Args:
        positions: 1D numpy array of SNP positions (must be sorted for best performance)
        bed_starts: Sorted 1D numpy array of BED region start positions (0-based)
        bed_ends: Sorted 1D numpy array of BED region end positions (0-based, exclusive)
        mode: "keep" = keep SNPs inside BED regions, "remove" = keep SNPs outside BED regions
        
    Returns:
        Boolean mask (True = keep this SNP)
    """
    if len(positions) == 0:
        return np.array([], dtype=bool)
    
    if len(bed_starts) == 0:
        # No BED regions: if mode=keep, keep nothing; if mode=remove, keep all
        if mode == "keep":
            return np.zeros(len(positions), dtype=bool)
        else:
            return np.ones(len(positions), dtype=bool)
    
    # For each SNP position, find the index of the first BED region whose end > position
    # This is the candidate region that might contain the SNP
    candidate_idx = np.searchsorted(bed_ends, positions, side='right')
    
    # Check if the SNP is within the candidate BED region
    # A SNP at position p is in region i if: bed_starts[i] <= p < bed_ends[i]
    valid_idx = candidate_idx < len(bed_starts)
    
    # Initialize mask
    in_bed = np.zeros(len(positions), dtype=bool)
    
    # For valid indices, check if position >= bed_starts[candidate_idx]
    in_bed[valid_idx] = positions[valid_idx] >= bed_starts[candidate_idx[valid_idx]]
    
    # Return based on mode
    if mode == "keep":
        return in_bed
    else:  # mode == "remove"
        return ~in_bed


def classify_mutations_numpy(refs: np.ndarray, alts: np.ndarray) -> np.ndarray:
    """
    Vectorized mutation type classification (transition vs transversion).
    
    Args:
        refs: 1D numpy array of reference alleles (single characters)
        alts: 1D numpy array of alternate alleles (single characters)
        
    Returns:
        1D numpy array of mutation types ('ts' or 'tv')
    """
    # Transition pairs: A<->G, C<->T
    # Convert to uppercase for comparison
    refs_upper = np.char.upper(refs.astype(str))
    alts_upper = np.char.upper(alts.astype(str))
    
    # Check for transitions
    is_transition = (
        ((refs_upper == 'A') & (alts_upper == 'G')) |
        ((refs_upper == 'G') & (alts_upper == 'A')) |
        ((refs_upper == 'C') & (alts_upper == 'T')) |
        ((refs_upper == 'T') & (alts_upper == 'C'))
    )
    
    # Create result array
    result = np.where(is_transition, 'ts', 'tv')
    
    return result


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
    Uses vectorized NumPy operations for filtering and classification.
    
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
    total_in_vcf = 0
    multi_allelic = 0
    indels_skipped = 0
    
    try:
        vcf = VCF(str(vcf_path))
        
        # ========== STEP 1: Batch read all variants into arrays ==========
        # Pre-allocate lists for batch collection (more efficient than individual appends)
        raw_positions = []
        raw_refs = []
        raw_alts = []
        raw_alt_counts = []  # Track multi-allelic
        
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="no intervals found")
            
            for variant in vcf(chrom):
                total_in_vcf += 1
                
                # Collect raw data - filtering done vectorized later
                ref = variant.REF
                alt_list = variant.ALT if variant.ALT else []
                first_alt = alt_list[0] if alt_list else ''
                
                raw_positions.append(variant.POS)
                raw_refs.append(ref)
                raw_alts.append(first_alt)
                raw_alt_counts.append(len(alt_list))
        
        vcf.close()
        
        if total_in_vcf == 0:
            return ChromosomeResult(
                chrom=chrom, snps=[], total_in_vcf=0, after_bed=0,
                after_cluster=0, multi_allelic=0, indels_skipped=0,
            )
        
        # Convert to NumPy arrays for vectorized operations
        positions = np.array(raw_positions, dtype=np.int64)
        refs = np.array(raw_refs, dtype=object)
        alts = np.array(raw_alts, dtype=object)
        alt_counts = np.array(raw_alt_counts, dtype=np.int32)
        
        # ========== STEP 2: Vectorized SNP filtering ==========
        # Filter 1: REF must be single nucleotide
        ref_lens = np.array([len(r) for r in refs], dtype=np.int32)
        ref_valid = ref_lens == 1
        
        # Filter 2: ALT must be single nucleotide and valid base
        valid_bases = {'A', 'C', 'G', 'T'}
        alt_valid = np.array([
            len(a) == 1 and a in valid_bases for a in alts
        ], dtype=bool)
        
        # Combined filter mask
        snp_valid = ref_valid & alt_valid
        
        # Count statistics
        indels_skipped = int(np.sum(~snp_valid))
        multi_allelic = int(np.sum(alt_counts[snp_valid] > 1))
        
        # Apply SNP filter
        positions = positions[snp_valid]
        refs = refs[snp_valid]
        alts = alts[snp_valid]
        
        # ========== STEP 3: Vectorized BED filtering ==========
        if bed_regions is not None:
            if len(bed_regions) == 0:
                # Empty BED regions list = no regions to keep
                positions = np.array([], dtype=np.int64)
                refs = np.array([], dtype=object)
                alts = np.array([], dtype=object)
            else:
                # Convert BED regions to sorted NumPy arrays
                bed_starts = np.array([start for start, end in bed_regions], dtype=np.int64)
                bed_ends = np.array([end for start, end in bed_regions], dtype=np.int64)
                
                # Sort BED regions by start position (required for searchsorted)
                sort_idx = np.argsort(bed_starts)
                bed_starts = bed_starts[sort_idx]
                bed_ends = bed_ends[sort_idx]
                
                # Use vectorized BED filtering - mode="keep" because bed_regions 
                # are already inverted (for remove_bed) or direct (for keep_bed)
                bed_mask = filter_by_bed_regions_numpy(positions, bed_starts, bed_ends, mode="keep")
                
                # Apply BED filter
                positions = positions[bed_mask]
                refs = refs[bed_mask]
                alts = alts[bed_mask]
        
        after_bed = len(positions)
        
        # ========== STEP 4: Vectorized mutation type classification ==========
        if len(positions) > 0:
            mut_types = classify_mutations_numpy(refs, alts)
        else:
            mut_types = np.array([], dtype=object)
        
        # ========== STEP 5: Vectorized cluster filtering ==========
        if cluster_enabled and len(positions) > 0:
            # Sort by position for cluster detection
            sort_idx = np.argsort(positions)
            sorted_positions = positions[sort_idx]
            
            # Detect clusters
            keep_mask = detect_clusters_numpy(sorted_positions, cluster_flank, max_cluster_snp)
            
            # Apply cluster filter (keep sorted order for output)
            keep_indices = sort_idx[keep_mask]
            keep_indices = np.sort(keep_indices)  # Restore original order
            
            positions = positions[keep_indices]
            refs = refs[keep_indices]
            alts = alts[keep_indices]
            mut_types = mut_types[keep_indices]
        
        after_cluster = len(positions)
        
        # ========== STEP 6: Build output SNP list ==========
        snps = [
            {
                'chr': chrom,
                'pos': int(pos),
                'ref': str(ref),
                'alt': str(alt),
                'type': str(mt),
            }
            for pos, ref, alt, mt in zip(positions, refs, alts, mut_types)
        ]
        
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
        import traceback
        logger.error(traceback.format_exc())
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
