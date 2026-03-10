"""
SNP selection module.

Implements strategies for selecting SNPs from filtered candidates:
  - Uniform distribution: equal spacing across genome
  - Random sampling: random selection with optional stratification
  - Weighted: selection based on biophysical scores
  - Priority: prioritize SNPs in specific genomic regions (BED)

This module corresponds to the original SNP_subsampler.py functionality.
"""

import logging
import random
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple, Set, Union
from dataclasses import dataclass, field
from collections import defaultdict
from enum import Enum

import numpy as np
import pandas as pd

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNP, SNPDataFrame

logger = logging.getLogger(__name__)


class SelectionStrategy(Enum):
    """SNP selection strategies."""
    UNIFORM = "uniform"
    RANDOM = "random"
    WEIGHTED = "weighted"
    PRIORITY = "priority"


# Biophysical columns that can be used for weighted selection
BIOPHYSICAL_COLUMNS = ['gc', 'tm', 'complexity', 'hairpin', 'dimer']


def _load_file_worker(path: Path):
    """Worker that loads a single SNP TSV file. Returns (path, snp_list, error)."""
    snp_df_result = SNPDataFrame.from_tsv(path)
    if snp_df_result.is_err():
        return path, None, snp_df_result.unwrap_err()
    snp_list = snp_df_result.unwrap().to_snps()
    return path, snp_list, None


def merge_snp_files(
    input_paths: List[Path],
    merge_mode: str = "intersection",
    threads: int = 1,
) -> Result[List[SNP], str]:
    """
    Merge SNPs from multiple input TSV files.
    
    Args:
        input_paths: List of SNP TSV file paths
        merge_mode: Merge operation ("intersection", "union", "difference", "symmetric_diff")
        threads: Number of threads for parallel file loading (default: 1)
        
    Returns:
        Result containing merged SNP list
        
    Merge modes:
        - intersection: Keep SNPs present in ALL files
        - union: Combine all unique SNPs
        - difference: Keep SNPs in first file but NOT in others
        - symmetric_diff: Keep SNPs present in exactly one file
    """
    if not input_paths:
        return Err("No input files provided")
    
    if len(input_paths) == 1:
        # Single input, just load and return
        snp_df_result = SNPDataFrame.from_tsv(input_paths[0])
        if snp_df_result.is_err():
            return Err(f"Failed to load {input_paths[0]}: {snp_df_result.unwrap_err()}")
        return Ok(snp_df_result.unwrap().to_snps())
    
    # Load all input files (parallel when threads > 1)
    all_snps_list = [None] * len(input_paths)  # preserve order
    
    n_workers = min(threads, len(input_paths)) if threads > 1 else 1
    if n_workers > 1:
        logger.info(f"Loading {len(input_paths)} files with {n_workers} threads")
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            future_to_idx = {
                executor.submit(_load_file_worker, path): idx
                for idx, path in enumerate(input_paths)
            }
            for future in as_completed(future_to_idx):
                idx = future_to_idx[future]
                path, snp_list, err = future.result()
                if err:
                    return Err(f"Failed to load {path}: {err}")
                logger.info(f"Loaded {len(snp_list)} SNPs from {path.name}")
                all_snps_list[idx] = snp_list
    else:
        for idx, path in enumerate(input_paths):
            snp_df_result = SNPDataFrame.from_tsv(path)
            if snp_df_result.is_err():
                return Err(f"Failed to load {path}: {snp_df_result.unwrap_err()}")
            snp_list = snp_df_result.unwrap().to_snps()
            logger.info(f"Loaded {len(snp_list)} SNPs from {path.name}")
            all_snps_list[idx] = snp_list
    
    # Build SNP ID sets for set operations
    all_snp_sets = [
        set(f"{snp.chrom}:{snp.pos}_{snp.ref}_{snp.alt}" for snp in snp_list)
        for snp_list in all_snps_list
    ]
    
    # Perform merge operation
    if merge_mode == "intersection":
        # Keep SNPs in ALL files
        result_set = all_snp_sets[0]
        for snp_set in all_snp_sets[1:]:
            result_set = result_set.intersection(snp_set)
    elif merge_mode == "union":
        # Combine all unique SNPs
        result_set = set()
        for snp_set in all_snp_sets:
            result_set = result_set.union(snp_set)
    elif merge_mode == "difference":
        # Only in first file
        result_set = all_snp_sets[0]
        for snp_set in all_snp_sets[1:]:
            result_set = result_set.difference(snp_set)
    elif merge_mode == "symmetric_diff":
        # In exactly one file
        result_set = all_snp_sets[0]
        for snp_set in all_snp_sets[1:]:
            result_set = result_set.symmetric_difference(snp_set)
    else:
        return Err(f"Unknown merge mode: {merge_mode}")
    
    # Filter original SNPs to keep only those in result_set
    result_snps = []
    for snps in all_snps_list:
        for snp in snps:
            snp_id = f"{snp.chrom}:{snp.pos}_{snp.ref}_{snp.alt}"
            if snp_id in result_set:
                result_snps.append(snp)
                result_set.discard(snp_id)  # Avoid duplicates
    
    logger.info(f"After {merge_mode} merge: {len(result_snps)} SNPs")
    return Ok(result_snps)


@dataclass
class SelectionConfig:
    """Configuration for SNP selection."""
    target_count: int = 10000
    window_size: int = 10000
    strategy: str = "random"
    seed: int = 42
    # Weights for biophysical features (gc, tm, complexity, hairpin, dimer)
    weights: Optional[List[float]] = None
    # Priority regions from BED file
    priority_regions: Optional[Set[Tuple[str, int, int]]] = None
    # Chromosomes to include (None = all)
    chromosomes: Optional[List[str]] = None


def calculate_window_index(pos: int, window_size: int) -> int:
    """Calculate window index for a genomic position."""
    return pos // window_size


def parse_bed_file(bed_path: Path) -> Result[Set[Tuple[str, int, int]], str]:
    """
    Parse BED file to extract priority regions.
    
    Args:
        bed_path: Path to BED file
        
    Returns:
        Set of (chrom, start, end) tuples
    """
    if not bed_path.exists():
        return Err(f"BED file not found: {bed_path}")
    
    regions = set()
    try:
        with open(bed_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('track'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                
                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    regions.add((chrom, start, end))
                except ValueError:
                    logger.warning(f"Invalid BED line {line_num}: {line[:50]}")
                    continue
        
        logger.info(f"Loaded {len(regions)} priority regions from BED file")
        return Ok(regions)
    
    except Exception as e:
        return Err(f"Failed to parse BED file: {e}")


def snp_in_priority_region(
    snp: SNP, 
    regions: Set[Tuple[str, int, int]]
) -> bool:
    """Check if SNP falls within any priority region."""
    for chrom, start, end in regions:
        if snp.chrom == chrom and start <= snp.pos <= end:
            return True
    return False


_BIOPHYSICAL_LABELS = {
    'gc': 'GC (%)',
    'tm': 'Tm (°C)',
    'complexity': 'Complexity',
    'hairpin': 'Hairpin',
    'dimer': 'Dimer',
}


def _log_biophysical_comparison(
    input_df: pd.DataFrame,
    selected_df: pd.DataFrame,
) -> None:
    """Log before/after biophysical stats after SNP selection."""
    available = [c for c in BIOPHYSICAL_COLUMNS if c in input_df.columns and c in selected_df.columns]
    if not available:
        return

    col_w = 10
    divider = '-' * (18 + col_w * 3 + 3)
    logger.info("Biophysical metrics: input pool → selected set")
    logger.info(divider)
    logger.info(f"  {'Metric':<16} {'Input':>{col_w}} {'Selected':>{col_w}} {'Δ (mean)':>{col_w}}")
    logger.info(divider)
    for col in available:
        label = _BIOPHYSICAL_LABELS.get(col, col)
        b = float(input_df[col].dropna().mean())
        a = float(selected_df[col].dropna().mean())
        d = a - b
        sign = '+' if d >= 0 else ''
        logger.info(f"  {label:<16} {b:>{col_w}.3f} {a:>{col_w}.3f} {sign}{d:>{col_w - 1}.3f}")
    logger.info(divider)


def _write_biophysical_summary(
    input_df: pd.DataFrame,
    selected_df: pd.DataFrame,
    output_prefix: Path,
    stats: Dict[str, Any],
) -> Path:
    """Write before/after biophysical distribution summary to {output_prefix}.select_summary.txt."""
    summary_path = Path(str(output_prefix) + ".select_summary.txt")
    available = [c for c in BIOPHYSICAL_COLUMNS if c in input_df.columns and c in selected_df.columns]

    with open(summary_path, 'w') as f:
        f.write("eProbe SNP Selection Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Strategy:        {stats.get('strategy', '-')}\n")
        f.write(f"Window size:     {stats.get('window_size', 0):,} bp\n")
        f.write(f"Input SNPs:      {stats.get('initial_count', 0):,}\n")
        f.write(f"Selected SNPs:   {stats.get('selected', 0):,}\n")
        f.write(f"Target count:    {stats.get('target_count', 0):,}\n")
        coverage = stats.get('coverage', {})
        f.write(f"Windows covered: {coverage.get('windows_covered', '-')}\n")
        f.write(f"Chromosomes:     {coverage.get('chromosomes', '-')}\n")
        f.write("\n")

        if available:
            f.write("Biophysical Metrics: Before \u2192 After Selection\n")
            f.write("-" * 50 + "\n\n")
            for col in available:
                label = _BIOPHYSICAL_LABELS.get(col, col)
                in_s = input_df[col].dropna()
                sel_s = selected_df[col].dropna()
                d_mean = float(sel_s.mean()) - float(in_s.mean())
                sign = '+' if d_mean >= 0 else ''
                f.write(f"{label}:\n")
                f.write(f"  Input:    Mean={float(in_s.mean()):.3f}  Std={float(in_s.std()):.3f}  "
                        f"Median={float(in_s.median()):.3f}  Min={float(in_s.min()):.3f}  Max={float(in_s.max()):.3f}\n")
                f.write(f"  Selected: Mean={float(sel_s.mean()):.3f}  Std={float(sel_s.std()):.3f}  "
                        f"Median={float(sel_s.median()):.3f}  Min={float(sel_s.min()):.3f}  Max={float(sel_s.max()):.3f}\n")
                f.write(f"  \u0394 Mean:   {sign}{d_mean:.3f}\n\n")

    logger.info(f"Selection summary saved to {summary_path}")
    return summary_path


# Default target values for biophysical scoring
# [gc%, tm°C, complexity, hairpin, dimer]
DEFAULT_TARGETS = [50.0, 70.0, 0.0, 0.0, 0.0]

# Maximum plausible range for each metric (used for normalization)
# score = 1 - |value - target| / range  →  clamped to [0, 1]
_NORMALIZE_RANGE = [50.0, 30.0, 5.0, None, None]  # None = data-driven


def calculate_snp_score(
    snp: SNP,
    weights: List[float],
    targets: Optional[List[float]] = None,
    probe_sequences: Optional[Dict[str, str]] = None,
) -> float:
    """
    Calculate weighted biophysical score for a SNP.
    
    Uses tags from SNP if available, or uses default values.
    For actual biophysical scoring, tags should be added during filtering.
    
    Args:
        snp: SNP object
        weights: [gc_weight, tm_weight, complexity_weight, hairpin_weight, dimer_weight]
        targets: Target values [gc, tm, complexity, hairpin, dimer]. Default: [50,70,0,0,0]
        probe_sequences: Optional dict of SNP ID -> probe sequence
        
    Returns:
        Weighted score (higher = better)
    """
    if weights is None or len(weights) != 5:
        weights = [1.0, 1.0, 1.0, 1.0, 1.0]
    if targets is None:
        targets = list(DEFAULT_TARGETS)
    
    tags = getattr(snp, 'tags', None)
    
    if tags and isinstance(tags, dict):
        gc = tags.get('gc', targets[0])
        tm = tags.get('tm', targets[1])
        complexity = tags.get('complexity', targets[2])
        hairpin = tags.get('hairpin', targets[3])
        dimer = tags.get('dimer', targets[4])
    else:
        gc = targets[0]
        tm = targets[1]
        complexity = targets[2]
        hairpin = targets[3]
        dimer = targets[4]
    
    gc_score = max(0, 1.0 - abs(gc - targets[0]) / _NORMALIZE_RANGE[0])
    tm_score = max(0, 1.0 - abs(tm - targets[1]) / _NORMALIZE_RANGE[1])
    complexity_score = max(0, 1.0 - abs(complexity - targets[2]) / _NORMALIZE_RANGE[2])
    hairpin_score = max(0, 1.0 - hairpin / 50.0) if targets[3] == 0 else max(0, 1.0 - abs(hairpin - targets[3]) / 50.0)
    dimer_score = max(0, 1.0 - dimer) if targets[4] == 0 else max(0, 1.0 - abs(dimer - targets[4]))
    
    total_score = (
        weights[0] * gc_score +
        weights[1] * tm_score +
        weights[2] * complexity_score +
        weights[3] * hairpin_score +
        weights[4] * dimer_score
    )
    
    return total_score


def select_uniform(
    snps: List[SNP],
    target_count: int,
    window_size: int = 100000,
    seed: int = 42,
) -> Result[List[SNP], str]:
    """
    Select SNPs with uniform genomic distribution.
    
    Divides the genome into windows and selects SNPs proportionally
    from each window to maintain even coverage.
    
    Algorithm:
    1. Group SNPs by chromosome and window
    2. Calculate SNPs per window based on target count
    3. Random sample within each window
    
    Args:
        snps: Input SNP list
        target_count: Target number of SNPs to select
        window_size: Window size for uniform distribution
        seed: Random seed for reproducibility
        
    Returns:
        Result containing selected SNP list
    """
    if not snps:
        return Ok([])
    
    random.seed(seed)
    np.random.seed(seed)
    
    logger.info(f"Uniform selection: {len(snps)} → {target_count} SNPs")
    logger.info(f"Window size: {window_size:,} bp")
    
    # Group SNPs by chromosome and window
    window_snps: Dict[Tuple[str, int], List[SNP]] = defaultdict(list)
    
    for snp in snps:
        window_idx = calculate_window_index(snp.pos, window_size)
        key = (snp.chrom, window_idx)
        window_snps[key].append(snp)
    
    total_windows = len(window_snps)
    logger.info(f"Total windows with SNPs: {total_windows}")
    
    # Calculate base SNPs per window
    if total_windows == 0:
        return Ok([])
    
    base_per_window = target_count // total_windows
    remainder = target_count % total_windows
    
    # Sort windows for deterministic selection
    sorted_windows = sorted(window_snps.keys())
    
    # Select SNPs from each window
    selected_snps: List[SNP] = []
    
    for i, window_key in enumerate(sorted_windows):
        window_candidates = window_snps[window_key]
        
        # Calculate how many to select from this window
        n_select = base_per_window + (1 if i < remainder else 0)
        
        if n_select <= 0:
            continue
        
        # Sample from window
        if len(window_candidates) <= n_select:
            selected_snps.extend(window_candidates)
        else:
            sampled = random.sample(window_candidates, n_select)
            selected_snps.extend(sampled)
    
    # Sort by chromosome and position
    selected_snps.sort(key=lambda s: (s.chrom, s.pos))
    
    logger.info(f"Selected {len(selected_snps)} SNPs from {total_windows} windows")
    
    return Ok(selected_snps)


def select_random(
    snps: List[SNP],
    target_count: int,
    seed: int = 42,
) -> Result[List[SNP], str]:
    """
    Randomly sample SNPs.
    
    Simple random sampling without considering genomic distribution.
    
    Args:
        snps: Input SNP list
        target_count: Target number of SNPs to select
        seed: Random seed for reproducibility
        
    Returns:
        Result containing selected SNP list
    """
    if not snps:
        return Ok([])
    
    random.seed(seed)
    
    logger.info(f"Random selection: {len(snps)} → {target_count} SNPs")
    
    if len(snps) <= target_count:
        logger.warning(f"Input count ({len(snps)}) <= target ({target_count}), returning all")
        return Ok(snps.copy())
    
    selected = random.sample(snps, target_count)
    
    # Sort by chromosome and position
    selected.sort(key=lambda s: (s.chrom, s.pos))
    
    logger.info(f"Selected {len(selected)} SNPs randomly")
    
    return Ok(selected)


def select_by_priority(
    snps: List[SNP],
    target_count: int,
    priority_regions: Set[Tuple[str, int, int]],
    window_size: int = 10000,
    seed: int = 42,
    weights: Optional[List[float]] = None,
    snp_df: Optional[pd.DataFrame] = None,
    targets: Optional[List[float]] = None,
) -> Result[Union[List[SNP], Tuple[List[SNP], pd.DataFrame]], str]:
    """
    Select SNPs with priority for specific genomic regions.
    
    Prioritizes SNPs that fall within priority regions (e.g., exons).
    Uses window-based selection to maintain uniform distribution.
    If weights are provided and snp_df has biophysical columns, uses
    weighted selection within priority regions.
    
    Args:
        snps: Input SNP list
        target_count: Target number of SNPs to select
        priority_regions: Set of (chrom, start, end) priority regions
        window_size: Window size for distribution
        seed: Random seed for reproducibility
        weights: Optional biophysical weights [gc, tm, complexity, hairpin, dimer]
        snp_df: Optional DataFrame with biophysical columns for weighted selection
        targets: Target values [gc, tm, complexity, hairpin, dimer]. Default: [50,70,0,0,0]
        
    Returns:
        Result containing selected SNP list (or tuple with DataFrame if weighted)
    """
    if not snps:
        return Ok([])
    
    random.seed(seed)
    
    logger.info(f"Priority selection: {len(snps)} → {target_count} SNPs")
    logger.info(f"Priority regions: {len(priority_regions)}")
    
    # Check if we can use weighted selection
    use_weighted = False
    if weights is not None and snp_df is not None:
        missing_cols = [col for col in BIOPHYSICAL_COLUMNS if col not in snp_df.columns]
        if not missing_cols:
            use_weighted = True
            logger.info("Using weighted selection within priority regions")
        else:
            logger.warning(f"Cannot use weighted selection: missing columns {missing_cols}")
    
    # Separate SNPs into priority and non-priority
    priority_snps = []
    other_snps = []
    priority_indices = []
    other_indices = []
    
    for i, snp in enumerate(snps):
        if snp_in_priority_region(snp, priority_regions):
            priority_snps.append(snp)
            priority_indices.append(i)
        else:
            other_snps.append(snp)
            other_indices.append(i)
    
    logger.info(f"SNPs in priority regions: {len(priority_snps)}")
    logger.info(f"SNPs outside priority regions: {len(other_snps)}")
    
    # Select from priority regions first
    selected = []
    selected_df = None
    
    if len(priority_snps) >= target_count:
        # All from priority regions
        if use_weighted:
            # Filter DataFrame to priority SNPs only
            priority_df = snp_df.iloc[priority_indices].copy()
            result = select_weighted(priority_snps, target_count, weights, window_size, seed, snp_df=priority_df, targets=targets)
            if result.is_err():
                return result
            selected, selected_df = result.unwrap()
        else:
            result = select_uniform(priority_snps, target_count, window_size, seed)
            if result.is_err():
                return result
            selected = result.unwrap()
    else:
        # All priority + selection from others
        if use_weighted:
            # Use all priority SNPs with weighted selection
            priority_df = snp_df.iloc[priority_indices].copy()
            if len(priority_snps) > 0:
                # Select best from priority using weights
                result = select_weighted(priority_snps, len(priority_snps), weights, window_size, seed, snp_df=priority_df, targets=targets)
                if result.is_err():
                    return result
                selected, selected_df = result.unwrap()
            
            remaining = target_count - len(selected)
            
            if remaining > 0 and other_snps:
                # Also use weighted for remaining from other regions
                other_df = snp_df.iloc[other_indices].copy()
                result = select_weighted(other_snps, remaining, weights, window_size, seed, snp_df=other_df, targets=targets)
                if result.is_err():
                    return result
                other_selected, other_df_selected = result.unwrap()
                selected.extend(other_selected)
                if selected_df is not None:
                    selected_df = pd.concat([selected_df, other_df_selected], ignore_index=True)
                else:
                    selected_df = other_df_selected
        else:
            selected.extend(priority_snps)
            remaining = target_count - len(priority_snps)
            
            if remaining > 0 and other_snps:
                result = select_uniform(other_snps, remaining, window_size, seed)
                if result.is_err():
                    return result
                selected.extend(result.unwrap())
    
    # Sort by chromosome and position
    selected.sort(key=lambda s: (s.chrom, s.pos))
    if selected_df is not None:
        selected_df = selected_df.sort_values(['chr', 'pos']).reset_index(drop=True)
    
    logger.info(f"Selected {len(selected)} SNPs ({len([s for s in selected if snp_in_priority_region(s, priority_regions)])} in priority regions)")
    
    if use_weighted and selected_df is not None:
        return Ok((selected, selected_df))
    return Ok(selected)


def select_weighted(
    snps: List[SNP],
    target_count: int,
    weights: List[float],
    window_size: int = 10000,
    seed: int = 42,
    snp_df: Optional[pd.DataFrame] = None,
    threads: int = 1,
    targets: Optional[List[float]] = None,
) -> Result[List[SNP], str]:
    """
    Select SNPs based on weighted biophysical scores within windows.
    
    REQUIRES biophysical tags (gc, tm, complexity, hairpin, dimer) in the DataFrame.
    Divides genome into windows and selects SNPs with best biophysical scores.
    
    Algorithm:
    1. Validate biophysical columns exist
    2. Calculate weighted composite score for each SNP (distance to targets)
    3. Group by window (chrom + pos // window_size)
    4. Select top-scoring SNP(s) from each window
    
    Args:
        snps: Input SNP list (for output compatibility)
        target_count: Target number of SNPs to select
        weights: [gc_weight, tm_weight, complexity_weight, hairpin_weight, dimer_weight]
        window_size: Window size for selection (default: 10kb)
        seed: Random seed for tie-breaking
        snp_df: DataFrame with biophysical columns (required for weighted selection)
        targets: Target values [gc, tm, complexity, hairpin, dimer]. Default: [50,70,0,0,0]
        
    Returns:
        Result containing selected SNP list, or error if biophysical columns missing
    """
    if not snps:
        return Ok([])
    
    if snp_df is None or snp_df.empty:
        return Err("Weighted selection requires DataFrame with biophysical columns")
    
    # Validate required biophysical columns
    missing_cols = [col for col in BIOPHYSICAL_COLUMNS if col not in snp_df.columns]
    if missing_cols:
        return Err(f"Weighted selection requires biophysical columns: {missing_cols}. "
                   f"Run filter with biophysical filtering first.")
    
    random.seed(seed)
    np.random.seed(seed)
    
    logger.info(f"Weighted selection: {len(snp_df)} → {target_count} SNPs")
    logger.info(f"Weights: gc={weights[0]}, tm={weights[1]}, complexity={weights[2]}, "
                f"hairpin={weights[3]}, dimer={weights[4]}")
    
    if targets is None:
        targets = list(DEFAULT_TARGETS)
    logger.info(f"Targets: gc={targets[0]}, tm={targets[1]}, complexity={targets[2]}, "
                f"hairpin={targets[3]}, dimer={targets[4]}")
    
    df = snp_df.copy()
    
    # === WINDOW-BASED RANK SCORING ===
    # Compute |value - target| as raw distance for each metric.
    # Within each window, convert distances to ranks and normalise to [0,1].
    # This ensures every metric contributes equally regardless of its natural
    # scale (e.g. GC% vs complexity) — only the relative ordering within a
    # window matters.
    
    dist_cols = []
    for col, tgt in zip(BIOPHYSICAL_COLUMNS, targets):
        dcol = f'_dist_{col}'
        df[dcol] = np.abs(df[col] - tgt)
        dist_cols.append(dcol)
    
    # Window key
    df['window'] = df['pos'] // window_size
    df['window_key'] = df['chr'].astype(str) + '_' + df['window'].astype(str)
    
    total_windows = df['window_key'].nunique()
    logger.info(f"Total windows with SNPs: {total_windows}")
    
    if total_windows == 0:
        return Ok([])
    
    # Rank distances within each window (ascending: smallest distance = rank 1 = best).
    # Normalise ranks to [0,1]: rank_score = 1 - (rank - 1) / (n - 1).
    # When a window has only 1 SNP, rank_score = 1.0 for all metrics.
    score_cols = [f'{col}_score' for col in BIOPHYSICAL_COLUMNS]
    for dcol, scol in zip(dist_cols, score_cols):
        # rank with method='average' handles ties symmetrically
        df[scol] = df.groupby('window_key')[dcol].rank(method='average', ascending=True)
        window_counts = df.groupby('window_key')[dcol].transform('count')
        df[scol] = 1.0 - (df[scol] - 1.0) / (window_counts - 1.0).replace(0, np.nan).fillna(1.0)
    
    # Drop distance columns
    df.drop(columns=dist_cols, inplace=True)
    
    # Weighted composite score from rank scores
    df['weighted_score'] = sum(
        w * df[sc] for w, sc in zip(weights, score_cols)
    )
    
    # Calculate SNPs per window
    base_per_window = max(1, target_count // total_windows)
    
    # Select top N from each window based on weighted rank score.
    def _process_windows(chrom_df: pd.DataFrame) -> List:
        indices = []
        for _, group in chrom_df.groupby('window_key'):
            sorted_group = group.nlargest(base_per_window, 'weighted_score')
            indices.extend(sorted_group.index.tolist())
        return indices
    
    selected_indices: List = []
    n_workers = min(threads, df['chr'].nunique()) if threads > 1 else 1
    if n_workers > 1:
        logger.info(f"Window selection using {n_workers} threads (by chromosome)")
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = [
                executor.submit(_process_windows, chrom_df)
                for _, chrom_df in df.groupby('chr', sort=False)
            ]
            for future in futures:
                selected_indices.extend(future.result())
    else:
        selected_indices = _process_windows(df)
    
    selected_df = df.loc[selected_indices]
    
    # If we have more than target, randomly sample
    if len(selected_df) > target_count:
        selected_df = selected_df.sample(n=target_count, random_state=seed)
    
    # Sort by chromosome and position
    selected_df = selected_df.sort_values(['chr', 'pos'])
    
    # Remove temporary columns, keep biophysical columns
    cols_to_drop = ['gc_score', 'tm_score', 'complexity_score', 'hairpin_score', 
                    'dimer_score', 'weighted_score', 'window', 'window_key']
    selected_df = selected_df.drop(columns=cols_to_drop, errors='ignore')
    
    # Convert back to SNP list
    selected_snps = []
    for _, row in selected_df.iterrows():
        snp = SNP(
            chrom=row['chr'],
            pos=int(row['pos']),
            ref=row['ref'],
            alt=row['alt'],
            mutation_type=row['type']
        )
        selected_snps.append(snp)
    
    logger.info(f"Selected {len(selected_snps)} SNPs by weighted score")
    logger.info(f"Score range: {df['weighted_score'].min():.3f} - {df['weighted_score'].max():.3f}")
    
    # Return tuple: (SNPs, DataFrame with biophysical cols)
    return Ok((selected_snps, selected_df))


def select_chromosomes(
    snps: List[SNP],
    chromosomes: List[str],
) -> Result[List[SNP], str]:
    """
    Filter SNPs to specific chromosomes.
    
    Args:
        snps: Input SNP list
        chromosomes: List of chromosome names to keep
        
    Returns:
        Result containing filtered SNP list
    """
    if not snps:
        return Ok([])
    
    chrom_set = set(chromosomes)
    filtered = [snp for snp in snps if snp.chrom in chrom_set]
    
    logger.info(f"Chromosome filter: {len(snps)} → {len(filtered)} SNPs")
    logger.info(f"Kept chromosomes: {chromosomes}")
    
    return Ok(filtered)


def run_select(
    input_path: Path,
    output_prefix: Path,
    strategy: str = "random",
    target_count: Optional[int] = None,
    window_size: int = 10000,
    weights: Optional[List[float]] = None,
    targets: Optional[List[float]] = None,
    priority_bed: Optional[Path] = None,
    chromosomes: Optional[List[str]] = None,
    seed: int = 42,
    keep_biophysical: bool = False,
    merged_snps: Optional[List[SNP]] = None,
    merge_details: Optional[Dict[str, Any]] = None,
    threads: int = 1,
    verbose: bool = False,
    **kwargs,  # Accept extra kwargs for CLI compatibility
) -> Result[Dict[str, Any], str]:
    """
    Run SNP selection pipeline.
    
    Main entry point for the select command. Applies selection strategy
    to reduce SNP count to target while maintaining desired distribution.
    
    Args:
        input_path: Input SNP TSV file (used for reference format and biophysical data)
        output_prefix: Output file prefix
        strategy: Selection strategy (uniform, random, weighted, priority)
        target_count: Target number of SNPs to select (None = keep all)
        window_size: Window size for selection
        weights: Biophysical weights [gc, tm, complexity, hairpin, dimer]
        targets: Biophysical target values [gc, tm, complexity, hairpin, dimer] (default: 50,70,0,0,0)
        priority_bed: BED file with priority regions
        chromosomes: Optional list of chromosomes to include
        seed: Random seed for reproducibility
        keep_biophysical: Keep biophysical columns in output (default: False)
        merged_snps: Pre-merged SNP list (if None, load from input_path)
        merge_details: Details dictionary for merge operation (mode, file counts, etc.)
        threads: Number of threads for parallel operations (default: 1)
        verbose: Enable verbose logging
        
    Returns:
        Result containing selection statistics
    """
    # Configure logging
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        logger.addHandler(handler)
    
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    
    logger.info(f"Starting SNP selection from {input_path}")
    logger.info(f"Strategy: {strategy}, Window size: {window_size}")
    logger.info(f"Threads: {threads}")
    
    # Load input SNPs (or use pre-merged SNPs)
    if merged_snps is not None:
        # Use pre-merged SNPs
        snps = merged_snps
        snp_df_result = SNPDataFrame.from_tsv(input_path)
        if snp_df_result.is_err():
            return Err(f"Failed to load input for reference: {snp_df_result.unwrap_err()}")
        snp_df_obj = snp_df_result.unwrap()
        logger.info(f"Using pre-merged SNPs: {len(snps)} total")
    else:
        # Load from input file
        snp_df_result = SNPDataFrame.from_tsv(input_path)
        if snp_df_result.is_err():
            return Err(f"Failed to load input: {snp_df_result.unwrap_err()}")
        
        snp_df_obj = snp_df_result.unwrap()
        snps = snp_df_obj.to_snps()
    initial_count = len(snps)
    logger.info(f"Loaded {initial_count} SNPs")
    
    # Keep the raw DataFrame for weighted selection
    raw_df = snp_df_obj.df.copy()
    
    # Filter by chromosomes if specified
    if chromosomes:
        chrom_result = select_chromosomes(snps, chromosomes)
        if chrom_result.is_err():
            return Err(chrom_result.unwrap_err())
        snps = chrom_result.unwrap()
        # Also filter the DataFrame for weighted selection
        raw_df = raw_df[raw_df['chr'].isin(chromosomes)]
        logger.info(f"After chromosome filter: {len(snps)} SNPs")
    
    # Determine target count
    probe_number = target_count if target_count else len(snps)
    logger.info(f"Target count: {probe_number}")
    
    # Variable to hold DataFrame with biophysical columns (for weighted selection)
    selected_with_biophysical: Optional[pd.DataFrame] = None
    
    # Check if we already have fewer SNPs than target
    if len(snps) <= probe_number:
        logger.warning(f"Input count ({len(snps)}) <= target ({probe_number}), keeping all")
        selected = snps
    else:
        # Apply selection strategy
        strategy_lower = strategy.lower()
        
        if strategy_lower == "uniform":
            result = select_uniform(snps, probe_number, window_size, seed)
        elif strategy_lower == "random":
            result = select_random(snps, probe_number, seed)
        elif strategy_lower == "weighted":
            if weights is None:
                weights = [1.0, 1.0, 1.0, 1.0, 1.0]
            result = select_weighted(snps, probe_number, weights, window_size, seed, snp_df=raw_df, threads=threads, targets=targets)
            # Weighted returns (snps, df_with_biophysical)
            if result.is_ok():
                selected_snps, selected_with_biophysical = result.unwrap()
                selected = selected_snps
                # Calculate coverage statistics
                coverage_stats = calculate_coverage_stats(selected, window_size)
                
                # Save output
                output_path = Path(str(output_prefix) + ".selected.tsv")
                output_path.parent.mkdir(parents=True, exist_ok=True)
                
                if keep_biophysical:
                    # Use the DataFrame directly to preserve biophysical columns
                    selected_with_biophysical.to_csv(output_path, sep='\t', index=False)
                else:
                    # Remove biophysical columns from output
                    output_df = selected_with_biophysical.drop(
                        columns=BIOPHYSICAL_COLUMNS, errors='ignore'
                    )
                    output_df.to_csv(output_path, sep='\t', index=False)
                
                stats = {
                    "initial_count": initial_count,
                    "selected": len(selected),
                    "final_count": len(selected),
                    "target_count": probe_number,
                    "strategy": strategy,
                    "window_size": window_size,
                    "windows": coverage_stats["windows_covered"],
                    "coverage": coverage_stats,
                    "output_file": str(output_path),
                }
                
                # Add merge details if available
                if merge_details:
                    stats["merge_details"] = merge_details
                
                _log_biophysical_comparison(raw_df, selected_with_biophysical)
                summary_path = _write_biophysical_summary(raw_df, selected_with_biophysical, output_prefix, stats)
                stats["summary_file"] = str(summary_path)
                logger.info(f"Selection complete: {len(selected)} SNPs selected")
                logger.info(f"Windows covered: {coverage_stats['windows_covered']}")
                logger.info(f"Output saved to {output_path}")
                
                return Ok(stats)
            else:
                return Err(result.unwrap_err())
        elif strategy_lower == "priority":
            if priority_bed is None:
                return Err("Priority strategy requires --priority_bed")
            
            bed_result = parse_bed_file(priority_bed)
            if bed_result.is_err():
                return Err(bed_result.unwrap_err())
            
            priority_regions = bed_result.unwrap()
            
            # Priority strategy can optionally use weights for biophysical scoring
            result = select_by_priority(
                snps, probe_number, priority_regions, window_size, seed,
                weights=weights, snp_df=raw_df, targets=targets
            )
            
            # Check if result is tuple (weighted was used) or just list
            if result.is_ok():
                unwrapped = result.unwrap()
                if isinstance(unwrapped, tuple):
                    # Weighted selection was used within priority
                    selected_snps, selected_with_biophysical = unwrapped
                    selected = selected_snps
                    
                    # Calculate coverage statistics
                    coverage_stats = calculate_coverage_stats(selected, window_size)
                    
                    # Save output
                    output_path = Path(str(output_prefix) + ".selected.tsv")
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    
                    if keep_biophysical:
                        selected_with_biophysical.to_csv(output_path, sep='\t', index=False)
                    else:
                        output_df = selected_with_biophysical.drop(
                            columns=BIOPHYSICAL_COLUMNS, errors='ignore'
                        )
                        output_df.to_csv(output_path, sep='\t', index=False)
                    
                    stats = {
                        "initial_count": initial_count,
                        "selected": len(selected),
                        "final_count": len(selected),
                        "target_count": probe_number,
                        "strategy": f"{strategy}+weighted",
                        "window_size": window_size,
                        "windows": coverage_stats["windows_covered"],
                        "coverage": coverage_stats,
                        "output_file": str(output_path),
                    }
                    
                    # Add merge details if available
                    if merge_details:
                        stats["merge_details"] = merge_details
                    
                    _log_biophysical_comparison(raw_df, selected_with_biophysical)
                    summary_path = _write_biophysical_summary(raw_df, selected_with_biophysical, output_prefix, stats)
                    stats["summary_file"] = str(summary_path)
                    logger.info(f"Selection complete: {len(selected)} SNPs selected")
                    logger.info(f"Windows covered: {coverage_stats['windows_covered']}")
                    logger.info(f"Output saved to {output_path}")
                    
                    return Ok(stats)
                else:
                    selected = unwrapped
            else:
                return Err(result.unwrap_err())
        else:
            return Err(f"Unknown selection strategy: {strategy}. Use: uniform, random, weighted, priority")
        
        # For non-weighted strategies (uniform, random, priority without weights)
        # We need to unwrap the result and set selected
        if strategy_lower in ["uniform", "random"]:
            if result.is_err():
                return Err(result.unwrap_err())
            selected = result.unwrap()
    
    # Calculate coverage statistics
    coverage_stats = calculate_coverage_stats(selected, window_size)
    
    # Save output
    output_path = Path(str(output_prefix) + ".selected.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    selected_df = SNPDataFrame.from_snps(selected)
    save_result = selected_df.to_tsv(output_path)
    if save_result.is_err():
        return Err(f"Failed to save output: {save_result.unwrap_err()}")
    
    stats = {
        "initial_count": initial_count,
        "selected": len(selected),
        "final_count": len(selected),
        "target_count": probe_number,
        "strategy": strategy,
        "window_size": window_size,
        "windows": coverage_stats["windows_covered"],
        "coverage": coverage_stats,
        "output_file": str(output_path),
    }
    
    # Add merge details if available
    if merge_details:
        stats["merge_details"] = merge_details
    
    # Log biophysical comparison for non-weighted strategies (if biophysical cols available)
    sel_bio_df = pd.DataFrame()
    if any(c in raw_df.columns for c in BIOPHYSICAL_COLUMNS):
        sel_keys = set(zip(
            (str(s.chrom) for s in selected),
            (int(s.pos) for s in selected)
        ))
        sel_bio_df = raw_df[
            pd.Series(
                list(zip(raw_df['chr'].astype(str), raw_df['pos'].astype(int))),
                index=raw_df.index
            ).isin(sel_keys)
        ]
        _log_biophysical_comparison(raw_df, sel_bio_df)

    summary_path = _write_biophysical_summary(raw_df, sel_bio_df, output_prefix, stats)
    stats["summary_file"] = str(summary_path)
    logger.info(f"Selection complete: {len(selected)} SNPs selected")
    logger.info(f"Windows covered: {coverage_stats['windows_covered']}")
    logger.info(f"Output saved to {output_path}")
    
    return Ok(stats)


def calculate_coverage_stats(
    snps: List[SNP],
    window_size: int,
) -> Dict[str, Any]:
    """
    Calculate coverage statistics for selected SNPs.
    
    Args:
        snps: Selected SNP list
        window_size: Window size for calculation
        
    Returns:
        Dictionary of coverage statistics
    """
    if not snps:
        return {"windows_covered": 0, "chromosomes": 0}
    
    # Count unique windows
    windows = set()
    chroms = set()
    
    for snp in snps:
        chroms.add(snp.chrom)
        window_idx = calculate_window_index(snp.pos, window_size)
        windows.add((snp.chrom, window_idx))
    
    # Calculate per-chromosome stats
    chrom_counts = defaultdict(int)
    for snp in snps:
        chrom_counts[snp.chrom] += 1
    
    return {
        "windows_covered": len(windows),
        "chromosomes": len(chroms),
        "per_chromosome": dict(chrom_counts),
        "mean_per_chromosome": np.mean(list(chrom_counts.values())),
        "std_per_chromosome": np.std(list(chrom_counts.values())),
    }
