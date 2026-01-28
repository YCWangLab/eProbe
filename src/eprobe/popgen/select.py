"""
SNP selection module.

Implements strategies for selecting SNPs from filtered candidates:
  - Uniform distribution: equal spacing across genome
  - Random sampling: random selection with optional stratification
  - Priority-based: selection based on quality scores

This module corresponds to the original SNP_subsampler.py functionality.
"""

import logging
import random
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass
from collections import defaultdict
from enum import Enum

import numpy as np

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNP, SNPDataFrame

logger = logging.getLogger(__name__)


class SelectionStrategy(Enum):
    """SNP selection strategies."""
    UNIFORM = "uniform"
    RANDOM = "random"
    PRIORITY = "priority"


@dataclass
class WindowConfig:
    """Configuration for window-based selection."""
    size: int = 100000  # Window size in bp
    per_window: int = 1  # SNPs to select per window


def calculate_window_index(pos: int, window_size: int) -> int:
    """Calculate window index for a genomic position."""
    return pos // window_size


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
    priority_column: str = "score",
    ascending: bool = False,
) -> Result[List[SNP], str]:
    """
    Select SNPs based on priority score.
    
    Selects top-scoring SNPs based on a priority metric.
    Higher scores are selected by default.
    
    Args:
        snps: Input SNP list with scores
        target_count: Target number of SNPs to select
        priority_column: Column name containing scores
        ascending: If True, select lowest scores instead
        
    Returns:
        Result containing selected SNP list
    """
    if not snps:
        return Ok([])
    
    logger.info(f"Priority selection: {len(snps)} → {target_count} SNPs")
    logger.info(f"Priority column: {priority_column}, ascending: {ascending}")
    
    # Check if SNPs have the priority attribute
    # Note: This would require SNP metadata support
    # For now, use random as fallback
    logger.warning(f"Priority selection not yet implemented, falling back to random")
    return select_random(snps, target_count)


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
    strategy: str = "uniform",
    probe_number: int = 10000,
    window_size: int = 100000,
    chromosomes: Optional[List[str]] = None,
    seed: int = 42,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Run SNP selection pipeline.
    
    Main entry point for the select command. Applies selection strategy
    to reduce SNP count to target while maintaining desired distribution.
    
    Args:
        input_path: Input SNP TSV file
        output_prefix: Output file prefix
        strategy: Selection strategy (uniform, random, priority)
        probe_number: Target number of SNPs to select
        window_size: Window size for uniform strategy
        chromosomes: Optional list of chromosomes to include
        seed: Random seed for reproducibility
        verbose: Enable verbose logging
        
    Returns:
        Result containing selection statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Starting SNP selection from {input_path}")
    logger.info(f"Strategy: {strategy}, Target: {probe_number}")
    
    # Load input SNPs
    snp_df = SNPDataFrame.from_tsv(input_path)
    if snp_df.is_err():
        return Err(f"Failed to load input: {snp_df.unwrap_err()}")
    
    snps = snp_df.unwrap().to_snps()
    initial_count = len(snps)
    logger.info(f"Loaded {initial_count} SNPs")
    
    # Filter by chromosomes if specified
    if chromosomes:
        chrom_result = select_chromosomes(snps, chromosomes)
        if chrom_result.is_err():
            return Err(chrom_result.unwrap_err())
        snps = chrom_result.unwrap()
        logger.info(f"After chromosome filter: {len(snps)} SNPs")
    
    # Check if we already have fewer SNPs than target
    if len(snps) <= probe_number:
        logger.warning(f"Input count ({len(snps)}) <= target ({probe_number})")
        selected = snps
    else:
        # Apply selection strategy
        strategy_lower = strategy.lower()
        
        if strategy_lower == "uniform":
            result = select_uniform(snps, probe_number, window_size, seed)
        elif strategy_lower == "random":
            result = select_random(snps, probe_number, seed)
        elif strategy_lower == "priority":
            result = select_by_priority(snps, probe_number)
        else:
            return Err(f"Unknown selection strategy: {strategy}")
        
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
        "final_count": len(selected),
        "target_count": probe_number,
        "strategy": strategy,
        "window_size": window_size,
        "coverage": coverage_stats,
        "output_file": str(output_path),
    }
    
    logger.info(f"Selection complete: {len(selected)} SNPs selected")
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
