"""
Probe quality assessment module.

Calculates and reports quality metrics for probe sets:
  - Biophysical properties (GC, Tm, complexity, hairpin, dimer)
  - Genomic coverage statistics
  - Distribution plots

This module corresponds to the original SNP_assessor.py functionality.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from dataclasses import dataclass

import pandas as pd
import numpy as np

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNPDataFrame, ProbeSet
from eprobe.core.fasta import read_fasta
from eprobe.biophysics import (
    calculate_gc,
    calculate_tm,
    calculate_complexity,
    calculate_hairpin_score,
)
from eprobe.biophysics.dimer import DimerCalculator
from eprobe.biophysics.entropy import calculate_entropy

logger = logging.getLogger(__name__)


# Available assessment tags
AVAILABLE_TAGS = {
    "gc": "GC content (%)",
    "tm": "Melting temperature (°C)",
    "complexity": "DUST complexity score",
    "hairpin": "Self-complementarity score",
    "dimer": "Inter-probe complementarity score",
    "entropy": "Shannon entropy",
}


@dataclass
class AssessmentResult:
    """Container for assessment results."""
    probe_id: str
    sequence: str
    gc: Optional[float] = None
    tm: Optional[float] = None
    complexity: Optional[float] = None
    hairpin: Optional[float] = None
    dimer: Optional[float] = None
    entropy: Optional[float] = None


def assess_single_probe(
    probe_id: str,
    sequence: str,
    tags: List[str],
) -> AssessmentResult:
    """
    Calculate assessment metrics for a single probe.
    
    Args:
        probe_id: Probe identifier
        sequence: Probe sequence
        tags: List of metrics to calculate
        
    Returns:
        AssessmentResult with calculated values
    """
    result = AssessmentResult(probe_id=probe_id, sequence=sequence)
    
    if "gc" in tags:
        result.gc = calculate_gc(sequence)
    
    if "tm" in tags:
        result.tm = calculate_tm(sequence)
    
    if "complexity" in tags:
        result.complexity = calculate_complexity(sequence)
    
    if "hairpin" in tags:
        result.hairpin = calculate_hairpin_score(sequence)
    
    if "entropy" in tags:
        result.entropy = calculate_entropy(sequence)
    
    return result


def assess_probes(
    sequences: Dict[str, str],
    tags: List[str],
) -> List[AssessmentResult]:
    """
    Assess all probes for specified metrics.
    
    Args:
        sequences: Dictionary of probe_id -> sequence
        tags: List of metrics to calculate
        
    Returns:
        List of AssessmentResult objects
    """
    results = []
    
    for probe_id, sequence in sequences.items():
        result = assess_single_probe(probe_id, sequence, tags)
        results.append(result)
    
    return results


def calculate_dimer_scores(
    sequences: Dict[str, str],
    sample_size: int = 1000,
) -> Dict[str, float]:
    """
    Calculate dimer scores by sampling pairs.
    
    For large probe sets, calculating all pairwise dimer scores
    is computationally expensive. This function samples pairs
    to estimate dimer propensity.
    
    Args:
        sequences: Dictionary of probe_id -> sequence
        sample_size: Number of pairs to sample
        
    Returns:
        Dictionary of probe_id -> max dimer score
    """
    calculator = DimerCalculator()
    probe_ids = list(sequences.keys())
    n_probes = len(probe_ids)
    
    if n_probes <= 1:
        return {pid: 0.0 for pid in probe_ids}
    
    # Initialize scores
    max_scores: Dict[str, float] = {pid: 0.0 for pid in probe_ids}
    
    # Sample pairs
    n_pairs = min(sample_size, n_probes * (n_probes - 1) // 2)
    
    for _ in range(n_pairs):
        # Random pair selection
        i, j = np.random.choice(n_probes, 2, replace=False)
        
        seq1 = sequences[probe_ids[i]]
        seq2 = sequences[probe_ids[j]]
        
        score = calculator.calculate(seq1, seq2)
        
        max_scores[probe_ids[i]] = max(max_scores[probe_ids[i]], score)
        max_scores[probe_ids[j]] = max(max_scores[probe_ids[j]], score)
    
    return max_scores


def calculate_summary_statistics(
    results: List[AssessmentResult],
    tags: List[str],
) -> Dict[str, Dict[str, float]]:
    """
    Calculate summary statistics for each metric.
    
    Args:
        results: List of assessment results
        tags: List of metrics to summarize
        
    Returns:
        Dictionary of metric -> {mean, std, min, max, median}
    """
    summary = {}
    
    for tag in tags:
        values = [getattr(r, tag) for r in results if getattr(r, tag) is not None]
        
        if not values:
            continue
        
        summary[tag] = {
            "mean": float(np.mean(values)),
            "std": float(np.std(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "median": float(np.median(values)),
            "count": len(values),
        }
    
    return summary


def generate_assessment_plots(
    results: List[AssessmentResult],
    tags: List[str],
    output_prefix: Path,
) -> Result[List[Path], str]:
    """
    Generate distribution plots for assessment metrics.
    
    Args:
        results: List of assessment results
        tags: List of metrics to plot
        output_prefix: Output file prefix
        
    Returns:
        Result containing list of generated plot paths
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        return Err("Plotting requires matplotlib and seaborn. Install with: pip install eprobe[plot]")
    
    plot_paths = []
    
    for tag in tags:
        values = [getattr(r, tag) for r in results if getattr(r, tag) is not None]
        
        if not values:
            continue
        
        fig, ax = plt.subplots(figsize=(8, 5))
        
        # Histogram with KDE
        sns.histplot(values, kde=True, ax=ax)
        
        # Add statistics annotation
        mean_val = np.mean(values)
        std_val = np.std(values)
        ax.axvline(mean_val, color='red', linestyle='--', label=f'Mean: {mean_val:.2f}')
        
        ax.set_xlabel(AVAILABLE_TAGS.get(tag, tag))
        ax.set_ylabel("Count")
        ax.set_title(f"Distribution of {tag.upper()}")
        ax.legend()
        
        # Save plot
        plot_path = Path(str(output_prefix) + f".{tag}_dist.png")
        fig.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        
        plot_paths.append(plot_path)
        logger.info(f"Saved plot: {plot_path}")
    
    return Ok(plot_paths)


def run_assess(
    input_path: Path,
    output_prefix: Path,
    tags: Optional[List[str]] = None,
    generate_plots: bool = True,
    sample_dimer: int = 1000,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Assess quality of probe set.
    
    Main entry point for the assess command. Calculates biophysical
    metrics and generates summary statistics and plots.
    
    Args:
        input_path: Input file (FASTA or TSV)
        output_prefix: Output file prefix
        tags: List of metrics to calculate (default: all)
        generate_plots: Generate distribution plots
        sample_dimer: Number of pairs to sample for dimer calculation
        threads: Number of threads
        verbose: Enable verbose logging
        
    Returns:
        Result containing assessment statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    # Default to all tags except dimer (expensive)
    if tags is None:
        tags = ["gc", "tm", "complexity", "hairpin"]
    
    # Normalize tag names
    tags = [t.lower() for t in tags]
    
    # Validate tags
    invalid_tags = [t for t in tags if t not in AVAILABLE_TAGS]
    if invalid_tags:
        return Err(f"Invalid tags: {invalid_tags}. Available: {list(AVAILABLE_TAGS.keys())}")
    
    logger.info(f"Starting probe assessment from {input_path}")
    logger.info(f"Tags: {tags}")
    
    # Load input
    suffix = input_path.suffix.lower()
    
    if suffix in ['.fa', '.fasta', '.fna']:
        # Load from FASTA
        fasta_result = read_fasta(input_path)
        if fasta_result.is_err():
            return Err(f"Failed to read FASTA: {fasta_result.unwrap_err()}")
        sequences = fasta_result.unwrap()
    elif suffix == '.tsv':
        # Load from TSV (assume SNP format)
        snp_df = SNPDataFrame.from_tsv(input_path)
        if snp_df.is_err():
            return Err(f"Failed to read TSV: {snp_df.unwrap_err()}")
        
        snps = snp_df.unwrap().to_snps()
        sequences = {
            snp.snp_id: snp.left_flank + snp.ref + snp.right_flank
            for snp in snps
        }
    else:
        return Err(f"Unsupported input format: {suffix}")
    
    logger.info(f"Loaded {len(sequences)} sequences")
    
    # Calculate individual probe metrics
    logger.info("Calculating probe metrics...")
    tags_without_dimer = [t for t in tags if t != "dimer"]
    results = assess_probes(sequences, tags_without_dimer)
    
    # Calculate dimer scores separately (sampling-based)
    if "dimer" in tags:
        logger.info(f"Calculating dimer scores (sampling {sample_dimer} pairs)...")
        dimer_scores = calculate_dimer_scores(sequences, sample_dimer)
        
        for result in results:
            result.dimer = dimer_scores.get(result.probe_id, 0.0)
    
    # Calculate summary statistics
    summary = calculate_summary_statistics(results, tags)
    
    # Save detailed results to TSV
    output_path = Path(str(output_prefix) + ".stats.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Convert results to DataFrame
    data = []
    for r in results:
        row = {"probe_id": r.probe_id}
        for tag in tags:
            row[tag] = getattr(r, tag)
        data.append(row)
    
    df = pd.DataFrame(data)
    df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Saved detailed stats to {output_path}")
    
    # Save summary
    summary_path = Path(str(output_prefix) + ".summary.txt")
    with open(summary_path, 'w') as f:
        f.write("eProbe Assessment Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Input: {input_path}\n")
        f.write(f"Total probes: {len(sequences)}\n\n")
        
        for tag, stats in summary.items():
            f.write(f"{AVAILABLE_TAGS.get(tag, tag)}:\n")
            f.write(f"  Mean:   {stats['mean']:.3f}\n")
            f.write(f"  Std:    {stats['std']:.3f}\n")
            f.write(f"  Min:    {stats['min']:.3f}\n")
            f.write(f"  Max:    {stats['max']:.3f}\n")
            f.write(f"  Median: {stats['median']:.3f}\n\n")
    
    logger.info(f"Saved summary to {summary_path}")
    
    # Generate plots
    plot_paths = []
    if generate_plots:
        plot_result = generate_assessment_plots(results, tags, output_prefix)
        if plot_result.is_ok():
            plot_paths = plot_result.unwrap()
    
    stats = {
        "probe_count": len(sequences),
        "tags": tags,
        "summary": summary,
        "stats_file": str(output_path),
        "summary_file": str(summary_path),
        "plot_files": [str(p) for p in plot_paths],
    }
    
    logger.info(f"Assessment complete for {len(sequences)} probes")
    
    return Ok(stats)
