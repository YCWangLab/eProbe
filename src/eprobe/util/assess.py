"""
Biophysical assessment and filtering for probe FASTA files.

Modes:
  tags:   Calculate GC, Tm, complexity, hairpin, dimer → TSV + plots
  filter: Apply biophysical thresholds → filtered FASTA + rejected FASTA

Reuses the same fast biophysics calculators as popgen filter/assess.
"""

import logging
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass

import numpy as np

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta
from eprobe.biophysics.biophysics import (
    calculate_gc_fast,
    calculate_tm_fast,
    calculate_dust_fast,
    calculate_hairpin_fast,
    DimerCalculatorFast,
    calculate_percentile_threshold,
    compute_biophysical_parallel,
)

logger = logging.getLogger(__name__)

AVAILABLE_TAGS = ["gc", "tm", "complexity", "hairpin", "dimer"]


def _resolve_threshold_assess(value: float, scores: list) -> float:
    """
    Interpret dual-mode threshold: >=1 absolute, 0<x<1 percentile, 0 skip.
    
    Returns float('inf') if value <= 0 (effectively disabling the filter).
    """
    if value <= 0:
        return float('inf')
    elif value < 1.0:
        return calculate_percentile_threshold(
            scores, value * 100, higher_is_worse=True
        )
    else:
        return value


@dataclass
class FilterThresholds:
    """Biophysical filtering thresholds.
    
    Dual-mode for hairpin/dimer:
    - >=1: absolute threshold
    - 0<x<1: percentile (e.g. 0.95 keeps 95%)
    - 0: skip filter
    """
    gc_min: float = 35.0
    gc_max: float = 65.0
    tm_min: float = 55.0
    tm_max: float = 75.0
    complexity_max: float = 2.0
    hairpin: float = 15.0
    dimer: float = 0.95


def assess_fasta(
    sequences: Dict[str, str],
    tags: List[str],
    threads: int = 1,
) -> List[Dict[str, Any]]:
    """
    Calculate biophysical tags for all sequences.

    Args:
        sequences: {probe_id: sequence}
        tags: List of tags to compute
        threads: Number of worker processes for parallel computation

    Returns:
        List of dicts with probe_id and tag values
    """
    seq_ids = list(sequences.keys())
    seq_list = [seq.upper() for seq in sequences.values()]

    # Parallel computation for gc, tm, complexity, hairpin
    need_parallel = any(t in tags for t in ("gc", "tm", "complexity", "hairpin"))
    gc_vals, tm_vals, cplx_vals, hp_vals = [], [], [], []
    if need_parallel:
        gc_vals, tm_vals, cplx_vals, hp_vals = compute_biophysical_parallel(
            seq_list, threads=threads,
        )

    # Pre-compute dimer scores if needed (requires pool-level index)
    dimer_scores: Dict[str, float] = {}
    if "dimer" in tags and len(seq_list) > 1:
        calc = DimerCalculatorFast(k=11)
        calc.build_index(seq_list)
        scores = calc.calculate_all_scores()
        dimer_scores = dict(zip(seq_ids, scores))

    results = []
    for i, (sid, seq) in enumerate(zip(seq_ids, seq_list)):
        row: Dict[str, Any] = {"probe_id": sid, "length": len(seq)}

        if "gc" in tags:
            row["gc"] = round(gc_vals[i], 2) if gc_vals else None
        if "tm" in tags:
            row["tm"] = round(tm_vals[i], 2) if tm_vals else None
        if "complexity" in tags:
            row["complexity"] = round(cplx_vals[i], 4) if cplx_vals else None
        if "hairpin" in tags:
            row["hairpin"] = round(hp_vals[i], 2) if hp_vals else None
        if "dimer" in tags:
            row["dimer"] = round(dimer_scores.get(sid, 0.0), 4)

        results.append(row)

    return results


def filter_fasta(
    sequences: Dict[str, str],
    thresholds: FilterThresholds,
    threads: int = 1,
) -> tuple:
    """
    Filter probes by biophysical thresholds.

    3-stage filtering:
      Stage 1: GC, Tm, DUST (absolute thresholds)
      Stage 2: Hairpin (percentile-based)
      Stage 3: Dimer (percentile-based)

    Returns:
        (passed_dict, failed_dict, filter_stats)
    """
    from collections import OrderedDict

    passed: Dict[str, str] = OrderedDict()
    failed: Dict[str, str] = OrderedDict()
    stats = {
        "gc_failed": 0, "tm_failed": 0, "complexity_failed": 0,
        "hairpin_failed": 0, "dimer_failed": 0,
    }

    # Parallel computation of all 4 metrics for all sequences
    seq_ids = list(sequences.keys())
    seq_list = [seq.upper() for seq in sequences.values()]
    gc_all, tm_all, dust_all, hp_all = compute_biophysical_parallel(
        seq_list, threads=threads,
    )

    # Stage 1: GC, Tm, DUST
    stage1_passed: Dict[str, str] = OrderedDict()
    stage1_hairpin: Dict[str, float] = {}  # reuse hairpin from parallel computation
    for sid, seq, gc, tm, dust, hp in zip(seq_ids, list(sequences.values()), gc_all, tm_all, dust_all, hp_all):
        if not (thresholds.gc_min <= gc <= thresholds.gc_max):
            stats["gc_failed"] += 1
            failed[sid] = seq
            continue

        if not (thresholds.tm_min <= tm <= thresholds.tm_max):
            stats["tm_failed"] += 1
            failed[sid] = seq
            continue

        if dust > thresholds.complexity_max:
            stats["complexity_failed"] += 1
            failed[sid] = seq
            continue

        stage1_passed[sid] = seq
        stage1_hairpin[sid] = hp

    # Stage 2: Hairpin (dual-mode) — reuse pre-computed values
    if stage1_passed and thresholds.hairpin > 0:
        hairpin_scores = stage1_hairpin
        
        threshold = _resolve_threshold_assess(
            thresholds.hairpin, list(hairpin_scores.values())
        )
        
        stage2_passed: Dict[str, str] = OrderedDict()
        for sid, seq in stage1_passed.items():
            if hairpin_scores[sid] <= threshold:
                stage2_passed[sid] = seq
            else:
                stats["hairpin_failed"] += 1
                failed[sid] = seq
    else:
        stage2_passed = stage1_passed
    
    # Stage 3: Dimer (dual-mode)
    if len(stage2_passed) > 1 and thresholds.dimer > 0:
        calc = DimerCalculatorFast(k=11)
        seq_list = list(stage2_passed.values())
        calc.build_index(seq_list)
        dimer_scores_list = calc.calculate_all_scores()
        
        ids = list(stage2_passed.keys())
        dimer_scores = dict(zip(ids, dimer_scores_list))
        
        threshold = _resolve_threshold_assess(
            thresholds.dimer, dimer_scores_list
        )
        
        for sid, seq in stage2_passed.items():
            if dimer_scores[sid] <= threshold:
                passed[sid] = seq
            else:
                stats["dimer_failed"] += 1
                failed[sid] = seq
    else:
        passed = stage2_passed
    
    stats["input_count"] = len(sequences)
    stats["passed_count"] = len(passed)
    stats["failed_count"] = len(failed)
    
    return passed, failed, stats


def _generate_tag_plots(results, tags, output_prefix):
    """Generate distribution plots for biophysical tags."""
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib/seaborn not available, skipping plots")
        return []
    
    colors = ['#83639F', '#EA7827', '#C22f2F', '#449945', '#1F70A9']
    plot_paths = []
    
    for idx, tag in enumerate(tags):
        values = [r[tag] for r in results if r.get(tag) is not None]
        if not values:
            continue
        
        color = colors[idx % len(colors)]
        sns.set(style='darkgrid')
        fig, ax = plt.subplots(figsize=(8, 6))
        
        sns.histplot(values, alpha=0.7, stat='percent', edgecolor=color,
                     color=color, element='step', line_kws={'linewidth': 2}, ax=ax)
        
        ax.set_title(f'{tag.upper()} Distribution', fontsize=18)
        ax.set_xlabel(tag, fontsize=14)
        ax.set_ylabel('Percent (%)', fontsize=14)
        
        plt.tight_layout()
        plot_path = Path(str(output_prefix) + f".{tag}_dist.jpg")
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)
        
        plot_paths.append(str(plot_path))
    
    return plot_paths


def run_assess(
    input_path: Path,
    output_prefix: Path,
    mode: str = "tags",
    tags: Optional[List[str]] = None,
    gc_min: float = 35.0,
    gc_max: float = 65.0,
    tm_min: float = 55.0,
    tm_max: float = 75.0,
    complexity_max: float = 2.0,
    hairpin: float = 15.0,
    dimer: float = 0.95,
    generate_plots: bool = True,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Assess or filter probe FASTA by biophysical properties.
    
    Modes:
      tags: Calculate biophysical metrics → TSV + plots
      filter: Apply thresholds → filtered FASTA
      
    Args:
        input_path: Input FASTA file
        output_prefix: Output prefix
        mode: "tags" or "filter"
        tags: List of tags (for tags mode)
        gc_min/gc_max: GC range (for filter mode)
        tm_min/tm_max: Tm range (for filter mode)
        complexity_max: DUST max (for filter mode)
        hairpin: Hairpin threshold (>=1 absolute, <1 percentile)
        dimer: Dimer threshold (>=1 absolute, <1 percentile)
        generate_plots: Generate distribution plots
        verbose: Verbose logging
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    fasta_result = read_fasta(input_path)
    if fasta_result.is_err():
        return Err(f"Failed to read input: {fasta_result.unwrap_err()}")
    
    sequences = fasta_result.unwrap()
    n_input = len(sequences)
    logger.info(f"Loaded {n_input} probes")
    
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    if mode == "tags":
        if tags is None:
            tags = ["gc", "tm", "complexity", "hairpin"]
        
        invalid = [t for t in tags if t not in AVAILABLE_TAGS]
        if invalid:
            return Err(f"Invalid tags: {invalid}. Available: {AVAILABLE_TAGS}")
        
        logger.info(f"Computing tags: {tags}")
        results = assess_fasta(sequences, tags, threads=threads)
        
        import pandas as pd
        df = pd.DataFrame(results)
        tsv_path = Path(str(output_prefix) + ".tags_stats.tsv")
        df.to_csv(tsv_path, sep='\t', index=False)
        
        summary: Dict[str, Dict[str, float]] = {}
        for tag in tags:
            values = [r[tag] for r in results if r.get(tag) is not None]
            if values:
                summary[tag] = {
                    "mean": round(float(np.mean(values)), 3),
                    "std": round(float(np.std(values)), 3),
                    "min": round(float(np.min(values)), 3),
                    "max": round(float(np.max(values)), 3),
                    "median": round(float(np.median(values)), 3),
                }
        
        summary_path = Path(str(output_prefix) + ".tags_summary.txt")
        with open(summary_path, 'w') as f:
            f.write("Biophysical Assessment Summary\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Input: {input_path}\n")
            f.write(f"Probes: {n_input}\n\n")
            for tag, s in summary.items():
                f.write(f"{tag}:\n")
                for k, v in s.items():
                    f.write(f"  {k}: {v}\n")
                f.write("\n")
        
        plot_paths: List[str] = []
        if generate_plots:
            plot_paths = _generate_tag_plots(results, tags, output_prefix)
        
        return Ok({
            "probe_count": n_input,
            "tags": tags,
            "summary": summary,
            "stats_file": str(tsv_path),
            "summary_file": str(summary_path),
            "plot_files": plot_paths,
        })
    
    elif mode == "filter":
        thresholds = FilterThresholds(
            gc_min=gc_min, gc_max=gc_max,
            tm_min=tm_min, tm_max=tm_max,
            complexity_max=complexity_max,
            hairpin=hairpin,
            dimer=dimer,
        )
        
        logger.info(
            f"Filtering: GC={gc_min}-{gc_max}%, Tm={tm_min}-{tm_max}°C, "
            f"DUST≤{complexity_max}"
        )
        
        passed, failed_seqs, filter_stats = filter_fasta(sequences, thresholds, threads=threads)
        
        passed_path = Path(str(output_prefix) + ".filtered.fa")
        write_result = write_fasta(passed, passed_path)
        if write_result.is_err():
            return Err(f"Write failed: {write_result.unwrap_err()}")
        
        if failed_seqs:
            failed_path = Path(str(output_prefix) + ".filtered_out.fa")
            write_fasta(failed_seqs, failed_path)
        
        return Ok({
            "input_count": n_input,
            "passed_count": len(passed),
            "failed_count": len(failed_seqs),
            "filter_stats": filter_stats,
            "fasta_file": str(passed_path),
        })
    
    else:
        return Err(f"Unknown mode: {mode}. Use 'tags' or 'filter'.")
