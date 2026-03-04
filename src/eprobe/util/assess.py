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
from eprobe.biophysics.fast_biophysics import (
    calculate_gc_fast,
    calculate_tm_fast,
    calculate_dust_fast,
    calculate_hairpin_fast,
    DimerCalculatorFast,
    calculate_percentile_threshold,
)

logger = logging.getLogger(__name__)

AVAILABLE_TAGS = ["gc", "tm", "complexity", "hairpin", "dimer"]


@dataclass
class FilterThresholds:
    """Biophysical filtering thresholds."""
    gc_min: float = 35.0
    gc_max: float = 65.0
    tm_min: float = 55.0
    tm_max: float = 75.0
    complexity_max: float = 2.0
    hairpin_percentile: float = 90.0
    dimer_percentile: float = 90.0


def assess_fasta(
    sequences: Dict[str, str],
    tags: List[str],
) -> List[Dict[str, Any]]:
    """
    Calculate biophysical tags for all sequences.
    
    Args:
        sequences: {probe_id: sequence}
        tags: List of tags to compute
        
    Returns:
        List of dicts with probe_id and tag values
    """
    results = []
    seq_ids = list(sequences.keys())
    seq_list = list(sequences.values())
    
    # Pre-compute dimer scores if needed (requires pool-level index)
    dimer_scores: Dict[str, float] = {}
    if "dimer" in tags and len(seq_list) > 1:
        calc = DimerCalculatorFast(k=11, include_revcomp=True)
        calc.build_index(seq_list)
        scores = calc.calculate_all_scores()
        dimer_scores = dict(zip(seq_ids, scores))
    
    for sid, seq in sequences.items():
        row: Dict[str, Any] = {"probe_id": sid, "length": len(seq)}
        seq_upper = seq.upper()
        
        if "gc" in tags:
            try:
                row["gc"] = round(calculate_gc_fast(seq_upper), 2)
            except Exception:
                row["gc"] = None
        
        if "tm" in tags:
            try:
                row["tm"] = round(calculate_tm_fast(seq_upper), 2)
            except Exception:
                row["tm"] = None
        
        if "complexity" in tags:
            try:
                row["complexity"] = round(calculate_dust_fast(seq_upper), 4)
            except Exception:
                row["complexity"] = None
        
        if "hairpin" in tags:
            try:
                row["hairpin"] = round(calculate_hairpin_fast(seq_upper), 2)
            except Exception:
                row["hairpin"] = None
        
        if "dimer" in tags:
            row["dimer"] = round(dimer_scores.get(sid, 0.0), 4)
        
        results.append(row)
    
    return results


def filter_fasta(
    sequences: Dict[str, str],
    thresholds: FilterThresholds,
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
    
    # Stage 1: GC, Tm, DUST
    stage1_passed: Dict[str, str] = OrderedDict()
    for sid, seq in sequences.items():
        seq_upper = seq.upper()
        try:
            gc = calculate_gc_fast(seq_upper)
            if not (thresholds.gc_min <= gc <= thresholds.gc_max):
                stats["gc_failed"] += 1
                failed[sid] = seq
                continue
            
            tm = calculate_tm_fast(seq_upper)
            if not (thresholds.tm_min <= tm <= thresholds.tm_max):
                stats["tm_failed"] += 1
                failed[sid] = seq
                continue
            
            dust = calculate_dust_fast(seq_upper)
            if dust > thresholds.complexity_max:
                stats["complexity_failed"] += 1
                failed[sid] = seq
                continue
            
            stage1_passed[sid] = seq
        except Exception:
            failed[sid] = seq
            continue
    
    # Stage 2: Hairpin (percentile)
    if stage1_passed:
        hairpin_scores = {
            sid: calculate_hairpin_fast(seq.upper())
            for sid, seq in stage1_passed.items()
        }
        
        threshold = calculate_percentile_threshold(
            list(hairpin_scores.values()),
            thresholds.hairpin_percentile,
            higher_is_worse=True,
        )
        
        stage2_passed: Dict[str, str] = OrderedDict()
        for sid, seq in stage1_passed.items():
            if hairpin_scores[sid] <= threshold:
                stage2_passed[sid] = seq
            else:
                stats["hairpin_failed"] += 1
                failed[sid] = seq
    else:
        stage2_passed = OrderedDict()
    
    # Stage 3: Dimer (percentile)
    if len(stage2_passed) > 1:
        calc = DimerCalculatorFast(k=11, include_revcomp=True)
        seq_list = list(stage2_passed.values())
        calc.build_index(seq_list)
        dimer_scores_list = calc.calculate_all_scores()
        
        ids = list(stage2_passed.keys())
        dimer_scores = dict(zip(ids, dimer_scores_list))
        
        threshold = calculate_percentile_threshold(
            dimer_scores_list,
            thresholds.dimer_percentile,
            higher_is_worse=True,
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
    hairpin_percentile: float = 90.0,
    dimer_percentile: float = 90.0,
    generate_plots: bool = True,
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
        hairpin_percentile: Hairpin threshold percentile
        dimer_percentile: Dimer threshold percentile
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
        results = assess_fasta(sequences, tags)
        
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
            hairpin_percentile=hairpin_percentile,
            dimer_percentile=dimer_percentile,
        )
        
        logger.info(
            f"Filtering: GC={gc_min}-{gc_max}%, Tm={tm_min}-{tm_max}°C, "
            f"DUST≤{complexity_max}"
        )
        
        passed, failed_seqs, filter_stats = filter_fasta(sequences, thresholds)
        
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
