#!/usr/bin/env python3
"""
Benchmark eProbe biophysical metrics against open-source tools.

Compares:
  - Tm:      eProbe (Biopython Tm_NN, R_DNA_NN1)  vs  primer3 (SantaLucia 1998)
  - Hairpin:  eProbe (PairwiseAligner self-RC)      vs  primer3 calc_hairpin (ΔG)
                                                    vs  ViennaRNA MFE (Zuker, C)
                                                    vs  seqfold MFE (Zuker, Python)
  - Dimer:    eProbe (k-mer frequency)              vs  primer3 calc_homodimer (ΔG)
  - GC:       eProbe  vs  simple count (should be identical, sanity check)
  - Complexity: eProbe DUST  vs  ViennaRNA MFE (secondary structure propensity)

Usage:
    python scripts/benchmark_biophysics.py -f probes.fasta -o benchmark_output [-n 2000]
    python scripts/benchmark_biophysics.py -f selected.tsv -r ref.fa -o benchmark_output

Requires: pip install primer3-py seqfold ViennaRNA
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_sequences(
    fasta_path: Optional[Path] = None,
    tsv_path: Optional[Path] = None,
    ref_path: Optional[Path] = None,
    max_n: Optional[int] = None,
) -> Dict[str, str]:
    """Load probe sequences from FASTA or TSV + reference."""
    seqs: Dict[str, str] = OrderedDict()

    if fasta_path and fasta_path.exists():
        from Bio import SeqIO
        for rec in SeqIO.parse(str(fasta_path), "fasta"):
            seqs[rec.id] = str(rec.seq).upper()
    elif tsv_path and tsv_path.exists():
        if ref_path is None or not ref_path.exists():
            raise FileNotFoundError("TSV input requires --ref")
        from eprobe.core.fasta import read_fasta
        from eprobe.popgen.assess import _extract_probe_sequences
        ref = read_fasta(str(ref_path))
        df = pd.read_csv(tsv_path, sep="\t")
        seqs = _extract_probe_sequences(df, ref)
    else:
        raise FileNotFoundError("Provide --fasta or --tsv (+ --ref)")

    if max_n and len(seqs) > max_n:
        keys = list(seqs.keys())[:max_n]
        seqs = OrderedDict((k, seqs[k]) for k in keys)
        logger.info(f"Subsampled to {max_n} probes")

    logger.info(f"Loaded {len(seqs)} probes")
    return seqs


# ---------------------------------------------------------------------------
# eProbe calculations
# ---------------------------------------------------------------------------

def eprobe_tm(seqs: Dict[str, str]) -> np.ndarray:
    from eprobe.biophysics.tm import calculate_tm
    vals = []
    for seq in seqs.values():
        r = calculate_tm(seq)
        vals.append(r.unwrap() if r.is_ok() else np.nan)
    return np.array(vals)


def eprobe_gc(seqs: Dict[str, str]) -> np.ndarray:
    from eprobe.biophysics.gc import calculate_gc
    vals = []
    for seq in seqs.values():
        r = calculate_gc(seq)
        vals.append(r.unwrap() if r.is_ok() else np.nan)
    return np.array(vals)


def eprobe_hairpin(seqs: Dict[str, str]) -> np.ndarray:
    from eprobe.biophysics.hairpin import calculate_hairpin_score
    vals = []
    for seq in seqs.values():
        r = calculate_hairpin_score(seq)
        vals.append(r.unwrap() if r.is_ok() else np.nan)
    return np.array(vals)


def eprobe_dust(seqs: Dict[str, str]) -> np.ndarray:
    from eprobe.biophysics.complexity import calculate_dust_score
    vals = []
    for seq in seqs.values():
        r = calculate_dust_score(seq)
        vals.append(r.unwrap() if r.is_ok() else np.nan)
    return np.array(vals)


def eprobe_dimer(seqs: Dict[str, str]) -> np.ndarray:
    from eprobe.biophysics.dimer import DimerCalculator
    calc = DimerCalculator(k=11, min_freq=2)
    calc.build_kmer_index(OrderedDict(seqs))
    vals = []
    for seq in seqs.values():
        r = calc.calculate_score(seq)
        vals.append(r.unwrap() if r.is_ok() else np.nan)
    return np.array(vals)


# ---------------------------------------------------------------------------
# primer3 calculations
# ---------------------------------------------------------------------------

def primer3_tm(seqs: Dict[str, str]) -> np.ndarray:
    import primer3
    vals = []
    for seq in seqs.values():
        try:
            vals.append(primer3.calc_tm(seq))
        except Exception:
            vals.append(np.nan)
    return np.array(vals)


def _primer3_sliding_min_dg(seq: str, calc_func, window: int = 60) -> float:
    """
    primer3 structural calcs have a 60bp limit.
    Slide a window across longer sequences, return the minimum (most stable) ΔG.
    """
    import primer3
    if len(seq) <= window:
        res = calc_func(seq)
        return res.dg / 1000.0

    min_dg = 0.0
    for i in range(len(seq) - window + 1):
        sub = seq[i:i + window]
        res = calc_func(sub)
        if res.structure_found and res.dg / 1000.0 < min_dg:
            min_dg = res.dg / 1000.0
    return min_dg


def primer3_hairpin_dg(seqs: Dict[str, str]) -> np.ndarray:
    """Return hairpin ΔG in kcal/mol (most stable window)."""
    import primer3
    vals = []
    for seq in seqs.values():
        try:
            vals.append(_primer3_sliding_min_dg(seq, primer3.calc_hairpin))
        except Exception:
            vals.append(np.nan)
    return np.array(vals)


def primer3_homodimer_dg(seqs: Dict[str, str]) -> np.ndarray:
    """Return homodimer ΔG in kcal/mol (most stable window)."""
    import primer3
    vals = []
    for seq in seqs.values():
        try:
            vals.append(_primer3_sliding_min_dg(seq, primer3.calc_homodimer))
        except Exception:
            vals.append(np.nan)
    return np.array(vals)


# ---------------------------------------------------------------------------
# ViennaRNA calculations (Zuker MFE, C implementation — fast, no length limit)
# ---------------------------------------------------------------------------

def vienna_mfe(seqs: Dict[str, str]) -> np.ndarray:
    import RNA
    vals = []
    for seq in seqs.values():
        try:
            _ss, mfe = RNA.fold(seq)
            vals.append(mfe)
        except Exception:
            vals.append(np.nan)
    return np.array(vals)


# ---------------------------------------------------------------------------
# seqfold calculations (Zuker MFE, pure Python — slow, no length limit)
# ---------------------------------------------------------------------------

def seqfold_mfe(seqs: Dict[str, str]) -> np.ndarray:
    from seqfold import dg
    vals = []
    for seq in seqs.values():
        try:
            vals.append(dg(seq))
        except Exception:
            vals.append(np.nan)
    return np.array(vals)


# ---------------------------------------------------------------------------
# Correlation & plotting
# ---------------------------------------------------------------------------

def pearson_r(a: np.ndarray, b: np.ndarray) -> Tuple[float, float]:
    """Return (r, p-value) for non-nan pairs."""
    from scipy import stats as sp_stats
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return (np.nan, np.nan)
    return sp_stats.pearsonr(a[mask], b[mask])


def spearman_r(a: np.ndarray, b: np.ndarray) -> Tuple[float, float]:
    from scipy import stats as sp_stats
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return (np.nan, np.nan)
    return sp_stats.spearmanr(a[mask], b[mask])


def make_scatter(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
) -> None:
    mask = np.isfinite(x) & np.isfinite(y)
    xm, ym = x[mask], y[mask]
    ax.scatter(xm, ym, s=4, alpha=0.3, edgecolors="none")
    r_p, _ = pearson_r(x, y)
    rho, _ = spearman_r(x, y)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(f"{title}\nr={r_p:.3f}  ρ={rho:.3f}", fontsize=9)
    # Trend line
    if len(xm) > 2:
        z = np.polyfit(xm, ym, 1)
        p = np.poly1d(z)
        xs = np.linspace(xm.min(), xm.max(), 100)
        ax.plot(xs, p(xs), "r-", linewidth=1, alpha=0.7)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_benchmark(
    seqs: Dict[str, str],
    output_prefix: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(seqs)
    results: Dict[str, np.ndarray] = {}

    # --- Compute all metrics ---
    computations = [
        ("eprobe_tm",           eprobe_tm),
        ("eprobe_gc",           eprobe_gc),
        ("eprobe_hairpin",      eprobe_hairpin),
        ("eprobe_dust",         eprobe_dust),
        ("eprobe_dimer",        eprobe_dimer),
        ("primer3_tm",          primer3_tm),
        ("primer3_hairpin_dg",  primer3_hairpin_dg),
        ("primer3_homodimer_dg", primer3_homodimer_dg),
        ("vienna_mfe",         vienna_mfe),
        ("seqfold_mfe",        seqfold_mfe),
    ]

    for name, func in computations:
        logger.info(f"Computing {name} ({n} probes)...")
        t0 = time.time()
        results[name] = func(seqs)
        dt = time.time() - t0
        logger.info(f"  {name}: {dt:.1f}s")

    # --- Save raw data ---
    df = pd.DataFrame({"probe_id": list(seqs.keys())})
    for name, arr in results.items():
        df[name] = arr
    tsv_path = Path(str(output_prefix) + ".benchmark.tsv")
    df.to_csv(tsv_path, sep="\t", index=False)
    logger.info(f"Saved raw data to {tsv_path}")

    # --- Define comparison pairs ---
    comparisons = [
        # (x_key, y_key, xlabel, ylabel, title)
        ("eprobe_tm", "primer3_tm",
         "eProbe Tm (°C)", "primer3 Tm (°C)", "Tm: eProbe vs primer3"),
        ("eprobe_hairpin", "vienna_mfe",
         "eProbe Hairpin Score", "ViennaRNA MFE (kcal/mol)", "Hairpin: eProbe vs ViennaRNA"),
        ("eprobe_hairpin", "primer3_hairpin_dg",
         "eProbe Hairpin Score", "primer3 Hairpin ΔG (kcal/mol)", "Hairpin: eProbe vs primer3"),
        ("vienna_mfe", "primer3_hairpin_dg",
         "ViennaRNA MFE (kcal/mol)", "primer3 Hairpin ΔG (kcal/mol)", "ViennaRNA vs primer3 hairpin"),
        ("vienna_mfe", "seqfold_mfe",
         "ViennaRNA MFE (kcal/mol)", "seqfold MFE (kcal/mol)", "ViennaRNA vs seqfold (both Zuker)"),
        ("eprobe_dimer", "primer3_homodimer_dg",
         "eProbe Dimer Score", "primer3 Homodimer ΔG (kcal/mol)", "Dimer: eProbe vs primer3"),
        ("eprobe_dust", "vienna_mfe",
         "eProbe DUST Score", "ViennaRNA MFE (kcal/mol)", "Complexity vs secondary structure"),
    ]

    # --- Plot ---
    n_plots = len(comparisons)
    ncols = 3
    nrows = (n_plots + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows))
    axes_flat = axes.flatten() if n_plots > 1 else [axes]

    for idx, (xk, yk, xl, yl, title) in enumerate(comparisons):
        make_scatter(axes_flat[idx], results[xk], results[yk], xl, yl, title)

    # Hide unused axes
    for ax in axes_flat[n_plots:]:
        ax.set_visible(False)

    fig.suptitle(f"eProbe Biophysics Benchmark  (n={n})", fontsize=12, y=1.01)
    fig.tight_layout()
    plot_path = Path(str(output_prefix) + ".benchmark.png")
    fig.savefig(plot_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved scatter plots to {plot_path}")

    # --- Summary text ---
    summary_path = Path(str(output_prefix) + ".benchmark_summary.txt")
    with open(summary_path, "w") as f:
        f.write("eProbe Biophysics Benchmark Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Probes evaluated: {n}\n\n")
        f.write(f"{'Comparison':<44} {'Pearson r':>10} {'Spearman ρ':>10}\n")
        f.write("-" * 66 + "\n")
        for xk, yk, xl, yl, title in comparisons:
            r_p, _ = pearson_r(results[xk], results[yk])
            rho, _ = spearman_r(results[xk], results[yk])
            f.write(f"{title:<44} {r_p:>10.4f} {rho:>10.4f}\n")
        f.write("-" * 66 + "\n")

        # Per-metric descriptive stats
        f.write("\nDescriptive Statistics\n")
        f.write("-" * 60 + "\n")
        for name, arr in results.items():
            valid = arr[np.isfinite(arr)]
            if len(valid) == 0:
                continue
            f.write(f"\n{name}:\n")
            f.write(f"  Mean={np.mean(valid):.3f}  Std={np.std(valid):.3f}  "
                    f"Median={np.median(valid):.3f}  "
                    f"Min={np.min(valid):.3f}  Max={np.max(valid):.3f}\n")

    logger.info(f"Saved summary to {summary_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark eProbe biophysical metrics against primer3 & seqfold"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fasta", type=Path, help="Probe FASTA file")
    group.add_argument("--tsv", type=Path, help="eProbe selected TSV (needs --ref)")
    parser.add_argument("-r", "--ref", type=Path, default=None, help="Reference FASTA (for TSV input)")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output prefix")
    parser.add_argument("-n", "--max_probes", type=int, default=None,
                        help="Max probes to benchmark (random first N)")
    args = parser.parse_args()

    seqs = load_sequences(
        fasta_path=args.fasta,
        tsv_path=args.tsv,
        ref_path=args.ref,
        max_n=args.max_probes,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    run_benchmark(seqs, args.output)
    logger.info("Done.")


if __name__ == "__main__":
    main()
