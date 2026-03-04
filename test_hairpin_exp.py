#!/usr/bin/env python3
"""
Comprehensive test for exponential bonus hairpin scoring algorithm.

Tests across 52, 80, 100, 120bp probes with:
1. Random DNA baseline statistics (N=2000)
2. Designed hairpin stems (4-12bp) with non-complementary filler
3. Cross-length stability verification
4. Threshold analysis with stem-length interpretation

Output format: tables for academic paper.
"""

import random
import math
import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from eprobe.biophysics.fast_biophysics import (
    calculate_hairpin_fast,
    _calculate_hairpin_exponential,
)


def revcomp(seq):
    return seq[::-1].translate(str.maketrans('ATCG', 'TAGC'))


def make_hairpin_probe(stem_len, loop_len, total_len):
    """
    Build a probe with a KNOWN hairpin structure using non-self-complementary
    filler to avoid accidental extra stems.
    
    Structure: [filler]--STEM--LOOP--revcomp(STEM)--[filler]
    
    Filler uses "AACC" repeat which has minimal self-complementarity
    (revcomp(AACC) = GGTT, cannot match AACC).
    """
    # Fixed stem per stem_len (reproducible, non-palindromic)
    rng = random.Random(stem_len * 31 + 7)
    stem = ''.join(rng.choices('ATCG', k=stem_len))
    loop = 'T' * loop_len  # poly-T loop
    rc_stem = revcomp(stem)
    
    hairpin_core = stem + loop + rc_stem
    
    if len(hairpin_core) >= total_len:
        return hairpin_core[:total_len]
    
    # Non-self-complementary filler
    filler_needed = total_len - len(hairpin_core)
    left_pad = filler_needed // 2
    right_pad = filler_needed - left_pad
    
    filler_unit = "AACC"
    left = (filler_unit * (left_pad // 4 + 1))[:left_pad]
    right = (filler_unit * (right_pad // 4 + 1))[:right_pad]
    
    probe = left + hairpin_core + right
    assert len(probe) == total_len
    return probe


def raw_to_stem_bp(raw_score, k=4):
    """Convert raw exponential score to effective stem length in bp."""
    if raw_score <= 0:
        return 0
    n = math.log(raw_score, 4) + 1
    return round(n + k - 1, 1)


# ═══════════════════════════════════════════════════════════════════════
def print_scoring_reference():
    print("=" * 72)
    print("REFERENCE: Exponential Scoring Formula")
    print("  k=4, base=4. n consecutive k-mer matches → stem = n+3 bp")
    print("  Raw score = 4^(n-1), reported as MAX across all stems")
    print("=" * 72)
    
    print(f"\n{'Stem (bp)':>10}  {'n consec':>9}  {'Raw Score':>10}  "
          f"{'Norm(L=80)':>11}  {'Norm(L=120)':>12}")
    print("-" * 58)
    for stem_bp in [4, 5, 6, 7, 8, 9, 10, 11, 12]:
        n = stem_bp - 3
        raw = 4 ** (n - 1)
        norm80 = raw / math.log(80, 4)
        norm120 = raw / math.log(120, 4)
        print(f"{stem_bp:>10}  {n:>9}  {raw:>10,}  {norm80:>11.2f}  {norm120:>12.2f}")


# ═══════════════════════════════════════════════════════════════════════
def run_random_baseline(lengths, n_samples=2000, seed=42):
    print("\n" + "=" * 72)
    print(f"TABLE 1: Random DNA Baseline (N={n_samples} per length)")
    print("  Score = max(4^(n-1)) / log4(L)")
    print("=" * 72)
    
    print(f"\n{'Length':>8}  {'Mean':>7}  {'SD':>7}  {'P50':>7}  "
          f"{'P90':>7}  {'P95':>7}  {'P99':>7}  {'Max':>8}")
    print("-" * 72)
    
    results = {}
    for length in lengths:
        random.seed(seed)
        scores = []
        raws = []
        for _ in range(n_samples):
            seq = ''.join(random.choices('ATCG', k=length))
            raw = _calculate_hairpin_exponential(seq, k=4, min_loop=3)
            norm = calculate_hairpin_fast(seq, method="exp")
            scores.append(norm)
            raws.append(raw)
        
        scores.sort()
        raws.sort()
        mean = sum(scores) / len(scores)
        std = (sum((s - mean) ** 2 for s in scores) / len(scores)) ** 0.5
        p50 = scores[len(scores) // 2]
        p90 = scores[int(len(scores) * 0.90)]
        p95 = scores[int(len(scores) * 0.95)]
        p99 = scores[int(len(scores) * 0.99)]
        mx = scores[-1]
        
        raw_p50 = raws[len(raws) // 2]
        raw_p90 = raws[int(len(raws) * 0.90)]
        raw_p95 = raws[int(len(raws) * 0.95)]
        raw_p99 = raws[int(len(raws) * 0.99)]
        
        print(f"{length:>6}bp  {mean:>7.2f}  {std:>7.2f}  {p50:>7.2f}  "
              f"{p90:>7.2f}  {p95:>7.2f}  {p99:>7.2f}  {mx:>8.2f}")
        
        results[length] = {
            "mean": mean, "std": std, "p50": p50,
            "p90": p90, "p95": p95, "p99": p99, "max": mx,
            "raw_p50": raw_p50, "raw_p90": raw_p90,
            "raw_p95": raw_p95, "raw_p99": raw_p99,
        }
    
    print(f"\n  Raw score → effective stem length:")
    for length in lengths:
        r = results[length]
        print(f"    {length}bp: P50 raw={r['raw_p50']:.0f}→{raw_to_stem_bp(r['raw_p50'])}bp, "
              f"P90 raw={r['raw_p90']:.0f}→{raw_to_stem_bp(r['raw_p90'])}bp, "
              f"P99 raw={r['raw_p99']:.0f}→{raw_to_stem_bp(r['raw_p99'])}bp stem")
    
    return results


# ═══════════════════════════════════════════════════════════════════════
def run_hairpin_detection(lengths, stem_lengths, loop_len=4):
    print("\n" + "=" * 72)
    print(f"TABLE 2: Known Hairpin Detection (loop={loop_len}bp, AACC filler)")
    print("=" * 72)
    
    header = f"\n{'Stem':>6}"
    for L in lengths:
        header += f"  {'Raw':>8}  {'Norm':>8}"
    print(header)
    
    sub = f"{'(bp)':>6}"
    for L in lengths:
        sub += f"  {f'({L}bp)':>8}  {f'({L}bp)':>8}"
    print(sub)
    print("-" * (6 + len(lengths) * 18))
    
    results = {}
    for stem_len in stem_lengths:
        row = f"{stem_len:>4}bp"
        results[stem_len] = {}
        
        for L in lengths:
            probe = make_hairpin_probe(stem_len, loop_len, L)
            raw = _calculate_hairpin_exponential(probe, k=4, min_loop=3)
            norm = calculate_hairpin_fast(probe, method="exp")
            row += f"  {raw:>8.0f}  {norm:>8.2f}"
            results[stem_len][L] = {"raw": raw, "norm": norm}
        
        print(row)
    
    # Verify
    print(f"\n  Expected raw = 4^(stem-4). Verification:")
    for stem_len in stem_lengths:
        n = max(stem_len - 3, 1)
        expected = 4 ** (n - 1)
        actuals = [results[stem_len][L]["raw"] for L in lengths]
        match = all(a == expected for a in actuals)
        status = "✓ exact" if match else f"≠ got {actuals}"
        print(f"    {stem_len}bp stem: expected {expected}, {status}")
    
    return results


# ═══════════════════════════════════════════════════════════════════════
def run_stability_test(lengths):
    print("\n" + "=" * 72)
    print("TABLE 3: Cross-length Stability")
    print("  Same hairpin (10bp stem) in different probe lengths")
    print("  Old method: raw_stem/log4(L) → DILUTED in longer probes")
    print("  New method: 4^(n-1)/log4(L) → STABLE, always caught")
    print("=" * 72)
    
    stem_len = 10
    print(f"\n{'Length':>8}  {'Old Stem':>10}  {'Old Caught':>11}  "
          f"{'New Exp':>10}  {'New Caught':>11}")
    print("-" * 58)
    
    for L in lengths:
        probe = make_hairpin_probe(stem_len, 4, L)
        old = calculate_hairpin_fast(probe, method="stem")
        new = calculate_hairpin_fast(probe, method="exp")
        old_caught = "YES" if old > 3.0 else "NO ← miss!"
        new_caught = "YES" if new > 20.0 else "NO ← miss!"
        print(f"{L:>6}bp  {old:>10.2f}  {old_caught:>11}  "
              f"{new:>10.2f}  {new_caught:>11}")


# ═══════════════════════════════════════════════════════════════════════
def run_threshold_recommendation(baseline, hairpin, lengths):
    print("\n" + "=" * 72)
    print("TABLE 4: Threshold Recommendation")
    print("=" * 72)
    
    print(f"\n  At threshold = 20.0 (recommended):")
    print(f"  {'Length':>8}  {'Pass random':>12}  "
          f"{'5bp':>8}  {'6bp':>8}  {'7bp':>8}  {'8bp':>8}  {'10bp':>8}")
    print(f"  {'-' * 62}")
    
    for L in lengths:
        b = baseline[L]
        pr = ">99%" if 20 > b["p99"] else ">95%" if 20 > b["p95"] else \
             ">90%" if 20 > b["p90"] else "<90%"
        
        cells = []
        for s in [5, 6, 7, 8, 10]:
            v = hairpin.get(s, {}).get(L, {}).get("norm", 0)
            mark = "✗" if v > 20 else "✓"
            cells.append(f"{v:.1f}{mark}")
        
        print(f"  {L:>6}bp  {pr:>12}  " + "  ".join(f"{c:>8}" for c in cells))
    
    print("""
  ✓ = passes (kept), ✗ = filtered (removed)
  
  Physical interpretation:
  │ Threshold │ Raw equiv │ Min stem caught │ Biology                    │
  │───────────│───────────│─────────────────│────────────────────────────│
  │    5.0    │   ~16     │   6bp           │ Marginally stable hairpin  │
  │   20.0    │   ~64     │   7bp           │ Thermodynamically stable   │
  │   80.0    │   ~256    │   8bp           │ Very stable hairpin only   │
  │  300.0    │   ~1024   │   9bp           │ Extreme hairpins only      │
  
  Default = 20.0 → filters stems ≥ 7bp (stable hairpins that compete
  with target hybridization at typical capture temperatures of 55-65°C)""")


# ═══════════════════════════════════════════════════════════════════════
def run_benchmark(n_probes=10000, length=81):
    print(f"\n{'=' * 72}")
    print(f"BENCHMARK: {n_probes:,} probes x {length}bp")
    print(f"{'=' * 72}")
    
    random.seed(42)
    probes = [''.join(random.choices('ATCG', k=length)) for _ in range(n_probes)]
    
    start = time.perf_counter()
    for p in probes:
        calculate_hairpin_fast(p, method="exp")
    exp_time = time.perf_counter() - start
    
    start = time.perf_counter()
    for p in probes:
        calculate_hairpin_fast(p, method="exp", early_stop_raw=64)
    exp_es_time = time.perf_counter() - start
    
    start = time.perf_counter()
    for p in probes:
        calculate_hairpin_fast(p, method="stem")
    stem_time = time.perf_counter() - start
    
    print(f"  Exp method:           {exp_time:.3f}s ({n_probes/exp_time:,.0f} probes/sec)")
    print(f"  Exp + early stop:     {exp_es_time:.3f}s ({n_probes/exp_es_time:,.0f} probes/sec)")
    print(f"  Stem method (legacy): {stem_time:.3f}s ({n_probes/stem_time:,.0f} probes/sec)")
    
    print(f"\n  Time complexity: O(L + M)")
    print(f"  L={length}, M (expected matches) ≈ L²/4^k = {(length-3)**2 // 256}")
    print(f"  Early stopping avoids full scan for obviously bad probes")


if __name__ == "__main__":
    lengths = [52, 80, 100, 120]
    stem_lengths = [4, 5, 6, 7, 8, 10, 12]
    
    print_scoring_reference()
    baseline = run_random_baseline(lengths, n_samples=2000)
    hairpin = run_hairpin_detection(lengths, stem_lengths)
    run_stability_test(lengths)
    run_threshold_recommendation(baseline, hairpin, lengths)
    run_benchmark()
