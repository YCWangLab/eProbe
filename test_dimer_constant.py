#!/usr/bin/env python3
"""
Test the dimer constant (100) to check if it produces reasonable score distributions.

Formula: D(g) = (100/N) * Σ c_i(g)
Where:
  - c_i(g) = frequency of kmer i from probe g in the pool
  - N = total number of probes
  - Constant 100 normalizes to percentage scale
"""

import numpy as np
from pathlib import Path
from typing import Dict
from collections import Counter

# Test on real probe set
test_fasta = Path("/Users/zh384/Desktop/scripts_dev/vs_code/developed_package/eProbe/tests/fixtures/test_seqs.fa")

def read_simple_fasta(filepath):
    """Simple FASTA parser."""
    seqs = {}
    current_id = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id is not None:
            seqs[current_id] = ''.join(current_seq)
    
    return seqs


def calc_dimer_legacy(sequences: Dict[str, str], k: int = 11):
    """Calculate dimer scores with different constants."""
    
    if len(sequences) <= 1:
        return {pid: 0.0 for pid in sequences.keys()}
    
    probe_ids = list(sequences.keys())
    probe_seqs = list(sequences.values())
    N = len(probe_seqs)
    
    # Build k-mer frequency dictionary
    kmer_freq: Dict[str, int] = {}
    total_kmers = 0
    for seq in probe_seqs:
        seq_upper = seq.upper()
        for i in range(len(seq_upper) - k + 1):
            kmer = seq_upper[i : i + k]
            if "N" not in kmer:
                kmer_freq[kmer] = kmer_freq.get(kmer, 0) + 1
                total_kmers += 1
    
    # Calculate scores with constant = 100
    scores_100 = {}
    scores_1 = {}
    scores_no_norm = {}
    
    for pid, seq in zip(probe_ids, probe_seqs):
        seq_upper = seq.upper()
        kmer_sum = 0
        
        for i in range(len(seq_upper) - k + 1):
            kmer = seq_upper[i : i + k]
            if "N" not in kmer and kmer in kmer_freq:
                kmer_sum += kmer_freq[kmer]
        
        # Three different normalizations
        score_100 = (kmer_sum / N) * 100.0  # Current: (100/N) * Σ
        score_1 = (kmer_sum / N) * 1.0      # Alternative: (1/N) * Σ
        score_no_norm = kmer_sum / (len(seq_upper) - k + 1)  # Per-kmer average
        
        scores_100[pid] = round(score_100, 2)
        scores_1[pid] = round(score_1, 2)
        scores_no_norm[pid] = round(score_no_norm, 2)
    
    return scores_100, scores_1, scores_no_norm, N, len(kmer_freq), total_kmers


def analyze_distribution(scores: Dict[str, float], label: str):
    """Analyze score distribution."""
    vals = list(scores.values())
    
    print(f"\n=== {label} ===")
    print(f"  Min:    {min(vals):.2f}")
    print(f"  Max:    {max(vals):.2f}")
    print(f"  Mean:   {np.mean(vals):.2f}")
    print(f"  Median: {np.median(vals):.2f}")
    print(f"  Std:    {np.std(vals):.2f}")
    print(f"  Range:  [{min(vals):.2f}, {max(vals):.2f}]")
    print(f"  P25:    {np.percentile(vals, 25):.2f}")
    print(f"  P75:    {np.percentile(vals, 75):.2f}")
    print(f"  P95:    {np.percentile(vals, 95):.2f}")
    
    # Check if values are in reasonable ranges
    if all(0 <= v <= 100 for v in vals):
        print(f"  ✓ All values in [0, 100] range (percentage scale)")
    elif all(0 <= v <= 1 for v in vals):
        print(f"  ✓ All values in [0, 1] range (normalized)")
    else:
        print(f"  ⚠ Values outside typical ranges")


if __name__ == '__main__':
    print("Testing dimer constant...")
    
    if not test_fasta.exists():
        print(f"❌ Test file not found: {test_fasta}")
        # Use synthetic data instead
        print("\nUsing synthetic data for testing...")
        sequences = {
            f"probe_{i}": "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
            for i in range(50)
        }
        # Add some variation
        for i in range(10, 20):
            sequences[f"probe_{i}"] = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"
    else:
        sequences = read_simple_fasta(test_fasta)
    
    print(f"Loaded {len(sequences)} sequences")
    
    scores_100, scores_1, scores_no_norm, N, n_kmers, total_kmers = calc_dimer_legacy(sequences)
    
    print(f"\nPool statistics:")
    print(f"  Total probes (N):    {N}")
    print(f"  Unique k-mers (k=11): {n_kmers}")
    print(f"  Total k-mer instances: {total_kmers}")
    print(f"  Avg k-mers per probe: {total_kmers / N:.1f}")
    print(f"  Avg unique k-mers per probe: {total_kmers / (N * (len(list(sequences.values())[0]) - 11 + 1)):.1f}")
    
    # Analyze all three approaches
    analyze_distribution(scores_100, "CONSTANT = 100  (D(g) = (100/N)*Σc_i)")
    analyze_distribution(scores_1, "CONSTANT = 1    (D(g) = (1/N)*Σc_i)")
    analyze_distribution(scores_no_norm, "NO DIVISION BY N (D(g) = Σc_i / n_kmers_in_probe)")
    
    print("\n" + "="*60)
    print("ANALYSIS:")
    print("="*60)
    
    vals_100 = list(scores_100.values())
    vals_1 = list(scores_1.values())
    
    print(f"\nWith CONSTANT=100:")
    print(f"  - Range is approximately [0, {max(vals_100):.1f}]")
    print(f"  - Interpretation: percentage of probe k-mers × total probe count")
    print(f"  - Scales with pool size N")
    
    print(f"\nWith CONSTANT=1:")
    print(f"  - Range is approximately [0, {max(vals_1):.1f}]")
    print(f"  - Interpretation: direct average k-mer frequency per probe")
    print(f"  - Independent of pool size N")
    
    print("\n" + "="*60)
    print("RECOMMENDATION:")
    print("="*60)
    
    # Check which is more reasonable
    if max(vals_100) > 1000:
        print("\n⚠ CONSTANT=100 produces very large values (>1000)")
        print("  → Consider using CONSTANT=1 or CONSTANT=10")
    elif max(vals_100) > 100:
        print("\n✓ CONSTANT=100 produces values in range [0-100] (percentage-like)")
        print("  → This is reasonable for interpretation")
    elif max(vals_100) > 10:
        print("\n✓ CONSTANT=100 produces values in range [0-10]")
        print("  → Small values, but reasonable")
    else:
        print("\n✓ CONSTANT=100 produces small normalized values [0-1]")
        print("  → Use as dimensionless score")
    
    print("\nKey insight:")
    print("  The constant controls the scale of the output.")
    print("  Choose based on:")
    print("  1. Interpretability: 100→percentage, 1→fraction")
    print("  2. Consistency with other metrics (GC%, Tm°C)")
