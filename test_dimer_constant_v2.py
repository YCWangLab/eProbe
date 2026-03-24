#!/usr/bin/env python3
"""
Test dimer constant with realistic probe data.

The issue: with constant=100, scores can become very large (>10000).

This test explores:
1. How the constant scales with pool size
2. What would be reasonable constants for different scenarios
3. Whether original formula needs adjustment
"""

import numpy as np
from typing import Dict
from collections import defaultdict


def generate_realistic_probes(n_probes: int = 1000, probe_length: int = 81):
    """Generate realistic probe sequences (semi-random, some repetition)."""
    bases = ['A', 'T', 'C', 'G']
    np.random.seed(42)
    
    sequences = {}
    for i in range(n_probes):
        # Add some correlation between probes using shared blocks
        if i < n_probes // 3:
            # Group 1: mostly A-T rich
            seq = "".join(np.random.choice(bases, size=probe_length, p=[0.4, 0.4, 0.1, 0.1]))
        elif i < 2 * n_probes // 3:
            # Group 2: mostly G-C rich
            seq = "".join(np.random.choice(bases, size=probe_length, p=[0.1, 0.1, 0.4, 0.4]))
        else:
            # Group 3: random
            seq = "".join(np.random.choice(bases, size=probe_length, p=[0.25, 0.25, 0.25, 0.25]))
        
        sequences[f"probe_{i:05d}"] = seq
    
    return sequences


def calc_dimer_scores(sequences: Dict[str, str], k: int = 11, constant: float = 100):
    """Calculate dimer scores with given constant."""
    
    if len(sequences) <= 1:
        return {pid: 0.0 for pid in sequences.keys()}
    
    probe_ids = list(sequences.keys())
    probe_seqs = list(sequences.values())
    N = len(probe_seqs)
    probe_len = len(probe_seqs[0])
    
    # Build k-mer frequency dictionary
    kmer_freq: Dict[str, int] = {}
    for seq in probe_seqs:
        seq_upper = seq.upper()
        for i in range(len(seq_upper) - k + 1):
            kmer = seq_upper[i : i + k]
            if "N" not in kmer:
                kmer_freq[kmer] = kmer_freq.get(kmer, 0) + 1
    
    n_kmers_in_pool = len(kmer_freq)
    
    # Calculate scores
    scores = {}
    for pid, seq in zip(probe_ids, probe_seqs):
        seq_upper = seq.upper()
        kmer_sum = 0
        n_kmers_in_probe = 0
        
        for i in range(len(seq_upper) - k + 1):
            kmer = seq_upper[i : i + k]
            if "N" not in kmer:
                n_kmers_in_probe += 1
                if kmer in kmer_freq:
                    kmer_sum += kmer_freq[kmer]
        
        # Original formula: (constant/N) * Σ c_i(g)
        score = (constant / N) * kmer_sum
        scores[pid] = round(score, 2)
    
    return scores, N, n_kmers_in_pool, probe_len


def analyze(scores: Dict[str, float], label: str):
    """Analyze score distribution."""
    vals = list(scores.values())
    print(f"\n{label}:")
    print(f"  Range:   [{min(vals):.2f}, {max(vals):.2f}]")
    print(f"  Mean:    {np.mean(vals):.2f}")
    print(f"  Median:  {np.median(vals):.2f}")
    print(f"  Std:     {np.std(vals):.2f}")
    print(f"  P5:      {np.percentile(vals, 5):.2f}")
    print(f"  P95:     {np.percentile(vals, 95):.2f}")


if __name__ == '__main__':
    print("="*70)
    print("Testing dimer constant with realistic probe pools")
    print("="*70)
    
    # Test with different pool sizes
    for n_probes in [100, 1000, 5000]:
        print(f"\n{'='*70}")
        print(f"SCENARIO: {n_probes} probes of 81bp (realistic eProbe set)")
        print(f"{'='*70}")
        
        sequences = generate_realistic_probes(n_probes, probe_length=81)
        
        # Test different constants
        print(f"\nFormula: D(g) = (C/N) × Σ c_i(g)")
        print(f"Where: N={n_probes}, C=constant, c_i(g)=k-mer frequency")
        
        scores_100, N, n_kmers, probe_len = calc_dimer_scores(sequences, constant=100)
        scores_1, _, _, _ = calc_dimer_scores(sequences, constant=1)
        scores_10, _, _, _ = calc_dimer_scores(sequences, constant=10)
        
        avg_kmer_in_probe = probe_len - 11 + 1  # ~71 for 81bp probe
        
        print(f"\nPool statistics:")
        print(f"  Total probes:           {N}")
        print(f"  Unique k-mers:          {n_kmers}")
        print(f"  K-mers per probe:       ~{avg_kmer_in_probe}")
        print(f"  Avg kmer frequency:     {sum(scores_100.values()) / len(scores_100) * N / 100:.1f}")
        
        print(f"\nConstant = 100 (D(g) = (100/N)*Σc_i)")
        analyze(scores_100, "  Scores")
        max_100 = max(scores_100.values())
        if max_100 > 1000:
            print(f"  ⚠️ MAX={max_100:.0f} - TOO LARGE! Not interpretable as percentage")
        elif max_100 > 100:
            print(f"  ⚠️ MAX={max_100:.0f} - Quite large, resembles percentages")
        else:
            print(f"  ✓ MAX={max_100:.0f} - Reasonable percentage-like scale")
        
        print(f"\nConstant = 10 (D(g) = (10/N)*Σc_i)")
        analyze(scores_10, "  Scores")
        
        print(f"\nConstant = 1 (D(g) = (1/N)*Σc_i)")
        analyze(scores_1, "  Scores")
    
    print(f"\n{'='*70}")
    print("ANALYSIS & RECOMMENDATIONS")
    print(f"{'='*70}")
    
    print("""
The constant 100 produces different scale behaviors depending on pool size:

1. Small pool (100 probes):
   - Scores will be 10-100x larger than more constrained pools
   - This reflects real high complementarity risk with small sets

2. Medium pool (1000 probes):
   - Scores scale as ~(100/1000) = 0.1× the sum of frequencies
   - Results in scores that CAN exceed 100 for repetitive probes

3. Large pool (10000+ probes):
   - Constant 100 produces small values (similar to using constant=1 on small pools)
   - Interpretation becomes more stable

ISSUE: The constant does NOT normalize against common reference points.

SOLUTIONS:
1. ✓ Keep constant=100 but accept it scales with pool size
   (This is what original code does - it's intentional!)
   
2. ✗ Normalize by probe length or expected k-mer overlap
   (Would change the semantic meaning)
   
3. ✗ Use a smaller constant (e.g., 1)
   (Loses the percentage interpretation)

ORIGINAL DESIGN LIKELY INTENDED:
- Larger pools → smaller scores (more diluted complementarity)
- Smaller pools → larger scores (higher complementarity risk)
- This is actually biologically sensible for multiplexing!

CONCLUSION:
=========
The constant 100 IS REASONABLE if you accept that:
- Dimer score SCALES with pool size
- Smaller probe sets have naturally higher dimer risk
- The metric reflects "frequency-based complementarity risk"

Alternative: If you want constant scale across pool sizes,
use constant=N (divide by N-1 instead).
Formula: D(g) = (N / (N-1)) × Σ c_i(g) / N = Σ c_i(g) / (N-1)
This gives "average k-mer frequency excluding self"
    """)
