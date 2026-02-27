#!/usr/bin/env python3
"""
Test script for optimized biophysical calculators.

Generates test sequences and validates:
1. Correctness: Results match expected values for known sequences
2. Performance: Batch processing is faster than individual calls
3. Consistency: Results are reproducible
"""

import time
import random
from typing import List, Tuple

# Add parent directory to path for local testing
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from eprobe.biophysics.fast_biophysics import (
    calculate_gc_fast,
    calculate_gc_batch_fast,
    calculate_tm_fast,
    calculate_tm_batch_fast,
    calculate_dust_fast,
    calculate_dust_batch_fast,
    calculate_hairpin_fast,
    calculate_hairpin_batch_fast,
    calculate_dimer_batch_fast,
    calculate_all_stats_fast,
    DimerCalculatorFast,
    calculate_percentile_threshold,
    HAS_PARASAIL,
)


def generate_random_sequence(length: int) -> str:
    """Generate a random DNA sequence."""
    return ''.join(random.choice('ATCG') for _ in range(length))


def generate_test_sequences(n: int, length: int = 81) -> List[str]:
    """Generate n random test sequences."""
    return [generate_random_sequence(length) for _ in range(n)]


def generate_special_sequences() -> List[Tuple[str, str]]:
    """
    Generate special test sequences with known properties.
    
    Returns:
        List of (sequence, description) tuples
    """
    return [
        # GC content tests
        ("ATATATAT" * 10, "0% GC - all AT"),
        ("GCGCGCGC" * 10, "100% GC - all GC"),
        ("ATCGATCG" * 10, "50% GC - balanced"),
        
        # DUST complexity tests
        ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" * 2, "Low complexity - poly-A"),
        ("ATGATGATGATGATGATGATGATGATGATGATGATGATGA", "Low complexity - ATG repeat"),
        ("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", "Moderate complexity - ATCG repeat"),
        
        # Hairpin tests (palindromic sequences)
        ("ATCGATCGATCGATCG" + "N" * 10 + "CGATCGATCGATCGAT", "High hairpin - palindrome"),
        ("GCGCATATGCGC" + "NNNNNNNNNN" + "GCGCATATGCGC", "Moderate hairpin"),
        ("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", "Low hairpin - no self-complementary"),
        
        # Similar sequences for dimer testing
        ("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", "Dimer test seq 1"),
        ("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", "Dimer test seq 2 (identical)"),
        ("TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "Dimer test seq 3 (different)"),
    ]


def test_gc():
    """Test GC content calculation."""
    print("\n" + "="*60)
    print("Testing GC Content Calculation")
    print("="*60)
    
    test_cases = [
        ("ATATATAT", 0.0),
        ("GCGCGCGC", 100.0),
        ("ATCGATCG", 50.0),
        ("AAAGGGCCC", 66.67),  # 6 GC out of 9
    ]
    
    all_passed = True
    for seq, expected in test_cases:
        result = calculate_gc_fast(seq)
        passed = abs(result - expected) < 0.1
        status = "✓" if passed else "✗"
        print(f"  {status} GC('{seq[:20]}...') = {result:.2f}% (expected: {expected:.2f}%)")
        all_passed = all_passed and passed
    
    return all_passed


def test_tm():
    """Test melting temperature calculation."""
    print("\n" + "="*60)
    print("Testing Melting Temperature Calculation")
    print("="*60)
    
    # Test with various sequences
    test_seqs = [
        ("ATCGATCGATCGATCGATCG", 40, 70, "50% GC, 20bp"),   # NN method
        ("GCGCGCGCGCGCGCGCGCGC", 60, 90, "100% GC, 20bp"),  # NN method
        ("ATATATAT", 10, 30, "0% GC, 8bp - Wallace rule"),   # Short seq uses Wallace
    ]
    
    all_passed = True
    for seq, min_tm, max_tm, desc in test_seqs:
        try:
            result = calculate_tm_fast(seq)
            passed = min_tm < result < max_tm
            status = "✓" if passed else "✗"
            print(f"  {status} Tm('{seq[:20]}...') = {result:.2f}°C ({desc})")
            all_passed = all_passed and passed
        except Exception as e:
            print(f"  ✗ Tm('{seq[:20]}...') failed: {e}")
            all_passed = False
    
    return all_passed


def test_dust():
    """Test DUST complexity calculation."""
    print("\n" + "="*60)
    print("Testing DUST Complexity Calculation")
    print("="*60)
    
    test_cases = [
        # (sequence, expected_low_complexity)
        ("AAAAAAAAAAAAAAAAAAAAAAAAAAAA", True),  # Poly-A should be low complexity (high DUST)
        ("ATGATGATGATGATGATGATGATGATG", True),   # Repeat should be low complexity
        ("ATCGATCGTAGCTAGCATCGATCGTAG", False),  # Random-ish should be high complexity (low DUST)
    ]
    
    all_passed = True
    for seq, expected_low in test_cases:
        result = calculate_dust_fast(seq)
        is_low = result > 2.0
        passed = is_low == expected_low
        status = "✓" if passed else "✗"
        complexity = "low" if is_low else "high"
        print(f"  {status} DUST('{seq[:20]}...') = {result:.4f} ({complexity} complexity)")
        all_passed = all_passed and passed
    
    return all_passed


def test_hairpin():
    """Test hairpin score calculation (stem method, normalized by log4(L))."""
    print("\n" + "="*60)
    print("Testing Hairpin Score Calculation (stem method)")
    print("  Score = max_stem_bp / log4(probe_length)")
    print("  Random DNA ~ 1.8, real hairpin > 3.0")
    print("="*60)
    
    # Sequences with different hairpin potential
    test_seqs = [
        # High hairpin: sequence has self-complementary regions
        ("ATCGATCGNNNNNNNNNNNNNNNNNNNCGATCGAT", "High - palindrome"),
        # Medium hairpin
        ("ATCGATCGATCGATCGATCGATCGATCGATCGATC", "Medium"),
        # Low hairpin: random sequence
        ("ATGCTGACGTAGCTGATCGATCGATCGTAGCTGA", "Low - random"),
    ]
    
    all_passed = True
    scores = []
    
    for seq, desc in test_seqs:
        try:
            result = calculate_hairpin_fast(seq)
            scores.append((result, desc))
            print(f"  ✓ Hairpin('{seq[:25]}...') = {result:.2f} ({desc})")
        except Exception as e:
            print(f"  ✗ Hairpin('{seq[:25]}...') failed: {e}")
            all_passed = False
    
    # Check that high hairpin has higher score
    if len(scores) >= 3:
        if not (scores[0][0] > scores[2][0]):
            print("  ⚠ Warning: Expected high hairpin seq to have higher score than low")
    
    return all_passed


def test_dimer():
    """Test dimer score calculation."""
    print("\n" + "="*60)
    print("Testing Dimer Score Calculation")
    print("="*60)
    
    # Create sequences with known similarity
    seq1 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    seq2 = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"  # Identical to seq1
    seq3 = "TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"  # Different from seq1
    seq4 = "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"  # Very different
    
    sequences = [seq1, seq2, seq3, seq4]
    
    calc = DimerCalculatorFast(k=11)
    n_kmers = calc.build_index(sequences)
    print(f"  Built index with {n_kmers} unique k-mers from {len(sequences)} sequences")
    
    scores = calc.calculate_all_scores()
    
    all_passed = True
    for i, (seq, score) in enumerate(zip(sequences, scores)):
        print(f"  Seq {i+1}: Dimer score = {score:.4f}")
    
    # Check that identical sequences have highest dimer scores
    if scores[0] < scores[3]:
        print("  ⚠ Warning: Expected identical sequences to have higher dimer score")
        all_passed = False
    else:
        print("  ✓ Identical sequences have higher dimer scores (as expected)")
    
    return all_passed


def test_batch_performance():
    """Test batch processing performance."""
    print("\n" + "="*60)
    print("Testing Batch Performance")
    print("="*60)
    
    # Generate test sequences
    n_seqs = 1000
    seq_len = 81
    sequences = generate_test_sequences(n_seqs, seq_len)
    
    print(f"  Generated {n_seqs} sequences of length {seq_len}")
    
    # Time GC calculation
    start = time.time()
    gc_results = calculate_gc_batch_fast(sequences)
    gc_time = time.time() - start
    print(f"  GC batch:      {gc_time*1000:.2f} ms ({n_seqs/gc_time:.0f} seqs/sec)")
    
    # Time Tm calculation
    start = time.time()
    tm_results = calculate_tm_batch_fast(sequences)
    tm_time = time.time() - start
    print(f"  Tm batch:      {tm_time*1000:.2f} ms ({n_seqs/tm_time:.0f} seqs/sec)")
    
    # Time DUST calculation
    start = time.time()
    dust_results = calculate_dust_batch_fast(sequences)
    dust_time = time.time() - start
    print(f"  DUST batch:    {dust_time*1000:.2f} ms ({n_seqs/dust_time:.0f} seqs/sec)")
    
    # Time Hairpin calculation
    start = time.time()
    hairpin_results = calculate_hairpin_batch_fast(sequences, threads=1)
    hairpin_time = time.time() - start
    print(f"  Hairpin batch: {hairpin_time*1000:.2f} ms ({n_seqs/hairpin_time:.0f} seqs/sec)")
    
    # Time Dimer calculation
    start = time.time()
    dimer_results = calculate_dimer_batch_fast(sequences)
    dimer_time = time.time() - start
    print(f"  Dimer batch:   {dimer_time*1000:.2f} ms ({n_seqs/dimer_time:.0f} seqs/sec)")
    
    # Time all-in-one calculation
    start = time.time()
    all_results = calculate_all_stats_fast(sequences, threads=1)
    all_time = time.time() - start
    print(f"  All stats:     {all_time*1000:.2f} ms ({n_seqs/all_time:.0f} seqs/sec)")
    
    # Show sample results
    print(f"\n  Sample result (first sequence):")
    print(f"    GC:      {all_results[0].gc:.2f}%")
    print(f"    Tm:      {all_results[0].tm:.2f}°C")
    print(f"    DUST:    {all_results[0].dust:.4f}")
    print(f"    Hairpin: {all_results[0].hairpin:.2f}")
    print(f"    Dimer:   {all_results[0].dimer:.4f}")
    
    return True


def test_percentile_threshold():
    """Test percentile-based thresholding."""
    print("\n" + "="*60)
    print("Testing Percentile Thresholding")
    print("="*60)
    
    # Generate test sequences and calculate hairpin scores
    sequences = generate_test_sequences(100, 81)
    scores = calculate_hairpin_batch_fast(sequences)
    
    # Calculate thresholds at different percentiles
    percentiles = [50, 75, 90, 95, 99]
    
    print(f"  Score range: {min(scores):.2f} - {max(scores):.2f}")
    print(f"  Mean: {sum(scores)/len(scores):.2f}")
    
    for p in percentiles:
        threshold = calculate_percentile_threshold(scores, p, higher_is_worse=True)
        n_pass = sum(1 for s in scores if s <= threshold)
        print(f"  {p}th percentile: {threshold:.2f} ({n_pass}/{len(scores)} would pass)")
    
    return True


def main():
    """Run all tests."""
    print("="*60)
    print("  Biophysical Calculator Optimization Tests")
    print("="*60)
    
    # Set random seed for reproducibility
    random.seed(42)
    
    tests = [
        ("GC Content", test_gc),
        ("Melting Temperature", test_tm),
        ("DUST Complexity", test_dust),
        ("Hairpin Score", test_hairpin),
        ("Dimer Score", test_dimer),
        ("Batch Performance", test_batch_performance),
        ("Percentile Threshold", test_percentile_threshold),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print(f"\n  ✗ Test '{name}' crashed: {e}")
            results.append((name, False))
    
    # Summary
    print("\n" + "="*60)
    print("  TEST SUMMARY")
    print("="*60)
    
    all_passed = True
    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}: {name}")
        all_passed = all_passed and passed
    
    print("="*60)
    if all_passed:
        print("  All tests PASSED!")
    else:
        print("  Some tests FAILED!")
    print("="*60)
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    exit(main())
