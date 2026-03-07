#!/usr/bin/env python3
"""
Integration test for the optimized filter_biophysical function.

Creates synthetic SNP data and runs through the full biophysical filter pipeline.
"""

import sys
import random
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from dataclasses import dataclass
from typing import Optional, Dict, List

# Import the SNP model and filter
from eprobe.core.models import SNP
from eprobe.popgen.filter import filter_biophysical, BiophysicalThresholds


def generate_random_sequence(length: int = 81, gc_bias: float = 0.5) -> str:
    """
    Generate a random DNA sequence with controlled GC content.
    
    Args:
        length: Sequence length
        gc_bias: Probability of G or C (0.5 = 50% GC)
    """
    seq = []
    for _ in range(length):
        if random.random() < gc_bias:
            seq.append(random.choice('GC'))
        else:
            seq.append(random.choice('AT'))
    return ''.join(seq)


def create_mock_snps(n: int = 100) -> tuple:
    """
    Create mock SNP objects and probe sequences for testing.
    
    Returns:
        (list of SNP objects, dict of probe sequences)
    """
    snps = []
    probe_sequences: Dict[str, str] = {}
    
    for i in range(n):
        # Vary GC content to test filtering
        if i < n * 0.1:
            gc_bias = 0.2  # Low GC - should fail
        elif i < n * 0.2:
            gc_bias = 0.85  # High GC - should fail
        else:
            gc_bias = 0.5  # Normal GC - should pass GC filter
        
        # Generate probe sequence (81bp)
        probe_seq = generate_random_sequence(81, gc_bias)
        
        # Create SNP object (ref/alt are single bases)
        chrom = f"chr{i % 22 + 1}"
        pos = i * 1000 + 500
        ref_base = "A"
        alt_base = "G"
        
        snp = SNP.from_vcf_record(chrom=chrom, pos=pos, ref=ref_base, alt=alt_base)
        snps.append(snp)
        
        # Store probe sequence in dictionary
        probe_sequences[snp.id] = probe_seq
    
    # Add some sequences with low complexity (should fail DUST filter)
    for i in range(5):
        low_complexity_seq = "ATATATATAT" * 8 + "A"  # Very repetitive
        
        snp = SNP.from_vcf_record(
            chrom="chrX",
            pos=1000 + i * 100,
            ref="A",
            alt="T",
        )
        snps.append(snp)
        probe_sequences[snp.id] = low_complexity_seq
    
    # Add some palindromic sequences (should have high hairpin scores)
    # Create a hairpin-forming sequence
    stem = "ATCGATCGAT"
    loop = "NNNNN"
    rc_stem = "ATCGATCGAT"[::-1].translate(str.maketrans("ATCG", "TAGC"))
    palindrome = stem + loop + rc_stem
    # Pad to 81bp
    padding = generate_random_sequence(81 - len(palindrome), 0.5)
    palindrome = palindrome + padding
    
    for i in range(5):
        snp = SNP.from_vcf_record(
            chrom="chrY",
            pos=2000 + i * 100,
            ref="C",
            alt="G",
        )
        snps.append(snp)
        probe_sequences[snp.id] = palindrome
    
    return snps, probe_sequences


def test_filter_biophysical():
    """Test the filter_biophysical function."""
    print("=" * 60)
    print("  Integration Test: filter_biophysical")
    print("=" * 60)
    
    # Create test data
    random.seed(42)
    snps, probe_sequences = create_mock_snps(100)
    print(f"\nCreated {len(snps)} mock SNPs for testing")
    
    # Count expected failures
    low_gc_expected = 10
    high_gc_expected = 10
    low_complexity_expected = 5
    palindrome_expected = 5
    
    # Test 1: Default thresholds
    print("\n" + "-" * 40)
    print("Test 1: Default thresholds")
    print("-" * 40)
    
    thresholds = BiophysicalThresholds()
    print(f"GC: {thresholds.gc_min}-{thresholds.gc_max}%")
    print(f"Tm: {thresholds.tm_min}-{thresholds.tm_max}°C")
    print(f"DUST: ≤{thresholds.complexity_max}")
    print(f"Hairpin: {thresholds.hairpin}")
    print(f"Dimer: {thresholds.dimer}")
    
    result = filter_biophysical(snps, thresholds, probe_sequences)
    
    if result.is_ok():
        passed_snps, stats = result.value
        print(f"\nResult: {len(passed_snps)}/{len(snps)} SNPs passed")
        print(f"  GC failed: {stats['gc_failed']}")
        print(f"  Tm failed: {stats['tm_failed']}")
        print(f"  Complexity failed: {stats['complexity_failed']}")
        print(f"  Hairpin failed: {stats['hairpin_failed']}")
        print(f"  Dimer failed: {stats['dimer_failed']}")
        
        # Verify GC filtering worked
        assert stats['gc_failed'] >= low_gc_expected, "Low GC sequences should fail"
        print("  ✓ GC filtering works as expected")
    else:
        print(f"ERROR: {result.error}")
        return False
    
    # Test 2: Strict thresholds
    print("\n" + "-" * 40)
    print("Test 2: Strict thresholds (no hairpin/dimer)")
    print("-" * 40)
    
    strict_thresholds = BiophysicalThresholds(
        gc_min=40.0,
        gc_max=60.0,
        tm_min=60.0,
        tm_max=80.0,
        complexity_max=1.0,
        hairpin=0,  # Disable hairpin filter
        dimer=0,    # Disable dimer filter
    )
    
    result = filter_biophysical(snps, strict_thresholds, probe_sequences)
    
    if result.is_ok():
        passed_snps, stats = result.value
        print(f"\nResult: {len(passed_snps)}/{len(snps)} SNPs passed")
        print(f"  GC failed: {stats['gc_failed']}")
        print(f"  Tm failed: {stats['tm_failed']}")
        print(f"  Complexity failed: {stats['complexity_failed']}")
        print(f"  Hairpin failed: {stats['hairpin_failed']} (disabled)")
        print(f"  Dimer failed: {stats['dimer_failed']} (disabled)")
        
        # With stricter thresholds, more should fail
        assert len(passed_snps) < len(snps), "Strict thresholds should filter more"
        assert stats['hairpin_failed'] == 0, "Hairpin should be disabled"
        assert stats['dimer_failed'] == 0, "Dimer should be disabled"
        print("  ✓ Strict filtering works as expected")
    else:
        print(f"ERROR: {result.error}")
        return False
    
    # Test 3: Only GC filter
    print("\n" + "-" * 40)
    print("Test 3: GC-only filter (wide Tm/complexity)")
    print("-" * 40)
    
    gc_only = BiophysicalThresholds(
        gc_min=35.0,
        gc_max=65.0,
        tm_min=0.0,
        tm_max=200.0,
        complexity_max=100.0,
        hairpin=0,
        dimer=0,
    )
    
    result = filter_biophysical(snps, gc_only, probe_sequences)
    
    if result.is_ok():
        passed_snps, stats = result.value
        print(f"\nResult: {len(passed_snps)}/{len(snps)} SNPs passed")
        print(f"  Only GC failed: {stats['gc_failed']}")
        
        # Only GC should fail, not Tm/complexity
        expected_gc_failures = low_gc_expected + high_gc_expected
        assert stats['tm_failed'] == 0, "Tm should not fail with wide range"
        print("  ✓ GC-only filtering works as expected")
    else:
        print(f"ERROR: {result.error}")
        return False
    
    # Test 4: Empty input
    print("\n" + "-" * 40)
    print("Test 4: Empty input")
    print("-" * 40)
    
    result = filter_biophysical([], thresholds, {})
    
    if result.is_ok():
        passed_snps, stats = result.value
        assert len(passed_snps) == 0, "Empty input should return empty list"
        print("  ✓ Empty input handled correctly")
    else:
        print(f"ERROR: {result.error}")
        return False
    
    # Test 5: Large batch performance
    print("\n" + "-" * 40)
    print("Test 5: Performance test (1000 SNPs)")
    print("-" * 40)
    
    import time
    
    large_snps, large_probe_sequences = create_mock_snps(1000)
    
    start = time.time()
    result = filter_biophysical(large_snps, thresholds, large_probe_sequences)
    elapsed = time.time() - start
    
    if result.is_ok():
        passed_snps, stats = result.value
        snps_per_sec = len(large_snps) / elapsed
        print(f"\nProcessed {len(large_snps)} SNPs in {elapsed:.2f}s")
        print(f"  Throughput: {snps_per_sec:.0f} SNPs/sec")
        print(f"  Result: {len(passed_snps)} passed")
        
        # Should be faster than 100 SNPs/sec
        assert snps_per_sec > 100, f"Too slow: {snps_per_sec:.0f} SNPs/sec"
        print("  ✓ Performance is acceptable")
    else:
        print(f"ERROR: {result.error}")
        return False
    
    print("\n" + "=" * 60)
    print("  All integration tests PASSED!")
    print("=" * 60)
    return True


if __name__ == "__main__":
    success = test_filter_biophysical()
    sys.exit(0 if success else 1)
