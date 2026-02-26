#!/usr/bin/env python3
"""
Comprehensive tests for SNP selection module.

Tests all selection strategies:
  - random: Simple random sampling
  - uniform: Window-based uniform distribution
  - weighted: Biophysical score-based selection
  - priority: BED region prioritization
"""

import tempfile
import shutil
from pathlib import Path
from typing import List, Dict
import random

import numpy as np

from eprobe.core.models import SNP, SNPDataFrame
from eprobe.popgen.select import (
    SelectionStrategy,
    SelectionConfig,
    calculate_window_index,
    parse_bed_file,
    snp_in_priority_region,
    calculate_snp_score,
    select_uniform,
    select_random,
    select_weighted,
    select_by_priority,
    select_chromosomes,
    run_select,
    calculate_coverage_stats,
)


# =============================================================================
# Test Data Generation
# =============================================================================

def generate_test_snps(
    n_snps: int = 1000,
    n_chroms: int = 5,
    chrom_length: int = 1000000,
    seed: int = 42,
) -> List[SNP]:
    """Generate test SNP data."""
    random.seed(seed)
    np.random.seed(seed)
    
    snps = []
    bases = ['A', 'C', 'G', 'T']
    
    for i in range(n_snps):
        chrom = f"chr{(i % n_chroms) + 1}"
        pos = random.randint(1, chrom_length)
        ref = random.choice(bases)
        alt = random.choice([b for b in bases if b != ref])
        
        snp = SNP.from_vcf_record(chrom=chrom, pos=pos, ref=ref, alt=alt)
        snps.append(snp)
    
    # Sort by chrom and pos
    snps.sort(key=lambda s: (s.chrom, s.pos))
    return snps


def create_test_bed_file(
    path: Path,
    snps: List[SNP],
    coverage_fraction: float = 0.3,
    region_size: int = 5000,
    seed: int = 42,
) -> None:
    """Create a BED file covering some SNP positions."""
    random.seed(seed)
    
    # Select some SNPs to create regions around
    n_regions = int(len(snps) * coverage_fraction)
    selected_snps = random.sample(snps, min(n_regions, len(snps)))
    
    regions = []
    for snp in selected_snps:
        start = max(0, snp.pos - region_size // 2)
        end = snp.pos + region_size // 2
        regions.append((snp.chrom, start, end))
    
    # Sort and write
    regions.sort()
    
    with open(path, 'w') as f:
        for chrom, start, end in regions:
            f.write(f"{chrom}\t{start}\t{end}\n")


def create_test_tsv(path: Path, snps: List[SNP]) -> None:
    """Create test SNP TSV file."""
    df = SNPDataFrame.from_snps(snps)
    df.to_tsv(path)


def create_test_tsv_with_biophysical(path: Path, snps: List[SNP], seed: int = 42) -> None:
    """Create test SNP TSV file with biophysical columns."""
    np.random.seed(seed)
    
    data = {
        'chr': [s.chrom for s in snps],
        'pos': [s.pos for s in snps],
        'ref': [s.ref for s in snps],
        'alt': [s.alt for s in snps],
        'type': [s.mutation_type for s in snps],
        'gc': [40.0 + np.random.uniform(-10, 10) for _ in snps],
        'tm': [65.0 + np.random.uniform(-10, 10) for _ in snps],
        'complexity': [np.random.uniform(0, 3) for _ in snps],
        'hairpin': [np.random.uniform(0, 0.5) for _ in snps],
        'dimer': [np.random.uniform(0, 0.5) for _ in snps],
    }
    
    import pandas as pd
    df = pd.DataFrame(data)
    df.to_csv(path, sep='\t', index=False)


# =============================================================================
# Unit Tests
# =============================================================================

def test_calculate_window_index():
    """Test window index calculation."""
    print("Testing calculate_window_index...")
    
    assert calculate_window_index(0, 10000) == 0
    assert calculate_window_index(9999, 10000) == 0
    assert calculate_window_index(10000, 10000) == 1
    assert calculate_window_index(25000, 10000) == 2
    assert calculate_window_index(100000, 10000) == 10
    
    print("  ✓ Window index calculation correct")


def test_parse_bed_file():
    """Test BED file parsing."""
    print("Testing parse_bed_file...")
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write("chr1\t100\t200\n")
        f.write("chr1\t500\t600\n")
        f.write("chr2\t1000\t2000\n")
        f.write("# comment line\n")
        f.write("track name=test\n")
        bed_path = Path(f.name)
    
    try:
        result = parse_bed_file(bed_path)
        assert result.is_ok()
        regions = result.unwrap()
        
        assert len(regions) == 3
        assert ("chr1", 100, 200) in regions
        assert ("chr1", 500, 600) in regions
        assert ("chr2", 1000, 2000) in regions
        
        print(f"  ✓ Parsed {len(regions)} regions correctly")
    finally:
        bed_path.unlink()


def test_snp_in_priority_region():
    """Test SNP-in-region checking."""
    print("Testing snp_in_priority_region...")
    
    regions = {
        ("chr1", 100, 200),
        ("chr1", 500, 600),
        ("chr2", 1000, 2000),
    }
    
    # SNP inside region
    snp1 = SNP.from_vcf_record(chrom="chr1", pos=150, ref="A", alt="G")
    assert snp_in_priority_region(snp1, regions) == True
    
    # SNP at region boundary
    snp2 = SNP.from_vcf_record(chrom="chr1", pos=100, ref="A", alt="G")
    assert snp_in_priority_region(snp2, regions) == True
    
    # SNP outside region
    snp3 = SNP.from_vcf_record(chrom="chr1", pos=300, ref="A", alt="G")
    assert snp_in_priority_region(snp3, regions) == False
    
    # SNP on different chromosome
    snp4 = SNP.from_vcf_record(chrom="chr3", pos=150, ref="A", alt="G")
    assert snp_in_priority_region(snp4, regions) == False
    
    print("  ✓ Region checking works correctly")


def test_calculate_snp_score():
    """Test biophysical score calculation."""
    print("Testing calculate_snp_score...")
    
    # SNP without tags (should use defaults)
    snp = SNP.from_vcf_record(chrom="chr1", pos=100, ref="A", alt="G")
    
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
    score = calculate_snp_score(snp, weights)
    
    # Should get a reasonable score with defaults
    assert 0 <= score <= 5.0, f"Score {score} out of expected range"
    
    print(f"  ✓ Default SNP score: {score:.3f}")


def test_select_random():
    """Test random selection."""
    print("Testing select_random...")
    
    snps = generate_test_snps(n_snps=1000, seed=42)
    
    # Select 100 SNPs
    result = select_random(snps, target_count=100, seed=42)
    assert result.is_ok()
    selected = result.unwrap()
    
    assert len(selected) == 100
    assert len(set(s.id for s in selected)) == 100  # All unique
    
    # Reproducibility test
    result2 = select_random(snps, target_count=100, seed=42)
    selected2 = result2.unwrap()
    assert [s.id for s in selected] == [s.id for s in selected2]
    
    print(f"  ✓ Selected {len(selected)} SNPs randomly")
    print(f"  ✓ Reproducibility verified")


def test_select_uniform():
    """Test uniform (window-based) selection."""
    print("Testing select_uniform...")
    
    snps = generate_test_snps(n_snps=1000, n_chroms=5, seed=42)
    
    # Select with 10kb windows
    result = select_uniform(snps, target_count=100, window_size=10000, seed=42)
    assert result.is_ok()
    selected = result.unwrap()
    
    # Check approximate target
    assert 50 <= len(selected) <= 150, f"Expected ~100 SNPs, got {len(selected)}"
    
    # Check distribution across chromosomes
    chrom_counts = {}
    for snp in selected:
        chrom_counts[snp.chrom] = chrom_counts.get(snp.chrom, 0) + 1
    
    print(f"  ✓ Selected {len(selected)} SNPs uniformly")
    print(f"  ✓ Distribution: {chrom_counts}")
    
    # At least some chromosomes should be represented
    assert len(chrom_counts) >= 1


def test_select_weighted():
    """Test weighted (biophysical score) selection."""
    print("Testing select_weighted...")
    
    snps = generate_test_snps(n_snps=1000, seed=42)
    
    # Create DataFrame with biophysical columns
    data = {
        'chr': [s.chrom for s in snps],
        'pos': [s.pos for s in snps],
        'ref': [s.ref for s in snps],
        'alt': [s.alt for s in snps],
        'type': [s.mutation_type for s in snps],
        'gc': [40.0 + np.random.uniform(-10, 10) for _ in snps],
        'tm': [65.0 + np.random.uniform(-10, 10) for _ in snps],
        'complexity': [np.random.uniform(0, 3) for _ in snps],
        'hairpin': [np.random.uniform(0, 0.5) for _ in snps],
        'dimer': [np.random.uniform(0, 0.5) for _ in snps],
    }
    import pandas as pd
    snp_df = pd.DataFrame(data)
    
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
    
    result = select_weighted(snps, target_count=100, weights=weights, window_size=10000, seed=42, snp_df=snp_df)
    assert result.is_ok()
    selected_snps, selected_df = result.unwrap()
    
    assert len(selected_snps) > 0
    # Verify biophysical columns are preserved
    assert 'gc' in selected_df.columns
    assert 'tm' in selected_df.columns
    
    print(f"  ✓ Selected {len(selected_snps)} SNPs by score")
    print(f"  ✓ Biophysical columns preserved in output")


def test_select_weighted_requires_biophysical():
    """Test that weighted selection fails without biophysical columns."""
    print("Testing select_weighted requires biophysical columns...")
    
    snps = generate_test_snps(n_snps=100, seed=42)
    
    # No DataFrame provided
    result = select_weighted(snps, target_count=50, weights=[1.0]*5, window_size=10000, seed=42)
    assert result.is_err()
    assert "biophysical" in result.unwrap_err().lower()
    
    # DataFrame without biophysical columns
    import pandas as pd
    snp_df = pd.DataFrame({
        'chr': [s.chrom for s in snps],
        'pos': [s.pos for s in snps],
        'ref': [s.ref for s in snps],
        'alt': [s.alt for s in snps],
        'type': [s.mutation_type for s in snps],
    })
    
    result = select_weighted(snps, target_count=50, weights=[1.0]*5, window_size=10000, seed=42, snp_df=snp_df)
    assert result.is_err()
    assert "gc" in result.unwrap_err() or "biophysical" in result.unwrap_err().lower()
    
    print("  ✓ Correctly requires biophysical columns")


def test_select_by_priority():
    """Test priority (BED region) selection."""
    print("Testing select_by_priority...")
    
    snps = generate_test_snps(n_snps=1000, n_chroms=5, seed=42)
    
    # Create priority regions covering ~30% of SNPs
    regions = set()
    for snp in snps[:300]:  # First 300 SNPs get regions
        start = max(0, snp.pos - 100)
        end = snp.pos + 100
        regions.add((snp.chrom, start, end))
    
    result = select_by_priority(snps, target_count=100, priority_regions=regions, window_size=10000, seed=42)
    assert result.is_ok()
    selected = result.unwrap()
    
    # Count how many selected are in priority regions
    in_priority = sum(1 for s in selected if snp_in_priority_region(s, regions))
    
    print(f"  ✓ Selected {len(selected)} SNPs")
    print(f"  ✓ In priority regions: {in_priority} ({in_priority/len(selected)*100:.1f}%)")
    
    # Should prioritize SNPs in regions
    assert in_priority > 0


def test_select_priority_with_weights():
    """Test priority selection combined with biophysical weights."""
    print("Testing select_by_priority with weights...")
    
    snps = generate_test_snps(n_snps=1000, n_chroms=5, seed=42)
    
    # Create DataFrame with biophysical columns
    import pandas as pd
    np.random.seed(42)
    snp_df = pd.DataFrame({
        'chr': [s.chrom for s in snps],
        'pos': [s.pos for s in snps],
        'ref': [s.ref for s in snps],
        'alt': [s.alt for s in snps],
        'type': [s.mutation_type for s in snps],
        'gc': [40.0 + np.random.uniform(-10, 10) for _ in snps],
        'tm': [65.0 + np.random.uniform(-10, 10) for _ in snps],
        'complexity': [np.random.uniform(0, 3) for _ in snps],
        'hairpin': [np.random.uniform(0, 0.5) for _ in snps],
        'dimer': [np.random.uniform(0, 0.5) for _ in snps],
    })
    
    # Create priority regions
    regions = set()
    for snp in snps[:300]:
        start = max(0, snp.pos - 100)
        end = snp.pos + 100
        regions.add((snp.chrom, start, end))
    
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
    
    result = select_by_priority(
        snps, target_count=100, priority_regions=regions,
        window_size=10000, seed=42, weights=weights, snp_df=snp_df
    )
    assert result.is_ok()
    
    # Should return tuple when weights are used
    unwrapped = result.unwrap()
    assert isinstance(unwrapped, tuple), "Should return (snps, df) when using weights"
    selected_snps, selected_df = unwrapped
    
    assert len(selected_snps) > 0
    assert 'gc' in selected_df.columns
    
    in_priority = sum(1 for s in selected_snps if snp_in_priority_region(s, regions))
    print(f"  ✓ Selected {len(selected_snps)} SNPs with weighted priority")
    print(f"  ✓ In priority regions: {in_priority} ({in_priority/len(selected_snps)*100:.1f}%)")
    print("  ✓ Biophysical columns preserved in output DataFrame")


def test_select_chromosomes():
    """Test chromosome filtering."""
    print("Testing select_chromosomes...")
    
    snps = generate_test_snps(n_snps=1000, n_chroms=5, seed=42)
    
    # Filter to chr1 and chr3
    result = select_chromosomes(snps, chromosomes=["chr1", "chr3"])
    assert result.is_ok()
    filtered = result.unwrap()
    
    chroms = set(s.chrom for s in filtered)
    assert chroms == {"chr1", "chr3"}
    
    print(f"  ✓ Filtered to {len(filtered)} SNPs on chr1, chr3")


def test_run_select_full_pipeline():
    """Test the full selection pipeline."""
    print("\nTesting run_select full pipeline...")
    
    # Create temporary directory
    temp_dir = Path(tempfile.mkdtemp())
    
    try:
        # Generate test data
        snps = generate_test_snps(n_snps=500, n_chroms=3, seed=42)
        input_path = temp_dir / "test_snps.tsv"
        create_test_tsv(input_path, snps)
        
        # Create a separate input with biophysical columns for weighted test
        input_path_bio = temp_dir / "test_snps_biophysical.tsv"
        create_test_tsv_with_biophysical(input_path_bio, snps, seed=42)
        
        output_prefix = temp_dir / "output"
        
        # Test random strategy
        print("\n  Testing random strategy...")
        result = run_select(
            input_path=input_path,
            output_prefix=output_prefix,
            strategy="random",
            target_count=100,
            window_size=10000,
            seed=42,
        )
        assert result.is_ok(), f"Random selection failed: {result.unwrap_err()}"
        stats = result.unwrap()
        print(f"    ✓ Random: {stats['selected']} SNPs selected")
        
        # Test uniform strategy
        print("\n  Testing uniform strategy...")
        result = run_select(
            input_path=input_path,
            output_prefix=temp_dir / "output_uniform",
            strategy="uniform",
            target_count=100,
            window_size=10000,
            seed=42,
        )
        assert result.is_ok(), f"Uniform selection failed: {result.unwrap_err()}"
        stats = result.unwrap()
        print(f"    ✓ Uniform: {stats['selected']} SNPs, {stats['windows']} windows")
        
        # Test weighted strategy (requires biophysical columns)
        print("\n  Testing weighted strategy...")
        result = run_select(
            input_path=input_path_bio,  # Use file WITH biophysical columns
            output_prefix=temp_dir / "output_weighted",
            strategy="weighted",
            target_count=100,
            window_size=10000,
            weights=[1.0, 1.0, 1.0, 1.0, 1.0],
            seed=42,
            keep_biophysical=False,  # Default: don't keep biophysical columns
        )
        assert result.is_ok(), f"Weighted selection failed: {result.unwrap_err()}"
        stats = result.unwrap()
        print(f"    ✓ Weighted: {stats['selected']} SNPs")
        
        # Verify biophysical columns are NOT preserved by default
        import pandas as pd
        output_weighted = temp_dir / "output_weighted.selected.tsv"
        df_out = pd.read_csv(output_weighted, sep='\t')
        assert 'gc' not in df_out.columns, "Biophysical columns should NOT be in output by default"
        print("    ✓ Biophysical columns correctly removed (default)")
        
        # Test with keep_biophysical=True
        result = run_select(
            input_path=input_path_bio,
            output_prefix=temp_dir / "output_weighted_bio",
            strategy="weighted",
            target_count=100,
            window_size=10000,
            weights=[1.0, 1.0, 1.0, 1.0, 1.0],
            seed=42,
            keep_biophysical=True,  # Explicitly keep biophysical columns
        )
        assert result.is_ok()
        output_weighted_bio = temp_dir / "output_weighted_bio.selected.tsv"
        df_out_bio = pd.read_csv(output_weighted_bio, sep='\t')
        assert 'gc' in df_out_bio.columns, "Biophysical columns should be preserved with keep_biophysical=True"
        print("    ✓ Biophysical columns preserved with --keep_biophysical")
        
        # Test priority strategy
        print("\n  Testing priority strategy...")
        bed_path = temp_dir / "priority.bed"
        create_test_bed_file(bed_path, snps, coverage_fraction=0.3, seed=42)
        
        result = run_select(
            input_path=input_path,
            output_prefix=temp_dir / "output_priority",
            strategy="priority",
            target_count=100,
            window_size=10000,
            priority_bed=bed_path,
            seed=42,
        )
        assert result.is_ok(), f"Priority selection failed: {result.unwrap_err()}"
        stats = result.unwrap()
        print(f"    ✓ Priority: {stats['selected']} SNPs")
        
        # Verify output file exists
        assert Path(stats['output_file']).exists()
        print(f"    ✓ Output file created: {stats['output_file']}")
        
        # Test priority + weights combination
        print("\n  Testing priority + weights combination...")
        result = run_select(
            input_path=input_path_bio,  # Use biophysical data
            output_prefix=temp_dir / "output_priority_weighted",
            strategy="priority",
            target_count=100,
            window_size=10000,
            priority_bed=bed_path,
            weights=[1.0, 1.0, 1.0, 1.0, 1.0],
            seed=42,
            keep_biophysical=True,
        )
        assert result.is_ok(), f"Priority+weighted failed: {result.unwrap_err()}"
        stats = result.unwrap()
        assert "weighted" in stats['strategy'].lower(), "Strategy should indicate weighted was used"
        print(f"    ✓ Priority+weighted: {stats['selected']} SNPs, strategy={stats['strategy']}")
        
        # Verify biophysical columns preserved
        output_pw = temp_dir / "output_priority_weighted.selected.tsv"
        df_pw = pd.read_csv(output_pw, sep='\t')
        assert 'gc' in df_pw.columns, "Biophysical columns should be preserved"
        print("    ✓ Biophysical columns preserved in priority+weighted output")
        
    finally:
        # Cleanup
        shutil.rmtree(temp_dir)
        print("\n  ✓ Cleanup complete")


def test_coverage_stats():
    """Test coverage statistics calculation."""
    print("Testing calculate_coverage_stats...")
    
    snps = generate_test_snps(n_snps=100, n_chroms=3, seed=42)
    
    stats = calculate_coverage_stats(snps, window_size=10000)
    
    assert 'windows_covered' in stats
    assert 'chromosomes' in stats
    assert 'per_chromosome' in stats
    assert 'mean_per_chromosome' in stats
    assert 'std_per_chromosome' in stats
    
    print(f"  ✓ Windows covered: {stats['windows_covered']}")
    print(f"  ✓ Chromosomes: {stats['chromosomes']}")
    print(f"  ✓ Mean per chromosome: {stats['mean_per_chromosome']:.1f}")


# =============================================================================
# Main Test Runner
# =============================================================================

def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("SNP Selection Module Tests")
    print("=" * 60)
    
    # Unit tests
    test_calculate_window_index()
    test_parse_bed_file()
    test_snp_in_priority_region()
    test_calculate_snp_score()
    
    # Selection strategy tests
    print("\n" + "-" * 40)
    print("Selection Strategy Tests")
    print("-" * 40)
    
    test_select_random()
    test_select_uniform()
    test_select_weighted()
    test_select_weighted_requires_biophysical()
    test_select_by_priority()
    test_select_priority_with_weights()
    test_select_chromosomes()
    
    # Integration test
    print("\n" + "-" * 40)
    print("Integration Tests")
    print("-" * 40)
    
    test_run_select_full_pipeline()
    test_coverage_stats()
    
    print("\n" + "=" * 60)
    print("All tests passed! ✓")
    print("=" * 60)


if __name__ == '__main__':
    run_all_tests()
