"""
Tests for eprobe.popgen modules.
"""

import pytest
from pathlib import Path
from eprobe.popgen.extract import (
    detect_clusters,
    apply_cluster_filter,
    ClusterConfig,
)
from eprobe.popgen.select import (
    select_uniform,
    select_random,
    calculate_window_index,
)
from eprobe.popgen.build import (
    calculate_probe_coordinates,
    ProbeConfig,
)
from eprobe.core.models import SNP


# Helper to create SNPs easily
def make_snp(chrom: str, pos: int, ref: str = "A", alt: str = "G") -> SNP:
    """Create SNP using from_vcf_record (auto-determines mutation_type)."""
    return SNP.from_vcf_record(chrom, pos, ref, alt)


class TestClusterDetection:
    """Tests for SNP cluster detection."""
    
    def test_no_clusters(self):
        """Test when no clusters exist."""
        snps = [
            make_snp("chr1", 100, "A", "G"),
            make_snp("chr1", 500, "C", "T"),
            make_snp("chr1", 1000, "G", "A"),
        ]
        
        config = ClusterConfig(flank=60, max_snp=3)
        result = apply_cluster_filter(snps, config)
        
        assert result.is_ok()
        assert len(result.unwrap()) == 3  # No SNPs removed
    
    def test_cluster_detected(self):
        """Test cluster is detected and removed."""
        # Create a cluster of 5 SNPs within 60bp window
        snps = [
            make_snp("chr1", 100, "A", "G"),
            make_snp("chr1", 110, "C", "T"),
            make_snp("chr1", 120, "G", "A"),
            make_snp("chr1", 130, "T", "C"),
            make_snp("chr1", 140, "A", "T"),
            make_snp("chr1", 1000, "G", "C"),  # Far away, not in cluster
        ]
        
        config = ClusterConfig(flank=60, max_snp=3)
        result = apply_cluster_filter(snps, config)
        
        assert result.is_ok()
        # Clustered SNPs should be removed
        assert len(result.unwrap()) < 6
    
    def test_disabled_filter(self):
        """Test cluster filter when disabled."""
        snps = [
            make_snp("chr1", 100, "A", "G"),
            make_snp("chr1", 110, "C", "T"),  # Would be in cluster
        ]
        
        config = ClusterConfig(flank=60, max_snp=1, enabled=False)
        result = apply_cluster_filter(snps, config)
        
        assert result.is_ok()
        assert len(result.unwrap()) == 2  # None removed


class TestWindowCalculation:
    """Tests for window-based calculations."""
    
    def test_window_index(self):
        """Test window index calculation."""
        assert calculate_window_index(50000, 100000) == 0
        assert calculate_window_index(150000, 100000) == 1
        assert calculate_window_index(250000, 100000) == 2
    
    def test_window_boundary(self):
        """Test positions at window boundaries."""
        assert calculate_window_index(100000, 100000) == 1
        assert calculate_window_index(99999, 100000) == 0


class TestSNPSelection:
    """Tests for SNP selection strategies."""
    
    def test_uniform_selection(self):
        """Test uniform distribution selection."""
        # Create SNPs across different windows
        snps = [
            make_snp("chr1", 50000 + i * 100, "A", "G")
            for i in range(10)
        ] + [
            make_snp("chr1", 150000 + i * 100, "C", "T")
            for i in range(10)
        ]
        
        result = select_uniform(snps, target_count=5, window_size=100000)
        
        assert result.is_ok()
        selected = result.unwrap()
        assert len(selected) == 5
    
    def test_random_selection(self):
        """Test random selection."""
        snps = [make_snp("chr1", i * 1000, "A", "G") for i in range(1, 101)]
        
        result = select_random(snps, target_count=10, seed=42)
        
        assert result.is_ok()
        selected = result.unwrap()
        assert len(selected) == 10
    
    def test_selection_preserves_order(self):
        """Test that selection maintains sorted order."""
        snps = [make_snp("chr1", i * 1000, "A", "G") for i in range(1, 51)]
        
        result = select_random(snps, target_count=20, seed=42)
        selected = result.unwrap()
        
        # Should be sorted by position
        positions = [s.pos for s in selected]
        assert positions == sorted(positions)
    
    def test_target_exceeds_input(self):
        """Test when target exceeds input count."""
        snps = [make_snp("chr1", 100, "A", "G"), make_snp("chr1", 200, "C", "T")]
        
        result = select_random(snps, target_count=100)
        
        assert result.is_ok()
        assert len(result.unwrap()) == 2  # Returns all SNPs


class TestProbeCoordinates:
    """Tests for probe coordinate calculation."""
    
    def test_centered_probe(self):
        """Test probe with SNP at center."""
        start, end, offset = calculate_probe_coordinates(
            snp_pos=100,
            probe_length=81,
            shift=0,
        )
        
        # SNP should be at center (position 40 in 0-based)
        assert offset == 40
        assert end - start == 81
    
    def test_shifted_probe(self):
        """Test probe with SNP shifted."""
        start, end, offset = calculate_probe_coordinates(
            snp_pos=100,
            probe_length=81,
            shift=20,  # Shift right, SNP moves left
        )
        
        # SNP should be offset from center
        assert offset == 20  # 40 - 20
    
    def test_probe_length_odd(self):
        """Test odd probe length."""
        start, end, offset = calculate_probe_coordinates(
            snp_pos=100,
            probe_length=81,
            shift=0,
        )
        
        assert end - start == 81
