"""
Tests for eprobe.funcgen modules.
"""

import pytest
from pathlib import Path
from eprobe.funcgen.from_fasta import (
    tile_sequence,
    TilingConfig,
    parse_haplotype_id,
    group_haplotypes,
    HaplotypeConfig,
)
from eprobe.funcgen.from_bed import (
    parse_bed_file,
    BedRegion,
)


class TestTileSequence:
    """Tests for sequence tiling."""
    
    def test_basic_tiling(self):
        """Test basic tiling with step size."""
        config = TilingConfig(probe_length=10, step_size=5)
        sequence = "A" * 30
        
        probes = tile_sequence("test", sequence, config)
        
        assert len(probes) > 0
        assert all(len(p.sequence) == 10 for p in probes)
    
    def test_short_sequence(self):
        """Test sequence shorter than probe length."""
        config = TilingConfig(probe_length=20, step_size=10)
        sequence = "ATCGATCGATCG"  # 12 bp
        
        probes = tile_sequence("test", sequence, config)
        
        # Should return empty or single truncated probe
        assert len(probes) == 0
    
    def test_exact_probe_length(self):
        """Test sequence exactly probe length."""
        config = TilingConfig(probe_length=10, step_size=5)
        sequence = "ATCGATCGAT"  # Exactly 10 bp
        
        probes = tile_sequence("test", sequence, config)
        
        assert len(probes) == 1
        assert probes[0].sequence == "ATCGATCGAT"
    
    def test_probe_ids_sequential(self):
        """Test that probe IDs are sequential."""
        config = TilingConfig(probe_length=10, step_size=5)
        sequence = "A" * 50
        
        probes = tile_sequence("gene1", sequence, config)
        
        for i, probe in enumerate(probes):
            assert f"P{i+1:04d}" in probe.id


class TestHaplotypeProcessing:
    """Tests for haplotype ID parsing."""
    
    def test_parse_with_allele(self):
        """Test parsing ID with allele suffix."""
        gene_id, allele_id = parse_haplotype_id("BRCA1_1", "_")
        assert gene_id == "BRCA1"
        assert allele_id == "1"
    
    def test_parse_without_allele(self):
        """Test parsing ID without allele suffix."""
        gene_id, allele_id = parse_haplotype_id("BRCA1", "_")
        assert gene_id == "BRCA1"
        assert allele_id is None
    
    def test_parse_multiple_separators(self):
        """Test ID with multiple separators."""
        gene_id, allele_id = parse_haplotype_id("Gene_Name_1", "_")
        assert gene_id == "Gene_Name"
        assert allele_id == "1"
    
    def test_group_haplotypes(self):
        """Test grouping sequences by gene."""
        sequences = {
            "GENE1_1": "ATCG",
            "GENE1_2": "ATCG",
            "GENE2_1": "GCTA",
        }
        
        config = HaplotypeConfig(enabled=True, separator="_")
        groups = group_haplotypes(sequences, config)
        
        assert "GENE1" in groups
        assert "GENE2" in groups
        assert len(groups["GENE1"]) == 2
        assert len(groups["GENE2"]) == 1


class TestBedParsing:
    """Tests for BED file parsing."""
    
    def test_parse_valid_bed(self, temp_dir):
        """Test parsing valid BED file."""
        bed_path = temp_dir / "test.bed"
        bed_path.write_text(
            "chr1\t100\t200\tregion1\n"
            "chr1\t300\t400\tregion2\n"
            "chr2\t100\t500\tregion3\n"
        )
        
        result = parse_bed_file(bed_path)
        
        assert result.is_ok()
        regions = result.unwrap()
        assert len(regions) == 3
    
    def test_parse_minimal_bed(self, temp_dir):
        """Test parsing minimal 3-column BED."""
        bed_path = temp_dir / "test.bed"
        bed_path.write_text(
            "chr1\t100\t200\n"
            "chr1\t300\t400\n"
        )
        
        result = parse_bed_file(bed_path)
        
        assert result.is_ok()
        regions = result.unwrap()
        assert len(regions) == 2
    
    def test_bed_region_length(self):
        """Test BedRegion length property."""
        region = BedRegion("chr1", 100, 200)
        assert region.length == 100
    
    def test_skip_comments(self, temp_dir):
        """Test that comments are skipped."""
        bed_path = temp_dir / "test.bed"
        bed_path.write_text(
            "# Header comment\n"
            "chr1\t100\t200\n"
            "# Another comment\n"
            "chr1\t300\t400\n"
        )
        
        result = parse_bed_file(bed_path)
        
        assert result.is_ok()
        regions = result.unwrap()
        assert len(regions) == 2
