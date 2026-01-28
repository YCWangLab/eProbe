"""
Tests for eprobe.util modules.
"""

import pytest
from pathlib import Path
from eprobe.util.tiling import parse_positions, TilePosition
from eprobe.util.dedup import deduplicate_exact
from eprobe.util.validate import (
    validate_snp_dataframe,
    validate_probe_fasta,
)


class TestPositionParsing:
    """Tests for tiling position parsing."""
    
    def test_named_positions(self):
        """Test parsing named positions."""
        result = parse_positions("center,left,right", default_offset=20)
        
        assert result.is_ok()
        positions = result.unwrap()
        
        assert len(positions) == 3
        assert positions[0].name == "C"
        assert positions[0].offset == 0
        assert positions[1].name == "L20"
        assert positions[1].offset == -20
        assert positions[2].name == "R20"
        assert positions[2].offset == 20
    
    def test_numeric_positions(self):
        """Test parsing numeric positions."""
        result = parse_positions("-30,0,+30", default_offset=20)
        
        assert result.is_ok()
        positions = result.unwrap()
        
        assert len(positions) == 3
        assert positions[0].offset == -30
        assert positions[1].offset == 0
        assert positions[2].offset == 30
    
    def test_mixed_positions(self):
        """Test parsing mixed named and numeric."""
        result = parse_positions("center,-20,+20", default_offset=25)
        
        assert result.is_ok()
        positions = result.unwrap()
        assert len(positions) == 3
    
    def test_invalid_position(self):
        """Test invalid position returns error."""
        result = parse_positions("invalid_position")
        assert result.is_err()


class TestDeduplication:
    """Tests for sequence deduplication."""
    
    def test_no_duplicates(self):
        """Test when no duplicates exist."""
        sequences = {
            "seq1": "ATCGATCG",
            "seq2": "GCTAGCTA",
            "seq3": "TGCATGCA",
        }
        
        result = deduplicate_exact(sequences)
        
        assert len(result) == 3
    
    def test_exact_duplicates(self):
        """Test removing exact duplicates."""
        sequences = {
            "seq1": "ATCGATCG",
            "seq2": "ATCGATCG",  # Duplicate
            "seq3": "GCTAGCTA",
        }
        
        result = deduplicate_exact(sequences, keep="first")
        
        assert len(result) == 2
        assert "ATCGATCG" in result.values()
    
    def test_case_insensitive(self):
        """Test that deduplication is case insensitive."""
        sequences = {
            "seq1": "ATCGATCG",
            "seq2": "atcgatcg",  # Lowercase duplicate
        }
        
        result = deduplicate_exact(sequences)
        
        assert len(result) == 1


class TestValidation:
    """Tests for file validation."""
    
    def test_valid_snp_dataframe(self, temp_dir):
        """Test validating correct SNP TSV."""
        tsv_path = temp_dir / "valid.tsv"
        tsv_path.write_text(
            "chrom\tpos\tref\talt\n"
            "chr1\t100\tA\tG\n"
            "chr1\t200\tC\tT\n"
        )
        
        result = validate_snp_dataframe(tsv_path)
        
        assert result.is_ok()
        validation = result.unwrap()
        assert validation["valid"] is True
    
    def test_missing_columns(self, temp_dir):
        """Test validation fails with missing columns."""
        tsv_path = temp_dir / "invalid.tsv"
        tsv_path.write_text(
            "chrom\tpos\n"  # Missing ref and alt
            "chr1\t100\n"
        )
        
        result = validate_snp_dataframe(tsv_path)
        
        assert result.is_ok()  # Function returns Ok with valid=False
        validation = result.unwrap()
        assert validation["valid"] is False
    
    def test_valid_fasta(self, sample_fasta):
        """Test validating correct FASTA."""
        result = validate_probe_fasta(sample_fasta)
        
        assert result.is_ok()
        validation = result.unwrap()
        assert validation["valid"] is True


class TestUtilSubset:
    """Tests for subset utility (edge cases)."""
    
    def test_sample_exceeds_input(self, sample_sequences):
        """Test sampling more than available."""
        from eprobe.util.subset import sample_random
        
        result = sample_random(sample_sequences, n=100)
        
        # Should return all sequences
        assert len(result) == len(sample_sequences)
    
    def test_sample_exact_count(self, sample_sequences):
        """Test sampling exact count."""
        from eprobe.util.subset import sample_random
        
        result = sample_random(sample_sequences, n=2, seed=42)
        
        assert len(result) == 2
    
    def test_reproducible_sampling(self, sample_sequences):
        """Test that sampling is reproducible with seed."""
        from eprobe.util.subset import sample_random
        
        result1 = sample_random(sample_sequences, n=2, seed=42)
        result2 = sample_random(sample_sequences, n=2, seed=42)
        
        assert set(result1.keys()) == set(result2.keys())
