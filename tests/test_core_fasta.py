"""
Tests for eprobe.core.fasta module.
"""

import pytest
from pathlib import Path
from eprobe.core.fasta import (
    read_fasta,
    write_fasta,
    extract_sequence,
    reverse_complement,
    get_chromosome_sizes,
)


class TestReadFasta:
    """Tests for read_fasta function."""
    
    def test_read_valid_fasta(self, sample_fasta):
        """Test reading a valid FASTA file."""
        result = read_fasta(sample_fasta)
        assert result.is_ok()
        
        sequences = result.unwrap()
        assert "chr1" in sequences
        assert "chr2" in sequences
        assert len(sequences["chr1"]) > 0
    
    def test_read_nonexistent_file(self, temp_dir):
        """Test reading nonexistent file returns Err."""
        result = read_fasta(temp_dir / "nonexistent.fa")
        assert result.is_err()
    
    def test_sequences_uppercase(self, sample_fasta):
        """Test that sequences are returned in uppercase."""
        result = read_fasta(sample_fasta)
        sequences = result.unwrap()
        
        for seq in sequences.values():
            assert seq == seq.upper()


class TestWriteFasta:
    """Tests for write_fasta function."""
    
    def test_write_and_read_back(self, temp_dir):
        """Test writing FASTA and reading it back."""
        sequences = {
            "seq1": "ATCGATCG",
            "seq2": "GCTAGCTA",
        }
        
        output_path = temp_dir / "output.fa"
        write_result = write_fasta(sequences, output_path)
        assert write_result.is_ok()
        
        read_result = read_fasta(output_path)
        assert read_result.is_ok()
        
        read_sequences = read_result.unwrap()
        assert read_sequences == sequences
    
    def test_write_creates_parent_dirs(self, temp_dir):
        """Test that write_fasta creates parent directories."""
        sequences = {"seq1": "ATCG"}
        output_path = temp_dir / "subdir" / "nested" / "output.fa"
        
        result = write_fasta(sequences, output_path)
        assert result.is_ok()
        assert output_path.exists()


class TestReverseComplement:
    """Tests for reverse_complement function."""
    
    def test_simple_reverse_complement(self):
        """Test basic reverse complement."""
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("GCGC") == "GCGC"
    
    def test_lowercase_input(self):
        """Test that lowercase input is handled."""
        assert reverse_complement("atcg") == "CGAT"
    
    def test_empty_string(self):
        """Test empty string."""
        assert reverse_complement("") == ""
    
    def test_palindrome(self):
        """Test palindromic sequence."""
        assert reverse_complement("ATAT") == "ATAT"


class TestExtractSequence:
    """Tests for extract_sequence function."""
    
    def test_extract_middle(self):
        """Test extracting from middle of sequence."""
        sequences = {"chr1": "ATCGATCGATCG"}
        result = extract_sequence(sequences, "chr1", 4, 9)
        assert result.is_ok()
        assert result.unwrap() == "GATCGA"
    
    def test_extract_start(self):
        """Test extracting from start."""
        sequences = {"chr1": "ATCGATCG"}
        result = extract_sequence(sequences, "chr1", 1, 4)
        assert result.is_ok()
        assert result.unwrap() == "ATCG"
    
    def test_invalid_chromosome(self):
        """Test extracting from nonexistent chromosome."""
        sequences = {"chr1": "ATCG"}
        result = extract_sequence(sequences, "chr99", 1, 4)
        assert result.is_err()
    
    def test_out_of_bounds(self):
        """Test extracting beyond sequence end."""
        sequences = {"chr1": "ATCG"}
        result = extract_sequence(sequences, "chr1", 1, 100)
        assert result.is_err()


class TestGetChromosomeSizes:
    """Tests for get_chromosome_sizes function."""
    
    def test_get_sizes(self, sample_fasta):
        """Test getting chromosome sizes."""
        read_result = read_fasta(sample_fasta)
        sequences = read_result.unwrap()
        
        sizes = get_chromosome_sizes(sequences)
        
        assert "chr1" in sizes
        assert "chr2" in sizes
        assert sizes["chr1"] > 0
        assert sizes["chr2"] > 0
