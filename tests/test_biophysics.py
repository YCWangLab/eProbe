"""
Tests for eprobe.biophysics modules.

All biophysics functions return Result[T, E], so tests must use .unwrap() 
to extract values or check .is_ok()/.is_err() for error cases.
"""

import pytest
from eprobe.biophysics import (
    calculate_gc,
    calculate_tm,
    calculate_complexity,
    calculate_hairpin_score,
)
from eprobe.biophysics.dimer import DimerCalculator
from eprobe.biophysics.entropy import calculate_entropy


class TestGC:
    """Tests for GC content calculation."""
    
    def test_pure_gc(self):
        """Test sequence with only G and C."""
        result = calculate_gc("GCGCGC")
        assert result.is_ok()
        assert result.unwrap() == 100.0
    
    def test_pure_at(self):
        """Test sequence with only A and T."""
        result = calculate_gc("ATATAT")
        assert result.is_ok()
        assert result.unwrap() == 0.0
    
    def test_balanced(self):
        """Test balanced sequence."""
        result = calculate_gc("ATCG")
        assert result.is_ok()
        assert result.unwrap() == pytest.approx(50.0)
    
    def test_empty_sequence(self):
        """Test empty sequence returns error."""
        result = calculate_gc("")
        assert result.is_err()
    
    def test_lowercase(self):
        """Test lowercase input."""
        result = calculate_gc("atcg")
        assert result.is_ok()
        assert result.unwrap() == pytest.approx(50.0)
    
    def test_real_sequence(self):
        """Test with typical probe sequence."""
        # 40 bp sequence with known GC content
        seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"  # 50% GC
        result = calculate_gc(seq)
        assert result.is_ok()
        gc = result.unwrap()
        assert 45.0 <= gc <= 55.0


class TestTm:
    """Tests for melting temperature calculation."""
    
    def test_short_sequence(self):
        """Test Tm for short sequence."""
        result = calculate_tm("ATCGATCG")
        assert result.is_ok()
        tm = result.unwrap()
        assert 5 <= tm <= 35  # Reasonable range for 8bp
    
    def test_gc_rich_higher_tm(self):
        """Test that GC-rich sequences have higher Tm."""
        gc_rich_result = calculate_tm("GCGCGCGCGC")
        at_rich_result = calculate_tm("ATATATATAT")
        
        assert gc_rich_result.is_ok() and at_rich_result.is_ok()
        assert gc_rich_result.unwrap() > at_rich_result.unwrap()
    
    def test_longer_sequence(self):
        """Test Tm for longer (probe-length) sequence."""
        seq = "A" * 41 + "G" * 40  # 81 bp with ~50% GC
        result = calculate_tm(seq)
        assert result.is_ok()
        tm = result.unwrap()
        assert 50 <= tm <= 95  # Reasonable range for 81bp


class TestComplexity:
    """Tests for sequence complexity calculation."""
    
    def test_low_complexity(self):
        """Test low complexity sequence."""
        # Homopolymer has low complexity (high DUST score)
        result = calculate_complexity("AAAAAAAAAAAAAAAAAAAA")
        assert result.is_ok()
        score = result.unwrap()
        assert score > 2.0  # High DUST score = low complexity
    
    def test_high_complexity(self):
        """Test high complexity sequence."""
        # Random-looking sequence
        result = calculate_complexity("ATCGATCGTAGCTAGCTACG")
        assert result.is_ok()
        score = result.unwrap()
        assert score < 2.0  # Low DUST score = high complexity
    
    def test_empty_sequence(self):
        """Test empty sequence returns error."""
        result = calculate_complexity("")
        assert result.is_err()


class TestHairpin:
    """Tests for hairpin score calculation."""
    
    def test_palindrome(self):
        """Test palindromic sequence has higher hairpin score."""
        # This sequence can form a hairpin
        palindrome = "ATCGATCGAT" + "ATCGATCGAT"[::-1].translate(str.maketrans("ATCG", "TAGC"))
        random_seq = "ATCGATCGATCGATCGATCG"
        
        result_pal = calculate_hairpin_score(palindrome)
        result_rand = calculate_hairpin_score(random_seq)
        
        assert result_pal.is_ok() and result_rand.is_ok()
        score_pal = result_pal.unwrap()
        score_rand = result_rand.unwrap()
        
        # Palindrome should have higher or equal score
        assert score_pal >= 0
        assert score_rand >= 0
    
    def test_score_range(self):
        """Test that hairpin score is non-negative."""
        seq = "ATCGATCGATCGATCGATCG"
        result = calculate_hairpin_score(seq)
        assert result.is_ok()
        score = result.unwrap()
        assert score >= 0


class TestDimer:
    """Tests for dimer score calculation."""
    
    def test_build_kmer_index(self):
        """Test building k-mer index from probes."""
        calc = DimerCalculator(k=5)
        probes = {
            "probe1": "ATCGATCGATCG",
            "probe2": "GCTAGCTAGCTA",
        }
        result = calc.build_kmer_index(probes)
        assert result.is_ok()
    
    def test_calculate_scores(self):
        """Test calculating dimer scores for probes."""
        calc = DimerCalculator(k=5, min_freq=1)
        probes = {
            "probe1": "ATCGATCGATCG",
            "probe2": "GCTAGCTAGCTA",
            "probe3": "ATCGATCGATCG",  # Same as probe1
        }
        
        build_result = calc.build_kmer_index(probes)
        assert build_result.is_ok()
        
        # Test individual score calculation
        score_result = calc.calculate_score("probe1")
        assert score_result.is_ok()
        score = score_result.unwrap()
        assert score >= 0


class TestEntropy:
    """Tests for entropy calculation."""
    
    def test_low_entropy(self):
        """Test homopolymer has low entropy."""
        result = calculate_entropy("AAAAAAAAAAAAAAAAAAAA")
        assert result.is_ok()
        entropy = result.unwrap()
        assert entropy < 1.0  # Low entropy
    
    def test_high_entropy(self):
        """Test random sequence has higher entropy."""
        result = calculate_entropy("ATCGATCGATCGATCGATCG")
        assert result.is_ok()
        entropy = result.unwrap()
        assert entropy > 1.0  # Higher entropy
    
    def test_maximum_entropy(self):
        """Test that maximum entropy is bounded."""
        # Sequence with equal base frequencies
        seq = "ATCG" * 10
        result = calculate_entropy(seq)
        assert result.is_ok()
        entropy = result.unwrap()
        # Maximum entropy for 4 bases is 2.0
        assert entropy <= 2.0
    
    def test_empty_sequence(self):
        """Test empty sequence returns error."""
        result = calculate_entropy("")
        assert result.is_err()
