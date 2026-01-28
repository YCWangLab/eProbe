"""
Dimer (inter-probe complementarity) score calculator.

Dimer formation occurs when two probes hybridize to each other
instead of their target sequences. This is problematic in multiplex
capture experiments where many probes are pooled together.

The score is based on k-mer sharing between probes:
- Higher scores indicate greater dimer formation potential
- Probes with high shared k-mer frequency should be avoided in the same pool

Design consideration:
    The current approach uses k-mer frequency analysis which is:
    1. Computationally efficient for large probe sets
    2. Provides relative ranking for probe selection
    3. Does not calculate absolute thermodynamic stability
    
    For precise dimer ΔG calculations, tools like Primer3 or
    OligoAnalyzer could be integrated in the future.

Alternative approaches:
    1. Full pairwise alignment: O(n²) complexity, expensive for large sets
    2. Thermodynamic ΔG calculation: Most accurate but slow
    3. K-mer frequency (current): Good balance of speed and accuracy
"""

from __future__ import annotations
from collections import Counter, OrderedDict
from typing import Sequence
from dataclasses import dataclass, field

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import FastaDict


@dataclass
class DimerCalculator:
    """
    Calculator for dimer scores across a probe set.
    
    Pre-computes k-mer frequencies from all probes, then calculates
    individual probe scores based on how often their k-mers appear
    in the overall pool.
    
    Usage:
        >>> calc = DimerCalculator(k=11)
        >>> result = calc.build_kmer_index(probes_dict)
        >>> if result.is_ok():
        ...     scores = calc.calculate_scores()
    """
    
    k: int = 11  # K-mer size
    min_freq: int = 2  # Minimum k-mer frequency to consider
    
    _kmer_freq: Counter = field(default_factory=Counter)
    _probes: FastaDict = field(default_factory=OrderedDict)
    _total_probes: int = 0
    
    def __post_init__(self) -> None:
        """Initialize with empty state."""
        self._kmer_freq = Counter()
        self._probes = OrderedDict()
        self._total_probes = 0
    
    def build_kmer_index(self, probes: FastaDict) -> Result[int, str]:
        """
        Build k-mer frequency index from probe sequences.
        
        Also includes reverse complements to detect sense/antisense dimers.
        
        Args:
            probes: Dictionary mapping probe IDs to sequences
            
        Returns:
            Ok(total_kmers) on success, Err(message) on failure
        """
        if not probes:
            return Err("Empty probe dictionary provided")
        
        self._probes = OrderedDict(probes)
        self._total_probes = len(probes)
        self._kmer_freq = Counter()
        
        try:
            for seq in probes.values():
                seq = seq.upper()
                
                # Count k-mers in forward sequence
                for i in range(len(seq) - self.k + 1):
                    kmer = seq[i:i + self.k]
                    self._kmer_freq[kmer] += 1
                
                # Count k-mers in reverse complement
                rc_seq = self._reverse_complement(seq)
                for i in range(len(rc_seq) - self.k + 1):
                    kmer = rc_seq[i:i + self.k]
                    self._kmer_freq[kmer] += 1
            
            # Filter low-frequency k-mers
            self._kmer_freq = Counter({
                k: v for k, v in self._kmer_freq.items() 
                if v >= self.min_freq
            })
            
            return Ok(len(self._kmer_freq))
        
        except Exception as e:
            return Err(f"Failed to build k-mer index: {e}")
    
    def calculate_score(self, sequence: str) -> Result[float, str]:
        """
        Calculate dimer score for a single sequence.
        
        Score represents the average k-mer frequency normalized
        by the total number of probes. Higher scores indicate
        more shared sequence content with other probes.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Ok(dimer_score) on success, Err(message) on failure
        """
        if self._total_probes == 0:
            return Err("No probes indexed. Call build_kmer_index first.")
        
        sequence = sequence.upper()
        
        if len(sequence) < self.k:
            return Err(f"Sequence shorter than k-mer size ({self.k})")
        
        # Sum k-mer frequencies for this sequence
        total_freq = 0
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i + self.k]
            total_freq += self._kmer_freq.get(kmer, 0)
        
        # Normalize by probe count and scale to percentage
        score = (total_freq / self._total_probes) * 100
        
        return Ok(round(score, 4))
    
    def calculate_all_scores(self) -> Result[dict[str, float], str]:
        """
        Calculate dimer scores for all indexed probes.
        
        Returns:
            Ok(dict mapping probe_id to score) on success
            Err(message) on failure
        """
        if not self._probes:
            return Err("No probes indexed")
        
        scores: dict[str, float] = {}
        
        for probe_id, sequence in self._probes.items():
            result = self.calculate_score(sequence)
            if result.is_err():
                return Err(f"Failed to score probe {probe_id}: {result.unwrap_err()}")
            scores[probe_id] = result.unwrap()
        
        return Ok(scores)
    
    @staticmethod
    def _reverse_complement(sequence: str) -> str:
        """Get reverse complement of sequence."""
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        return "".join(complement.get(base, "N") for base in reversed(sequence))


def calculate_dimer_score(
    sequence: str,
    all_sequences: Sequence[str],
    k: int = 11,
    min_freq: int = 2,
) -> Result[float, str]:
    """
    Calculate dimer score for a sequence against a pool.
    
    Convenience function for single-sequence scoring.
    For multiple sequences, use DimerCalculator class for efficiency.
    
    Args:
        sequence: Target sequence to score
        all_sequences: Pool of sequences to check against
        k: K-mer size
        min_freq: Minimum k-mer frequency threshold
        
    Returns:
        Ok(dimer_score) on success, Err(message) on failure
    """
    # Build temporary probe dictionary
    probes: FastaDict = OrderedDict()
    for i, seq in enumerate(all_sequences):
        probes[f"probe_{i}"] = seq
    
    calc = DimerCalculator(k=k, min_freq=min_freq)
    result = calc.build_kmer_index(probes)
    
    if result.is_err():
        return Err(result.unwrap_err())
    
    return calc.calculate_score(sequence)
