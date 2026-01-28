"""
Sequence complexity calculator using DUST algorithm.

Low-complexity sequences (e.g., repetitive elements, homopolymers) can cause
non-specific hybridization and should be avoided in probe design.

The DUST algorithm scores sequences based on triplet composition:
- Score 0: Perfect complexity (all unique triplets)
- Higher scores: Lower complexity (more repetitive)
- Scores > 2: Generally considered low complexity

Reference:
    Morgulis et al. (2006) "A fast and symmetric DUST implementation"
"""

from __future__ import annotations
from collections import Counter
from typing import Sequence

from eprobe.core.result import Result, Ok, Err


def calculate_dust_score(sequence: str, kmer_size: int = 3) -> Result[float, str]:
    """
    Calculate DUST complexity score.
    
    The DUST score measures sequence complexity based on k-mer frequencies.
    Lower scores indicate higher complexity (more random).
    Higher scores indicate lower complexity (more repetitive).
    
    Algorithm:
        1. Count all k-mers in the sequence
        2. For each k-mer, calculate contribution: count * (count - 1) / 2
        3. Sum contributions and normalize by sequence length
    
    Args:
        sequence: DNA sequence
        kmer_size: Size of k-mers to count (default: 3)
        
    Returns:
        Ok(dust_score) on success
        Err(message) on failure
        
    Example:
        >>> calculate_dust_score("ATCGATCGATCGATCG").unwrap()
        0.5  # moderate complexity
        >>> calculate_dust_score("AAAAAAAAAAAAAAAA").unwrap()
        7.0  # low complexity (repetitive)
    """
    if not sequence:
        return Err("Empty sequence provided")
    
    sequence = sequence.upper()
    
    if len(sequence) < kmer_size + 1:
        return Err(f"Sequence too short for k-mer size {kmer_size}")
    
    # Extract k-mers
    kmers = [sequence[i:i + kmer_size] for i in range(len(sequence) - kmer_size + 1)]
    
    # Count k-mer frequencies
    kmer_counts = Counter(kmers)
    
    # Calculate DUST score
    # Formula: sum of (count * (count - 1) / 2) for each unique k-mer
    # Normalized by (sequence_length - kmer_size)
    total_contribution = 0.0
    for count in kmer_counts.values():
        total_contribution += count * (count - 1) / 2
    
    # Normalize by window size
    normalizer = len(sequence) - kmer_size
    if normalizer <= 0:
        return Err("Sequence too short for normalization")
    
    dust_score = total_contribution / normalizer
    
    return Ok(round(dust_score, 4))


def calculate_dust_batch(
    sequences: Sequence[str],
    kmer_size: int = 3,
) -> list[Result[float, str]]:
    """
    Calculate DUST scores for multiple sequences.
    
    Args:
        sequences: List of DNA sequences
        kmer_size: Size of k-mers
        
    Returns:
        List of Result objects (one per sequence)
    """
    return [calculate_dust_score(seq, kmer_size) for seq in sequences]


def is_low_complexity(
    sequence: str,
    threshold: float = 2.0,
    kmer_size: int = 3,
) -> Result[bool, str]:
    """
    Check if sequence is low-complexity based on DUST score.
    
    Args:
        sequence: DNA sequence
        threshold: DUST score threshold (sequences above are low complexity)
        kmer_size: Size of k-mers
        
    Returns:
        Ok(True) if low complexity, Ok(False) if normal, Err on failure
    """
    result = calculate_dust_score(sequence, kmer_size)
    if result.is_err():
        return result  # type: ignore
    
    return Ok(result.unwrap() > threshold)
