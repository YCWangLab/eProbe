"""
GC content calculator.

GC content affects probe stability and hybridization efficiency.
Optimal GC content for capture probes is typically 35-65%.
"""

from __future__ import annotations
from typing import Sequence

from eprobe.core.result import Result, Ok, Err


def calculate_gc(sequence: str) -> Result[float, str]:
    """
    Calculate GC content as percentage.
    
    Args:
        sequence: DNA sequence (A/T/C/G/N)
        
    Returns:
        Ok(gc_percentage) on success (0-100 scale)
        Err(message) on failure
        
    Example:
        >>> calculate_gc("ATCGATCG").unwrap()
        50.0
        >>> calculate_gc("GCGCGC").unwrap()
        100.0
    """
    if not sequence:
        return Err("Empty sequence provided")
    
    sequence = sequence.upper()
    
    # Count valid bases only (exclude N and other ambiguous)
    gc_count = sequence.count("G") + sequence.count("C")
    total_valid = sum(sequence.count(base) for base in "ATCG")
    
    if total_valid == 0:
        return Err("No valid bases (A/T/C/G) found in sequence")
    
    gc_percent = (gc_count / total_valid) * 100
    
    return Ok(round(gc_percent, 4))


def calculate_gc_batch(sequences: Sequence[str]) -> list[Result[float, str]]:
    """
    Calculate GC content for multiple sequences.
    
    Args:
        sequences: List of DNA sequences
        
    Returns:
        List of Result objects (one per sequence)
    """
    return [calculate_gc(seq) for seq in sequences]


def gc_in_range(sequence: str, min_gc: float, max_gc: float) -> Result[bool, str]:
    """
    Check if GC content falls within specified range.
    
    Args:
        sequence: DNA sequence
        min_gc: Minimum GC percentage (inclusive)
        max_gc: Maximum GC percentage (inclusive)
        
    Returns:
        Ok(True) if in range, Ok(False) if not, Err on calculation failure
    """
    result = calculate_gc(sequence)
    if result.is_err():
        return result  # type: ignore
    
    gc = result.unwrap()
    return Ok(min_gc <= gc <= max_gc)
