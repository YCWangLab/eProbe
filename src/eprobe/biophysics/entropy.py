"""
Sequence entropy calculator.

Shannon entropy measures the information content of a sequence.
Lower entropy indicates more repetitive/predictable sequences,
while higher entropy indicates more random/complex sequences.

This is an ALTERNATIVE complexity metric to DUST score:
- Entropy: Based on information theory, considers base frequencies
- DUST: Based on triplet frequencies, more sensitive to short repeats

Both metrics can be used together for comprehensive complexity filtering.

NOTE: This module is a placeholder for future implementation.
The interface is defined but full functionality is not yet implemented.
"""

from __future__ import annotations
import math
from collections import Counter
from typing import Sequence

from eprobe.core.result import Result, Ok, Err


def calculate_entropy(sequence: str, base: int = 2) -> Result[float, str]:
    """
    Calculate Shannon entropy of a DNA sequence.
    
    Formula: H = -sum(p(x) * log(p(x))) for each base x
    
    Theoretical values:
        - Maximum entropy (random sequence): 2.0 bits (base 2)
        - Minimum entropy (single base): 0.0 bits
        
    Args:
        sequence: DNA sequence
        base: Logarithm base (2 for bits, e for nats)
        
    Returns:
        Ok(entropy_value) on success
        Err(message) on failure
        
    Example:
        >>> calculate_entropy("ATCGATCG").unwrap()
        2.0  # maximum entropy (equal base frequencies)
        >>> calculate_entropy("AAAAAAAA").unwrap()
        0.0  # minimum entropy (single base)
        
    Note:
        This is a basic implementation. Future versions may include:
        - Sliding window entropy for local complexity
        - Di/tri-nucleotide entropy
        - Conditional entropy
    """
    if not sequence:
        return Err("Empty sequence provided")
    
    sequence = sequence.upper()
    
    # Count base frequencies
    base_counts = Counter(sequence)
    total_bases = len(sequence)
    
    # Calculate entropy
    entropy = 0.0
    for count in base_counts.values():
        if count > 0:
            probability = count / total_bases
            entropy -= probability * math.log(probability, base)
    
    return Ok(round(entropy, 4))


def calculate_entropy_batch(
    sequences: Sequence[str],
    base: int = 2,
) -> list[Result[float, str]]:
    """
    Calculate entropy for multiple sequences.
    
    Args:
        sequences: List of DNA sequences
        base: Logarithm base
        
    Returns:
        List of Result objects (one per sequence)
    """
    return [calculate_entropy(seq, base) for seq in sequences]


def calculate_window_entropy(
    sequence: str,
    window_size: int = 20,
    step: int = 1,
    base: int = 2,
) -> Result[list[float], str]:
    """
    Calculate sliding window entropy across a sequence.
    
    Useful for identifying local low-complexity regions within
    an otherwise normal-complexity sequence.
    
    Args:
        sequence: DNA sequence
        window_size: Size of sliding window
        step: Step size for window movement
        base: Logarithm base
        
    Returns:
        Ok(list of entropy values) on success
        Err(message) on failure
        
    Note:
        PLACEHOLDER - Full implementation pending.
        Returns basic single-value entropy for now.
    """
    if not sequence:
        return Err("Empty sequence provided")
    
    if len(sequence) < window_size:
        return Err(f"Sequence shorter than window size ({window_size})")
    
    entropies: list[float] = []
    
    for i in range(0, len(sequence) - window_size + 1, step):
        window = sequence[i:i + window_size]
        result = calculate_entropy(window, base)
        if result.is_err():
            return Err(f"Failed at position {i}: {result.unwrap_err()}")
        entropies.append(result.unwrap())
    
    return Ok(entropies)


def is_low_entropy(
    sequence: str,
    threshold: float = 1.5,
    base: int = 2,
) -> Result[bool, str]:
    """
    Check if sequence has low entropy (potentially repetitive).
    
    Args:
        sequence: DNA sequence
        threshold: Entropy threshold below which sequence is low entropy
        base: Logarithm base
        
    Returns:
        Ok(True) if low entropy, Ok(False) if normal, Err on failure
        
    Note:
        Default threshold of 1.5 bits is conservative.
        Random DNA has ~2.0 bits, so 1.5 catches moderately biased sequences.
    """
    result = calculate_entropy(sequence, base)
    if result.is_err():
        return result  # type: ignore
    
    return Ok(result.unwrap() < threshold)
