"""
Hairpin (secondary structure) score calculator using ViennaRNA.

Hairpin structures form when a single-stranded DNA probe folds back
on itself due to internal complementary regions. This reduces probe
availability for target hybridization.

The score is calculated using ViennaRNA's Zuker MFE algorithm with
DNA parameters (Mathews 2004). The returned value is the absolute
value of the minimum free energy (|ΔG| in kcal/mol):

    - 0.0: No secondary structure predicted
    - Higher values: More stable self-folding (worse for probes)
    - Typical threshold: 5.0 kcal/mol (reject |ΔG| > 5)

ViennaRNA (C implementation) is fast and has no sequence length limit,
making it suitable for probes of any length (40–200+ bp).

Requires: pip install ViennaRNA  (or: pip install eprobe[thermo])
"""

from __future__ import annotations
from typing import Sequence

from eprobe.core.result import Result, Ok, Err

try:
    import RNA as _RNA
    # Load DNA parameters (Mathews 2004) once at import time
    _RNA.params_load_DNA_Mathews2004()
    _HAS_VIENNA = True
except ImportError:
    _HAS_VIENNA = False


def calculate_hairpin_score(sequence: str) -> Result[float, str]:
    """
    Calculate hairpin formation potential as |MFE| (kcal/mol) via ViennaRNA.

    Uses the Zuker MFE algorithm with DNA thermodynamic parameters
    (Mathews 2004). Returns the absolute value of the minimum free
    energy so that higher values indicate stronger self-folding.

    Args:
        sequence: DNA sequence (A/T/C/G)

    Returns:
        Ok(|MFE|) on success (0.0 = no structure, higher = more hairpin)
        Err(message) on failure
    """
    if not _HAS_VIENNA:
        return Err(
            "ViennaRNA is not installed. "
            "Install with: pip install ViennaRNA"
        )

    if not sequence:
        return Err("Empty sequence provided")

    sequence = sequence.upper()

    valid_bases = set("ATCGN")
    invalid = set(sequence) - valid_bases
    if invalid:
        return Err(f"Invalid bases in sequence: {invalid}")

    # ViennaRNA cannot handle N; replace with A (conservative)
    clean = sequence.replace("N", "A")

    try:
        _ss, mfe = _RNA.fold(clean)
        return Ok(round(abs(mfe), 4))
    except Exception as e:
        return Err(f"ViennaRNA fold failed: {e}")


def calculate_hairpin_batch(
    sequences: Sequence[str],
) -> list[Result[float, str]]:
    """
    Calculate hairpin scores for multiple sequences.

    Args:
        sequences: List of DNA sequences

    Returns:
        List of Result objects (one per sequence)
    """
    return [calculate_hairpin_score(seq) for seq in sequences]


def has_strong_hairpin(
    sequence: str,
    threshold: float = 5.0,
) -> Result[bool, str]:
    """
    Check if sequence has strong hairpin potential.

    Args:
        sequence: DNA sequence
        threshold: |MFE| threshold (kcal/mol). Default 5.0.

    Returns:
        Ok(True) if |MFE| > threshold, Ok(False) if acceptable
    """
    result = calculate_hairpin_score(sequence)
    if result.is_err():
        return result  # type: ignore
    return Ok(result.unwrap() > threshold)
