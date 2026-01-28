"""
Hairpin (self-complementarity) score calculator.

Hairpin structures form when a probe folds back on itself due to
internal complementary regions. This reduces probe availability
for target hybridization.

The score is calculated using local alignment of the sequence
against its reverse complement. Higher scores indicate greater
potential for hairpin formation.

Design consideration:
    The current approach uses relative scores based on alignment.
    For absolute thermodynamic hairpin stability, consider using
    tools like UNAFold or NUPACK. However, relative scores are
    sufficient for ranking and filtering probes.

Alternative approaches considered:
    1. Delta G calculation: More accurate but computationally expensive
    2. Simple palindrome detection: Fast but misses bulged hairpins
    3. Current approach: Balanced accuracy and speed using alignment
"""

from __future__ import annotations
from typing import Sequence, Optional

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from eprobe.core.result import Result, Ok, Err


def calculate_hairpin_score(
    sequence: str,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_open: int = -5,
    gap_extend: int = -5,
) -> Result[float, str]:
    """
    Calculate hairpin formation potential score.
    
    Uses local alignment of sequence against its reverse complement
    to detect potential hairpin-forming regions.
    
    Scoring rationale:
        - Higher scores indicate more stable hairpin structures
        - Score components: matches (+3), mismatches (-2), gaps (-5)
        - The final score reflects both strength and extent of self-complementarity
    
    Args:
        sequence: DNA sequence (A/T/C/G)
        match_score: Score for matching bases (default: 1)
        mismatch_score: Penalty for mismatches (default: -1)
        gap_open: Gap opening penalty (default: -5)
        gap_extend: Gap extension penalty (default: -5)
        
    Returns:
        Ok(hairpin_score) on success
        Err(message) on failure
        
    Note:
        The score is relative and useful for comparing/ranking probes.
        It does not directly correspond to thermodynamic stability (ΔG).
        For capture probe design, lower scores are preferred.
    """
    if not sequence:
        return Err("Empty sequence provided")
    
    sequence = sequence.upper()
    
    # Validate sequence
    valid_bases = set("ATCGN")
    invalid = set(sequence) - valid_bases
    if invalid:
        return Err(f"Invalid bases in sequence: {invalid}")
    
    try:
        # Create aligner with local alignment mode
        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        aligner.open_gap_score = gap_open
        aligner.extend_gap_score = gap_extend
        
        # Get reverse complement
        reverse_comp = str(Seq(sequence).reverse_complement())
        
        # Perform alignment
        alignments = list(aligner.align(sequence, reverse_comp))
        
        if not alignments:
            return Ok(0.0)
        
        # Get best alignment
        best_alignment = max(alignments, key=lambda x: x.score)
        
        # Parse alignment to get detailed scoring
        # Extract alignment string representation
        alignment_str = str(best_alignment).split("\n")
        
        if len(alignment_str) >= 2:
            # Second line contains match indicators
            match_line = alignment_str[1].strip()
            
            matches = match_line.count("|")
            mismatches = match_line.count(".")
            gaps = match_line.count("-")
            
            # Calculate weighted score
            # Rationale: Matches are most important (x3), mismatches
            # destabilize (x2), gaps severely destabilize (x5)
            score = 3 * matches - 2 * mismatches - 5 * gaps
        else:
            # Fallback to raw alignment score
            score = best_alignment.score
        
        return Ok(float(score))
    
    except Exception as e:
        return Err(f"Hairpin calculation failed: {e}")


def calculate_hairpin_batch(
    sequences: Sequence[str],
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_open: int = -5,
    gap_extend: int = -5,
) -> list[Result[float, str]]:
    """
    Calculate hairpin scores for multiple sequences.
    
    Args:
        sequences: List of DNA sequences
        match_score: Score for matching bases
        mismatch_score: Penalty for mismatches
        gap_open: Gap opening penalty
        gap_extend: Gap extension penalty
        
    Returns:
        List of Result objects (one per sequence)
    """
    return [
        calculate_hairpin_score(seq, match_score, mismatch_score, gap_open, gap_extend)
        for seq in sequences
    ]


def has_strong_hairpin(
    sequence: str,
    threshold: float = 30.0,
) -> Result[bool, str]:
    """
    Check if sequence has strong hairpin potential.
    
    Args:
        sequence: DNA sequence
        threshold: Score threshold above which hairpin is considered strong
        
    Returns:
        Ok(True) if strong hairpin potential, Ok(False) if acceptable, Err on failure
    """
    result = calculate_hairpin_score(sequence)
    if result.is_err():
        return result  # type: ignore
    
    return Ok(result.unwrap() > threshold)
