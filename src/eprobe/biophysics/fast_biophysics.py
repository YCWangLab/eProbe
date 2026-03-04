"""
Optimized biophysical property calculators for probe sequences.

This module provides fast implementations of biophysical calculations:
- GC content: Vectorized numpy calculation
- Tm: Nearest-neighbor with lookup table optimization
- DUST: Optimized k-mer counting
- Hairpin: Fast local alignment using parasail (if available) or optimized fallback
- Dimer: K-mer frequency sharing with pre-built index

All batch functions are optimized for processing large numbers of sequences.
"""

from __future__ import annotations
import math
from collections import Counter
from typing import List, Dict, Tuple, Optional, Sequence
from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as BioTm

from eprobe.core.result import Result, Ok, Err

# Try to import parasail for fast alignment, fall back to Bio.Align
try:
    import parasail
    HAS_PARASAIL = True
except ImportError:
    HAS_PARASAIL = False
    from Bio.Align import PairwiseAligner


# =============================================================================
# GC Content - Optimized
# =============================================================================

def calculate_gc_fast(sequence: str) -> float:
    """
    Calculate GC content as percentage (optimized version).
    
    Args:
        sequence: DNA sequence (A/T/C/G/N)
        
    Returns:
        GC percentage (0-100 scale)
        
    Raises:
        ValueError: If sequence is empty or has no valid bases
    """
    if not sequence:
        raise ValueError("Empty sequence provided")
    
    seq = sequence.upper()
    gc_count = seq.count('G') + seq.count('C')
    total_valid = seq.count('A') + seq.count('T') + seq.count('G') + seq.count('C')
    
    if total_valid == 0:
        raise ValueError("No valid bases (A/T/C/G) found")
    
    return (gc_count / total_valid) * 100


def calculate_gc_batch_fast(sequences: List[str]) -> List[float]:
    """
    Calculate GC content for multiple sequences (optimized).
    
    Args:
        sequences: List of DNA sequences
        
    Returns:
        List of GC percentages
    """
    return [calculate_gc_fast(seq) for seq in sequences]


# =============================================================================
# Melting Temperature - Using Biopython NN tables
# =============================================================================

# Available NN tables from Biopython for user selection
# Users should choose based on their experimental setup (DNA/DNA vs RNA/DNA hybrid)
NN_TABLE_OPTIONS = {
    # DNA/DNA hybridization tables
    "DNA_NN1": BioTm.DNA_NN1,  # Breslauer 1986 - older, less accurate
    "DNA_NN2": BioTm.DNA_NN2,  # Sugimoto 1996 - DNA/DNA
    "DNA_NN3": BioTm.DNA_NN3,  # Allawi 1997 - improved parameters
    "DNA_NN4": BioTm.DNA_NN4,  # SantaLucia 1998 - most widely used, recommended
    
    # RNA/DNA hybrid table (for RNA baits capturing DNA targets)
    "R_DNA_NN1": BioTm.R_DNA_NN1,  # Sugimoto 1995 - RNA probe / DNA target
}


def calculate_tm_fast(
    sequence: str, 
    na_conc: float = 50.0,
    nn_table: str = "DNA_NN4",
    dnac1: float = 250.0,
    dnac2: float = 0.0,
) -> float:
    """
    Calculate melting temperature using nearest-neighbor method.
    
    Uses Biopython's Tm_NN with user-selectable NN parameter tables.
    For short sequences (<14bp), uses Wallace rule as fallback.
    
    Args:
        sequence: DNA sequence
        na_conc: Na+ concentration in mM (default: 50)
        nn_table: NN parameter table name. Options:
            DNA/DNA hybridization:
            - "DNA_NN1": Breslauer 1986 (older parameters)
            - "DNA_NN2": Sugimoto 1996
            - "DNA_NN3": Allawi 1997
            - "DNA_NN4": SantaLucia 1998 (recommended, default)
            
            RNA/DNA hybrid (for RNA baits):
            - "R_DNA_NN1": Sugimoto 1995
            
        dnac1: Concentration of higher concentrated strand (nM, default: 250)
        dnac2: Concentration of lower concentrated strand (nM, default: 0 = excess)
        
    Returns:
        Melting temperature in °C
        
    Raises:
        ValueError: If sequence is invalid or nn_table not recognized
        
    Note:
        For capture probe design:
        - DNA probes + DNA targets: use "DNA_NN4" (default)
        - RNA probes + DNA targets: use "R_DNA_NN1"
        
        References:
        - Breslauer et al. (1986) PNAS 83:3746-3750
        - Sugimoto et al. (1995) Biochemistry 34:11211-11216
        - Sugimoto et al. (1996) Nucleic Acids Res 24:4501-4505
        - Allawi & SantaLucia (1997) Biochemistry 36:10581-10594
        - SantaLucia (1998) PNAS 95:1460-1465
    """
    if not sequence:
        raise ValueError("Empty sequence provided")
    
    seq = sequence.upper()
    
    if len(seq) < 2:
        raise ValueError("Sequence must be at least 2 bases")
    
    # Validate nn_table selection
    if nn_table not in NN_TABLE_OPTIONS:
        valid_tables = ", ".join(NN_TABLE_OPTIONS.keys())
        raise ValueError(f"Unknown nn_table '{nn_table}'. Valid options: {valid_tables}")
    
    # For short sequences, use Wallace rule (more reliable)
    if len(seq) < 14:
        gc_count = seq.count('G') + seq.count('C')
        at_count = seq.count('A') + seq.count('T')
        return 4.0 * gc_count + 2.0 * at_count
    
    # Check for invalid characters
    valid = set('ATCGN')
    if not all(c in valid for c in seq):
        # Use Biopython fallback with default table
        return float(BioTm.Tm_NN(seq))
    
    # Use Biopython's NN calculation with user-selected table
    try:
        tm = BioTm.Tm_NN(
            seq,
            nn_table=NN_TABLE_OPTIONS[nn_table],
            Na=na_conc,
            dnac1=dnac1,
            dnac2=dnac2,
        )
        return round(float(tm), 2)
    except Exception:
        # Fallback to simpler GC-based calculation (Owczarzy formula)
        gc_content = (seq.count('G') + seq.count('C')) / len(seq)
        tm = 81.5 + 16.6 * math.log10(na_conc / 1000) + 41 * gc_content - 500 / len(seq)
        return round(tm, 2)


def calculate_tm_batch_fast(
    sequences: List[str], 
    na_conc: float = 50.0,
    nn_table: str = "DNA_NN4",
) -> List[float]:
    """
    Calculate Tm for multiple sequences (optimized).
    
    Args:
        sequences: List of DNA sequences
        na_conc: Na+ concentration in mM
        nn_table: NN parameter table name (see calculate_tm_fast for options)
        
    Returns:
        List of Tm values in °C
    """
    return [calculate_tm_fast(seq, na_conc, nn_table) for seq in sequences]


# =============================================================================
# DUST Complexity - Optimized
# =============================================================================

def calculate_dust_fast(sequence: str, k: int = 3) -> float:
    """
    Calculate DUST complexity score (optimized).
    
    Lower scores = higher complexity (more random)
    Higher scores = lower complexity (more repetitive)
    
    Args:
        sequence: DNA sequence
        k: K-mer size (default: 3)
        
    Returns:
        DUST score
        
    Raises:
        ValueError: If sequence is too short
    """
    if not sequence:
        raise ValueError("Empty sequence provided")
    
    seq = sequence.upper()
    
    if len(seq) < k + 1:
        raise ValueError(f"Sequence too short for k-mer size {k}")
    
    # Count k-mers using Counter (optimized)
    kmer_counts = Counter(seq[i:i+k] for i in range(len(seq) - k + 1))
    
    # Calculate DUST score: sum(count * (count-1) / 2) / (len - k)
    total = sum(c * (c - 1) / 2 for c in kmer_counts.values())
    normalizer = len(seq) - k
    
    return round(total / normalizer, 4) if normalizer > 0 else 0.0


def calculate_dust_batch_fast(sequences: List[str], k: int = 3) -> List[float]:
    """
    Calculate DUST scores for multiple sequences.
    
    Args:
        sequences: List of DNA sequences
        k: K-mer size
        
    Returns:
        List of DUST scores
    """
    return [calculate_dust_fast(seq, k) for seq in sequences]


# =============================================================================
# Hairpin Score - Fast Local Alignment
# =============================================================================

# Complement table for fast reverse complement
_COMPLEMENT = str.maketrans('ATCGN', 'TAGCN')


def _reverse_complement(sequence: str) -> str:
    """Get reverse complement of sequence."""
    return sequence.upper().translate(_COMPLEMENT)[::-1]


def _calculate_max_stem(seq: str, rev_comp: str, min_loop: int = 3) -> int:
    """
    Find the longest consecutive complementary stem between a sequence
    and its reverse complement (i.e., the longest possible hairpin stem).
    
    At each offset k, seq[i] is aligned with rev_comp[k+i], meaning
    seq[i] pairs with complement(seq[n-1-k-i]) in the original sequence.
    The offset represents the approximate loop + stem region.
    
    Args:
        seq: Uppercase DNA sequence
        rev_comp: Reverse complement of seq
        min_loop: Minimum hairpin loop size (default: 3)
        
    Returns:
        Length of the longest consecutive complementary stem (in bp)
    """
    n = len(seq)
    if n < min_loop + 4:
        return 0
    
    best_stem = 0
    
    for offset in range(min_loop, n - 3):
        overlap_len = n - offset
        if overlap_len < 4:
            break
        
        current_run = 0
        max_run = 0
        
        for i in range(overlap_len):
            if seq[i] == rev_comp[offset + i]:
                current_run += 1
                if current_run > max_run:
                    max_run = current_run
            else:
                current_run = 0
        
        if max_run > best_stem:
            best_stem = max_run
    
    return best_stem


# =============================================================================
# Hairpin: Exponential Bonus Algorithm (k-mer + spatial continuity)
# =============================================================================

# 2-bit encoding lookup: A=00, C=01, G=10, T=11
# Complement: XOR with 0b11 → A(00)↔T(11), C(01)↔G(10)
_BASE_TO_2BIT = [4] * 256  # 4 = invalid sentinel
_BASE_TO_2BIT[ord('A')] = 0b00
_BASE_TO_2BIT[ord('a')] = 0b00
_BASE_TO_2BIT[ord('C')] = 0b01
_BASE_TO_2BIT[ord('c')] = 0b01
_BASE_TO_2BIT[ord('G')] = 0b10
_BASE_TO_2BIT[ord('g')] = 0b10
_BASE_TO_2BIT[ord('T')] = 0b11
_BASE_TO_2BIT[ord('t')] = 0b11


def _rc_hash(h: int, k: int) -> int:
    """
    Compute reverse-complement k-mer hash via bit manipulation.
    
    Each base is 2-bit encoded. Complement = XOR 0b11 per base,
    then reverse the order of 2-bit chunks.
    
    Time: O(k). Called once per k-mer position.
    """
    rc = 0
    for _ in range(k):
        rc = (rc << 2) | ((h & 0b11) ^ 0b11)
        h >>= 2
    return rc


def _calculate_hairpin_exponential(
    seq: str,
    k: int = 4,
    min_loop: int = 3,
    base: float = 4.0,
    early_stop_raw: float = None,
) -> float:
    """
    Calculate hairpin raw score using exponential bonus for consecutive
    k-mer reverse-complement matches — reports the MAXIMUM single-run
    score across all potential stem positions.
    
    Algorithm overview:
    ──────────────────
    1. Encode the probe sequence as 2-bit integers (A=00, C=01, G=10, T=11).
    2. Compute forward k-mer hashes using rolling bit-shift: O(L).
    3. Compute reverse-complement hashes via XOR + bit reversal: O(L).
    4. Build a hash→positions index for forward k-mers: O(L).
    5. For each position j, look up rc_hash[j] in the forward index.
       Each hit (i, j) means kmer[i] = revcomp(kmer[j]), indicating
       a potential stem pairing. Group by anti-diagonal d = i + j
       (consecutive stem extension keeps d constant): O(L + M).
    6. For each anti-diagonal, find the longest consecutive run of
       positions and apply exponential scoring: n consecutive → base^(n−1).
       Report the MAXIMUM single-run score across all diagonals.
    
    Why MAX instead of SUM:
    ───────────────────────
    A probe is problematic if it contains even ONE strong stem. Summing
    all matches (including scattered isolated ones) adds O(L²/4^k) noise
    that grows quadratically with probe length. Reporting the MAX run
    score captures only the strongest local structure — exactly what
    determines thermodynamic stability of the worst-case hairpin fold.
    
    Scoring:
    ────────
    n consecutive k-mer matches → stem of (n + k − 1) bases → base^(n−1) points.
    
    - Stem 4bp (n=1) → 4^0 =     1
    - Stem 5bp (n=2) → 4^1 =     4
    - Stem 6bp (n=3) → 4^2 =    16
    - Stem 7bp (n=4) → 4^3 =    64
    - Stem 8bp (n=5) → 4^4 =   256
    - Stem 10bp(n=7) → 4^6 = 4,096
    
    The exponential scaling ensures the score is DOMINATED by the longest
    contiguous stem. A 7bp stem (score 64) dwarfs any background from
    random 4-mer matches (score 1 each).
    
    Time complexity:
    ────────────────
    O(L + M), where L = sequence length, M = total k-mer matches.
    For random DNA with k=4: M ≈ L²/256 ≈ 50 for 120bp probe.
    Early stopping further reduces computation for bad probes.
    
    Args:
        seq: DNA sequence (uppercase, no ambiguous bases expected)
        k: k-mer size (default: 4, giving 256 possible k-mers)
        min_loop: minimum loop size between stem arms (default: 3bp)
        base: exponential base for consecutive bonus (default: 4.0)
        early_stop_raw: if set, return immediately when max >= this value
        
    Returns:
        Max exponential hairpin score (largest single stem contribution)
    """
    L = len(seq)
    n_kmers = L - k + 1
    
    if n_kmers <= 0:
        return 0.0
    
    # ─── Step 1: 2-bit encode ─────────────────────────────────────────
    encoded = [_BASE_TO_2BIT[ord(c)] for c in seq]
    
    # ─── Step 2: Rolling forward k-mer hashes via bit-shift ──────────
    mask = (1 << (2 * k)) - 1  # e.g. k=4 → 0xFF
    fwd_hashes = [0] * n_kmers
    
    h = 0
    for i in range(k):
        h = (h << 2) | encoded[i]
    fwd_hashes[0] = h & mask
    
    for i in range(1, n_kmers):
        h = ((h << 2) | encoded[i + k - 1]) & mask
        fwd_hashes[i] = h
    
    # ─── Step 3: Reverse-complement hashes ────────────────────────────
    rc_hashes = [_rc_hash(fh, k) for fh in fwd_hashes]
    
    # ─── Step 4: Build forward position index ─────────────────────────
    fwd_index: Dict[int, List[int]] = {}
    for i, fh in enumerate(fwd_hashes):
        if fh in fwd_index:
            fwd_index[fh].append(i)
        else:
            fwd_index[fh] = [i]
    
    # ─── Step 5: Match k-mer pairs, group by anti-diagonal ───────────
    #
    # A hairpin stem pairs position i (5′ arm) with position j (3′ arm):
    #   kmer[i] == revcomp(kmer[j])
    # Consecutive stem extension: i→i+1, j→j−1, so d = i + j is constant.
    #
    # min_gap ensures stem arms don't overlap and loop is ≥ min_loop:
    #   |i − j| ≥ k + min_loop
    #
    min_gap = k + min_loop
    
    diagonal_positions: Dict[int, List[int]] = {}
    
    for j in range(n_kmers):
        rc_h = rc_hashes[j]
        if rc_h not in fwd_index:
            continue
        for i in fwd_index[rc_h]:
            if abs(i - j) >= min_gap:
                d = i + j
                if d in diagonal_positions:
                    diagonal_positions[d].append(i)
                else:
                    diagonal_positions[d] = [i]
    
    # ─── Step 6: Find MAX run score across all diagonals ──────────────
    #
    # For each anti-diagonal (= one potential stem), find the longest
    # consecutive run of k-mer matches and compute its exponential score.
    # Return the MAXIMUM across all diagonals — this is the strongest
    # potential hairpin stem in the probe.
    #
    max_run_score = 0.0
    
    for d in diagonal_positions:
        # Deduplicate & sort (usually very short lists for real probes)
        positions = sorted(set(diagonal_positions[d]))
        
        run_length = 1
        for idx in range(1, len(positions)):
            if positions[idx] == positions[idx - 1] + 1:
                run_length += 1
            else:
                # Settle current run: n consecutive k-mer matches → base^(n−1)
                run_score = base ** (run_length - 1)
                if run_score > max_run_score:
                    max_run_score = run_score
                run_length = 1
        
        # Settle last run on this diagonal
        run_score = base ** (run_length - 1)
        if run_score > max_run_score:
            max_run_score = run_score
        
        # Early stopping: if max is already above threshold, no need to check more
        if early_stop_raw is not None and max_run_score >= early_stop_raw:
            return max_run_score
    
    return max_run_score


def _find_best_local_match(seq: str, rev_comp: str, min_match: int = 4) -> Tuple[int, int, int]:
    """
    Find the best local matching region between seq and its reverse complement.
    
    Uses a sliding window approach to find complementary regions.
    
    Returns:
        Tuple of (score, best_start, best_len)
    """
    n = len(seq)
    best_score = 0
    best_start = 0
    best_len = 0
    
    # Try different offsets (allowing for hairpin loop)
    for offset in range(min_match, n - min_match):
        # Compare seq[0:n-offset] with rev_comp[offset:n]
        match_len = n - offset
        if match_len < min_match:
            continue
        
        # Count matches in this alignment
        matches = 0
        current_run = 0
        max_run = 0
        
        for i in range(match_len):
            if seq[i] == rev_comp[offset + i]:
                matches += 1
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 0
        
        # Score based on matches and longest consecutive run
        score = matches * 2 + max_run * 3
        
        if score > best_score:
            best_score = score
            best_start = offset
            best_len = match_len
    
    return best_score, best_start, best_len


def calculate_hairpin_fast(
    sequence: str,
    method: str = "exp",
    match: int = 2,
    mismatch: int = -1,
    gap_open: int = 3,
    gap_extend: int = 1,
    kmer_size: int = 6,
    min_loop: int = 3,
    early_stop_raw: float = None,
) -> float:
    """
    Calculate hairpin formation potential of a probe sequence.
    
    Methods available:
    - "exp": Exponential bonus for consecutive k-mer matches (default, recommended)
    - "stem": Normalized max stem length (legacy simple method)
    - "kmer": K-mer sharing + local match composite score (legacy)
    - "align": Full local alignment (most accurate but slowest)
    
    The "exp" method (RECOMMENDED):
    ───────────────────────────────
    Scans all 4-mer reverse-complement matches across the probe, then
    rewards spatially CONSECUTIVE matches exponentially (4^(n-1) for n
    consecutive k-mer matches, corresponding to an (n+3)bp stem).
    
    This captures the biophysics of hairpin stability: a continuous stem
    has much higher binding energy than scattered random matches. The
    exponential scoring ensures that even in long probes, a real stem
    produces a score orders of magnitude above background noise.
    
    Score interpretation (exp method):
    - Random DNA: normalized ~1-6 (max chance stem rarely exceeds 6bp)
    - 6bp stem: normalized ~5 (marginally stable, usually passes)
    - 7bp stem: normalized ~18-22 (thermodynamically stable)
    - 8bp+ stem: normalized ≥74 (always caught)
    
    Default threshold: 18.0 (normalized). This catches stems ≥ 7bp
    (biologically stable hairpins) at probe lengths 40-120bp.
    
    Args:
        sequence: DNA sequence
        method: "exp" (default), "stem" (legacy), "kmer" (legacy), or "align"
        match: Match score for alignment method
        mismatch: Mismatch penalty for alignment method
        gap_open: Gap opening penalty for alignment method
        gap_extend: Gap extension penalty for alignment method
        kmer_size: K-mer size for kmer method (default: 6)
        min_loop: Minimum hairpin loop size (default: 3)
        early_stop_raw: For "exp" method: stop immediately when raw score
                       exceeds this value (saves time on bad probes)
        
    Returns:
        Hairpin score (higher = more likely to form hairpin).
        For "exp": normalized score = raw_exponential / log4(L).
        For "stem": normalized stem score = raw_stem / log4(L).
        For "kmer"/"align": composite score (legacy).
    """
    if not sequence:
        raise ValueError("Empty sequence provided")
    
    seq = sequence.upper()
    
    if method == "exp":
        # Exponential bonus method: k-mer scan with spatial continuity weighting.
        # Consecutive k-mer matches (stem structure) get exponential bonus:
        #   n consecutive 4-mer matches → 4^(n-1) points
        # Normalized by log4(L) for cross-length comparability.
        raw = _calculate_hairpin_exponential(
            seq, k=4, min_loop=min_loop, base=4.0,
            early_stop_raw=early_stop_raw,
        )
        log4_len = math.log(len(seq), 4) if len(seq) > 1 else 1.0
        return round(raw / log4_len, 2)
    
    elif method == "stem":
        # Legacy simple method: max consecutive stem length / log4(L)
        rev_comp = _reverse_complement(seq)
        raw_stem = _calculate_max_stem(seq, rev_comp, min_loop)
        log4_len = math.log(len(seq), 4) if len(seq) > 1 else 1.0
        return round(raw_stem / log4_len, 2)
    
    elif method == "kmer":
        # Legacy k-mer based method (NOT recommended for long sequences)
        rev_comp = _reverse_complement(seq)
        if len(seq) < kmer_size * 2:
            return 0.0
        
        seq_kmers = set(seq[i:i+kmer_size] for i in range(len(seq) - kmer_size + 1))
        rc_kmers = set(rev_comp[i:i+kmer_size] for i in range(len(rev_comp) - kmer_size + 1))
        
        shared = len(seq_kmers & rc_kmers)
        total = len(seq_kmers)
        
        score, _, _ = _find_best_local_match(seq, rev_comp, min_match=kmer_size)
        
        kmer_score = (shared / total * 50) if total > 0 else 0
        hairpin_score = kmer_score + score
        
        return round(hairpin_score, 2)
    
    elif method == "align":
        rev_comp = _reverse_complement(seq)
        if HAS_PARASAIL:
            matrix = parasail.matrix_create("ACGTN", match, mismatch)
            result = parasail.sw_striped_16(seq, rev_comp, gap_open, gap_extend, matrix)
            return float(result.score)
        else:
            aligner = PairwiseAligner()
            aligner.mode = "local"
            aligner.match_score = match
            aligner.mismatch_score = mismatch
            aligner.open_gap_score = -gap_open
            aligner.extend_gap_score = -gap_extend
            
            alignments = list(aligner.align(seq, rev_comp))
            
            if not alignments:
                return 0.0
            
            return float(max(a.score for a in alignments))
    
    else:
        raise ValueError(f"Unknown method: {method}. Use 'exp', 'stem', 'kmer', or 'align'.")


def calculate_hairpin_batch_fast(
    sequences: List[str],
    method: str = "stem",
    match: int = 2,
    mismatch: int = -1,
    gap_open: int = 3,
    gap_extend: int = 1,
    kmer_size: int = 6,
    threads: int = 1,
) -> List[float]:
    """
    Calculate hairpin scores for multiple sequences.
    
    Uses threading for parallel computation when method="align".
    The "stem" and "kmer" methods are already fast enough for single-threaded use.
    
    Args:
        sequences: List of DNA sequences
        method: "stem" (default), "kmer" (legacy), or "align"
        match: Match score for alignment method
        mismatch: Mismatch penalty for alignment method
        gap_open: Gap opening penalty for alignment method
        gap_extend: Gap extension penalty for alignment method
        kmer_size: K-mer size for kmer method
        threads: Number of threads for parallel processing (align method only)
        
    Returns:
        List of hairpin scores
    """
    if method in ("stem", "kmer") or threads <= 1:
        return [calculate_hairpin_fast(seq, method, match, mismatch, gap_open, gap_extend, kmer_size) 
                for seq in sequences]
    
    # Parallel processing with ThreadPoolExecutor for align method
    results = [None] * len(sequences)
    
    def compute(idx: int, seq: str) -> Tuple[int, float]:
        return idx, calculate_hairpin_fast(seq, method, match, mismatch, gap_open, gap_extend, kmer_size)
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(compute, i, seq) for i, seq in enumerate(sequences)]
        for future in as_completed(futures):
            idx, score = future.result()
            results[idx] = score
    
    return results


# =============================================================================
# Dimer Score - K-mer Frequency Sharing
# =============================================================================

@dataclass
class DimerCalculatorFast:
    """
    Fast calculator for dimer (inter-probe complementarity) scores.
    
    Pre-builds a k-mer index from all probes, then calculates each probe's
    dimer score based on how often its k-mers appear in the overall pool.
    
    Higher scores indicate greater potential for dimer formation with other probes.
    
    Usage:
        >>> calc = DimerCalculatorFast(k=11)
        >>> calc.build_index(sequences)
        >>> scores = calc.calculate_all_scores()
    """
    
    k: int = 11  # K-mer size
    include_revcomp: bool = True  # Include reverse complement k-mers
    
    _kmer_freq: Counter = field(default_factory=Counter, repr=False)
    _sequences: List[str] = field(default_factory=list, repr=False)
    _total_sequences: int = 0
    
    def __post_init__(self) -> None:
        """Initialize with empty state."""
        self._kmer_freq = Counter()
        self._sequences = []
        self._total_sequences = 0
    
    def build_index(self, sequences: List[str]) -> int:
        """
        Build k-mer frequency index from all probe sequences.
        
        Args:
            sequences: List of probe sequences
            
        Returns:
            Number of unique k-mers in index
        """
        self._sequences = [seq.upper() for seq in sequences]
        self._total_sequences = len(sequences)
        self._kmer_freq = Counter()
        
        for seq in self._sequences:
            # Count k-mers in forward sequence
            for i in range(len(seq) - self.k + 1):
                kmer = seq[i:i + self.k]
                self._kmer_freq[kmer] += 1
            
            # Also count reverse complement k-mers (for detecting sense/antisense dimers)
            if self.include_revcomp:
                rc_seq = _reverse_complement(seq)
                for i in range(len(rc_seq) - self.k + 1):
                    kmer = rc_seq[i:i + self.k]
                    self._kmer_freq[kmer] += 1
        
        return len(self._kmer_freq)
    
    def calculate_score(self, sequence: str) -> float:
        """
        Calculate dimer score for a single sequence.
        
        Score = sum of k-mer frequencies / total_sequences
        Higher scores mean more k-mers are shared with other probes.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dimer score (normalized)
        """
        if self._total_sequences == 0:
            raise ValueError("No sequences indexed. Call build_index first.")
        
        seq = sequence.upper()
        
        if len(seq) < self.k:
            return 0.0
        
        # Sum frequencies of k-mers in this sequence
        total_freq = 0
        n_kmers = 0
        
        for i in range(len(seq) - self.k + 1):
            kmer = seq[i:i + self.k]
            # Subtract 1 to not count the sequence's own contribution
            freq = self._kmer_freq.get(kmer, 0)
            if freq > 0:
                total_freq += freq - 1  # -1 to exclude self
            n_kmers += 1
        
        if n_kmers == 0:
            return 0.0
        
        # Normalize by number of k-mers and total sequences
        # This gives average k-mer sharing per probe
        score = total_freq / (n_kmers * self._total_sequences) * 100
        
        return round(score, 4)
    
    def calculate_all_scores(self) -> List[float]:
        """
        Calculate dimer scores for all indexed sequences.
        
        Returns:
            List of dimer scores (one per sequence)
        """
        return [self.calculate_score(seq) for seq in self._sequences]


def calculate_dimer_batch_fast(
    sequences: List[str],
    k: int = 11,
    include_revcomp: bool = True,
) -> List[float]:
    """
    Calculate dimer scores for a set of sequences.
    
    Convenience function that creates calculator, builds index,
    and returns all scores.
    
    Args:
        sequences: List of DNA sequences
        k: K-mer size (default: 11)
        include_revcomp: Whether to include reverse complement k-mers
        
    Returns:
        List of dimer scores
    """
    calc = DimerCalculatorFast(k=k, include_revcomp=include_revcomp)
    calc.build_index(sequences)
    return calc.calculate_all_scores()


# =============================================================================
# Unified Batch Calculator
# =============================================================================

@dataclass
class BiophysicalStats:
    """Statistics for a single probe sequence."""
    gc: float
    tm: float
    dust: float
    hairpin: float
    dimer: float = 0.0  # Requires pool-level calculation


def calculate_all_stats_fast(
    sequences: List[str],
    threads: int = 1,
    na_conc: float = 50.0,
    dimer_k: int = 11,
) -> List[BiophysicalStats]:
    """
    Calculate all biophysical statistics for a set of sequences.
    
    This is the main entry point for batch biophysical calculation.
    Optimized for processing large numbers of sequences.
    
    Args:
        sequences: List of DNA sequences
        threads: Number of threads for hairpin calculation
        na_conc: Na+ concentration for Tm calculation
        dimer_k: K-mer size for dimer calculation
        
    Returns:
        List of BiophysicalStats objects
    """
    n = len(sequences)
    
    if n == 0:
        return []
    
    # Calculate all metrics
    gc_scores = calculate_gc_batch_fast(sequences)
    tm_scores = calculate_tm_batch_fast(sequences, na_conc)
    dust_scores = calculate_dust_batch_fast(sequences)
    hairpin_scores = calculate_hairpin_batch_fast(sequences, threads=threads)
    dimer_scores = calculate_dimer_batch_fast(sequences, k=dimer_k)
    
    # Combine into stats objects
    return [
        BiophysicalStats(
            gc=gc_scores[i],
            tm=tm_scores[i],
            dust=dust_scores[i],
            hairpin=hairpin_scores[i],
            dimer=dimer_scores[i],
        )
        for i in range(n)
    ]


# =============================================================================
# Percentile-based Thresholding
# =============================================================================

def calculate_percentile_threshold(
    scores: List[float],
    percentile: float,
    higher_is_worse: bool = True,
) -> float:
    """
    Calculate threshold value at given percentile.
    
    Args:
        scores: List of scores
        percentile: Percentile (0-100)
        higher_is_worse: If True, threshold is upper bound; if False, lower bound
        
    Returns:
        Threshold value
    """
    if not scores:
        return 0.0
    
    sorted_scores = sorted(scores)
    n = len(sorted_scores)
    
    # Calculate percentile index
    idx = int((percentile / 100.0) * n)
    idx = min(idx, n - 1)
    
    return sorted_scores[idx]
