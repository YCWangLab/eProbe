"""
Optimized biophysical property calculators for probe sequences.

This module provides fast implementations of biophysical calculations:
- GC content: Vectorized numpy calculation
- Tm: Nearest-neighbor with lookup table optimization
- DUST: Optimized k-mer counting
- Hairpin: ViennaRNA MFE (DNA Mathews 2004 parameters)
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



# =============================================================================
# ViennaRNA for hairpin MFE calculation
# =============================================================================

try:
    import RNA as _RNA
    _RNA.params_load_DNA_Mathews2004()
    _HAS_VIENNA = True
except ImportError:
    _HAS_VIENNA = False


def calculate_hairpin_fast(
    sequence: str,
    **_kwargs,
) -> float:
    """
    Calculate hairpin formation potential as |MFE| (kcal/mol) via ViennaRNA.

    Uses the Zuker MFE algorithm with DNA thermodynamic parameters
    (Mathews 2004). Returns the absolute value of the minimum free energy
    so that higher values indicate stronger self-folding (worse for probes).

    No sequence length limit. ~1.8 ms per 81 bp probe (C implementation).

    Legacy keyword arguments (method, match, mismatch, etc.) are accepted
    but silently ignored for backward compatibility.

    Args:
        sequence: DNA sequence

    Returns:
        |MFE| in kcal/mol (0.0 = no structure, higher = more hairpin)

    Raises:
        RuntimeError: If ViennaRNA is not installed.
    """
    if not _HAS_VIENNA:
        raise RuntimeError(
            "ViennaRNA is not installed. Install with: pip install ViennaRNA"
        )
    if not sequence:
        raise ValueError("Empty sequence provided")
    seq = sequence.upper().replace("N", "A")
    _ss, mfe = _RNA.fold(seq)
    return round(abs(mfe), 2)


def calculate_hairpin_batch_fast(
    sequences: List[str],
    threads: int = 1,
    **_kwargs,
) -> List[float]:
    """
    Calculate hairpin |MFE| scores for multiple sequences via ViennaRNA.

    Args:
        sequences: List of DNA sequences
        threads: Ignored (ViennaRNA is already fast in single-thread).

    Returns:
        List of |MFE| scores (kcal/mol)
    """
    return [calculate_hairpin_fast(seq) for seq in sequences]


# =============================================================================
# Dimer Score - K-mer Frequency Sharing
# =============================================================================

@dataclass
class DimerCalculatorFast:
    """
    Fast calculator for dimer (inter-probe complementarity) scores.
    
    Indexes only the reverse-complement k-mers of all probes, then scores
    each probe by querying its forward k-mers against this RC index.
    A hit means another probe's RC contains that k-mer → potential
    Watson-Crick hybridization (dimer formation).
    
    Self-complementarity (own RC hits) is subtracted so the score
    reflects only inter-probe dimer risk (hairpin is a separate metric).
    
    Higher scores indicate greater potential for dimer formation with other probes.
    
    Usage:
        >>> calc = DimerCalculatorFast(k=11)
        >>> calc.build_index(sequences)
        >>> scores = calc.calculate_all_scores()
    """
    
    k: int = 11  # K-mer size
    
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
        Build k-mer frequency index from reverse complements of all probes.

        Only RC k-mers are indexed so that querying a probe's forward k-mers
        against this index directly counts complementary matches (i.e. which
        other probes could hybridize to this one via Watson-Crick pairing).

        Args:
            sequences: List of probe sequences

        Returns:
            Number of unique k-mers in index
        """
        self._sequences = [seq.upper() for seq in sequences]
        self._total_sequences = len(sequences)
        self._kmer_freq = Counter()

        for seq in self._sequences:
            rc_seq = _reverse_complement(seq)
            for i in range(len(rc_seq) - self.k + 1):
                kmer = rc_seq[i:i + self.k]
                self._kmer_freq[kmer] += 1

        return len(self._kmer_freq)

    def calculate_score(self, sequence: str) -> float:
        """
        Calculate dimer score for a single sequence.

        Queries the probe's forward k-mers against the RC index.
        Subtracts the probe's own RC contribution (self-complementarity
        is hairpin, not inter-probe dimer).

        Score = (sum of other-probe RC hits) / (n_kmers × N_probes) × 10000

        Args:
            sequence: DNA sequence

        Returns:
            Dimer score (normalized, higher = more complementary to pool)
        """
        if self._total_sequences == 0:
            raise ValueError("No sequences indexed. Call build_index first.")

        seq = sequence.upper()

        if len(seq) < self.k:
            return 0.0

        # Count this probe's own RC k-mers so we can subtract self-contribution
        rc_seq = _reverse_complement(seq)
        rc_self: Dict[str, int] = {}
        for i in range(len(rc_seq) - self.k + 1):
            kmer = rc_seq[i:i + self.k]
            rc_self[kmer] = rc_self.get(kmer, 0) + 1

        # Query forward k-mers against RC index
        total_freq = 0
        n_kmers = 0

        for i in range(len(seq) - self.k + 1):
            kmer = seq[i:i + self.k]
            freq = self._kmer_freq.get(kmer, 0)
            # Subtract own RC contribution (self-complementarity = hairpin, not dimer)
            self_rc_count = rc_self.get(kmer, 0)
            other_freq = max(freq - self_rc_count, 0)
            total_freq += other_freq
            n_kmers += 1

        if n_kmers == 0:
            return 0.0

        # Normalize by number of k-mers and total sequences
        score = total_freq / (n_kmers * self._total_sequences) * 10000

        return round(score, 2)
    
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
) -> List[float]:
    """
    Calculate dimer scores for a set of sequences.
    
    Convenience function that creates calculator, builds index,
    and returns all scores.
    
    Args:
        sequences: List of DNA sequences
        k: K-mer size (default: 11)
        
    Returns:
        List of dimer scores
    """
    calc = DimerCalculatorFast(k=k)
    calc.build_index(sequences)
    return calc.calculate_all_scores()


# =============================================================================
# Smart Dimer Filter - Graph-Based Deduplication
# =============================================================================

class SmartDimerFilter:
    """
    Graph-based dimer filter: identifies groups of complementary probes
    and keeps only the best probe from each dimer-risk group.
    
    Instead of removing all probes above a score threshold, this approach:
    1. Builds a canonical k-mer index (canonical = min(kmer, rc(kmer)))
    2. Finds probe pairs sharing significant k-mer content (potential dimers)
    3. Groups connected probes via Union-Find
    4. Keeps one representative (the least-connected) from each group
    
    This maximizes retained probe count while eliminating dimer risk clusters.
    
    Args:
        k: K-mer size for fingerprinting (default: 11)
        min_shared_fraction: Minimum fraction of shared canonical k-mers
            to consider two probes "dimer partners" (default: 0.15 = 15%).
            For 81bp probes with ~60 unique canonical 11-mers, 0.15 means
            ~9 shared k-mers ≈ a 19bp complementary stretch.
        max_kmer_freq: Skip canonical k-mers appearing in more probes than
            this (likely repetitive, not informative). Default: 100.
    """
    
    def __init__(self, k: int = 11, min_shared_fraction: float = 0.15,
                 max_kmer_freq: int = 100):
        self.k = k
        self.min_shared_fraction = min_shared_fraction
        self.max_kmer_freq = max_kmer_freq
    
    def filter(self, sequences: List[str]) -> Tuple[List[int], Dict]:
        """
        Identify dimer-risk groups and return indices of probes to keep.
        
        Args:
            sequences: List of probe sequences
            
        Returns:
            Tuple of:
            - Sorted list of indices to keep
            - Stats dict with dimer group information
        """
        from collections import defaultdict
        
        n = len(sequences)
        if n == 0:
            return [], {
                "dimer_failed": 0, "dimer_groups": 0,
                "dimer_edges": 0, "max_group_size": 0,
                "mean_group_size": 0,
                "mode": f"smart(k={self.k}, shared≥{self.min_shared_fraction:.0%})",
            }
        
        # Step 1: Build canonical k-mer inverted index
        # canonical(K) = min(K, rc(K)) — groups forward and RC together
        kmer_to_probes: Dict[str, set] = defaultdict(set)
        probe_n_kmers: List[int] = []  # unique canonical k-mers per probe
        
        for i, seq in enumerate(sequences):
            seq_upper = seq.upper()
            seen: set = set()
            for j in range(len(seq_upper) - self.k + 1):
                kmer = seq_upper[j:j + self.k]
                if 'N' in kmer:
                    continue
                rc = _reverse_complement(kmer)
                canonical = kmer if kmer <= rc else rc
                if canonical not in seen:
                    kmer_to_probes[canonical].add(i)
                    seen.add(canonical)
            probe_n_kmers.append(max(len(seen), 1))
        
        # Step 2: Count pairwise shared canonical k-mers
        pair_sharing: Counter = Counter()
        
        for _kmer, probe_set in kmer_to_probes.items():
            ps = len(probe_set)
            if ps < 2 or ps > self.max_kmer_freq:
                continue
            probes = sorted(probe_set)
            for a_idx in range(ps):
                for b_idx in range(a_idx + 1, ps):
                    pair_sharing[(probes[a_idx], probes[b_idx])] += 1
        
        # Step 3: Build adjacency from significant pairs
        adj: Dict[int, set] = defaultdict(set)
        n_edges = 0
        
        for (i, j), count in pair_sharing.items():
            frac = count / min(probe_n_kmers[i], probe_n_kmers[j])
            if frac >= self.min_shared_fraction:
                adj[i].add(j)
                adj[j].add(i)
                n_edges += 1
        
        # Free memory — pair_sharing can be large
        del pair_sharing
        
        # Step 4: Union-Find for connected components
        parent = list(range(n))
        uf_rank = [0] * n
        
        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]  # path compression
                x = parent[x]
            return x
        
        def union(x: int, y: int) -> None:
            px, py = find(x), find(y)
            if px == py:
                return
            if uf_rank[px] < uf_rank[py]:
                px, py = py, px
            parent[py] = px
            if uf_rank[px] == uf_rank[py]:
                uf_rank[px] += 1
        
        for i, neighbors in adj.items():
            for j in neighbors:
                union(i, j)
        
        # Step 5: Group by connected component
        components: Dict[int, List[int]] = defaultdict(list)
        for i in range(n):
            components[find(i)].append(i)
        
        # Step 6: Keep best from each group
        # "Best" = fewest dimer edges (least cross-hybridization risk)
        keep_indices: set = set()
        n_groups = 0
        group_sizes: List[int] = []
        
        for _root, members in components.items():
            if len(members) == 1:
                keep_indices.add(members[0])
            else:
                n_groups += 1
                group_sizes.append(len(members))
                best = min(members, key=lambda idx: len(adj.get(idx, set())))
                keep_indices.add(best)
        
        keep_sorted = sorted(keep_indices)
        n_removed = n - len(keep_sorted)
        
        stats = {
            "dimer_failed": n_removed,
            "dimer_groups": n_groups,
            "dimer_edges": n_edges,
            "max_group_size": max(group_sizes) if group_sizes else 0,
            "mean_group_size": round(sum(group_sizes) / len(group_sizes), 1) if group_sizes else 0,
            "mode": f"smart(k={self.k}, shared≥{self.min_shared_fraction:.0%})",
        }
        
        return keep_sorted, stats


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
