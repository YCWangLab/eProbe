"""
SNP filtering module.

Implements multi-stage filtering for SNP probe candidates:
  - Background noise filtering (Kraken2)
  - Accessibility filtering (Bowtie2)
  - Taxonomic filtering (NCBI taxonomy)
  - Biophysical filtering (GC, Tm, complexity, hairpin, dimer)

This module corresponds to the original SNP_filter.py functionality.
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Dict, Any, List, Set, Callable, Tuple
from dataclasses import dataclass, field
from enum import Enum

import pandas as pd
import numpy as np

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNP, SNPDataFrame
from eprobe.core.fasta import write_fasta
from eprobe.biophysics import (
    calculate_gc,
    calculate_tm,
    calculate_complexity,
    calculate_hairpin_score,
)
from eprobe.biophysics.dimer import DimerCalculator

logger = logging.getLogger(__name__)


class FilterType(Enum):
    """Available filter types."""
    BG = "bg"              # Background noise (Kraken2)
    AC = "ac"              # Accessibility (Bowtie2)
    TX = "tx"              # Taxonomy
    BIOPHYSICAL = "biophysical"


@dataclass
class BiophysicalThresholds:
    """Thresholds for biophysical filtering."""
    gc_min: float = 35.0
    gc_max: float = 65.0
    tm_min: float = 55.0
    tm_max: float = 75.0
    complexity_max: float = 2.0
    hairpin_max: Optional[float] = None  # Optional hairpin threshold
    dimer_max: Optional[float] = None    # Optional dimer threshold


@dataclass
class FilterConfig:
    """Configuration for all filter stages."""
    enabled_filters: List[str] = field(default_factory=lambda: ["biophysical"])
    bg_db: Optional[Path] = None
    ac_db: Optional[Path] = None
    tx_db: Optional[Path] = None
    tx_ids: Optional[List[int]] = None
    biophysical: BiophysicalThresholds = field(default_factory=BiophysicalThresholds)
    threads: int = 1


# =============================================================================
# Background Noise Filter (Kraken2)
# =============================================================================

def run_kraken2(
    fasta_path: Path,
    db_path: Path,
    output_path: Path,
    threads: int = 1,
) -> Result[None, str]:
    """
    Run Kraken2 classification on probe sequences and save all output files.
    
    Args:
        fasta_path: Input FASTA file
        db_path: Kraken2 database path
        output_path: Output classification file
        threads: Number of threads
        
    Returns:
        Result indicating success/failure
    """
    # Create additional output files
    report_path = output_path.parent / "kraken2.report"
    classified_path = output_path.parent / "kraken2.classified.fasta"
    unclassified_path = output_path.parent / "kraken2.unclassified.fasta"
    
    cmd = [
        "kraken2",
        "--db", str(db_path),
        "--threads", str(threads),
        "--output", str(output_path),
        "--report", str(report_path),
        "--classified-out", str(classified_path),
        "--unclassified-out", str(unclassified_path),
        str(fasta_path),
    ]
    
    logger.info(f"Running Kraken2: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        
        if not output_path.exists():
            return Err(f"Kraken2 output file not created: {output_path}")
        
        logger.info(f"Kraken2 completed successfully. Output files:")
        logger.info(f"  - Main output: {output_path}")
        logger.info(f"  - Report: {report_path}")
        logger.info(f"  - Classified: {classified_path}")
        logger.info(f"  - Unclassified: {unclassified_path}")
        
        # Print first few lines of main output for debugging
        logger.info("First 5 lines of kraken2 output:")
        try:
            with open(output_path, 'r') as f:
                for i, line in enumerate(f):
                    if i < 5:
                        logger.info(f"  {line.strip()}")
                    else:
                        break
        except Exception as e:
            logger.warning(f"Could not read output file for preview: {e}")
        
        return Ok(None)
        
    except subprocess.CalledProcessError as e:
        return Err(f"Kraken2 failed: {e}\nStdout: {e.stdout}\nStderr: {e.stderr}")
    except Exception as e:
        return Err(f"Kraken2 execution error: {e}")


def parse_kraken_kmer_matches(output_path: Path, match_threshold: int = 10) -> Result[tuple[Set[str], Dict[str, int]], str]:
    """
    Parse Kraken2 output as simple kmer matcher (ignore taxonomy classification).
    
    Counts how many kmers match to non-zero taxid in column 5.
    Column 5 format: "taxid:count taxid:count ..."
    Examples:
        - "0:35 586:1 0:29"  -> 1 kmer matched (586:1)
        - "562:15 0:50 585:8" -> 23 kmers matched (15+8)
        - "0:65"              -> 0 kmers matched
        - "A:35 0:30"         -> 0 kmers matched (A=ambiguous, not real match)
    
    Args:
        output_path: Kraken2 main output file
        match_threshold: Remove sequences with >= this many matched kmers
        
    Returns:
        Tuple of (sequences_to_remove, kmer_match_counts)
        - sequences_to_remove: Set of sequence IDs exceeding threshold
        - kmer_match_counts: Dict {seq_id: matched_kmer_count}
    """
    try:
        sequences_to_remove = set()
        kmer_match_counts = {}
        
        with open(output_path, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                
                seq_id = parts[1]
                lca_mappings = parts[4]  # Column 5: kmer match details
                
                # Count kmers matched to database (non-zero, non-ambiguous taxid)
                matched_kmers = 0
                
                for mapping in lca_mappings.split():
                    if ':' not in mapping:
                        continue
                    
                    try:
                        taxid_str, count_str = mapping.split(':', 1)
                        
                        # Skip special markers
                        if taxid_str in ['0', 'A']:  # 0=unassigned, A=ambiguous
                            continue
                        
                        # Count real database matches (taxid > 0)
                        taxid = int(taxid_str)
                        count = int(count_str)
                        
                        if taxid > 0:
                            matched_kmers += count
                            
                    except (ValueError, IndexError):
                        continue
                
                kmer_match_counts[seq_id] = matched_kmers
                
                # Mark for removal if exceeds threshold
                if matched_kmers >= match_threshold:
                    sequences_to_remove.add(seq_id)
        
        logger.info(f"Kmer match analysis: {len(kmer_match_counts)} sequences analyzed")
        logger.info(f"Match threshold: {match_threshold} kmers")
        logger.info(f"Sequences to remove: {len(sequences_to_remove)} (>= {match_threshold} matched kmers)")
        
        return Ok((sequences_to_remove, kmer_match_counts))
        
    except Exception as e:
        return Err(f"Failed to parse kraken kmer matches: {e}")


def display_kraken_summary(kmer_match_counts: Dict[str, int], match_threshold: int) -> None:
    """
    Display kmer match statistics summary.
    
    Args:
        kmer_match_counts: Dict {seq_id: matched_kmer_count}
        match_threshold: Threshold used for filtering
    """
    if not kmer_match_counts:
        return
    
    matched_counts = list(kmer_match_counts.values())
    total_seqs = len(matched_counts)
    
    logger.info("=" * 60)
    logger.info("=== Kraken2 Kmer Match Analysis ===")
    logger.info(f"Total sequences: {total_seqs}")
    logger.info(f"Filtering threshold: >= {match_threshold} matched kmers")
    
    logger.info(f"\nMatched kmer statistics:")
    logger.info(f"  Min:    {min(matched_counts)}")
    logger.info(f"  Max:    {max(matched_counts)}")
    logger.info(f"  Mean:   {sum(matched_counts)/len(matched_counts):.1f}")
    logger.info(f"  Median: {sorted(matched_counts)[len(matched_counts)//2]}")
    
    # Distribution
    bins = [(0, 0), (1, 5), (6, 10), (11, 20), (21, 50), (51, float('inf'))]
    logger.info(f"\nMatch distribution:")
    for low, high in bins:
        if high == float('inf'):
            count = sum(1 for c in matched_counts if c > low)
            label = f"  >{low:>3} kmers"
        elif low == high:
            count = sum(1 for c in matched_counts if c == low)
            label = f"  {low:>3} kmers "
        else:
            count = sum(1 for c in matched_counts if low <= c <= high)
            label = f"  {low:>3}-{high:<3}   "
        pct = count / total_seqs * 100
        logger.info(f"{label}: {count:>6} seqs ({pct:>5.1f}%)")
    
    logger.info("=" * 60)


def filter_background_noise(
    snps: List[SNP],
    db_path: Path,
    fasta_path: Optional[Path] = None,
    threads: int = 1,
    kmer_threshold: int = 1,
) -> Result[List[SNP], str]:
    """
    Filter SNPs using Kraken2 as simple kmer matcher (ignore taxonomy).
    
    Removes sequences with >= kmer_threshold kmers matching background database.
    Kraken2 column 5 analysis: counts kmers matching non-zero taxid.
    
    Args:
        snps: Input SNP list
        db_path: Kraken2 database path
        fasta_path: Optional pre-computed FASTA file
        threads: Number of threads
        kmer_threshold: Remove sequences with >= this many matched kmers (default: 1)
        
    Returns:
        Result containing filtered SNP list
    """
    if not snps:
        return Ok([])
    
    logger.info(f"Running background noise filter (Kraken2 kmer matching)")
    logger.info(f"Filtering threshold: >= {kmer_threshold} matched kmers")
    logger.info(f"Processing {len(snps)} SNPs")
    
    # Use shared FASTA or create temporary one
    cleanup_fasta = False
    if fasta_path is None:
        # Create temporary FASTA
        import tempfile
        temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
        fasta_path = Path(temp_fasta.name)
        cleanup_fasta = True
        
        sequences = {snp.id: probe_sequences[snp.id] for snp in snps}
        write_result = write_fasta(sequences, fasta_path)
        if write_result.is_err():
            return Err(f"Failed to write temp FASTA: {write_result.unwrap_err()}")
    
    # Create temporary output file
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.kraken.out', delete=False) as temp_out:
        kraken_output_path = Path(temp_out.name)
    
    try:
        # Run Kraken2
        kraken_result = run_kraken2(fasta_path, db_path, kraken_output_path, threads)
        if kraken_result.is_err():
            return Err(kraken_result.unwrap_err())
        
        # Parse kmer matches
        parse_result = parse_kraken_kmer_matches(kraken_output_path, kmer_threshold)
        if parse_result.is_err():
            return Err(parse_result.unwrap_err())
        
        sequences_to_remove, kmer_match_counts = parse_result.unwrap()
        
        # Display statistics
        display_kraken_summary(kmer_match_counts, kmer_threshold)
        
        # Filter SNPs
        filtered_snps = [snp for snp in snps if snp.id not in sequences_to_remove]
        removed = len(snps) - len(filtered_snps)
        
        logger.info(f"Background filter: removed {removed} SNPs, {len(filtered_snps)} remaining")
        
        return Ok(filtered_snps)
        
    finally:
        # Clean up temporary files
        try:
            if cleanup_fasta and fasta_path.exists():
                fasta_path.unlink()
                logger.debug(f"Cleaned up temp FASTA: {fasta_path}")
            
            if kraken_output_path.exists():
                kraken_output_path.unlink()
                logger.debug(f"Cleaned up kraken output: {kraken_output_path}")
            
            # Clean up additional kraken files
            report_path = kraken_output_path.parent / "kraken2.report"
            classified_path = kraken_output_path.parent / "kraken2.classified.fasta"
            unclassified_path = kraken_output_path.parent / "kraken2.unclassified.fasta"
            
            for cleanup_file in [report_path, classified_path, unclassified_path]:
                if cleanup_file.exists():
                    cleanup_file.unlink()
                    logger.debug(f"Cleaned up: {cleanup_file}")
                    
        except Exception as e:
            logger.warning(f"Failed to clean up temporary files: {e}")


# =============================================================================
# Accessibility Filter (Bowtie2)
# =============================================================================
# Accessibility Filter (Bowtie2) - Multi-genome Universality + Mappability
# =============================================================================

def run_bowtie2(
    fasta_path: Path,
    index_path: Path,
    output_path: Path,
    threads: int = 1,
    report_all: bool = False,
) -> Result[Path, str]:
    """
    Run Bowtie2 alignment for accessibility check.
    
    Args:
        fasta_path: Input FASTA file
        index_path: Bowtie2 index prefix
        output_path: Output SAM file
        threads: Number of threads
        report_all: If True, report all alignments (-a); else report up to 10 (-k 10)
        
    Returns:
        Result containing path to output file
    """
    cmd = [
        "bowtie2",
        "-x", str(index_path),
        "-f", str(fasta_path),
        "-S", str(output_path),
        "-p", str(threads),
    ]
    
    if report_all:
        cmd.append("-a")  # Report all alignments
    else:
        cmd.extend(["-k", "10"])  # Report up to 10 alignments
    
    logger.debug(f"Running Bowtie2: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        return Ok(output_path)
    except subprocess.CalledProcessError as e:
        return Err(f"Bowtie2 failed: {e.stderr}")
    except FileNotFoundError:
        return Err("Bowtie2 not found. Please install Bowtie2 and ensure it's in PATH.")


def parse_bowtie2_accessibility(
    sam_path: Path,
    mode: str = "strict",
    score_diff_threshold: int = 10,
) -> Result[Dict[str, Dict], str]:
    """
    Parse Bowtie2 SAM output for accessibility analysis.
    
    Returns detailed alignment info for each sequence:
    - aligned: whether the sequence aligned to this genome
    - hit_count: number of alignment positions
    - best_score: best alignment score (AS tag)
    - second_score: second-best alignment score (XS tag)
    - score_diff: difference between best and second-best
    - mappable: whether the sequence passes mappability check
    
    Args:
        sam_path: Path to SAM file
        mode: "strict" (unique only) or "relaxed" (allow multi if distinguishable)
        score_diff_threshold: Minimum score difference for relaxed mode
        
    Returns:
        Dict {seq_id: alignment_info}
    """
    results: Dict[str, Dict] = {}
    
    try:
        with open(sam_path) as f:
            for line in f:
                if line.startswith("@"):  # Skip header
                    continue
                    
                parts = line.strip().split("\t")
                if len(parts) < 11:
                    continue
                
                seq_id = parts[0]
                flag = int(parts[1])
                
                # Initialize result for this sequence
                if seq_id not in results:
                    results[seq_id] = {
                        "aligned": False,
                        "hit_count": 0,
                        "best_score": None,
                        "second_score": None,
                        "score_diff": None,
                        "mappable": False,
                    }
                
                # Check if mapped (flag bit 4 = unmapped)
                if flag & 4:  # Unmapped
                    continue
                
                results[seq_id]["aligned"] = True
                results[seq_id]["hit_count"] += 1
                
                # Parse optional tags for alignment scores
                tags = {tag.split(":")[0]: tag.split(":")[-1] for tag in parts[11:] if ":" in tag}
                
                # AS = alignment score of current alignment
                # XS = alignment score of second-best alignment
                if "AS" in tags:
                    current_score = int(tags["AS"])
                    
                    if results[seq_id]["best_score"] is None:
                        results[seq_id]["best_score"] = current_score
                    elif current_score > results[seq_id]["best_score"]:
                        # This is now the best, old best becomes second
                        results[seq_id]["second_score"] = results[seq_id]["best_score"]
                        results[seq_id]["best_score"] = current_score
                    elif results[seq_id]["second_score"] is None or current_score > results[seq_id]["second_score"]:
                        results[seq_id]["second_score"] = current_score
                
                # XS tag directly gives second-best score (if available)
                if "XS" in tags and results[seq_id]["second_score"] is None:
                    results[seq_id]["second_score"] = int(tags["XS"])
        
        # Calculate mappability for each sequence
        for seq_id, info in results.items():
            if not info["aligned"]:
                info["mappable"] = False
                continue
            
            # Calculate score difference
            if info["best_score"] is not None and info["second_score"] is not None:
                info["score_diff"] = info["best_score"] - info["second_score"]
            else:
                info["score_diff"] = None  # Only one alignment, perfectly unique
            
            # Determine mappability based on mode
            if mode == "strict":
                # Strict: must have exactly 1 alignment position
                info["mappable"] = info["hit_count"] == 1
            elif mode == "relaxed":
                # Relaxed: unique OR (multiple with sufficient score difference)
                if info["hit_count"] == 1:
                    info["mappable"] = True
                elif info["score_diff"] is not None and info["score_diff"] >= score_diff_threshold:
                    info["mappable"] = True
                else:
                    info["mappable"] = False
            else:
                info["mappable"] = info["hit_count"] >= 1  # Any alignment counts
        
        return Ok(results)
        
    except Exception as e:
        return Err(f"Failed to parse Bowtie2 output: {e}")


def filter_accessibility(
    snps: List[SNP],
    index_paths: List[Path],
    threads: int = 1,
    mode: str = "strict",
    score_diff_threshold: int = 10,
    min_genomes: Optional[int] = None,
    fasta_path: Optional[Path] = None,
    probe_sequences: Optional[Dict[str, str]] = None,
) -> Result[List[SNP], str]:
    """
    Filter SNPs based on genomic accessibility across multiple genomes.
    
    Two-step filtering:
    1. Universality: Probe must align to at least min_genomes (default: ALL)
    2. Mappability: Probe must pass mappability check in ALL aligned genomes
    
    Args:
        snps: Input SNP list
        index_paths: List of Bowtie2 index prefix paths (multiple genomes)
        threads: Number of threads
        mode: "strict" (unique alignment only) or "relaxed" (allow distinguishable multi-mapping)
        score_diff_threshold: Minimum AS-XS score difference for relaxed mode
        min_genomes: Minimum number of genomes to align to (default: len(index_paths))
        fasta_path: Path to existing FASTA file with probe sequences (optional)
        probe_sequences: Dict mapping SNP ID to probe sequence (optional, used if fasta_path not provided)
        
    Returns:
        Result containing filtered SNP list
    """
    if not snps:
        return Ok([])
    
    if not index_paths:
        return Err("No genome indices provided for accessibility filter")
    
    if fasta_path is None and probe_sequences is None:
        return Err("Either fasta_path or probe_sequences must be provided")
    
    # Set default min_genomes to all genomes
    if min_genomes is None:
        min_genomes = len(index_paths)
    elif min_genomes > len(index_paths):
        return Err(f"min_genomes ({min_genomes}) cannot exceed total genomes ({len(index_paths)})")
    elif min_genomes < 1:
        return Err(f"min_genomes must be at least 1")
    
    logger.info(f"Running accessibility filter on {len(snps)} SNPs")
    logger.info(f"Mode: {mode}, Score diff threshold: {score_diff_threshold}")
    logger.info(f"Genomes to check: {len(index_paths)}")
    logger.info(f"Minimum genomes required: {min_genomes}")
    
    # Distribute threads across genomes
    threads_per_genome = max(1, threads // len(index_paths))
    logger.info(f"Threads per genome: {threads_per_genome} (total: {threads})")
    
    # Track which SNPs pass each genome
    # snp_id -> genome_index -> alignment_info
    all_results: Dict[str, Dict[int, Dict]] = {snp.id: {} for snp in snps}
    
    # Use existing FASTA or create a temp file from probe_sequences
    use_temp_fasta = fasta_path is None
    
    if use_temp_fasta:
        import tempfile as tf
        temp_dir = tf.mkdtemp()
        fasta_path = Path(temp_dir) / "probes.fa"
        
        # Write probe sequences
        sequences = {snp.id: probe_sequences[snp.id] for snp in snps}
        write_result = write_fasta(sequences, fasta_path)
        if write_result.is_err():
            return Err(f"Failed to write temp FASTA: {write_result.unwrap_err()}")
    
    try:
        # Run Bowtie2 against each genome
        for genome_idx, index_path in enumerate(index_paths):
            genome_name = index_path.stem if hasattr(index_path, 'stem') else str(index_path)
            logger.info(f"  Checking genome {genome_idx + 1}/{len(index_paths)}: {genome_name}")
            
            sam_path = fasta_path.parent / f"genome_{genome_idx}.sam"
            
            # Run Bowtie2 with distributed threads
            bt2_result = run_bowtie2(fasta_path, index_path, sam_path, threads_per_genome)
            if bt2_result.is_err():
                return Err(bt2_result.unwrap_err())
            
            # Parse results
            parse_result = parse_bowtie2_accessibility(sam_path, mode, score_diff_threshold)
            if parse_result.is_err():
                return Err(parse_result.unwrap_err())
            
            genome_results = parse_result.unwrap()
            
            # Store results for this genome
            for seq_id in all_results:
                if seq_id in genome_results:
                    all_results[seq_id][genome_idx] = genome_results[seq_id]
                else:
                    # Sequence not in results = not aligned
                    all_results[seq_id][genome_idx] = {
                        "aligned": False,
                        "hit_count": 0,
                        "mappable": False,
                    }
            
            # Clean up SAM file
            sam_path.unlink()
    finally:
        # Clean up temp FASTA if we created it
        if use_temp_fasta:
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)
    
    # Filter SNPs based on results across all genomes
    n_genomes = len(index_paths)
    passed_snps = []
    
    # Detailed statistics
    stats = {
        "passed": 0,
        "insufficient_genomes": 0,  # Aligned to < min_genomes
        "multi_mapped": 0,           # Multi-mapped in some genome (not mappable)
    }
    
    # Per-genome statistics
    genome_stats = {
        "aligned_count": [0] * n_genomes,
        "mappable_count": [0] * n_genomes,
        "multi_mapped_count": [0] * n_genomes,  # Count of multi-mapped probes per genome
    }
    
    for snp in snps:
        snp_results = all_results[snp.id]
        
        # Count aligned and mappable genomes
        aligned_genomes = []
        mappable_genomes = []
        multi_mapped_genomes = []
        
        for i in range(n_genomes):
            result = snp_results.get(i, {})
            if result.get("aligned", False):
                aligned_genomes.append(i)
                genome_stats["aligned_count"][i] += 1
                
                if result.get("mappable", False):
                    mappable_genomes.append(i)
                    genome_stats["mappable_count"][i] += 1
                else:
                    # Aligned but not mappable = multi-mapped
                    multi_mapped_genomes.append(i)
                    genome_stats["multi_mapped_count"][i] += 1
        
        # Check 1: Must align to at least min_genomes
        if len(aligned_genomes) < min_genomes:
            stats["insufficient_genomes"] += 1
            continue
        
        # Check 2: All aligned genomes must be mappable (no multi-mapping)
        if len(multi_mapped_genomes) > 0:
            stats["multi_mapped"] += 1
            continue
        
        # Passed both checks
        stats["passed"] += 1
        passed_snps.append(snp)
    
    # Log detailed statistics
    logger.info(f"Accessibility filter results:")
    logger.info(f"  Input: {len(snps)} SNPs")
    logger.info(f"")
    logger.info(f"  Per-genome alignment statistics:")
    for i in range(n_genomes):
        genome_name = index_paths[i].stem if hasattr(index_paths[i], 'stem') else f"genome_{i}"
        aligned_pct = genome_stats["aligned_count"][i] / len(snps) * 100
        mappable_pct = genome_stats["mappable_count"][i] / len(snps) * 100
        multi_pct = genome_stats["multi_mapped_count"][i] / len(snps) * 100
        logger.info(f"    Genome {i+1} ({genome_name}):")
        logger.info(f"      - Aligned: {genome_stats['aligned_count'][i]} ({aligned_pct:.1f}%)")
        logger.info(f"      - Mappable (unique): {genome_stats['mappable_count'][i]} ({mappable_pct:.1f}%)")
        logger.info(f"      - Multi-mapped: {genome_stats['multi_mapped_count'][i]} ({multi_pct:.1f}%)")
    logger.info(f"")
    logger.info(f"  Filtering summary:")
    logger.info(f"    - Insufficient genomes (< {min_genomes}): {stats['insufficient_genomes']}")
    logger.info(f"    - Multi-mapped in some genome: {stats['multi_mapped']}")
    logger.info(f"    - Passed: {stats['passed']} ({stats['passed']/len(snps)*100:.1f}%)")
    
    return Ok(passed_snps)


# =============================================================================
# Biophysical Filter
# =============================================================================

def calculate_probe_stats(
    snp: SNP,
    probe_sequence: str = None,
) -> Dict[str, float]:
    """
    Calculate biophysical statistics for a SNP probe sequence.
    
    Args:
        snp: SNP object
        probe_sequence: Complete probe sequence (if None, uses snp.ref as fallback)
        
    Returns:
        Dictionary of calculated statistics
    """
    if probe_sequence is None:
        probe_seq = snp.ref  # Fallback for compatibility
    else:
        probe_seq = probe_sequence
    
    return {
        "gc": calculate_gc(probe_seq),
        "tm": calculate_tm(probe_seq),
        "complexity": calculate_complexity(probe_seq),
        "hairpin": calculate_hairpin_score(probe_seq),
    }


def _calculate_batch_stats(snps: List[SNP]) -> List[Dict[str, float]]:
    """
    Calculate biophysical stats for a batch of SNPs (for parallel processing).
    
    This function is designed to be called by ProcessPoolExecutor.
    
    Args:
        snps: List of SNPs to process
        
    Returns:
        List of statistics dictionaries
    """
    return [calculate_probe_stats(snp) for snp in snps]


def filter_biophysical(
    snps: List[SNP],
    thresholds: BiophysicalThresholds,
    threads: int = 1,
) -> Result[Tuple[List[SNP], Dict[str, Any]], str]:
    """
    Filter SNPs based on biophysical properties (parallelized).
    
    Calculates GC content, melting temperature, sequence complexity,
    and optionally hairpin/dimer scores. Filters based on thresholds.
    Uses multiprocessing for faster computation on large datasets.
    
    Args:
        snps: Input SNP list
        thresholds: Biophysical threshold configuration
        threads: Number of parallel workers (default: 1)
        
    Returns:
        Result containing (filtered SNPs, statistics dict)
    """
    if not snps:
        return Ok(([], {}))
    
    logger.info(f"Running biophysical filter on {len(snps)} SNPs (threads={threads})")
    logger.info(f"Thresholds: GC={thresholds.gc_min}-{thresholds.gc_max}%, "
                f"Tm={thresholds.tm_min}-{thresholds.tm_max}°C, "
                f"Complexity≤{thresholds.complexity_max}")
    
    filter_stats = {
        "gc_failed": 0,
        "tm_failed": 0,
        "complexity_failed": 0,
        "hairpin_failed": 0,
        "dimer_failed": 0,
    }
    
    # Parallel computation of biophysical properties
    if threads > 1:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        
        # Split into batches to reduce process creation overhead
        batch_size = max(100, len(snps) // (threads * 4))
        batches = [snps[i:i + batch_size] for i in range(0, len(snps), batch_size)]
        
        logger.info(f"Processing {len(batches)} batches with {threads} workers")
        
        all_stats = []
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(_calculate_batch_stats, batch): batch
                for batch in batches
            }
            
            for future in as_completed(futures):
                try:
                    batch_stats = future.result()
                    all_stats.extend(batch_stats)
                except Exception as e:
                    logger.error(f"Batch processing failed: {e}")
                    return Err(f"Biophysical calculation failed: {e}")
    else:
        # Serial computation
        all_stats = [calculate_probe_stats(snp) for snp in snps]
    
    # Filter based on thresholds
    passed_snps = []
    for snp, stats in zip(snps, all_stats):
        # GC check
        if not (thresholds.gc_min <= stats["gc"] <= thresholds.gc_max):
            filter_stats["gc_failed"] += 1
            continue
        
        # Tm check
        if not (thresholds.tm_min <= stats["tm"] <= thresholds.tm_max):
            filter_stats["tm_failed"] += 1
            continue
        
        # Complexity check
        if stats["complexity"] > thresholds.complexity_max:
            filter_stats["complexity_failed"] += 1
            continue
        
        # Hairpin check (optional)
        if thresholds.hairpin_max is not None:
            if stats["hairpin"] > thresholds.hairpin_max:
                filter_stats["hairpin_failed"] += 1
                continue
        
        passed_snps.append(snp)
    
    removed = len(snps) - len(passed_snps)
    logger.info(f"Biophysical filter: removed {removed} SNPs")
    logger.info(f"  GC failed: {filter_stats['gc_failed']}")
    logger.info(f"  Tm failed: {filter_stats['tm_failed']}")
    logger.info(f"  Complexity failed: {filter_stats['complexity_failed']}")
    if thresholds.hairpin_max is not None:
        logger.info(f"  Hairpin failed: {filter_stats['hairpin_failed']}")
    
    return Ok((passed_snps, filter_stats))


# =============================================================================
# Taxonomic Filter (ngsLCA)
# =============================================================================

def run_taxonomic_mapping(
    fasta_path: Path,
    index_path: Path,
    output_prefix: Path,
    keep_hits: int = 100,
    threads: int = 1,
) -> Result[Path, str]:
    """
    Run Bowtie2 alignment for taxonomic assignment.
    
    Args:
        fasta_path: Input FASTA file
        index_path: Bowtie2 index prefix
        output_prefix: Output file prefix
        keep_hits: Number of hits to report (for ngsLCA, default: 100)
        threads: Number of threads
        
    Returns:
        Result containing path to output SAM file
    """
    db_name = Path(index_path).stem
    sam_path = Path(str(output_prefix) + f".{db_name}.sam")
    
    cmd = [
        "bowtie2",
        "-f",
        "-k", str(keep_hits),
        "--threads", str(threads),
        "-x", str(index_path),
        "-U", str(fasta_path),
        "--no-unal",
        "-S", str(sam_path),
    ]
    
    logger.debug(f"Running Bowtie2 for taxonomy: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        return Ok(sam_path)
    except subprocess.CalledProcessError as e:
        return Err(f"Bowtie2 mapping failed: {e.stderr}")
    except FileNotFoundError:
        return Err("Bowtie2 not found. Please install Bowtie2 and ensure it's in PATH.")


def process_sam_to_bam(
    sam_path: Path,
    threads: int = 1,
) -> Result[Path, str]:
    """
    Convert SAM to sorted BAM file.
    
    Args:
        sam_path: Input SAM file
        threads: Number of threads
        
    Returns:
        Result containing path to output BAM file
    """
    bam_path = sam_path.with_suffix(".bam")
    
    # Check if SAM has alignments
    check_cmd = f"samtools view {sam_path} | head -n 1"
    result = subprocess.run(check_cmd, shell=True, capture_output=True, text=True)
    
    if not result.stdout.strip():
        logger.warning(f"No alignments found in {sam_path}")
        return Err("No alignments found")
    
    # Convert to BAM
    cmd = [
        "samtools", "view",
        "-@", str(threads),
        "-b",
        "-o", str(bam_path),
        str(sam_path),
    ]
    
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        sam_path.unlink()  # Remove SAM file
        return Ok(bam_path)
    except subprocess.CalledProcessError as e:
        return Err(f"SAM to BAM conversion failed: {e.stderr}")
    except FileNotFoundError:
        return Err("samtools not found. Please install samtools and ensure it's in PATH.")


def merge_and_sort_bams(
    bam_files: List[Path],
    output_path: Path,
    threads: int = 1,
) -> Result[Path, str]:
    """
    Merge multiple BAM files and sort by name (required for ngsLCA).
    
    Args:
        bam_files: List of BAM files to merge
        output_path: Output sorted BAM path
        threads: Number of threads
        
    Returns:
        Result containing path to sorted BAM file
    """
    if len(bam_files) == 0:
        return Err("No BAM files to merge")
    
    if len(bam_files) == 1:
        # Single file, just sort
        sorted_bam = output_path
        sort_cmd = [
            "samtools", "sort",
            "-n",  # Sort by read name (required for ngsLCA)
            "-@", str(threads),
            "-O", "bam",
            "-o", str(sorted_bam),
            str(bam_files[0]),
        ]
    else:
        # Multiple files, merge then sort
        merged_bam = output_path.with_suffix(".merged.bam")
        merge_cmd = [
            "samtools", "merge",
            "-n",  # Assume input sorted by name
            "-@", str(threads),
            "-o", str(merged_bam),
        ] + [str(f) for f in bam_files]
        
        try:
            subprocess.run(merge_cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            return Err(f"BAM merge failed: {e.stderr}")
        
        sorted_bam = output_path
        sort_cmd = [
            "samtools", "sort",
            "-n",
            "-@", str(threads),
            "-O", "bam",
            "-o", str(sorted_bam),
            str(merged_bam),
        ]
        merged_bam.unlink()  # Clean up merged file
    
    try:
        subprocess.run(sort_cmd, capture_output=True, text=True, check=True)
        return Ok(sorted_bam)
    except subprocess.CalledProcessError as e:
        return Err(f"BAM sort failed: {e.stderr}")


def run_ngslca(
    bam_path: Path,
    names_dmp: Path,
    nodes_dmp: Path,
    acc2tax: Path,
    min_edit: int = 0,
    max_edit: int = 2,
    output_prefix: Path = None,
) -> Result[Path, str]:
    """
    Run ngsLCA taxonomic classification.
    
    Args:
        bam_path: Input sorted BAM file
        names_dmp: NCBI taxonomy names.dmp file
        nodes_dmp: NCBI taxonomy nodes.dmp file
        acc2tax: Accession to taxid mapping file
        min_edit: Minimum edit distance (default: 0)
        max_edit: Maximum edit distance (default: 2)
        output_prefix: Output prefix (default: bam_path.stem)
        
    Returns:
        Result containing path to .lca output file
    """
    if output_prefix is None:
        output_prefix = bam_path.parent / bam_path.stem
    
    lca_output = Path(str(output_prefix) + f".min{min_edit}max{max_edit}.lca")
    
    cmd = [
        "ngsLCA",
        "-names", str(names_dmp),
        "-nodes", str(nodes_dmp),
        "-acc2tax", str(acc2tax),
        "-editdistmin", str(min_edit),
        "-editdistmax", str(max_edit),
        "-bam", str(bam_path),
        "-outnames", str(output_prefix) + f".min{min_edit}max{max_edit}",
    ]
    
    logger.debug(f"Running ngsLCA: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return Ok(lca_output)
    except subprocess.CalledProcessError as e:
        return Err(f"ngsLCA failed: {e.stderr}")
    except FileNotFoundError:
        return Err("ngsLCA not found. Please install ngsLCA and ensure it's in PATH.")


def parse_lca_results(
    lca_path: Path,
    target_taxids: List[int],
) -> Result[Set[str], str]:
    """
    Parse ngsLCA output to extract sequences assigned to target taxa.
    
    Args:
        lca_path: Path to .lca file
        target_taxids: List of target taxonomy IDs
        
    Returns:
        Result containing set of sequence IDs assigned to target taxa
    """
    import re
    
    assigned_ids: Set[str] = set()
    
    try:
        with open(lca_path) as f:
            for line in f:
                # Check if line contains any target taxid
                for taxid in target_taxids:
                    if re.search(rf'\s+{taxid}:', line):
                        # Extract sequence ID (before last 3 colons)
                        parts = line.split('\t')[0].split(':')
                        seq_id = ':'.join(parts[:-3])
                        assigned_ids.add(seq_id.strip())
                        break
        
        return Ok(assigned_ids)
    except Exception as e:
        return Err(f"Failed to parse LCA results: {e}")


def filter_taxonomy(
    snps: List[SNP],
    index_paths: List[Path],
    names_dmp: Path,
    nodes_dmp: Path,
    acc2tax: Path,
    target_taxids: List[int],
    fasta_path: Optional[Path] = None,
    min_edit: int = 0,
    max_edit: int = 2,
    keep_hits: int = 100,
    threads: int = 1,
) -> Result[List[SNP], str]:
    """
    Filter SNPs based on taxonomic assignment using ngsLCA.
    
    Workflow:
    1. Map probe sequences to reference database(s) with Bowtie2
    2. Convert SAM to BAM and merge if multiple databases
    3. Run ngsLCA for taxonomic classification
    4. Extract SNPs assigned to target taxa
    
    Args:
        snps: Input SNP list
        index_paths: List of Bowtie2 index paths
        names_dmp: NCBI taxonomy names.dmp file
        nodes_dmp: NCBI taxonomy nodes.dmp file
        acc2tax: Accession to taxid mapping file
        target_taxids: List of target taxonomy IDs to keep
        fasta_path: Optional pre-computed FASTA file (optimization)
        min_edit: Minimum edit distance for ngsLCA (default: 0)
        max_edit: Maximum edit distance for ngsLCA (default: 2)
        keep_hits: Number of alignments to report per sequence (default: 100)
        threads: Number of threads
        
    Returns:
        Result containing filtered SNP list
    """
    if not snps:
        return Ok([])
    
    logger.info(f"Running taxonomic filter (ngsLCA) on {len(snps)} SNPs")
    logger.info(f"Target taxonomy IDs: {target_taxids}")
    logger.info(f"Using {len(index_paths)} database(s)")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        
        # Use shared FASTA or create temporary one
        if fasta_path is not None:
            # Use shared FASTA (optimization)
            working_fasta = fasta_path
            logger.debug(f"Using shared FASTA: {fasta_path}")
        else:
            # Fallback: create temporary FASTA
            working_fasta = tmpdir_path / "probes.fa"
            sequences = {
                snp.id: probe_sequences[snp.id]
                for snp in snps
            }
            
            write_result = write_fasta(sequences, working_fasta)
            if write_result.is_err():
                return Err(f"Failed to write temp FASTA: {write_result.unwrap_err()}")
        
        # Map to all databases
        bam_files = []
        for idx, index_path in enumerate(index_paths, 1):
            logger.info(f"Mapping to database {idx}/{len(index_paths)}: {index_path.name}")
            output_prefix = tmpdir_path / f"mapping_{idx}"
            
            sam_result = run_taxonomic_mapping(
                working_fasta, index_path, output_prefix, keep_hits, threads
            )
            if sam_result.is_err():
                logger.warning(f"Mapping to {index_path} failed: {sam_result.unwrap_err()}")
                continue
            
            # Convert to BAM
            bam_result = process_sam_to_bam(sam_result.unwrap(), threads)
            if bam_result.is_err():
                logger.warning(f"SAM to BAM conversion failed: {bam_result.unwrap_err()}")
                continue
            
            bam_files.append(bam_result.unwrap())
        
        if not bam_files:
            return Err("No successful alignments from any database")
        
        # Merge and sort BAMs
        logger.info(f"Merging and sorting {len(bam_files)} BAM file(s)")
        sorted_bam_path = tmpdir_path / "merged.sorted.bam"
        merge_result = merge_and_sort_bams(bam_files, sorted_bam_path, threads)
        if merge_result.is_err():
            return Err(merge_result.unwrap_err())
        
        sorted_bam = merge_result.unwrap()
        
        # Run ngsLCA
        logger.info(f"Running ngsLCA (min_edit={min_edit}, max_edit={max_edit})")
        lca_result = run_ngslca(
            sorted_bam, names_dmp, nodes_dmp, acc2tax,
            min_edit, max_edit, tmpdir_path / "taxonomy"
        )
        if lca_result.is_err():
            return Err(lca_result.unwrap_err())
        
        lca_path = lca_result.unwrap()
        
        # Parse LCA results
        logger.info("Parsing taxonomic assignments")
        parse_result = parse_lca_results(lca_path, target_taxids)
        if parse_result.is_err():
            return Err(parse_result.unwrap_err())
        
        assigned_ids = parse_result.unwrap()
    
    # Filter SNPs based on assignments
    filtered_snps = [snp for snp in snps if snp.id in assigned_ids]
    removed = len(snps) - len(filtered_snps)
    
    logger.info(f"Taxonomy filter: {len(assigned_ids)} sequences assigned to target taxa")
    logger.info(f"Taxonomy filter: removed {removed} SNPs, {len(filtered_snps)} remaining")
    
    return Ok(filtered_snps)


# =============================================================================
# Main Filter Pipeline
# =============================================================================

def run_filter(
    input_path: Path,
    reference_path: Path,
    output_prefix: Path,
    filters: List[str],
    bg_db: Optional[str] = None,
    ac_db: Optional[str] = None,
    ac_mode: str = "strict",
    ac_score_diff: int = 10,
    ac_min_genomes: Optional[int] = None,
    tx_db: Optional[str] = None,
    tx_ids: Optional[List[int]] = None,
    names_dmp: Optional[str] = None,
    nodes_dmp: Optional[str] = None,
    acc2tax: Optional[str] = None,
    tx_min_edit: int = 0,
    tx_max_edit: int = 2,
    tx_keep_hits: int = 100,
    gc_range: Tuple[float, float] = (35.0, 65.0),
    tm_range: Tuple[float, float] = (55.0, 75.0),
    max_complexity: float = 2.0,
    max_hairpin: Optional[float] = None,
    max_dimer: Optional[float] = None,
    probe_length: int = 100,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Run multi-stage SNP filtering pipeline.
    
    Main entry point for the filter command. Applies filters in sequence:
    1. Background noise (BG) - if enabled
    2. Accessibility (AC) - if enabled
    3. Taxonomy (TX) - if enabled
    4. Biophysical - if enabled
    
    Args:
        input_path: Input SNP TSV file
        reference_path: Reference genome FASTA
        output_prefix: Output file prefix
        filters: List of filter names to apply
        bg_db: Kraken2 database(s) for background filter (comma-separated)
        ac_db: Bowtie2 index/indices for accessibility filter (comma-separated)
        ac_mode: Accessibility filter mode ('strict' or 'relaxed')
        ac_score_diff: Minimum score difference for relaxed mode (default: 10)
        ac_min_genomes: Minimum number of genomes to align to (default: all)
        tx_db: Bowtie2 index/indices for taxonomic filter (comma-separated)
        tx_ids: Taxonomy IDs to keep (target taxa)
        names_dmp: NCBI taxonomy names.dmp file path
        nodes_dmp: NCBI taxonomy nodes.dmp file path
        acc2tax: Accession to taxid mapping file path
        tx_min_edit: Minimum edit distance for ngsLCA (default: 0)
        tx_max_edit: Maximum edit distance for ngsLCA (default: 2)
        tx_keep_hits: Number of alignments to report (default: 100)
        gc_range: (min, max) GC content range
        tm_range: (min, max) melting temperature range
        max_complexity: Maximum DUST complexity score
        max_hairpin: Maximum hairpin score (optional)
        max_dimer: Maximum dimer score (optional)
        probe_length: Length of probe sequences (for flanking calculation)
        threads: Number of threads
        verbose: Enable verbose logging
        
    Returns:
        Result containing filtering statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Starting SNP filtering from {input_path}")
    logger.info(f"Enabled filters: {filters}")
    
    # Load input SNPs
    snp_df = SNPDataFrame.from_tsv(input_path)
    if snp_df.is_err():
        return Err(f"Failed to load input: {snp_df.unwrap_err()}")
    
    snps = snp_df.unwrap().to_snps()
    initial_count = len(snps)
    logger.info(f"Loaded {initial_count} SNPs")
    
    # Add flanking sequences to SNPs
    logger.info("Adding flanking sequences from reference genome")
    from eprobe.core.fasta import read_fasta
    
    ref_result = read_fasta(reference_path)
    if ref_result.is_err():
        return Err(f"Failed to read reference: {ref_result.unwrap_err()}")
    
    ref_dict = ref_result.unwrap()
    
    # Create probe sequences dictionary (SNP ID -> probe sequence)
    flank_size = (probe_length - 1) // 2
    probe_sequences = {}
    
    for snp in snps:
        if snp.chrom not in ref_dict:
            # Set empty probe for unknown chromosomes
            probe_sequences[snp.id] = snp.ref
            continue
        
        chrom_seq = ref_dict[snp.chrom]
        start_pos = max(0, snp.pos - 1 - flank_size)  # Convert to 0-based
        end_pos = min(len(chrom_seq), snp.pos + flank_size)
        
        # Extract flanking regions
        left_end = snp.pos - 1  # Position before SNP (0-based)
        right_start = snp.pos  # Position after SNP (0-based)
        
        left_flank = chrom_seq[start_pos:left_end] if left_end > start_pos else ""
        right_flank = chrom_seq[right_start:end_pos] if end_pos > right_start else ""
        
        # Build complete probe sequence
        probe_sequences[snp.id] = left_flank + snp.ref + right_flank

    stats = {
        "initial_count": initial_count,
        "filters_applied": [],
    }
    
    # Normalize filter names
    filters_normalized = [f.lower() for f in filters]
    
    # Write shared FASTA file for external filters (BG/AC/TX)
    # Only create if any external filter is enabled
    needs_fasta = any(f in filters_normalized for f in ["bg", "ac", "tx"])
    fasta_path = None
    
    if needs_fasta:
        logger.info("Writing probe sequences to temporary FASTA (shared by external filters)")
        import tempfile
        temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
        fasta_path = Path(temp_fasta.name)
        
        write_result = write_fasta(probe_sequences, fasta_path)
        if write_result.is_err():
            return Err(f"Failed to write shared FASTA: {write_result.unwrap_err()}")
        
        logger.info(f"Shared FASTA created: {fasta_path} ({len(probe_sequences)} sequences)")
    
    # Helper function to update shared FASTA after each filter
    def update_shared_fasta(current_snps: List[SNP], path: Path, sequences: Dict[str, str]) -> Result[None, str]:
        """Update shared FASTA file with current SNPs only."""
        if path is None:
            return Ok(None)
        
        updated_sequences = {
            snp.id: sequences[snp.id]
            for snp in current_snps
        }
        
        write_result = write_fasta(updated_sequences, path)
        if write_result.is_err():
            return Err(f"Failed to update shared FASTA: {write_result.unwrap_err()}")
        
        logger.debug(f"Updated shared FASTA: {len(updated_sequences)} sequences")
        return Ok(None)
    
    try:
        # Apply Background filter (support multiple databases)
        if "bg" in filters_normalized:
            if bg_db is None:
                return Err("BG filter requires --bg_db")
            
            # Parse comma-separated databases
            bg_databases = [Path(db.strip()) for db in bg_db.split(",")]
            logger.info(f"Running background filter with {len(bg_databases)} database(s)")
            
            for idx, db_path in enumerate(bg_databases, 1):
                logger.info(f"Background filter {idx}/{len(bg_databases)}: {db_path.name}")
                result = filter_background_noise(snps, db_path, fasta_path, threads)
                if result.is_err():
                    return Err(result.unwrap_err())
                snps = result.unwrap()
            
            # Update shared FASTA with remaining SNPs for subsequent filters
            update_result = update_shared_fasta(snps, fasta_path, probe_sequences)
            if update_result.is_err():
                return Err(update_result.unwrap_err())
            
            stats["filters_applied"].append("BG")
            stats["bg_remaining"] = len(snps)
        
        # Apply Accessibility filter (support multiple databases)
        if "ac" in filters_normalized:
            if ac_db is None:
                return Err("AC filter requires --ac_db")
            
            # Parse comma-separated databases - all genomes checked together
            ac_databases = [Path(db.strip()) for db in ac_db.split(",")]
            logger.info(f"Running accessibility filter with {len(ac_databases)} genome(s)")
            logger.info(f"Accessibility mode: {ac_mode}, score_diff_threshold: {ac_score_diff}")
            
            # All genomes are checked together for universality + mappability
            result = filter_accessibility(
                snps, ac_databases, threads, 
                mode=ac_mode, 
                score_diff_threshold=ac_score_diff,
                min_genomes=ac_min_genomes,
                fasta_path=fasta_path,
            )
            if result.is_err():
                return Err(result.unwrap_err())
            snps = result.unwrap()
            
            # Update shared FASTA with remaining SNPs for subsequent filters
            update_result = update_shared_fasta(snps, fasta_path, probe_sequences)
            if update_result.is_err():
                return Err(update_result.unwrap_err())
            
            stats["filters_applied"].append("AC")
            stats["ac_remaining"] = len(snps)
        
        # Apply Taxonomy filter
        if "tx" in filters_normalized:
            if tx_db is None:
                return Err("TX filter requires --tx_db")
            if tx_ids is None or len(tx_ids) == 0:
                return Err("TX filter requires --tx_ids (target taxonomy IDs)")
            
            # Parse database paths
            index_paths = [Path(db.strip()) for db in tx_db.split(",")]
            
            # Check for required NCBI taxonomy files
            if names_dmp is None:
                return Err("TX filter requires --tx_names (NCBI names.dmp)")
            if nodes_dmp is None:
                return Err("TX filter requires --tx_nodes (NCBI nodes.dmp)")
            if acc2tax is None:
                return Err("TX filter requires --tx_acc2tax (accession to taxid mapping)")
            
            result = filter_taxonomy(
                snps=snps,
                index_paths=index_paths,
                names_dmp=Path(names_dmp),
                nodes_dmp=Path(nodes_dmp),
                acc2tax=Path(acc2tax),
                target_taxids=tx_ids,
                fasta_path=fasta_path,
                min_edit=tx_min_edit if tx_min_edit is not None else 0,
                max_edit=tx_max_edit if tx_max_edit is not None else 2,
                keep_hits=tx_keep_hits if tx_keep_hits is not None else 100,
                threads=threads,
            )
            if result.is_err():
                return Err(result.unwrap_err())
            
            snps = result.unwrap()
            stats["filters_applied"].append("TX")
            stats["tx_remaining"] = len(snps)
        
        # Apply Biophysical filter
        if "biophysical" in filters_normalized:
            thresholds = BiophysicalThresholds(
            gc_min=gc_range[0],
            gc_max=gc_range[1],
            tm_min=tm_range[0],
            tm_max=tm_range[1],
            complexity_max=max_complexity,
            hairpin_max=max_hairpin,
            dimer_max=max_dimer,
        )
        
            result = filter_biophysical(snps, thresholds, threads)
            if result.is_err():
                return Err(result.unwrap_err())
            
            snps, biophys_stats = result.unwrap()
            stats["filters_applied"].append("biophysical")
            stats["biophysical_remaining"] = len(snps)
            stats["biophysical_details"] = biophys_stats
        
        # Save output
        output_path = Path(str(output_prefix) + ".filtered.tsv")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        filtered_df = SNPDataFrame.from_snps(snps)
        save_result = filtered_df.to_tsv(output_path)
        if save_result.is_err():
            return Err(f"Failed to save output: {save_result.unwrap_err()}")
        
        final_count = len(snps)
        removed = initial_count - final_count
        
        stats["final_count"] = final_count
        stats["total_removed"] = removed
        stats["output_file"] = str(output_path)
        
        logger.info(f"Filtering complete: {final_count} SNPs remaining (removed {removed})")
        logger.info(f"Output saved to {output_path}")
        
        return Ok(stats)
    
    finally:
        # Clean up shared FASTA file
        if fasta_path is not None and fasta_path.exists():
            fasta_path.unlink()
            logger.debug(f"Cleaned up temporary FASTA: {fasta_path}")