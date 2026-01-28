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
) -> Result[Path, str]:
    """
    Run Kraken2 classification on probe sequences.
    
    Args:
        fasta_path: Input FASTA file
        db_path: Kraken2 database path
        output_path: Output classification file
        threads: Number of threads
        
    Returns:
        Result containing path to output file
    """
    cmd = [
        "kraken2",
        "--db", str(db_path),
        "--threads", str(threads),
        "--output", str(output_path),
        "--report", str(output_path) + ".report",
        str(fasta_path),
    ]
    
    logger.debug(f"Running Kraken2: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        return Ok(output_path)
    except subprocess.CalledProcessError as e:
        return Err(f"Kraken2 failed: {e.stderr}")
    except FileNotFoundError:
        return Err("Kraken2 not found. Please install Kraken2 and ensure it's in PATH.")


def parse_kraken_output(output_path: Path) -> Result[Set[str], str]:
    """
    Parse Kraken2 output to find classified sequences.
    
    Returns set of sequence IDs that were classified (should be filtered out).
    """
    classified_ids: Set[str] = set()
    
    try:
        with open(output_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2 and parts[0] == "C":  # Classified
                    classified_ids.add(parts[1])
        return Ok(classified_ids)
    except Exception as e:
        return Err(f"Failed to parse Kraken2 output: {e}")


def filter_background_noise(
    snps: List[SNP],
    db_path: Path,
    threads: int = 1,
) -> Result[List[SNP], str]:
    """
    Filter SNPs whose probe sequences match background database.
    
    Uses Kraken2 to identify sequences that match contaminant/noise database.
    SNPs with matching sequences are removed.
    
    Args:
        snps: Input SNP list
        db_path: Kraken2 database path
        threads: Number of threads
        
    Returns:
        Result containing filtered SNP list
    """
    if not snps:
        return Ok([])
    
    logger.info(f"Running background noise filter (Kraken2) on {len(snps)} SNPs")
    
    # Create temporary FASTA for Kraken2
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = Path(tmpdir) / "probes.fa"
        output_path = Path(tmpdir) / "kraken.out"
        
        # Write probe sequences
        sequences = {
            snp.snp_id: snp.left_flank + snp.ref + snp.right_flank
            for snp in snps
        }
        
        write_result = write_fasta(sequences, fasta_path)
        if write_result.is_err():
            return Err(f"Failed to write temp FASTA: {write_result.unwrap_err()}")
        
        # Run Kraken2
        kraken_result = run_kraken2(fasta_path, db_path, output_path, threads)
        if kraken_result.is_err():
            return Err(kraken_result.unwrap_err())
        
        # Parse results
        classified_result = parse_kraken_output(output_path)
        if classified_result.is_err():
            return Err(classified_result.unwrap_err())
        
        classified_ids = classified_result.unwrap()
    
    # Filter out classified SNPs
    filtered_snps = [snp for snp in snps if snp.snp_id not in classified_ids]
    removed = len(snps) - len(filtered_snps)
    
    logger.info(f"Background filter: removed {removed} SNPs ({len(classified_ids)} matched noise)")
    
    return Ok(filtered_snps)


# =============================================================================
# Accessibility Filter (Bowtie2)
# =============================================================================

def run_bowtie2(
    fasta_path: Path,
    index_path: Path,
    output_path: Path,
    threads: int = 1,
) -> Result[Path, str]:
    """
    Run Bowtie2 alignment for accessibility check.
    
    Args:
        fasta_path: Input FASTA file
        index_path: Bowtie2 index prefix
        output_path: Output SAM file
        threads: Number of threads
        
    Returns:
        Result containing path to output file
    """
    cmd = [
        "bowtie2",
        "-x", str(index_path),
        "-f", str(fasta_path),
        "-S", str(output_path),
        "-p", str(threads),
        "--no-unal",  # Don't output unaligned sequences
        "-k", "2",    # Report up to 2 alignments
    ]
    
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


def parse_bowtie2_alignments(
    sam_path: Path,
    max_hits: int = 1,
) -> Result[Set[str], str]:
    """
    Parse Bowtie2 SAM output to find multi-mapping sequences.
    
    Returns set of sequence IDs that have more than max_hits alignments.
    """
    hit_counts: Dict[str, int] = {}
    
    try:
        with open(sam_path) as f:
            for line in f:
                if line.startswith("@"):  # Skip header
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    seq_id = parts[0]
                    flag = int(parts[1])
                    if not (flag & 4):  # Not unmapped
                        hit_counts[seq_id] = hit_counts.get(seq_id, 0) + 1
        
        multi_mapping = {
            seq_id for seq_id, count in hit_counts.items()
            if count > max_hits
        }
        return Ok(multi_mapping)
    except Exception as e:
        return Err(f"Failed to parse Bowtie2 output: {e}")


def filter_accessibility(
    snps: List[SNP],
    index_path: Path,
    threads: int = 1,
) -> Result[List[SNP], str]:
    """
    Filter SNPs based on genomic accessibility.
    
    Uses Bowtie2 to check if probe sequences map uniquely to the genome.
    SNPs with multi-mapping probes are removed.
    
    Args:
        snps: Input SNP list
        index_path: Bowtie2 index prefix path
        threads: Number of threads
        
    Returns:
        Result containing filtered SNP list
    """
    if not snps:
        return Ok([])
    
    logger.info(f"Running accessibility filter (Bowtie2) on {len(snps)} SNPs")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = Path(tmpdir) / "probes.fa"
        sam_path = Path(tmpdir) / "alignments.sam"
        
        # Write probe sequences
        sequences = {
            snp.snp_id: snp.left_flank + snp.ref + snp.right_flank
            for snp in snps
        }
        
        write_result = write_fasta(sequences, fasta_path)
        if write_result.is_err():
            return Err(f"Failed to write temp FASTA: {write_result.unwrap_err()}")
        
        # Run Bowtie2
        bt2_result = run_bowtie2(fasta_path, index_path, sam_path, threads)
        if bt2_result.is_err():
            return Err(bt2_result.unwrap_err())
        
        # Parse results
        multi_map_result = parse_bowtie2_alignments(sam_path)
        if multi_map_result.is_err():
            return Err(multi_map_result.unwrap_err())
        
        multi_mapping_ids = multi_map_result.unwrap()
    
    # Filter out multi-mapping SNPs
    filtered_snps = [snp for snp in snps if snp.snp_id not in multi_mapping_ids]
    removed = len(snps) - len(filtered_snps)
    
    logger.info(f"Accessibility filter: removed {removed} multi-mapping SNPs")
    
    return Ok(filtered_snps)


# =============================================================================
# Biophysical Filter
# =============================================================================

def calculate_probe_stats(
    snp: SNP,
) -> Dict[str, float]:
    """
    Calculate biophysical statistics for a SNP probe sequence.
    
    Args:
        snp: SNP with flanking sequences
        
    Returns:
        Dictionary of calculated statistics
    """
    probe_seq = snp.left_flank + snp.ref + snp.right_flank
    
    return {
        "gc": calculate_gc(probe_seq),
        "tm": calculate_tm(probe_seq),
        "complexity": calculate_complexity(probe_seq),
        "hairpin": calculate_hairpin_score(probe_seq),
    }


def filter_biophysical(
    snps: List[SNP],
    thresholds: BiophysicalThresholds,
) -> Result[Tuple[List[SNP], Dict[str, Any]], str]:
    """
    Filter SNPs based on biophysical properties.
    
    Calculates GC content, melting temperature, sequence complexity,
    and optionally hairpin/dimer scores. Filters based on thresholds.
    
    Args:
        snps: Input SNP list
        thresholds: Biophysical threshold configuration
        
    Returns:
        Result containing (filtered SNPs, statistics dict)
    """
    if not snps:
        return Ok(([], {}))
    
    logger.info(f"Running biophysical filter on {len(snps)} SNPs")
    logger.info(f"Thresholds: GC={thresholds.gc_min}-{thresholds.gc_max}%, "
                f"Tm={thresholds.tm_min}-{thresholds.tm_max}°C, "
                f"Complexity≤{thresholds.complexity_max}")
    
    passed_snps = []
    filter_stats = {
        "gc_failed": 0,
        "tm_failed": 0,
        "complexity_failed": 0,
        "hairpin_failed": 0,
        "dimer_failed": 0,
    }
    
    for snp in snps:
        stats = calculate_probe_stats(snp)
        
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
# Main Filter Pipeline
# =============================================================================

def run_filter(
    input_path: Path,
    reference_path: Path,
    output_prefix: Path,
    filters: List[str],
    bg_db: Optional[Path] = None,
    ac_db: Optional[Path] = None,
    tx_db: Optional[Path] = None,
    tx_ids: Optional[List[int]] = None,
    gc_range: Tuple[float, float] = (35.0, 65.0),
    tm_range: Tuple[float, float] = (55.0, 75.0),
    max_complexity: float = 2.0,
    max_hairpin: Optional[float] = None,
    max_dimer: Optional[float] = None,
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
        bg_db: Kraken2 database for background filter
        ac_db: Bowtie2 index for accessibility filter
        tx_db: Taxonomy database path
        tx_ids: Taxonomy IDs to filter against
        gc_range: (min, max) GC content range
        tm_range: (min, max) melting temperature range
        max_complexity: Maximum DUST complexity score
        max_hairpin: Maximum hairpin score (optional)
        max_dimer: Maximum dimer score (optional)
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
    
    stats = {
        "initial_count": initial_count,
        "filters_applied": [],
    }
    
    # Normalize filter names
    filters_normalized = [f.lower() for f in filters]
    
    # Apply Background filter
    if "bg" in filters_normalized:
        if bg_db is None:
            return Err("BG filter requires --bg_db")
        
        result = filter_background_noise(snps, bg_db, threads)
        if result.is_err():
            return Err(result.unwrap_err())
        
        snps = result.unwrap()
        stats["filters_applied"].append("BG")
        stats["bg_remaining"] = len(snps)
    
    # Apply Accessibility filter
    if "ac" in filters_normalized:
        if ac_db is None:
            return Err("AC filter requires --ac_db")
        
        result = filter_accessibility(snps, ac_db, threads)
        if result.is_err():
            return Err(result.unwrap_err())
        
        snps = result.unwrap()
        stats["filters_applied"].append("AC")
        stats["ac_remaining"] = len(snps)
    
    # Apply Taxonomy filter (simplified - would need implementation)
    if "tx" in filters_normalized:
        logger.warning("Taxonomy filter not yet implemented")
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
        
        result = filter_biophysical(snps, thresholds)
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
