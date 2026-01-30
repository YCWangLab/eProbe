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


def parse_kraken_output(output_path: Path) -> Result[Set[str], str]:
    """
    Parse Kraken2 output and return classified sequence IDs.
    
    Args:
        output_path: Kraken2 output file
        
    Returns:
        Result containing set of classified sequence IDs
    """
    try:
        classified_ids = set()
        
        with open(output_path, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        classification = parts[0]  # 'C' or 'U'
                        seq_id = parts[1]
                        
                        if classification == 'C':  # Classified
                            classified_ids.add(seq_id)
        
        logger.info(f"Parsed Kraken2 output: {len(classified_ids)} classified sequences")
        return Ok(classified_ids)
        
    except Exception as e:
        return Err(f"Failed to parse Kraken output: {e}")


def display_kraken_summary(kraken_output_file: Path, max_examples: int = 10) -> None:
    """
    Display summary statistics from kraken2 main output file.
    
    Args:
        kraken_output_file: Path to kraken2 main output
        max_examples: Maximum number of example sequences to show
    """
    try:
        classified_count = 0
        unclassified_count = 0
        taxon_counts = {}
        classified_examples = []
        unclassified_examples = []
        
        logger.info("=== Kraken2 Analysis Summary ===")
        
        with open(kraken_output_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if not line.strip():
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                    
                status = parts[0]  # 'C' or 'U'
                seq_id = parts[1]
                taxon_id = parts[2] if len(parts) > 2 else "0"
                
                if status == 'C':
                    classified_count += 1
                    if taxon_id in taxon_counts:
                        taxon_counts[taxon_id] += 1
                    else:
                        taxon_counts[taxon_id] = 1
                    
                    if len(classified_examples) < max_examples:
                        seq_len = parts[3] if len(parts) > 3 else "unknown"
                        class_path = parts[4] if len(parts) > 4 else ""
                        classified_examples.append((seq_id, taxon_id, seq_len, class_path))
                        
                elif status == 'U':
                    unclassified_count += 1
                    if len(unclassified_examples) < max_examples:
                        unclassified_examples.append(seq_id)
        
        total = classified_count + unclassified_count
        logger.info(f"Total sequences analyzed: {total}")
        logger.info(f"  - Classified: {classified_count} ({classified_count/total*100:.1f}%)")
        logger.info(f"  - Unclassified: {unclassified_count} ({unclassified_count/total*100:.1f}%)")
        
        if taxon_counts:
            logger.info(f"Top taxa (by sequence count):")
            sorted_taxa = sorted(taxon_counts.items(), key=lambda x: x[1], reverse=True)[:10]
            for taxon_id, count in sorted_taxa:
                logger.info(f"  - Taxon {taxon_id}: {count} sequences")
        
        if classified_examples:
            logger.info(f"Classified examples (first {len(classified_examples)}):")
            for seq_id, taxon_id, seq_len, class_path in classified_examples:
                path_preview = class_path[:50] + "..." if len(class_path) > 50 else class_path
                logger.info(f"  - {seq_id}: taxon={taxon_id}, len={seq_len}, path={path_preview}")
        
        if unclassified_examples:
            logger.info(f"Unclassified examples: {', '.join(unclassified_examples[:5])}")
        
        logger.info("=" * 40)
        
    except Exception as e:
        logger.warning(f"Could not parse kraken output for summary: {e}")


def filter_background_noise(
    snps: List[SNP],
    db_path: Path,
    fasta_path: Optional[Path] = None,
    threads: int = 1,
    output_dir: Optional[Path] = None,
) -> Result[List[SNP], str]:
    """
    Filter SNPs whose probe sequences match background database.
    
    Uses Kraken2 to identify sequences that match contaminant/noise database.
    SNPs with matching sequences are removed.
    
    Args:
        snps: Input SNP list
        db_path: Kraken2 database path
        fasta_path: Optional pre-computed FASTA file (optimization)
        threads: Number of threads
        output_dir: Directory to save kraken2 output files (if None, uses cwd)
        
    Returns:
        Result containing filtered SNP list
    """
    if not snps:
        return Ok([])
    
    logger.info(f"Running background noise filter (Kraken2) on {len(snps)} SNPs")
    
    # Determine output directory for kraken2 files
    if output_dir is None:
        output_dir = Path.cwd()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Always save kraken2 output to the specified output directory
    kraken_output_path = output_dir / "kraken2_main.out"
    
    # Use shared FASTA or create temporary one
    if fasta_path is not None:
        # Use shared FASTA
        kraken_result = run_kraken2(fasta_path, db_path, kraken_output_path, threads)
        if kraken_result.is_err():
            return Err(kraken_result.unwrap_err())
        
        classified_result = parse_kraken_output(kraken_output_path)
        if classified_result.is_err():
            return Err(classified_result.unwrap_err())
        
        classified_ids = classified_result.unwrap()
        
    else:
        # Fallback: create temporary FASTA
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_fasta_path = Path(tmpdir) / "probes.fa"
            
            # Write probe sequences
            sequences = {
                snp.id: probe_sequences[snp.id]
                for snp in snps
            }
            
            write_result = write_fasta(sequences, temp_fasta_path)
            if write_result.is_err():
                return Err(f"Failed to write temp FASTA: {write_result.unwrap_err()}")
            
            # Run Kraken2 - output to permanent location
            kraken_result = run_kraken2(temp_fasta_path, db_path, kraken_output_path, threads)
            if kraken_result.is_err():
                return Err(kraken_result.unwrap_err())
            
            # Parse results
            classified_result = parse_kraken_output(kraken_output_path)
            if classified_result.is_err():
                return Err(classified_result.unwrap_err())
            
            classified_ids = classified_result.unwrap()
    
    # Log saved kraken2 files
    logger.info(f"Kraken2 output files saved to: {output_dir}")
    logger.info(f"  - Main output: {kraken_output_path}")
    logger.info(f"  - Report: {output_dir / 'kraken2.report'}")
    logger.info(f"  - Classified: {output_dir / 'kraken2.classified.fasta'}")
    logger.info(f"  - Unclassified: {output_dir / 'kraken2.unclassified.fasta'}")
    
    # Filter out classified SNPs
    filtered_snps = [snp for snp in snps if snp.id not in classified_ids]
    removed = len(snps) - len(filtered_snps)
    
    logger.info(f"Background filter: removed {removed} SNPs ({len(classified_ids)} matched noise)")
    
    # Generate and display kraken2 analysis summary
    if kraken_output_path.exists():
        display_kraken_summary(kraken_output_path)
    
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
            snp.id: probe_sequences[snp.id]
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
    filtered_snps = [snp for snp in snps if snp.id not in multi_mapping_ids]
    removed = len(snps) - len(filtered_snps)
    
    logger.info(f"Accessibility filter: removed {removed} multi-mapping SNPs")
    
    return Ok(filtered_snps)


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
                result = filter_background_noise(snps, db_path, fasta_path, threads, output_path.parent)
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
            
            # Parse comma-separated databases
            ac_databases = [Path(db.strip()) for db in ac_db.split(",")]
            logger.info(f"Running accessibility filter with {len(ac_databases)} database(s)")
            
            for idx, db_path in enumerate(ac_databases, 1):
                logger.info(f"Accessibility filter {idx}/{len(ac_databases)}: {db_path.name}")
                result = filter_accessibility(snps, db_path, threads)
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