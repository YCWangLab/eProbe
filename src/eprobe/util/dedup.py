"""
Probe deduplication.

Removes duplicate or highly similar sequences from probe sets.
Supports exact match and similarity-based clustering (CD-HIT).
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Dict, Any, List
import random

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def deduplicate_exact(
    sequences: Dict[str, str],
    keep: str = "first",
) -> Dict[str, str]:
    """
    Remove exact duplicate sequences.
    
    Args:
        sequences: Input sequences
        keep: Which to keep (first, longest, random)
        
    Returns:
        Deduplicated sequences
    """
    # Group by sequence
    seq_to_ids: Dict[str, List[str]] = {}
    
    for seq_id, seq in sequences.items():
        seq_upper = seq.upper()
        if seq_upper not in seq_to_ids:
            seq_to_ids[seq_upper] = []
        seq_to_ids[seq_upper].append(seq_id)
    
    # Select representative for each unique sequence
    result = {}
    
    for seq, ids in seq_to_ids.items():
        if keep == "first":
            selected_id = ids[0]
        elif keep == "longest":
            # This doesn't make sense for exact duplicates, use first
            selected_id = ids[0]
        elif keep == "random":
            selected_id = random.choice(ids)
        else:
            selected_id = ids[0]
        
        result[selected_id] = seq
    
    return result


def run_cdhit(
    fasta_path: Path,
    output_path: Path,
    identity: float = 0.95,
    coverage: float = 0.9,
    threads: int = 1,
) -> Result[Path, str]:
    """
    Run CD-HIT for sequence clustering.
    
    Args:
        fasta_path: Input FASTA file
        output_path: Output FASTA file
        identity: Sequence identity threshold
        coverage: Coverage threshold
        threads: Number of threads
        
    Returns:
        Result containing output path
    """
    cmd = [
        "cd-hit",
        "-i", str(fasta_path),
        "-o", str(output_path),
        "-c", str(identity),
        "-aS", str(coverage),
        "-T", str(threads),
        "-M", "0",  # Use all available memory
        "-d", "0",  # Full sequence description
    ]
    
    logger.debug(f"Running CD-HIT: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        return Ok(output_path)
    except subprocess.CalledProcessError as e:
        return Err(f"CD-HIT failed: {e.stderr}")
    except FileNotFoundError:
        return Err("CD-HIT not found. Please install CD-HIT and ensure it's in PATH.")


def parse_cdhit_clusters(cluster_path: Path) -> Dict[str, List[str]]:
    """
    Parse CD-HIT cluster file.
    
    Args:
        cluster_path: Path to .clstr file
        
    Returns:
        Dictionary of representative_id -> [member_ids]
    """
    clusters: Dict[str, List[str]] = {}
    current_cluster: List[str] = []
    representative: Optional[str] = None
    
    with open(cluster_path) as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('>Cluster'):
                # Save previous cluster
                if representative and current_cluster:
                    clusters[representative] = current_cluster
                current_cluster = []
                representative = None
            else:
                # Parse member line
                # Format: 0       81aa, >SeqID... *
                parts = line.split('>')
                if len(parts) >= 2:
                    seq_id = parts[1].split('...')[0]
                    current_cluster.append(seq_id)
                    
                    if line.endswith('*'):
                        representative = seq_id
        
        # Save last cluster
        if representative and current_cluster:
            clusters[representative] = current_cluster
    
    return clusters


def run_dedup(
    input_path: Path,
    output_path: Path,
    method: str = "simple",
    identity: float = 0.95,
    coverage: float = 0.9,
    keep: str = "first",
    clusters_path: Optional[Path] = None,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Deduplicate probe sequences.
    
    Main entry point for the dedup command.
    
    Args:
        input_path: Input FASTA file
        output_path: Output FASTA file
        method: Deduplication method (simple or cdhit)
        identity: Identity threshold for CD-HIT
        coverage: Coverage threshold for CD-HIT
        keep: Which duplicate to keep
        clusters_path: Optional path to save cluster information
        threads: Number of threads for CD-HIT
        verbose: Enable verbose logging
        
    Returns:
        Result containing deduplication statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Deduplicating {input_path}")
    logger.info(f"Method: {method}")
    
    # Read input
    fasta_result = read_fasta(input_path)
    if fasta_result.is_err():
        return Err(f"Failed to read input: {fasta_result.unwrap_err()}")
    
    sequences = fasta_result.unwrap()
    input_count = len(sequences)
    logger.info(f"Loaded {input_count} sequences")
    
    if method == "simple":
        # Simple exact deduplication
        deduped = deduplicate_exact(sequences, keep)
        
        # Write output
        output_path.parent.mkdir(parents=True, exist_ok=True)
        write_result = write_fasta(deduped, output_path)
        if write_result.is_err():
            return Err(f"Failed to write output: {write_result.unwrap_err()}")
        
        output_count = len(deduped)
        
    elif method == "cdhit":
        # CD-HIT clustering
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        cdhit_result = run_cdhit(
            input_path,
            output_path,
            identity=identity,
            coverage=coverage,
            threads=threads,
        )
        
        if cdhit_result.is_err():
            return Err(cdhit_result.unwrap_err())
        
        # Read output to count
        output_result = read_fasta(output_path)
        if output_result.is_err():
            return Err(f"Failed to read CD-HIT output: {output_result.unwrap_err()}")
        
        output_count = len(output_result.unwrap())
        
        # Parse clusters if requested
        if clusters_path:
            cluster_file = Path(str(output_path) + ".clstr")
            if cluster_file.exists():
                clusters = parse_cdhit_clusters(cluster_file)
                
                with open(clusters_path, 'w') as f:
                    f.write("representative\tmembers\tcluster_size\n")
                    for rep, members in clusters.items():
                        f.write(f"{rep}\t{','.join(members)}\t{len(members)}\n")
    
    else:
        return Err(f"Unknown method: {method}")
    
    removed = input_count - output_count
    logger.info(f"Removed {removed} duplicates, {output_count} remaining")
    
    stats = {
        "input_count": input_count,
        "output_count": output_count,
        "removed_count": removed,
        "method": method,
        "output_file": str(output_path),
    }
    
    if method == "cdhit":
        stats["identity"] = identity
        stats["coverage"] = coverage
    
    return Ok(stats)
