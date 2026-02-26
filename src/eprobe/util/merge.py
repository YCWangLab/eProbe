"""
Merge multiple probe sets with sequence-level deduplication.

Combines FASTA files and removes exact duplicate sequences (by content).
Reports which IDs were removed and which were kept as representative.
"""

import logging
from pathlib import Path
from typing import Dict, Any, List, Tuple
from collections import OrderedDict

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def merge_and_dedup(
    fasta_files: List[Path],
    keep: str = "first",
) -> Tuple[Dict[str, str], List[Tuple[str, str]]]:
    """
    Merge FASTA files and remove exact duplicate sequences.
    
    Duplicates are identified by identical sequence content (case-insensitive).
    Among duplicates, one representative is kept and the rest are removed.
    
    Args:
        fasta_files: List of FASTA file paths
        keep: Strategy for selecting representative ("first")
        
    Returns:
        (unique_sequences, removed_list)
        removed_list: [(removed_id, kept_id), ...]
    """
    all_seqs: Dict[str, str] = OrderedDict()
    
    for fpath in fasta_files:
        result = read_fasta(fpath)
        if result.is_err():
            logger.warning(f"Failed to read {fpath}: {result.unwrap_err()}")
            continue
        seqs = result.unwrap()
        for sid, seq in seqs.items():
            if sid in all_seqs:
                sid = f"{sid}_{fpath.stem}"
            all_seqs[sid] = seq
    
    # Group by sequence content (case-insensitive)
    seq_to_ids: Dict[str, List[str]] = OrderedDict()
    for sid, seq in all_seqs.items():
        key = seq.upper()
        if key not in seq_to_ids:
            seq_to_ids[key] = []
        seq_to_ids[key].append(sid)
    
    unique: Dict[str, str] = OrderedDict()
    removed: List[Tuple[str, str]] = []
    
    for seq_upper, ids in seq_to_ids.items():
        kept = ids[0]
        unique[kept] = all_seqs[kept]
        for rid in ids[1:]:
            removed.append((rid, kept))
    
    return unique, removed


def run_merge(
    input_files: List[Path],
    output_prefix: Path,
    keep: str = "first",
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Merge probe FASTA files with sequence-level deduplication.
    
    Args:
        input_files: List of input FASTA paths
        output_prefix: Output prefix
        keep: Which duplicate to keep ("first")
        verbose: Verbose logging
        
    Returns:
        Result with merge statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    if not input_files:
        return Err("No input files provided")
    
    for f in input_files:
        if not f.exists():
            return Err(f"File not found: {f}")
    
    logger.info(f"Merging {len(input_files)} FASTA files")
    
    input_counts: Dict[str, int] = {}
    total_input = 0
    for f in input_files:
        result = read_fasta(f)
        if result.is_ok():
            count = len(result.unwrap())
            input_counts[f.name] = count
            total_input += count
    
    unique, removed = merge_and_dedup(input_files, keep=keep)
    
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fasta_path = Path(str(output_prefix) + ".merged.fa")
    write_result = write_fasta(unique, fasta_path)
    if write_result.is_err():
        return Err(f"Failed to write output: {write_result.unwrap_err()}")
    
    removed_path = Path(str(output_prefix) + ".removed_duplicates.tsv")
    with open(removed_path, 'w') as f:
        f.write("removed_id\tkept_id\n")
        for rid, kid in removed:
            f.write(f"{rid}\t{kid}\n")
    
    summary_path = Path(str(output_prefix) + ".merge_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("Merge Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write("Input files:\n")
        for fname, cnt in input_counts.items():
            f.write(f"  {fname}: {cnt} probes\n")
        f.write(f"\nTotal input: {total_input}\n")
        f.write(f"Unique probes: {len(unique)}\n")
        f.write(f"Duplicates removed: {len(removed)}\n")
    
    stats = {
        "total_input": total_input,
        "unique_count": len(unique),
        "removed_count": len(removed),
        "input_counts": input_counts,
        "fasta_file": str(fasta_path),
        "removed_file": str(removed_path),
        "summary_file": str(summary_path),
    }
    
    return Ok(stats)
