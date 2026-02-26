"""
Simple FASTA tiling for probe generation.

Generates probes by sliding window across input sequences.
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional
from collections import OrderedDict

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def tile_sequences(
    sequences: Dict[str, str],
    probe_length: int = 81,
    step: int = 30,
) -> Dict[str, str]:
    """
    Generate probes by sliding window tiling.
    
    Args:
        sequences: Input sequences {id: seq}
        probe_length: Probe length
        step: Step size between probes
    
    Returns:
        Dict of probe_id -> probe_sequence
    """
    probes: Dict[str, str] = OrderedDict()
    
    for seq_id, seq in sequences.items():
        seq = seq.upper()
        
        if len(seq) < probe_length:
            probes[f"{seq_id}_tile1"] = seq
            continue
        
        tile_idx = 0
        for start in range(0, len(seq) - probe_length + 1, step):
            tile_idx += 1
            probe_seq = seq[start:start + probe_length]
            end = start + probe_length
            probes[f"{seq_id}:{start + 1}-{end}_tile{tile_idx}"] = probe_seq
        
        # Add a final tile reaching the end if not already covered
        last_start = len(seq) - probe_length
        if last_start > 0 and last_start % step != 0:
            tile_idx += 1
            probes[f"{seq_id}:{last_start + 1}-{len(seq)}_tile{tile_idx}"] = seq[last_start:]
    
    return probes


def run_tile(
    input_path: Path,
    output_prefix: Path,
    probe_length: int = 81,
    step: int = 30,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Tile FASTA sequences into fixed-length probes.
    
    Args:
        input_path: Input FASTA file
        output_prefix: Output prefix
        probe_length: Probe length (default: 81)
        step: Step size (default: 30)
        verbose: Verbose logging
    
    Returns:
        Result with tiling statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Tiling {input_path}: length={probe_length}, step={step}")
    
    fasta_result = read_fasta(input_path)
    if fasta_result.is_err():
        return Err(f"Failed to read input: {fasta_result.unwrap_err()}")
    
    sequences = fasta_result.unwrap()
    logger.info(f"Loaded {len(sequences)} sequences")
    
    probes = tile_sequences(sequences, probe_length, step)
    
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fasta_path = Path(str(output_prefix) + ".tiled.fa")
    write_result = write_fasta(probes, fasta_path)
    if write_result.is_err():
        return Err(f"Failed to write output: {write_result.unwrap_err()}")
    
    stats = {
        "input_count": len(sequences),
        "probe_count": len(probes),
        "probe_length": probe_length,
        "step": step,
        "fasta_file": str(fasta_path),
    }
    
    return Ok(stats)
