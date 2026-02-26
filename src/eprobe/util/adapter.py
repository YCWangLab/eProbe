"""
Add adapter sequences to probes.

Batch append adapter sequences to the 5' or 3' end of all probes
in a FASTA file.
"""

import logging
from pathlib import Path
from typing import Dict, Any
from collections import OrderedDict

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def add_adapters(
    sequences: Dict[str, str],
    adapter_seq: str,
    end: str = "5prime",
) -> Dict[str, str]:
    """
    Add adapter sequences to probes.
    
    Args:
        sequences: Input probe sequences
        adapter_seq: Adapter sequence to add
        end: "5prime" or "3prime"
    
    Returns:
        Adapted sequences
    """
    adapted: Dict[str, str] = OrderedDict()
    adapter = adapter_seq.upper()
    
    for sid, seq in sequences.items():
        if end == "5prime":
            adapted[sid] = adapter + seq
        else:
            adapted[sid] = seq + adapter
    
    return adapted


def run_adapter(
    input_path: Path,
    output_prefix: Path,
    adapter_seq: str,
    end: str = "5prime",
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Add adapter sequences to a probe FASTA.
    
    Args:
        input_path: Input FASTA file
        output_prefix: Output prefix
        adapter_seq: Adapter sequence
        end: "5prime" or "3prime"
        verbose: Verbose logging
        
    Returns:
        Result with statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    if not adapter_seq:
        return Err("Adapter sequence cannot be empty")
    
    valid_bases = set("ATCGN")
    if not set(adapter_seq.upper()).issubset(valid_bases):
        return Err(f"Adapter contains invalid bases: {adapter_seq}")
    
    logger.info(f"Adding {end} adapter ({len(adapter_seq)}bp): {adapter_seq}")
    
    fasta_result = read_fasta(input_path)
    if fasta_result.is_err():
        return Err(f"Failed to read input: {fasta_result.unwrap_err()}")
    
    sequences = fasta_result.unwrap()
    adapted = add_adapters(sequences, adapter_seq, end)
    
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fasta_path = Path(str(output_prefix) + ".adapted.fa")
    write_result = write_fasta(adapted, fasta_path)
    if write_result.is_err():
        return Err(f"Failed to write output: {write_result.unwrap_err()}")
    
    original_len = len(next(iter(sequences.values()))) if sequences else 0
    new_len = len(next(iter(adapted.values()))) if adapted else 0
    
    stats = {
        "probe_count": len(adapted),
        "adapter_seq": adapter_seq,
        "adapter_length": len(adapter_seq),
        "end": end,
        "original_length": original_len,
        "new_length": new_len,
        "fasta_file": str(fasta_path),
    }
    
    return Ok(stats)
