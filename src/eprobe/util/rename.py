"""
Batch rename probe IDs.

Two modes:
  prefix: Rename all to {prefix}_0001, {prefix}_0002, ...
  id_map: Rename according to old_id -> new_id mapping file (TSV)

Always outputs a correspondence map (new_id -> old_id).
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional, Tuple, List
from collections import OrderedDict

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def rename_by_prefix(
    sequences: Dict[str, str],
    prefix: str,
    start_index: int = 1,
    zero_pad: int = 0,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Rename sequences with prefix + sequential index.
    
    Returns:
        (renamed_sequences, id_map: {new_id: old_id})
    """
    renamed: Dict[str, str] = OrderedDict()
    id_map: Dict[str, str] = OrderedDict()
    
    if zero_pad == 0:
        zero_pad = max(len(str(len(sequences) + start_index - 1)), 4)
    
    for idx, (old_id, seq) in enumerate(sequences.items(), start=start_index):
        new_id = f"{prefix}_{str(idx).zfill(zero_pad)}"
        renamed[new_id] = seq
        id_map[new_id] = old_id
    
    return renamed, id_map


def rename_by_map(
    sequences: Dict[str, str],
    mapping: Dict[str, str],
) -> Tuple[Dict[str, str], Dict[str, str], List[str]]:
    """
    Rename sequences according to a provided mapping.
    
    Args:
        sequences: Input sequences
        mapping: {old_id: new_id}
    
    Returns:
        (renamed_sequences, id_map: {new_id: old_id}, unmapped_ids)
    """
    renamed: Dict[str, str] = OrderedDict()
    id_map: Dict[str, str] = OrderedDict()
    unmapped: List[str] = []
    
    for old_id, seq in sequences.items():
        if old_id in mapping:
            new_id = mapping[old_id]
            renamed[new_id] = seq
            id_map[new_id] = old_id
        else:
            renamed[old_id] = seq
            unmapped.append(old_id)
    
    return renamed, id_map, unmapped


def run_rename(
    input_path: Path,
    output_prefix: Path,
    prefix: Optional[str] = None,
    id_map_path: Optional[Path] = None,
    start_index: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Batch rename probe IDs.
    
    Mode 1 (prefix): Rename all to {prefix}_0001, _0002, ...
    Mode 2 (id_map): Rename according to mapping file (TSV: old_id<tab>new_id)
    
    Args:
        input_path: Input FASTA
        output_prefix: Output prefix
        prefix: Prefix for mode 1
        id_map_path: Mapping file for mode 2
        start_index: Starting index for prefix mode
        verbose: Verbose logging
        
    Returns:
        Result with rename statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    if prefix is None and id_map_path is None:
        return Err("Either --prefix or --id_map must be provided")
    
    fasta_result = read_fasta(input_path)
    if fasta_result.is_err():
        return Err(f"Failed to read input: {fasta_result.unwrap_err()}")
    
    sequences = fasta_result.unwrap()
    logger.info(f"Loaded {len(sequences)} sequences")
    
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    unmapped_count = 0
    
    if prefix is not None:
        renamed, id_map = rename_by_prefix(sequences, prefix, start_index)
        logger.info(f"Renamed {len(renamed)} probes with prefix '{prefix}'")
    else:
        import pandas as pd
        try:
            map_df = pd.read_csv(id_map_path, sep='\t', header=None)
            mapping = dict(zip(
                map_df.iloc[:, 0].astype(str),
                map_df.iloc[:, 1].astype(str),
            ))
        except Exception as e:
            return Err(f"Failed to read id_map: {e}")
        
        renamed, id_map, unmapped = rename_by_map(sequences, mapping)
        unmapped_count = len(unmapped)
        
        if unmapped:
            logger.warning(f"{unmapped_count} IDs not found in map, kept original names")
    
    # Write renamed FASTA
    fasta_path = Path(str(output_prefix) + ".renamed.fa")
    write_result = write_fasta(renamed, fasta_path)
    if write_result.is_err():
        return Err(f"Failed to write output: {write_result.unwrap_err()}")
    
    # Write ID map (new_id -> old_id)
    map_out_path = Path(str(output_prefix) + ".id_map.tsv")
    with open(map_out_path, 'w') as f:
        f.write("new_id\told_id\n")
        for new_id, old_id in id_map.items():
            f.write(f"{new_id}\t{old_id}\n")
    
    stats = {
        "probe_count": len(renamed),
        "renamed_count": len(id_map),
        "unmapped_count": unmapped_count,
        "mode": "prefix" if prefix else "id_map",
        "fasta_file": str(fasta_path),
        "map_file": str(map_out_path),
    }
    
    return Ok(stats)
