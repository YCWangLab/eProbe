"""
Merge multiple probe sets.

Combines FASTA or TSV files from different sources into a single file.
Optionally tags sequences with their source.
"""

import logging
from pathlib import Path
from typing import Dict, Any, List, Tuple

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def run_merge(
    input_files: List[Tuple[str, Path]],
    output_path: Path,
    tag_source: bool = False,
    file_format: str = "fasta",
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Merge multiple probe files.
    
    Args:
        input_files: List of (tag, path) tuples
        output_path: Output file path
        tag_source: Add source prefix to sequence IDs
        file_format: Input/output format (fasta or tsv)
        verbose: Enable verbose logging
        
    Returns:
        Result containing merge statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Merging {len(input_files)} files")
    
    merged_sequences: Dict[str, str] = {}
    source_counts: Dict[str, int] = {}
    duplicates = 0
    
    if file_format == "fasta":
        for tag, path in input_files:
            logger.info(f"Reading {tag}: {path}")
            
            fasta_result = read_fasta(path)
            if fasta_result.is_err():
                return Err(f"Failed to read {path}: {fasta_result.unwrap_err()}")
            
            sequences = fasta_result.unwrap()
            source_counts[tag] = len(sequences)
            
            for seq_id, seq in sequences.items():
                # Apply source tag if requested
                if tag_source:
                    new_id = f"{tag}_{seq_id}"
                else:
                    new_id = seq_id
                
                # Handle duplicates
                if new_id in merged_sequences:
                    duplicates += 1
                    new_id = f"{new_id}_{duplicates}"
                
                merged_sequences[new_id] = seq
        
        # Write merged output
        output_path.parent.mkdir(parents=True, exist_ok=True)
        write_result = write_fasta(merged_sequences, output_path)
        if write_result.is_err():
            return Err(f"Failed to write output: {write_result.unwrap_err()}")
    
    elif file_format == "tsv":
        # TSV merge
        import pandas as pd
        
        dfs = []
        for tag, path in input_files:
            logger.info(f"Reading {tag}: {path}")
            
            try:
                df = pd.read_csv(path, sep='\t')
                if tag_source:
                    df['source'] = tag
                dfs.append(df)
                source_counts[tag] = len(df)
            except Exception as e:
                return Err(f"Failed to read {path}: {e}")
        
        merged_df = pd.concat(dfs, ignore_index=True)
        merged_df.to_csv(output_path, sep='\t', index=False)
        
        merged_sequences = {str(i): "" for i in range(len(merged_df))}
    
    else:
        return Err(f"Unsupported format: {file_format}")
    
    logger.info(f"Merged {len(merged_sequences)} sequences to {output_path}")
    
    stats = {
        "total_probes": len(merged_sequences),
        "sources": source_counts,
        "duplicates_renamed": duplicates,
        "output_file": str(output_path),
    }
    
    return Ok(stats)
