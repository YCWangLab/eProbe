"""
Probe subsetting and filtering.

Filter or sample probes based on various criteria.
"""

import logging
import random
from pathlib import Path
from typing import Optional, Dict, Any, List
from collections import defaultdict

import numpy as np
import pandas as pd

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def filter_by_expression(
    df: pd.DataFrame,
    expression: str,
) -> pd.DataFrame:
    """
    Filter DataFrame using expression string.
    
    Args:
        df: Input DataFrame
        expression: Filter expression
        
    Returns:
        Filtered DataFrame
    """
    try:
        return df.query(expression)
    except Exception as e:
        logger.error(f"Filter expression error: {e}")
        return df


def sample_random(
    sequences: Dict[str, str],
    n: int,
    seed: int = 42,
) -> Dict[str, str]:
    """
    Random sample sequences.
    
    Args:
        sequences: Input sequences
        n: Number to sample
        seed: Random seed
        
    Returns:
        Sampled sequences
    """
    random.seed(seed)
    
    if len(sequences) <= n:
        return sequences
    
    selected_ids = random.sample(list(sequences.keys()), n)
    return {k: sequences[k] for k in selected_ids}


def sample_uniform(
    sequences: Dict[str, str],
    metadata: pd.DataFrame,
    n: int,
    window_size: int = 100000,
    seed: int = 42,
) -> Dict[str, str]:
    """
    Sample sequences with uniform genomic distribution.
    
    Requires metadata with 'chrom' and 'start' columns.
    
    Args:
        sequences: Input sequences
        metadata: DataFrame with genomic coordinates
        n: Target number to sample
        window_size: Window size for uniform distribution
        seed: Random seed
        
    Returns:
        Sampled sequences
    """
    random.seed(seed)
    np.random.seed(seed)
    
    if len(sequences) <= n:
        return sequences
    
    # Group by window
    if 'chrom' not in metadata.columns or 'start' not in metadata.columns:
        logger.warning("Missing chrom/start columns, falling back to random")
        return sample_random(sequences, n, seed)
    
    # Add window column
    metadata = metadata.copy()
    metadata['window'] = metadata['start'] // window_size
    
    # Group by chrom and window
    window_groups: Dict[tuple, List[str]] = defaultdict(list)
    
    for _, row in metadata.iterrows():
        seq_id = row.get('probe_id') or row.get('snp_id') or str(row.name)
        if seq_id in sequences:
            key = (row['chrom'], row['window'])
            window_groups[key].append(seq_id)
    
    # Distribute samples across windows
    n_windows = len(window_groups)
    base_per_window = n // n_windows
    remainder = n % n_windows
    
    selected_ids = []
    sorted_windows = sorted(window_groups.keys())
    
    for i, window_key in enumerate(sorted_windows):
        candidates = window_groups[window_key]
        n_select = base_per_window + (1 if i < remainder else 0)
        
        if len(candidates) <= n_select:
            selected_ids.extend(candidates)
        else:
            selected_ids.extend(random.sample(candidates, n_select))
    
    return {k: sequences[k] for k in selected_ids}


def run_subset(
    input_path: Path,
    output_path: Path,
    filter_expr: Optional[str] = None,
    sample_count: Optional[int] = None,
    chromosomes: Optional[List[str]] = None,
    bed_path: Optional[Path] = None,
    reference_path: Optional[Path] = None,
    strategy: str = "random",
    window_size: int = 100000,
    seed: int = 42,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Subset probes by filtering or sampling.
    
    Args:
        input_path: Input file (FASTA or TSV)
        output_path: Output file
        filter_expr: Filter expression for TSV
        sample_count: Target sample size
        chromosomes: Chromosomes to keep
        bed_path: BED file for region filtering
        reference_path: Reference genome (for uniform sampling)
        strategy: Sampling strategy (random or uniform)
        window_size: Window size for uniform sampling
        seed: Random seed
        verbose: Enable verbose logging
        
    Returns:
        Result containing subset statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Subsetting {input_path}")
    
    suffix = input_path.suffix.lower()
    is_fasta = suffix in ['.fa', '.fasta', '.fna']
    
    if is_fasta:
        # FASTA input
        fasta_result = read_fasta(input_path)
        if fasta_result.is_err():
            return Err(f"Failed to read input: {fasta_result.unwrap_err()}")
        
        sequences = fasta_result.unwrap()
        input_count = len(sequences)
        metadata = None
        
    else:
        # TSV input
        try:
            df = pd.read_csv(input_path, sep='\t')
            input_count = len(df)
            
            # Apply filter expression
            if filter_expr:
                df = filter_by_expression(df, filter_expr)
                logger.info(f"After filter: {len(df)} rows")
            
            # Filter by chromosomes
            if chromosomes and 'chrom' in df.columns:
                df = df[df['chrom'].isin(chromosomes)]
                logger.info(f"After chromosome filter: {len(df)} rows")
            
            # Build sequence dict from TSV
            if 'sequence' in df.columns:
                id_col = 'probe_id' if 'probe_id' in df.columns else df.columns[0]
                sequences = dict(zip(df[id_col], df['sequence']))
            else:
                # No sequences - just filter TSV
                output_path.parent.mkdir(parents=True, exist_ok=True)
                
                if sample_count and len(df) > sample_count:
                    df = df.sample(n=sample_count, random_state=seed)
                
                df.to_csv(output_path, sep='\t', index=False)
                
                return Ok({
                    "input_count": input_count,
                    "output_count": len(df),
                    "output_file": str(output_path),
                })
            
            metadata = df
            
        except Exception as e:
            return Err(f"Failed to read input: {e}")
    
    logger.info(f"Loaded {input_count} sequences")
    
    # Apply sampling
    if sample_count and len(sequences) > sample_count:
        if strategy == "uniform" and metadata is not None:
            sequences = sample_uniform(sequences, metadata, sample_count, window_size, seed)
        else:
            sequences = sample_random(sequences, sample_count, seed)
        
        logger.info(f"Sampled to {len(sequences)} sequences")
    
    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if is_fasta or output_path.suffix.lower() in ['.fa', '.fasta', '.fna']:
        write_result = write_fasta(sequences, output_path)
        if write_result.is_err():
            return Err(f"Failed to write output: {write_result.unwrap_err()}")
    else:
        # Write as TSV
        data = [{"id": k, "sequence": v} for k, v in sequences.items()]
        pd.DataFrame(data).to_csv(output_path, sep='\t', index=False)
    
    stats = {
        "input_count": input_count,
        "output_count": len(sequences),
        "strategy": strategy,
        "output_file": str(output_path),
    }
    
    return Ok(stats)
