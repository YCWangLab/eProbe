"""
Format conversion utilities.

Convert between different file formats used in probe design.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any

import pandas as pd

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta
from eprobe.biophysics import calculate_gc, calculate_tm, calculate_complexity

logger = logging.getLogger(__name__)


def vcf_to_snp_tsv(
    input_path: Path,
    output_path: Path,
) -> Result[int, str]:
    """
    Convert VCF to SNP TSV format.
    
    Note: This is a simplified conversion. For full VCF processing
    with flanking sequences, use 'eprobe popgen extract'.
    """
    from eprobe.core.vcf import extract_snps_from_vcf
    from eprobe.core.models import SNPDataFrame
    
    result = extract_snps_from_vcf(input_path)
    if result.is_err():
        return Err(result.unwrap_err())
    
    snps = result.unwrap()
    snp_df = SNPDataFrame.from_snps(snps)
    save_result = snp_df.to_tsv(output_path)
    
    if save_result.is_err():
        return Err(save_result.unwrap_err())
    
    return Ok(len(snps))


def snp_tsv_to_bed(
    input_path: Path,
    output_path: Path,
) -> Result[int, str]:
    """
    Convert SNP TSV to BED format.
    """
    try:
        df = pd.read_csv(input_path, sep='\t')
        
        if 'chrom' not in df.columns or 'pos' not in df.columns:
            return Err("TSV must have 'chrom' and 'pos' columns")
        
        bed_df = pd.DataFrame({
            'chrom': df['chrom'],
            'start': df['pos'] - 1,  # 0-based
            'end': df['pos'],        # 1-based exclusive
            'name': df.get('snp_id', df.index.astype(str)),
        })
        
        bed_df.to_csv(output_path, sep='\t', index=False, header=False)
        return Ok(len(bed_df))
        
    except Exception as e:
        return Err(f"Conversion failed: {e}")


def fasta_to_probe_tsv(
    input_path: Path,
    output_path: Path,
    include_stats: bool = False,
) -> Result[int, str]:
    """
    Convert probe FASTA to TSV with metadata.
    """
    fasta_result = read_fasta(input_path)
    if fasta_result.is_err():
        return Err(fasta_result.unwrap_err())
    
    sequences = fasta_result.unwrap()
    
    data = []
    for seq_id, seq in sequences.items():
        row = {
            'probe_id': seq_id,
            'sequence': seq,
            'length': len(seq),
        }
        
        if include_stats:
            row['gc'] = round(calculate_gc(seq), 2)
            row['tm'] = round(calculate_tm(seq), 2)
            row['complexity'] = round(calculate_complexity(seq), 3)
        
        data.append(row)
    
    df = pd.DataFrame(data)
    df.to_csv(output_path, sep='\t', index=False)
    
    return Ok(len(data))


def run_convert(
    input_path: Path,
    output_path: Path,
    from_format: str,
    to_format: str,
    include_stats: bool = False,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Convert between file formats.
    
    Args:
        input_path: Input file
        output_path: Output file
        from_format: Input format
        to_format: Output format
        include_stats: Include sequence statistics
        verbose: Enable verbose logging
        
    Returns:
        Result containing conversion statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Converting {from_format} → {to_format}")
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    conversion_key = f"{from_format}_{to_format}"
    
    if conversion_key == "vcf_snp_tsv":
        result = vcf_to_snp_tsv(input_path, output_path)
    elif conversion_key == "snp_tsv_bed":
        result = snp_tsv_to_bed(input_path, output_path)
    elif conversion_key == "fasta_probe_tsv":
        result = fasta_to_probe_tsv(input_path, output_path, include_stats)
    else:
        return Err(f"Unsupported conversion: {from_format} → {to_format}")
    
    if result.is_err():
        return Err(result.unwrap_err())
    
    record_count = result.unwrap()
    
    stats = {
        "record_count": record_count,
        "from_format": from_format,
        "to_format": to_format,
        "output_file": str(output_path),
    }
    
    return Ok(stats)
