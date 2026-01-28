"""
File and database validation utilities.

Validate file formats and external database availability.
"""

import logging
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any, List

import pandas as pd

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta

logger = logging.getLogger(__name__)


def validate_snp_dataframe(file_path: Path) -> Result[Dict[str, Any], str]:
    """
    Validate SNP TSV file format.
    
    Required columns: chrom, pos, ref, alt
    Optional columns: snp_id, left_flank, right_flank
    """
    try:
        df = pd.read_csv(file_path, sep='\t', nrows=100)  # Sample first 100 rows
        
        required_cols = {'chrom', 'pos', 'ref', 'alt'}
        missing = required_cols - set(df.columns)
        
        if missing:
            return Ok({
                "valid": False,
                "errors": [f"Missing required columns: {missing}"],
            })
        
        messages = []
        
        # Check data types
        if not pd.api.types.is_integer_dtype(df['pos']):
            messages.append("Warning: 'pos' column should be integer")
        
        # Check ref/alt are single characters (for SNPs)
        if (df['ref'].str.len() > 1).any():
            messages.append("Warning: some ref alleles have length > 1 (indels)")
        
        if (df['alt'].str.len() > 1).any():
            messages.append("Warning: some alt alleles have length > 1 (indels)")
        
        # Check for duplicates
        dup_count = df.duplicated(subset=['chrom', 'pos']).sum()
        if dup_count > 0:
            messages.append(f"Warning: {dup_count} duplicate positions found")
        
        return Ok({
            "valid": True,
            "messages": messages,
            "row_count": len(df),
            "columns": list(df.columns),
        })
        
    except Exception as e:
        return Ok({
            "valid": False,
            "errors": [f"Failed to read file: {e}"],
        })


def validate_probe_fasta(file_path: Path) -> Result[Dict[str, Any], str]:
    """
    Validate probe FASTA file format.
    """
    fasta_result = read_fasta(file_path)
    
    if fasta_result.is_err():
        return Ok({
            "valid": False,
            "errors": [fasta_result.unwrap_err()],
        })
    
    sequences = fasta_result.unwrap()
    
    if not sequences:
        return Ok({
            "valid": False,
            "errors": ["No sequences found in FASTA file"],
        })
    
    messages = []
    
    # Check sequence properties
    lengths = [len(seq) for seq in sequences.values()]
    unique_lengths = set(lengths)
    
    if len(unique_lengths) > 1:
        messages.append(f"Variable probe lengths: {min(lengths)}-{max(lengths)}")
    else:
        messages.append(f"Uniform probe length: {lengths[0]}")
    
    # Check for invalid characters
    valid_bases = set('ACGTN')
    invalid_seqs = []
    
    for seq_id, seq in list(sequences.items())[:100]:  # Sample first 100
        if not set(seq.upper()).issubset(valid_bases):
            invalid_seqs.append(seq_id)
    
    if invalid_seqs:
        messages.append(f"Sequences with invalid bases: {len(invalid_seqs)}")
    
    return Ok({
        "valid": True,
        "messages": messages,
        "sequence_count": len(sequences),
        "length_range": [min(lengths), max(lengths)],
    })


def validate_bed_file(file_path: Path) -> Result[Dict[str, Any], str]:
    """
    Validate BED file format.
    """
    try:
        with open(file_path) as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
        
        if not lines:
            return Ok({
                "valid": False,
                "errors": ["No data lines found in BED file"],
            })
        
        messages = []
        errors = []
        
        for i, line in enumerate(lines[:100]):
            parts = line.split('\t')
            
            if len(parts) < 3:
                errors.append(f"Line {i+1}: insufficient columns")
                continue
            
            try:
                start = int(parts[1])
                end = int(parts[2])
                
                if start < 0:
                    errors.append(f"Line {i+1}: negative start coordinate")
                if end < start:
                    errors.append(f"Line {i+1}: end < start")
            except ValueError:
                errors.append(f"Line {i+1}: invalid coordinates")
        
        if errors:
            return Ok({
                "valid": False,
                "errors": errors[:10],  # Limit error messages
            })
        
        messages.append(f"Total regions: {len(lines)}")
        
        return Ok({
            "valid": True,
            "messages": messages,
            "region_count": len(lines),
        })
        
    except Exception as e:
        return Ok({
            "valid": False,
            "errors": [f"Failed to read file: {e}"],
        })


def validate_kraken2_db(db_path: Path) -> Result[Dict[str, Any], str]:
    """
    Validate Kraken2 database.
    """
    # Check required files
    required_files = ['hash.k2d', 'opts.k2d', 'taxo.k2d']
    missing = []
    
    for fname in required_files:
        if not (db_path / fname).exists():
            missing.append(fname)
    
    if missing:
        return Ok({
            "valid": False,
            "errors": [f"Missing database files: {missing}"],
        })
    
    # Try running kraken2 --version to check installation
    try:
        result = subprocess.run(
            ["kraken2", "--version"],
            capture_output=True,
            text=True,
        )
        version = result.stdout.split('\n')[0] if result.returncode == 0 else "Unknown"
    except FileNotFoundError:
        return Ok({
            "valid": False,
            "errors": ["Kraken2 not found in PATH"],
        })
    
    return Ok({
        "valid": True,
        "messages": [f"Kraken2 version: {version}"],
        "database_path": str(db_path),
    })


def validate_bowtie2_index(index_path: Path) -> Result[Dict[str, Any], str]:
    """
    Validate Bowtie2 index.
    """
    # Check for index files
    index_extensions = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    
    found_files = []
    missing_files = []
    
    for ext in index_extensions:
        full_path = Path(str(index_path) + ext)
        if full_path.exists():
            found_files.append(ext)
        else:
            missing_files.append(ext)
    
    if missing_files:
        return Ok({
            "valid": False,
            "errors": [f"Missing index files: {missing_files}"],
        })
    
    # Check bowtie2 installation
    try:
        result = subprocess.run(
            ["bowtie2", "--version"],
            capture_output=True,
            text=True,
        )
        version = result.stdout.split('\n')[0] if result.returncode == 0 else "Unknown"
    except FileNotFoundError:
        return Ok({
            "valid": False,
            "errors": ["Bowtie2 not found in PATH"],
        })
    
    return Ok({
        "valid": True,
        "messages": [f"Bowtie2 version: {version}"],
        "index_prefix": str(index_path),
    })


def run_validate(
    file_path: Optional[Path] = None,
    schema: Optional[str] = None,
    database_path: Optional[Path] = None,
    database_type: Optional[str] = None,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Validate files or databases.
    
    Args:
        file_path: Path to file to validate
        schema: Schema type for file validation
        database_path: Path to database
        database_type: Database type
        verbose: Enable verbose logging
        
    Returns:
        Result containing validation results
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    if file_path and schema:
        # File validation
        logger.info(f"Validating {file_path} as {schema}")
        
        if schema == "snp_dataframe":
            return validate_snp_dataframe(file_path)
        elif schema == "probe_fasta":
            return validate_probe_fasta(file_path)
        elif schema == "bed":
            return validate_bed_file(file_path)
        else:
            return Err(f"Unknown schema: {schema}")
    
    elif database_path and database_type:
        # Database validation
        logger.info(f"Validating {database_type} database at {database_path}")
        
        if database_type == "kraken2":
            return validate_kraken2_db(database_path)
        elif database_type == "bowtie2":
            return validate_bowtie2_index(database_path)
        else:
            return Err(f"Unknown database type: {database_type}")
    
    else:
        return Err("Invalid validation parameters")
