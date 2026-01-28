"""
FASTA file operations.

Provides pure functions for reading and writing FASTA files with
proper error handling using Result types.
"""

from __future__ import annotations
from collections import OrderedDict
from pathlib import Path
from typing import Iterator

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from eprobe.core.result import Result, Ok, Err


# Type alias for FASTA dictionary
FastaDict = OrderedDict[str, str]


def read_fasta(path: Path | str) -> Result[FastaDict, str]:
    """
    Read FASTA file into ordered dictionary.
    
    Args:
        path: Path to FASTA file
        
    Returns:
        Ok(FastaDict) mapping sequence IDs to sequences (uppercase)
        Err(message) on failure
        
    Example:
        >>> result = read_fasta("reference.fa")
        >>> if result.is_ok():
        ...     genome = result.unwrap()
        ...     print(genome["chr1"][:100])
    """
    try:
        path = Path(path)
        if not path.exists():
            return Err(f"FASTA file not found: {path}")
        
        fasta_dict: FastaDict = OrderedDict()
        
        for record in SeqIO.parse(path, "fasta"):
            # Store uppercase sequence, strip any whitespace
            fasta_dict[record.id] = str(record.seq).upper().replace("\n", "")
        
        if not fasta_dict:
            return Err(f"No sequences found in FASTA file: {path}")
        
        return Ok(fasta_dict)
    
    except Exception as e:
        return Err(f"Failed to read FASTA file: {e}")


def write_fasta(
    sequences: FastaDict | dict[str, str],
    path: Path | str,
    line_width: int = 80,
) -> Result[Path, str]:
    """
    Write sequences to FASTA file.
    
    Args:
        sequences: Dictionary mapping IDs to sequences
        path: Output file path
        line_width: Number of characters per line (0 for no wrapping)
        
    Returns:
        Ok(path) on success, Err(message) on failure
    """
    try:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(path, "w") as f:
            for seq_id, sequence in sequences.items():
                f.write(f">{seq_id}\n")
                
                if line_width > 0:
                    # Wrap sequence at specified width
                    for i in range(0, len(sequence), line_width):
                        f.write(f"{sequence[i:i + line_width]}\n")
                else:
                    f.write(f"{sequence}\n")
        
        return Ok(path)
    
    except Exception as e:
        return Err(f"Failed to write FASTA file: {e}")


def get_chromosome_sizes(fasta_dict: FastaDict) -> dict[str, int]:
    """
    Get chromosome/contig sizes from FASTA dictionary.
    
    Args:
        fasta_dict: Dictionary mapping IDs to sequences
        
    Returns:
        Dictionary mapping chromosome IDs to lengths
    """
    return {chrom: len(seq) for chrom, seq in fasta_dict.items()}


def write_chromosome_sizes(
    fasta_dict: FastaDict,
    path: Path | str,
) -> Result[Path, str]:
    """
    Write chromosome sizes to TSV file.
    
    Args:
        fasta_dict: FASTA dictionary
        path: Output path
        
    Returns:
        Ok(path) on success, Err(message) on failure
    """
    try:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(path, "w") as f:
            for chrom, seq in fasta_dict.items():
                f.write(f"{chrom}\t{len(seq)}\n")
        
        return Ok(path)
    
    except Exception as e:
        return Err(f"Failed to write chromosome sizes: {e}")


def read_chromosome_sizes(path: Path | str) -> Result[dict[str, int], str]:
    """
    Read chromosome sizes from TSV file.
    
    Args:
        path: Path to TSV file (chrom<tab>size format)
        
    Returns:
        Ok(dict) mapping chromosome to size, Err(message) on failure
    """
    try:
        path = Path(path)
        if not path.exists():
            return Err(f"Chromosome size file not found: {path}")
        
        sizes: dict[str, int] = {}
        
        with open(path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 2:
                    return Err(f"Invalid format at line {line_num}: expected 'chrom<tab>size'")
                
                chrom = parts[0]
                try:
                    size = int(parts[1])
                except ValueError:
                    return Err(f"Invalid size at line {line_num}: '{parts[1]}'")
                
                sizes[chrom] = size
        
        return Ok(sizes)
    
    except Exception as e:
        return Err(f"Failed to read chromosome sizes: {e}")


def extract_sequence(
    fasta_dict: FastaDict,
    chrom: str,
    start: int,
    end: int,
) -> Result[str, str]:
    """
    Extract sequence from genomic region.
    
    Args:
        fasta_dict: FASTA dictionary
        chrom: Chromosome name
        start: Start position (1-based, inclusive)
        end: End position (1-based, inclusive)
        
    Returns:
        Ok(sequence) on success, Err(message) on failure
        
    Note:
        Coordinates are 1-based inclusive, following genomic convention.
        Internally converts to 0-based for Python slicing.
    """
    if chrom not in fasta_dict:
        return Err(f"Chromosome not found: {chrom}")
    
    chrom_seq = fasta_dict[chrom]
    chrom_len = len(chrom_seq)
    
    if start < 1:
        return Err(f"Start position must be >= 1, got {start}")
    
    if end > chrom_len:
        return Err(f"End position {end} exceeds chromosome length {chrom_len}")
    
    if start > end:
        return Err(f"Start ({start}) must be <= end ({end})")
    
    # Convert to 0-based indexing for Python slice
    sequence = chrom_seq[start - 1:end]
    
    return Ok(sequence)


def split_fasta_dict(
    fasta_dict: FastaDict,
    num_parts: int,
) -> list[FastaDict]:
    """
    Split FASTA dictionary into multiple parts for parallel processing.
    
    Distributes sequences round-robin across parts to balance workload.
    
    Args:
        fasta_dict: Input FASTA dictionary
        num_parts: Number of parts to split into
        
    Returns:
        List of FastaDict objects
    """
    if num_parts <= 0:
        raise ValueError("num_parts must be positive")
    
    parts: list[FastaDict] = [OrderedDict() for _ in range(num_parts)]
    
    for i, (seq_id, sequence) in enumerate(fasta_dict.items()):
        part_idx = i % num_parts
        parts[part_idx][seq_id] = sequence
    
    return parts


def reverse_complement(sequence: str) -> str:
    """
    Get reverse complement of DNA sequence.
    
    Args:
        sequence: DNA sequence (A/T/C/G/N)
        
    Returns:
        Reverse complement sequence
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(base, "N") for base in reversed(sequence.upper()))
