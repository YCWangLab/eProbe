"""
Melting temperature (Tm) calculator.

Tm indicates the temperature at which 50% of probe-target duplexes dissociate.
Important for designing probes that hybridize efficiently at capture temperature.

Supports multiple calculation methods:
- Nearest-neighbor (NN): Most accurate for oligonucleotides
- Wallace rule: Quick estimation (4°C per GC, 2°C per AT)
- GC-based: Simple percentage-based calculation

For capture probes, optimal Tm is typically 65-85°C depending on the
hybridization conditions.
"""

from __future__ import annotations
from enum import Enum
from typing import Sequence

from Bio.SeqUtils import MeltingTemp as Tm
from Bio.Seq import Seq

from eprobe.core.result import Result, Ok, Err


class TmMethod(str, Enum):
    """Available Tm calculation methods."""
    
    NEAREST_NEIGHBOR = "tm_NN"  # Nearest-neighbor thermodynamic model
    GC_BASED = "tm_GC"  # Simple GC-based formula
    WALLACE = "wallace"  # Wallace rule: 4*GC + 2*AT


class TmTable(str, Enum):
    """
    Thermodynamic parameter tables for nearest-neighbor method.
    
    Different tables are optimized for different experimental conditions:
    - DNA_NN1-4: Various DNA/DNA hybridization parameters
    - R_DNA_NN1: RNA probe / DNA target (recommended for capture)
    """
    
    # DNA/DNA tables
    DNA_NN1 = "DNA_NN1"  # Breslauer et al. (1986)
    DNA_NN2 = "DNA_NN2"  # Sugimoto et al. (1996)
    DNA_NN3 = "DNA_NN3"  # Allawi & SantaLucia (1997)
    DNA_NN4 = "DNA_NN4"  # SantaLucia & Hicks (2004)
    
    # RNA/DNA table (recommended for RNA baits capturing DNA)
    R_DNA_NN1 = "R_DNA_NN1"  # Sugimoto et al. (1995)


# Mapping from TmTable enum to Biopython table objects
_TM_TABLE_MAP = {
    TmTable.DNA_NN1: Tm.DNA_NN1,
    TmTable.DNA_NN2: Tm.DNA_NN2,
    TmTable.DNA_NN3: Tm.DNA_NN3,
    TmTable.DNA_NN4: Tm.DNA_NN4,
    TmTable.R_DNA_NN1: Tm.R_DNA_NN1,
}


def calculate_tm(
    sequence: str,
    method: TmMethod = TmMethod.NEAREST_NEIGHBOR,
    table: TmTable = TmTable.R_DNA_NN1,
) -> Result[float, str]:
    """
    Calculate melting temperature of DNA sequence.
    
    Args:
        sequence: DNA sequence (A/T/C/G)
        method: Calculation method (default: nearest-neighbor)
        table: NN parameter table (only used with NEAREST_NEIGHBOR method)
        
    Returns:
        Ok(tm_celsius) on success
        Err(message) on failure
        
    Example:
        >>> calculate_tm("GCGCATCGATCGATCG").unwrap()
        52.5  # approximate value
    """
    if not sequence:
        return Err("Empty sequence provided")
    
    sequence = sequence.upper()
    
    # Validate sequence
    valid_bases = set("ATCG")
    invalid = set(sequence) - valid_bases
    if invalid:
        return Err(f"Invalid bases in sequence: {invalid}")
    
    try:
        if method == TmMethod.NEAREST_NEIGHBOR:
            # For RNA/DNA table, need to transcribe sequence
            if table == TmTable.R_DNA_NN1:
                rna_seq = Seq(sequence).transcribe()
                tm = Tm.Tm_NN(rna_seq, nn_table=_TM_TABLE_MAP[table])
            else:
                tm = Tm.Tm_NN(sequence, nn_table=_TM_TABLE_MAP[table])
        
        elif method == TmMethod.GC_BASED:
            tm = Tm.Tm_GC(sequence)
        
        elif method == TmMethod.WALLACE:
            # Wallace rule: 4°C per G/C, 2°C per A/T
            gc_count = sequence.count("G") + sequence.count("C")
            at_count = sequence.count("A") + sequence.count("T")
            tm = 4 * gc_count + 2 * at_count
        
        else:
            return Err(f"Unknown Tm method: {method}")
        
        return Ok(round(float(tm), 4))
    
    except Exception as e:
        return Err(f"Tm calculation failed: {e}")


def calculate_tm_batch(
    sequences: Sequence[str],
    method: TmMethod = TmMethod.NEAREST_NEIGHBOR,
    table: TmTable = TmTable.R_DNA_NN1,
) -> list[Result[float, str]]:
    """
    Calculate Tm for multiple sequences.
    
    Args:
        sequences: List of DNA sequences
        method: Calculation method
        table: NN parameter table
        
    Returns:
        List of Result objects (one per sequence)
    """
    return [calculate_tm(seq, method, table) for seq in sequences]


def tm_in_range(
    sequence: str,
    min_tm: float,
    max_tm: float,
    method: TmMethod = TmMethod.NEAREST_NEIGHBOR,
    table: TmTable = TmTable.R_DNA_NN1,
) -> Result[bool, str]:
    """
    Check if Tm falls within specified range.
    
    Args:
        sequence: DNA sequence
        min_tm: Minimum Tm (inclusive)
        max_tm: Maximum Tm (inclusive)
        method: Calculation method
        table: NN parameter table
        
    Returns:
        Ok(True) if in range, Ok(False) if not, Err on failure
    """
    result = calculate_tm(sequence, method, table)
    if result.is_err():
        return result  # type: ignore
    
    tm = result.unwrap()
    return Ok(min_tm <= tm <= max_tm)
