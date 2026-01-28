"""
Thermodynamic entropy (ΔS°) calculator for DNA hybridization.

This module calculates the standard entropy change (ΔS°) for probe-target
hybridization using the nearest-neighbor (NN) thermodynamic model.

The entropy values are from Sugimoto et al. (1996) for DNA/DNA duplexes,
consistent with the parameters used in the eProbe manuscript.

Reference:
    Sugimoto N, Nakano S, Yoneyama M, Honda K. (1996)
    Improved thermodynamic parameters and helix initiation factor to predict
    stability of DNA duplexes.
    Nucleic Acids Res. 24(22):4501-4505.
    DOI: 10.1093/nar/24.22.4501

Units:
    ΔS° is reported in cal·mol⁻¹·K⁻¹ (entropy units, e.u.)
    
Note:
    More negative ΔS° indicates greater loss of conformational freedom upon
    hybridization, typically correlating with stronger base stacking (GC-rich).
"""

from __future__ import annotations
from enum import Enum
from typing import Sequence

from eprobe.core.result import Result, Ok, Err


class NNTable(str, Enum):
    """
    Available nearest-neighbor parameter tables.
    
    SUGIMOTO_1996 is the default, matching the eProbe manuscript.
    """
    
    SUGIMOTO_1996 = "sugimoto_1996"  # DNA/DNA, Sugimoto et al. (1996)
    SANTALUCIA_1998 = "santalucia_1998"  # DNA/DNA, unified parameters


# Sugimoto et al. (1996) - DNA/DNA thermodynamic parameters
# Format: 'XY/X'Y'': (ΔH° in kcal/mol, ΔS° in cal/(mol·K))
# Where XY is 5'->3' on one strand, X'Y' is 3'->5' on complementary strand
SUGIMOTO_1996_PARAMS = {
    # Nearest-neighbor dinucleotides (ΔH°, ΔS°)
    "AA/TT": (-8.0, -21.9),
    "AT/TA": (-5.6, -15.2),
    "TA/AT": (-6.6, -18.4),
    "CA/GT": (-8.2, -21.0),
    "GT/CA": (-9.4, -25.5),
    "CT/GA": (-6.6, -16.4),
    "GA/CT": (-8.8, -23.5),
    "CG/GC": (-11.8, -29.0),
    "GC/CG": (-10.5, -26.4),
    "GG/CC": (-10.9, -28.4),
    # Initiation parameter
    "init": (0.6, -9.0),
    # Symmetry correction (for self-complementary sequences)
    "sym": (0.0, -1.4),
}

# SantaLucia (1998) unified parameters (alternative)
# Commonly used, slightly different from Sugimoto
SANTALUCIA_1998_PARAMS = {
    "AA/TT": (-7.9, -22.2),
    "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
    # Initiation with terminal G/C
    "init_G/C": (0.1, -2.8),
    # Initiation with terminal A/T
    "init_A/T": (2.3, 4.1),
    # Symmetry correction
    "sym": (0.0, -1.4),
}

# Mapping from base to complement
_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def _get_nn_params(base1: str, base2: str, params: dict) -> tuple[str, tuple[float, float]] | tuple[None, None]:
    """
    Get the nearest-neighbor parameters for a dinucleotide step.
    
    NN parameter tables only store 10 unique dinucleotide pairs (considering
    Watson-Crick symmetry). This function looks up the correct key by trying
    both the forward and reverse orientations.
    
    For sequence 5'-XY-3':
    - Forward key: XY/X'Y' (X' is complement of X)
    - Reverse key: Y'X'/YX (same pair from complement strand perspective)
    
    Args:
        base1: First base (5' position)
        base2: Second base (3' position)
        params: Parameter dictionary
        
    Returns:
        Tuple of (key, (ΔH°, ΔS°)) if found, (None, None) otherwise
    """
    comp1 = _COMPLEMENT[base1]
    comp2 = _COMPLEMENT[base2]
    
    # Forward key: XY/X'Y'
    forward_key = f"{base1}{base2}/{comp1}{comp2}"
    if forward_key in params:
        return forward_key, params[forward_key]
    
    # Reverse key: Y'X'/YX (same base-pair step from other strand)
    reverse_key = f"{comp2}{comp1}/{base2}{base1}"
    if reverse_key in params:
        return reverse_key, params[reverse_key]
    
    return None, None


def _is_self_complementary(sequence: str) -> bool:
    """Check if sequence is self-complementary (palindromic)."""
    complement = "".join(_COMPLEMENT[b] for b in sequence)
    return sequence == complement[::-1]


def calculate_thermo_entropy(
    sequence: str,
    table: NNTable = NNTable.SUGIMOTO_1996,
) -> Result[float, str]:
    """
    Calculate thermodynamic entropy (ΔS°) for DNA hybridization.
    
    Uses nearest-neighbor model to sum entropy contributions from each
    dinucleotide step, plus initiation and symmetry corrections.
    
    Args:
        sequence: DNA sequence (5'->3', probe sequence)
        table: Thermodynamic parameter table to use
        
    Returns:
        Ok(entropy) in cal·mol⁻¹·K⁻¹ on success
        Err(message) on failure
        
    Example:
        >>> calculate_thermo_entropy("GCGATCGATCGC").unwrap()
        -295.3  # approximate value
        
    Note:
        More negative values indicate:
        - Greater entropy loss upon hybridization
        - Stronger base stacking (GC-rich sequences)
        - Higher thermal stability
    """
    if not sequence:
        return Err("Empty sequence provided")
    
    sequence = sequence.upper()
    
    # Validate sequence
    valid_bases = set("ATCG")
    invalid = set(sequence) - valid_bases
    if invalid:
        return Err(f"Invalid bases in sequence: {invalid}")
    
    if len(sequence) < 2:
        return Err("Sequence must be at least 2 bases")
    
    # Select parameter table
    if table == NNTable.SUGIMOTO_1996:
        params = SUGIMOTO_1996_PARAMS
    elif table == NNTable.SANTALUCIA_1998:
        params = SANTALUCIA_1998_PARAMS
    else:
        return Err(f"Unknown parameter table: {table}")
    
    try:
        # Sum nearest-neighbor entropy contributions
        total_entropy = 0.0
        
        for i in range(len(sequence) - 1):
            base1 = sequence[i]
            base2 = sequence[i + 1]
            nn_key, nn_params = _get_nn_params(base1, base2, params)
            
            if nn_key is None:
                return Err(f"Missing NN parameter for {base1}{base2}")
            
            _, delta_s = nn_params  # (ΔH°, ΔS°)
            total_entropy += delta_s
        
        # Add initiation parameter
        if table == NNTable.SUGIMOTO_1996:
            _, init_s = params["init"]
            total_entropy += init_s
        elif table == NNTable.SANTALUCIA_1998:
            # SantaLucia uses different init based on terminal bases
            first_base = sequence[0]
            last_base = sequence[-1]
            
            if first_base in "GC":
                _, init_s = params["init_G/C"]
                total_entropy += init_s
            else:
                _, init_s = params["init_A/T"]
                total_entropy += init_s
            
            if last_base in "GC":
                _, init_s = params["init_G/C"]
                total_entropy += init_s
            else:
                _, init_s = params["init_A/T"]
                total_entropy += init_s
        
        # Symmetry correction for self-complementary sequences
        if _is_self_complementary(sequence):
            _, sym_s = params["sym"]
            total_entropy += sym_s
        
        return Ok(round(total_entropy, 2))
    
    except Exception as e:
        return Err(f"Entropy calculation failed: {e}")


def calculate_thermo_entropy_batch(
    sequences: Sequence[str],
    table: NNTable = NNTable.SUGIMOTO_1996,
) -> list[Result[float, str]]:
    """
    Calculate thermodynamic entropy for multiple sequences.
    
    Args:
        sequences: List of DNA sequences
        table: Thermodynamic parameter table
        
    Returns:
        List of Result objects (one per sequence)
    """
    return [calculate_thermo_entropy(seq, table) for seq in sequences]


def calculate_thermo_enthalpy(
    sequence: str,
    table: NNTable = NNTable.SUGIMOTO_1996,
) -> Result[float, str]:
    """
    Calculate thermodynamic enthalpy (ΔH°) for DNA hybridization.
    
    Complementary function to entropy calculation, using the same NN model.
    
    Args:
        sequence: DNA sequence (5'->3', probe sequence)
        table: Thermodynamic parameter table to use
        
    Returns:
        Ok(enthalpy) in kcal·mol⁻¹ on success
        Err(message) on failure
        
    Note:
        More negative ΔH° indicates stronger hydrogen bonding (more stable duplex).
    """
    if not sequence:
        return Err("Empty sequence provided")
    
    sequence = sequence.upper()
    
    valid_bases = set("ATCG")
    invalid = set(sequence) - valid_bases
    if invalid:
        return Err(f"Invalid bases in sequence: {invalid}")
    
    if len(sequence) < 2:
        return Err("Sequence must be at least 2 bases")
    
    if table == NNTable.SUGIMOTO_1996:
        params = SUGIMOTO_1996_PARAMS
    elif table == NNTable.SANTALUCIA_1998:
        params = SANTALUCIA_1998_PARAMS
    else:
        return Err(f"Unknown parameter table: {table}")
    
    try:
        total_enthalpy = 0.0
        
        for i in range(len(sequence) - 1):
            base1 = sequence[i]
            base2 = sequence[i + 1]
            nn_key, nn_params = _get_nn_params(base1, base2, params)
            
            if nn_key is None:
                return Err(f"Missing NN parameter for {base1}{base2}")
            
            delta_h, _ = nn_params
            total_enthalpy += delta_h
        
        # Initiation
        if table == NNTable.SUGIMOTO_1996:
            init_h, _ = params["init"]
            total_enthalpy += init_h
        elif table == NNTable.SANTALUCIA_1998:
            first_base = sequence[0]
            last_base = sequence[-1]
            
            if first_base in "GC":
                init_h, _ = params["init_G/C"]
                total_enthalpy += init_h
            else:
                init_h, _ = params["init_A/T"]
                total_enthalpy += init_h
            
            if last_base in "GC":
                init_h, _ = params["init_G/C"]
                total_enthalpy += init_h
            else:
                init_h, _ = params["init_A/T"]
                total_enthalpy += init_h
        
        # Symmetry correction
        if _is_self_complementary(sequence):
            sym_h, _ = params["sym"]
            total_enthalpy += sym_h
        
        return Ok(round(total_enthalpy, 2))
    
    except Exception as e:
        return Err(f"Enthalpy calculation failed: {e}")


def calculate_gibbs_energy(
    sequence: str,
    temperature: float = 310.15,  # 37°C in Kelvin
    table: NNTable = NNTable.SUGIMOTO_1996,
) -> Result[float, str]:
    """
    Calculate Gibbs free energy (ΔG°) for DNA hybridization.
    
    ΔG° = ΔH° - T·ΔS°
    
    Args:
        sequence: DNA sequence
        temperature: Temperature in Kelvin (default 37°C = 310.15K)
        table: Thermodynamic parameter table
        
    Returns:
        Ok(gibbs_energy) in kcal·mol⁻¹ on success
        Err(message) on failure
        
    Note:
        More negative ΔG° indicates thermodynamically favorable hybridization.
    """
    enthalpy_result = calculate_thermo_enthalpy(sequence, table)
    if enthalpy_result.is_err():
        return enthalpy_result
    
    entropy_result = calculate_thermo_entropy(sequence, table)
    if entropy_result.is_err():
        return Err(entropy_result.unwrap_err())
    
    delta_h = enthalpy_result.unwrap()  # kcal/mol
    delta_s = entropy_result.unwrap()  # cal/(mol·K)
    
    # Convert ΔS to kcal/(mol·K) for calculation
    delta_s_kcal = delta_s / 1000.0
    
    # ΔG = ΔH - TΔS
    delta_g = delta_h - (temperature * delta_s_kcal)
    
    return Ok(round(delta_g, 2))
