"""
Biophysical property calculators for probe sequences.

This module provides pure functions for calculating various biophysical
and biochemical properties of DNA sequences that affect probe performance.

Available calculators:
- GC content: Percentage of G and C bases
- Melting temperature (Tm): Thermodynamic stability estimation
- Complexity (DUST): Low-complexity region detection
- Hairpin: Self-complementarity score
- Dimer: Inter-probe complementarity score
- Entropy: Shannon entropy (placeholder for future implementation)

All functions are designed to be pure (no side effects) and thread-safe.
"""

from eprobe.biophysics.gc import calculate_gc, calculate_gc_batch
from eprobe.biophysics.tm import calculate_tm, calculate_tm_batch, TmMethod, TmTable
from eprobe.biophysics.complexity import calculate_dust_score, calculate_dust_batch
from eprobe.biophysics.hairpin import calculate_hairpin_score, calculate_hairpin_batch
from eprobe.biophysics.dimer import calculate_dimer_score, DimerCalculator
from eprobe.biophysics.entropy import calculate_entropy, calculate_entropy_batch
from eprobe.biophysics.thermo_entropy import (
    calculate_thermo_entropy,
    calculate_thermo_entropy_batch,
    calculate_thermo_enthalpy,
    calculate_gibbs_energy,
    NNTable,
)

# Aliases for convenience
calculate_complexity = calculate_dust_score
calculate_complexity_batch = calculate_dust_batch

__all__ = [
    # GC content
    "calculate_gc",
    "calculate_gc_batch",
    # Melting temperature
    "calculate_tm",
    "calculate_tm_batch",
    "TmMethod",
    "TmTable",
    # Complexity (DUST algorithm)
    "calculate_dust_score",
    "calculate_dust_batch",
    "calculate_complexity",  # alias
    "calculate_complexity_batch",  # alias
    # Hairpin
    "calculate_hairpin_score",
    "calculate_hairpin_batch",
    # Dimer
    "calculate_dimer_score",
    "DimerCalculator",
    # Shannon entropy
    "calculate_entropy",
    "calculate_entropy_batch",
    # Thermodynamic entropy (ΔS°, ΔH°, ΔG°)
    "calculate_thermo_entropy",
    "calculate_thermo_entropy_batch",
    "calculate_thermo_enthalpy",
    "calculate_gibbs_energy",
    "NNTable",
]
