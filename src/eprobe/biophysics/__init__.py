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
- Entropy: Shannon entropy (optional, for advanced use)

All functions are designed to be pure (no side effects) and thread-safe.

Fast versions (recommended for large datasets):
- calculate_gc_fast, calculate_tm_fast, calculate_dust_fast
- calculate_hairpin_fast, calculate_hairpin_batch_fast
- DimerCalculatorFast, calculate_dimer_batch_fast
- calculate_all_stats_fast (unified batch calculation)
"""

from eprobe.biophysics.gc import calculate_gc, calculate_gc_batch
from eprobe.biophysics.tm import calculate_tm, calculate_tm_batch, TmMethod, TmTable
from eprobe.biophysics.complexity import calculate_dust_score, calculate_dust_batch
from eprobe.biophysics.hairpin import calculate_hairpin_score, calculate_hairpin_batch
from eprobe.biophysics.entropy import calculate_entropy, calculate_entropy_batch
from eprobe.biophysics.thermo_entropy import (
    calculate_thermo_entropy,
    calculate_thermo_entropy_batch,
    calculate_thermo_enthalpy,
    calculate_gibbs_energy,
    NNTable,
)

# Fast optimized versions
from eprobe.biophysics.biophysics import (
    calculate_gc_fast,
    calculate_gc_batch_fast,
    calculate_tm_fast,
    calculate_tm_batch_fast,
    NN_TABLE_OPTIONS,  # Available NN tables for Tm calculation
    calculate_dust_fast,
    calculate_dust_batch_fast,
    calculate_hairpin_fast,
    calculate_hairpin_batch_fast,
    DimerCalculatorFast,
    calculate_dimer_batch_fast,
    compute_biophysical_parallel,
    calculate_all_stats_fast,
    calculate_percentile_threshold,
    BiophysicalStats,
)

# Aliases for convenience
calculate_complexity = calculate_dust_score
calculate_complexity_batch = calculate_dust_batch

__all__ = [
    # GC content
    "calculate_gc",
    "calculate_gc_batch",
    "calculate_gc_fast",
    "calculate_gc_batch_fast",
    # Melting temperature
    "calculate_tm",
    "calculate_tm_batch",
    "calculate_tm_fast",
    "calculate_tm_batch_fast",
    "TmMethod",
    "TmTable",
    "NN_TABLE_OPTIONS",  # Available NN tables: DNA_NN1-4, R_DNA_NN1-4
    # Complexity (DUST algorithm)
    "calculate_dust_score",
    "calculate_dust_batch",
    "calculate_dust_fast",
    "calculate_dust_batch_fast",
    "calculate_complexity",  # alias
    "calculate_complexity_batch",  # alias
    # Hairpin
    "calculate_hairpin_score",
    "calculate_hairpin_batch",
    "calculate_hairpin_fast",
    "calculate_hairpin_batch_fast",
    # Dimer
    "calculate_dimer_score",
    "DimerCalculator",
    "DimerCalculatorFast",
    "calculate_dimer_batch_fast",
    # Shannon entropy (optional)
    "calculate_entropy",
    "calculate_entropy_batch",
    # Thermodynamic entropy (ΔS°, ΔH°, ΔG°) - optional
    "calculate_thermo_entropy",
    "calculate_thermo_entropy_batch",
    "calculate_thermo_enthalpy",
    "calculate_gibbs_energy",
    "NNTable",
    # Parallel computation
    "compute_biophysical_parallel",
    # Unified fast calculation
    "calculate_all_stats_fast",
    "calculate_percentile_threshold",
    "BiophysicalStats",
]
