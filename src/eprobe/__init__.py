"""
eProbe: A one-stop capture probe design toolkit for genetic diversity
reconstructions from ancient environmental DNA.

This package provides tools for designing capture probes for both:
- POPGEN panel: SNP-based probes for population genetics
- FUNCGEN panel: Sequence-based probes for functional genetics
"""

__version__ = "1.0.0"
__author__ = "Zi-Hao Huang"
__email__ = "zh384@cam.ac.uk"

from eprobe.core.result import Result, Ok, Err

__all__ = [
    "__version__",
    "Result",
    "Ok",
    "Err",
]
