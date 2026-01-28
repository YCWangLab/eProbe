"""
FUNCGEN panel implementation modules.

This package contains the implementation for the FUNCGEN panel commands:
  from_fasta - Generate probes from FASTA sequences
  from_bed   - Generate probes from BED regions
  from_gff   - Generate probes from GFF annotations
  assess     - Quality assessment
"""

from eprobe.funcgen.from_fasta import run_from_fasta
from eprobe.funcgen.from_bed import run_from_bed
from eprobe.funcgen.from_gff import run_from_gff
from eprobe.funcgen.assess import run_assess

__all__ = [
    "run_from_fasta",
    "run_from_bed",
    "run_from_gff",
    "run_assess",
]
