"""
Core module for eProbe.

Contains fundamental data structures, result types, and shared utilities.
"""

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNP, SNPDataFrame, Probe, ProbeSet
from eprobe.core.fasta import read_fasta, write_fasta, FastaDict

__all__ = [
    "Result",
    "Ok", 
    "Err",
    "SNP",
    "SNPDataFrame",
    "Probe",
    "ProbeSet",
    "read_fasta",
    "write_fasta",
    "FastaDict",
]
