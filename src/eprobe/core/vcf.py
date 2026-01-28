"""
VCF file operations.

Provides functions for reading and processing VCF files to extract SNPs.
Uses pysam for efficient indexed access to compressed VCF files.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional
import warnings

import pysam
import pandas as pd

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNP, SNPDataFrame


@dataclass
class VCFReader:
    """
    Reader for VCF files with SNP extraction capabilities.
    
    Handles compressed (vcf.gz) and indexed VCF files using pysam.
    Extracts only biallelic SNPs by default.
    """
    
    path: Path
    _vcf: Optional[pysam.VariantFile] = None
    
    def __post_init__(self) -> None:
        """Convert path to Path object."""
        self.path = Path(self.path)
    
    def open(self) -> Result[None, str]:
        """
        Open VCF file for reading.
        
        Returns:
            Ok(None) on success, Err(message) on failure
        """
        try:
            if not self.path.exists():
                return Err(f"VCF file not found: {self.path}")
            
            # Check for index file
            index_path = Path(f"{self.path}.tbi")
            if not index_path.exists():
                return Err(f"VCF index not found: {index_path}. Run 'tabix -p vcf {self.path}'")
            
            self._vcf = pysam.VariantFile(str(self.path))
            return Ok(None)
        
        except Exception as e:
            return Err(f"Failed to open VCF: {e}")
    
    def close(self) -> None:
        """Close VCF file."""
        if self._vcf is not None:
            self._vcf.close()
            self._vcf = None
    
    def __enter__(self) -> VCFReader:
        result = self.open()
        if result.is_err():
            raise IOError(result.unwrap_err())
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()
    
    @staticmethod
    def is_biallelic_snp(record: pysam.VariantRecord) -> bool:
        """
        Check if variant record is a biallelic SNP.
        
        Args:
            record: pysam VariantRecord
            
        Returns:
            True if biallelic SNP (single-base ref and alt)
        """
        # Check reference is single base
        if len(record.ref) != 1:
            return False
        
        # Check exactly one alternative allele
        if len(record.alts) != 1:
            return False
        
        # Check alternative is single base
        alt = record.alts[0]
        if alt not in {"A", "C", "G", "T"}:
            return False
        
        return True
    
    def fetch_snps(
        self,
        chrom: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
    ) -> Iterator[SNP]:
        """
        Fetch SNPs from VCF file.
        
        Args:
            chrom: Chromosome to fetch from (None for all)
            start: Start position (1-based, optional)
            end: End position (1-based, optional)
            
        Yields:
            SNP objects for each biallelic SNP
        """
        if self._vcf is None:
            raise RuntimeError("VCF file not opened. Call open() first.")
        
        # Determine fetch region
        if chrom is not None:
            records = self._vcf.fetch(chrom, start, end)
        else:
            records = self._vcf.fetch()
        
        for record in records:
            if self.is_biallelic_snp(record):
                yield SNP.from_vcf_record(
                    chrom=record.chrom,
                    pos=record.pos,
                    ref=record.ref,
                    alt=record.alts[0],
                )
    
    def extract_snps_by_chromosome(
        self,
        chrom: str,
    ) -> Result[list[SNP], str]:
        """
        Extract all SNPs from a specific chromosome.
        
        Args:
            chrom: Chromosome name
            
        Returns:
            Ok(list of SNPs) on success, Err(message) on failure
        """
        try:
            snps = list(self.fetch_snps(chrom=chrom))
            return Ok(snps)
        except Exception as e:
            return Err(f"Failed to extract SNPs from {chrom}: {e}")
    
    @property
    def chromosomes(self) -> list[str]:
        """Get list of chromosomes in VCF."""
        if self._vcf is None:
            raise RuntimeError("VCF file not opened")
        return list(self._vcf.header.contigs)


def extract_snps_from_vcf(
    vcf_path: Path | str,
    chromosomes: Optional[list[str]] = None,
) -> Result[SNPDataFrame, str]:
    """
    Extract all biallelic SNPs from VCF file.
    
    Args:
        vcf_path: Path to VCF file (vcf.gz)
        chromosomes: Optional list of chromosomes to extract from
        
    Returns:
        Ok(SNPDataFrame) on success, Err(message) on failure
        
    Note:
        For large VCF files, consider using multiprocessing with
        extract_snps_parallel() instead.
    """
    try:
        reader = VCFReader(Path(vcf_path))
        result = reader.open()
        if result.is_err():
            return Err(result.unwrap_err())
        
        all_snps: list[SNP] = []
        target_chroms = chromosomes or reader.chromosomes
        
        non_biallelic_count = 0
        
        for chrom in target_chroms:
            for snp in reader.fetch_snps(chrom=chrom):
                all_snps.append(snp)
        
        reader.close()
        
        if non_biallelic_count > 0:
            warnings.warn(f"Skipped {non_biallelic_count} non-biallelic variants")
        
        return Ok(SNPDataFrame.from_snps(all_snps))
    
    except Exception as e:
        return Err(f"Failed to extract SNPs: {e}")


def process_chromosome_snps(args: tuple[str, str]) -> Result[list[dict], str]:
    """
    Process SNPs for a single chromosome (for multiprocessing).
    
    Args:
        args: Tuple of (chromosome, vcf_path)
        
    Returns:
        Ok(list of SNP dicts) on success, Err(message) on failure
    """
    chrom, vcf_path = args
    
    try:
        snp_dicts: list[dict] = []
        
        with pysam.VariantFile(vcf_path) as vcf:
            for record in vcf.fetch(chrom):
                if VCFReader.is_biallelic_snp(record):
                    snp = SNP.from_vcf_record(
                        chrom=record.chrom,
                        pos=record.pos,
                        ref=record.ref,
                        alt=record.alts[0],
                    )
                    snp_dicts.append(snp.to_dict())
        
        return Ok(snp_dicts)
    
    except Exception as e:
        return Err(f"Failed to process chromosome {chrom}: {e}")
