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
    
    @staticmethod
    def _can_force_biallelic(record: pysam.VariantRecord) -> bool:
        """
        Quick check if variant can be forced to biallelic.
        
        Checks:
        1. REF is single base
        2. At least one ALT exists
        3. First ALT is single base
        
        Args:
            record: pysam VariantRecord
            
        Returns:
            True if can be forced to biallelic
        """
        # Check ref is single base
        if len(record.ref) != 1:
            return False
        
        # Check has at least one alt
        if not record.alts or len(record.alts) == 0:
            return False
        
        # Check first alt is single base and valid
        first_alt = record.alts[0]
        if first_alt in {"A", "C", "G", "T"}:
            return True
        
        return False
    
    @staticmethod
    def _record_to_snp(record: pysam.VariantRecord, force_first_alt: bool = False) -> SNP:
        """
        Convert pysam VariantRecord to SNP object.
        
        Args:
            record: pysam VariantRecord
            force_first_alt: If True, use first alt even if multiple exist
            
        Returns:
            SNP object
        """
        alt_allele = record.alts[0]
        
        return SNP.from_vcf_record(
            chrom=record.chrom,
            pos=record.pos,
            ref=record.ref,
            alt=alt_allele,
            snp_id=record.id if record.id else None,
        )
    
    def fetch_snps(
        self,
        chrom: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        bed_path: Optional[Path] = None,
        force_biallelic: bool = False,
    ) -> Iterator[SNP]:
        """
        Fetch biallelic SNPs from VCF.
        
        Args:
            chrom: Chromosome to fetch (None = all)
            start: Start position (0-based)
            end: End position
            bed_path: BED file to restrict regions
            force_biallelic: Force biallelic by taking first alt
            
        Yields:
            SNP objects for biallelic variants
        """
        if self._vcf is None:
            raise RuntimeError("VCF not opened")
        
        # Read BED regions if provided
        regions = []
        if bed_path is not None:
            with open(bed_path) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        regions.append((parts[0], int(parts[1]), int(parts[2])))
        
        # Fetch variants
        if regions:
            # Fetch from specified regions
            for region_chrom, region_start, region_end in regions:
                for record in self._vcf.fetch(region_chrom, region_start, region_end):
                    if self.is_biallelic_snp(record):
                        yield self._record_to_snp(record)
                    elif force_biallelic and self._can_force_biallelic(record):
                        yield self._record_to_snp(record, force_first_alt=True)
        elif chrom is not None:
            # Fetch from specified chromosome/region
            for record in self._vcf.fetch(chrom, start, end):
                if self.is_biallelic_snp(record):
                    yield self._record_to_snp(record)
                elif force_biallelic and self._can_force_biallelic(record):
                    yield self._record_to_snp(record, force_first_alt=True)
        else:
            # Fetch all
            for record in self._vcf:
                if self.is_biallelic_snp(record):
                    yield self._record_to_snp(record)
                elif force_biallelic and self._can_force_biallelic(record):
                    yield self._record_to_snp(record, force_first_alt=True)
    
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
    vcf_path: Path,
    bed_path: Optional[Path] = None,
    force_biallelic: bool = False,
) -> Result[list[SNP], str]:
    """
    Extract biallelic SNPs from VCF file.
    
    Args:
        vcf_path: Path to VCF file (.vcf.gz with .tbi index)
        bed_path: Optional BED file to filter regions
        force_biallelic: If True, take first alt allele when multiple exist
        
    Returns:
        Result containing list of SNPs
    """
    reader = VCFReader(vcf_path)
    open_result = reader.open()
    if open_result.is_err():
        return Err(open_result.unwrap_err())
    
    try:
        snps = list(reader.fetch_snps(bed_path=bed_path, force_biallelic=force_biallelic))
        reader.close()
        return Ok(snps)
    except Exception as e:
        reader.close()
        return Err(f"SNP extraction failed: {e}")
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
