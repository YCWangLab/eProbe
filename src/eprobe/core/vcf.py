"""
VCF file operations.

Provides functions for reading and processing VCF files to extract SNPs.
Uses pysam for efficient indexed access to compressed VCF files.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional, List, Tuple
import warnings
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

import pysam
import pandas as pd

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNP, SNPDataFrame

logger = logging.getLogger(__name__)


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


def _extract_snps_for_chromosome(
    vcf_path: Path,
    chrom: str,
    bed_path: Optional[Path],
    force_biallelic: bool,
) -> Tuple[str, List[SNP], int]:
    """
    Extract SNPs from a single chromosome (for multiprocessing).
    
    Args:
        vcf_path: Path to VCF file
        chrom: Chromosome name
        bed_path: Optional BED file to filter regions
        force_biallelic: Force biallelic mode
        
    Returns:
        Tuple of (chromosome, snp_list, non_biallelic_count)
    """
    logger.info(f"→ Processing chromosome {chrom}...")
    
    reader = VCFReader(vcf_path)
    open_result = reader.open()
    if open_result.is_err():
        logger.error(f"Failed to open VCF for {chrom}: {open_result.unwrap_err()}")
        return (chrom, [], 0)
    
    snps = []
    non_biallelic = 0
    
    try:
        for snp in reader.fetch_snps(chrom=chrom, bed_path=bed_path, force_biallelic=force_biallelic):
            snps.append(snp)
        
        # Count non-biallelic variants if not forcing
        if not force_biallelic:
            for record in reader._vcf.fetch(chrom):
                if not reader.is_biallelic_snp(record):
                    non_biallelic += 1
        
        reader.close()
        logger.info(f"✓ Chromosome {chrom}: {len(snps)} SNPs extracted")
        
        return (chrom, snps, non_biallelic)
        
    except Exception as e:
        logger.error(f"Error processing {chrom}: {e}")
        reader.close()
        return (chrom, [], 0)


def extract_snps_from_vcf(
    vcf_path: Path,
    bed_path: Optional[Path] = None,
    force_biallelic: bool = False,
    threads: int = 1,
) -> Result[list[SNP], str]:
    """
    Extract biallelic SNPs from VCF file (with multiprocessing).
    
    Args:
        vcf_path: Path to VCF file (.vcf.gz with .tbi index)
        bed_path: Optional BED file to filter regions
        force_biallelic: If True, take first alt allele when multiple exist
        threads: Number of parallel processes (default: 1)
        
    Returns:
        Result containing list of SNPs
    """
    # Get chromosome list from VCF header
    try:
        reader = VCFReader(vcf_path)
        open_result = reader.open()
        if open_result.is_err():
            return Err(open_result.unwrap_err())
        
        chromosomes = reader.chromosomes
        reader.close()
        
        logger.info(f"Found {len(chromosomes)} chromosomes in VCF")
        logger.info(f"Using {threads} threads for parallel extraction")
        
    except Exception as e:
        return Err(f"Failed to read VCF header: {e}")
    
    all_snps = []
    total_non_biallelic = 0
    
    if threads == 1:
        # Serial processing
        for chrom in chromosomes:
            chrom, snps, non_biallelic = _extract_snps_for_chromosome(
                vcf_path, chrom, bed_path, force_biallelic
            )
            all_snps.extend(snps)
            total_non_biallelic += non_biallelic
    else:
        # Parallel processing by chromosome
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {}
            for chrom in chromosomes:
                future = executor.submit(
                    _extract_snps_for_chromosome,
                    vcf_path, chrom, bed_path, force_biallelic
                )
                futures[future] = chrom
            
            for future in as_completed(futures):
                chrom_name = futures[future]
                try:
                    chrom, snps, non_biallelic = future.result()
                    all_snps.extend(snps)
                    total_non_biallelic += non_biallelic
                except Exception as e:
                    logger.error(f"Failed to process {chrom_name}: {e}")
    
    if total_non_biallelic > 0 and not force_biallelic:
        logger.warning(f"Skipped {total_non_biallelic} non-biallelic variants")
    
    return Ok(all_snps)
