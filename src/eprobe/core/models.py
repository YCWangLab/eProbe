"""
Core data models for eProbe.

Defines immutable data classes for SNPs, Probes, and related structures.
All models use dataclasses with slots for memory efficiency.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Iterator, Literal
from collections import OrderedDict
from pathlib import Path
import pandas as pd

from eprobe.core.result import Result, Ok, Err


# Type definitions
MutationType = Literal["ts", "tv"]  # transition or transversion


@dataclass(frozen=True, slots=True)
class SNP:
    """
    Immutable representation of a Single Nucleotide Polymorphism.
    
    Attributes:
        chrom: Chromosome name
        pos: 1-based position on chromosome
        ref: Reference allele (A/T/C/G)
        alt: Alternative allele (A/T/C/G)
        mutation_type: 'ts' for transition, 'tv' for transversion
    """
    chrom: str
    pos: int
    ref: str
    alt: str
    mutation_type: MutationType
    
    def __post_init__(self) -> None:
        """Validate SNP data on creation."""
        valid_bases = {"A", "T", "C", "G"}
        if self.ref not in valid_bases:
            raise ValueError(f"Invalid reference allele: {self.ref}")
        if self.alt not in valid_bases:
            raise ValueError(f"Invalid alternative allele: {self.alt}")
        if self.pos < 1:
            raise ValueError(f"Position must be >= 1, got {self.pos}")
    
    @staticmethod
    def determine_mutation_type(ref: str, alt: str) -> MutationType:
        """
        Determine if mutation is transition (ts) or transversion (tv).
        
        Transitions: A<->G, C<->T (purine<->purine or pyrimidine<->pyrimidine)
        Transversions: All other substitutions
        """
        transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
        return "ts" if (ref, alt) in transitions else "tv"
    
    @classmethod
    def from_vcf_record(cls, chrom: str, pos: int, ref: str, alt: str) -> SNP:
        """Create SNP from VCF-style record."""
        mutation_type = cls.determine_mutation_type(ref, alt)
        return cls(chrom=chrom, pos=pos, ref=ref, alt=alt, mutation_type=mutation_type)
    
    def to_dict(self) -> dict:
        """Convert to dictionary for DataFrame creation."""
        return {
            "chr": self.chrom,
            "pos": self.pos,
            "type": self.mutation_type,
            "ref": self.ref,
            "alt": self.alt,
        }
    
    @property
    def id(self) -> str:
        """Generate unique identifier for this SNP.
        
        Format: {chrom}:{pos}_{ref}_to_{alt}
        Note: Uses '_to_' instead of '>' to avoid conflicts with FASTA format,
        where '>' is the sequence header delimiter.
        """
        return f"{self.chrom}:{self.pos}_{self.ref}_to_{self.alt}"


@dataclass
class SNPDataFrame:
    """
    Container for SNP data with DataFrame backend.
    
    Provides a pandas DataFrame interface with domain-specific methods.
    The internal DataFrame has columns: chr, pos, type, ref, alt
    Additional columns may be added for biophysical tags.
    """
    
    _df: pd.DataFrame
    
    # Required columns in SNP DataFrame
    REQUIRED_COLUMNS = ["chr", "pos", "type", "ref", "alt"]
    
    # Column order for consistent output
    COLUMN_ORDER = ["chr", "start", "end", "pos", "type", "ref", "alt", 
                    "gc", "tm", "complexity", "hairpin", "dimer", "entropy"]
    
    def __init__(self, df: pd.DataFrame | None = None) -> None:
        """
        Initialize SNPDataFrame.
        
        Args:
            df: Optional pandas DataFrame. If None, creates empty DataFrame.
        """
        if df is None:
            self._df = pd.DataFrame(columns=self.REQUIRED_COLUMNS)
        else:
            self._df = df.copy()
    
    @classmethod
    def from_snps(cls, snps: list[SNP]) -> SNPDataFrame:
        """Create SNPDataFrame from list of SNP objects."""
        if not snps:
            return cls()
        data = [snp.to_dict() for snp in snps]
        return cls(pd.DataFrame(data))
    
    @classmethod
    def from_tsv(cls, path: Path | str, has_header: bool = True) -> Result[SNPDataFrame, str]:
        """
        Load SNPDataFrame from TSV file.
        
        Args:
            path: Path to TSV file
            has_header: Whether file has header row
            
        Returns:
            Ok(SNPDataFrame) on success, Err(message) on failure
        """
        try:
            path = Path(path)
            if not path.exists():
                return Err(f"File not found: {path}")
            
            if has_header:
                df = pd.read_csv(path, sep="\t", header=0)
            else:
                df = pd.read_csv(path, sep="\t", header=None, 
                                names=cls.REQUIRED_COLUMNS)
            
            # Validate required columns
            missing = set(cls.REQUIRED_COLUMNS) - set(df.columns)
            if missing:
                return Err(f"Missing required columns: {missing}")
            
            return Ok(cls(df))
        except Exception as e:
            return Err(f"Failed to read TSV: {e}")
    
    def to_tsv(self, path: Path | str, include_header: bool = True) -> Result[Path, str]:
        """
        Save SNPDataFrame to TSV file.
        
        Args:
            path: Output path
            include_header: Whether to include header row
            
        Returns:
            Ok(path) on success, Err(message) on failure
        """
        try:
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            
            # Reorder columns to standard order, keeping only existing columns
            cols = [c for c in self.COLUMN_ORDER if c in self._df.columns]
            # Add any extra columns not in standard order
            extra_cols = [c for c in self._df.columns if c not in self.COLUMN_ORDER]
            cols.extend(extra_cols)
            
            self._df[cols].to_csv(path, sep="\t", index=False, header=include_header)
            return Ok(path)
        except Exception as e:
            return Err(f"Failed to write TSV: {e}")
    
    def __len__(self) -> int:
        return len(self._df)
    
    def __iter__(self) -> Iterator[pd.Series]:
        return (row for _, row in self._df.iterrows())
    
    @property
    def df(self) -> pd.DataFrame:
        """Access underlying DataFrame (read-only copy)."""
        return self._df.copy()
    
    @property
    def chromosomes(self) -> list[str]:
        """Get list of unique chromosomes."""
        return self._df["chr"].unique().tolist()
    
    def filter_by_chromosome(self, chroms: list[str]) -> SNPDataFrame:
        """Return new SNPDataFrame with only specified chromosomes."""
        filtered = self._df[self._df["chr"].isin(chroms)]
        return SNPDataFrame(filtered)
    
    def filter_by_position(self, chrom: str, start: int, end: int) -> SNPDataFrame:
        """Return SNPs within genomic region."""
        mask = (self._df["chr"] == chrom) & \
               (self._df["pos"] >= start) & \
               (self._df["pos"] <= end)
        return SNPDataFrame(self._df[mask])
    
    def filter_by_mutation_type(self, mutation_type: MutationType) -> SNPDataFrame:
        """Return SNPs of specified mutation type."""
        filtered = self._df[self._df["type"] == mutation_type]
        return SNPDataFrame(filtered)
    
    def add_column(self, name: str, values: pd.Series | list) -> SNPDataFrame:
        """
        Add or update a column (returns new SNPDataFrame, immutable pattern).
        
        Args:
            name: Column name
            values: Column values (must match DataFrame length)
            
        Returns:
            New SNPDataFrame with added column
        """
        new_df = self._df.copy()
        new_df[name] = values
        return SNPDataFrame(new_df)
    
    def sample(self, n: int, seed: int = 42) -> SNPDataFrame:
        """Return random sample of n SNPs."""
        if n >= len(self._df):
            return SNPDataFrame(self._df.copy())
        sampled = self._df.sample(n=n, random_state=seed)
        return SNPDataFrame(sampled.sort_values(["chr", "pos"]))
    
    def sort(self) -> SNPDataFrame:
        """Return sorted SNPDataFrame by chromosome and position."""
        sorted_df = self._df.sort_values(["chr", "pos"]).reset_index(drop=True)
        return SNPDataFrame(sorted_df)
    
    def to_snps(self) -> list[SNP]:
        """Convert to list of SNP objects."""
        snps = []
        for _, row in self._df.iterrows():
            snp = SNP(
                chrom=row["chr"],
                pos=int(row["pos"]),
                ref=row["ref"],
                alt=row["alt"],
                mutation_type=row["type"]
            )
            snps.append(snp)
        return snps


@dataclass(frozen=True, slots=True)
class Probe:
    """
    Immutable representation of a capture probe sequence.
    
    Attributes:
        id: Unique identifier for the probe
        sequence: DNA sequence (uppercase)
        source_chrom: Source chromosome (for SNP-based probes)
        source_start: Start position on source (1-based)
        source_end: End position on source (1-based)
        snp_pos: SNP position if SNP-based probe
        snp_ref: Reference allele if SNP-based
        snp_alt: Alternative allele if SNP-based
    """
    id: str
    sequence: str
    source_chrom: Optional[str] = None
    source_start: Optional[int] = None
    source_end: Optional[int] = None
    snp_pos: Optional[int] = None
    snp_ref: Optional[str] = None
    snp_alt: Optional[str] = None
    mutation_type: Optional[MutationType] = None
    
    def __post_init__(self) -> None:
        """Validate and normalize sequence."""
        # Use object.__setattr__ for frozen dataclass
        object.__setattr__(self, "sequence", self.sequence.upper())
        
        valid_bases = set("ATCGN")
        if not all(b in valid_bases for b in self.sequence):
            invalid = set(self.sequence) - valid_bases
            raise ValueError(f"Invalid bases in sequence: {invalid}")
    
    @property
    def length(self) -> int:
        return len(self.sequence)
    
    @property
    def gc_content(self) -> float:
        """Calculate GC content as percentage."""
        if not self.sequence:
            return 0.0
        gc_count = self.sequence.count("G") + self.sequence.count("C")
        return (gc_count / len(self.sequence)) * 100
    
    def to_fasta_record(self) -> str:
        """Format as FASTA record."""
        return f">{self.id}\n{self.sequence}"
    
    @classmethod
    def from_snp_region(
        cls,
        snp: SNP,
        sequence: str,
        start: int,
        end: int,
    ) -> Probe:
        """
        Create probe from SNP and surrounding sequence.
        
        Args:
            snp: Source SNP
            sequence: Probe sequence
            start: Start position (1-based)
            end: End position (1-based)
        """
        probe_id = f"{snp.chrom}:{start}-{end}_{snp.pos}_{snp.mutation_type}_{snp.ref}_{snp.alt}"
        return cls(
            id=probe_id,
            sequence=sequence,
            source_chrom=snp.chrom,
            source_start=start,
            source_end=end,
            snp_pos=snp.pos,
            snp_ref=snp.ref,
            snp_alt=snp.alt,
            mutation_type=snp.mutation_type,
        )


@dataclass
class ProbeSet:
    """
    Collection of probes with associated metadata.
    
    Uses OrderedDict to maintain insertion order and enable
    efficient lookup by probe ID.
    """
    
    _probes: OrderedDict[str, Probe] = field(default_factory=OrderedDict)
    name: str = "probe_set"
    
    def add(self, probe: Probe) -> None:
        """Add a probe to the set."""
        self._probes[probe.id] = probe
    
    def get(self, probe_id: str) -> Optional[Probe]:
        """Get probe by ID."""
        return self._probes.get(probe_id)
    
    def __len__(self) -> int:
        return len(self._probes)
    
    def __iter__(self) -> Iterator[Probe]:
        return iter(self._probes.values())
    
    def __contains__(self, probe_id: str) -> bool:
        return probe_id in self._probes
    
    @property
    def ids(self) -> list[str]:
        """Get all probe IDs."""
        return list(self._probes.keys())
    
    @property
    def sequences(self) -> list[str]:
        """Get all probe sequences."""
        return [p.sequence for p in self._probes.values()]
    
    def to_fasta(self, path: Path | str) -> Result[Path, str]:
        """
        Write probes to FASTA file.
        
        Args:
            path: Output file path
            
        Returns:
            Ok(path) on success, Err(message) on failure
        """
        try:
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            
            with open(path, "w") as f:
                for probe in self._probes.values():
                    f.write(f">{probe.id}\n{probe.sequence}\n")
            
            return Ok(path)
        except Exception as e:
            return Err(f"Failed to write FASTA: {e}")
    
    @classmethod
    def from_fasta(cls, path: Path | str, name: Optional[str] = None) -> Result[ProbeSet, str]:
        """
        Load ProbeSet from FASTA file.
        
        Args:
            path: Path to FASTA file
            name: Optional name for the probe set
            
        Returns:
            Ok(ProbeSet) on success, Err(message) on failure
        """
        try:
            from Bio import SeqIO
            
            path = Path(path)
            if not path.exists():
                return Err(f"File not found: {path}")
            
            probe_set = cls(name=name or path.stem)
            
            for record in SeqIO.parse(path, "fasta"):
                probe = Probe(
                    id=record.id,
                    sequence=str(record.seq),
                )
                probe_set.add(probe)
            
            return Ok(probe_set)
        except Exception as e:
            return Err(f"Failed to read FASTA: {e}")
    
    def filter_by_gc(self, min_gc: float, max_gc: float) -> ProbeSet:
        """Return new ProbeSet with probes in GC range."""
        filtered = ProbeSet(name=f"{self.name}_gc_filtered")
        for probe in self._probes.values():
            if min_gc <= probe.gc_content <= max_gc:
                filtered.add(probe)
        return filtered
    
    def deduplicate(self) -> ProbeSet:
        """
        Remove probes with duplicate sequences.
        Keeps first occurrence of each unique sequence.
        """
        seen_seqs: set[str] = set()
        deduped = ProbeSet(name=f"{self.name}_dedup")
        
        for probe in self._probes.values():
            if probe.sequence not in seen_seqs:
                seen_seqs.add(probe.sequence)
                deduped.add(probe)
        
        return deduped
    
    def merge(self, other: ProbeSet) -> ProbeSet:
        """
        Merge with another ProbeSet.
        
        On ID collision, keeps probe from self.
        """
        merged = ProbeSet(name=f"{self.name}_merged")
        
        for probe in self._probes.values():
            merged.add(probe)
        
        for probe in other._probes.values():
            if probe.id not in merged:
                merged.add(probe)
        
        return merged
