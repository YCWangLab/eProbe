"""
Probe sequence generation module.

Generates final probe sequences from selected SNPs:
  - Extracts sequences with configurable length and offset
  - Supports alternate allele probe generation
  - Handles edge cases near chromosome boundaries

This module corresponds to the original SNP_generator.py functionality.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNP, SNPDataFrame, Probe, ProbeSet
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


@dataclass
class ProbeConfig:
    """Configuration for probe generation."""
    length: int = 81
    shift: int = 0  # SNP offset from center (0 = centered)
    replace_mode: bool = True  # Generate alternate allele probes
    include_ref: bool = True   # Include reference allele probes
    include_alt: bool = True   # Include alternate allele probes


def calculate_probe_coordinates(
    snp_pos: int,
    probe_length: int,
    shift: int = 0,
) -> Tuple[int, int, int]:
    """
    Calculate probe start, end, and SNP offset within probe.
    
    Args:
        snp_pos: 1-based SNP position
        probe_length: Desired probe length
        shift: Offset from center (positive = shift right, SNP moves left in probe)
        
    Returns:
        Tuple of (start_pos, end_pos, snp_offset)
        - start_pos: 1-based start position
        - end_pos: 1-based end position (exclusive)
        - snp_offset: 0-based offset of SNP within probe
    """
    # Center position in probe (0-based)
    center_offset = (probe_length - 1) // 2
    
    # Apply shift: positive shift moves SNP left (toward start)
    snp_offset = center_offset - shift
    
    # Ensure SNP offset is within probe
    snp_offset = max(0, min(probe_length - 1, snp_offset))
    
    # Calculate genomic positions (1-based)
    start_pos = snp_pos - snp_offset
    end_pos = start_pos + probe_length
    
    return start_pos, end_pos, snp_offset


def generate_probe_sequence(
    snp: SNP,
    reference_seq: str,
    chrom_length: int,
    config: ProbeConfig,
) -> Result[List[Probe], str]:
    """
    Generate probe sequences for a single SNP.
    
    Args:
        snp: SNP with position and alleles
        reference_seq: Reference chromosome sequence
        chrom_length: Length of chromosome
        config: Probe generation configuration
        
    Returns:
        Result containing list of Probe objects
    """
    start, end, snp_offset = calculate_probe_coordinates(
        snp.pos, config.length, config.shift
    )
    
    # Check boundaries
    if start < 1:
        return Err(f"Probe start ({start}) before chromosome start for {snp.snp_id}")
    if end > chrom_length + 1:
        return Err(f"Probe end ({end}) beyond chromosome end for {snp.snp_id}")
    
    # Extract reference sequence (convert to 0-based for Python string)
    probe_ref_seq = reference_seq[start - 1:end - 1]
    
    # Validate sequence length
    if len(probe_ref_seq) != config.length:
        return Err(f"Extracted sequence length mismatch for {snp.snp_id}")
    
    # Validate SNP position matches reference
    ref_at_snp = probe_ref_seq[snp_offset]
    if ref_at_snp.upper() != snp.ref.upper():
        return Err(
            f"Reference mismatch at {snp.snp_id}: "
            f"expected {snp.ref}, found {ref_at_snp}"
        )
    
    probes: List[Probe] = []
    
    # Generate reference allele probe
    if config.include_ref:
        ref_probe = Probe(
            probe_id=f"{snp.snp_id}_REF",
            sequence=probe_ref_seq.upper(),
            snp_id=snp.snp_id,
            chrom=snp.chrom,
            start=start,
            end=end - 1,
            snp_offset=snp_offset,
            allele=snp.ref,
            is_ref=True,
        )
        probes.append(ref_probe)
    
    # Generate alternate allele probe
    if config.include_alt and config.replace_mode:
        # Replace reference allele with alternate
        seq_list = list(probe_ref_seq.upper())
        seq_list[snp_offset] = snp.alt.upper()
        probe_alt_seq = "".join(seq_list)
        
        alt_probe = Probe(
            probe_id=f"{snp.snp_id}_ALT",
            sequence=probe_alt_seq,
            snp_id=snp.snp_id,
            chrom=snp.chrom,
            start=start,
            end=end - 1,
            snp_offset=snp_offset,
            allele=snp.alt,
            is_ref=False,
        )
        probes.append(alt_probe)
    
    return Ok(probes)


def generate_probes_batch(
    snps: List[SNP],
    reference_sequences: Dict[str, str],
    config: ProbeConfig,
) -> Result[ProbeSet, str]:
    """
    Generate probes for a batch of SNPs.
    
    Args:
        snps: List of SNPs
        reference_sequences: Reference sequences by chromosome
        config: Probe generation configuration
        
    Returns:
        Result containing ProbeSet with all generated probes
    """
    all_probes: List[Probe] = []
    errors: List[str] = []
    
    for snp in snps:
        ref_seq = reference_sequences.get(snp.chrom)
        if ref_seq is None:
            errors.append(f"Chromosome {snp.chrom} not in reference")
            continue
        
        result = generate_probe_sequence(snp, ref_seq, len(ref_seq), config)
        
        if result.is_err():
            errors.append(result.unwrap_err())
            logger.debug(f"Failed to generate probe: {result.unwrap_err()}")
            continue
        
        all_probes.extend(result.unwrap())
    
    if errors:
        logger.warning(f"Failed to generate {len(errors)} probes")
        for err in errors[:5]:  # Show first 5 errors
            logger.debug(f"  {err}")
        if len(errors) > 5:
            logger.debug(f"  ... and {len(errors) - 5} more")
    
    probe_set = ProbeSet(probes=all_probes)
    return Ok(probe_set)


def run_build(
    input_path: Path,
    reference_path: Path,
    output_prefix: Path,
    length: int = 81,
    shift: int = 0,
    replace_mode: bool = True,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Build probe sequences from selected SNPs.
    
    Main entry point for the build command. Generates final probe
    sequences ready for synthesis.
    
    Args:
        input_path: Input SNP TSV file
        reference_path: Reference genome FASTA
        output_prefix: Output file prefix
        length: Probe length
        shift: SNP position offset from center
        replace_mode: Generate alternate allele probes
        verbose: Enable verbose logging
        
    Returns:
        Result containing build statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Starting probe generation from {input_path}")
    logger.info(f"Probe length: {length}, Shift: {shift}, Replace mode: {replace_mode}")
    
    # Load input SNPs
    snp_df = SNPDataFrame.from_tsv(input_path)
    if snp_df.is_err():
        return Err(f"Failed to load input: {snp_df.unwrap_err()}")
    
    snps = snp_df.unwrap().to_snps()
    logger.info(f"Loaded {len(snps)} SNPs")
    
    # Load reference genome
    logger.info(f"Loading reference genome from {reference_path}...")
    ref_result = read_fasta(reference_path)
    if ref_result.is_err():
        return Err(f"Failed to read reference: {ref_result.unwrap_err()}")
    
    reference_sequences = ref_result.unwrap()
    logger.info(f"Loaded {len(reference_sequences)} chromosomes")
    
    # Configure probe generation
    config = ProbeConfig(
        length=length,
        shift=shift,
        replace_mode=replace_mode,
    )
    
    # Generate probes
    logger.info("Generating probe sequences...")
    probe_result = generate_probes_batch(snps, reference_sequences, config)
    if probe_result.is_err():
        return Err(probe_result.unwrap_err())
    
    probe_set = probe_result.unwrap()
    
    # Save FASTA output
    fasta_path = Path(str(output_prefix) + ".probes.fa")
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    
    sequences = {probe.probe_id: probe.sequence for probe in probe_set.probes}
    write_result = write_fasta(sequences, fasta_path)
    if write_result.is_err():
        return Err(f"Failed to save FASTA: {write_result.unwrap_err()}")
    
    logger.info(f"Saved {len(sequences)} probes to {fasta_path}")
    
    # Save TSV metadata
    tsv_path = Path(str(output_prefix) + ".probes.tsv")
    probe_set.to_tsv(tsv_path)
    
    # Calculate statistics
    ref_probes = sum(1 for p in probe_set.probes if p.is_ref)
    alt_probes = sum(1 for p in probe_set.probes if not p.is_ref)
    
    stats = {
        "snp_count": len(snps),
        "probe_count": len(probe_set.probes),
        "ref_probes": ref_probes,
        "alt_probes": alt_probes,
        "probe_length": length,
        "shift": shift,
        "fasta_file": str(fasta_path),
        "tsv_file": str(tsv_path),
    }
    
    logger.info(f"Probe generation complete: {stats['probe_count']} probes from {stats['snp_count']} SNPs")
    
    return Ok(stats)
