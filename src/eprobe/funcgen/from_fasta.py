"""
Probe generation from FASTA sequences.

Generates capture probes by tiling across input sequences using
a sliding window approach. Supports haplotype-aware probe design.

This module corresponds to part of the original Seq_generator.py functionality.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass
from collections import defaultdict
import re

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta
from eprobe.core.models import Probe

logger = logging.getLogger(__name__)


@dataclass
class TilingConfig:
    """Configuration for tiling probe generation."""
    probe_length: int = 81
    step_size: int = 30
    min_seq_length: int = 0  # Minimum input sequence length (0 = use probe_length)


@dataclass
class HaplotypeConfig:
    """Configuration for haplotype-aware processing."""
    enabled: bool = False
    separator: str = "_"  # e.g., Gene_1, Gene_2
    min_freq: float = 0.05  # Minimum haplotype frequency


def parse_haplotype_id(
    seq_id: str,
    separator: str = "_",
) -> Tuple[str, Optional[str]]:
    """
    Parse sequence ID to extract gene and allele components.
    
    Args:
        seq_id: Sequence identifier (e.g., "BRCA1_1", "BRCA1_2")
        separator: Character separating gene from allele
        
    Returns:
        Tuple of (gene_id, allele_id or None)
    """
    # Try to split by separator from the right
    if separator in seq_id:
        parts = seq_id.rsplit(separator, 1)
        if len(parts) == 2 and parts[1].isdigit():
            return parts[0], parts[1]
    
    # No allele suffix found
    return seq_id, None


def tile_sequence(
    seq_id: str,
    sequence: str,
    config: TilingConfig,
) -> List[Probe]:
    """
    Generate tiled probes from a sequence using sliding window.
    
    Args:
        seq_id: Sequence identifier
        sequence: Input sequence
        config: Tiling configuration
        
    Returns:
        List of generated Probe objects
    """
    probes = []
    seq_len = len(sequence)
    
    # Minimum length check
    min_len = config.min_seq_length or config.probe_length
    if seq_len < min_len:
        logger.debug(f"Sequence {seq_id} too short ({seq_len} < {min_len}), skipping")
        return []
    
    # Calculate number of probes
    if seq_len <= config.probe_length:
        # Single probe if sequence equals or less than probe length
        n_probes = 1
    else:
        n_probes = (seq_len - config.probe_length) // config.step_size + 1
    
    for i in range(n_probes):
        start = i * config.step_size
        end = start + config.probe_length
        
        # Handle last probe - ensure we don't go past sequence end
        if end > seq_len:
            if seq_len >= config.probe_length:
                start = seq_len - config.probe_length
                end = seq_len
            else:
                continue
        
        probe_seq = sequence[start:end].upper()
        
        # Skip probes with N bases
        if 'N' in probe_seq:
            logger.debug(f"Skipping probe {seq_id}_{i+1}: contains N")
            continue
        
        probe = Probe(
            id=f"{seq_id}_P{i+1:04d}",
            sequence=probe_seq,
            source_chrom=seq_id,
            source_start=start + 1,  # 1-based
            source_end=end,
        )
        probes.append(probe)
    
    return probes


def group_haplotypes(
    sequences: Dict[str, str],
    config: HaplotypeConfig,
) -> Dict[str, Dict[str, str]]:
    """
    Group sequences by gene ID for haplotype processing.
    
    Args:
        sequences: All input sequences
        config: Haplotype configuration
        
    Returns:
        Dictionary of gene_id -> {allele_id: sequence}
    """
    groups: Dict[str, Dict[str, str]] = defaultdict(dict)
    
    for seq_id, sequence in sequences.items():
        gene_id, allele_id = parse_haplotype_id(seq_id, config.separator)
        
        if allele_id is None:
            allele_id = "1"  # Default allele
        
        groups[gene_id][allele_id] = sequence
    
    return dict(groups)


def process_haplotypes(
    groups: Dict[str, Dict[str, str]],
    tiling_config: TilingConfig,
    haplo_config: HaplotypeConfig,
) -> List[Probe]:
    """
    Process haplotype groups to generate non-redundant probes.
    
    For each gene, identifies variant positions across haplotypes
    and generates probes that capture the variation.
    
    Args:
        groups: Grouped sequences by gene
        tiling_config: Tiling configuration
        haplo_config: Haplotype configuration
        
    Returns:
        List of generated probes
    """
    all_probes = []
    
    for gene_id, alleles in groups.items():
        if len(alleles) == 1:
            # Single allele - just tile normally
            allele_id, seq = list(alleles.items())[0]
            probes = tile_sequence(gene_id, seq, tiling_config)
            all_probes.extend(probes)
        else:
            # Multiple alleles - process each uniquely
            logger.info(f"Gene {gene_id}: {len(alleles)} haplotypes")
            
            seen_probes: set = set()
            
            for allele_id, seq in alleles.items():
                seq_id = f"{gene_id}{haplo_config.separator}{allele_id}"
                probes = tile_sequence(seq_id, seq, tiling_config)
                
                # Deduplicate by sequence
                for probe in probes:
                    if probe.sequence not in seen_probes:
                        seen_probes.add(probe.sequence)
                        all_probes.append(probe)
    
    return all_probes


def run_from_fasta(
    fasta_path: Path,
    output_prefix: Path,
    probe_length: int = 81,
    step_size: int = 30,
    haplotype_sep: str = "_",
    min_freq: float = 0.05,
    haplotyping: bool = False,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Generate probes from FASTA sequences.
    
    Main entry point for the from_fasta command. Tiles probes across
    each input sequence using a sliding window.
    
    Args:
        fasta_path: Input FASTA file
        output_prefix: Output file prefix
        probe_length: Probe length
        step_size: Sliding window step
        haplotype_sep: Separator for haplotype IDs
        min_freq: Minimum haplotype frequency
        haplotyping: Enable haplotype-aware processing
        verbose: Enable verbose logging
        
    Returns:
        Result containing generation statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Starting probe generation from {fasta_path}")
    logger.info(f"Probe length: {probe_length}, Step: {step_size}")
    
    # Load input sequences
    fasta_result = read_fasta(fasta_path)
    if fasta_result.is_err():
        return Err(f"Failed to read FASTA: {fasta_result.unwrap_err()}")
    
    sequences = fasta_result.unwrap()
    logger.info(f"Loaded {len(sequences)} sequences")
    
    if not sequences:
        return Err("No sequences found in input FASTA")
    
    # Configure tiling
    tiling_config = TilingConfig(
        probe_length=probe_length,
        step_size=step_size,
    )
    
    # Generate probes
    if haplotyping:
        # Haplotype-aware processing
        haplo_config = HaplotypeConfig(
            enabled=True,
            separator=haplotype_sep,
            min_freq=min_freq,
        )
        
        groups = group_haplotypes(sequences, haplo_config)
        logger.info(f"Grouped into {len(groups)} genes")
        
        probes = process_haplotypes(groups, tiling_config, haplo_config)
    else:
        # Simple tiling
        probes = []
        for seq_id, sequence in sequences.items():
            seq_probes = tile_sequence(seq_id, sequence, tiling_config)
            probes.extend(seq_probes)
    
    logger.info(f"Generated {len(probes)} probes")
    
    # Save FASTA output
    fasta_output = Path(str(output_prefix) + ".probes.fa")
    fasta_output.parent.mkdir(parents=True, exist_ok=True)
    
    probe_sequences = {p.id: p.sequence for p in probes}
    write_result = write_fasta(probe_sequences, fasta_output)
    if write_result.is_err():
        return Err(f"Failed to save FASTA: {write_result.unwrap_err()}")
    
    logger.info(f"Saved probes to {fasta_output}")
    
    # Save TSV metadata
    tsv_output = Path(str(output_prefix) + ".probes.tsv")
    
    with open(tsv_output, 'w') as f:
        f.write("probe_id\tsequence\tgene_id\tstart\tend\tlength\n")
        for p in probes:
            f.write(f"{p.id}\t{p.sequence}\t{p.source_chrom}\t{p.source_start}\t{p.source_end}\t{len(p.sequence)}\n")
    
    stats = {
        "gene_count": len(sequences),
        "probe_count": len(probes),
        "probe_length": probe_length,
        "step_size": step_size,
        "haplotyping": haplotyping,
        "fasta_file": str(fasta_output),
        "tsv_file": str(tsv_output),
    }
    
    return Ok(stats)
