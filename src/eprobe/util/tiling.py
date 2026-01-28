"""
Multi-position probe tiling.

Generates multiple probes per SNP with the SNP at different positions:
  - Center (default position)
  - Left-shifted (SNP toward 3' end)
  - Right-shifted (SNP toward 5' end)

This creates redundancy for more robust capture across different
hybridization conditions.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNP, SNPDataFrame, Probe
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


@dataclass
class TilePosition:
    """Definition of a tiling position."""
    name: str
    offset: int  # Offset from center (negative = SNP moves left, positive = SNP moves right)


def parse_positions(
    positions_str: str,
    default_offset: int = 20,
) -> Result[List[TilePosition], str]:
    """
    Parse position specification string.
    
    Supports two formats:
    1. Named: "center,left,right" (uses default_offset)
    2. Numeric: "-30,-15,0,+15,+30" (explicit offsets)
    
    Args:
        positions_str: Position specification
        default_offset: Offset for named positions
        
    Returns:
        Result containing list of TilePosition objects
    """
    positions = []
    parts = [p.strip() for p in positions_str.split(',')]
    
    for part in parts:
        part_lower = part.lower()
        
        if part_lower == "center" or part_lower == "c":
            positions.append(TilePosition(name="C", offset=0))
        elif part_lower == "left" or part_lower == "l":
            positions.append(TilePosition(name=f"L{default_offset}", offset=-default_offset))
        elif part_lower == "right" or part_lower == "r":
            positions.append(TilePosition(name=f"R{default_offset}", offset=default_offset))
        else:
            # Try to parse as numeric offset
            try:
                offset = int(part)
                if offset < 0:
                    name = f"L{abs(offset)}"
                elif offset > 0:
                    name = f"R{offset}"
                else:
                    name = "C"
                positions.append(TilePosition(name=name, offset=offset))
            except ValueError:
                return Err(f"Invalid position specification: {part}")
    
    if not positions:
        return Err("No valid positions specified")
    
    return Ok(positions)


def generate_tiled_probes_snp(
    snp: SNP,
    reference_seq: str,
    chrom_length: int,
    probe_length: int,
    positions: List[TilePosition],
    suffix_style: str = "position",
) -> List[Probe]:
    """
    Generate tiled probes for a single SNP.
    
    Args:
        snp: SNP with position information
        reference_seq: Reference chromosome sequence
        chrom_length: Length of chromosome
        probe_length: Probe length
        positions: List of tiling positions
        suffix_style: "position" (_L20, _C, _R20) or "index" (_1, _2, _3)
        
    Returns:
        List of generated Probe objects
    """
    probes = []
    center_offset = (probe_length - 1) // 2
    
    for idx, tile_pos in enumerate(positions):
        # Calculate SNP position in probe
        snp_in_probe = center_offset - tile_pos.offset
        
        # Ensure SNP is within probe
        if snp_in_probe < 0 or snp_in_probe >= probe_length:
            logger.debug(f"Position {tile_pos.name} puts SNP outside probe for {snp.snp_id}")
            continue
        
        # Calculate genomic coordinates (1-based)
        start = snp.pos - snp_in_probe
        end = start + probe_length
        
        # Check boundaries
        if start < 1 or end > chrom_length + 1:
            logger.debug(f"Position {tile_pos.name} exceeds boundaries for {snp.snp_id}")
            continue
        
        # Extract sequence (0-based indexing)
        probe_seq = reference_seq[start - 1:end - 1].upper()
        
        if len(probe_seq) != probe_length:
            continue
        
        # Generate suffix
        if suffix_style == "position":
            suffix = f"_{tile_pos.name}"
        else:
            suffix = f"_{idx + 1}"
        
        # Create reference allele probe
        probe = Probe(
            probe_id=f"{snp.snp_id}{suffix}_REF",
            sequence=probe_seq,
            snp_id=snp.snp_id,
            chrom=snp.chrom,
            start=start,
            end=end - 1,
            snp_offset=snp_in_probe,
            allele=snp.ref,
            is_ref=True,
        )
        probes.append(probe)
        
        # Create alternate allele probe
        seq_list = list(probe_seq)
        seq_list[snp_in_probe] = snp.alt.upper()
        alt_seq = "".join(seq_list)
        
        alt_probe = Probe(
            probe_id=f"{snp.snp_id}{suffix}_ALT",
            sequence=alt_seq,
            snp_id=snp.snp_id,
            chrom=snp.chrom,
            start=start,
            end=end - 1,
            snp_offset=snp_in_probe,
            allele=snp.alt,
            is_ref=False,
        )
        probes.append(alt_probe)
    
    return probes


def run_tiling(
    input_path: Path,
    reference_path: Optional[Path],
    output_prefix: Path,
    probe_length: int = 81,
    positions: str = "center,left,right",
    offset: int = 20,
    suffix_style: str = "position",
    step: Optional[int] = None,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Generate multi-position tiled probes.
    
    Main entry point for the tiling command.
    
    For SNP input: Creates multiple probes per SNP at different positions.
    For FASTA input: Uses sliding window tiling.
    
    Args:
        input_path: Input file (SNP TSV or FASTA)
        reference_path: Reference genome (required for SNP input)
        output_prefix: Output file prefix
        probe_length: Probe length
        positions: Position specification string
        offset: Default offset for named positions
        suffix_style: Suffix style (position or index)
        step: Step size for FASTA tiling
        verbose: Enable verbose logging
        
    Returns:
        Result containing tiling statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Starting tiled probe generation from {input_path}")
    
    # Parse positions
    pos_result = parse_positions(positions, offset)
    if pos_result.is_err():
        return Err(pos_result.unwrap_err())
    
    tile_positions = pos_result.unwrap()
    logger.info(f"Tiling positions: {[p.name for p in tile_positions]}")
    
    # Determine input type
    suffix = input_path.suffix.lower()
    is_fasta = suffix in ['.fa', '.fasta', '.fna']
    
    if is_fasta:
        # FASTA tiling - use funcgen from_fasta
        from eprobe.funcgen.from_fasta import run_from_fasta
        
        step_size = step or 30
        return run_from_fasta(
            fasta_path=input_path,
            output_prefix=output_prefix,
            probe_length=probe_length,
            step_size=step_size,
            verbose=verbose,
        )
    
    # SNP-based tiling
    if reference_path is None:
        return Err("Reference genome required for SNP input")
    
    # Load SNPs
    snp_df = SNPDataFrame.from_tsv(input_path)
    if snp_df.is_err():
        return Err(f"Failed to load SNPs: {snp_df.unwrap_err()}")
    
    snps = snp_df.unwrap().to_snps()
    logger.info(f"Loaded {len(snps)} SNPs")
    
    # Load reference
    ref_result = read_fasta(reference_path)
    if ref_result.is_err():
        return Err(f"Failed to read reference: {ref_result.unwrap_err()}")
    
    reference_sequences = ref_result.unwrap()
    
    # Generate tiled probes
    all_probes: List[Probe] = []
    
    for snp in snps:
        ref_seq = reference_sequences.get(snp.chrom)
        if ref_seq is None:
            logger.warning(f"Chromosome {snp.chrom} not in reference")
            continue
        
        probes = generate_tiled_probes_snp(
            snp=snp,
            reference_seq=ref_seq,
            chrom_length=len(ref_seq),
            probe_length=probe_length,
            positions=tile_positions,
            suffix_style=suffix_style,
        )
        all_probes.extend(probes)
    
    logger.info(f"Generated {len(all_probes)} tiled probes")
    
    # Save outputs
    fasta_output = Path(str(output_prefix) + ".tiled.fa")
    fasta_output.parent.mkdir(parents=True, exist_ok=True)
    
    probe_sequences = {p.probe_id: p.sequence for p in all_probes}
    write_result = write_fasta(probe_sequences, fasta_output)
    if write_result.is_err():
        return Err(f"Failed to save FASTA: {write_result.unwrap_err()}")
    
    # Calculate tiling factor
    tiling_factor = len(all_probes) / (2 * len(snps)) if snps else 0  # /2 for ref+alt
    
    stats = {
        "input_count": len(snps),
        "probe_count": len(all_probes),
        "tiling_factor": round(tiling_factor, 1),
        "positions": [p.name for p in tile_positions],
        "probe_length": probe_length,
        "fasta_file": str(fasta_output),
    }
    
    return Ok(stats)
