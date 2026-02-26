"""
Probe sequence generation module.

Generates final probe sequences from selected SNPs:
  - Extracts sequences with configurable length and offset
  - Supports SNP replacement (third-base) to reduce allele capture bias
  - Supports mutation type filtering (transition/transversion)
  - Handles edge cases near chromosome boundaries

This module corresponds to the original SNP_generator.py functionality.

FASTA header format (legacy-compatible):
  >{chrom}:{start}-{end}_{pos}_{type}_{ref}_{alt}
"""

import logging
import random
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from collections import OrderedDict
from dataclasses import dataclass

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNP, SNPDataFrame, Probe, ProbeSet
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class ProbeConfig:
    """Configuration for probe generation."""
    length: int = 81
    shift: int = 0        # SNP offset from center (0 = centered)
    replace_snp: bool = False  # Replace SNP with third base (neither REF nor ALT)
    mutation_type: Optional[str] = None  # Filter: "ts", "tv", or None (both)
    tiling: bool = False   # Generate 3 probes per SNP (center + left + right)
    tiling_shift: int = 15 # Offset for tiling left/right probes


# =============================================================================
# Coordinate Calculation
# =============================================================================

def calculate_probe_coordinates(
    snp_pos: int,
    probe_length: int,
    shift: int = 0,
) -> Tuple[int, int, int, int]:
    """
    Calculate probe start, end, and SNP offset within probe.
    
    For odd-length probes (recommended), left_flank == right_flank.
    For even-length probes, right_flank = left_flank + 1.
    
    Shift follows the legacy convention:
      - Positive shift: probe slides right (SNP moves toward 5' end)
      - Negative shift: probe slides left (SNP moves toward 3' end)
    
    Args:
        snp_pos: 1-based SNP position
        probe_length: Desired probe length
        shift: Offset from center
        
    Returns:
        Tuple of (start_pos, end_pos, snp_offset, actual_length)
        - start_pos: 0-based start position (for Python slicing)
        - end_pos: 0-based end position (exclusive, for Python slicing)
        - snp_offset: 0-based offset of SNP within probe
        - actual_length: actual probe length
    """
    if probe_length % 2 == 1:
        left_flank = right_flank = probe_length // 2
    else:
        left_flank = (probe_length // 2) - 1
        right_flank = left_flank + 1
    
    # 1-based to 0-based; apply shift (positive = slide right)
    start = snp_pos - left_flank - 1 + shift
    end = snp_pos + right_flank - 1 + shift
    
    snp_offset = snp_pos - 1 - start  # 0-based offset within probe
    
    return start, end + 1, snp_offset, end - start + 1


def validate_shift(shift: int, probe_length: int) -> Result[None, str]:
    """
    Validate that shift does not move SNP outside the probe.
    
    Legacy rule: abs(shift) must be < min(left_flank, right_flank).
    """
    if probe_length % 2 == 1:
        max_shift = probe_length // 2
    else:
        max_shift = (probe_length // 2) - 1
    
    if abs(shift) >= max_shift:
        return Err(
            f"Shift {shift} is too large for probe length {probe_length}. "
            f"Maximum allowed: {max_shift - 1} (to keep SNP inside probe)."
        )
    return Ok(None)


# =============================================================================
# Probe Generation
# =============================================================================

def make_probe_id(snp: SNP, start_1based: int, end_1based: int) -> str:
    """
    Generate legacy-compatible FASTA header / probe ID.
    
    Format: {chrom}:{start}-{end}_{pos}_{type}_{ref}_{alt}
    """
    return (
        f"{snp.chrom}:{start_1based}-{end_1based}"
        f"_{snp.pos}_{snp.mutation_type}_{snp.ref}_{snp.alt}"
    )


def replace_with_third_base(ref: str, alt: str, seed: Optional[int] = None) -> str:
    """
    Pick a random base that is neither REF nor ALT.
    
    Used to reduce allele capture bias in hybridization.
    """
    bases = {"A", "T", "C", "G"}
    candidates = list(bases - {ref.upper(), alt.upper()})
    if seed is not None:
        rng = random.Random(seed)
        return rng.choice(candidates)
    return random.choice(candidates)


def generate_probe_for_snp(
    snp: SNP,
    reference_seq: str,
    chrom_length: int,
    config: ProbeConfig,
) -> Result[Probe, str]:
    """
    Generate one probe for a single SNP.
    
    Args:
        snp: SNP object
        reference_seq: Full chromosome sequence (0-based indexing)
        chrom_length: Length of chromosome
        config: Probe generation configuration
        
    Returns:
        Result containing a Probe object
    """
    start_0, end_0, snp_offset, length = calculate_probe_coordinates(
        snp.pos, config.length, config.shift
    )
    
    # Boundary check
    if start_0 < 0:
        return Err(
            f"Probe extends before chromosome start for "
            f"{snp.chrom}:{snp.pos} (start={start_0})"
        )
    if end_0 > chrom_length:
        return Err(
            f"Probe extends beyond chromosome end for "
            f"{snp.chrom}:{snp.pos} (end={end_0}, chrom_len={chrom_length})"
        )
    
    # Extract sequence
    probe_seq = reference_seq[start_0:end_0].upper()
    
    if len(probe_seq) != config.length:
        return Err(
            f"Extracted sequence length {len(probe_seq)} != "
            f"expected {config.length} for {snp.chrom}:{snp.pos}"
        )
    
    # Validate REF base matches reference genome
    ref_at_snp = probe_seq[snp_offset]
    if ref_at_snp != snp.ref.upper():
        return Err(
            f"Reference mismatch at {snp.chrom}:{snp.pos}: "
            f"expected '{snp.ref}', found '{ref_at_snp}' in genome"
        )
    
    # Apply SNP replacement if requested
    if config.replace_snp:
        new_base = replace_with_third_base(snp.ref, snp.alt)
        seq_list = list(probe_seq)
        seq_list[snp_offset] = new_base
        probe_seq = "".join(seq_list)
    
    # Build probe ID (legacy format: chr:start-end_pos_type_ref_alt)
    start_1based = start_0 + 1
    end_1based = end_0  # end_0 is exclusive in 0-based, so equals 1-based end
    probe_id = make_probe_id(snp, start_1based, end_1based)
    
    probe = Probe(
        id=probe_id,
        sequence=probe_seq,
        source_chrom=snp.chrom,
        source_start=start_1based,
        source_end=end_1based,
        snp_pos=snp.pos,
        snp_ref=snp.ref,
        snp_alt=snp.alt,
        mutation_type=snp.mutation_type,
    )
    
    return Ok(probe)


def generate_probes_batch(
    snps: List[SNP],
    reference_sequences: Dict[str, str],
    config: ProbeConfig,
) -> Result[ProbeSet, str]:
    """
    Generate probes for a batch of SNPs.
    
    In tiling mode, generates 3 probes per SNP:
      - Center probe (shift=0)
      - Left probe (shift=+tiling_shift, SNP toward 5' end)
      - Right probe (shift=-tiling_shift, SNP toward 3' end)
    
    Args:
        snps: List of SNPs
        reference_sequences: Reference sequences by chromosome
        config: Probe generation configuration
        
    Returns:
        Result containing ProbeSet
    """
    probe_set = ProbeSet(name="popgen_probes")
    errors: List[str] = []
    
    if config.tiling:
        # Tiling mode: 3 probes per SNP
        shifts = [
            (0, "C"),                        # center
            (config.tiling_shift, "L"),       # left-shifted (SNP toward 5')
            (-config.tiling_shift, "R"),      # right-shifted (SNP toward 3')
        ]
        for snp in snps:
            ref_seq = reference_sequences.get(snp.chrom)
            if ref_seq is None:
                errors.append(f"Chromosome '{snp.chrom}' not in reference")
                continue
            
            for shift_val, tag in shifts:
                tile_config = ProbeConfig(
                    length=config.length,
                    shift=shift_val,
                    replace_snp=config.replace_snp,
                    mutation_type=config.mutation_type,
                )
                result = generate_probe_for_snp(
                    snp, ref_seq, len(ref_seq), tile_config
                )
                if result.is_err():
                    errors.append(f"{tag}:{result.unwrap_err()}")
                    logger.debug(f"Skipped tiling {tag}: {result.unwrap_err()}")
                    continue
                
                probe = result.unwrap()
                # Append tiling tag to probe ID for uniqueness
                tagged_probe = Probe(
                    id=f"{probe.id}_{tag}",
                    sequence=probe.sequence,
                    source_chrom=probe.source_chrom,
                    source_start=probe.source_start,
                    source_end=probe.source_end,
                    snp_pos=probe.snp_pos,
                    snp_ref=probe.snp_ref,
                    snp_alt=probe.snp_alt,
                    mutation_type=probe.mutation_type,
                )
                probe_set.add(tagged_probe)
    else:
        # Standard mode: 1 probe per SNP
        for snp in snps:
            ref_seq = reference_sequences.get(snp.chrom)
            if ref_seq is None:
                errors.append(f"Chromosome '{snp.chrom}' not in reference")
                continue
            
            result = generate_probe_for_snp(snp, ref_seq, len(ref_seq), config)
            
            if result.is_err():
                errors.append(result.unwrap_err())
                logger.debug(f"Skipped: {result.unwrap_err()}")
                continue
            
            probe_set.add(result.unwrap())
    
    if errors:
        logger.warning(f"Failed to generate {len(errors)} probes")
        for err in errors[:5]:
            logger.warning(f"  {err}")
        if len(errors) > 5:
            logger.warning(f"  ... and {len(errors) - 5} more")
    
    return Ok(probe_set)


# =============================================================================
# Main Entry Point
# =============================================================================

def run_build(
    input_path: Path,
    reference_path: Path,
    output_prefix: Path,
    probe_length: int = 81,
    snp_position: str = "center",
    offset: int = 0,
    replace_snp: bool = False,
    mutation_type: Optional[str] = None,
    tiling_offset: int = 15,
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
        probe_length: Probe length (odd number recommended)
        snp_position: "center", "left", "right", or "tiling"
        offset: Additional offset from center position
        replace_snp: Replace SNP with third base (neither REF nor ALT)
        mutation_type: Filter by "ts" or "tv", None = keep both
        tiling_offset: Shift for tiling left/right probes (default: 15)
        verbose: Enable verbose logging
        
    Returns:
        Result containing build statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    # Tiling mode: 3 probes per SNP (center, left, right)
    use_tiling = snp_position == "tiling"
    
    if use_tiling:
        shift = 0  # center probe is at 0
        # Validate tiling_offset
        shift_check = validate_shift(tiling_offset, probe_length)
        if shift_check.is_err():
            return Err(
                f"Tiling offset {tiling_offset} is too large. "
                f"{shift_check.unwrap_err()}"
            )
    else:
        # Convert snp_position to shift value
        # Legacy convention: positive shift slides probe right
        if snp_position == "center":
            shift = offset
        elif snp_position == "left":
            # SNP toward left (5') end → probe slides right → positive shift
            shift = abs(offset) if offset != 0 else probe_length // 4
        elif snp_position == "right":
            # SNP toward right (3') end → probe slides left → negative shift
            shift = -(abs(offset) if offset != 0 else probe_length // 4)
        else:
            return Err(f"Invalid snp_position: {snp_position}")
        
        # Validate shift range
        shift_check = validate_shift(shift, probe_length)
        if shift_check.is_err():
            return Err(shift_check.unwrap_err())
    
    logger.info(f"Starting probe generation from {input_path}")
    logger.info(f"Probe length: {probe_length}, SNP position: {snp_position}, "
                f"Shift: {shift}, Replace SNP: {replace_snp}")
    
    # Load input SNPs
    snp_df_result = SNPDataFrame.from_tsv(input_path)
    if snp_df_result.is_err():
        return Err(f"Failed to load input: {snp_df_result.unwrap_err()}")
    
    snp_df = snp_df_result.unwrap()
    
    # Filter by mutation type if specified
    if mutation_type in ("ts", "tv"):
        before = len(snp_df)
        snp_df = snp_df.filter_by_mutation_type(mutation_type)
        logger.info(f"Filtered by mutation type '{mutation_type}': "
                    f"{before} → {len(snp_df)} SNPs")
    
    snps = snp_df.to_snps()
    if not snps:
        return Err("No SNPs remaining after filtering")
    
    logger.info(f"Loaded {len(snps)} SNPs")
    
    # Load reference genome
    logger.info(f"Loading reference genome from {reference_path}...")
    ref_result = read_fasta(reference_path)
    if ref_result.is_err():
        return Err(f"Failed to read reference: {ref_result.unwrap_err()}")
    
    reference_sequences = ref_result.unwrap()
    logger.info(f"Loaded {len(reference_sequences)} chromosomes")
    
    # Configure and generate
    config = ProbeConfig(
        length=probe_length,
        shift=shift,
        replace_snp=replace_snp,
        mutation_type=mutation_type,
        tiling=use_tiling,
        tiling_shift=tiling_offset,
    )
    
    logger.info("Generating probe sequences...")
    probe_result = generate_probes_batch(snps, reference_sequences, config)
    if probe_result.is_err():
        return Err(probe_result.unwrap_err())
    
    probe_set = probe_result.unwrap()
    
    if len(probe_set) == 0:
        return Err("No probes were generated. Check reference genome and SNP coordinates.")
    
    # Save FASTA output
    fasta_path = Path(str(output_prefix) + ".probes.fa")
    write_result = probe_set.to_fasta(fasta_path)
    if write_result.is_err():
        return Err(f"Failed to save FASTA: {write_result.unwrap_err()}")
    
    logger.info(f"Saved {len(probe_set)} probes to {fasta_path}")
    
    # Save summary
    summary_path = Path(str(output_prefix) + ".build_summary.txt")
    ts_count = sum(1 for p in probe_set if p.mutation_type == "ts")
    tv_count = sum(1 for p in probe_set if p.mutation_type == "tv")
    
    with open(summary_path, "w") as f:
        f.write(f"Probe Generation Summary\n")
        f.write(f"========================\n")
        f.write(f"Input SNPs: {len(snps)}\n")
        f.write(f"Probes generated: {len(probe_set)}\n")
        f.write(f"Probe length: {probe_length}\n")
        f.write(f"SNP position: {snp_position} (shift={shift})\n")
        if use_tiling:
            f.write(f"Tiling mode: 3 probes per SNP (center, L±{tiling_offset}, R±{tiling_offset})\n")
        f.write(f"Replace SNP (3rd base): {replace_snp}\n")
        if mutation_type:
            f.write(f"Mutation type filter: {mutation_type}\n")
        f.write(f"Transitions (ts): {ts_count}\n")
        f.write(f"Transversions (tv): {tv_count}\n")
        f.write(f"\nOutput: {fasta_path}\n")
    
    stats = {
        "snp_count": len(snps),
        "probe_count": len(probe_set),
        "ts_count": ts_count,
        "tv_count": tv_count,
        "probe_length": probe_length,
        "shift": shift,
        "replace_snp": replace_snp,
        "fasta_file": str(fasta_path),
        "summary_file": str(summary_path),
    }
    
    logger.info(f"Probe generation complete: {len(probe_set)} probes "
                f"from {len(snps)} SNPs (ts={ts_count}, tv={tv_count})")
    
    return Ok(stats)
