"""
Probe generation from BED regions.

Extracts sequences from reference genome at BED-specified regions,
then generates tiled probes. Supports haplotype inference via VCF.

This module corresponds to part of the original Seq_generator.py functionality.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta, extract_sequence
from eprobe.core.models import Probe
from eprobe.funcgen.from_fasta import tile_sequence, TilingConfig

logger = logging.getLogger(__name__)


@dataclass
class BedRegion:
    """Representation of a BED region."""
    chrom: str
    start: int  # 0-based
    end: int    # 0-based, exclusive
    name: Optional[str] = None
    score: Optional[float] = None
    strand: Optional[str] = None
    
    @property
    def length(self) -> int:
        return self.end - self.start


def parse_bed_file(bed_path: Path) -> Result[List[BedRegion], str]:
    """
    Parse BED file into list of BedRegion objects.
    
    Args:
        bed_path: Path to BED file
        
    Returns:
        Result containing list of BedRegion objects
    """
    regions = []
    
    try:
        with open(bed_path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # Skip empty lines and headers
                if not line or line.startswith('#') or line.startswith('track'):
                    continue
                
                parts = line.split('\t')
                
                if len(parts) < 3:
                    logger.warning(f"Line {line_num}: insufficient columns, skipping")
                    continue
                
                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    
                    name = parts[3] if len(parts) > 3 else f"region_{line_num}"
                    score = float(parts[4]) if len(parts) > 4 and parts[4] != '.' else None
                    strand = parts[5] if len(parts) > 5 else None
                    
                    region = BedRegion(
                        chrom=chrom,
                        start=start,
                        end=end,
                        name=name,
                        score=score,
                        strand=strand,
                    )
                    regions.append(region)
                    
                except (ValueError, IndexError) as e:
                    logger.warning(f"Line {line_num}: parse error ({e}), skipping")
                    continue
        
        return Ok(regions)
        
    except Exception as e:
        return Err(f"Failed to parse BED file: {e}")


def extract_region_sequence(
    region: BedRegion,
    reference_sequences: Dict[str, str],
) -> Result[str, str]:
    """
    Extract sequence for a BED region from reference.
    
    Args:
        region: BED region
        reference_sequences: Reference sequences by chromosome
        
    Returns:
        Result containing extracted sequence
    """
    chrom_seq = reference_sequences.get(region.chrom)
    
    if chrom_seq is None:
        return Err(f"Chromosome {region.chrom} not in reference")
    
    if region.start < 0 or region.end > len(chrom_seq):
        return Err(f"Region {region.name} ({region.chrom}:{region.start}-{region.end}) "
                   f"outside chromosome bounds (length: {len(chrom_seq)})")
    
    sequence = chrom_seq[region.start:region.end].upper()
    
    # Handle strand
    if region.strand == '-':
        from eprobe.core.fasta import reverse_complement
        sequence = reverse_complement(sequence)
    
    return Ok(sequence)


def run_from_bed(
    bed_path: Path,
    reference_path: Path,
    output_prefix: Path,
    probe_length: int = 81,
    step_size: int = 30,
    vcf_path: Optional[Path] = None,
    phase: bool = False,
    min_freq: float = 0.05,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Generate probes from BED regions.
    
    Main entry point for the from_bed command. Extracts sequences
    from reference at BED regions, then tiles probes.
    
    Args:
        bed_path: Input BED file
        reference_path: Reference genome FASTA
        output_prefix: Output file prefix
        probe_length: Probe length
        step_size: Sliding window step
        vcf_path: VCF for haplotype inference (optional)
        phase: Enable haplotype phasing
        min_freq: Minimum haplotype frequency
        threads: Number of threads
        verbose: Enable verbose logging
        
    Returns:
        Result containing generation statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Starting probe generation from {bed_path}")
    logger.info(f"Reference: {reference_path}")
    logger.info(f"Probe length: {probe_length}, Step: {step_size}")
    
    # Parse BED file
    bed_result = parse_bed_file(bed_path)
    if bed_result.is_err():
        return Err(bed_result.unwrap_err())
    
    regions = bed_result.unwrap()
    logger.info(f"Loaded {len(regions)} regions from BED")
    
    if not regions:
        return Err("No valid regions found in BED file")
    
    # Load reference genome
    logger.info("Loading reference genome...")
    ref_result = read_fasta(reference_path)
    if ref_result.is_err():
        return Err(f"Failed to read reference: {ref_result.unwrap_err()}")
    
    reference_sequences = ref_result.unwrap()
    logger.info(f"Loaded {len(reference_sequences)} chromosomes")
    
    # Configure tiling
    tiling_config = TilingConfig(
        probe_length=probe_length,
        step_size=step_size,
    )
    
    # Extract sequences and generate probes
    all_probes: List[Probe] = []
    extraction_errors = 0
    
    for region in regions:
        # Extract sequence
        seq_result = extract_region_sequence(region, reference_sequences)
        
        if seq_result.is_err():
            logger.warning(f"Failed to extract {region.name}: {seq_result.unwrap_err()}")
            extraction_errors += 1
            continue
        
        sequence = seq_result.unwrap()
        
        # Tile probes
        probes = tile_sequence(region.name, sequence, tiling_config)
        
        # Update coordinates to genomic (for reference strand regions)
        for probe in probes:
            if region.strand != '-':
                probe.start = region.start + probe.start
                probe.end = region.start + probe.end
            probe.chrom = region.chrom
        
        all_probes.extend(probes)
    
    logger.info(f"Generated {len(all_probes)} probes from {len(regions) - extraction_errors} regions")
    
    if extraction_errors > 0:
        logger.warning(f"Failed to extract {extraction_errors} regions")
    
    # Handle haplotype phasing (placeholder for future implementation)
    if phase and vcf_path:
        logger.warning("Haplotype phasing not yet implemented for BED mode")
    
    # Save FASTA output
    fasta_output = Path(str(output_prefix) + ".probes.fa")
    fasta_output.parent.mkdir(parents=True, exist_ok=True)
    
    probe_sequences = {p.probe_id: p.sequence for p in all_probes}
    write_result = write_fasta(probe_sequences, fasta_output)
    if write_result.is_err():
        return Err(f"Failed to save FASTA: {write_result.unwrap_err()}")
    
    logger.info(f"Saved probes to {fasta_output}")
    
    # Save TSV metadata
    tsv_output = Path(str(output_prefix) + ".probes.tsv")
    
    with open(tsv_output, 'w') as f:
        f.write("probe_id\tsequence\tchrom\tstart\tend\tregion\tlength\n")
        for p in all_probes:
            f.write(f"{p.probe_id}\t{p.sequence}\t{p.chrom}\t{p.start}\t{p.end}\t{p.snp_id}\t{len(p.sequence)}\n")
    
    stats = {
        "region_count": len(regions) - extraction_errors,
        "probe_count": len(all_probes),
        "extraction_errors": extraction_errors,
        "probe_length": probe_length,
        "step_size": step_size,
        "fasta_file": str(fasta_output),
        "tsv_file": str(tsv_output),
    }
    
    return Ok(stats)
