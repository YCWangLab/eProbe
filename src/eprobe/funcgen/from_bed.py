"""
Probe generation from BED regions.

Entry points:
  - Simple mode: extract sequences from reference at BED regions, tile probes
  - Haplotype mode: use phased VCF to construct per-region haplotype sequences,
    then generate haplotype-aware probes

The shared function ``process_regions()`` is also used by from_gff.

Legacy equivalent: Get_probe_from_bed.py + Get_allele_from_vcf.py + Get_probe_from_allele.py (bed mode)
"""

import logging
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pysam

from eprobe.core.fasta import write_fasta
from eprobe.core.models import Probe
from eprobe.core.result import Err, Ok, Result
from eprobe.funcgen.haplotype import (
    HaplotypeSet,
    extract_haplotypes,
    find_phase_common,
    generate_haplotype_probes,
    phase_vcf_region,
    tile_sequence,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# BED region data model
# ---------------------------------------------------------------------------

@dataclass
class BedRegion:
    """Representation of a BED region (0-based, half-open)."""
    chrom: str
    start: int    # 0-based
    end: int      # 0-based exclusive
    name: Optional[str] = None
    score: Optional[float] = None
    strand: Optional[str] = None

    @property
    def length(self) -> int:
        return self.end - self.start


# ---------------------------------------------------------------------------
# BED parsing
# ---------------------------------------------------------------------------

def parse_bed_file(bed_path: Path) -> Result[List[BedRegion], str]:
    """
    Parse BED file into list of BedRegion objects.

    Supports 3-6 column BED format. Skips comment lines and headers.

    Args:
        bed_path: Path to BED file

    Returns:
        Ok(list of BedRegion), Err(message) on failure
    """
    regions: List[BedRegion] = []

    try:
        with open(bed_path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue

                parts = line.split("\t")
                if len(parts) < 3:
                    logger.warning(f"BED line {line_num}: insufficient columns, skipping")
                    continue

                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
                    score = float(parts[4]) if len(parts) > 4 and parts[4] != "." else None
                    strand = parts[5] if len(parts) > 5 else None

                    regions.append(BedRegion(
                        chrom=chrom, start=start, end=end,
                        name=name, score=score, strand=strand,
                    ))
                except (ValueError, IndexError) as e:
                    logger.warning(f"BED line {line_num}: parse error ({e}), skipping")

        return Ok(regions)

    except Exception as e:
        return Err(f"Failed to parse BED file: {e}")


# ---------------------------------------------------------------------------
# Per-region worker (used by parallel haplotype mode)
# ---------------------------------------------------------------------------

def _process_one_region_haplo(
    region: BedRegion,
    reference_path: Path,
    vcf_path: Path,
    temp_dir: Path,
    phase: bool,
    min_freq: float,
    variant_only: bool,
    probe_length: int,
    step_size: int,
    aligner: str,
) -> Tuple[str, List[Probe], bool, int]:
    """Process a single region in haplotype mode.

    Returns:
        (status, probes, has_multi_haplo, n_haplo)
        status: "ok", "err", or "ref_only"
    """
    rid = region.name or f"{region.chrom}:{region.start}-{region.end}"

    # Step 1: Phase if requested
    if phase:
        phase_result = phase_vcf_region(
            vcf_path, region.chrom, region.start, region.end,
            temp_dir, rid,
        )
        if phase_result.is_err():
            err_msg = phase_result.unwrap_err()
            # No variants in region → fall back to reference sequence
            if err_msg.startswith("no_variants:"):
                logger.debug(f"  {rid}: no variants in VCF, using reference")
            else:
                logger.warning(f"  {rid}: {err_msg}")
            if not variant_only:
                ref_fasta = pysam.FastaFile(str(reference_path))
                ref_seq = ref_fasta.fetch(region.chrom, region.start, region.end)
                ref_fasta.close()
                probes = tile_sequence(
                    rid, ref_seq, probe_length, step_size,
                    source_chrom=region.chrom,
                    source_offset=region.start,
                )
                return ("ref_only", probes, False, 0)
            return ("ref_only", [], False, 0)
        vcf_to_use = phase_result.unwrap()
    else:
        vcf_to_use = vcf_path

    # Step 2: Extract haplotypes
    haplo_result = extract_haplotypes(
        phased_vcf_path=vcf_to_use,
        reference_fasta_path=reference_path,
        chrom=region.chrom,
        start=region.start,
        end=region.end,
        min_freq=min_freq,
        region_id=rid,
    )

    if haplo_result.is_err():
        return ("err", [], False, 0)

    haplo_set = haplo_result.unwrap()
    has_multi = len(haplo_set.alleles) > 1

    # Step 3: Generate haplotype-aware probes
    probe_result = generate_haplotype_probes(
        region_id=rid,
        alleles=haplo_set.alleles,
        probe_length=probe_length,
        step_size=step_size,
        variant_only=variant_only,
        source_chrom=region.chrom,
        source_offset=region.start,
        aligner=aligner,
    )

    if probe_result.is_err():
        return ("err", [], False, 0)

    probes, _stats = probe_result.unwrap()
    return ("ok", probes, has_multi, len(haplo_set.alleles))


# ---------------------------------------------------------------------------
# Core region processing (shared with from_gff)
# ---------------------------------------------------------------------------

def process_regions(
    regions: List[BedRegion],
    reference_path: Path,
    output_prefix: Path,
    probe_length: int,
    step_size: int,
    vcf_path: Optional[Path] = None,
    phase: bool = False,
    min_freq: float = 0.05,
    variant_only: bool = False,
    aligner: str = "muscle",
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Generate probes for a list of genomic regions.

    This is the shared engine used by both from_bed and from_gff.

    Simple mode (no VCF):
      Extracts reference sequences at each BED region, tiles probes.

    Haplotype mode (with VCF):
      For each region:
        1. Phase VCF if requested (shapeit5 per region)
        2. Extract per-sample haplotype sequences from phased genotypes
        3. Frequency-filter unique haplotypes
        4. Generate haplotype-aware probes (1 probe at non-variant,
           N probes at variant positions)

    Args:
        regions: List of BedRegion objects
        reference_path: Reference genome FASTA (must have .fai index)
        output_prefix: Output file prefix
        probe_length: Probe length (bp)
        step_size: Sliding window step (bp)
        vcf_path: VCF for haplotype inference (optional)
        phase: Run shapeit5 phasing before haplotype extraction
        min_freq: Minimum haplotype frequency to retain
        variant_only: Only generate probes at variant positions
        aligner: MSA aligner for indel-containing haplotypes
        threads: Number of threads (for parallel region processing)
        verbose: Enable verbose logging

    Returns:
        Ok(stats dict) on success, Err(message) on failure
    """
    if not regions:
        return Err("No regions to process")

    # Early check: if phasing requested, verify external tools exist
    if phase:
        if find_phase_common() is None:
            return Err(
                "phase_common / SHAPEIT5_phase_common not found in PATH. "
                "Install shapeit5 ('conda install -c bioconda shapeit5') "
                "or use --no_phase if VCF is already phased."
            )
        for tool in ("bcftools", "bgzip", "tabix"):
            if not shutil.which(tool):
                return Err(f"'{tool}' not found in PATH. Required for --phase.")

    all_probes: List[Probe] = []
    n_regions_ok = 0
    n_regions_err = 0
    n_haplo_regions = 0  # regions with >1 haplotype

    if vcf_path is not None:
        # ======================= HAPLOTYPE MODE =======================
        logger.info(f"Haplotype mode: VCF={vcf_path}, phase={phase}, min_freq={min_freq}")

        temp_dir = Path(str(output_prefix) + "_eprobe_temp")
        temp_dir.mkdir(parents=True, exist_ok=True)

        try:
            n_workers = min(threads, len(regions))
            if n_workers > 1:
                logger.info(f"Processing {len(regions)} regions with {n_workers} workers")
                futures = {}
                with ThreadPoolExecutor(max_workers=n_workers) as executor:
                    for region in regions:
                        rid = region.name or f"{region.chrom}:{region.start}-{region.end}"
                        fut = executor.submit(
                            _process_one_region_haplo,
                            region, reference_path, vcf_path, temp_dir,
                            phase, min_freq, variant_only,
                            probe_length, step_size, aligner,
                        )
                        futures[fut] = rid

                    for fut in as_completed(futures):
                        rid = futures[fut]
                        try:
                            status, probes, has_multi, n_haplo = fut.result()
                        except Exception as exc:
                            logger.warning(f"  {rid}: worker exception: {exc}")
                            n_regions_err += 1
                            continue
                        if status == "err":
                            n_regions_err += 1
                        else:
                            all_probes.extend(probes)
                            n_regions_ok += 1
                            if has_multi:
                                n_haplo_regions += 1
            else:
                # Single-threaded: process regions sequentially
                for region in regions:
                    rid = region.name or f"{region.chrom}:{region.start}-{region.end}"
                    status, probes, has_multi, n_haplo = _process_one_region_haplo(
                        region, reference_path, vcf_path, temp_dir,
                        phase, min_freq, variant_only,
                        probe_length, step_size, aligner,
                    )
                    if status == "err":
                        n_regions_err += 1
                    else:
                        all_probes.extend(probes)
                        n_regions_ok += 1
                        if has_multi:
                            n_haplo_regions += 1

        finally:
            # Cleanup temp directory
            if temp_dir.exists():
                shutil.rmtree(temp_dir, ignore_errors=True)
                logger.debug(f"Cleaned up temp dir: {temp_dir}")

    else:
        # ======================= SIMPLE MODE ==========================
        logger.info("Simple mode: extracting reference sequences")

        ref_fasta = pysam.FastaFile(str(reference_path))

        for region in regions:
            rid = region.name or f"{region.chrom}:{region.start}-{region.end}"

            try:
                seq = ref_fasta.fetch(region.chrom, region.start, region.end)
            except Exception as e:
                logger.warning(f"  {rid}: failed to extract sequence: {e}")
                n_regions_err += 1
                continue

            if not seq:
                logger.warning(f"  {rid}: empty sequence")
                n_regions_err += 1
                continue

            # Handle reverse strand
            if region.strand == "-":
                from eprobe.core.fasta import reverse_complement
                seq = reverse_complement(seq)

            probes = tile_sequence(
                rid, seq, probe_length, step_size,
                source_chrom=region.chrom,
                source_offset=region.start,
            )
            all_probes.extend(probes)
            n_regions_ok += 1

        ref_fasta.close()

    logger.info(
        f"Processed {n_regions_ok} regions OK, {n_regions_err} errors, "
        f"{n_haplo_regions} with multiple haplotypes"
    )
    logger.info(f"Generated {len(all_probes)} probes total")

    if not all_probes:
        return Err("No probes generated")

    # --- Save FASTA output ---
    fasta_output = Path(str(output_prefix) + ".probes.fasta")
    fasta_output.parent.mkdir(parents=True, exist_ok=True)

    probe_seqs = {p.id: p.sequence for p in all_probes}
    wr = write_fasta(probe_seqs, fasta_output)
    if wr.is_err():
        return Err(f"Failed to write FASTA: {wr.unwrap_err()}")

    logger.info(f"Saved probes to {fasta_output}")

    stats = {
        "region_count": n_regions_ok,
        "region_errors": n_regions_err,
        "haplo_regions": n_haplo_regions,
        "probe_count": len(all_probes),
        "probe_length": probe_length,
        "step_size": step_size,
        "haplotyping": vcf_path is not None,
        "variant_only": variant_only,
        "fasta_file": str(fasta_output),
    }

    return Ok(stats)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_from_bed(
    bed_path: Path,
    reference_path: Path,
    output_prefix: Path,
    probe_length: int = 81,
    step_size: int = 30,
    vcf_path: Optional[Path] = None,
    phase: bool = False,
    min_freq: float = 0.05,
    variant_only: bool = False,
    aligner: str = "muscle",
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Generate probes from BED regions.

    Args:
        bed_path: Input BED file
        reference_path: Reference genome FASTA (must have .fai index)
        output_prefix: Output file prefix
        probe_length: Probe length (bp)
        step_size: Sliding window step (bp)
        vcf_path: VCF for haplotype inference (optional)
        phase: Phase VCF with shapeit5 before haplotype extraction
        min_freq: Minimum haplotype frequency to retain
        variant_only: Only generate probes at variant positions
        aligner: MSA aligner ('muscle', 'clustalo', 'mafft')
        threads: Number of threads
        verbose: Enable verbose logging

    Returns:
        Ok(stats dict) on success, Err(message) on failure

    Output files:
        {output_prefix}.probes.fasta - Probe sequences
        {output_prefix}.probes.tsv  - Probe metadata
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

    # Delegate to shared region processor
    return process_regions(
        regions=regions,
        reference_path=reference_path,
        output_prefix=output_prefix,
        probe_length=probe_length,
        step_size=step_size,
        vcf_path=vcf_path,
        phase=phase,
        min_freq=min_freq,
        variant_only=variant_only,
        aligner=aligner,
        threads=threads,
        verbose=verbose,
    )
