"""
Haplotype processing for funcgen probe design.

Core module providing:
- VCF phasing with shapeit5 (phase_common)
- Haplotype sequence extraction from phased VCF
- Multiple sequence alignment (muscle/clustalo/mafft)
- Haplotype-aware probe generation with variant/non-variant detection

The haplotype-aware probe generation strategy:
1. For single allele: simple sliding-window tiling
2. For multiple alleles of equal length (SNP-only):
   - Tile all alleles at same positions, deduplicate per window
   - 1 probe at non-variant windows, N probes at variant windows
3. For multiple alleles of different length (indels present):
   - MSA to establish coordinate mapping
   - Slide window across alignment, extract original sequences
   - Deduplicate per window

This avoids creating artificial haplotypes: all sequences come from
actual phased data or user-provided allele sequences.
"""

import logging
import os
import shutil
import subprocess
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pysam

from eprobe.core.fasta import write_fasta as _write_fasta_to_file
from eprobe.core.models import Probe
from eprobe.core.result import Err, Ok, Result

logger = logging.getLogger(__name__)


def find_phase_common() -> Optional[str]:
    """Find the shapeit5 phase_common executable.

    Returns the command name if found, None otherwise.
    Searches for both 'phase_common' and 'SHAPEIT5_phase_common'.
    """
    for name in ("phase_common", "SHAPEIT5_phase_common"):
        if shutil.which(name):
            return name
    return None


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class HaplotypeSet:
    """
    Collection of unique haplotype sequences for a genomic region.

    Attributes:
        region_id: Human-readable region identifier
        alleles: allele_label -> DNA sequence (e.g. {"H1": "ACGT...", "H2": "ACGT..."})
        frequencies: allele_label -> population frequency (0.0-1.0)
        n_samples: Number of diploid samples in the VCF
        n_variants: Number of variant sites in the region
    """
    region_id: str
    alleles: Dict[str, str]
    frequencies: Dict[str, float]
    n_samples: int = 0
    n_variants: int = 0


# ---------------------------------------------------------------------------
# Tiling: sliding-window probe generation
# ---------------------------------------------------------------------------

def tile_sequence(
    seq_id: str,
    sequence: str,
    probe_length: int,
    step_size: int,
    source_chrom: Optional[str] = None,
    source_offset: int = 0,
) -> List[Probe]:
    """
    Generate tiled probes from a single sequence using sliding window.

    Args:
        seq_id: Base identifier for probe naming
        sequence: Input DNA sequence
        probe_length: Desired probe length (bp)
        step_size: Sliding window step (bp)
        source_chrom: Chromosome name for genomic coordinates
        source_offset: Genomic start offset (0-based) for coordinate mapping

    Returns:
        List of Probe objects
    """
    probes = []
    seq_len = len(sequence)

    if seq_len < probe_length:
        logger.debug(f"{seq_id}: too short ({seq_len} < {probe_length}), skipping")
        return []

    positions = list(range(0, seq_len - probe_length + 1, step_size))

    # Ensure the last window is included (covers sequence end)
    last_start = seq_len - probe_length
    if positions and positions[-1] < last_start:
        positions.append(last_start)

    for idx, start in enumerate(positions):
        end = start + probe_length
        probe_seq = sequence[start:end].upper()

        if "N" in probe_seq:
            continue

        chrom = source_chrom or seq_id
        genomic_start = source_offset + start + 1  # 1-based
        genomic_end = source_offset + end

        probe = Probe(
            id=f"{seq_id}_P{idx + 1:04d}",
            sequence=probe_seq,
            source_chrom=chrom,
            source_start=genomic_start,
            source_end=genomic_end,
        )
        probes.append(probe)

    return probes


# ---------------------------------------------------------------------------
# VCF phasing with shapeit5
# ---------------------------------------------------------------------------

def phase_vcf_region(
    vcf_path: Path,
    chrom: str,
    start: int,
    end: int,
    output_dir: Path,
    region_name: str,
) -> Result[Path, str]:
    """
    Phase a VCF region using shapeit5 phase_common.

    Runs shapeit5 on a single region, converts output to VCF.gz + tabix index.

    Args:
        vcf_path: Input VCF file (bgzipped + tabix indexed)
        chrom: Chromosome
        start: Start position (BED 0-based)
        end: End position (BED 0-based exclusive)
        output_dir: Directory for intermediate files
        region_name: Unique name for this region (used for filenames)

    Returns:
        Ok(Path) to phased VCF.gz, Err(message) on failure
    """
    output_bcf = output_dir / f"{region_name}.bcf"
    output_vcf = output_dir / f"{region_name}.vcf"
    output_vcf_gz = output_dir / f"{region_name}.vcf.gz"

    # Reuse existing result
    if output_bcf.exists() and output_vcf_gz.exists():
        logger.debug(f"Region {region_name} already phased, reusing")
        return Ok(output_vcf_gz)

    try:
        region_str = f"{chrom}:{start}-{end}"

        # shapeit5 phase_common (supports both phase_common and SHAPEIT5_phase_common)
        phase_cmd = find_phase_common()
        if phase_cmd is None:
            return Err("phase_common / SHAPEIT5_phase_common not found in PATH")
        subprocess.run(
            f"{phase_cmd} --input {vcf_path} "
            f"--output-format bcf --output {output_bcf} "
            f"--region {region_str} --thread 1",
            shell=True, check=True, capture_output=True, text=True,
        )
        # BCF → VCF → VCF.gz + tabix
        subprocess.run(
            f"bcftools view {output_bcf} -Ov -o {output_vcf}",
            shell=True, check=True, capture_output=True, text=True,
        )
        subprocess.run(
            f"bgzip -c {output_vcf} > {output_vcf_gz}",
            shell=True, check=True, capture_output=True, text=True,
        )
        subprocess.run(
            f"tabix -f {output_vcf_gz}",
            shell=True, check=True, capture_output=True, text=True,
        )

        return Ok(output_vcf_gz)

    except subprocess.CalledProcessError as e:
        logger.warning(
            f"Phasing failed for {region_name} ({chrom}:{start}-{end}): "
            f"{e.stderr or 'possibly no variants in region'}"
        )
        return Err(f"Phasing failed for {region_name}")


# ---------------------------------------------------------------------------
# Haplotype extraction from phased VCF
# ---------------------------------------------------------------------------

def _apply_variant(
    template: list,
    pos_rel: int,
    ref: str,
    allele: str,
) -> None:
    """
    Apply a variant to a template sequence (list of characters) IN PLACE.

    Handles SNPs, insertions, and deletions.

    Args:
        template: Mutable list of single characters
        pos_rel: 0-based position relative to the template start
        ref: Reference allele string
        allele: Alternative allele string
    """
    ref_len = len(ref)
    allele_len = len(allele)

    if ref_len == allele_len == 1:
        # SNP
        template[pos_rel] = allele
    elif ref_len > allele_len:
        # Deletion: blank ref positions, put allele at start
        for i in range(pos_rel, min(pos_rel + ref_len, len(template))):
            template[i] = ""
        template[pos_rel] = allele
    elif allele_len > ref_len:
        # Insertion: blank ref positions, insert full allele
        for i in range(pos_rel, min(pos_rel + ref_len, len(template))):
            template[i] = ""
        template[pos_rel] = allele


def extract_haplotypes(
    phased_vcf_path: Path,
    reference_fasta_path: Path,
    chrom: str,
    start: int,
    end: int,
    min_freq: float = 0.0,
    region_id: Optional[str] = None,
) -> Result[HaplotypeSet, str]:
    """
    Extract unique haplotype sequences from a phased VCF for a genomic region.

    For each sample, constructs two haplotype sequences (hap0 / hap1)
    by applying phased genotypes to the reference. Then collects all unique
    haplotypes with their population frequencies.

    This produces REAL haplotypes from phased data — no artificial recombination.

    Args:
        phased_vcf_path: Path to phased VCF (bgzipped + tabix indexed)
        reference_fasta_path: Path to reference genome FASTA (with .fai)
        chrom: Chromosome name
        start: Region start (0-based)
        end: Region end (0-based exclusive)
        min_freq: Minimum frequency to retain a haplotype (0.0 = keep all)
        region_id: Optional label for this region

    Returns:
        HaplotypeSet with unique alleles and their frequencies
    """
    try:
        ref_fasta = pysam.FastaFile(str(reference_fasta_path))
    except Exception as e:
        return Err(f"Failed to open reference FASTA: {e}")

    try:
        vcf = pysam.VariantFile(str(phased_vcf_path))
    except Exception as e:
        ref_fasta.close()
        return Err(f"Failed to open VCF: {e}")

    rid = region_id or f"{chrom}:{start}-{end}"

    try:
        # Reference sequence for this region
        ref_seq = ref_fasta.fetch(chrom, start, end)

        sample_list = list(vcf.header.samples)
        if not sample_list:
            return Err("No samples found in VCF")

        # Per-sample haplotype templates (diploid)
        sample_haplotypes: Dict[str, Dict[int, list]] = {
            sample: {0: list(ref_seq), 1: list(ref_seq)}
            for sample in sample_list
        }

        # Apply variants
        n_variants = 0
        for record in vcf.fetch(chrom, start, end):
            n_variants += 1
            pos_rel = record.pos - 1 - start  # 0-based relative

            ref = record.ref
            alts = record.alts
            if alts is None:
                continue
            alleles = [ref] + list(alts)

            for sample in sample_list:
                gt = record.samples[sample]["GT"]
                if gt is None or len(gt) < 2:
                    continue
                if gt[0] is None or gt[1] is None:
                    continue

                gt0, gt1 = gt[0], gt[1]

                if gt0 != 0:
                    _apply_variant(
                        sample_haplotypes[sample][0],
                        pos_rel, ref, alleles[gt0],
                    )
                if gt1 != 0:
                    _apply_variant(
                        sample_haplotypes[sample][1],
                        pos_rel, ref, alleles[gt1],
                    )

        # Join character lists → sequences
        all_seqs: List[str] = []
        for sample in sample_list:
            for hap_idx in (0, 1):
                seq = "".join(sample_haplotypes[sample][hap_idx]).upper()
                all_seqs.append(seq)

        # Count unique haplotypes
        counter = Counter(all_seqs)
        total = len(all_seqs)

        alleles_out: Dict[str, str] = {}
        frequencies: Dict[str, float] = {}
        allele_num = 1

        for seq, count in counter.most_common():
            freq = count / total
            if freq >= min_freq:
                label = f"H{allele_num}"
                alleles_out[label] = seq
                frequencies[label] = freq
                allele_num += 1

        return Ok(HaplotypeSet(
            region_id=rid,
            alleles=alleles_out,
            frequencies=frequencies,
            n_samples=len(sample_list),
            n_variants=n_variants,
        ))

    finally:
        vcf.close()
        ref_fasta.close()


# ---------------------------------------------------------------------------
# Parallel worker for VCF haplotype extraction per BED region
# ---------------------------------------------------------------------------

def _worker_extract_region(
    region_tuple: Tuple[str, int, int, str],
    reference_path: str,
    vcf_path: str,
    phase: bool,
    temp_dir: str,
    min_freq: float,
) -> Tuple[str, Optional[HaplotypeSet]]:
    """
    Worker function for parallel haplotype extraction.

    Each worker opens its own file handles (pysam objects are not picklable).

    Args:
        region_tuple: (chrom, start, end, region_name)
        reference_path: Path to reference FASTA
        vcf_path: Path to input VCF
        phase: Whether to run shapeit5 first
        temp_dir: Temporary directory for phasing output
        min_freq: Minimum haplotype frequency

    Returns:
        (region_name, HaplotypeSet or None)
    """
    chrom, start, end, region_name = region_tuple

    # Phase if requested
    if phase:
        phase_result = phase_vcf_region(
            Path(vcf_path), chrom, start, end,
            Path(temp_dir), region_name,
        )
        if phase_result.is_err():
            logger.warning(f"No variants in {region_name}, using reference")
            return region_name, None
        vcf_to_use = phase_result.unwrap()
    else:
        vcf_to_use = Path(vcf_path)

    # Extract haplotypes
    result = extract_haplotypes(
        phased_vcf_path=vcf_to_use,
        reference_fasta_path=Path(reference_path),
        chrom=chrom, start=start, end=end,
        min_freq=min_freq,
        region_id=region_name,
    )

    if result.is_ok():
        return region_name, result.unwrap()
    else:
        logger.warning(f"Haplotype extraction failed for {region_name}: {result.unwrap_err()}")
        return region_name, None


# ---------------------------------------------------------------------------
# Multiple sequence alignment
# ---------------------------------------------------------------------------

def run_msa(
    sequences: Dict[str, str],
    aligner: str = "muscle",
    temp_dir: Optional[Path] = None,
) -> Result[Dict[str, str], str]:
    """
    Run multiple sequence alignment on a set of sequences.

    Supports muscle (v5), clustalo (Clustal Omega), and mafft.

    Args:
        sequences: sequence_id -> DNA sequence
        aligner: Aligner command ('muscle', 'clustalo', 'mafft')
        temp_dir: Temporary directory (created if None)

    Returns:
        Dict of sequence_id -> aligned sequence (with '-' gaps)
    """
    if len(sequences) <= 1:
        return Ok(dict(sequences))

    cleanup_dir = False
    if temp_dir is None:
        temp_dir = Path(tempfile.mkdtemp(prefix="eprobe_msa_"))
        cleanup_dir = True

    input_fa = temp_dir / "msa_input.fa"
    output_fa = temp_dir / "msa_output.fa"

    try:
        # Write input FASTA
        wr = _write_fasta_to_file(sequences, input_fa, line_width=0)
        if wr.is_err():
            return Err(f"Failed to write MSA input: {wr.unwrap_err()}")

        # Build command
        if aligner == "muscle":
            cmd = f"muscle -align {input_fa} -output {output_fa}"
        elif aligner == "clustalo":
            cmd = f"clustalo -i {input_fa} -o {output_fa} --force"
        elif aligner == "mafft":
            cmd = f"mafft --auto --quiet {input_fa} > {output_fa}"
        else:
            return Err(f"Unknown aligner: {aligner}. Use 'muscle', 'clustalo', or 'mafft'.")

        subprocess.run(
            cmd, shell=True, check=True,
            capture_output=True, text=True,
        )

        # Parse aligned output
        from Bio import SeqIO

        aligned: Dict[str, str] = {}
        for record in SeqIO.parse(str(output_fa), "fasta"):
            aligned[record.id] = str(record.seq)

        if not aligned:
            return Err("MSA produced no output sequences")

        return Ok(aligned)

    except subprocess.CalledProcessError as e:
        return Err(f"MSA failed ({aligner}): {e.stderr or str(e)}")
    except Exception as e:
        return Err(f"MSA error: {e}")
    finally:
        for f in [input_fa, output_fa]:
            if f.exists():
                f.unlink()
        if cleanup_dir and temp_dir.exists():
            try:
                temp_dir.rmdir()
            except OSError:
                pass


# ---------------------------------------------------------------------------
# Haplotype-aware probe generation
# ---------------------------------------------------------------------------

def generate_haplotype_probes(
    region_id: str,
    alleles: Dict[str, str],
    probe_length: int,
    step_size: int,
    variant_only: bool = False,
    source_chrom: Optional[str] = None,
    source_offset: int = 0,
    aligner: str = "mafft",
) -> Result[Tuple[List[Probe], Dict[str, int]], str]:
    """
    Generate probes with haplotype awareness.

    Strategy:
    - 1 allele → simple tiling (or empty if variant_only)
    - N alleles, all same length → fast path: tile each at same positions, dedup per window
    - N alleles, different lengths → MSA path: align, slide window, extract originals

    At non-variant windows: 1 probe (from any allele)
    At variant windows: 1 probe per unique haplotype sequence
    Variant probe IDs include mutation annotation, e.g.:
      GENE_H1_snp87:A_P0002   (SNP at gene pos 87, this haplotype has A)
      GENE_H2_snp87:G_P0003   (same window, H2 has G)

    Args:
        region_id: Identifier for this region/gene
        alleles: allele_label -> DNA sequence
        probe_length: Probe length (bp)
        step_size: Sliding window step (bp)
        variant_only: If True, only generate probes at variant positions
        source_chrom: Chromosome for coordinate annotation
        source_offset: Genomic offset (0-based) for coordinate mapping
        aligner: MSA aligner when needed ('mafft', 'muscle', 'clustalo')

    Returns:
        Ok((List[Probe], stats_dict)) where stats_dict has:
          n_variant_windows:   windows with >=2 unique probe sequences
          n_invariant_windows: windows with exactly 1 unique probe sequence
          n_variant_probes:    total probes from variant windows
          n_invariant_probes:  total probes from invariant windows
    """
    empty_stats: Dict[str, int] = {
        "n_variant_windows": 0, "n_invariant_windows": 0,
        "n_variant_probes": 0, "n_invariant_probes": 0,
    }

    if not alleles:
        return Ok(([], empty_stats))

    # --- Single allele: simple tiling ---
    if len(alleles) == 1:
        if variant_only:
            return Ok(([], empty_stats))
        seq = list(alleles.values())[0]
        probes = tile_sequence(
            region_id, seq, probe_length, step_size,
            source_chrom, source_offset,
        )
        stats = {
            "n_variant_windows": 0, "n_invariant_windows": len(probes),
            "n_variant_probes": 0, "n_invariant_probes": len(probes),
        }
        return Ok((probes, stats))

    # --- Multiple alleles ---
    lengths = set(len(s) for s in alleles.values())

    if len(lengths) == 1:
        return _probes_equal_length(
            region_id, alleles, probe_length, step_size,
            variant_only, source_chrom, source_offset,
        )
    else:
        return _probes_with_msa(
            region_id, alleles, probe_length, step_size,
            variant_only, source_chrom, source_offset, aligner,
        )


def _probes_equal_length(
    region_id: str,
    alleles: Dict[str, str],
    probe_length: int,
    step_size: int,
    variant_only: bool,
    source_chrom: Optional[str],
    source_offset: int,
) -> Result[Tuple[List[Probe], Dict[str, int]], str]:
    """
    Fast path: all alleles have equal length (SNPs only).

    Tiles all alleles at identical positions and deduplicates per window.
    Variant probe IDs encode the variable positions and bases:
      {gene}_{allele}_snp{pos1}:{base}[,{pos2}:{base}...]_P{n}
    where pos is 1-based gene-relative.
    """
    seq_len = len(next(iter(alleles.values())))

    if seq_len < probe_length:
        return Ok([])

    positions = list(range(0, seq_len - probe_length + 1, step_size))
    last_start = seq_len - probe_length
    if positions and positions[-1] < last_start:
        positions.append(last_start)

    allele_items = list(alleles.items())
    all_probes: List[Probe] = []
    probe_counter = 0
    n_variant_windows = 0
    n_invariant_windows = 0

    for start in positions:
        end = start + probe_length

        # Collect unique sequences at this window
        seen: Dict[str, str] = {}  # sequence -> first allele_label
        for allele_label, seq in allele_items:
            window_seq = seq[start:end].upper()
            if "N" in window_seq:
                continue
            if window_seq not in seen:
                seen[window_seq] = allele_label

        if not seen:
            continue

        is_variant = len(seen) > 1

        if variant_only and not is_variant:
            continue

        # For variant windows: find variable positions in gene coordinates
        # and build per-allele mutation annotation for the probe ID
        if is_variant:
            n_variant_windows += 1
            # Find positions (0-based within window) that differ across alleles
            unique_seqs = list(seen.keys())
            var_positions = [
                p for p in range(probe_length)
                if len({sq[p] for sq in unique_seqs}) > 1
            ]
        else:
            n_invariant_windows += 1
            var_positions = []

        chrom = source_chrom or region_id
        genomic_start = source_offset + start + 1
        genomic_end = source_offset + end

        for probe_seq, allele_label in seen.items():
            probe_counter += 1
            if is_variant:
                # Encode variable positions as snp{gene_pos1_1based}:{base},...
                parts = [
                    f"{start + p + 1}:{probe_seq[p]}"
                    for p in var_positions[:5]  # cap at 5 to keep ID readable
                ]
                extra = f"+{len(var_positions) - 5}more" if len(var_positions) > 5 else ""
                ann = "snp" + ",".join(parts) + extra
                probe_id = f"{region_id}_{allele_label}_{ann}_P{probe_counter:04d}"
            else:
                probe_id = f"{region_id}_P{probe_counter:04d}"

            probe = Probe(
                id=probe_id,
                sequence=probe_seq,
                source_chrom=chrom,
                source_start=genomic_start,
                source_end=genomic_end,
            )
            all_probes.append(probe)

    n_variant_probes = sum(1 for p in all_probes if "_snp" in p.id)
    n_invariant_probes = len(all_probes) - n_variant_probes
    stats = {
        "n_variant_windows": n_variant_windows,
        "n_invariant_windows": n_invariant_windows,
        "n_variant_probes": n_variant_probes,
        "n_invariant_probes": n_invariant_probes,
    }
    return Ok((all_probes, stats))


def _find_ungapped_position(aligned_window: str, original_seq: str) -> int:
    """Map aligned window back to position in the original (ungapped) sequence."""
    ungapped = aligned_window.replace("-", "").upper()
    if not ungapped:
        return -1
    return original_seq.upper().find(ungapped)


def _probes_with_msa(
    region_id: str,
    alleles: Dict[str, str],
    probe_length: int,
    step_size: int,
    variant_only: bool,
    source_chrom: Optional[str],
    source_offset: int,
    aligner: str,
) -> Result[Tuple[List[Probe], Dict[str, int]], str]:
    """
    MSA path: alleles have different lengths (indels present).

    1. Align all alleles
    2. Slide window across alignment coordinates
    3. For each window, extract original (ungapped) sequences per allele
    4. Deduplicate per window

    Variant probe IDs encode the allele + 'indel' type:
      {gene}_{allele}_indel_P{n}
    """
    msa_result = run_msa(alleles, aligner=aligner)
    if msa_result.is_err():
        return Err(f"MSA failed for {region_id}: {msa_result.unwrap_err()}")

    aligned = msa_result.unwrap()
    aln_length = len(next(iter(aligned.values())))

    if aln_length < probe_length:
        return Ok(([], {"n_variant_windows": 0, "n_invariant_windows": 0,
                        "n_variant_probes": 0, "n_invariant_probes": 0}))

    positions = list(range(0, aln_length - probe_length + 1, step_size))
    last_start = aln_length - probe_length
    if positions and positions[-1] < last_start:
        positions.append(last_start)

    all_probes: List[Probe] = []
    probe_counter = 0
    n_variant_windows = 0
    n_invariant_windows = 0

    for aln_start in positions:
        aln_end = aln_start + probe_length

        # Check variation in this alignment window
        has_variation = False
        for col in range(aln_start, aln_end):
            bases = set()
            for aln_seq in aligned.values():
                b = aln_seq[col]
                if b != "-":
                    bases.add(b.upper())
            if len(bases) > 1:
                has_variation = True
                break

        if variant_only and not has_variation:
            continue

        if has_variation:
            n_variant_windows += 1
        else:
            n_invariant_windows += 1

        # Extract original (ungapped) probe sequences per allele
        seen: Dict[str, str] = {}
        for allele_label in aligned:
            aln_window = aligned[allele_label][aln_start:aln_end]
            original_seq = alleles[allele_label]

            # Map back to original sequence
            orig_pos = _find_ungapped_position(aln_window, original_seq)

            if orig_pos >= 0 and orig_pos + probe_length <= len(original_seq):
                probe_seq = original_seq[orig_pos:orig_pos + probe_length].upper()
            else:
                # Fallback: use ungapped aligned window
                probe_seq = aln_window.replace("-", "").upper()

            if len(probe_seq) != probe_length:
                continue
            if "N" in probe_seq:
                continue

            if probe_seq not in seen:
                seen[probe_seq] = allele_label

        if not seen:
            continue

        is_variant = len(seen) > 1

        for probe_seq, allele_label in seen.items():
            probe_counter += 1
            if is_variant:
                probe_id = f"{region_id}_{allele_label}_indel_P{probe_counter:04d}"
            else:
                probe_id = f"{region_id}_P{probe_counter:04d}"
            probe = Probe(
                id=probe_id,
                sequence=probe_seq,
                source_chrom=source_chrom or region_id,
                source_start=source_offset + aln_start + 1,
                source_end=source_offset + aln_end,
            )
            all_probes.append(probe)

    n_variant_probes = sum(1 for p in all_probes if "_indel_" in p.id)
    n_invariant_probes = len(all_probes) - n_variant_probes
    stats = {
        "n_variant_windows": n_variant_windows,
        "n_invariant_windows": n_invariant_windows,
        "n_variant_probes": n_variant_probes,
        "n_invariant_probes": n_invariant_probes,
    }
    return Ok((all_probes, stats))
