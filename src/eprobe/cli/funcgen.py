"""
FUNCGEN panel CLI commands.

Commands for sequence-based / haplotype-aware probe design:
  from-fasta  - Generate probes from FASTA sequences (simple or haplotype mode)
  from-bed    - Generate probes from BED regions + reference (optional VCF haplotyping)
  from-gff    - Generate probes from GFF annotations + reference (optional VCF haplotyping)
  assess      - Quality assessment of generated probes
"""

import click
from pathlib import Path
from typing import Optional

from eprobe.cli.utils import (
    echo_success,
    echo_error,
    echo_info,
)
from eprobe.cli.main import AliasedGroup


# ---------------------------------------------------------------------------
# Group
# ---------------------------------------------------------------------------

@click.group(cls=AliasedGroup)
@click.pass_context
def funcgen(ctx: click.Context) -> None:
    """
    Functional genetics panel for sequence-based probe design.

    Generate probes from gene sequences, BED regions, or GFF annotations.
    Supports haplotype-aware probe design for capturing multiple alleles.

    \b
    Workflows:
      from-fasta  Multi-FASTA input (simple tiling or haplotype-aware)
      from-bed    BED regions + reference (+ optional VCF haplotyping)
      from-gff    GFF annotations → BED → probe generation
      assess      Quality assessment of generated probes

    \b
    Examples:
      eprobe funcgen from-fasta -f genes.fa -l 81 -s 30 -o probes
      eprobe funcgen from-fasta -f alleles.fa -l 81 -s 30 --haplotyping -o probes
      eprobe funcgen from-bed -b regions.bed -r ref.fa -l 81 -s 30 -o probes
      eprobe funcgen from-bed -b regions.bed -r ref.fa -v phased.vcf.gz -l 81 -s 30 -o probes
      eprobe funcgen from-gff -g genes.gff -r ref.fa --feature_type CDS -l 81 -s 30 -o probes
    """
    pass


# ---------------------------------------------------------------------------
# from-fasta
# ---------------------------------------------------------------------------

@funcgen.command(name="from-fasta")
@click.option(
    "-f", "--fasta",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input FASTA file with gene/region sequences.",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix.",
)
@click.option(
    "-l", "--length",
    required=True,
    type=int,
    help="Probe length (bp).",
)
@click.option(
    "-s", "--step",
    required=True,
    type=int,
    help="Sliding window step size (bp).",
)
@click.option(
    "--haplotyping/--no_haplotyping",
    default=False,
    help="Enable haplotype-aware probe design. Requires allele naming "
         "convention: GeneID{sep}AlleleNum (e.g. BRCA1_1, BRCA1_2). Default: off.",
)
@click.option(
    "--haplotype_sep",
    type=str,
    default="_",
    help="Separator for haplotype IDs (e.g. Gene_1, Gene_2). Default: '_'.",
)
@click.option(
    "--variant_only",
    is_flag=True,
    default=False,
    help="Only generate probes at variant positions (skip shared/non-variant windows).",
)
@click.option(
    "--aligner",
    type=click.Choice(["muscle", "clustalo", "mafft"]),
    default="mafft",
    help="MSA aligner for haplotype comparison when indels present. Default: mafft.",
)
@click.pass_context
def from_fasta(
    ctx: click.Context,
    fasta: Path,
    output: Path,
    length: int,
    step: int,
    haplotyping: bool,
    haplotype_sep: str,
    variant_only: bool,
    aligner: str,
) -> None:
    """
    Generate probes from FASTA sequences.

    \b
    Simple mode (default):
      Tiles probes across each input sequence with a sliding window.

    \b
    Haplotype mode (--haplotyping):
      Sequences are grouped by gene ID using the naming convention
      GeneID{sep}AlleleNum (e.g., BRCA1_1, BRCA1_2).
      For each gene:
        - Non-variant positions → 1 shared probe
        - Variant positions → 1 probe per unique haplotype
      Use --variant_only to skip non-variant positions.

    \b
    Output files:
      {output}.probes.fa    - Probe sequences in FASTA format
      {output}.probes.tsv   - Probe metadata table
    """
    from eprobe.funcgen.from_fasta import run_from_fasta

    verbose = ctx.obj.get("verbose", False)

    echo_info(f"Generating probes from {fasta}")
    echo_info(f"Length: {length}bp, Step: {step}bp")

    if haplotyping:
        echo_info(f"Haplotype mode: separator='{haplotype_sep}', aligner={aligner}")
        if variant_only:
            echo_info("  → variant-only probes")

    result = run_from_fasta(
        fasta_path=fasta,
        output_prefix=output,
        probe_length=length,
        step_size=step,
        haplotyping=haplotyping,
        haplotype_sep=haplotype_sep,
        variant_only=variant_only,
        aligner=aligner,
        verbose=verbose,
    )

    if result.is_err():
        echo_error(f"Generation failed: {result.unwrap_err()}")
        raise SystemExit(1)

    stats = result.unwrap()
    echo_success(
        f"Generated {stats['probe_count']} probes from {stats['gene_count']} sequences"
    )
    if stats.get("multi_allele_count", 0) > 0:
        echo_info(f"  Multi-allele genes: {stats['multi_allele_count']}")
    if stats.get("haplotyping") and stats.get("multi_allele_count", 0) > 0:
        v_win = stats.get("n_variant_windows_total", 0)
        i_win = stats.get("n_invariant_windows_total", 0)
        v_pr = stats.get("n_variant_probes_total", 0)
        i_pr = stats.get("n_invariant_probes_total", 0)
        if v_win + i_win > 0:
            echo_info(f"  Variant windows   (mutation-covering): {v_win} → {v_pr} probes")
            echo_info(f"  Invariant windows (shared sequence):   {i_win} → {i_pr} probes")


# ---------------------------------------------------------------------------
# from-bed
# ---------------------------------------------------------------------------

@funcgen.command(name="from-bed")
@click.option(
    "-b", "--bed",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input BED file with target regions.",
)
@click.option(
    "-r", "--reference",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Reference genome FASTA file (must have .fai index).",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix.",
)
@click.option(
    "-l", "--length",
    required=True,
    type=int,
    help="Probe length (bp).",
)
@click.option(
    "-s", "--step",
    required=True,
    type=int,
    help="Sliding window step size (bp).",
)
@click.option(
    "-v", "--vcf",
    type=click.Path(exists=True, path_type=Path),
    help="VCF file for haplotype inference (bgzipped + tabix indexed).",
)
@click.option(
    "--phase/--no_phase",
    default=False,
    help="Phase VCF with shapeit5 before haplotype extraction. "
         "Use --no_phase if VCF is already phased. Default: off.",
)
@click.option(
    "--min_freq",
    type=float,
    default=0.05,
    help="Minimum haplotype frequency to retain (default: 0.05).",
)
@click.option(
    "--variant_only",
    is_flag=True,
    default=False,
    help="Only generate probes at variant positions.",
)
@click.option(
    "--aligner",
    type=click.Choice(["muscle", "clustalo", "mafft"]),
    default="muscle",
    help="MSA aligner for indel-containing haplotypes. Default: muscle.",
)
@click.option(
    "-t", "--threads",
    default=1,
    type=int,
    help="Number of threads (default: 1).",
)
@click.pass_context
def from_bed(
    ctx: click.Context,
    bed: Path,
    reference: Path,
    output: Path,
    length: int,
    step: int,
    vcf: Optional[Path],
    phase: bool,
    min_freq: float,
    variant_only: bool,
    aligner: str,
    threads: int,
) -> None:
    """
    Generate probes from BED regions.

    \b
    Simple mode (no --vcf):
      Extracts reference sequences at each BED region, tiles probes.

    \b
    Haplotype mode (--vcf):
      For each BED region:
        1. Phase VCF with shapeit5 (if --phase)
        2. Extract per-sample haplotype sequences from phased genotypes
        3. Filter by frequency (--min_freq)
        4. Generate probes: 1 at shared positions, N at variant positions
      Only REAL haplotypes from phased data are used (no artificial recombination).

    \b
    Output files:
      {output}.probes.fa    - Probe sequences
      {output}.probes.tsv   - Probe metadata
    """
    from eprobe.funcgen.from_bed import run_from_bed

    verbose = ctx.obj.get("verbose", False)

    if phase and vcf is None:
        echo_error("--vcf is required when using --phase")
        raise SystemExit(1)

    echo_info(f"Generating probes from {bed}")
    echo_info(f"Reference: {reference}")
    echo_info(f"Length: {length}bp, Step: {step}bp")

    if vcf:
        echo_info(f"Haplotype mode: VCF={vcf}, phase={phase}, min_freq={min_freq}")
        if variant_only:
            echo_info("  → variant-only probes")

    result = run_from_bed(
        bed_path=bed,
        reference_path=reference,
        output_prefix=output,
        probe_length=length,
        step_size=step,
        vcf_path=vcf,
        phase=phase,
        min_freq=min_freq,
        variant_only=variant_only,
        aligner=aligner,
        threads=threads,
        verbose=verbose,
    )

    if result.is_err():
        echo_error(f"Generation failed: {result.unwrap_err()}")
        raise SystemExit(1)

    stats = result.unwrap()
    echo_success(
        f"Generated {stats['probe_count']} probes from {stats['region_count']} regions"
    )
    if stats.get("haplo_regions", 0) > 0:
        echo_info(f"  Regions with multiple haplotypes: {stats['haplo_regions']}")


# ---------------------------------------------------------------------------
# from-gff
# ---------------------------------------------------------------------------

@funcgen.command(name="from-gff")
@click.option(
    "-g", "--gff",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input GFF/GFF3/GTF annotation file.",
)
@click.option(
    "-r", "--reference",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Reference genome FASTA file (must have .fai index).",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix.",
)
@click.option(
    "-l", "--length",
    required=True,
    type=int,
    help="Probe length (bp).",
)
@click.option(
    "-s", "--step",
    required=True,
    type=int,
    help="Sliding window step size (bp).",
)
@click.option(
    "--feature_type",
    type=str,
    default="CDS",
    help="GFF feature type(s) to extract. Comma-separated for multiple "
         "(e.g. 'CDS,exon'). Default: CDS.",
)
@click.option(
    "--gene_list",
    type=click.Path(exists=True, path_type=Path),
    help="File with gene IDs to extract (one per line).",
)
@click.option(
    "--gene_id",
    type=str,
    help="Gene ID(s) to extract (comma-separated for multiple).",
)
@click.option(
    "-v", "--vcf",
    type=click.Path(exists=True, path_type=Path),
    help="VCF file for haplotype inference (bgzipped + tabix indexed).",
)
@click.option(
    "--phase/--no_phase",
    default=False,
    help="Phase VCF with shapeit5. Default: off.",
)
@click.option(
    "--min_freq",
    type=float,
    default=0.05,
    help="Minimum haplotype frequency (default: 0.05).",
)
@click.option(
    "--variant_only",
    is_flag=True,
    default=False,
    help="Only generate probes at variant positions.",
)
@click.option(
    "--aligner",
    type=click.Choice(["muscle", "clustalo", "mafft"]),
    default="muscle",
    help="MSA aligner for indel-containing haplotypes. Default: muscle.",
)
@click.option(
    "-t", "--threads",
    default=1,
    type=int,
    help="Number of threads (default: 1).",
)
@click.pass_context
def from_gff(
    ctx: click.Context,
    gff: Path,
    reference: Path,
    output: Path,
    length: int,
    step: int,
    feature_type: str,
    gene_list: Optional[Path],
    gene_id: Optional[str],
    vcf: Optional[Path],
    phase: bool,
    min_freq: float,
    variant_only: bool,
    aligner: str,
    threads: int,
) -> None:
    """
    Generate probes from GFF annotations.

    \b
    Workflow:
      1. Parse GFF, filter by --feature_type
      2. Optionally filter by gene IDs (--gene_list or --gene_id)
      3. Merge overlapping features by gene → BED regions
      4. Extract sequences + tile probes (or haplotype-aware with --vcf)

    \b
    Feature types:
      gene   Full gene regions including introns
      mRNA   Transcribed regions
      CDS    Coding sequences (recommended for protein-coding genes)
      exon   Individual exons

    \b
    Output files:
      {output}.probes.fa    - Probe sequences
      {output}.probes.tsv   - Probe metadata
    """
    from eprobe.funcgen.from_gff import run_from_gff

    verbose = ctx.obj.get("verbose", False)

    # Parse gene IDs
    gene_ids = None
    if gene_list:
        with open(gene_list) as f:
            gene_ids = [line.strip() for line in f if line.strip()]
    elif gene_id:
        gene_ids = [g.strip() for g in gene_id.split(",")]

    echo_info(f"Generating probes from {gff}")
    echo_info(f"Feature type: {feature_type}")
    echo_info(f"Length: {length}bp, Step: {step}bp")

    if gene_ids:
        echo_info(f"Target genes: {len(gene_ids)}")
    if vcf:
        echo_info(f"Haplotype mode: VCF={vcf}, phase={phase}, min_freq={min_freq}")

    result = run_from_gff(
        gff_path=gff,
        reference_path=reference,
        output_prefix=output,
        probe_length=length,
        step_size=step,
        feature_type=feature_type,
        gene_ids=gene_ids,
        vcf_path=vcf,
        phase=phase,
        min_freq=min_freq,
        variant_only=variant_only,
        aligner=aligner,
        threads=threads,
        verbose=verbose,
    )

    if result.is_err():
        echo_error(f"Generation failed: {result.unwrap_err()}")
        raise SystemExit(1)

    stats = result.unwrap()
    echo_success(
        f"Generated {stats['probe_count']} probes from {stats.get('gene_count', stats['region_count'])} genes"
    )
    if stats.get("haplo_regions", 0) > 0:
        echo_info(f"  Regions with multiple haplotypes: {stats['haplo_regions']}")


# ---------------------------------------------------------------------------
# assess
# ---------------------------------------------------------------------------

@funcgen.command(name="assess")
@click.option(
    "-f", "--fasta",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input probe FASTA file.",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix.",
)
@click.option(
    "--tags",
    type=str,
    default="gc,tm,complexity,hairpin,dimer",
    help="Tags to analyze (comma-separated). Default: gc,tm,complexity,hairpin,dimer.",
)
@click.option(
    "--plot/--no_plot",
    default=True,
    help="Generate distribution plots (default: yes).",
)
@click.option(
    "-t", "--threads",
    default=1,
    type=int,
    help="Number of threads (default: 1).",
)
@click.pass_context
def assess(
    ctx: click.Context,
    fasta: Path,
    output: Path,
    tags: str,
    plot: bool,
    threads: int,
) -> None:
    """
    Assess quality of FUNCGEN probe set.

    Calculates biophysical properties and generates summary statistics.

    \b
    Available tags:
      gc         - GC content (%)
      tm         - Melting temperature (°C)
      complexity - DUST complexity score
      hairpin    - Self-complementarity score
      dimer      - Inter-probe complementarity score

    \b
    Output files:
      {output}.stats.tsv    - Per-probe statistics
      {output}.summary.txt  - Summary statistics
      {output}.*.png        - Distribution plots (if --plot)
    """
    from eprobe.funcgen.assess import run_assess

    verbose = ctx.obj.get("verbose", False)

    echo_info(f"Assessing probe set: {fasta}")
    echo_info(f"Tags: {tags}")

    result = run_assess(
        fasta_path=fasta,
        output_prefix=output,
        tags=tags.split(","),
        generate_plots=plot,
        threads=threads,
        verbose=verbose,
    )

    if result.is_err():
        echo_error(f"Assessment failed: {result.unwrap_err()}")
        raise SystemExit(1)

    stats = result.unwrap()
    echo_success(f"Assessed {stats['probe_count']} probes")

    # Print summary
    echo_info("Summary statistics:")
    for tag, values in stats.get("summary", {}).items():
        echo_info(f"  {tag}: mean={values['mean']:.2f}, std={values['std']:.2f}")
