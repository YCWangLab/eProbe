"""
FUNCGEN panel CLI commands.

Commands for sequence-based probe design workflow:
  from_fasta  - Generate probes from FASTA sequences
  from_bed    - Generate probes from BED regions
  from_gff    - Generate probes from GFF annotations
  assess      - Quality assessment
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


@click.group(cls=AliasedGroup)
@click.pass_context
def funcgen(ctx: click.Context) -> None:
    """
    Functional genetics panel for sequence-based probe design.
    
    Generate probes from gene sequences, BED regions, or GFF annotations.
    Supports haplotype-aware probe design for capturing multiple alleles.
    
    \b
    Input options:
      from_fasta  - Multi-FASTA with gene/allele sequences
      from_bed    - BED regions + reference genome
      from_gff    - GFF annotations + reference genome
    
    \b
    Example:
      eprobe funcgen from_fasta -f genes.fa -l 81 -s 30 -o probes
      eprobe funcgen from_bed -b regions.bed -r ref.fa -l 81 -s 30 -o probes
    """
    pass


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
    "--haplotype_sep",
    type=str,
    default="_",
    help="Separator for haplotype IDs (e.g., Gene_1, Gene_2). Default: '_'",
)
@click.option(
    "--min_freq",
    type=float,
    default=0.05,
    help="Minimum allele frequency to retain (default: 0.05).",
)
@click.option(
    "--haplotyping/--no_haplotyping",
    default=False,
    help="Enable haplotype-aware probe design (default: no).",
)
@click.pass_context
def from_fasta(
    ctx: click.Context,
    fasta: Path,
    output: Path,
    length: int,
    step: int,
    haplotype_sep: str,
    min_freq: float,
    haplotyping: bool,
) -> None:
    """
    Generate probes from FASTA sequences.
    
    Uses a sliding window approach to tile probes across each input sequence.
    
    \b
    Haplotype-aware mode (--haplotyping):
      When sequences represent different alleles of the same gene,
      use naming convention: GeneID{sep}AlleleNum (e.g., BRCA1_1, BRCA1_2)
      The tool will generate probes considering variant positions.
    
    \b
    Output files:
      {output}.probes.fa    - Probe sequences
      {output}.probes.tsv   - Probe metadata
    """
    from eprobe.funcgen.from_fasta import run_from_fasta
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Generating probes from {fasta}")
    echo_info(f"Length: {length}bp, Step: {step}bp")
    
    if haplotyping:
        echo_info(f"Haplotype mode enabled (separator: '{haplotype_sep}')")
    
    result = run_from_fasta(
        fasta_path=fasta,
        output_prefix=output,
        probe_length=length,
        step_size=step,
        haplotype_sep=haplotype_sep,
        min_freq=min_freq,
        haplotyping=haplotyping,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Generation failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Generated {stats['probe_count']} probes from {stats['gene_count']} sequences")


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
    help="Reference genome FASTA file.",
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
    help="VCF file for haplotype inference (optional).",
)
@click.option(
    "--phase/--no_phase",
    default=False,
    help="Phase variants to infer haplotypes (requires --vcf).",
)
@click.option(
    "--min_freq",
    type=float,
    default=0.05,
    help="Minimum haplotype frequency to retain (default: 0.05).",
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
    threads: int,
) -> None:
    """
    Generate probes from BED regions.
    
    Extracts sequences from reference genome at BED-specified regions,
    then tiles probes across each region.
    
    \b
    Haplotype mode (--phase --vcf):
      Uses VCF to identify variants within each region.
      Phases variants using SHAPEIT to infer haplotypes.
      Generates probes for each common haplotype.
    
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
    
    if phase:
        echo_info(f"Haplotype phasing enabled with VCF: {vcf}")
    
    result = run_from_bed(
        bed_path=bed,
        reference_path=reference,
        output_prefix=output,
        probe_length=length,
        step_size=step,
        vcf_path=vcf,
        phase=phase,
        min_freq=min_freq,
        threads=threads,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Generation failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Generated {stats['probe_count']} probes from {stats['region_count']} regions")


@funcgen.command(name="from-gff")
@click.option(
    "-g", "--gff",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input GFF/GFF3 annotation file.",
)
@click.option(
    "-r", "--reference",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Reference genome FASTA file.",
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
    type=click.Choice(["gene", "mRNA", "CDS", "exon"]),
    default="CDS",
    help="GFF feature type to extract (default: CDS).",
)
@click.option(
    "--gene_list",
    type=click.Path(exists=True, path_type=Path),
    help="File with gene IDs to extract (one per line).",
)
@click.option(
    "--gene_id",
    type=str,
    help="Single gene ID to extract (comma-separated for multiple).",
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
    threads: int,
) -> None:
    """
    Generate probes from GFF annotations.
    
    Extracts specified feature types from GFF, retrieves sequences
    from reference genome, and tiles probes.
    
    \b
    Feature types:
      gene  - Full gene regions including introns
      mRNA  - Transcribed regions
      CDS   - Coding sequences (recommended for protein-coding genes)
      exon  - Individual exons
    
    \b
    Gene selection:
      --gene_list  - File with gene IDs (one per line)
      --gene_id    - Comma-separated gene IDs
      (if neither specified, processes all genes)
    
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
    
    result = run_from_gff(
        gff_path=gff,
        reference_path=reference,
        output_prefix=output,
        probe_length=length,
        step_size=step,
        feature_type=feature_type,
        gene_ids=gene_ids,
        threads=threads,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Generation failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Generated {stats['probe_count']} probes from {stats['gene_count']} genes")


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
    help="Tags to analyze (comma-separated).",
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
    
    Calculates biophysical properties for all probes and generates
    summary statistics and distribution plots.
    
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
