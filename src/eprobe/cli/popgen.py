"""
POPGEN panel CLI commands.

Commands for SNP-based probe design workflow:
  extract  - Extract SNPs from VCF file
  filter   - Apply multi-stage filtering
  select   - Window-based SNP selection
  build    - Generate probe sequences
  assess   - Quality assessment
"""

import click
from pathlib import Path
from typing import Optional

from eprobe.cli.utils import (
    validate_file_exists,
    validate_output_path,
    echo_success,
    echo_error,
    echo_info,
    echo_warning,
)
from eprobe.cli.main import AliasedGroup


@click.group(cls=AliasedGroup)
@click.pass_context
def popgen(ctx: click.Context) -> None:
    """
    Population genetics panel for SNP-based probe design.
    
    \b
    Workflow:
      1. extract  - Extract SNPs from VCF
      2. filter   - Apply quality filters
      3. select   - Select optimal SNPs per window
      4. build    - Generate probe sequences
      5. assess   - Evaluate probe set quality
    
    \b
    Example workflow:
      eprobe popgen extract -v input.vcf.gz -r ref.fa -o project/step1
      eprobe popgen filter -i project/step1.snps.tsv -r ref.fa -o project/step2
      eprobe popgen select -i project/step2.snps.tsv -r ref.fa -o project/step3
      eprobe popgen build -i project/step3.snps.tsv -r ref.fa -o project/probes
    """
    pass


@popgen.command()
@click.option(
    "-v", "--vcf",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input VCF file (compressed .vcf.gz with .tbi index).",
)
@click.option(
    "-r", "--reference",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Reference genome FASTA file (same as used for VCF calling).",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix for generated files.",
)
@click.option(
    "-t", "--threads",
    default=1,
    type=int,
    help="Number of threads (default: 1, recommended: number of chromosomes).",
)
@click.option(
    "--cluster_flank",
    default=60,
    type=int,
    help="Flanking distance (bp) for cluster detection (default: 60).",
)
@click.option(
    "--max_cluster_snp",
    default=3,
    type=int,
    help="Maximum SNPs allowed in flanking region (default: 3).",
)
@click.option(
    "--min_cluster_snp",
    default=1,
    type=int,
    help="Minimum SNPs required in flanking region (default: 1).",
)
@click.option(
    "--keep_bed",
    type=click.Path(exists=True, path_type=Path),
    help="BED file of regions to keep SNPs in.",
)
@click.option(
    "--remove_bed",
    type=click.Path(exists=True, path_type=Path),
    help="BED file of regions to remove SNPs from.",
)
@click.option(
    "--cluster_filter/--no_cluster_filter",
    default=True,
    help="Enable/disable cluster filtering (default: enabled).",
)
@click.pass_context
def extract(
    ctx: click.Context,
    vcf: Path,
    reference: Path,
    output: Path,
    threads: int,
    cluster_flank: int,
    max_cluster_snp: int,
    min_cluster_snp: int,
    keep_bed: Optional[Path],
    remove_bed: Optional[Path],
    cluster_filter: bool,
) -> None:
    """
    Extract SNPs from VCF file.
    
    Reads a compressed, indexed VCF file and extracts biallelic SNPs.
    Optionally filters out SNP clusters (regions with too many SNPs close together)
    as these indicate hypervariable regions that may cause probe design issues.
    
    \b
    Why filter SNP clusters?
      - Dense SNP regions indicate high variability
      - Probes from these regions may not hybridize consistently
      - Single SNP probes work best in conserved flanking contexts
    
    \b
    Output files:
      {output}.snps.tsv       - Extracted SNPs in TSV format
      {output}.chr_sizes.tsv  - Chromosome sizes
      {output}.extract.log    - Processing log
    """
    from eprobe.popgen.extract import run_extract
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Extracting SNPs from {vcf}")
    
    result = run_extract(
        vcf_path=vcf,
        reference_path=reference,
        output_prefix=output,
        threads=threads,
        cluster_flank=cluster_flank,
        max_cluster_snp=max_cluster_snp,
        min_cluster_snp=min_cluster_snp,
        keep_bed=keep_bed,
        remove_bed=remove_bed,
        enable_cluster_filter=cluster_filter,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Extraction failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Extracted {stats['total_snps']} SNPs")
    
    if cluster_filter:
        echo_info(f"  Filtered {stats.get('cluster_filtered', 0)} SNPs in clusters")
    
    echo_info(f"Output: {output}.snps.tsv")


@popgen.command()
@click.option(
    "-i", "--input",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input SNPs TSV file (from extract step).",
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
    "-t", "--threads",
    default=1,
    type=int,
    help="Number of threads (default: 1).",
)
@click.option(
    "-l", "--length",
    default=81,
    type=int,
    help="Probe length for filtering (default: 81).",
)
@click.option(
    "--bg_db",
    type=str,
    help="Kraken2 database(s) for background filtering (comma-separated).",
)
@click.option(
    "--ac_db",
    type=str,
    help="Bowtie2 database(s) for accessibility filtering (comma-separated).",
)
@click.option(
    "--tx_db",
    type=str,
    help="Bowtie2 database for taxonomic filtering.",
)
@click.option(
    "--tx_taxid",
    type=int,
    help="Target taxon ID for taxonomic filtering.",
)
@click.option(
    "--gc",
    type=str,
    default="35,65",
    help="GC content range (min,max percentage, default: 35,65).",
)
@click.option(
    "--tm",
    type=str,
    default="65,85",
    help="Melting temperature range (min,max Celsius, default: 65,85).",
)
@click.option(
    "--complexity",
    type=float,
    default=2.0,
    help="Maximum DUST complexity score (default: 2.0).",
)
@click.option(
    "--hairpin",
    type=float,
    default=30.0,
    help="Maximum hairpin score (default: 30.0).",
)
@click.option(
    "--dimer",
    type=float,
    default=50.0,
    help="Maximum dimer score (default: 50.0).",
)
@click.pass_context
def filter(
    ctx: click.Context,
    input: Path,
    reference: Path,
    output: Path,
    threads: int,
    length: int,
    bg_db: Optional[str],
    ac_db: Optional[str],
    tx_db: Optional[str],
    tx_taxid: Optional[int],
    gc: str,
    tm: str,
    complexity: float,
    hairpin: float,
    dimer: float,
) -> None:
    """
    Apply multi-stage filtering to SNPs.
    
    \b
    Filtering stages:
      1. Background (--bg_db): Remove sequences matching background organisms
      2. Accessibility (--ac_db): Remove sequences in repetitive/inaccessible regions
      3. Taxonomic (--tx_db): Keep only sequences assignable to target taxon
      4. Biophysical: Filter by GC, Tm, complexity, hairpin, dimer scores
    
    \b
    Output files:
      {output}.snps.tsv     - Filtered SNPs with biophysical tags
      {output}.filter.log   - Filtering statistics
    """
    from eprobe.popgen.filter import run_filter
    
    verbose = ctx.obj.get("verbose", False)
    
    # Parse GC and Tm ranges
    gc_min, gc_max = map(float, gc.split(","))
    tm_min, tm_max = map(float, tm.split(","))
    
    echo_info(f"Filtering SNPs from {input}")
    
    result = run_filter(
        input_path=input,
        reference_path=reference,
        output_prefix=output,
        threads=threads,
        probe_length=length,
        bg_databases=bg_db.split(",") if bg_db else None,
        ac_databases=ac_db.split(",") if ac_db else None,
        tx_database=tx_db,
        tx_taxid=tx_taxid,
        gc_range=(gc_min, gc_max),
        tm_range=(tm_min, tm_max),
        max_complexity=complexity,
        max_hairpin=hairpin,
        max_dimer=dimer,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Filtering failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Retained {stats['passed']} of {stats['total']} SNPs ({stats['pass_rate']:.1f}%)")


@popgen.command()
@click.option(
    "-i", "--input",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input filtered SNPs TSV file.",
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
    "-t", "--threads",
    default=1,
    type=int,
    help="Number of threads (default: 1).",
)
@click.option(
    "--window_size",
    default=10000,
    type=int,
    help="Window size for SNP selection (default: 10000).",
)
@click.option(
    "--strategy",
    type=click.Choice(["random", "weighted", "priority"]),
    default="random",
    help="Selection strategy within windows (default: random).",
)
@click.option(
    "--weights",
    type=str,
    help="Weights for biophysical features (gc,tm,complexity,hairpin,dimer).",
)
@click.option(
    "--priority_bed",
    type=click.Path(exists=True, path_type=Path),
    help="BED file of priority regions (e.g., exons).",
)
@click.option(
    "--target_count",
    type=int,
    help="Target number of SNPs (random subsample if exceeded).",
)
@click.option(
    "--seed",
    default=42,
    type=int,
    help="Random seed for reproducibility (default: 42).",
)
@click.pass_context
def select(
    ctx: click.Context,
    input: Path,
    reference: Path,
    output: Path,
    threads: int,
    window_size: int,
    strategy: str,
    weights: Optional[str],
    priority_bed: Optional[Path],
    target_count: Optional[int],
    seed: int,
) -> None:
    """
    Select optimal SNPs using window-based sampling.
    
    Ensures even distribution of probes across the genome by selecting
    one SNP per non-overlapping window.
    
    \b
    Selection strategies:
      random   - Random selection within each window
      weighted - Select SNP with best biophysical score (use --weights)
      priority - Prioritize SNPs in BED regions (use --priority_bed)
    
    \b
    Output files:
      {output}.snps.tsv     - Selected SNPs
      {output}.select.log   - Selection statistics
    """
    from eprobe.popgen.select import run_select
    
    verbose = ctx.obj.get("verbose", False)
    
    # Parse weights if provided
    weight_values = None
    if weights:
        weight_values = list(map(float, weights.split(",")))
    
    echo_info(f"Selecting SNPs from {input}")
    echo_info(f"Window size: {window_size}bp, Strategy: {strategy}")
    
    result = run_select(
        input_path=input,
        reference_path=reference,
        output_prefix=output,
        threads=threads,
        window_size=window_size,
        strategy=strategy,
        weights=weight_values,
        priority_bed=priority_bed,
        target_count=target_count,
        seed=seed,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Selection failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Selected {stats['selected']} SNPs from {stats['windows']} windows")


@popgen.command()
@click.option(
    "-i", "--input",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input selected SNPs TSV file.",
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
    default=81,
    type=int,
    help="Probe length (default: 81, odd number recommended).",
)
@click.option(
    "-s", "--snp_position",
    type=click.Choice(["center", "left", "right"]),
    default="center",
    help="SNP position within probe (default: center).",
)
@click.option(
    "--offset",
    default=0,
    type=int,
    help="Position offset from center (-/+ for left/right shift).",
)
@click.option(
    "--replace_snp/--no_replace_snp",
    default=False,
    help="Replace SNP base with non-ref/alt base (default: no).",
)
@click.option(
    "--mutation_type",
    type=click.Choice(["ts", "tv", "both"]),
    default="both",
    help="Filter by mutation type (default: both).",
)
@click.pass_context
def build(
    ctx: click.Context,
    input: Path,
    reference: Path,
    output: Path,
    length: int,
    snp_position: str,
    offset: int,
    replace_snp: bool,
    mutation_type: str,
) -> None:
    """
    Generate probe sequences from selected SNPs.
    
    Extracts genomic sequences around each SNP to create probe sequences.
    
    \b
    SNP position options:
      center - SNP at probe center (recommended for hybridization)
      left   - SNP shifted toward 5' end
      right  - SNP shifted toward 3' end
    
    \b
    Replace SNP option:
      When enabled, replaces the SNP base with a third base (not REF or ALT).
      This can reduce allele bias in capture efficiency.
    
    \b
    Output files:
      {output}.probes.fa    - Probe sequences in FASTA format
      {output}.probes.tsv   - Probe metadata
      {output}.build.log    - Build statistics
    """
    from eprobe.popgen.build import run_build
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Building probes from {input}")
    echo_info(f"Length: {length}bp, SNP position: {snp_position}")
    
    result = run_build(
        input_path=input,
        reference_path=reference,
        output_prefix=output,
        probe_length=length,
        snp_position=snp_position,
        offset=offset,
        replace_snp=replace_snp,
        mutation_type=mutation_type if mutation_type != "both" else None,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Build failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Generated {stats['probe_count']} probes")
    echo_info(f"Output: {output}.probes.fa")


@popgen.command()
@click.option(
    "-i", "--input",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input SNPs TSV file (with biophysical tags).",
)
@click.option(
    "-r", "--reference",
    type=click.Path(exists=True, path_type=Path),
    help="Reference genome FASTA (required for --mode tags).",
)
@click.option(
    "-v", "--vcf",
    type=click.Path(exists=True, path_type=Path),
    help="Original VCF file (required for --mode distance/missing).",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix.",
)
@click.option(
    "--mode",
    type=click.Choice(["tags", "distance", "missing", "all"]),
    default="tags",
    help="Assessment mode (default: tags).",
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
@click.pass_context
def assess(
    ctx: click.Context,
    input: Path,
    reference: Optional[Path],
    vcf: Optional[Path],
    output: Path,
    mode: str,
    tags: str,
    plot: bool,
) -> None:
    """
    Assess quality of probe set.
    
    \b
    Assessment modes:
      tags     - Plot biophysical property distributions
      distance - Compare IBS distance matrices (probe vs genome-wide)
      missing  - Simulate missing data tolerance
      all      - Run all assessments
    
    \b
    Output files (depends on mode):
      {output}.tags.png         - Biophysical distributions
      {output}.distance.png     - Distance matrix heatmaps
      {output}.missing.png      - Missing rate PCA plots
      {output}.assess.log       - Assessment statistics
    """
    from eprobe.popgen.assess import run_assess
    
    verbose = ctx.obj.get("verbose", False)
    
    # Validate required inputs for each mode
    if mode in ["tags", "all"] and reference is None:
        echo_error("--reference is required for tags assessment")
        raise SystemExit(1)
    
    if mode in ["distance", "missing", "all"] and vcf is None:
        echo_error("--vcf is required for distance/missing assessment")
        raise SystemExit(1)
    
    echo_info(f"Assessing probe set: {input}")
    echo_info(f"Mode: {mode}")
    
    result = run_assess(
        input_path=input,
        reference_path=reference,
        vcf_path=vcf,
        output_prefix=output,
        mode=mode,
        tags=tags.split(","),
        generate_plots=plot,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Assessment failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success("Assessment complete")
    
    if plot:
        echo_info(f"Plots saved to {output}.*")
