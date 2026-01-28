"""
Utility CLI commands.

General-purpose tools for probe manipulation:
  tiling   - Generate multi-position probes
  merge    - Combine multiple probe sets
  dedup    - Remove duplicate probes
  subset   - Filter/sample probes
  convert  - Format conversion
  validate - File validation
"""

import click
from pathlib import Path
from typing import Optional, List

from eprobe.cli.utils import (
    echo_success,
    echo_error,
    echo_info,
    echo_warning,
)
from eprobe.cli.main import AliasedGroup


@click.group(cls=AliasedGroup)
@click.pass_context
def util(ctx: click.Context) -> None:
    """
    Utility tools for probe manipulation.
    
    General-purpose commands for working with probe sets,
    applicable to both POPGEN and FUNCGEN panels.
    
    \b
    Available commands:
      tiling   - Generate probes at multiple positions per SNP
      merge    - Combine multiple probe FASTA files
      dedup    - Remove duplicate sequences
      subset   - Filter or sample probes
      convert  - Convert between formats
      validate - Validate file formats
    """
    pass


@util.command()
@click.option(
    "-i", "--input",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input SNPs TSV or FASTA file.",
)
@click.option(
    "-r", "--reference",
    type=click.Path(exists=True, path_type=Path),
    help="Reference genome FASTA (required for SNP input).",
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
    help="Probe length (default: 81).",
)
@click.option(
    "--positions",
    type=str,
    default="center,left,right",
    help="Positions: 'center,left,right' or offsets like '-20,0,+20' (default: center,left,right).",
)
@click.option(
    "--offset",
    default=20,
    type=int,
    help="Offset distance for left/right positions (default: 20).",
)
@click.option(
    "--suffix_style",
    type=click.Choice(["position", "index"]),
    default="position",
    help="Suffix style: 'position' (_L20, _C, _R20) or 'index' (_1, _2, _3).",
)
@click.option(
    "--step",
    type=int,
    help="Step size for FASTA tiling (sliding window mode).",
)
@click.pass_context
def tiling(
    ctx: click.Context,
    input: Path,
    reference: Optional[Path],
    output: Path,
    length: int,
    positions: str,
    offset: int,
    suffix_style: str,
    step: Optional[int],
) -> None:
    """
    Generate probes at multiple positions.
    
    For SNP-based input, creates multiple probes per SNP with the SNP
    at different positions (center, shifted left, shifted right).
    
    For FASTA input, uses sliding window tiling with specified step.
    
    \b
    Position syntax:
      Named:  center,left,right (uses --offset for shift distance)
      Numeric: -30,-15,0,+15,+30 (explicit offsets from center)
    
    \b
    Use cases:
      - Triple tiling for robust SNP capture
      - Dense coverage of target regions
      - Reducing position-specific capture bias
    
    \b
    Example (triple tiling):
      eprobe util tiling -i snps.tsv -r ref.fa -o tiled --positions center,left,right --offset 25
      # Generates 3 probes per SNP: _C (center), _L25 (left), _R25 (right)
    """
    from eprobe.util.tiling import run_tiling
    
    verbose = ctx.obj.get("verbose", False)
    
    # Determine input type
    is_fasta = input.suffix.lower() in [".fa", ".fasta", ".fna"]
    
    if not is_fasta and reference is None:
        echo_error("--reference is required for SNP TSV input")
        raise SystemExit(1)
    
    echo_info(f"Generating tiled probes from {input}")
    echo_info(f"Positions: {positions}, Length: {length}bp")
    
    result = run_tiling(
        input_path=input,
        reference_path=reference,
        output_prefix=output,
        probe_length=length,
        positions=positions,
        offset=offset,
        suffix_style=suffix_style,
        step=step,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Tiling failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Generated {stats['probe_count']} probes from {stats['input_count']} inputs")
    echo_info(f"Tiling factor: {stats['tiling_factor']}x")


@util.command()
@click.option(
    "-i", "--input",
    required=True,
    multiple=True,
    type=str,
    help="Input FASTA files (can specify multiple: -i file1.fa -i file2.fa or prefix:file.fa).",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output FASTA file.",
)
@click.option(
    "--tag_source/--no_tag_source",
    default=False,
    help="Add source prefix to sequence IDs (default: no).",
)
@click.option(
    "--format",
    type=click.Choice(["fasta", "tsv"]),
    default="fasta",
    help="Input/output format (default: fasta).",
)
@click.pass_context
def merge(
    ctx: click.Context,
    input: tuple,
    output: Path,
    tag_source: bool,
    format: str,
) -> None:
    """
    Merge multiple probe sets.
    
    Combines multiple FASTA or TSV files into a single file.
    Optionally adds source tags to distinguish origins.
    
    \b
    Input syntax:
      Simple:    -i file1.fa -i file2.fa
      With tag:  -i popgen:popgen.fa -i funcgen:funcgen.fa
    
    \b
    Output with --tag_source:
      >popgen_SNP001
      >funcgen_GENE001
    
    \b
    Note: Does not automatically deduplicate. Use 'eprobe util dedup' after merging.
    """
    from eprobe.util.merge import run_merge
    
    verbose = ctx.obj.get("verbose", False)
    
    # Parse input files (handle prefix:path syntax)
    input_files: List[tuple[str, Path]] = []
    for inp in input:
        if ":" in inp and not inp.startswith("/"):
            parts = inp.split(":", 1)
            tag = parts[0]
            path = Path(parts[1])
        else:
            path = Path(inp)
            tag = path.stem
        
        if not path.exists():
            echo_error(f"Input file not found: {path}")
            raise SystemExit(1)
        
        input_files.append((tag, path))
    
    echo_info(f"Merging {len(input_files)} probe sets")
    
    result = run_merge(
        input_files=input_files,
        output_path=output,
        tag_source=tag_source,
        file_format=format,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Merge failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Merged {stats['total_probes']} probes into {output}")


@util.command()
@click.option(
    "-i", "--input",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input FASTA file.",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output FASTA file.",
)
@click.option(
    "--method",
    type=click.Choice(["simple", "cdhit"]),
    default="simple",
    help="Deduplication method (default: simple).",
)
@click.option(
    "--identity",
    default=0.95,
    type=float,
    help="Sequence identity threshold for cd-hit (default: 0.95).",
)
@click.option(
    "--coverage",
    default=0.9,
    type=float,
    help="Coverage threshold for cd-hit (default: 0.9).",
)
@click.option(
    "--keep",
    type=click.Choice(["first", "longest", "random"]),
    default="first",
    help="Which duplicate to keep (default: first).",
)
@click.option(
    "--output_clusters",
    type=click.Path(path_type=Path),
    help="Output cluster information file (cd-hit mode).",
)
@click.option(
    "-t", "--threads",
    default=1,
    type=int,
    help="Number of threads for cd-hit (default: 1).",
)
@click.pass_context
def dedup(
    ctx: click.Context,
    input: Path,
    output: Path,
    method: str,
    identity: float,
    coverage: float,
    keep: str,
    output_clusters: Optional[Path],
    threads: int,
) -> None:
    """
    Remove duplicate probe sequences.
    
    \b
    Methods:
      simple - Exact sequence match (fast, 100% identity only)
      cdhit  - Clustering by similarity (flexible threshold)
    
    \b
    Simple mode:
      Removes exact duplicates, keeping first/longest/random occurrence.
      Fast and memory-efficient.
    
    \b
    CD-HIT mode:
      Clusters similar sequences and keeps representative.
      Useful for reducing redundancy while maintaining coverage.
      Requires cd-hit to be installed.
    
    \b
    Example (remove 95% similar sequences):
      eprobe util dedup -i probes.fa -o dedup.fa --method cdhit --identity 0.95
    """
    from eprobe.util.dedup import run_dedup
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Deduplicating {input}")
    echo_info(f"Method: {method}")
    
    if method == "cdhit":
        echo_info(f"Identity threshold: {identity}")
    
    result = run_dedup(
        input_path=input,
        output_path=output,
        method=method,
        identity=identity,
        coverage=coverage,
        keep=keep,
        clusters_path=output_clusters,
        threads=threads,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Deduplication failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    removed = stats['input_count'] - stats['output_count']
    echo_success(f"Removed {removed} duplicates ({stats['output_count']} remaining)")


@util.command()
@click.option(
    "-i", "--input",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input file (FASTA or TSV).",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output file.",
)
@click.option(
    "--filter",
    "filter_expr",
    type=str,
    help="Filter expression (e.g., 'gc >= 40 and gc <= 60').",
)
@click.option(
    "--sample",
    type=int,
    help="Random sample to this many sequences.",
)
@click.option(
    "--chrom",
    type=str,
    help="Keep only these chromosomes (comma-separated).",
)
@click.option(
    "--bed",
    type=click.Path(exists=True, path_type=Path),
    help="Keep only sequences in BED regions.",
)
@click.option(
    "-r", "--reference",
    type=click.Path(exists=True, path_type=Path),
    help="Reference genome (required for --strategy uniform).",
)
@click.option(
    "--strategy",
    type=click.Choice(["random", "uniform"]),
    default="random",
    help="Sampling strategy (default: random).",
)
@click.option(
    "--window_size",
    default=100000,
    type=int,
    help="Window size for uniform sampling (default: 100000).",
)
@click.option(
    "--seed",
    default=42,
    type=int,
    help="Random seed (default: 42).",
)
@click.pass_context
def subset(
    ctx: click.Context,
    input: Path,
    output: Path,
    filter_expr: Optional[str],
    sample: Optional[int],
    chrom: Optional[str],
    bed: Optional[Path],
    reference: Optional[Path],
    strategy: str,
    window_size: int,
    seed: int,
) -> None:
    """
    Filter or sample probes.
    
    \b
    Filter options:
      --filter  - Expression filter on TSV columns
      --chrom   - Keep specific chromosomes
      --bed     - Keep sequences in BED regions
    
    \b
    Sample options:
      --sample    - Target number of sequences
      --strategy  - random or uniform (maintains distribution)
    
    \b
    Filter expression syntax:
      'gc >= 40 and gc <= 60'
      'tm > 70 and complexity < 2'
      'chrom == "chr1"'
    
    \b
    Uniform sampling:
      Maintains even genomic distribution by sampling within windows.
      Useful for reducing probe count while keeping coverage even.
    """
    from eprobe.util.subset import run_subset
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Subsetting {input}")
    
    if filter_expr:
        echo_info(f"Filter: {filter_expr}")
    if sample:
        echo_info(f"Target count: {sample} ({strategy} sampling)")
    
    result = run_subset(
        input_path=input,
        output_path=output,
        filter_expr=filter_expr,
        sample_count=sample,
        chromosomes=chrom.split(",") if chrom else None,
        bed_path=bed,
        reference_path=reference,
        strategy=strategy,
        window_size=window_size,
        seed=seed,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Subsetting failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Retained {stats['output_count']} of {stats['input_count']} sequences")


@util.command()
@click.option(
    "-i", "--input",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input file.",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output file.",
)
@click.option(
    "--from",
    "from_fmt",
    required=True,
    type=click.Choice(["vcf", "snp_tsv", "fasta", "bed"]),
    help="Input format.",
)
@click.option(
    "--to",
    "to_fmt",
    required=True,
    type=click.Choice(["snp_tsv", "fasta", "bed", "probe_tsv"]),
    help="Output format.",
)
@click.option(
    "--include_stats/--no_include_stats",
    default=False,
    help="Include sequence statistics in output (for probe_tsv).",
)
@click.pass_context
def convert(
    ctx: click.Context,
    input: Path,
    output: Path,
    from_fmt: str,
    to_fmt: str,
    include_stats: bool,
) -> None:
    """
    Convert between file formats.
    
    \b
    Supported conversions:
      vcf → snp_tsv      Extract SNPs from VCF
      snp_tsv → bed      Convert SNP positions to BED
      fasta → probe_tsv  Extract probe info with stats
      bed → snp_tsv      (requires additional info)
    
    \b
    Example:
      eprobe util convert -i snps.tsv -o snps.bed --from snp_tsv --to bed
    """
    from eprobe.util.convert import run_convert
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Converting {from_fmt} → {to_fmt}")
    
    result = run_convert(
        input_path=input,
        output_path=output,
        from_format=from_fmt,
        to_format=to_fmt,
        include_stats=include_stats,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Conversion failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Converted {stats['record_count']} records to {output}")


@util.command()
@click.option(
    "-i", "--input",
    type=click.Path(exists=True, path_type=Path),
    help="Input file to validate.",
)
@click.option(
    "--db",
    type=click.Path(exists=True, path_type=Path),
    help="Database path to validate.",
)
@click.option(
    "--schema",
    type=click.Choice(["snp_dataframe", "probe_fasta", "bed"]),
    help="Schema for file validation.",
)
@click.option(
    "--type",
    "db_type",
    type=click.Choice(["kraken2", "bowtie2"]),
    help="Database type for database validation.",
)
@click.pass_context
def validate(
    ctx: click.Context,
    input: Optional[Path],
    db: Optional[Path],
    schema: Optional[str],
    db_type: Optional[str],
) -> None:
    """
    Validate file formats or databases.
    
    \b
    File validation (--input --schema):
      snp_dataframe - Validate SNP TSV format
      probe_fasta   - Validate probe FASTA format
      bed           - Validate BED format
    
    \b
    Database validation (--db --type):
      kraken2 - Validate Kraken2 database
      bowtie2 - Validate Bowtie2 index
    
    \b
    Example:
      eprobe util validate -i snps.tsv --schema snp_dataframe
      eprobe util validate --db /path/to/kraken_db --type kraken2
    """
    from eprobe.util.validate import run_validate
    
    verbose = ctx.obj.get("verbose", False)
    
    if input and schema:
        echo_info(f"Validating {input} as {schema}")
        result = run_validate(
            file_path=input,
            schema=schema,
            verbose=verbose,
        )
    elif db and db_type:
        echo_info(f"Validating {db_type} database at {db}")
        result = run_validate(
            database_path=db,
            database_type=db_type,
            verbose=verbose,
        )
    else:
        echo_error("Specify either (--input --schema) or (--db --type)")
        raise SystemExit(1)
    
    if result.is_err():
        echo_error(f"Validation failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    validation = result.unwrap()
    
    if validation["valid"]:
        echo_success("Validation passed")
        for msg in validation.get("messages", []):
            echo_info(f"  {msg}")
    else:
        echo_error("Validation failed")
        for error in validation.get("errors", []):
            echo_error(f"  {error}")
        raise SystemExit(1)
