"""
Utility CLI commands.

Post-design probe manipulation toolkit:
  merge    - Merge probe sets with sequence-level dedup
  tile     - Simple FASTA sliding-window tiling
  adapter  - Add 5'/3' adapter sequences
  rename   - Batch rename with prefix or ID map
  assess   - Biophysical tags assessment & filtering
  sample   - Random subsampling & CD-HIT clustering
  target   - Bowtie2 mapping → target BED generation

Legacy tools:
  tiling   - Multi-position SNP-based probe tiling
  dedup    - Exact/CD-HIT deduplication
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
    
    Post-design toolkit for working with probe sets:
    merge, tile, add adapters, rename, assess quality,
    sample/cluster, filter, and generate target regions.
    
    \b
    Main commands:
      merge    - Merge probe sets with deduplication
      tile     - Tile FASTA sequences into probes
      adapter  - Add adapter sequences
      rename   - Batch rename probe IDs
      assess   - Biophysical quality assessment/filtering
      sample   - Random or cluster-based sampling
      target   - Map probes → target BED regions
      filter   - Apply popgen filter stages to FASTA
    
    \b
    Legacy commands:
      tiling   - Multi-position SNP probe tiling
      dedup    - Remove duplicate sequences
      subset   - Filter/sample probes
      convert  - Format conversion
      validate - File validation
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
    type=click.Path(exists=True, path_type=Path),
    help="Input FASTA files (specify multiple: -i file1.fa -i file2.fa).",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix.",
)
@click.option(
    "--keep",
    type=click.Choice(["first"]),
    default="first",
    help="Which duplicate to keep (default: first).",
)
@click.pass_context
def merge(
    ctx: click.Context,
    input: tuple,
    output: Path,
    keep: str,
) -> None:
    """
    Merge probe sets with sequence-level deduplication.
    
    Combines multiple FASTA files, detects exact duplicate sequences
    (case-insensitive), keeps one representative per unique sequence,
    and reports all removed duplicates.
    
    \b
    Outputs:
      {output}.merged.fa              - Unique probes
      {output}.removed_duplicates.tsv - Removed/kept ID pairs
      {output}.merge_summary.txt      - Statistics
    
    \b
    Example:
      eprobe util merge -i popgen.fa -i funcgen.fa -o combined
    """
    from eprobe.util.merge import run_merge
    
    verbose = ctx.obj.get("verbose", False)
    input_files = [Path(p) for p in input]
    
    echo_info(f"Merging {len(input_files)} probe sets")
    
    result = run_merge(
        input_files=input_files,
        output_prefix=output,
        keep=keep,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Merge failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(
        f"Merged {stats['total_input']} → {stats['unique_count']} unique probes "
        f"({stats['removed_count']} duplicates removed)"
    )


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


# ===================================================================
#  New commands
# ===================================================================


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
    help="Output prefix.",
)
@click.option(
    "-l", "--length",
    default=81,
    type=int,
    help="Probe length (default: 81).",
)
@click.option(
    "--step",
    default=30,
    type=int,
    help="Step size between probes (default: 30).",
)
@click.pass_context
def tile(
    ctx: click.Context,
    input: Path,
    output: Path,
    length: int,
    step: int,
) -> None:
    """
    Tile FASTA sequences into fixed-length probes.
    
    Generates probes by sliding window across input sequences.
    Useful for converting gene/region sequences into probe panels.
    
    \b
    Outputs:
      {output}.tiled.fa - Tiled probes
    
    \b
    Probe IDs: {original_id}:{start}-{end}_tile{n}
    
    \b
    Example:
      eprobe util tile -i genes.fa -o tiled -l 81 --step 30
    """
    from eprobe.util.tile import run_tile
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Tiling {input}: length={length}bp, step={step}bp")
    
    result = run_tile(
        input_path=input,
        output_prefix=output,
        probe_length=length,
        step=step,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Tiling failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(
        f"Generated {stats['probe_count']} probes from "
        f"{stats['input_count']} sequences"
    )


@util.command()
@click.option(
    "-i", "--input",
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
    "--adapter",
    "adapter_seq",
    required=True,
    type=str,
    help="Adapter sequence (DNA bases only, e.g. AGATCGGAAGAGC).",
)
@click.option(
    "--end",
    type=click.Choice(["5prime", "3prime"]),
    default="5prime",
    help="Which end to add adapter (default: 5prime).",
)
@click.pass_context
def adapter(
    ctx: click.Context,
    input: Path,
    output: Path,
    adapter_seq: str,
    end: str,
) -> None:
    """
    Add adapter sequences to probes.
    
    Batch prepend (5') or append (3') an adapter sequence to
    all probes in a FASTA file.
    
    \b
    Outputs:
      {output}.adapted.fa - Probes with adapter added
    
    \b
    Example:
      eprobe util adapter -i probes.fa -o adapted --adapter AGATCGGAAGAGC --end 5prime
    """
    from eprobe.util.adapter import run_adapter
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Adding {end} adapter: {adapter_seq}")
    
    result = run_adapter(
        input_path=input,
        output_prefix=output,
        adapter_seq=adapter_seq,
        end=end,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Adapter addition failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(
        f"Added {stats['adapter_length']}bp adapter to {stats['probe_count']} probes "
        f"({stats['original_length']}bp → {stats['new_length']}bp)"
    )


@util.command()
@click.option(
    "-i", "--input",
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
    "--prefix",
    type=str,
    help="Rename with prefix: {prefix}_0001, _0002, ...",
)
@click.option(
    "--id_map",
    type=click.Path(exists=True, path_type=Path),
    help="TSV mapping file (old_id<TAB>new_id).",
)
@click.option(
    "--start_index",
    default=1,
    type=int,
    help="Starting index for prefix mode (default: 1).",
)
@click.pass_context
def rename(
    ctx: click.Context,
    input: Path,
    output: Path,
    prefix: Optional[str],
    id_map: Optional[Path],
    start_index: int,
) -> None:
    """
    Batch rename probe IDs.
    
    \b
    Mode 1 (--prefix):
      Rename all IDs to {prefix}_0001, {prefix}_0002, ...
      Auto zero-padded to fit total count.
    
    \b
    Mode 2 (--id_map):
      Rename using TSV mapping file (old_id<TAB>new_id).
      Unmapped IDs keep their original names.
    
    \b
    Outputs:
      {output}.renamed.fa  - Renamed probes
      {output}.id_map.tsv  - Correspondence map (new_id → old_id)
    
    \b
    Example:
      eprobe util rename -i probes.fa -o renamed --prefix PANEL_A
      eprobe util rename -i probes.fa -o renamed --id_map mapping.tsv
    """
    from eprobe.util.rename import run_rename
    
    verbose = ctx.obj.get("verbose", False)
    
    if prefix is None and id_map is None:
        echo_error("Either --prefix or --id_map must be provided")
        raise SystemExit(1)
    
    mode_str = f"prefix='{prefix}'" if prefix else f"id_map={id_map}"
    echo_info(f"Renaming probes ({mode_str})")
    
    result = run_rename(
        input_path=input,
        output_prefix=output,
        prefix=prefix,
        id_map_path=id_map,
        start_index=start_index,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Rename failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    msg = f"Renamed {stats['renamed_count']} probes"
    if stats.get("unmapped_count", 0) > 0:
        msg += f" ({stats['unmapped_count']} unmapped, kept original)"
    echo_success(msg)


@util.command()
@click.option(
    "-i", "--input",
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
    "-m", "--mode",
    type=click.Choice(["tags", "filter"]),
    default="tags",
    help="Mode: 'tags' (compute metrics) or 'filter' (apply thresholds).",
)
@click.option(
    "--tag",
    type=str,
    default="gc,tm,complexity,hairpin",
    help="Tags to compute (comma-separated: gc,tm,complexity,hairpin,dimer).",
)
@click.option(
    "--gc_min", default=35.0, type=float, help="Min GC%% (filter mode, default: 35). Set -1 to disable GC.",
)
@click.option(
    "--gc_max", default=65.0, type=float, help="Max GC%% (filter mode, default: 65). Set -1 to disable GC.",
)
@click.option(
    "--tm_min", default=55.0, type=float, help="Min Tm (filter mode, default: 55). Set -1 to disable Tm.",
)
@click.option(
    "--tm_max", default=75.0, type=float, help="Max Tm (filter mode, default: 75). Set -1 to disable Tm.",
)
@click.option(
    "--complexity_max", default=2.0, type=float,
    help="Max DUST complexity (filter mode, default: 2.0). Set -1 to disable.",
)
@click.option(
    "--hairpin", default=0.95, type=float,
    help="Hairpin threshold (filter mode, default: 0.95 = top 5%% removed). "
         ">=1: absolute; <1: percentile (e.g. 0.95 keeps 95%%). Set -1 to disable.",
)
@click.option(
    "--dimer", default=0.15, type=float,
    help="Smart dimer filter sensitivity (default: 0.15). "
         "Graph-based: identifies cross-hybridizing groups, keeps one per group. "
         "Value = min fraction of shared k-mers for dimer edge. Set -1 to disable.",
)
@click.option(
    "--no_plots",
    is_flag=True,
    default=False,
    help="Skip generating distribution plots (tags mode).",
)
@click.pass_context
def assess(
    ctx: click.Context,
    input: Path,
    output: Path,
    mode: str,
    tag: str,
    gc_min: float,
    gc_max: float,
    tm_min: float,
    tm_max: float,
    complexity_max: float,
    hairpin: float,
    dimer: float,
    no_plots: bool,
) -> None:
    """
    Biophysical assessment or filtering of probes.
    
    \b
    Tags mode (default):
      Computes biophysical metrics for each probe → TSV + plots.
      Available tags: gc, tm, complexity, hairpin, dimer.
    
    \b
    Filter mode:
      3-stage biophysical filter (same as popgen filter):
        Stage 1: GC, Tm, DUST (absolute thresholds)
        Stage 2: Hairpin (percentile)
        Stage 3: Dimer (percentile)
    
    \b
    Output (tags mode):
      {output}.tags_stats.tsv     - Per-probe metrics
      {output}.tags_summary.txt   - Summary statistics
      {output}.{tag}_dist.jpg     - Distribution plots
    
    \b
    Output (filter mode):
      {output}.filtered.fa        - Passed probes
      {output}.filtered_out.fa    - Rejected probes
    
    \b
    Example:
      eprobe util assess -i probes.fa -o qc -m tags --tag gc,tm,complexity
      eprobe util assess -i probes.fa -o qc -m filter --gc_min 40 --gc_max 60
    """
    from eprobe.util.assess import run_assess
    
    verbose = ctx.obj.get("verbose", False)
    tags_list = [t.strip() for t in tag.split(",")]
    
    if mode == "tags":
        echo_info(f"Assessing {input} — tags: {', '.join(tags_list)}")
    else:
        echo_info(f"Filtering {input} — GC={gc_min}-{gc_max}%, Tm={tm_min}-{tm_max}°C")
    
    result = run_assess(
        input_path=input,
        output_prefix=output,
        mode=mode,
        tags=tags_list,
        gc_min=gc_min,
        gc_max=gc_max,
        tm_min=tm_min,
        tm_max=tm_max,
        complexity_max=complexity_max,
        hairpin=hairpin,
        dimer=dimer,
        generate_plots=not no_plots,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Assess failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    if mode == "tags":
        echo_success(f"Assessed {stats['probe_count']} probes → {stats['stats_file']}")
        for t, s in stats.get("summary", {}).items():
            echo_info(f"  {t}: mean={s['mean']}, median={s['median']}, "
                       f"range=[{s['min']}, {s['max']}]")
    else:
        echo_success(
            f"Filtered: {stats['passed_count']}/{stats['input_count']} passed "
            f"({stats['failed_count']} removed)"
        )


@util.command()
@click.option(
    "-i", "--input",
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
    "-m", "--mode",
    type=click.Choice(["random", "cluster"]),
    default="random",
    help="Sampling mode (default: random).",
)
@click.option(
    "-n", "--n",
    type=int,
    help="Target number of probes (required for random, optional for cluster).",
)
@click.option(
    "--seed",
    default=42,
    type=int,
    help="Random seed (default: 42).",
)
@click.option(
    "--identity",
    default=0.95,
    type=float,
    help="CD-HIT identity threshold (cluster mode, default: 0.95).",
)
@click.option(
    "--coverage",
    default=0.9,
    type=float,
    help="CD-HIT coverage threshold (cluster mode, default: 0.9).",
)
@click.option(
    "--min_members",
    type=int,
    help="Keep families with ≥ N members (cluster mode).",
)
@click.option(
    "--max_members",
    type=int,
    help="Keep families with ≤ N members (cluster mode).",
)
@click.option(
    "--min_proportion",
    type=float,
    help="Keep families with ≥ proportion (0–1, cluster mode).",
)
@click.option(
    "--max_proportion",
    type=float,
    help="Keep families with ≤ proportion (0–1, cluster mode).",
)
@click.option(
    "-t", "--threads",
    default=1,
    type=int,
    help="Number of threads (cluster mode, default: 1).",
)
@click.pass_context
def sample(
    ctx: click.Context,
    input: Path,
    output: Path,
    mode: str,
    n: Optional[int],
    seed: int,
    identity: float,
    coverage: float,
    min_members: Optional[int],
    max_members: Optional[int],
    min_proportion: Optional[float],
    max_proportion: Optional[float],
    threads: int,
) -> None:
    """
    Sample probes: random or cluster-based.
    
    \b
    Random mode:
      Randomly sample --n probes from input.
    
    \b
    Cluster mode:
      1. Run cd-hit-est clustering
      2. Report families (cluster report TSV)
      3. Select one representative per family
      4. Optionally filter families by member count/proportion
      5. Optionally further downsample with --n
    
    \b
    Outputs:
      {output}.sampled.fa          - Selected probes
      {output}.cluster_report.tsv  - Family details (cluster mode)
      {output}.sample_summary.txt  - Statistics (cluster mode)
    
    \b
    Example:
      eprobe util sample -i probes.fa -o sub -m random -n 1000
      eprobe util sample -i probes.fa -o sub -m cluster --identity 0.90
      eprobe util sample -i probes.fa -o sub -m cluster --max_members 5
    """
    from eprobe.util.sample import run_sample
    
    verbose = ctx.obj.get("verbose", False)
    
    if mode == "random":
        if n is None:
            echo_error("--n is required for random sampling")
            raise SystemExit(1)
        echo_info(f"Random sampling {n} probes from {input}")
    else:
        echo_info(f"Cluster-based sampling: identity={identity}")
    
    result = run_sample(
        input_path=input,
        output_prefix=output,
        mode=mode,
        n=n,
        seed=seed,
        identity=identity,
        coverage=coverage,
        min_members=min_members,
        max_members=max_members,
        min_proportion=min_proportion,
        max_proportion=max_proportion,
        threads=threads,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Sampling failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    if mode == "random":
        echo_success(
            f"Sampled {stats['output_count']}/{stats['input_count']} probes"
        )
    else:
        echo_success(
            f"Selected {stats['output_count']} representatives from "
            f"{stats['total_families']} families "
            f"({stats['input_count']} input probes)"
        )


@util.command()
@click.option(
    "-i", "--input",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input probe FASTA file.",
)
@click.option(
    "-g", "--genome",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Reference genome FASTA.",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix.",
)
@click.option(
    "--index",
    "index_prefix",
    type=click.Path(path_type=Path),
    help="Pre-built bowtie2 index prefix (skips index building).",
)
@click.option(
    "-t", "--threads",
    default=1,
    type=int,
    help="Number of threads (default: 1).",
)
@click.option(
    "--bt2_mode",
    type=click.Choice([
        "very-sensitive-local", "sensitive-local",
        "very-sensitive", "sensitive",
    ]),
    default="very-sensitive-local",
    help="Bowtie2 alignment preset (default: very-sensitive-local).",
)
@click.option(
    "--merge_distance",
    default=0,
    type=int,
    help="Max distance for merging adjacent regions (default: 0).",
)
@click.option(
    "--min_score",
    type=int,
    help="Minimum alignment score.",
)
@click.pass_context
def target(
    ctx: click.Context,
    input: Path,
    genome: Path,
    output: Path,
    index_prefix: Optional[Path],
    threads: int,
    bt2_mode: str,
    merge_distance: int,
    min_score: Optional[int],
) -> None:
    """
    Generate target BED regions from probe-to-genome mapping.
    
    Maps probes to a genome with bowtie2, merges overlapping
    alignment regions, and outputs non-redundant target regions.
    
    \b
    Use case:
      After switching reference genomes, define capture target
      regions based on where existing probes hybridize.
    
    \b
    Pipeline:
      1. Build bowtie2 index (or use --index)
      2. Map probes (local alignment)
      3. BAM → BED → merge overlaps
      4. Output target BED
    
    \b
    Outputs:
      {output}.target.bed          - Target regions
      {output}.target_summary.txt  - Mapping statistics
    
    \b
    Example:
      eprobe util target -i probes.fa -g genome.fa -o targets -t 8
      eprobe util target -i probes.fa -g genome.fa -o targets --index /path/to/bt2_idx
    """
    from eprobe.util.target import run_target
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Mapping {input} → {genome}")
    echo_info(f"Mode: {bt2_mode}, threads: {threads}")
    
    result = run_target(
        probe_fasta=input,
        genome_path=genome,
        output_prefix=output,
        index_prefix=index_prefix,
        threads=threads,
        mode=bt2_mode,
        merge_distance=merge_distance,
        min_score=min_score,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Target generation failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(
        f"Generated {stats['n_regions']} target regions "
        f"({stats['total_bp']:,}bp) from {stats['n_mapped']}/{stats['n_probes']} "
        f"mapped probes ({stats['mapping_rate']}%)"
    )


# =====================================================================
# filter — FASTA-level popgen filter stages
# =====================================================================

@util.command("filter")
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
    help="Output prefix.",
)
@click.option(
    "-t", "--threads",
    default=1,
    type=int,
    help="Threads for external tools (default: 1).",
)
@click.option(
    "--step",
    "steps",
    required=True,
    multiple=True,
    type=click.Choice(["bg", "ac", "tx", "biophysical"], case_sensitive=False),
    help="Filter step(s) to apply (can specify multiple, applied in order).",
)
# --- BG options ---
@click.option(
    "--bg_db",
    type=str,
    help="Kraken2 database path(s) for background filtering (comma-separated).",
)
# --- AC options ---
@click.option(
    "--ac_db",
    type=str,
    help="Bowtie2 index path(s) for accessibility filtering (comma-separated).",
)
@click.option(
    "--ac_mode",
    type=click.Choice(["strict", "relaxed"]),
    default="strict",
    help="Accessibility mode: 'strict' (unique only) or 'relaxed' (default: strict).",
)
@click.option(
    "--ac_score_diff",
    type=int,
    default=10,
    help="Min AS-XS score difference for relaxed mode (default: 10).",
)
@click.option(
    "--ac_min_genomes",
    type=int,
    help="Min genomes a probe must map to (default: all).",
)
# --- TX options ---
@click.option(
    "--tx_db",
    type=str,
    help="Bowtie2 index path(s) for taxonomic filtering (comma-separated).",
)
@click.option(
    "--tx_taxid",
    type=str,
    help="Target taxon ID(s) (comma-separated).",
)
@click.option(
    "--tx_target_names",
    type=str,
    help="Target taxon name(s) (comma-separated, requires --tx_names).",
)
@click.option(
    "--tx_names",
    type=click.Path(exists=True, path_type=Path),
    help="NCBI names.dmp file.",
)
@click.option(
    "--tx_nodes",
    type=click.Path(exists=True, path_type=Path),
    help="NCBI nodes.dmp file.",
)
@click.option(
    "--tx_acc2tax",
    type=click.Path(exists=True, path_type=Path),
    help="Accession-to-taxid mapping file.",
)
@click.option(
    "--tx_minedit",
    type=int,
    default=0,
    help="Min edit distance for ngsLCA (default: 0).",
)
@click.option(
    "--tx_maxedit",
    type=int,
    default=2,
    help="Max edit distance for ngsLCA (default: 2).",
)
@click.option(
    "--tx_keep_hits",
    type=int,
    default=100,
    help="Alignments to report per sequence (default: 100).",
)
@click.option(
    "--tx_mode",
    type=click.Choice(["lca", "besthit"]),
    default="lca",
    help="TX classification mode: 'lca' or 'besthit' (default: lca).",
)
@click.option(
    "--tx_outgroup_ids",
    type=str,
    help="Outgroup taxon ID(s) for best-hit mode (comma-separated).",
)
@click.option(
    "--tx_outgroup_names",
    type=str,
    help="Outgroup taxon name(s) for best-hit mode (comma-separated).",
)
@click.option(
    "--tx_exclude_mode",
    type=click.Choice(["sim", "nm", "diff"]),
    default="diff",
    help="Best-hit exclusion strategy (default: diff).",
)
@click.option(
    "--tx_exclude_val",
    type=float,
    default=2,
    help="Threshold for best-hit exclusion (default: 2).",
)
@click.option(
    "--tx_target_min_sim",
    type=float,
    default=90.0,
    help="Min target similarity %% for best-hit (default: 90).",
)
@click.option(
    "--tx_target_max_nm",
    type=int,
    default=5,
    help="Max target NM for best-hit (default: 5).",
)
# --- Biophysical options ---
@click.option(
    "--gc",
    type=str,
    default="35,65",
    help="GC content range (min,max, default: 35,65). Set -1 to disable.",
)
@click.option(
    "--tm",
    type=str,
    default="65,85",
    help="Melting temperature range (min,max, default: 65,85). Set -1 to disable.",
)
@click.option(
    "--complexity",
    type=float,
    default=2.0,
    help="Max DUST complexity score (default: 2.0). Set -1 to disable.",
)
@click.option(
    "--hairpin",
    type=float,
    default=0.95,
    help="Hairpin filter threshold (default: 0.95 = top 5%% removed). "
         ">=1: absolute score threshold; <1: percentile (e.g. 0.95 keeps 95%%). "
         "Uses exponential k-mer continuity scoring. Set -1 to disable.",
)
@click.option(
    "--dimer",
    type=float,
    default=0.15,
    help="Smart dimer filter sensitivity (default: 0.15). "
         "Graph-based: identifies cross-hybridizing groups, keeps one per group. "
         "Value = min fraction of shared k-mers for dimer edge. Set -1 to disable.",
)
@click.option(
    "--nn-table",
    type=click.Choice(["DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4", "R_DNA_NN1"]),
    default="DNA_NN4",
    help="Nearest-neighbor table for Tm (default: DNA_NN4).",
)
@click.option(
    "--na-conc",
    type=float,
    default=50.0,
    help="Sodium concentration in mM (default: 50.0).",
)
@click.option(
    "--keep-temp",
    is_flag=True,
    help="Keep temporary files.",
)
@click.pass_context
def filter_cmd(
    ctx: click.Context,
    input: Path,
    output: Path,
    threads: int,
    steps: tuple,
    bg_db: Optional[str],
    ac_db: Optional[str],
    ac_mode: str,
    ac_score_diff: int,
    ac_min_genomes: Optional[int],
    tx_db: Optional[str],
    tx_taxid: Optional[str],
    tx_target_names: Optional[str],
    tx_names: Optional[Path],
    tx_nodes: Optional[Path],
    tx_acc2tax: Optional[Path],
    tx_minedit: int,
    tx_maxedit: int,
    tx_keep_hits: int,
    tx_mode: str,
    tx_outgroup_ids: Optional[str],
    tx_outgroup_names: Optional[str],
    tx_exclude_mode: str,
    tx_exclude_val: float,
    tx_target_min_sim: float,
    tx_target_max_nm: int,
    gc: str,
    tm: str,
    complexity: float,
    hairpin: float,
    dimer: float,
    nn_table: str,
    na_conc: float,
    keep_temp: bool,
) -> None:
    """
    Filter FASTA sequences through popgen filter stages.

    \b
    Applies one or more filter steps directly to a FASTA file:
      bg          - Background noise removal (Kraken2)
      ac          - Accessibility / mappability check (Bowtie2)
      tx          - Taxonomic assignment (LCA or best-hit)
      biophysical - GC, Tm, complexity, hairpin, dimer

    \b
    Examples:
      eprobe util filter -i probes.fa -o out --step biophysical
      eprobe util filter -i probes.fa -o out --step bg --bg_db kraken_db
      eprobe util filter -i probes.fa -o out --step tx --tx_db bt2_idx \\
          --tx_taxid 9606 --tx_names names.dmp --tx_nodes nodes.dmp \\
          --tx_acc2tax acc2taxid --tx_mode besthit

    \b
    Output:
      {output}.filtered.fa   - Sequences passing all steps
      {output}.rejected.fa   - Sequences removed
      {output}.filter.log    - Filtering log
    """
    from eprobe.util.filter import run_util_filter

    verbose = ctx.obj.get("verbose", False)

    # Parse ranges (support -1 to disable)
    if gc.strip() == "-1":
        gc_min, gc_max = -1.0, -1.0
    else:
        gc_min, gc_max = map(float, gc.split(","))
    
    if tm.strip() == "-1":
        tm_min, tm_max = -1.0, -1.0
    else:
        tm_min, tm_max = map(float, tm.split(","))

    # Parse target taxonomy IDs
    tx_ids = None
    if tx_taxid:
        tx_ids = [int(t.strip()) for t in tx_taxid.split(",")]

    # Convert target names → taxids
    if tx_target_names:
        if not tx_names:
            echo_error("--tx_target_names requires --tx_names")
            raise SystemExit(1)
        from eprobe.popgen.filter import parse_taxonomy_names_to_taxids
        name_list = [n.strip() for n in tx_target_names.split(",")]
        name_result = parse_taxonomy_names_to_taxids(tx_names, name_list)
        if name_result.is_err():
            echo_error(f"Failed to parse target names: {name_result.unwrap_err()}")
            raise SystemExit(1)
        name_taxids = name_result.unwrap()
        echo_info(f"Resolved target names → taxids: {name_taxids}")
        if len(name_taxids) != len(name_list):
            echo_warning(
                f"Note: {len(name_list)} name(s) resolved to {len(name_taxids)} taxid(s). "
                f"Use --tx_taxid with explicit IDs for precise control."
            )
        tx_ids = list(dict.fromkeys((tx_ids or []) + name_taxids))

    # Parse outgroup IDs
    parsed_outgroup_ids = None
    if tx_outgroup_ids:
        parsed_outgroup_ids = [int(t.strip()) for t in tx_outgroup_ids.split(",")]

    # Convert outgroup names → taxids
    if tx_outgroup_names:
        if not tx_names:
            echo_error("--tx_outgroup_names requires --tx_names")
            raise SystemExit(1)
        from eprobe.popgen.filter import parse_taxonomy_names_to_taxids
        og_list = [n.strip() for n in tx_outgroup_names.split(",")]
        og_result = parse_taxonomy_names_to_taxids(tx_names, og_list)
        if og_result.is_err():
            echo_error(f"Failed to parse outgroup names: {og_result.unwrap_err()}")
            raise SystemExit(1)
        og_taxids = og_result.unwrap()
        echo_info(f"Resolved outgroup names → taxids: {og_taxids}")
        if len(og_taxids) != len(og_list):
            echo_warning(
                f"Note: {len(og_list)} outgroup name(s) resolved to {len(og_taxids)} taxid(s). "
                f"Use --tx_outgroup_ids with explicit IDs for precise control."
            )
        parsed_outgroup_ids = list(dict.fromkeys((parsed_outgroup_ids or []) + og_taxids))

    step_list = list(steps)
    echo_info(f"Filtering {input}")
    echo_info(f"Steps: {' → '.join(s.upper() for s in step_list)}")

    result = run_util_filter(
        input_fasta=input,
        output_prefix=output,
        steps=step_list,
        threads=threads,
        bg_db=bg_db,
        ac_db=ac_db,
        ac_mode=ac_mode,
        ac_score_diff=ac_score_diff,
        ac_min_genomes=ac_min_genomes,
        tx_db=tx_db,
        tx_ids=tx_ids,
        tx_mode=tx_mode,
        names_dmp=str(tx_names) if tx_names else None,
        nodes_dmp=str(tx_nodes) if tx_nodes else None,
        acc2tax=str(tx_acc2tax) if tx_acc2tax else None,
        tx_min_edit=tx_minedit,
        tx_max_edit=tx_maxedit,
        tx_keep_hits=tx_keep_hits,
        tx_outgroup_ids=parsed_outgroup_ids,
        tx_exclude_mode=tx_exclude_mode,
        tx_exclude_val=tx_exclude_val,
        tx_target_min_sim=tx_target_min_sim,
        tx_target_max_nm=tx_target_max_nm,
        gc_range=(gc_min, gc_max),
        tm_range=(tm_min, tm_max),
        max_complexity=complexity,
        max_hairpin=hairpin,
        max_dimer=dimer,
        nn_table=nn_table,
        na_conc=na_conc,
        keep_temp=keep_temp,
        verbose=verbose,
    )

    if result.is_err():
        echo_error(f"Filtering failed: {result.unwrap_err()}")
        raise SystemExit(1)

    stats = result.unwrap()
    initial = stats["initial_count"]
    final = stats["final_count"]
    removed = stats["total_removed"]
    rate = (final / initial * 100) if initial > 0 else 0

    echo_success(f"\n→ Filter Summary:")
    echo_info(f"  Input:  {initial} sequences")
    for step_name, step_stats in stats["steps"].items():
        echo_info(f"  {step_name}: {step_stats['before']} → {step_stats['after']} "
                  f"(-{step_stats['removed']})")
    echo_info(f"  Output: {final} sequences ({rate:.1f}% pass rate)")
    echo_success(f"  → {output}.filtered.fa")
    if removed > 0:
        echo_info(f"  → {output}.rejected.fa ({removed} rejected)")
