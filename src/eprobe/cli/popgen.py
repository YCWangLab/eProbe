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
      4. assess   - Evaluate selected SNPs quality
      5. build    - Generate probe sequences
    
    \b
    Example workflow:
      eprobe popgen extract -v input.vcf.gz -r ref.fa -o project/step1
      eprobe popgen filter -i project/step1.snps.tsv -r ref.fa -o project/step2
      eprobe popgen select -i project/step2.filtered.tsv -o project/step3
      eprobe popgen assess -i project/step3.selected.tsv -o project/step4
      eprobe popgen build -i project/step4.snps.tsv -r ref.fa -o project/probes
    """
    pass


@popgen.command()
@click.option(
    "-v", "--vcf",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input VCF file (MUST be bgzip compressed .vcf.gz with tabix index .tbi).",
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
    
    Reads a compressed, indexed VCF file and extracts SNPs.
    For multi-allelic sites, always extracts the first alt allele (alt[0]).
    Optionally filters out SNP clusters (regions with too many SNPs close together)
    as these indicate hypervariable regions that may cause probe design issues.
    
    \b
    Multi-allelic handling:
      - Multi-allelic sites: extracts first alt allele only
      - Example: REF=A, ALT=G,T → extracts A/G SNP
      - Indels (REF or ALT length > 1) are automatically skipped
    
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
    
    # Check VCF index exists
    vcf_index = Path(str(vcf) + ".tbi")
    if not vcf_index.exists():
        echo_error(f"VCF index not found: {vcf_index}")
        echo_error("Please create index with: tabix -p vcf {vcf}")
        raise SystemExit(1)
    
    # Report BED filtering step
    if keep_bed and remove_bed:
        echo_error("Cannot use both --keep_bed and --remove_bed simultaneously")
        raise SystemExit(1)
    
    bed_filter_type = None
    bed_path_used = None
    if keep_bed:
        echo_info(f"→ Step 1: Filtering VCF with --keep_bed: {keep_bed}")
        bed_filter_type = "keep"
        bed_path_used = keep_bed
    elif remove_bed:
        echo_info(f"→ Step 1: Filtering VCF with --remove_bed: {remove_bed}")
        bed_filter_type = "remove"
        bed_path_used = remove_bed
    else:
        echo_info(f"→ Step 1: Extracting SNPs from entire VCF")
    
    echo_info(f"Extracting SNPs from {vcf}")
    
    result = run_extract(
        vcf_path=vcf,
        reference_path=reference,
        output_prefix=output,
        threads=threads,
        cluster_flank=cluster_flank,
        max_cluster_snp=max_cluster_snp,
        cluster_mode=cluster_filter,
        bed_path=bed_path_used,
        bed_mode=bed_filter_type if bed_filter_type else "keep",
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Extraction failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    
    # Step 2: Report BED filtering (if applied)
    if stats.get('bed_applied'):
        bed_mode = stats.get('bed_mode')
        total_vcf = stats.get('total_in_vcf')
        kept = stats.get('bed_kept')
        removed = stats.get('bed_removed')
        
        if bed_mode == 'keep':
            echo_success(f"→ Step 2: BED filter (--keep_bed)")
            echo_info(f"  ├─ Total SNPs in VCF: {total_vcf:,}")
            echo_info(f"  ├─ Kept (in BED regions): {kept:,}")
            echo_info(f"  └─ Removed (outside BED): {removed:,}")
        elif bed_mode == 'remove':
            echo_success(f"→ Step 2: BED filter (--remove_bed)")
            echo_info(f"  ├─ Total SNPs in VCF: {total_vcf:,}")
            echo_info(f"  ├─ Removed (in BED regions): {removed:,}")
            echo_info(f"  └─ Kept (outside BED): {kept:,}")
    else:
        echo_success(f"→ Step 2: Extracted {stats['total_in_vcf']:,} SNPs from VCF")
    
    # Step 3: Report cluster filtering
    if cluster_filter:
        filtered_count = stats.get('cluster_removed', 0)
        final_count = stats['final_snp_count']
        echo_success(f"→ Step 3: Cluster filter")
        echo_info(f"  ├─ Removed (in clusters): {filtered_count:,}")
        echo_info(f"  └─ Final output: {final_count:,} SNPs")
    else:
        final_count = stats['final_snp_count']
        echo_success(f"→ Step 3: Final output: {final_count:,} SNPs")
        echo_info(f"  └─ Cluster filtering disabled")
    
    echo_info(f"Output: {stats['output_file']}")


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
    "--ac_mode",
    type=click.Choice(["strict", "relaxed"]),
    default="strict",
    help="Accessibility filter mode: 'strict' requires unique alignment, 'relaxed' allows multi-map with score difference (default: strict).",
)
@click.option(
    "--ac_score_diff",
    type=int,
    default=10,
    help="Minimum alignment score difference for relaxed mode (default: 10). "
         "This threshold determines how much better the best alignment must be compared to the second-best. "
         "With Bowtie2's default scoring (~6 points per mismatch), a value of 10 requires ~2 mismatches difference. "
         "Higher values (15+) = stricter filtering, lower values (5) = more permissive.",
)
@click.option(
    "--ac_min_genomes",
    type=int,
    default=None,
    help="Minimum number of genomes a probe must align to (default: all genomes). "
         "Set to lower value to allow probes missing from some genomes. "
         "Example: with 3 genomes, --ac_min_genomes 2 requires alignment to at least 2 genomes.",
)
@click.option(
    "--tx_db",
    type=str,
    help="Bowtie2 database(s) for taxonomic filtering (comma-separated).",
)
@click.option(
    "--tx_taxid",
    type=str,
    help="Target taxon ID(s) for taxonomic filtering (comma-separated, e.g., '4530,4531').",
)
@click.option(
    "--tx_target_names",
    type=str,
    help="Target taxonomy name(s) for filtering (comma-separated, e.g., 'Tacamahaca,Populus').",
)
@click.option(
    "--tx_names",
    type=click.Path(exists=True, path_type=Path),
    help="NCBI taxonomy names.dmp file (required for taxonomic filter).",
)
@click.option(
    "--tx_nodes",
    type=click.Path(exists=True, path_type=Path),
    help="NCBI taxonomy nodes.dmp file (required for taxonomic filter).",
)
@click.option(
    "--tx_acc2tax",
    type=click.Path(exists=True, path_type=Path),
    help="Accession to taxid mapping file (required for taxonomic filter).",
)
@click.option(
    "--tx_minedit",
    type=int,
    default=0,
    help="Minimum edit distance for ngsLCA (default: 0).",
)
@click.option(
    "--tx_maxedit",
    type=int,
    default=2,
    help="Maximum edit distance for ngsLCA (default: 2).",
)
@click.option(
    "--tx_keep_hits",
    type=int,
    default=100,
    help="Number of alignments to report per sequence (default: 100).",
)
@click.option(
    "--tx_mode",
    type=click.Choice(["lca", "besthit"]),
    default="lca",
    help="Taxonomic classification mode: 'lca' (ngsLCA, default) or 'besthit'.",
)
@click.option(
    "--tx_outgroup_ids",
    type=str,
    help="Outgroup taxon ID(s) for best-hit mode (comma-separated). "
         "If omitted, ALL non-target taxa are treated as outgroups.",
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
    help="Best-hit exclusion strategy: "
         "'sim' (reject if outgroup sim >= val), "
         "'nm' (reject if outgroup NM <= val), "
         "'diff' (reject if outgroup_NM - target_NM < val, default).",
)
@click.option(
    "--tx_exclude_val",
    type=float,
    default=2,
    help="Threshold for best-hit exclusion (default: 2). "
         "For 'sim': percentage, 'nm': edit distance, 'diff': NM delta.",
)
@click.option(
    "--tx_target_min_sim",
    type=float,
    default=90.0,
    help="Min target similarity %% for best-hit acceptance (default: 90).",
)
@click.option(
    "--tx_target_max_nm",
    type=int,
    default=5,
    help="Max target NM for best-hit acceptance (default: 5).",
)
@click.option(
    "--gc",
    type=str,
    default="35,65",
    help="GC content range (min,max percentage, default: 35,65). Set -1 to disable.",
)
@click.option(
    "--tm",
    type=str,
    default="65,85",
    help="Melting temperature range (min,max Celsius, default: 65,85). Set -1 to disable.",
)
@click.option(
    "--complexity",
    type=float,
    default=2.0,
    help="Maximum DUST complexity score (default: 2.0). Set -1 to disable.",
)
@click.option(
    "--hairpin",
    type=float,
    default=50.0,
    help="Hairpin filter threshold (default: 50.0, allows stems up to 7bp). "
         ">=1: absolute score threshold; <1: percentile (e.g. 0.95 keeps 95%%). "
         "Uses exponential k-mer continuity scoring: 4^(n-1) bonus, normalized by log4(L). "
         "Stem levels: 4bp~0.34, 5bp~1.35, 6bp~5.42, 7bp~21.67, 8bp~86.68. Set -1 to disable.",
)
@click.option(
    "--dimer",
    type=float,
    default=0.50,
    help="Smart dimer filter sensitivity (default: 0.50). "
         "Identifies groups of cross-hybridizing probes via canonical k-mer sharing "
         "and keeps only one representative per risk group. "
         "Value = minimum fraction of shared k-mers to form a dimer edge: "
         "0.05 = very sensitive (remove more), 0.30 = balanced, 0.50 = conservative. "
         "Set -1 to disable.",
)
@click.option(
    "--nn-table",
    type=click.Choice([
        "DNA_NN1", "DNA_NN2", "DNA_NN3", "DNA_NN4", "R_DNA_NN1"
    ]),
    default="R_DNA_NN1",
    help="Nearest-neighbor thermodynamic table for Tm calculation. "
         "R_DNA_NN1 (Sugimoto 1995, default) for RNA probes hybridizing to DNA targets. "
         "DNA_NN tables for DNA/DNA hybridization: DNA_NN1 (Breslauer 1986), "
         "DNA_NN2 (Sugimoto 1996), DNA_NN3 (Allawi 1997), DNA_NN4 (SantaLucia 1998).",
)
@click.option(
    "--na-conc",
    type=float,
    default=50.0,
    help="Sodium ion concentration in mM for Tm calculation (default: 50.0).",
)
@click.option(
    "--no-biophysical",
    is_flag=True,
    help="Skip biophysical filtering (GC, Tm, complexity, hairpin, dimer).",
)
@click.option(
    "--max_ambiguous",
    type=int,
    default=0,
    help="Max ambiguous bases (N, Y, R, W, etc.) per probe. "
         "0 (default): remove all probes with any ambiguous base. "
         ">0: replace ambiguous bases randomly if count <= N, else remove.",
)
@click.option(
    "--no-prefilter",
    is_flag=True,
    help="Skip the ambiguous base prefilter step entirely.",
)
@click.option(
    "--keep-temp",
    is_flag=True,
    help="Keep temporary files for debugging (SAM, BAM, FASTA, etc.).",
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
    no_biophysical: bool,
    max_ambiguous: int,
    no_prefilter: bool,
    keep_temp: bool,
) -> None:
    """
    Apply multi-stage filtering to SNPs.
    
    \b
    Filtering stages:
      0. Prefilter: Remove/replace probes with ambiguous bases (N, Y, R, etc.)
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
    
    # Parse GC and Tm ranges (support -1 to disable)
    if gc.strip() == "-1":
        gc_min, gc_max = -1.0, -1.0
    else:
        gc_min, gc_max = map(float, gc.split(","))
    
    if tm.strip() == "-1":
        tm_min, tm_max = -1.0, -1.0
    else:
        tm_min, tm_max = map(float, tm.split(","))
    
    # Parse taxonomy IDs (from both taxid and target names)
    tx_ids = None
    if tx_taxid:
        tx_ids = [int(t.strip()) for t in tx_taxid.split(",")]
    
    # If target names provided, convert to taxids
    if tx_target_names:
        if not tx_names:
            echo_error("--tx_target_names requires --tx_names (names.dmp file)")
            raise SystemExit(1)
        
        # Parse names.dmp and convert target names to taxids
        from eprobe.popgen.filter import parse_taxonomy_names_to_taxids
        target_name_list = [n.strip() for n in tx_target_names.split(",")]
        
        name_result = parse_taxonomy_names_to_taxids(tx_names, target_name_list)
        if name_result.is_err():
            echo_error(f"Failed to parse taxonomy names: {name_result.unwrap_err()}")
            raise SystemExit(1)
        
        name_taxids = name_result.unwrap()
        # Show resolved mapping (may be more taxids than names if duplicates exist)
        echo_info(f"Resolved target names → taxids: {name_taxids}")
        if len(name_taxids) != len(target_name_list):
            echo_warning(
                f"Note: {len(target_name_list)} name(s) resolved to {len(name_taxids)} taxid(s). "
                f"Use --tx_taxid with explicit IDs for precise control."
            )
        
        # Merge with taxids from --tx_taxid
        if tx_ids is None:
            tx_ids = name_taxids
        else:
            tx_ids.extend(name_taxids)
            tx_ids = list(dict.fromkeys(tx_ids))  # Deduplicate preserving order
            
    # Parse outgroup IDs for best-hit mode
    parsed_outgroup_ids = None
    if tx_outgroup_ids:
        parsed_outgroup_ids = [int(t.strip()) for t in tx_outgroup_ids.split(",")]
    
    # Convert outgroup names → taxids for best-hit mode
    if tx_outgroup_names:
        if not tx_names:
            echo_error("--tx_outgroup_names requires --tx_names (names.dmp file)")
            raise SystemExit(1)
        
        from eprobe.popgen.filter import parse_taxonomy_names_to_taxids
        outgroup_name_list = [n.strip() for n in tx_outgroup_names.split(",")]
        og_result = parse_taxonomy_names_to_taxids(tx_names, outgroup_name_list)
        if og_result.is_err():
            echo_error(f"Failed to parse outgroup names: {og_result.unwrap_err()}")
            raise SystemExit(1)
        
        og_taxids = og_result.unwrap()
        echo_info(f"Resolved outgroup names → taxids: {og_taxids}")
        if len(og_taxids) != len(outgroup_name_list):
            echo_warning(
                f"Note: {len(outgroup_name_list)} outgroup name(s) resolved to {len(og_taxids)} taxid(s). "
                f"Use --tx_outgroup_ids with explicit IDs for precise control."
            )
        
        if parsed_outgroup_ids is None:
            parsed_outgroup_ids = og_taxids
        else:
            parsed_outgroup_ids.extend(og_taxids)
            parsed_outgroup_ids = list(dict.fromkeys(parsed_outgroup_ids))
    
    # Validate best-hit mode requirements
    if tx_mode == "besthit" and tx_db:
        if not tx_acc2tax:
            echo_error("Best-hit mode requires --tx_acc2tax (accession→taxid mapping)")
            raise SystemExit(1)
        if not tx_names:
            echo_error("Best-hit mode requires --tx_names (names.dmp)")
            raise SystemExit(1)
        if not tx_nodes:
            echo_error("Best-hit mode requires --tx_nodes (nodes.dmp)")
            raise SystemExit(1)
        if tx_ids is None:
            echo_error("Best-hit mode requires --tx_taxid or --tx_target_names to identify target taxa")
            raise SystemExit(1)
    
    # Determine which filters to apply
    filters = []
    if bg_db:
        filters.append("bg")
    if ac_db:
        filters.append("ac")
    if tx_db:
        filters.append("tx")
    # Add biophysical filter unless explicitly disabled
    if not no_biophysical:
        filters.append("biophysical")
    
    echo_info(f"Filtering SNPs from {input}")
    echo_info(f"Active filters: {', '.join(filters)}")
    
    result = run_filter(
        input_path=input,
        reference_path=reference,
        output_prefix=output,
        filters=filters,
        threads=threads,
        probe_length=length,
        bg_db=bg_db,
        ac_db=ac_db,
        ac_mode=ac_mode,
        ac_score_diff=ac_score_diff,
        ac_min_genomes=ac_min_genomes,
        tx_db=tx_db,
        tx_ids=tx_ids,
        names_dmp=str(tx_names) if tx_names else None,
        nodes_dmp=str(tx_nodes) if tx_nodes else None,
        acc2tax=str(tx_acc2tax) if tx_acc2tax else None,
        tx_min_edit=tx_minedit,
        tx_max_edit=tx_maxedit,
        tx_keep_hits=tx_keep_hits,
        tx_mode=tx_mode,
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
        verbose=verbose,
        keep_temp=keep_temp,
        max_ambiguous=max_ambiguous,
        no_prefilter=no_prefilter,
    )
    
    if result.is_err():
        echo_error(f"Filtering failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    initial = stats['initial_count']
    final = stats['final_count']
    removed = stats['total_removed']
    pass_rate = (final / initial * 100) if initial > 0 else 0
    
    # Report with 3-step structure
    echo_success(f"\n→ Filtering Summary:")
    echo_info(f"  Input: {initial} SNPs")
    
    # Step 0: Prefilter
    if 'prefilter' in stats['filters_applied'] and 'prefilter_details' in stats:
        pf = stats['prefilter_details']
        echo_success(f"\n  Step 0: Prefilter (ambiguous bases)")
        echo_info(f"    ├─ Removed (ambiguous): {pf['removed']}")
        if pf['replaced'] > 0:
            echo_info(f"    ├─ Replaced: {pf['replaced']} probes ({pf['bases_replaced']} bases)")
        echo_info(f"    └─ Remaining: {pf['remaining']} SNPs")
    
    # Step 1: External filters
    external_filters = [f for f in ['BG', 'AC', 'TX'] if f in stats['filters_applied']]
    if external_filters:
        echo_success(f"\n  Step 1: External Filtering ({', '.join(external_filters)})")
        if 'bg_remaining' in stats:
            echo_info(f"    ├─ Background filter: {stats['bg_remaining']} SNPs remaining")
        if 'ac_remaining' in stats:
            echo_info(f"    ├─ Accessibility filter: {stats['ac_remaining']} SNPs remaining")
        if 'tx_remaining' in stats:
            echo_info(f"    ├─ Taxonomic filter: {stats['tx_remaining']} SNPs remaining")
            
            # Show per-taxid statistics if available
            if 'tx_taxid_counts' in stats and stats['tx_taxid_counts']:
                taxid_counts = stats['tx_taxid_counts']
                
                # Try to get taxon names if names.dmp available
                taxid_names = {}
                if tx_names:
                    from eprobe.popgen.filter import get_taxid_names
                    name_result = get_taxid_names(tx_names, list(taxid_counts.keys()))
                    if name_result.is_ok():
                        taxid_names = name_result.unwrap()
                
                echo_info(f"    │  Probes per taxonomy node:")
                for taxid in sorted(taxid_counts.keys()):
                    count = taxid_counts[taxid]
                    name = taxid_names.get(taxid, "")
                    if name:
                        echo_info(f"    │    ├─ {taxid} ({name}): {count} probes")
                    else:
                        echo_info(f"    │    ├─ {taxid}: {count} probes")
    
    # Step 2: Biophysical filter
    if 'biophysical' in stats['filters_applied'] and 'biophysical_details' in stats:
        bio_details = stats['biophysical_details']
        echo_success(f"\n  Step 2: Biophysical Filtering")
        
        # GC with distribution stats
        gc_stats = bio_details.get('gc_stats', {})
        gc_failed = bio_details.get('gc_failed', 0)
        if gc_stats:
            gc_range_str = f"{gc_min}-{gc_max}%" if gc_min >= 0 else "disabled"
            echo_info(f"    ├─ GC content failed: {gc_failed} "
                     f"(range={gc_range_str}, "
                     f"mean={gc_stats['mean']}%, median={gc_stats['median']}%, P95={gc_stats['p95']}%)")
        else:
            echo_info(f"    ├─ GC content failed: {gc_failed}")
        
        # Tm with distribution stats
        tm_stats = bio_details.get('tm_stats', {})
        tm_failed = bio_details.get('tm_failed', 0)
        if tm_stats:
            tm_range_str = f"{tm_min}-{tm_max}°C" if tm_min >= 0 else "disabled"
            echo_info(f"    ├─ Tm failed: {tm_failed} "
                     f"(range={tm_range_str}, "
                     f"mean={tm_stats['mean']}°C, median={tm_stats['median']}°C, P95={tm_stats['p95']}°C)")
        else:
            echo_info(f"    ├─ Tm failed: {tm_failed}")
        
        # Complexity with distribution stats
        dust_stats = bio_details.get('dust_stats', {})
        complexity_failed = bio_details.get('complexity_failed', 0)
        if dust_stats:
            dust_thr_str = f"max≤{complexity}" if complexity >= 0 else "disabled"
            echo_info(f"    ├─ Complexity (DUST) failed: {complexity_failed} "
                     f"({dust_thr_str}, "
                     f"mean={dust_stats['mean']}, median={dust_stats['median']}, P95={dust_stats['p95']})")
        else:
            echo_info(f"    ├─ Complexity (DUST) failed: {complexity_failed}")
        
        # Hairpin with distribution stats
        hp_stats = bio_details.get('hairpin_stats', {})
        hp_failed = bio_details.get('hairpin_failed', 0)
        if hp_stats:
            echo_info(f"    ├─ Hairpin failed: {hp_failed} "
                     f"(threshold={hp_stats['threshold']}, {hp_stats['mode']}, "
                     f"median={hp_stats['median']}, P95={hp_stats['p95']})")
        else:
            echo_info(f"    ├─ Hairpin failed: {hp_failed}")
        
        # Dimer with smart filter stats
        dm_stats = bio_details.get('dimer_stats', {})
        dm_failed = bio_details.get('dimer_failed', 0)
        if dm_stats:
            echo_info(f"    ├─ Dimer failed: {dm_failed} "
                     f"({dm_stats['mode']}, "
                     f"{dm_stats['dimer_groups']} risk groups, "
                     f"max_group={dm_stats['max_group_size']}, "
                     f"{dm_stats['dimer_edges']} edges)")
        else:
            echo_info(f"    ├─ Dimer failed: {dm_failed}")
        
        echo_info(f"    └─ Passed: {stats.get('biophysical_remaining', final)} SNPs")
    
    # Step 3: Output
    echo_success(f"\n  Step 3: Output Saved")
    echo_info(f"    ├─ Final SNPs: {final} ({pass_rate:.1f}% of input)")
    echo_info(f"    ├─ Total removed: {removed}")
    echo_info(f"    └─ File: {stats['output_file']}")


@popgen.command()
@click.option(
    "-i", "--input",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="Input filtered SNPs TSV file(s). Can specify multiple: -i file1 -i file2 ...",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix.",
)
@click.option(
    "--merge_mode",
    type=click.Choice(["intersection", "union", "difference", "symmetric_diff"]),
    default="intersection",
    help="Merge mode for multiple inputs: intersection (default, keep SNPs in ALL), "
         "union (combine all), difference (only in 1st), symmetric_diff (in one only).",
)
@click.option(
    "--window_size",
    default=10000,
    type=int,
    help="Window size for SNP selection in bp (default: 10000).",
)
@click.option(
    "--strategy",
    type=click.Choice(["random", "uniform", "weighted", "priority"]),
    default="random",
    help="Selection strategy: "
         "random (simple random sampling), "
         "uniform (even genomic distribution via windows), "
         "weighted (best biophysical score per window, requires filter output), "
         "priority (prioritize BED regions, requires --priority_bed).",
)
@click.option(
    "--weights",
    type=str,
    help="Weights for biophysical features: gc,tm,complexity,hairpin,dimer (default: 1,1,1,1,1).",
)
@click.option(
    "--targets",
    type=str,
    help="Target values for biophysical features: gc,tm,complexity,hairpin,dimer "
         "(default: 50,70,0,0,0). Scoring penalises distance from each target.",
)
@click.option(
    "--priority_bed",
    type=click.Path(exists=True, path_type=Path),
    help="BED file of priority regions (e.g., exons). Required for 'priority' strategy.",
)
@click.option(
    "--target_count",
    type=int,
    help="Target number of SNPs to select. If not specified, keeps all SNPs after merge.",
)
@click.option(
    "--seed",
    default=42,
    type=int,
    help="Random seed for reproducibility (default: 42).",
)
@click.option(
    "--threads", "-t",
    default=1,
    type=int,
    help="Number of threads for parallel operations: file loading and window selection (default: 1).",
)
@click.option(
    "--keep_biophysical",
    is_flag=True,
    default=False,
    help="Keep biophysical columns (gc, tm, complexity, hairpin, dimer) in output. Default: False.",
)
@click.pass_context
def select(
    ctx: click.Context,
    input: tuple,
    output: Path,
    merge_mode: str,
    window_size: int,
    strategy: str,
    weights: Optional[str],
    targets: Optional[str],
    priority_bed: Optional[Path],
    target_count: Optional[int],
    seed: int,
    threads: int,
    keep_biophysical: bool,
) -> None:
    """
    Select optimal SNPs using window-based sampling.

    Supports multiple input files with merge operations:
    - intersection (default): Keep SNPs in ALL files
    - union: Combine all SNPs
    - difference: Only in first file
    - symmetric_diff: In one file only

    \b
    Selection strategies:
      random   - Simple random sampling without genomic distribution constraints.
                 Randomly selects target_count SNPs from input pool.

      uniform  - Window-based uniform distribution across genome.
                 Divides genome into windows (--window_size), calculates
                 SNPs per window = target_count / total_windows, then
                 randomly samples within each window to ensure even coverage.

      weighted - Window-based selection using biophysical scores.
                 REQUIRES biophysical columns (gc, tm, complexity, hairpin, dimer).
                 Calculates composite score = Σ(weight_i × distance_to_target_i).
                 Default targets: gc=50%, tm=70°C, complexity=0, hairpin=0, dimer=0.
                 Selects top-scoring SNP(s) from each window.
                 Use --weights to adjust priorities, --targets to set optimal values.

      priority - Window-based selection with priority for specific genomic regions.
                 REQUIRES --priority_bed file defining priority regions (e.g., exons).
                 For EACH WINDOW: if SNPs exist in priority regions, select from those;
                 otherwise select from non-priority SNPs. This ensures uniform genomic
                 coverage while maximizing SNPs in priority regions.
                 Can combine with --weights for weighted selection within windows.

    \b
    Output files:
      {output}.selected.tsv  - Selected SNPs
    """
    from eprobe.popgen.select import run_select, merge_snp_files
    
    verbose = ctx.obj.get("verbose", False)
    
    if not input:
        echo_error("At least one --input file is required")
        raise SystemExit(1)
    
    # Parse weights if provided
    weight_values = None
    if weights:
        weight_values = list(map(float, weights.split(",")))
        if len(weight_values) != 5:
            echo_error("Weights must have 5 values: gc,tm,complexity,hairpin,dimer")
            raise SystemExit(1)

    # Parse targets if provided
    target_values = None
    if targets:
        target_values = list(map(float, targets.split(",")))
        if len(target_values) != 5:
            echo_error("Targets must have 5 values: gc,tm,complexity,hairpin,dimer")
            raise SystemExit(1)

    # Handle multiple input files
    merge_details = None
    biophysical_reference = None  # Path to file with biophysical columns

    if len(input) > 1:
        echo_info(f"Merging {len(input)} input files using {merge_mode} mode")
        for i, f in enumerate(input, 1):
            echo_info(f"  {i}. {f}")

        merge_result = merge_snp_files(list(input), merge_mode, threads=threads)
        if merge_result.is_err():
            echo_error(f"Failed to merge inputs: {merge_result.unwrap_err()}")
            raise SystemExit(1)

        merged_snps = merge_result.unwrap()
        echo_info(f"After merge ({merge_mode}): {len(merged_snps)} SNPs")

        # Find input file with biophysical columns for reference
        # Check each input file for biophysical columns (gc, tm, complexity, hairpin, dimer)
        biophysical_cols = {'gc', 'tm', 'complexity', 'hairpin', 'dimer'}
        for f in input:
            try:
                import pandas as pd
                df_head = pd.read_csv(f, sep='\t', nrows=1)
                if biophysical_cols.issubset(set(df_head.columns)):
                    biophysical_reference = f
                    break
            except Exception:
                continue

        input_path = biophysical_reference if biophysical_reference else input[0]
        merged_data = merged_snps

        # Create merge details for reporting
        merge_details = {
            "mode": merge_mode,
            "file_count": len(input),
            "files": [str(f) for f in input],
            "merged_count": len(merged_snps),
        }

        if biophysical_reference:
            echo_info(f"Using {biophysical_reference.name} for biophysical metrics")
    else:
        input_path = input[0]
        merged_data = None
    
    echo_info(f"Selecting SNPs from {input_path}")
    echo_info(f"Window size: {window_size}bp, Strategy: {strategy}, Threads: {threads}")
    if target_count:
        echo_info(f"Target count: {target_count}")
    
    result = run_select(
        input_path=input_path,
        output_prefix=output,
        window_size=window_size,
        strategy=strategy,
        weights=weight_values,
        targets=target_values,
        priority_bed=priority_bed,
        target_count=target_count,
        seed=seed,
        keep_biophysical=keep_biophysical,
        merged_snps=merged_data,
        merge_details=merge_details,
        threads=threads,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Selection failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"\n→ Selection Summary:")
    
    # Show merge info if available in stats
    if "merge_details" in stats:
        merge_info = stats["merge_details"]
        echo_info(f"  Merge (mode: {merge_info['mode']})")
        echo_info(f"    ├─ Input files: {merge_info['file_count']}")
        for i, f in enumerate(merge_info['files'], 1):
            echo_info(f"    │  {i}. {Path(f).name}")
        echo_info(f"    └─ After merge: {merge_info['merged_count']} SNPs")
    
    echo_info(f"  Input: {stats['initial_count']} SNPs")
    echo_info(f"  Selected: {stats['selected']} SNPs")
    echo_info(f"  Windows covered: {stats['windows']}")
    echo_info(f"  Strategy: {stats['strategy']}")
    echo_info(f"  Output: {stats['output_file']}")


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
    type=click.Choice(["center", "left", "right", "tiling"]),
    default="center",
    help="SNP position within probe (default: center). 'tiling' generates 3 probes per SNP.",
)
@click.option(
    "--offset",
    default=0,
    type=int,
    help="Position offset from center (-/+ for left/right shift).",
)
@click.option(
    "--tiling_offset",
    default=15,
    type=int,
    help="Shift for tiling left/right probes (default: 15bp). Only used with -s tiling.",
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
    tiling_offset: int,
    replace_snp: bool,
    mutation_type: str,
) -> None:
    """
    Generate probe sequences from selected SNPs.
    
    Extracts genomic sequences around each SNP to create probe sequences.
    
    \b
    SNP position options:
      center  - SNP at probe center (recommended for hybridization)
      left    - SNP shifted toward 5' end
      right   - SNP shifted toward 3' end
      tiling  - 3 probes per SNP: center + left-shifted + right-shifted
    
    \b
    Replace SNP option:
      When enabled, replaces the SNP base with a third base (not REF or ALT).
      This can reduce allele bias in capture efficiency.
    
    \b
    Output files:
      {output}.probes.fa          - Probe sequences in FASTA format
      {output}.build_summary.txt  - Build statistics
    """
    from eprobe.popgen.build import run_build
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Building probes from {input}")
    echo_info(f"Length: {length}bp, SNP position: {snp_position}")
    if snp_position == "tiling":
        echo_info(f"Tiling mode: 3 probes/SNP, offset=±{tiling_offset}bp")
    if replace_snp:
        echo_info("Replace SNP: ON (third-base replacement)")
    if mutation_type != "both":
        echo_info(f"Mutation type filter: {mutation_type}")
    
    result = run_build(
        input_path=input,
        reference_path=reference,
        output_prefix=output,
        probe_length=length,
        snp_position=snp_position,
        offset=offset,
        replace_snp=replace_snp,
        mutation_type=mutation_type if mutation_type != "both" else None,
        tiling_offset=tiling_offset,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Build failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Generated {stats['probe_count']} probes from {stats['snp_count']} SNPs")
    echo_info(f"  Transitions: {stats['ts_count']}, Transversions: {stats['tv_count']}")
    echo_info(f"Output: {stats['fasta_file']}")


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
    help="[tags mode] Reference genome FASTA (required for --mode tags).",
)
@click.option(
    "-v", "--vcf",
    type=click.Path(exists=True, path_type=Path),
    help="[distance/pca/sfs modes] Original VCF file.",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix.",
)
@click.option(
    "--mode",
    type=click.Choice(["tags", "distance", "pca", "sfs", "all"]),
    default="tags",
    help="Assessment mode (default: tags).",
)
@click.option(
    "--tags",
    type=str,
    default="gc,tm,complexity,hairpin,dimer",
    help="[tags mode] Tags to analyze (comma-separated).",
)
@click.option(
    "--plot/--no_plot",
    default=True,
    help="[tags/distance/pca/sfs modes] Generate distribution plots (default: yes).",
)
@click.option(
    "--max_vcf_sites",
    default=100000,
    type=int,
    help="[distance mode] Maximum sites to load from VCF (default: 100000).",
)
@click.option(
    "--pop_file",
    type=click.Path(exists=True, path_type=Path),
    help="[distance/pca/sfs modes] Population file (sample_id<tab>pop_id). "
         "Required for --mode sfs. Optional for --mode pca/distance (enables stratified subsampling and plot coloring).",
)
@click.option(
    "--n_samples_per_pop",
    default=100,
    type=int,
    help="[distance/pca modes] Max individuals to use per population (with --pop_file, stratified; "
         "without --pop_file, total random subsample). Default: 100. "
         "[sfs mode] Number of samples per population (default: 100, typical sfs usage: 5).",
)
@click.option(
    "--samples",
    type=str,
    default=None,
    help="[sfs/pca modes] Specific samples to use (comma-separated). Overrides random subsampling.",
)
@click.option(
    "--projection",
    type=str,
    default=None,
    help="[sfs mode] easySFS projection values (e.g., '10,10,10'). Auto-calculated if not specified.",
)
@click.option(
    "--pops",
    type=str,
    default=None,
    help="[sfs mode] Populations for SFS analysis (comma-separated, max 3). "
         "Samples labeled 'unknown' are excluded. Default: first 3 non-unknown populations.",
)
@click.option(
    "--seed",
    default=42,
    type=int,
    help="[distance/pca/sfs modes] Random seed for sample subsampling (default: 42).",
)
@click.option(
    "-t", "--threads",
    "threads",
    default=1,
    type=int,
    help="[distance/pca/sfs modes] Number of threads for PLINK and bcftools (default: 1).",
)
@click.option(
    "-c", "--compare",
    "compare",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="[tags mode] Second TSV (e.g. post-filter) to overlay distributions against the main input.",
)
@click.option(
    "--xlim",
    "xlim",
    multiple=True,
    metavar="TAG:MIN,MAX",
    help="[tags mode] Fix x-axis range for a tag (e.g. --xlim gc:0,1 --xlim tm:30,80). "
         "Default: auto (P99-clipped). Repeatable per tag.",
)
@click.option(
    "--ylim",
    "ylim",
    multiple=True,
    metavar="TAG:MIN,MAX",
    help="[tags mode] Fix y-axis (percent) range for a tag (e.g. --ylim gc:0,50). "
         "Default: auto. Repeatable per tag.",
)
@click.option(
    "--bins",
    "bins",
    multiple=True,
    metavar="TAG:N",
    help="[tags mode] Number of histogram bins for a tag (e.g. --bins gc:100 --bins tm:50). "
         "Default: 50 per tag. Repeatable per tag.",
)
@click.option(
    "-l", "--probe_length",
    default=81,
    type=int,
    help="[tags mode] Probe length when generating sequences from TSV (default: 81).",
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
    max_vcf_sites: int,
    pop_file: Optional[Path],
    n_samples_per_pop: int,
    samples: Optional[str],
    projection: Optional[str],
    pops: Optional[str],
    seed: int,
    threads: int,
    compare: Optional[Path],
    xlim: tuple,
    ylim: tuple,
    bins: tuple,
    probe_length: int,
) -> None:
    """
    Assess quality of probe set.
    
    \b
    Assessment modes:
      tags     - Plot biophysical property distributions (uses existing
                 columns from filter step, or calculates with --reference)
      distance - Compare 1-IBS distance matrices (probe vs full VCF)
                 Uses PLINK: 1-IBS = 1 - Identity By State (0=identical, 1=different)
                 Outputs: combined heatmap (lower triangle = probe, upper triangle = full)
                 Calculates: Pearson r, Spearman rho, Manhattan distance
      pca      - Compare Principal Component Analysis (probe vs full VCF)
                 Uses PLINK to compute PCA on both datasets
                 Supports population-based coloring with --pop_file
                 Calculates Procrustes similarity between PCAs
      sfs      - Compare multi-population joint Site Frequency Spectrum
                 using easySFS. Requires --pop_file with population assignments.
                 Generates heatmap plots for full VCF and probe set.
                 NOTE: Max 3 populations. 'unknown' samples are excluded.
                 Use --pops to select specific populations (e.g., --pops pop_A,pop_B).
      all      - Run all assessments
    
    \b
    Parameters by mode:
    
      ALWAYS REQUIRED:
        -i, --input        Input SNPs TSV file
        -o, --output       Output prefix
      
      TAGS MODE:
        -r, --reference    Reference genome FASTA (for sequence generation)
        --tags             Biophysical metrics to plot (default: gc,tm,complexity,hairpin,dimer)
        -c, --compare      Second TSV file to overlay for comparison
        --xlim, --ylim     Fix axis ranges (format: TAG:MIN,MAX, e.g., --xlim gc:0,1)
        --bins             Histogram bin count per tag (format: TAG:N, e.g., --bins gc:100)
        -l, --probe_length Probe length for sequence generation (default: 81)
      
      DISTANCE MODE (needs -v, --vcf):
        --max_vcf_sites    Max VCF sites to load (default: 100000)
        --pop_file         Population assignments (optional, for stratified subsampling)
        --n_samples_per_pop Max samples per population (default: 100)
      
      PCA MODE (needs -v, --vcf):
        --pop_file         Population assignments (optional, for coloring)
        --n_samples_per_pop Max samples per population (default: 100)
        --samples          Specific samples to use (comma-separated)
      
      SFS MODE (needs -v, --vcf and --pop_file):
        --pop_file         Population assignments (REQUIRED, format: sample_id<tab>pop_id)
        --n_samples_per_pop Samples per population (default: 100, typically 5 for SFS)
        --pops             Populations to analyze, max 3 (e.g., --pops pop_A,pop_B,pop_C)
        --projection       easySFS projection (e.g., --projection 10,10,10). If omitted, auto-calculated.
        --samples          Specific samples to use (comma-separated)
      
      COMMON TO DISTANCE/PCA/SFS:
        -v, --vcf          Original VCF file
        -t, --threads      Number of threads (default: 1)
        --seed             Random seed for subsampling (default: 42)
        --plot             Generate plots (default: yes, use --no_plot to disable)
    
    \b
    Output files (depends on mode):
      --mode tags:
        {output}.tags_summary.txt   - Biophysical statistics
        {output}.*_dist.jpg         - Distribution plots
      --mode distance:
        {output}.ibs_heatmap.png    - Combined heatmap
        {output}.ibs_full.tsv       - Full VCF 1-IBS matrix
        {output}.ibs_probe.tsv      - Probe set 1-IBS matrix
        {output}.distance_summary.txt - Correlation statistics
      --mode pca:
        {output}.pca_PC1_PC2.png    - Side-by-side PCA comparison (PC1 vs PC2)
        {output}.pca_full.tsv       - Full VCF eigenvectors
        {output}.pca_probe.tsv      - Probe set eigenvectors
        {output}.pca_summary.txt    - Variance explained and Procrustes similarity
      --mode sfs:
        {output}.sfs_full.png       - Full VCF multi-pop SFS heatmap
        {output}.sfs_probe.png      - Probe set multi-pop SFS heatmap
        {output}.sfs_summary.txt    - SFS correlations per population
        {output}_dadi/              - dadi-format SFS files
    
    \b
    SFS mode requires:
      --pop_file: Population assignment file (sample_id<tab>pop_id)
      --vcf: Original VCF file
      External tools: easySFS (set EASYSFS_PATH env var), bcftools
    
    \b
    PCA mode uses:
      --pop_file: (Optional) Population assignment for coloring
      --vcf: Original VCF file
      External tools: plink, bcftools
    """
    from eprobe.popgen.assess import run_assess
    
    verbose = ctx.obj.get("verbose", False)
    
    # Validate required inputs for each mode
    if mode in ["distance", "pca", "all"] and vcf is None:
        echo_error("--vcf is required for distance/pca assessment")
        raise SystemExit(1)
    
    if mode in ["sfs", "all"]:
        if vcf is None:
            echo_error("--vcf is required for SFS assessment")
            raise SystemExit(1)
        if pop_file is None:
            echo_error("--pop_file is required for SFS assessment")
            raise SystemExit(1)
    
    echo_info(f"Assessing probe set: {input}")
    echo_info(f"Mode: {mode}")
    
    # Parse specified samples
    specified_samples = None
    if samples is not None:
        specified_samples = [s.strip() for s in samples.split(",")]
    
    # Parse pops
    pops_list = None
    if pops is not None:
        pops_list = [p.strip() for p in pops.split(",")]

    # Parse --xlim / --ylim / --bins into dicts  (format: "tag:min,max" or "tag:n")
    def _parse_range(values: tuple) -> dict:
        out = {}
        for v in values:
            try:
                tag_part, rest = v.split(":", 1)
                lo, hi = rest.split(",", 1)
                out[tag_part.strip()] = (float(lo), float(hi))
            except Exception:
                echo_error(
                    f"Invalid format '{v}' — expected TAG:MIN,MAX (e.g. gc:0,1). Ignoring."
                )
        return out

    def _parse_bins(values: tuple) -> dict:
        out = {}
        for v in values:
            try:
                tag_part, n = v.split(":", 1)
                out[tag_part.strip()] = int(n.strip())
            except Exception:
                echo_error(
                    f"Invalid format '{v}' — expected TAG:N (e.g. gc:100). Ignoring."
                )
        return out

    plot_xlim = _parse_range(xlim) or None
    plot_ylim = _parse_range(ylim) or None
    plot_bins = _parse_bins(bins) or None

    result = run_assess(
        input_path=input,
        output_prefix=output,
        mode=mode,
        vcf_path=vcf,
        reference_path=reference,
        tags=tags.split(","),
        generate_plots=plot,
        probe_length=probe_length,
        max_vcf_sites=max_vcf_sites,
        pop_file=pop_file,
        n_samples_per_pop=n_samples_per_pop,
        specified_samples=specified_samples,
        projection=projection,
        pops=pops_list,
        seed=seed,
        threads=threads,
        verbose=verbose,
        compare_path=compare,
        plot_xlim=plot_xlim,
        plot_ylim=plot_ylim,
        plot_bins=plot_bins,
    )
    
    if result.is_err():
        echo_error(f"Assessment failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success("Assessment complete")
    
    # Print summary based on mode
    if "distance" in stats:
        dist_stats = stats["distance"]
        corr = dist_stats.get("correlation", {})
        echo_info(f"\n→ 1-IBS Distance Comparison:")
        echo_info(f"    ├─ Samples: {dist_stats.get('n_samples', 'N/A')}")
        echo_info(f"    ├─ Full VCF sites: {dist_stats.get('n_sites_full', 'N/A')}")
        echo_info(f"    ├─ Probe sites: {dist_stats.get('n_sites_probe', 'N/A')}")
        echo_info(f"    ├─ Pearson r: {corr.get('pearson_r', 'N/A'):.4f}")
        echo_info(f"    ├─ Spearman rho: {corr.get('spearman_r', 'N/A'):.4f}")
        echo_info(f"    └─ Manhattan: {corr.get('manhattan_distance', 'N/A'):.4f}")
    
    if "sfs" in stats:
        sfs_stats = stats["sfs"]
        echo_info(f"\n→ Site Frequency Spectrum (easySFS):")
        echo_info(f"    ├─ Populations: {', '.join(sfs_stats.get('populations', []))}")
        echo_info(f"    ├─ Samples per pop: {sfs_stats.get('n_samples_per_pop', 'N/A')}")
        echo_info(f"    ├─ Full VCF sites: {sfs_stats.get('n_sites_full', 'N/A'):,}")
        echo_info(f"    ├─ Probe sites: {sfs_stats.get('n_sites_probe', 'N/A'):,}")
        
        # Per-population correlations
        corrs = sfs_stats.get("sfs_correlations", {})
        if corrs:
            echo_info(f"    ├─ SFS correlations:")
            for pop, corr in corrs.items():
                echo_info(f"    │     {pop}: {corr:.4f}")
        
        avg_corr = sfs_stats.get("avg_correlation")
        if avg_corr is not None:
            echo_info(f"    └─ Average correlation: {avg_corr:.4f}")
        else:
            echo_info(f"    └─ Average correlation: N/A")
    
    if "pca" in stats:
        pca_stats = stats["pca"]
        echo_info(f"\n→ Principal Component Analysis (PLINK):")
        echo_info(f"    ├─ Samples: {pca_stats.get('n_samples', 'N/A')}")
        echo_info(f"    ├─ Probe positions: {pca_stats.get('n_probe_positions', 'N/A')}")
        
        # Variance explained for PC1-3
        var_full = pca_stats.get("variance_full", {})
        var_probe = pca_stats.get("variance_probe", {})
        if var_full and var_probe:
            echo_info(f"    ├─ Variance explained (PC1/PC2/PC3):")
            pc1_f = var_full.get("PC1", 0)
            pc2_f = var_full.get("PC2", 0)
            pc3_f = var_full.get("PC3", 0)
            pc1_p = var_probe.get("PC1", 0)
            pc2_p = var_probe.get("PC2", 0)
            pc3_p = var_probe.get("PC3", 0)
            echo_info(f"    │     Full VCF:  {pc1_f:.1f}% / {pc2_f:.1f}% / {pc3_f:.1f}%")
            echo_info(f"    │     Probe Set: {pc1_p:.1f}% / {pc2_p:.1f}% / {pc3_p:.1f}%")
        
        procrustes = pca_stats.get("procrustes_similarity")
        if procrustes is not None:
            echo_info(f"    └─ Procrustes similarity: {procrustes:.4f}")
        else:
            echo_info(f"    └─ Procrustes similarity: N/A")
    
    if "tags" in stats:
        tag_stats = stats["tags"]
        echo_info(f"\n→ Biophysical Tags:")
        echo_info(f"    ├─ Probes analyzed: {tag_stats.get('probe_count', 'N/A')}")
        echo_info(f"    └─ Summary: {tag_stats.get('summary_file', 'N/A')}")
    
    if plot:
        echo_info(f"\nPlots saved to {output}.*")
