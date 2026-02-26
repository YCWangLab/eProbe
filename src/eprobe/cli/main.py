"""
Main CLI entry point for eProbe.

Defines the root command group and registers all subcommands.
Uses Click framework for argument parsing and help generation.
"""

import click
from pathlib import Path
from typing import Optional

from eprobe import __version__


# Custom Click context settings for consistent behavior
CONTEXT_SETTINGS = dict(
    help_option_names=["-h", "--help"],
    max_content_width=120,
)


class AliasedGroup(click.Group):
    """
    Custom Click group that supports command aliases and underscore/hyphen interchangeability.
    
    Features:
    - Preserves command registration order (instead of alphabetical)
    - Allows users to type shortened versions of commands as long as unambiguous
    - Treats underscores and hyphens as equivalent (from_fasta == from-fasta)
    """

    def list_commands(self, ctx: click.Context) -> list:
        """Return commands in registration order."""
        return list(self.commands)

    def get_command(self, ctx: click.Context, cmd_name: str) -> Optional[click.Command]:
        # Normalize: convert underscores to hyphens for lookup
        normalized_name = cmd_name.replace("_", "-")
        
        # First try exact match with normalized name
        rv = click.Group.get_command(self, ctx, normalized_name)
        if rv is not None:
            return rv
        
        # Also try original name (in case command is registered with underscore)
        rv = click.Group.get_command(self, ctx, cmd_name)
        if rv is not None:
            return rv
        
        # Try prefix matching with normalized name
        matches = [x for x in self.list_commands(ctx) if x.startswith(normalized_name)]
        if not matches:
            # Also try with original name
            matches = [x for x in self.list_commands(ctx) if x.startswith(cmd_name)]
        
        if not matches:
            return None
        elif len(matches) == 1:
            return click.Group.get_command(self, ctx, matches[0])
        else:
            ctx.fail(f"Ambiguous command '{cmd_name}': could be {', '.join(sorted(matches))}")
            return None


@click.group(cls=AliasedGroup, context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__, "-V", "--version", prog_name="eprobe")
@click.option(
    "-v", "--verbose",
    is_flag=True,
    help="Enable verbose output with detailed logging.",
)
@click.option(
    "--quiet", "-q",
    is_flag=True,
    help="Suppress all output except errors.",
)
@click.pass_context
def cli(ctx: click.Context, verbose: bool, quiet: bool) -> None:
    """
    eProbe: Capture probe design toolkit for aDNA/eDNA.
    
    A comprehensive toolkit for designing capture probes for genetic diversity
    reconstructions from ancient and environmental DNA samples.
    
    \b
    Probe design panels:
      popgen   - SNP-based probes for population genetics
      funcgen  - Sequence-based probes for functional genetics
      barcode  - Kmer-based probes for taxonomic identification
                 with haploid genome (Mito/Plastid/Sex Chr) (developing)
    
    \b
    Post-design toolkit:
      util     - Merge, tile, adapt, rename, assess, sample, target
    
    \b
    Pipeline management:
      init     - Generate config templates for pipeline runner
      run      - Execute multi-step pipeline from config
      info     - Show dependencies and external tool versions
    
    \b
    Quick start:
      eprobe popgen extract -v input.vcf.gz -r ref.fa -o output
      eprobe popgen filter -i snps.tsv -r ref.fa -o filtered
      eprobe popgen build -i filtered.tsv -r ref.fa -o probes
    
    For detailed help on any command, use: eprobe <command> --help
    """
    # Ensure context object exists
    ctx.ensure_object(dict)
    
    # Store global options
    ctx.obj["verbose"] = verbose
    ctx.obj["quiet"] = quiet
    
    # Configure logging based on verbosity
    if verbose and not quiet:
        import logging
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        )
    elif quiet:
        import logging
        logging.basicConfig(level=logging.ERROR)


# Import and register subcommand groups (order defines help display order)
from eprobe.cli.popgen import popgen
from eprobe.cli.funcgen import funcgen
from eprobe.cli.util import util
from eprobe.cli.init import init
from eprobe.cli.run import run

cli.add_command(popgen)
cli.add_command(funcgen)
cli.add_command(util)
cli.add_command(init)
cli.add_command(run)


@cli.command()
@click.pass_context
def info(ctx: click.Context) -> None:
    """
    Show dependencies, external tools, and environment info.
    
    Reports versions of Python libraries and external bioinformatics
    tools required by various eProbe modules.
    """
    import sys
    import platform
    import subprocess
    import shutil
    
    click.echo(f"eProbe version: {__version__}")
    click.echo(f"Python: {sys.version.split()[0]}")
    click.echo(f"Platform: {platform.platform()}")
    
    # --- Python libraries ---
    click.echo("\nPython libraries:")
    
    libraries = [
        ("biopython", "Bio"),
        ("numpy", "numpy"),
        ("pandas", "pandas"),
        ("pysam", "pysam"),
        ("pybedtools", "pybedtools"),
        ("click", "click"),
        ("matplotlib", "matplotlib"),
        ("seaborn", "seaborn"),
        ("scipy", "scipy"),
        ("cyvcf2", "cyvcf2"),
        ("pyyaml", "yaml"),
    ]
    
    for name, import_name in libraries:
        try:
            module = __import__(import_name)
            version = getattr(module, "__version__", "unknown")
            click.echo(click.style(f"  ✓ {name}: {version}", fg="green"))
        except ImportError:
            click.echo(click.style(f"  ✗ {name}: not installed", fg="red"))
    
    # --- External tools ---
    click.echo("\nExternal tools:")
    
    tools = [
        ("bowtie2",     "bowtie2 --version",    0),
        ("samtools",    "samtools --version",    0),
        ("bedtools",    "bedtools --version",    0),
        ("bcftools",    "bcftools --version",    0),
        ("kraken2",     "kraken2 --version",     0),
        ("cd-hit-est",  "cd-hit-est -h",         None),  # exits non-zero
        ("plink",       "plink --version",       None),
        ("shapeit5",    "phase_common --help",   None),
        ("clustalo",    "clustalo --version",    0),
    ]
    
    for name, cmd, expected_rc in tools:
        exe = cmd.split()[0]
        if shutil.which(exe) is None:
            click.echo(click.style(f"  ✗ {name}: not found", fg="red"))
            continue
        try:
            r = subprocess.run(
                cmd.split(), capture_output=True, text=True, timeout=5,
            )
            # Extract version from first non-empty output line
            output = (r.stdout or r.stderr).strip()
            ver_line = ""
            for line in output.splitlines():
                line = line.strip()
                if line:
                    ver_line = line
                    break
            click.echo(click.style(f"  ✓ {name}: ", fg="green") + ver_line)
        except Exception:
            click.echo(click.style(f"  ? {name}: found but version unknown", fg="yellow"))


if __name__ == "__main__":
    cli()
