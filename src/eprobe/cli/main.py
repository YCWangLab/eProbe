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
    - Allows users to type shortened versions of commands as long as unambiguous
    - Treats underscores and hyphens as equivalent (from_fasta == from-fasta)
    """
    
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
    Two main panels are available:
      popgen   - SNP-based probes for population genetics
      funcgen  - Sequence-based probes for functional genetics
    
    \b
    Utility tools:
      util     - Data manipulation (tiling, merge, dedup, subset)
      run      - One-click pipeline execution
      init     - Initialize configuration templates
    
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


# Import and register subcommand groups
from eprobe.cli.popgen import popgen
from eprobe.cli.funcgen import funcgen
from eprobe.cli.util import util
from eprobe.cli.run import run
from eprobe.cli.init import init

cli.add_command(popgen)
cli.add_command(funcgen)
cli.add_command(util)
cli.add_command(run)
cli.add_command(init)


@cli.command()
@click.pass_context
def info(ctx: click.Context) -> None:
    """
    Display version and environment information.
    
    Shows eProbe version, Python version, and installed dependencies.
    """
    import sys
    import platform
    
    click.echo(f"eProbe version: {__version__}")
    click.echo(f"Python version: {sys.version}")
    click.echo(f"Platform: {platform.platform()}")
    
    click.echo("\nInstalled dependencies:")
    
    # Check key dependencies (map names to import names)
    dependencies = {
        "biopython": "Bio",
        "pandas": "pandas",
        "numpy": "numpy",
        "pysam": "pysam",
        "pybedtools": "pybedtools",
        "click": "click",
        "matplotlib": "matplotlib",
        "seaborn": "seaborn",
    }
    
    for name, import_name in dependencies.items():
        try:
            module = __import__(import_name)
            version = getattr(module, "__version__", "unknown")
            click.echo(f"  {name}: {version}")
        except ImportError:
            click.echo(f"  {name}: not installed")


if __name__ == "__main__":
    cli()
