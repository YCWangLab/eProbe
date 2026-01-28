"""
Pipeline runner CLI commands.

Commands for running multi-step pipelines:
  run      - Execute YAML-defined pipeline
  resume   - Resume interrupted pipeline
  status   - Check pipeline status
"""

import click
from pathlib import Path
from typing import Optional

from eprobe.cli.utils import (
    echo_success,
    echo_error,
    echo_info,
    echo_warning,
    echo_step,
)


@click.group()
@click.pass_context
def run(ctx: click.Context) -> None:
    """
    Pipeline runner for multi-step workflows.
    
    Execute complete probe design pipelines defined in YAML config files.
    Supports resuming interrupted pipelines and checking status.
    
    \b
    Example:
      eprobe run pipeline -c config.yaml -o results/
      eprobe run resume -d results/
    """
    pass


@run.command("pipeline")
@click.option(
    "-c", "--config",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Pipeline configuration YAML file.",
)
@click.option(
    "-o", "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output directory.",
)
@click.option(
    "--dry_run/--no_dry_run",
    default=False,
    help="Show what would be executed without running.",
)
@click.option(
    "--skip_completed/--no_skip_completed",
    default=True,
    help="Skip steps with existing outputs (default: yes).",
)
@click.option(
    "-t", "--threads",
    default=1,
    type=int,
    help="Number of threads (default: 1).",
)
@click.pass_context
def pipeline(
    ctx: click.Context,
    config: Path,
    output: Path,
    dry_run: bool,
    skip_completed: bool,
    threads: int,
) -> None:
    """
    Execute a pipeline from YAML configuration.
    
    \b
    Config file format:
      pipeline:
        name: "My Probe Design"
        panel: popgen  # or funcgen
        
      input:
        vcf: data/input.vcf.gz
        reference: data/ref.fa
        
      steps:
        - extract:
            flank: 60
            cluster_mode: on
        - filter:
            filters: [BG]
            bg_db: /path/to/kraken_db
        - select:
            strategy: uniform
            probe_number: 10000
        - build:
            length: 81
            shift: 0
    
    \b
    Output structure:
      {output}/
        01_extract/
        02_filter/
        03_select/
        04_build/
        pipeline.log
        status.json
    """
    from eprobe.pipeline.runner import run_pipeline
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Loading pipeline from {config}")
    
    if dry_run:
        echo_warning("DRY RUN - no commands will be executed")
    
    result = run_pipeline(
        config_path=config,
        output_dir=output,
        dry_run=dry_run,
        skip_completed=skip_completed,
        threads=threads,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Pipeline failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Pipeline completed: {stats['completed_steps']}/{stats['total_steps']} steps")


@run.command()
@click.option(
    "-d", "--directory",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Pipeline output directory to resume.",
)
@click.option(
    "--from_step",
    type=str,
    help="Restart from specific step (e.g., 'filter', '02').",
)
@click.option(
    "-t", "--threads",
    default=1,
    type=int,
    help="Number of threads (default: 1).",
)
@click.pass_context
def resume(
    ctx: click.Context,
    directory: Path,
    from_step: Optional[str],
    threads: int,
) -> None:
    """
    Resume an interrupted pipeline.
    
    Reads status.json from the pipeline directory and resumes
    from the last incomplete step, or from a specified step.
    
    \b
    Example:
      eprobe run resume -d results/          # Resume from last incomplete
      eprobe run resume -d results/ --from_step filter  # Restart from filter
    """
    from eprobe.pipeline.runner import resume_pipeline
    
    verbose = ctx.obj.get("verbose", False)
    
    echo_info(f"Resuming pipeline from {directory}")
    
    if from_step:
        echo_info(f"Restarting from step: {from_step}")
    
    result = resume_pipeline(
        output_dir=directory,
        from_step=from_step,
        threads=threads,
        verbose=verbose,
    )
    
    if result.is_err():
        echo_error(f"Resume failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    stats = result.unwrap()
    echo_success(f"Pipeline completed: {stats['completed_steps']}/{stats['total_steps']} steps")


@run.command()
@click.option(
    "-d", "--directory",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Pipeline output directory.",
)
@click.option(
    "--json/--no_json",
    default=False,
    help="Output as JSON.",
)
@click.pass_context
def status(
    ctx: click.Context,
    directory: Path,
    json: bool,
) -> None:
    """
    Check pipeline status.
    
    Shows which steps have completed, failed, or are pending.
    
    \b
    Example:
      eprobe run status -d results/
    """
    from eprobe.pipeline.runner import get_pipeline_status
    import json as json_lib
    
    result = get_pipeline_status(directory)
    
    if result.is_err():
        echo_error(f"Status check failed: {result.unwrap_err()}")
        raise SystemExit(1)
    
    status_info = result.unwrap()
    
    if json:
        click.echo(json_lib.dumps(status_info, indent=2))
        return
    
    # Pretty print status
    echo_info(f"Pipeline: {status_info.get('name', 'Unknown')}")
    echo_info(f"Panel: {status_info.get('panel', 'Unknown')}")
    echo_info(f"Started: {status_info.get('started', 'N/A')}")
    echo_info("")
    echo_info("Steps:")
    
    for step in status_info.get("steps", []):
        step_num = step.get("number", "?")
        step_name = step.get("name", "Unknown")
        step_status = step.get("status", "unknown")
        
        if step_status == "completed":
            status_icon = click.style("✓", fg="green")
        elif step_status == "failed":
            status_icon = click.style("✗", fg="red")
        elif step_status == "running":
            status_icon = click.style("●", fg="yellow")
        else:
            status_icon = click.style("○", fg="white")
        
        click.echo(f"  {status_icon} {step_num}. {step_name}: {step_status}")
        
        if step.get("error"):
            echo_error(f"      Error: {step['error']}")
    
    echo_info("")
    
    total = status_info.get("total_steps", 0)
    completed = status_info.get("completed_steps", 0)
    
    if completed == total:
        echo_success(f"All {total} steps completed")
    else:
        echo_info(f"Progress: {completed}/{total} steps")
