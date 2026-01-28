"""
Pipeline execution engine.

Executes probe design pipelines defined in YAML configuration files.
Supports step tracking, resumption, and dry-run mode.
"""

import logging
import json
from pathlib import Path
from typing import Optional, Dict, Any, List
from dataclasses import dataclass, field, asdict
from datetime import datetime
from enum import Enum

import yaml

from eprobe.core.result import Result, Ok, Err

logger = logging.getLogger(__name__)


class StepStatus(Enum):
    """Pipeline step status."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass
class PipelineStep:
    """Representation of a pipeline step."""
    number: int
    name: str
    config: Dict[str, Any]
    status: StepStatus = StepStatus.PENDING
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    error: Optional[str] = None
    output_files: List[str] = field(default_factory=list)


@dataclass
class PipelineState:
    """Pipeline execution state."""
    name: str
    panel: str
    config_path: str
    output_dir: str
    started_at: str
    steps: List[PipelineStep] = field(default_factory=list)
    current_step: Optional[int] = None
    completed: bool = False


def load_config(config_path: Path) -> Result[Dict[str, Any], str]:
    """
    Load pipeline configuration from YAML file.
    
    Args:
        config_path: Path to YAML config
        
    Returns:
        Result containing parsed config
    """
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
        return Ok(config)
    except Exception as e:
        return Err(f"Failed to load config: {e}")


def save_state(state: PipelineState, output_dir: Path) -> None:
    """Save pipeline state to status.json."""
    status_path = output_dir / "status.json"
    
    state_dict = {
        "name": state.name,
        "panel": state.panel,
        "config_path": state.config_path,
        "output_dir": state.output_dir,
        "started_at": state.started_at,
        "current_step": state.current_step,
        "completed": state.completed,
        "steps": [
            {
                "number": s.number,
                "name": s.name,
                "status": s.status.value,
                "started_at": s.started_at,
                "completed_at": s.completed_at,
                "error": s.error,
                "output_files": s.output_files,
            }
            for s in state.steps
        ],
    }
    
    with open(status_path, 'w') as f:
        json.dump(state_dict, f, indent=2)


def load_state(output_dir: Path) -> Result[PipelineState, str]:
    """Load pipeline state from status.json."""
    status_path = output_dir / "status.json"
    
    if not status_path.exists():
        return Err(f"Status file not found: {status_path}")
    
    try:
        with open(status_path) as f:
            data = json.load(f)
        
        steps = [
            PipelineStep(
                number=s["number"],
                name=s["name"],
                config={},
                status=StepStatus(s["status"]),
                started_at=s.get("started_at"),
                completed_at=s.get("completed_at"),
                error=s.get("error"),
                output_files=s.get("output_files", []),
            )
            for s in data.get("steps", [])
        ]
        
        state = PipelineState(
            name=data["name"],
            panel=data["panel"],
            config_path=data["config_path"],
            output_dir=data["output_dir"],
            started_at=data["started_at"],
            steps=steps,
            current_step=data.get("current_step"),
            completed=data.get("completed", False),
        )
        
        return Ok(state)
        
    except Exception as e:
        return Err(f"Failed to load state: {e}")


def execute_popgen_step(
    step: PipelineStep,
    config: Dict[str, Any],
    input_config: Dict[str, Any],
    output_dir: Path,
    prev_output: Optional[Path],
    threads: int,
    dry_run: bool,
) -> Result[List[str], str]:
    """
    Execute a POPGEN panel step.
    
    Args:
        step: Pipeline step to execute
        config: Step configuration
        input_config: Pipeline input configuration
        output_dir: Output directory
        prev_output: Previous step output file
        threads: Number of threads
        dry_run: If True, don't actually run
        
    Returns:
        Result containing output file paths
    """
    step_dir = output_dir / f"{step.number:02d}_{step.name}"
    step_dir.mkdir(parents=True, exist_ok=True)
    
    output_prefix = step_dir / "output"
    
    if dry_run:
        logger.info(f"[DRY RUN] Would execute: {step.name}")
        return Ok([str(output_prefix) + ".tsv"])
    
    # Import step modules
    if step.name == "extract":
        from eprobe.popgen.extract import run_extract
        
        result = run_extract(
            vcf_path=Path(input_config["vcf"]),
            reference_path=Path(input_config["reference"]),
            output_prefix=output_prefix,
            flank=config.get("flank", 60),
            cluster_mode=config.get("cluster_mode", "on") == "on",
            cluster_flank=config.get("cluster_flank", 60),
            max_cluster_snp=config.get("max_cluster_snp", 3),
            threads=threads,
        )
        
    elif step.name == "filter":
        from eprobe.popgen.filter import run_filter
        
        input_file = prev_output or Path(input_config.get("snp_tsv", ""))
        
        result = run_filter(
            input_path=input_file,
            reference_path=Path(input_config["reference"]),
            output_prefix=output_prefix,
            filters=config.get("filters", ["biophysical"]),
            bg_db=Path(config["bg_db"]) if config.get("bg_db") else None,
            ac_db=Path(config["ac_db"]) if config.get("ac_db") else None,
            gc_range=tuple(config.get("gc_range", [35, 65])),
            tm_range=tuple(config.get("tm_range", [55, 75])),
            max_complexity=config.get("max_complexity", 2.0),
            threads=threads,
        )
        
    elif step.name == "select":
        from eprobe.popgen.select import run_select
        
        input_file = prev_output or Path(input_config.get("snp_tsv", ""))
        
        result = run_select(
            input_path=input_file,
            output_prefix=output_prefix,
            strategy=config.get("strategy", "uniform"),
            probe_number=config.get("probe_number", 10000),
            window_size=config.get("window_size", 100000),
        )
        
    elif step.name == "build":
        from eprobe.popgen.build import run_build
        
        input_file = prev_output or Path(input_config.get("snp_tsv", ""))
        
        result = run_build(
            input_path=input_file,
            reference_path=Path(input_config["reference"]),
            output_prefix=output_prefix,
            length=config.get("length", 81),
            shift=config.get("shift", 0),
            replace_mode=config.get("replace_mode", "on") == "on",
        )
        
    else:
        return Err(f"Unknown step: {step.name}")
    
    if result.is_err():
        return Err(result.unwrap_err())
    
    stats = result.unwrap()
    output_file = stats.get("output_file") or stats.get("fasta_file", "")
    
    return Ok([output_file])


def run_pipeline(
    config_path: Path,
    output_dir: Path,
    dry_run: bool = False,
    skip_completed: bool = True,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Execute pipeline from configuration file.
    
    Args:
        config_path: Path to YAML config
        output_dir: Output directory
        dry_run: If True, show what would be done without running
        skip_completed: Skip steps with existing outputs
        threads: Number of threads
        verbose: Enable verbose logging
        
    Returns:
        Result containing execution statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    # Load configuration
    config_result = load_config(config_path)
    if config_result.is_err():
        return Err(config_result.unwrap_err())
    
    config = config_result.unwrap()
    
    # Extract pipeline info
    pipeline_config = config.get("pipeline", {})
    name = pipeline_config.get("name", "Unnamed Pipeline")
    panel = pipeline_config.get("panel", "popgen")
    
    input_config = config.get("input", {})
    steps_config = config.get("steps", [])
    
    logger.info(f"Running pipeline: {name}")
    logger.info(f"Panel: {panel}")
    logger.info(f"Steps: {len(steps_config)}")
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize state
    state = PipelineState(
        name=name,
        panel=panel,
        config_path=str(config_path),
        output_dir=str(output_dir),
        started_at=datetime.now().isoformat(),
        steps=[],
    )
    
    # Parse steps
    for i, step_config in enumerate(steps_config):
        if isinstance(step_config, dict):
            step_name = list(step_config.keys())[0]
            step_params = step_config[step_name] or {}
        else:
            step_name = step_config
            step_params = {}
        
        step = PipelineStep(
            number=i + 1,
            name=step_name,
            config=step_params,
        )
        state.steps.append(step)
    
    save_state(state, output_dir)
    
    # Execute steps
    prev_output: Optional[Path] = None
    completed_steps = 0
    
    for step in state.steps:
        state.current_step = step.number
        step.status = StepStatus.RUNNING
        step.started_at = datetime.now().isoformat()
        save_state(state, output_dir)
        
        logger.info(f"Step {step.number}/{len(state.steps)}: {step.name}")
        
        if panel == "popgen":
            result = execute_popgen_step(
                step=step,
                config=step.config,
                input_config=input_config,
                output_dir=output_dir,
                prev_output=prev_output,
                threads=threads,
                dry_run=dry_run,
            )
        else:
            result = Err(f"Unsupported panel: {panel}")
        
        if result.is_err():
            step.status = StepStatus.FAILED
            step.error = result.unwrap_err()
            save_state(state, output_dir)
            return Err(f"Step {step.name} failed: {step.error}")
        
        output_files = result.unwrap()
        step.output_files = output_files
        step.status = StepStatus.COMPLETED
        step.completed_at = datetime.now().isoformat()
        completed_steps += 1
        
        # Set output for next step
        if output_files:
            prev_output = Path(output_files[0])
        
        save_state(state, output_dir)
    
    state.completed = True
    save_state(state, output_dir)
    
    logger.info(f"Pipeline completed: {completed_steps}/{len(state.steps)} steps")
    
    return Ok({
        "name": name,
        "panel": panel,
        "total_steps": len(state.steps),
        "completed_steps": completed_steps,
        "output_dir": str(output_dir),
    })


def resume_pipeline(
    output_dir: Path,
    from_step: Optional[str] = None,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Resume an interrupted pipeline.
    
    Args:
        output_dir: Pipeline output directory
        from_step: Optional step name or number to restart from
        threads: Number of threads
        verbose: Enable verbose logging
        
    Returns:
        Result containing execution statistics
    """
    # Load state
    state_result = load_state(output_dir)
    if state_result.is_err():
        return Err(state_result.unwrap_err())
    
    state = state_result.unwrap()
    
    # Load original config
    config_result = load_config(Path(state.config_path))
    if config_result.is_err():
        return Err(config_result.unwrap_err())
    
    config = config_result.unwrap()
    
    # Find starting step
    start_index = 0
    
    if from_step:
        for i, step in enumerate(state.steps):
            if step.name == from_step or str(step.number) == from_step:
                start_index = i
                break
    else:
        # Find first incomplete step
        for i, step in enumerate(state.steps):
            if step.status != StepStatus.COMPLETED:
                start_index = i
                break
    
    logger.info(f"Resuming from step {start_index + 1}: {state.steps[start_index].name}")
    
    # Re-run pipeline with updated config
    return run_pipeline(
        config_path=Path(state.config_path),
        output_dir=output_dir,
        skip_completed=True,
        threads=threads,
        verbose=verbose,
    )


def get_pipeline_status(output_dir: Path) -> Result[Dict[str, Any], str]:
    """
    Get pipeline status.
    
    Args:
        output_dir: Pipeline output directory
        
    Returns:
        Result containing status information
    """
    state_result = load_state(output_dir)
    if state_result.is_err():
        return Err(state_result.unwrap_err())
    
    state = state_result.unwrap()
    
    completed_steps = sum(1 for s in state.steps if s.status == StepStatus.COMPLETED)
    
    return Ok({
        "name": state.name,
        "panel": state.panel,
        "started": state.started_at,
        "total_steps": len(state.steps),
        "completed_steps": completed_steps,
        "steps": [
            {
                "number": s.number,
                "name": s.name,
                "status": s.status.value,
                "error": s.error,
            }
            for s in state.steps
        ],
    })
