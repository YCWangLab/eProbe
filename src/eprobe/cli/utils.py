"""
CLI utility functions.

Common helpers for CLI commands including output formatting,
file validation, and error handling.
"""

import click
from pathlib import Path
from typing import Optional


# Color definitions for consistent styling
COLORS = {
    "success": "green",
    "error": "red",
    "warning": "yellow",
    "info": "blue",
    "highlight": "cyan",
}


def echo_success(message: str) -> None:
    """Print success message with green checkmark."""
    click.echo(click.style("✓ ", fg=COLORS["success"]) + message)


def echo_error(message: str) -> None:
    """Print error message with red X."""
    click.echo(click.style("✗ ", fg=COLORS["error"]) + message, err=True)


def echo_warning(message: str) -> None:
    """Print warning message with yellow exclamation."""
    click.echo(click.style("! ", fg=COLORS["warning"]) + message)


def echo_info(message: str) -> None:
    """Print info message with blue arrow."""
    click.echo(click.style("→ ", fg=COLORS["info"]) + message)


def echo_step(step: int, total: int, message: str) -> None:
    """Print step progress message."""
    progress = click.style(f"[{step}/{total}]", fg=COLORS["highlight"])
    click.echo(f"{progress} {message}")


def validate_file_exists(path: Path, description: str = "File") -> bool:
    """
    Validate that a file exists and print appropriate message.
    
    Args:
        path: Path to validate
        description: Description for error message
        
    Returns:
        True if file exists, False otherwise
    """
    if not path.exists():
        echo_error(f"{description} not found: {path}")
        return False
    return True


def validate_output_path(path: Path, overwrite: bool = False) -> bool:
    """
    Validate output path and create parent directories if needed.
    
    Args:
        path: Output path to validate
        overwrite: Whether to allow overwriting existing files
        
    Returns:
        True if path is valid for writing, False otherwise
    """
    if path.exists() and not overwrite:
        echo_error(f"Output file already exists: {path}")
        echo_info("Use --force to overwrite")
        return False
    
    # Create parent directories
    path.parent.mkdir(parents=True, exist_ok=True)
    return True


def format_number(n: int) -> str:
    """Format large numbers with thousand separators."""
    return f"{n:,}"


def format_percentage(value: float, decimals: int = 1) -> str:
    """Format value as percentage string."""
    return f"{value:.{decimals}f}%"


def format_duration(seconds: float) -> str:
    """Format duration in human-readable format."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}h"


class ProgressReporter:
    """
    Simple progress reporter for long-running operations.
    
    Usage:
        with ProgressReporter("Processing", total=100) as progress:
            for i in range(100):
                do_work()
                progress.update(i + 1)
    """
    
    def __init__(self, task: str, total: Optional[int] = None, quiet: bool = False):
        self.task = task
        self.total = total
        self.current = 0
        self.quiet = quiet
    
    def __enter__(self) -> "ProgressReporter":
        if not self.quiet:
            if self.total:
                click.echo(f"{self.task} (0/{self.total})...", nl=False)
            else:
                click.echo(f"{self.task}...", nl=False)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        if not self.quiet:
            click.echo(" done")
    
    def update(self, current: int) -> None:
        """Update progress counter."""
        self.current = current
        if not self.quiet and self.total:
            # Clear line and reprint
            click.echo(f"\r{self.task} ({current}/{self.total})...", nl=False)
    
    def increment(self) -> None:
        """Increment progress by 1."""
        self.update(self.current + 1)


def confirm_overwrite(path: Path) -> bool:
    """
    Ask user to confirm overwriting existing file.
    
    Args:
        path: Path that would be overwritten
        
    Returns:
        True if user confirms, False otherwise
    """
    if not path.exists():
        return True
    
    return click.confirm(
        click.style(f"File exists: {path}. Overwrite?", fg=COLORS["warning"])
    )
