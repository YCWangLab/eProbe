"""
FUNCGEN probe quality assessment.

Same functionality as popgen assess but tailored for sequence-based probes.
Re-exports the core assess function from popgen with funcgen-specific defaults.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List

from eprobe.core.result import Result, Ok, Err
from eprobe.popgen.assess import (
    run_assess as popgen_run_assess,
    assess_probes,
    calculate_summary_statistics,
    generate_assessment_plots,
    AVAILABLE_TAGS,
    AssessmentResult,
)

logger = logging.getLogger(__name__)


def run_assess(
    fasta_path: Path,
    output_prefix: Path,
    tags: Optional[List[str]] = None,
    generate_plots: bool = True,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Assess quality of FUNCGEN probe set.
    
    Wrapper around popgen assess with FUNCGEN-specific defaults.
    For FUNCGEN probes, we typically don't need dimer analysis
    since probes are designed for specific targets.
    
    Args:
        fasta_path: Input probe FASTA file
        output_prefix: Output file prefix
        tags: List of metrics to calculate
        generate_plots: Generate distribution plots
        threads: Number of threads
        verbose: Enable verbose logging
        
    Returns:
        Result containing assessment statistics
    """
    # Default tags for FUNCGEN (exclude dimer by default)
    if tags is None:
        tags = ["gc", "tm", "complexity", "hairpin"]
    
    return popgen_run_assess(
        input_path=fasta_path,
        output_prefix=output_prefix,
        tags=tags,
        generate_plots=generate_plots,
        sample_dimer=1000,
        threads=threads,
        verbose=verbose,
    )
