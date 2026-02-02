"""
BAM file processing utilities.

Provides optimized BAM compression and processing using:
  - compressbam (fast C++ BGZF compression, Linux only)
  - samtools fallback for cross-platform compatibility

This module integrates the fast BAM processing workflow from
the eDNA_taxonomic_identification pipeline into eProbe.
"""

import logging
import os
import platform
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Tuple

from eprobe.core.result import Result, Ok, Err

logger = logging.getLogger(__name__)


def get_compressbam_path() -> Optional[Path]:
    """
    Get path to compressbam binary if available.
    
    Checks:
    1. eprobe/bin/compressbam (bundled with package)
    2. PATH (system-installed)
    
    Returns:
        Path to compressbam binary, or None if not available
    """
    # Only available on Linux (ELF binary)
    if platform.system() != "Linux":
        logger.debug(f"compressbam not available on {platform.system()} (Linux only)")
        return None
    
    # Check bundled binary
    package_dir = Path(__file__).parent.parent
    bundled_path = package_dir / "bin" / "compressbam"
    
    if bundled_path.exists() and os.access(bundled_path, os.X_OK):
        logger.debug(f"Using bundled compressbam: {bundled_path}")
        return bundled_path
    
    # Check PATH
    compressbam = shutil.which("compressbam")
    if compressbam:
        logger.debug(f"Using system compressbam: {compressbam}")
        return Path(compressbam)
    
    logger.debug("compressbam not found, will use samtools fallback")
    return None


def compress_sam_with_compressbam(
    sam_path: Path,
    bam_path: Path,
    threads: int = 1,
) -> Result[Path, str]:
    """
    Convert SAM to BAM using compressbam (fast BGZF compression).
    
    Args:
        sam_path: Input SAM file
        bam_path: Output BAM file
        threads: Number of threads for compression
        
    Returns:
        Result containing path to output BAM file
    """
    compressbam = get_compressbam_path()
    if compressbam is None:
        return Err("compressbam not available")
    
    # compressbam needs an uncompressed BAM as input
    # First convert SAM to uncompressed BAM, then compress
    
    # Create temp uncompressed BAM
    temp_bam = sam_path.with_suffix(".uncomp.bam")
    
    try:
        # SAM -> uncompressed BAM using samtools
        sam_to_bam_cmd = [
            "samtools", "view",
            "-@", str(threads),
            "-b",
            "-u",  # Uncompressed BAM
            "-o", str(temp_bam),
            str(sam_path),
        ]
        
        subprocess.run(sam_to_bam_cmd, capture_output=True, text=True, check=True)
        
        # Compress with compressbam
        compress_cmd = [
            str(compressbam),
            "--threads", str(threads),
            "--input", str(temp_bam),
            "--output", str(bam_path),
        ]
        
        subprocess.run(compress_cmd, capture_output=True, text=True, check=True)
        
        # Verify output
        if not bam_path.exists() or bam_path.stat().st_size == 0:
            return Err(f"compressbam produced empty output: {bam_path}")
        
        return Ok(bam_path)
        
    except subprocess.CalledProcessError as e:
        return Err(f"BAM compression failed: {e.stderr}")
    except FileNotFoundError as e:
        return Err(f"Tool not found: {e}")
    finally:
        # Clean up temp file
        if temp_bam.exists():
            temp_bam.unlink()


def compress_sam_with_samtools(
    sam_path: Path,
    bam_path: Path,
    threads: int = 1,
) -> Result[Path, str]:
    """
    Convert SAM to BAM using samtools.
    
    Args:
        sam_path: Input SAM file
        bam_path: Output BAM file
        threads: Number of threads for compression
        
    Returns:
        Result containing path to output BAM file
    """
    cmd = [
        "samtools", "view",
        "-@", str(threads),
        "-b",
        "-o", str(bam_path),
        str(sam_path),
    ]
    
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if not bam_path.exists():
            return Err(f"samtools produced no output: {bam_path}")
        
        return Ok(bam_path)
        
    except subprocess.CalledProcessError as e:
        return Err(f"SAM to BAM conversion failed: {e.stderr}")
    except FileNotFoundError:
        return Err("samtools not found. Please install samtools and ensure it's in PATH.")


def sam_to_bam(
    sam_path: Path,
    bam_path: Optional[Path] = None,
    threads: int = 1,
    use_compressbam: bool = True,
    remove_sam: bool = True,
) -> Result[Path, str]:
    """
    Convert SAM to BAM with automatic tool selection.
    
    Uses compressbam for faster compression on Linux,
    falls back to samtools on other platforms or if compressbam is unavailable.
    
    Args:
        sam_path: Input SAM file
        bam_path: Output BAM file (default: sam_path with .bam extension)
        threads: Number of threads
        use_compressbam: Whether to try compressbam first (default: True)
        remove_sam: Whether to remove SAM file after conversion (default: True)
        
    Returns:
        Result containing path to output BAM file
    """
    sam_path = Path(sam_path)
    
    if bam_path is None:
        bam_path = sam_path.with_suffix(".bam")
    else:
        bam_path = Path(bam_path)
    
    # Check if SAM has alignments (skip all headers first)
    # Use awk to skip header lines and print first non-header line
    check_cmd = f"awk '/^[^@]/{{print; exit}}' {sam_path}"
    result = subprocess.run(check_cmd, shell=True, capture_output=True, text=True)
    
    if not result.stdout.strip():
        logger.warning(f"No alignments found in {sam_path}")
        return Err("No alignments found in SAM file")
    
    # Try compressbam first (Linux only)
    if use_compressbam and get_compressbam_path() is not None:
        logger.debug(f"Using compressbam for SAM->BAM conversion")
        result = compress_sam_with_compressbam(sam_path, bam_path, threads)
        
        if result.is_ok():
            if remove_sam:
                sam_path.unlink()
            return result
        else:
            logger.warning(f"compressbam failed, falling back to samtools: {result.unwrap_err()}")
    
    # Fallback to samtools
    logger.debug(f"Using samtools for SAM->BAM conversion")
    result = compress_sam_with_samtools(sam_path, bam_path, threads)
    
    if result.is_ok() and remove_sam:
        sam_path.unlink()
    
    return result


def _compress_worker(
    bam_path: Path,
    lib_name: str,
    out_base: Path,
    threads_each: int,
) -> Optional[Path]:
    """
    Worker function for parallel BAM compression.
    
    Args:
        bam_path: Input BAM file
        lib_name: Library name for output file naming
        out_base: Base path for output
        threads_each: Threads per compression job
        
    Returns:
        Output path or None if empty
    """
    out = Path(f"{out_base}.{lib_name}.bgzf.bam")
    
    compressbam = get_compressbam_path()
    if compressbam:
        cmd = [
            str(compressbam),
            "--threads", str(threads_each),
            "--input", str(bam_path),
            "--output", str(out),
        ]
    else:
        # Fallback: just copy with samtools (already compressed)
        cmd = [
            "samtools", "view",
            "-@", str(threads_each),
            "-b",
            "-o", str(out),
            str(bam_path),
        ]
    
    subprocess.run(cmd, check=True, capture_output=True)
    
    if not out.exists() or out.stat().st_size == 0:
        if out.exists():
            out.unlink()
        return None
    
    return out


def parallel_compress_bams(
    bam_files: List[Tuple[Path, str]],
    output_base: Path,
    total_threads: int = 1,
    max_threads_per_job: int = 8,
) -> Result[List[Path], str]:
    """
    Compress multiple BAM files in parallel.
    
    Distributes work across multiple processes, with each process
    using multiple threads for compression.
    
    Args:
        bam_files: List of (bam_path, library_name) tuples
        output_base: Base path for output files
        total_threads: Total available threads
        max_threads_per_job: Maximum threads per compression job
        
    Returns:
        Result containing list of compressed BAM paths
    """
    if not bam_files:
        return Ok([])
    
    # Calculate parallelism
    threads_each = min(max_threads_per_job, total_threads)
    workers = max(1, total_threads // threads_each)
    workers = min(workers, len(bam_files))
    
    # Adjust threads_each based on actual worker count
    threads_each = min(max_threads_per_job, max(1, total_threads // workers))
    
    logger.info(f"Parallel BAM compression: {len(bam_files)} files")
    logger.info(f"  Workers: {workers}, Threads/worker: {threads_each}")
    
    compressed: List[Path] = []
    
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_compress_worker, bam, name, output_base, threads_each): name
            for bam, name in bam_files
        }
        
        for future in as_completed(futures):
            lib = futures[future]
            try:
                path = future.result()
                if path:
                    compressed.append(path)
                    logger.debug(f"  Compressed: {lib}")
                else:
                    logger.warning(f"  Empty output: {lib}")
            except Exception as e:
                logger.error(f"  Failed: {lib} - {e}")
    
    return Ok(compressed)


def merge_and_namesort_bams(
    bam_files: List[Path],
    output_path: Path,
    threads: int = 1,
    cleanup: bool = True,
) -> Result[Path, str]:
    """
    Merge multiple BAMs and sort by read name (required for ngsLCA).
    
    Args:
        bam_files: List of BAM files to merge
        output_path: Output sorted BAM path
        threads: Number of threads
        cleanup: Whether to remove input BAMs after merging
        
    Returns:
        Result containing path to sorted BAM file
    """
    if len(bam_files) == 0:
        return Err("No BAM files to merge")
    
    temp_merged = output_path.with_suffix(".merged.tmp.bam")
    
    try:
        if len(bam_files) == 1:
            # Single file, just sort
            sort_cmd = [
                "samtools", "sort",
                "-n",  # Sort by read name (required for ngsLCA)
                "-@", str(threads),
                "-O", "bam",
                "-o", str(output_path),
                str(bam_files[0]),
            ]
            subprocess.run(sort_cmd, capture_output=True, text=True, check=True)
        else:
            # Multiple files: merge then sort
            merge_cmd = [
                "samtools", "merge",
                "-@", str(threads),
                "-o", str(temp_merged),
            ] + [str(f) for f in bam_files]
            
            subprocess.run(merge_cmd, capture_output=True, text=True, check=True)
            
            # Sort merged file
            sort_cmd = [
                "samtools", "sort",
                "-n",
                "-@", str(threads),
                "-O", "bam",
                "-o", str(output_path),
                str(temp_merged),
            ]
            subprocess.run(sort_cmd, capture_output=True, text=True, check=True)
            
            # Clean up temp merged file
            if temp_merged.exists():
                temp_merged.unlink()
        
        if cleanup:
            for bam in bam_files:
                if bam.exists():
                    bam.unlink()
        
        return Ok(output_path)
        
    except subprocess.CalledProcessError as e:
        return Err(f"BAM merge/sort failed: {e.stderr}")
    except FileNotFoundError:
        return Err("samtools not found. Please install samtools.")
