"""
Map probes to a genome and generate target BED regions.

Pipeline:
  1. Build bowtie2 index (if not provided)
  2. Map probes with bowtie2 (local alignment)
  3. Convert to BAM → BED
  4. Sort + merge overlapping regions
  5. Output target BED file

Use case: After switching reference genomes, define capture target
regions based on where existing probes would hybridize.
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Any, Optional

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta

logger = logging.getLogger(__name__)


def build_bowtie2_index(
    genome_path: Path,
    index_prefix: Path,
    threads: int = 1,
) -> Result[Path, str]:
    """Build bowtie2 index from genome FASTA."""
    cmd = [
        "bowtie2-build",
        "--threads", str(threads),
        str(genome_path),
        str(index_prefix),
    ]
    
    logger.info(f"Building bowtie2 index: {genome_path}")
    
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return Ok(index_prefix)
    except subprocess.CalledProcessError as e:
        return Err(f"bowtie2-build failed: {e.stderr}")
    except FileNotFoundError:
        return Err("bowtie2-build not found in PATH")


def map_probes(
    probe_fasta: Path,
    index_prefix: Path,
    output_bam: Path,
    threads: int = 1,
    mode: str = "very-sensitive-local",
    min_score: Optional[int] = None,
) -> Result[Path, str]:
    """
    Map probes to genome using bowtie2.
    
    Uses local alignment mode for partial probe matches.
    
    Args:
        probe_fasta: Probe FASTA file
        index_prefix: Bowtie2 index prefix
        output_bam: Output BAM path
        threads: Number of threads
        mode: Bowtie2 preset (default: very-sensitive-local)
        min_score: Minimum alignment score
        
    Returns:
        Result containing BAM path
    """
    cmd = [
        "bowtie2",
        "-x", str(index_prefix),
        "-f", str(probe_fasta),
        f"--{mode}",
        "--threads", str(threads),
        "--no-unal",
        "-k", "1",
    ]
    
    if min_score is not None:
        cmd.extend(["--score-min", f"L,{min_score},0"])
    
    logger.info(f"Mapping probes with bowtie2 ({mode})")
    
    try:
        # Pipe bowtie2 → samtools view → BAM
        bt2_proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        
        sam_cmd = [
            "samtools", "view", "-bS",
            "-@", str(threads),
            "-o", str(output_bam), "-",
        ]
        sam_proc = subprocess.Popen(
            sam_cmd, stdin=bt2_proc.stdout,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        
        bt2_proc.stdout.close()
        _, sam_err = sam_proc.communicate()
        _, bt2_err = bt2_proc.communicate()
        
        if bt2_proc.returncode != 0:
            return Err(f"bowtie2 failed: {bt2_err.decode()}")
        
        # Sort BAM
        sorted_bam = Path(str(output_bam) + ".sorted.bam")
        sort_cmd = [
            "samtools", "sort",
            "-@", str(threads),
            "-o", str(sorted_bam),
            str(output_bam),
        ]
        subprocess.run(sort_cmd, check=True, capture_output=True)
        sorted_bam.rename(output_bam)
        
        # Index BAM
        subprocess.run(
            ["samtools", "index", str(output_bam)],
            check=True, capture_output=True,
        )
        
        return Ok(output_bam)
        
    except subprocess.CalledProcessError as e:
        return Err(f"Mapping failed: {e}")
    except FileNotFoundError:
        return Err("bowtie2 or samtools not found in PATH")


def bam_to_bed(bam_path: Path, bed_path: Path) -> Result[Path, str]:
    """Convert BAM to BED using bedtools."""
    cmd = ["bedtools", "bamtobed", "-i", str(bam_path)]
    
    try:
        with open(bed_path, 'w') as f:
            subprocess.run(
                cmd, stdout=f, check=True,
                capture_output=False, stderr=subprocess.PIPE,
            )
        return Ok(bed_path)
    except subprocess.CalledProcessError as e:
        return Err(f"bedtools bamtobed failed: {e.stderr.decode()}")
    except FileNotFoundError:
        return Err("bedtools not found in PATH")


def merge_bed(
    bed_path: Path,
    merged_path: Path,
    merge_distance: int = 0,
) -> Result[Path, str]:
    """Sort and merge overlapping BED regions."""
    try:
        sort_cmd = ["sort", "-k1,1", "-k2,2n", str(bed_path)]
        sorted_result = subprocess.run(
            sort_cmd, capture_output=True, text=True, check=True,
        )
        
        merge_cmd = [
            "bedtools", "merge",
            "-d", str(merge_distance),
            "-i", "stdin",
        ]
        merge_result = subprocess.run(
            merge_cmd,
            input=sorted_result.stdout,
            capture_output=True, text=True, check=True,
        )
        
        with open(merged_path, 'w') as f:
            f.write(merge_result.stdout)
        
        return Ok(merged_path)
        
    except subprocess.CalledProcessError as e:
        return Err(f"BED merge failed: {e.stderr}")
    except FileNotFoundError:
        return Err("bedtools not found in PATH")


def run_target(
    probe_fasta: Path,
    genome_path: Path,
    output_prefix: Path,
    index_prefix: Optional[Path] = None,
    threads: int = 1,
    mode: str = "very-sensitive-local",
    merge_distance: int = 0,
    min_score: Optional[int] = None,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Map probes to genome and generate target BED regions.
    
    Pipeline:
      1. Build bowtie2 index (if not provided)
      2. Map probes with bowtie2 (local alignment)
      3. Convert BAM → BED
      4. Sort + merge overlapping regions
      5. Output target BED
      
    Args:
        probe_fasta: Input probe FASTA
        genome_path: Reference genome FASTA
        output_prefix: Output prefix
        index_prefix: Pre-built bowtie2 index (optional)
        threads: Number of threads
        mode: Bowtie2 alignment preset
        merge_distance: Max distance for merging regions
        min_score: Minimum alignment score
        verbose: Verbose logging
        
    Returns:
        Result with target region statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    # Count input probes
    fasta_r = read_fasta(probe_fasta)
    n_probes = len(fasta_r.unwrap()) if fasta_r.is_ok() else 0
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        
        # Step 1: Build or use existing index
        if index_prefix is None:
            idx = tmpdir_path / "genome_index"
            idx_result = build_bowtie2_index(genome_path, idx, threads)
            if idx_result.is_err():
                return Err(idx_result.unwrap_err())
            index_prefix = idx
        
        # Step 2-3: Map probes
        bam_path = tmpdir_path / "mapped.bam"
        map_result = map_probes(
            probe_fasta, index_prefix, bam_path,
            threads, mode, min_score,
        )
        if map_result.is_err():
            return Err(map_result.unwrap_err())
        
        # Count mapped reads
        try:
            count_cmd = ["samtools", "view", "-c", str(bam_path)]
            count_result = subprocess.run(
                count_cmd, capture_output=True, text=True, check=True,
            )
            n_mapped = int(count_result.stdout.strip())
        except Exception:
            n_mapped = 0
        
        # Step 4: BAM → BED
        raw_bed = tmpdir_path / "raw.bed"
        bed_result = bam_to_bed(bam_path, raw_bed)
        if bed_result.is_err():
            return Err(bed_result.unwrap_err())
        
        # Step 5: Merge overlapping regions
        merged_bed = Path(str(output_prefix) + ".target.bed")
        merge_result = merge_bed(raw_bed, merged_bed, merge_distance)
        if merge_result.is_err():
            return Err(merge_result.unwrap_err())
    
    # Count merged regions
    n_regions = 0
    total_bp = 0
    with open(merged_bed) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                n_regions += 1
                total_bp += int(parts[2]) - int(parts[1])
    
    # Write summary
    summary_path = Path(str(output_prefix) + ".target_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("Target Region Generation Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Probe FASTA: {probe_fasta}\n")
        f.write(f"Genome: {genome_path}\n")
        f.write(f"Alignment mode: {mode}\n")
        f.write(f"Merge distance: {merge_distance}bp\n\n")
        f.write(f"Input probes: {n_probes}\n")
        f.write(f"Mapped probes: {n_mapped}\n")
        if n_probes > 0:
            f.write(f"Mapping rate: {n_mapped / n_probes * 100:.1f}%\n")
        f.write(f"\nTarget regions: {n_regions}\n")
        f.write(f"Total target bp: {total_bp:,}\n")
    
    stats = {
        "n_probes": n_probes,
        "n_mapped": n_mapped,
        "mapping_rate": round(n_mapped / n_probes * 100, 1) if n_probes > 0 else 0,
        "n_regions": n_regions,
        "total_bp": total_bp,
        "bed_file": str(merged_bed),
        "summary_file": str(summary_path),
    }
    
    return Ok(stats)
