"""
Probe sampling: random subsampling and CD-HIT clustering.

Modes:
  random:  Randomly sample n probes
  cluster: CD-HIT clustering → family report → representative selection
           with optional family filtering by member count or proportion
"""

import logging
import random
import subprocess
from pathlib import Path
from typing import Dict, Any, Optional, List
from collections import OrderedDict

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


def _resolve_output_fasta(output_prefix: Path) -> Path:
    """Resolve FASTA output path from user -o argument.

    If the user provides a FASTA filename (e.g. *.fa, *.fasta, *.fna),
    write to that exact path. Otherwise, keep legacy prefix behavior by
    appending .sampled.fa.
    """
    suffixes = {s.lower() for s in output_prefix.suffixes}
    if {".fa", ".fasta", ".fna"} & suffixes:
        return output_prefix
    return Path(str(output_prefix) + ".sampled.fa")


def random_sample(
    sequences: Dict[str, str],
    n: int,
    seed: int = 42,
) -> Dict[str, str]:
    """Random sample n sequences."""
    random.seed(seed)
    if len(sequences) <= n:
        return sequences
    selected = random.sample(list(sequences.keys()), n)
    return OrderedDict((k, sequences[k]) for k in selected)


def parse_cdhit_clusters(cluster_file: Path) -> Dict[str, List[str]]:
    """
    Parse CD-HIT .clstr file.
    
    Returns:
        {representative_id: [member_ids]} — representative marked with '*'
    """
    clusters: Dict[str, List[str]] = OrderedDict()
    current_members: List[str] = []
    current_rep: Optional[str] = None
    
    with open(cluster_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                if current_rep is not None and current_members:
                    clusters[current_rep] = current_members
                current_members = []
                current_rep = None
            else:
                parts = line.split('>')
                if len(parts) >= 2:
                    seq_id = parts[1].split('...')[0].strip()
                    current_members.append(seq_id)
                    if line.rstrip().endswith('*'):
                        current_rep = seq_id
    
    if current_rep is not None and current_members:
        clusters[current_rep] = current_members
    
    return clusters


def run_cdhit(
    fasta_path: Path,
    output_path: Path,
    identity: float = 0.95,
    coverage: float = 0.9,
    threads: int = 1,
) -> Result[Path, str]:
    """Run CD-HIT-EST clustering for nucleotide sequences."""
    cmd = [
        "cd-hit-est",
        "-i", str(fasta_path),
        "-o", str(output_path),
        "-c", str(identity),
        "-aS", str(coverage),
        "-T", str(threads),
        "-M", "0",
        "-d", "0",
    ]
    
    # Word size based on identity threshold
    if identity >= 0.9:
        cmd.extend(["-n", "8"])
    elif identity >= 0.88:
        cmd.extend(["-n", "7"])
    elif identity >= 0.85:
        cmd.extend(["-n", "6"])
    elif identity >= 0.80:
        cmd.extend(["-n", "5"])
    else:
        cmd.extend(["-n", "4"])
    
    logger.debug(f"Running: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        return Ok(output_path)
    except subprocess.CalledProcessError as e:
        return Err(f"cd-hit-est failed: {e.stderr}")
    except FileNotFoundError:
        return Err("cd-hit-est not found. Install: conda install -c bioconda cd-hit")


def run_sample(
    input_path: Path,
    output_prefix: Path,
    mode: str = "random",
    n: Optional[int] = None,
    seed: int = 42,
    identity: float = 0.95,
    coverage: float = 0.9,
    min_members: Optional[int] = None,
    max_members: Optional[int] = None,
    min_proportion: Optional[float] = None,
    max_proportion: Optional[float] = None,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Sample probes: random or CD-HIT cluster-based.
    
    Modes:
      random:  Randomly sample n probes
      cluster: CD-HIT clustering → family report → representative selection
    
    Cluster mode family filtering:
      --min_members / --max_members: Filter by family member count
      --min_proportion / --max_proportion: Filter by family proportion (0-1)
      
    Args:
        input_path: Input FASTA
        output_prefix: Output prefix
        mode: "random" or "cluster"
        n: Target number of probes (required for random, optional for cluster)
        seed: Random seed
        identity: CD-HIT sequence identity threshold
        coverage: CD-HIT alignment coverage
        min_members: Min family members to keep
        max_members: Max family members to keep
        min_proportion: Min family proportion to keep
        max_proportion: Max family proportion to keep
        threads: Threads for CD-HIT
        verbose: Verbose logging
        
    Returns:
        Result with sampling statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    fasta_result = read_fasta(input_path)
    if fasta_result.is_err():
        return Err(f"Failed to read input: {fasta_result.unwrap_err()}")
    
    sequences = fasta_result.unwrap()
    n_input = len(sequences)
    logger.info(f"Loaded {n_input} probes")
    
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    if mode == "random":
        if n is None:
            return Err("--n is required for random sampling")
        
        sampled = random_sample(sequences, n, seed)
        
        fasta_path = _resolve_output_fasta(output_prefix)
        write_result = write_fasta(sampled, fasta_path)
        if write_result.is_err():
            return Err(f"Write failed: {write_result.unwrap_err()}")
        
        return Ok({
            "input_count": n_input,
            "output_count": len(sampled),
            "mode": "random",
            "fasta_file": str(fasta_path),
        })
    
    elif mode == "cluster":
        # Run CD-HIT
        cdhit_out = Path(str(output_prefix) + ".cdhit")
        cdhit_result = run_cdhit(input_path, cdhit_out, identity, coverage, threads)
        if cdhit_result.is_err():
            return Err(cdhit_result.unwrap_err())
        
        # Parse clusters
        cluster_file = Path(str(cdhit_out) + ".clstr")
        if not cluster_file.exists():
            return Err(f"Cluster file not found: {cluster_file}")
        
        clusters = parse_cdhit_clusters(cluster_file)
        total_families = len(clusters)
        logger.info(f"Found {total_families} cluster families")
        
        # Write full cluster report
        report_path = Path(str(output_prefix) + ".cluster_report.tsv")
        with open(report_path, 'w') as f:
            f.write("family_id\trepresentative\tn_members\tproportion\tmembers\n")
            for fam_idx, (rep, members) in enumerate(clusters.items(), 1):
                prop = len(members) / n_input
                member_str = ','.join(members)
                f.write(f"family_{fam_idx}\t{rep}\t{len(members)}\t{prop:.4f}\t{member_str}\n")
        
        # Apply family filters
        filtered_clusters: Dict[str, List[str]] = OrderedDict()
        for rep, members in clusters.items():
            n_mem = len(members)
            prop = n_mem / n_input
            
            if min_members is not None and n_mem < min_members:
                continue
            if max_members is not None and n_mem > max_members:
                continue
            if min_proportion is not None and prop < min_proportion:
                continue
            if max_proportion is not None and prop > max_proportion:
                continue
            
            filtered_clusters[rep] = members
        
        n_filtered = total_families - len(filtered_clusters)
        if n_filtered > 0:
            logger.info(f"Removed {n_filtered} families by filter, "
                        f"{len(filtered_clusters)} remaining")
        
        # Select representatives
        selected: Dict[str, str] = OrderedDict()
        for rep in filtered_clusters:
            if rep in sequences:
                selected[rep] = sequences[rep]
        
        # Additional random sampling if n is specified
        if n is not None and len(selected) > n:
            selected = random_sample(selected, n, seed)
        
        # Write selected probes
        fasta_path = _resolve_output_fasta(output_prefix)
        write_result = write_fasta(selected, fasta_path)
        if write_result.is_err():
            return Err(f"Write failed: {write_result.unwrap_err()}")
        
        # Write summary
        summary_path = Path(str(output_prefix) + ".sample_summary.txt")
        with open(summary_path, 'w') as f:
            f.write("Cluster Sampling Summary\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Input probes: {n_input}\n")
            f.write(f"CD-HIT identity: {identity}\n")
            f.write(f"CD-HIT coverage: {coverage}\n")
            f.write(f"Total families: {total_families}\n")
            if n_filtered > 0:
                f.write(f"Families removed by filter: {n_filtered}\n")
                if min_members is not None:
                    f.write(f"  min_members: {min_members}\n")
                if max_members is not None:
                    f.write(f"  max_members: {max_members}\n")
                if min_proportion is not None:
                    f.write(f"  min_proportion: {min_proportion}\n")
                if max_proportion is not None:
                    f.write(f"  max_proportion: {max_proportion}\n")
            f.write(f"Families kept: {len(filtered_clusters)}\n")
            f.write(f"Selected representatives: {len(selected)}\n")
        
        return Ok({
            "input_count": n_input,
            "total_families": total_families,
            "filtered_families": len(filtered_clusters),
            "families_removed": n_filtered,
            "output_count": len(selected),
            "mode": "cluster",
            "identity": identity,
            "fasta_file": str(fasta_path),
            "report_file": str(report_path),
            "summary_file": str(summary_path),
        })
    
    else:
        return Err(f"Unknown mode: {mode}. Use 'random' or 'cluster'.")
