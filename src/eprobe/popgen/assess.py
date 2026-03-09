"""
Probe quality assessment module.

Calculates and reports quality metrics for probe sets:
  - Biophysical properties (GC, Tm, complexity, hairpin, dimer)
  - Genomic coverage statistics
  - IBS distance comparison between full VCF and probe-covered SNPs
  - Site Frequency Spectrum comparison using easySFS
  - Distribution plots

This module corresponds to the original SNP_assessor.py functionality.
"""

import logging
import re
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict

import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

from eprobe.core.result import Result, Ok, Err
from eprobe.core.models import SNPDataFrame, ProbeSet
from eprobe.core.fasta import read_fasta
# Use fast biophysical calculators (same as filter.py)
from eprobe.biophysics.biophysics import (
    calculate_gc_fast,
    calculate_tm_fast,
    calculate_dust_fast,
    calculate_hairpin_fast,
    DimerCalculatorFast,
)
from eprobe.biophysics.entropy import calculate_entropy

logger = logging.getLogger(__name__)


# Available assessment tags
AVAILABLE_TAGS = {
    "gc": "GC content (%)",
    "tm": "Melting temperature (°C)",
    "complexity": "DUST complexity score",
    "hairpin": "Self-complementarity score",
    "dimer": "Inter-probe complementarity score",
    "entropy": "Shannon entropy",
}


@dataclass
class AssessmentResult:
    """Container for assessment results."""
    probe_id: str
    sequence: str
    gc: Optional[float] = None
    tm: Optional[float] = None
    complexity: Optional[float] = None
    hairpin: Optional[float] = None
    dimer: Optional[float] = None
    entropy: Optional[float] = None


def assess_single_probe(
    probe_id: str,
    sequence: str,
    tags: List[str],
) -> AssessmentResult:
    """
    Calculate assessment metrics for a single probe.
    
    Uses the same fast biophysics functions as the filter step.
    
    Args:
        probe_id: Probe identifier
        sequence: Probe sequence
        tags: List of metrics to calculate
        
    Returns:
        AssessmentResult with calculated values
    """
    result = AssessmentResult(probe_id=probe_id, sequence=sequence)
    
    if "gc" in tags:
        result.gc = calculate_gc_fast(sequence)
    
    if "tm" in tags:
        result.tm = calculate_tm_fast(sequence)
    
    if "complexity" in tags:
        result.complexity = calculate_dust_fast(sequence)
    
    if "hairpin" in tags:
        result.hairpin = calculate_hairpin_fast(sequence)
    
    if "entropy" in tags:
        result.entropy = calculate_entropy(sequence)
    
    return result


def assess_probes(
    sequences: Dict[str, str],
    tags: List[str],
) -> List[AssessmentResult]:
    """
    Assess all probes for specified metrics.
    
    Args:
        sequences: Dictionary of probe_id -> sequence
        tags: List of metrics to calculate
        
    Returns:
        List of AssessmentResult objects
    """
    results = []
    
    for probe_id, sequence in sequences.items():
        result = assess_single_probe(probe_id, sequence, tags)
        results.append(result)
    
    return results


def calculate_dimer_scores(
    sequences: Dict[str, str],
    sample_size: int = 1000,
) -> Dict[str, float]:
    """
    Calculate dimer scores for probe sequences.
    
    Uses the same DimerCalculatorFast as the filter step.
    Higher scores indicate more shared k-mers with other probes.
    
    Args:
        sequences: Dictionary of probe_id -> sequence
        sample_size: Not used (kept for API compatibility)
        
    Returns:
        Dictionary of probe_id -> dimer score
    """
    if len(sequences) <= 1:
        return {pid: 0.0 for pid in sequences.keys()}
    
    # Use same dimer calculator as filter.py
    dimer_calc = DimerCalculatorFast(k=11)
    
    # Build k-mer index from all probes
    seq_list = list(sequences.values())
    n_kmers = dimer_calc.build_index(seq_list)
    
    if n_kmers == 0:
        logger.warning("No k-mers found for dimer calculation")
        return {pid: 0.0 for pid in sequences.keys()}
    
    # Calculate scores for all probes
    scores = dimer_calc.calculate_all_scores()
    
    # Map scores back to probe IDs
    probe_ids = list(sequences.keys())
    return {pid: scores[i] for i, pid in enumerate(probe_ids)}


def calculate_summary_statistics(
    results: List[AssessmentResult],
    tags: List[str],
) -> Dict[str, Dict[str, float]]:
    """
    Calculate summary statistics for each metric.
    
    Args:
        results: List of assessment results
        tags: List of metrics to summarize
        
    Returns:
        Dictionary of metric -> {mean, std, min, max, median}
    """
    summary = {}
    
    for tag in tags:
        values = [getattr(r, tag) for r in results if getattr(r, tag) is not None]
        
        if not values:
            continue
        
        summary[tag] = {
            "mean": float(np.mean(values)),
            "std": float(np.std(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "median": float(np.median(values)),
            "count": len(values),
        }
    
    return summary


# =============================================================================
# IBS Distance Analysis Functions
# =============================================================================

def parse_vcf_genotypes(
    vcf_path: Path,
    positions: Optional[Set[Tuple[str, int]]] = None,
    max_sites: Optional[int] = None,
) -> Result[Tuple[List[str], np.ndarray, List[Tuple[str, int]]], str]:
    """
    Parse VCF file and extract genotype matrix.
    
    Args:
        vcf_path: Path to VCF file (can be .vcf or .vcf.gz)
        positions: Optional set of (chrom, pos) to filter to
        max_sites: Maximum number of sites to load (for memory)
        
    Returns:
        Result containing (sample_names, genotype_matrix, site_positions)
        Genotype matrix: 0=homref, 1=het, 2=homalt, -1=missing
    """
    import gzip
    
    # Determine if compressed
    open_func = gzip.open if str(vcf_path).endswith('.gz') else open
    mode = 'rt' if str(vcf_path).endswith('.gz') else 'r'
    
    samples = []
    genotypes = []
    site_pos = []
    
    try:
        with open_func(vcf_path, mode) as f:
            for line in f:
                if line.startswith('##'):
                    continue
                
                if line.startswith('#CHROM'):
                    # Header line with sample names
                    parts = line.strip().split('\t')
                    samples = parts[9:]  # Sample names start at column 9
                    continue
                
                # Data line
                parts = line.strip().split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                
                # Filter by positions if specified
                if positions is not None and (chrom, pos) not in positions:
                    continue
                
                # Parse genotypes
                gt_idx = None
                format_fields = parts[8].split(':')
                for i, field in enumerate(format_fields):
                    if field == 'GT':
                        gt_idx = i
                        break
                
                if gt_idx is None:
                    continue
                
                site_gts = []
                for sample_data in parts[9:]:
                    gt_field = sample_data.split(':')[gt_idx] if ':' in sample_data else sample_data
                    
                    # Parse genotype
                    if gt_field in ['.', './.', '.|.']:
                        site_gts.append(-1)  # Missing
                    elif '/' in gt_field:
                        alleles = gt_field.split('/')
                        try:
                            a1, a2 = int(alleles[0]), int(alleles[1])
                            site_gts.append(a1 + a2)  # 0+0=0, 0+1=1, 1+1=2
                        except ValueError:
                            site_gts.append(-1)
                    elif '|' in gt_field:
                        alleles = gt_field.split('|')
                        try:
                            a1, a2 = int(alleles[0]), int(alleles[1])
                            site_gts.append(a1 + a2)
                        except ValueError:
                            site_gts.append(-1)
                    else:
                        site_gts.append(-1)
                
                genotypes.append(site_gts)
                site_pos.append((chrom, pos))
                
                # Check max sites
                if max_sites and len(genotypes) >= max_sites:
                    break
        
        if not samples:
            return Err("No samples found in VCF")
        
        if not genotypes:
            return Err("No valid genotypes found in VCF")
        
        gt_matrix = np.array(genotypes, dtype=np.int8).T  # Shape: (n_samples, n_sites)
        
        return Ok((samples, gt_matrix, site_pos))
        
    except Exception as e:
        return Err(f"Error parsing VCF: {str(e)}")


def calculate_ibs_distance_matrix(
    genotype_matrix: np.ndarray,
) -> np.ndarray:
    """
    Calculate pairwise 1-IBS distance matrix.
    
    1-IBS (one minus Identity By State) measures genetic distance:
    Distance = 1 - IBS_similarity, where IBS_similarity is the proportion
    of alleles shared identical by state.
    
    - IBS2 (same genotype): distance = 0
    - IBS1 (one shared allele): distance = 0.5  
    - IBS0 (no shared alleles): distance = 1
    
    Args:
        genotype_matrix: Shape (n_samples, n_sites), values 0/1/2/-1
        
    Returns:
        Distance matrix of shape (n_samples, n_samples), values in [0, 1]
    """
    n_samples, n_sites = genotype_matrix.shape
    dist_matrix = np.zeros((n_samples, n_samples))
    
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            # Get genotypes for this pair
            gt_i = genotype_matrix[i]
            gt_j = genotype_matrix[j]
            
            # Mask missing data
            valid = (gt_i >= 0) & (gt_j >= 0)
            n_valid = np.sum(valid)
            
            if n_valid == 0:
                dist_matrix[i, j] = np.nan
                dist_matrix[j, i] = np.nan
                continue
            
            gt_i_valid = gt_i[valid]
            gt_j_valid = gt_j[valid]
            
            # Calculate 1-IBS distance
            diff = np.abs(gt_i_valid - gt_j_valid)
            # diff=0 -> IBS2 (1-IBS=0), diff=1 -> IBS1 (1-IBS=0.5), diff=2 -> IBS0 (1-IBS=1)
            ibs_dist = np.sum(diff) / (2 * n_valid)
            
            dist_matrix[i, j] = ibs_dist
            dist_matrix[j, i] = ibs_dist
    
    return dist_matrix


def calculate_distance_correlation(
    dist1: np.ndarray,
    dist2: np.ndarray,
) -> Dict[str, float]:
    """
    Calculate correlation and Manhattan distance between two distance matrices.
    
    Args:
        dist1: First distance matrix
        dist2: Second distance matrix
        
    Returns:
        Dictionary with pearson_r, spearman_r, manhattan_distance
    """
    # Get upper triangle (excluding diagonal)
    n = dist1.shape[0]
    triu_idx = np.triu_indices(n, k=1)
    
    vec1 = dist1[triu_idx]
    vec2 = dist2[triu_idx]
    
    # Remove NaN pairs
    valid = ~(np.isnan(vec1) | np.isnan(vec2))
    vec1_valid = vec1[valid]
    vec2_valid = vec2[valid]
    
    if len(vec1_valid) < 3:
        return {
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "spearman_r": np.nan,
            "spearman_p": np.nan,
            "manhattan_distance": np.nan,
            "n_pairs": len(vec1_valid),
        }
    
    # Pearson correlation
    pearson_r, pearson_p = stats.pearsonr(vec1_valid, vec2_valid)
    
    # Spearman correlation
    spearman_r, spearman_p = stats.spearmanr(vec1_valid, vec2_valid)
    
    # Manhattan distance (normalized)
    manhattan = np.sum(np.abs(vec1_valid - vec2_valid)) / len(vec1_valid)
    
    return {
        "pearson_r": float(pearson_r),
        "pearson_p": float(pearson_p),
        "spearman_r": float(spearman_r),
        "spearman_p": float(spearman_p),
        "manhattan_distance": float(manhattan),
        "n_pairs": len(vec1_valid),
    }


def generate_combined_heatmap(
    dist_full: np.ndarray,
    dist_probe: np.ndarray,
    sample_names: List[str],
    output_path: Path,
    correlation_stats: Dict[str, float],
    pop_file: Optional[Path] = None,
) -> Result[Path, str]:
    """
    Generate combined heatmap with full VCF (upper triangle) and probe set (lower triangle).
    
    Clustering is performed on probe set distances, and the same sample order
    is applied to both triangles for fair comparison.
    If pop_file is provided, a population color strip is drawn to the right of the heatmap.
    
    Args:
        dist_full: Distance matrix from full VCF
        dist_probe: Distance matrix from probe-covered SNPs
        sample_names: List of sample names
        output_path: Output file path
        correlation_stats: Correlation statistics to display
        pop_file: Optional population file for color strip annotation
        
    Returns:
        Result containing output path
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
        from mpl_toolkits.axes_grid1 import make_axes_locatable
    except ImportError:
        return Err("Plotting requires matplotlib. Install with: pip install matplotlib")
    
    n = dist_full.shape[0]
    
    # Hierarchical clustering on PROBE SET only
    dist_probe_clean = np.nan_to_num(dist_probe, nan=np.nanmean(dist_probe))
    np.fill_diagonal(dist_probe_clean, 0)
    
    try:
        condensed = squareform(dist_probe_clean)
        linkage_matrix = linkage(condensed, method='average')
        dendro = dendrogram(linkage_matrix, no_plot=True)
        order = dendro['leaves']
    except Exception:
        order = list(range(n))
    
    # Reorder both matrices using probe set clustering order
    dist_probe_ordered = dist_probe[np.ix_(order, order)]
    dist_full_ordered = dist_full[np.ix_(order, order)]
    sample_names_ordered = [sample_names[i] for i in order]
    
    # Combined matrix: lower triangle = probe, upper = full, diagonal = 1 (gap)
    combined = np.tril(dist_probe_ordered, -1) + np.triu(dist_full_ordered, 1) + np.eye(n)
    
    # --- Build population color strip if pop_file provided ---
    pop_color_strip = None
    pop_legend_handles = []
    if pop_file is not None and pop_file.exists():
        pop_dict_sample = {}
        try:
            with open(pop_file, 'r') as _pf:
                for _line in _pf:
                    _line = _line.strip()
                    if not _line or _line.startswith('#'):
                        continue
                    _parts = _line.split('\t') if '\t' in _line else _line.split()
                    if len(_parts) >= 2:
                        pop_dict_sample[_parts[0]] = _parts[1]
        except Exception as _e:
            logger.warning(f"Could not load pop file for heatmap coloring: {_e}")
        
        if pop_dict_sample:
            pops_ordered = [pop_dict_sample.get(s, 'Unknown') for s in sample_names_ordered]
            unique_pops = sorted(set(pops_ordered))
            cmap_pop = plt.get_cmap('tab10')
            pop_color_map = {p: cmap_pop(i % 10) for i, p in enumerate(unique_pops)}
            # Build (n, 1, 4) RGBA array
            pop_color_strip = np.array(
                [pop_color_map[p] for p in pops_ordered], dtype=float
            ).reshape(n, 1, 4)
            pop_legend_handles = [
                Patch(facecolor=pop_color_map[p], label=p) for p in unique_pops
            ]
    
    # --- Figure and axes ---
    fig, ax = plt.subplots(figsize=(9, 8))
    
    cmap = plt.cm.YlOrRd
    im = ax.imshow(combined, cmap=cmap, vmin=0, vmax=1.0, aspect='auto')
    
    # Build axes attached to the heatmap using make_axes_locatable
    divider = make_axes_locatable(ax)
    
    if pop_color_strip is not None:
        # Population color strip to the right of the heatmap
        ax_pop = divider.append_axes("right", size="2%", pad=0.06)
        ax_pop.imshow(pop_color_strip, aspect='auto', interpolation='nearest')
        ax_pop.set_xticks([])
        ax_pop.set_yticks([])
        for spine in ax_pop.spines.values():
            spine.set_linewidth(0.5)
        ax_pop.set_xlabel('Pop', fontsize=7, labelpad=2)
        cax = divider.append_axes("right", size="3%", pad=0.35)
    else:
        cax = divider.append_axes("right", size="3%", pad=0.1)
    
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label('1-IBS', fontsize=11)
    
    # Sample labels
    if n <= 50:
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(sample_names_ordered, rotation=90, fontsize=8)
        ax.set_yticklabels(sample_names_ordered, fontsize=8)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    
    # Title
    title = "1-IBS Distance: Upper triangle Full VCF | Lower triangle Probe Set\n"
    title += (
        f"Pearson = {correlation_stats['pearson_r']:.3f}, "
        f"Spearman = {correlation_stats['spearman_r']:.3f}, "
        f"Manhattan = {correlation_stats['manhattan_distance']:.3f}"
    )
    ax.set_title(title, fontsize=12, fontweight='bold')
    
    # Population legend below the figure
    if pop_legend_handles:
        ncol = min(len(pop_legend_handles), 6)
        fig.legend(
            handles=pop_legend_handles,
            title='Population',
            loc='lower center',
            bbox_to_anchor=(0.45, -0.02),
            ncol=ncol,
            fontsize=8,
            frameon=True,
        )
    
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    return Ok(output_path)


def generate_scatter_comparison(
    dist_full: np.ndarray,
    dist_probe: np.ndarray,
    output_path: Path,
    correlation_stats: Dict[str, float],
) -> Result[Path, str]:
    """
    Generate scatter plot comparing pairwise distances.
    
    Args:
        dist_full: Distance matrix from full VCF
        dist_probe: Distance matrix from probe-covered SNPs
        output_path: Output file path
        correlation_stats: Correlation statistics
        
    Returns:
        Result containing output path
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return Err("Plotting requires matplotlib")
    
    n = dist_full.shape[0]
    triu_idx = np.triu_indices(n, k=1)
    
    vec_full = dist_full[triu_idx]
    vec_probe = dist_probe[triu_idx]
    
    # Remove NaN
    valid = ~(np.isnan(vec_full) | np.isnan(vec_probe))
    vec_full_valid = vec_full[valid]
    vec_probe_valid = vec_probe[valid]
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    ax.scatter(vec_full_valid, vec_probe_valid, alpha=0.5, s=10, c='steelblue')
    
    # Add diagonal line (perfect correlation)
    lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1]),
    ]
    ax.plot(lims, lims, 'r--', alpha=0.75, zorder=0, label='y = x')
    
    # Add regression line
    if len(vec_full_valid) > 2:
        slope, intercept, _, _, _ = stats.linregress(vec_full_valid, vec_probe_valid)
        x_line = np.array(lims)
        y_line = slope * x_line + intercept
        ax.plot(x_line, y_line, 'g-', alpha=0.75, label=f'Regression (slope={slope:.3f})')
    
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect('equal')
    
    ax.set_xlabel('1-IBS (Full VCF)', fontsize=12)
    ax.set_ylabel('1-IBS (Probe Set)', fontsize=12)
    
    title = f"Pairwise 1-IBS Distance Comparison\n"
    title += f"r = {correlation_stats['pearson_r']:.4f}, ρ = {correlation_stats['spearman_r']:.4f}\n"
    title += f"n = {correlation_stats['n_pairs']} pairs"
    ax.set_title(title, fontsize=11)
    
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    return Ok(output_path)


# =============================================================================
# PLINK-based Analysis Functions (IBS Distance, PCA)
# =============================================================================

def check_plink_available() -> Result[str, str]:
    """
    Check if plink is available in PATH.
    
    Returns:
        Result containing path to plink executable
    """
    import shutil
    
    plink_cmd = shutil.which("plink")
    if plink_cmd is None:
        return Err("plink not found in PATH. Install with: conda install -c bioconda plink")
    return Ok(plink_cmd)


def run_plink_ibs_distance(
    vcf_path: Path,
    output_prefix: Path,
    threads: int = 1,
) -> Result[Tuple[np.ndarray, List[str]], str]:
    """
    Calculate 1-IBS distance matrix using PLINK.
    
    Args:
        vcf_path: Path to VCF file
        output_prefix: Output prefix for plink files
        
    Returns:
        Result containing (distance_matrix, sample_names)
    """
    import subprocess
    import tempfile
    
    plink_check = check_plink_available()
    if plink_check.is_err():
        return Err(plink_check.unwrap_err())
    
    plink_cmd = plink_check.unwrap()
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        tmp_prefix = tmpdir / "plink_ibs"
        
        # Run plink to calculate 1-IBS distance
        # --double-id: treats sample IDs as both FID and IID (avoids ID parsing issues)
        cmd = [
            plink_cmd,
            "--allow-extra-chr",
            "--double-id",
            "--vcf", str(vcf_path),
            "--distance", "square", "1-ibs",
            "--threads", str(threads),
            "--out", str(tmp_prefix)
        ]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
            )
            
            if result.returncode != 0:
                return Err(f"PLINK failed: {result.stderr}")
            
            # Read distance matrix - plink outputs .mdist and .mdist.id for 1-ibs
            dist_file = tmp_prefix.with_suffix(".mdist")
            id_file = Path(str(tmp_prefix) + ".mdist.id")
            
            if not dist_file.exists():
                return Err(f"PLINK distance file not found: {dist_file}")
            
            # Read sample IDs
            sample_ids = []
            with open(id_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        sample_ids.append(parts[1])  # IID is second column
                    elif len(parts) == 1:
                        sample_ids.append(parts[0])
            
            # Read distance matrix
            dist_matrix = np.loadtxt(dist_file)
            
            return Ok((dist_matrix, sample_ids))
            
        except Exception as e:
            return Err(f"Failed to run PLINK IBS: {e}")


def run_plink_pca(
    vcf_path: Path,
    output_prefix: Path,
    n_pcs: int = 10,
    threads: int = 1,
) -> Result[Tuple[pd.DataFrame, Dict[str, float], List[str]], str]:
    """
    Run PCA using PLINK.
    
    Args:
        vcf_path: Path to VCF file
        output_prefix: Output prefix for plink files
        n_pcs: Number of principal components to calculate
        
    Returns:
        Result containing (eigenvectors_df, variance_explained_dict, sample_names)
    """
    import subprocess
    import tempfile
    
    plink_check = check_plink_available()
    if plink_check.is_err():
        return Err(plink_check.unwrap_err())
    
    plink_cmd = plink_check.unwrap()
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        bed_prefix = tmpdir / "plink_bed"
        pca_prefix = tmpdir / "plink_pca"
        
        # Step 1: Convert VCF to bed format
        # --double-id: treats sample IDs as both FID and IID (avoids ID parsing issues)
        cmd_bed = [
            plink_cmd,
            "--allow-extra-chr",
            "--double-id",
            "--vcf", str(vcf_path),
            "--make-bed",
            "--threads", str(threads),
            "--out", str(bed_prefix)
        ]
        
        # Step 2: Run PCA
        cmd_pca = [
            plink_cmd,
            "--allow-extra-chr",
            "--bfile", str(bed_prefix),
            "--pca", str(n_pcs),
            "--threads", str(threads),
            "--out", str(pca_prefix)
        ]
        
        try:
            # Run bed conversion
            result = subprocess.run(cmd_bed, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                return Err(f"PLINK bed conversion failed: {result.stderr}")
            
            # Run PCA
            result = subprocess.run(cmd_pca, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                return Err(f"PLINK PCA failed: {result.stderr}")
            
            # Read eigenvectors
            eigenvec_file = pca_prefix.with_suffix(".eigenvec")
            eigenval_file = pca_prefix.with_suffix(".eigenval")
            
            if not eigenvec_file.exists():
                return Err(f"PLINK eigenvec file not found: {eigenvec_file}")
            
            # Parse eigenvectors
            eigenvec_df = pd.read_csv(eigenvec_file, sep=r'\s+', header=None)
            sample_ids = eigenvec_df.iloc[:, 1].tolist()  # IID is second column
            
            # Set column names
            n_cols = eigenvec_df.shape[1]
            col_names = ['FID', 'IID'] + [f'PC{i}' for i in range(1, n_cols - 1)]
            eigenvec_df.columns = col_names
            eigenvec_df = eigenvec_df.drop(columns=['FID'])
            eigenvec_df = eigenvec_df.rename(columns={'IID': 'ID'})
            
            # Parse eigenvalues and calculate variance explained
            variance_dict = {}
            if eigenval_file.exists():
                eigenvalues = pd.read_csv(eigenval_file, header=None)[0].values
                total_var = np.sum(eigenvalues)
                for i, val in enumerate(eigenvalues):
                    variance_dict[f'PC{i+1}'] = (val / total_var) * 100
            
            return Ok((eigenvec_df, variance_dict, sample_ids))
            
        except Exception as e:
            return Err(f"Failed to run PLINK PCA: {e}")


def generate_pca_comparison_plot(
    pca_full: pd.DataFrame,
    pca_probe: pd.DataFrame,
    var_full: Dict[str, float],
    var_probe: Dict[str, float],
    output_path: Path,
    pop_file: Optional[Path] = None,
    pc_x: int = 1,
    pc_y: int = 2,
) -> Result[Path, str]:
    """
    Generate side-by-side PCA comparison plots.
    
    Args:
        pca_full: Eigenvectors from full VCF
        pca_probe: Eigenvectors from probe VCF
        var_full: Variance explained for full VCF
        var_probe: Variance explained for probe VCF
        output_path: Output file path
        pop_file: Optional population file for coloring
        pc_x: X-axis PC number
        pc_y: Y-axis PC number
        
    Returns:
        Result containing output path
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return Err("Plotting requires matplotlib")
    
    # Load population info if provided
    pop_dict = {}
    if pop_file is not None and pop_file.exists():
        try:
            with open(pop_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t') if '\t' in line else line.split()
                    if len(parts) >= 2:
                        pop_dict[parts[0]] = parts[1]
        except Exception as e:
            logger.warning(f"Could not load pop file: {e}")
    
    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    pc_x_col = f'PC{pc_x}'
    pc_y_col = f'PC{pc_y}'
    
    # Plot function
    def plot_pca(ax, pca_df, var_dict, title):
        if pop_dict:
            # Add population column
            pca_df = pca_df.copy()
            pca_df['Population'] = pca_df['ID'].map(pop_dict)
            populations = pca_df['Population'].dropna().unique()
            
            # Color map
            cmap = plt.get_cmap('tab10')
            colors = {pop: cmap(i % 10) for i, pop in enumerate(sorted(populations))}
            
            for pop in sorted(populations):
                subset = pca_df[pca_df['Population'] == pop]
                ax.scatter(
                    subset[pc_x_col], subset[pc_y_col],
                    c=[colors[pop]], label=pop,
                    alpha=0.8, s=50, edgecolors='black', linewidths=0.5
                )
            
            # Handle samples without population
            no_pop = pca_df[pca_df['Population'].isna()]
            if len(no_pop) > 0:
                ax.scatter(
                    no_pop[pc_x_col], no_pop[pc_y_col],
                    c='gray', label='Unknown',
                    alpha=0.5, s=50, edgecolors='black', linewidths=0.5
                )
            
            ax.legend(title='Population', loc='best', fontsize=9, frameon=True)
        else:
            ax.scatter(
                pca_df[pc_x_col], pca_df[pc_y_col],
                c='steelblue', alpha=0.8, s=50, edgecolors='black', linewidths=0.5
            )
        
        # Labels
        var_x = var_dict.get(pc_x_col, 0)
        var_y = var_dict.get(pc_y_col, 0)
        ax.set_xlabel(f'{pc_x_col} ({var_x:.1f}%)', fontsize=12)
        ax.set_ylabel(f'{pc_y_col} ({var_y:.1f}%)', fontsize=12)
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.tick_params(labelsize=10)
    
    # Plot both
    plot_pca(axes[0], pca_full, var_full, 'Full VCF')
    plot_pca(axes[1], pca_probe, var_probe, 'Probe Set')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    return Ok(output_path)


def run_pca_assessment(
    snp_tsv_path: Path,
    vcf_path: Path,
    output_prefix: Path,
    pop_file: Optional[Path] = None,
    n_pcs: int = 10,
    n_samples: int = 100,
    seed: int = 42,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Compare PCA between full VCF and probe-covered SNPs using PLINK.

    Args:
        snp_tsv_path: Path to selected SNPs TSV file
        vcf_path: Path to original VCF file
        output_prefix: Output file prefix
        pop_file: Optional population file. If provided, stratified subsampling
            is used to ensure each population is represented (n_samples per pop).
            Also used for coloring PCA plots.
        n_pcs: Number of principal components
        n_samples: If pop_file is None: max total individuals (random subsample).
            If pop_file is provided: max individuals per population.
        seed: Random seed for subsampling
        threads: Number of threads for PLINK and bcftools
        verbose: Enable verbose logging

    Returns:
        Result containing PCA assessment statistics
    """
    import tempfile
    import shutil
    
    logger.info("Starting PCA assessment using PLINK")
    logger.info(f"SNP TSV: {snp_tsv_path}")
    logger.info(f"VCF: {vcf_path}")
    logger.info(f"Using {threads} thread(s) for PLINK and bcftools")
    
    # Check plink and bcftools
    plink_check = check_plink_available()
    if plink_check.is_err():
        return Err(plink_check.unwrap_err())
    
    if shutil.which("bcftools") is None:
        return Err("bcftools not found in PATH")
    
    # Validate input files
    if not snp_tsv_path.exists():
        return Err(f"SNP TSV file not found: {snp_tsv_path}")
    if not vcf_path.exists():
        return Err(f"VCF file not found: {vcf_path}")
    
    # Create output directory
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    # Load probe positions from TSV (only requires chr and pos columns)
    try:
        snp_df_pandas = pd.read_csv(snp_tsv_path, sep='\t')
        if 'chr' not in snp_df_pandas.columns or 'pos' not in snp_df_pandas.columns:
            return Err(f"SNP TSV must have 'chr' and 'pos' columns. Found: {list(snp_df_pandas.columns)}")
        
        probe_positions = set()
        for _, row in snp_df_pandas.iterrows():
            probe_positions.add((str(row['chr']), int(row['pos'])))
    except Exception as e:
        return Err(f"Failed to load SNP TSV: {e}")
    
    logger.info(f"Loaded {len(probe_positions)} probe positions")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # === Subsample individuals from VCF ===
        # If pop_file provided: stratified sampling (n_samples per population).
        # Otherwise: pure random subsample up to n_samples total.
        # Same individuals are used for both full-VCF and probe analyses.
        import subprocess, random as _random
        _random.seed(seed)
        try:
            header_result = subprocess.run(
                ["bcftools", "query", "-l", str(vcf_path)],
                capture_output=True, text=True, check=False
            )
            all_samples = [s for s in header_result.stdout.strip().splitlines() if s]
        except Exception as e:
            return Err(f"Failed to read sample IDs from VCF: {e}")

        if not all_samples:
            return Err("No samples found in VCF header")

        vcf_sample_set = set(all_samples)

        if pop_file is not None:
            # Stratified: sample n_samples per population so every pop is represented
            # Only consider members that actually exist in the VCF header
            pop_result = load_pop_file(pop_file)
            if pop_result.is_err():
                return Err(f"Failed to load pop file for stratified sampling: {pop_result.unwrap_err()}")
            pop_dict_strat = pop_result.unwrap()
            subsampled = []
            skipped_pops = []
            for pop_id, members in pop_dict_strat.items():
                if pop_id.lower() == 'unknown':
                    continue
                vcf_members = [m for m in members if m in vcf_sample_set]
                if not vcf_members:
                    skipped_pops.append(pop_id)
                    continue
                chosen = sorted(_random.sample(vcf_members, min(n_samples, len(vcf_members))))
                subsampled.extend(chosen)
            if skipped_pops:
                logger.warning(f"Populations with no matching VCF samples (skipped): {skipped_pops}")
            subsampled = sorted(set(subsampled))
            logger.info(
                f"Stratified subsampling: {len(subsampled)} individuals across "
                f"{sum(1 for p in pop_dict_strat if p.lower() != 'unknown')} populations "
                f"(≤{n_samples}/pop, seed={seed})"
            )
        elif len(all_samples) > n_samples:
            subsampled = sorted(_random.sample(all_samples, n_samples))
            logger.info(f"Subsampled {n_samples}/{len(all_samples)} individuals (seed={seed})")
        else:
            subsampled = all_samples
            logger.info(f"Using all {len(all_samples)} individuals")

        full_vcf_sub = tmpdir / "full_samples_subset.vcf.gz"
        sample_sub_result = subset_vcf_by_samples(vcf_path, subsampled, full_vcf_sub, threads=threads)
        if sample_sub_result.is_err():
            return Err(f"Failed to subset VCF by samples: {sample_sub_result.unwrap_err()}")

        # === Run PCA on sample-subsetted full VCF ===
        logger.info("Running PCA on full VCF...")
        full_pca_result = run_plink_pca(full_vcf_sub, tmpdir / "full", n_pcs=n_pcs, threads=threads)

        if full_pca_result.is_err():
            return Err(f"PCA on full VCF failed: {full_pca_result.unwrap_err()}")

        pca_full, var_full, samples_full = full_pca_result.unwrap()
        logger.info(f"Full VCF PCA: {len(samples_full)} samples, {n_pcs} PCs")

        # === Create probe subset VCF (positions, from sample-subsetted VCF) ===
        logger.info("Creating probe-subset VCF...")
        probe_vcf = tmpdir / "probe_subset.vcf.gz"
        subset_result = subset_vcf_by_positions(full_vcf_sub, probe_positions, probe_vcf, threads=threads)

        if subset_result.is_err():
            return Err(f"VCF subsetting failed: {subset_result.unwrap_err()}")

        # === Run PCA on probe VCF ===
        logger.info("Running PCA on probe-subset VCF...")
        probe_pca_result = run_plink_pca(probe_vcf, tmpdir / "probe", n_pcs=n_pcs, threads=threads)

        if probe_pca_result.is_err():
            return Err(f"PCA on probe VCF failed: {probe_pca_result.unwrap_err()}")

        pca_probe, var_probe, samples_probe = probe_pca_result.unwrap()
        logger.info(f"Probe VCF PCA: {len(samples_probe)} samples, {n_pcs} PCs")
    
    # === Generate comparison plots ===
    logger.info("Generating PCA comparison plots...")
    
    # PC1 vs PC2
    pca_plot_12 = Path(str(output_prefix) + ".pca_PC1_PC2.png")
    plot_result = generate_pca_comparison_plot(
        pca_full, pca_probe, var_full, var_probe,
        pca_plot_12, pop_file=pop_file, pc_x=1, pc_y=2
    )
    if plot_result.is_err():
        logger.warning(f"PCA plot failed: {plot_result.unwrap_err()}")
    
    # Calculate Procrustes similarity between PCAs
    procrustes_similarity = None
    try:
        from scipy.spatial import procrustes
        
        # Align on common samples
        common_samples = set(pca_full['ID']) & set(pca_probe['ID'])
        if len(common_samples) > 2:
            full_aligned = pca_full[pca_full['ID'].isin(common_samples)].sort_values('ID')
            probe_aligned = pca_probe[pca_probe['ID'].isin(common_samples)].sort_values('ID')
            
            # Use first 3 PCs for comparison
            n_compare = min(3, n_pcs)
            pc_cols = [f'PC{i}' for i in range(1, n_compare + 1)]
            
            mtx1 = full_aligned[pc_cols].values
            mtx2 = probe_aligned[pc_cols].values
            
            mtx1_std, mtx2_std, disparity = procrustes(mtx1, mtx2)
            procrustes_similarity = 1 - disparity  # Convert disparity to similarity
            logger.info(f"Procrustes similarity (PC1-{n_compare}): {procrustes_similarity:.4f}")
    except ImportError:
        logger.warning("scipy not available, skipping Procrustes analysis")
    except Exception as e:
        logger.warning(f"Procrustes analysis failed: {e}")
    
    # Save eigenvectors
    pca_full_path = Path(str(output_prefix) + ".pca_full.tsv")
    pca_probe_path = Path(str(output_prefix) + ".pca_probe.tsv")
    
    pca_full.to_csv(pca_full_path, sep='\t', index=False)
    pca_probe.to_csv(pca_probe_path, sep='\t', index=False)
    logger.info(f"Saved PCA results to {pca_full_path} and {pca_probe_path}")
    
    # Save summary
    summary_path = Path(str(output_prefix) + ".pca_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("PCA Assessment Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Full VCF: {vcf_path}\n")
        f.write(f"Probe SNPs: {snp_tsv_path}\n")
        if pop_file:
            f.write(f"Pop file: {pop_file}\n")
        f.write(f"\nSamples: {len(samples_full)}\n")
        f.write(f"Probe positions: {len(probe_positions)}\n\n")
        
        f.write("Variance Explained (Full VCF):\n")
        f.write("-" * 30 + "\n")
        for pc, var in sorted(var_full.items(), key=lambda x: int(x[0][2:])):
            f.write(f"  {pc}: {var:.2f}%\n")
        
        f.write("\nVariance Explained (Probe Set):\n")
        f.write("-" * 30 + "\n")
        for pc, var in sorted(var_probe.items(), key=lambda x: int(x[0][2:])):
            f.write(f"  {pc}: {var:.2f}%\n")
        
        if procrustes_similarity is not None:
            f.write(f"\nProcrustes Similarity: {procrustes_similarity:.4f}\n")
            f.write("(1.0 = identical, 0.0 = completely different)\n")
    
    logger.info(f"Saved summary to {summary_path}")
    
    stats_dict = {
        "n_samples": len(samples_full),
        "n_probe_positions": len(probe_positions),
        "variance_full": var_full,
        "variance_probe": var_probe,
        "procrustes_similarity": float(procrustes_similarity) if procrustes_similarity is not None else None,
        "pca_full_file": str(pca_full_path),
        "pca_probe_file": str(pca_probe_path),
        "plot_pc12": str(pca_plot_12),
        "summary_file": str(summary_path),
    }
    
    return Ok(stats_dict)


# =============================================================================
# SFS (Site Frequency Spectrum) Analysis Functions using easySFS
# =============================================================================

def check_easysfs_available() -> Result[str, str]:
    """
    Check if easySFS is available.
    
    Looks for easySFS in the following order:
    1. Bundled easySFS in eprobe.external (preferred)
    2. EASYSFS_PATH environment variable
    3. easySFS.py in PATH
    4. easySFS in PATH
    
    Returns:
        Result containing path to easySFS executable
    """
    import shutil
    import os
    
    # Check bundled easySFS first (preferred)
    try:
        from eprobe.external import get_easysfs_path
        bundled_path = get_easysfs_path()
        if bundled_path.exists():
            return Ok(str(bundled_path))
    except ImportError:
        pass
    
    # Check environment variable
    env_path = os.environ.get("EASYSFS_PATH")
    if env_path and Path(env_path).exists():
        return Ok(env_path)
    
    # Check PATH
    easysfs_cmd = shutil.which("easySFS.py") or shutil.which("easySFS")
    if easysfs_cmd is not None:
        return Ok(easysfs_cmd)
    
    return Err("easySFS not found. This should not happen as easySFS is bundled with eProbe.")


def load_pop_file(
    pop_file: Path,
) -> Result[Dict[str, List[str]], str]:
    """
    Load population assignment file.
    
    Args:
        pop_file: Path to pop file (sample_id<tab>pop_id)
        
    Returns:
        Result containing dict mapping pop_id -> list of sample_ids
    """
    # Validate file exists
    if not pop_file.exists():
        return Err(f"Population file not found: {pop_file}")
    
    try:
        pop_dict = defaultdict(list)
        line_num = 0
        with open(pop_file, 'r') as f:
            for line in f:
                line_num += 1
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                # Support both tab and space delimiters
                parts = line.split('\t') if '\t' in line else line.split()
                if len(parts) >= 2:
                    sample_id, pop_id = parts[0].strip(), parts[1].strip()
                    if sample_id and pop_id:
                        pop_dict[pop_id].append(sample_id)
                elif len(parts) == 1:
                    logger.warning(f"Line {line_num} in pop file has only one field, skipping: {line}")
        
        if not pop_dict:
            return Err("No valid population assignments found in pop file. "
                       "Expected format: sample_id<tab>pop_id")
        
        # Check for duplicate samples
        all_samples = []
        for samples in pop_dict.values():
            all_samples.extend(samples)
        if len(all_samples) != len(set(all_samples)):
            logger.warning("Duplicate sample IDs found in pop file")
        
        return Ok(dict(pop_dict))
        
    except Exception as e:
        return Err(f"Failed to load pop file: {e}")


MAX_SFS_POPULATIONS = 3


def select_sfs_populations(
    pop_dict: Dict[str, List[str]],
    pops: Optional[List[str]] = None,
) -> Result[tuple, str]:
    """
    Select populations for SFS analysis.

    - Samples labeled 'unknown' (case-insensitive) are always excluded.
    - If *pops* is given, validate and use those populations (max 3).
    - Otherwise, auto-select the first MAX_SFS_POPULATIONS sorted populations.

    Returns:
        Ok((filtered_pop_dict, selected_pop_ids))
    """
    # Exclude 'unknown' population
    valid = {k: v for k, v in pop_dict.items() if k.lower() != "unknown"}
    n_unknown = sum(len(v) for k, v in pop_dict.items() if k.lower() == "unknown")
    if n_unknown:
        logger.info(f"Excluded {n_unknown} samples labeled 'unknown' from SFS analysis")

    if not valid:
        return Err("No valid populations found after excluding 'unknown'")

    if pops is not None:
        if len(pops) > MAX_SFS_POPULATIONS:
            return Err(
                f"SFS supports at most {MAX_SFS_POPULATIONS} populations, got {len(pops)}: {pops}. "
                f"Use --pops to specify ≤{MAX_SFS_POPULATIONS} populations."
            )
        missing = [p for p in pops if p not in valid]
        if missing:
            return Err(
                f"Specified populations not found in pop_file: {missing}. "
                f"Available: {sorted(valid.keys())}"
            )
        selected = pops
        logger.info(f"Using user-specified populations for SFS: {selected}")
    else:
        all_pops = sorted(valid.keys())
        selected = all_pops[:MAX_SFS_POPULATIONS]
        if len(all_pops) > MAX_SFS_POPULATIONS:
            logger.warning(
                f"Found {len(all_pops)} populations ({all_pops}), "
                f"using first {MAX_SFS_POPULATIONS}: {selected}. "
                f"Use --pops to specify which populations to use."
            )
        logger.info(f"Selected populations for SFS: {selected}")

    filtered = {k: valid[k] for k in selected}
    return Ok((filtered, selected))


def subsample_populations(
    pop_dict: Dict[str, List[str]],
    n_samples_per_pop: int,
    specified_samples: Optional[List[str]] = None,
    seed: int = 42,
) -> Result[Dict[str, List[str]], str]:
    """
    Subsample populations to equal size.
    
    Args:
        pop_dict: Dict mapping pop_id -> list of sample_ids
        n_samples_per_pop: Number of samples to select per population
        specified_samples: Optional list of specific samples to use (overrides random)
        seed: Random seed for reproducibility
        
    Returns:
        Result containing subsampled pop dict
    """
    import random
    random.seed(seed)
    
    # Validate inputs
    if not pop_dict:
        return Err("Empty population dictionary")
    
    if n_samples_per_pop < 1:
        return Err(f"n_samples_per_pop must be >= 1, got {n_samples_per_pop}")
    
    if specified_samples is not None:
        # Use specified samples, assign to their populations
        specified_set = set(specified_samples)
        subsampled = defaultdict(list)
        
        for pop_id, samples in pop_dict.items():
            for s in samples:
                if s in specified_set:
                    subsampled[pop_id].append(s)
        
        if not subsampled:
            all_samples = set()
            for samples in pop_dict.values():
                all_samples.update(samples)
            return Err(f"None of the specified samples found in pop file. "
                       f"Available samples: {sorted(all_samples)[:10]}...")
        
        # Check if any populations have no samples
        empty_pops = [p for p in pop_dict.keys() if p not in subsampled]
        if empty_pops:
            logger.warning(f"Populations with no specified samples: {empty_pops}")
        
        return Ok(dict(subsampled))
    
    # Random subsampling
    subsampled = {}
    min_pop_size = min(len(samples) for samples in pop_dict.values())
    
    if n_samples_per_pop > min_pop_size:
        logger.warning(f"Requested {n_samples_per_pop} samples per pop, "
                       f"but smallest pop has {min_pop_size}. Using {min_pop_size}.")
        n_samples_per_pop = min_pop_size
    
    for pop_id, samples in pop_dict.items():
        subsampled[pop_id] = random.sample(samples, n_samples_per_pop)
    
    return Ok(subsampled)


def create_subsampled_pop_file(
    pop_dict: Dict[str, List[str]],
    output_path: Path,
) -> Result[Path, str]:
    """
    Create a population file for easySFS from subsampled populations.
    
    Args:
        pop_dict: Dict mapping pop_id -> list of sample_ids
        output_path: Output path for pops file
        
    Returns:
        Result containing path to pops file
    """
    try:
        with open(output_path, 'w') as f:
            for pop_id, samples in pop_dict.items():
                for sample in samples:
                    f.write(f"{sample}\t{pop_id}\n")
        return Ok(output_path)
    except Exception as e:
        return Err(f"Failed to create pop file: {e}")


def run_easysfs(
    vcf_path: Path,
    pops_file: Path,
    output_dir: Path,
    projection: Optional[str] = None,
    preview_only: bool = False,
) -> Result[Path, str]:
    """
    Run easySFS to calculate Site Frequency Spectrum from VCF.
    
    Args:
        vcf_path: Path to VCF file
        pops_file: Population assignment file (sample<tab>pop)
        output_dir: Output directory for SFS files
        projection: Projection values (e.g., "10,10,10" for 3 pops)
        preview_only: If True, only run preview mode
        
    Returns:
        Result containing path to output directory
    """
    import subprocess
    import sys

    easysfs_result = check_easysfs_available()
    if easysfs_result.is_err():
        return Err(easysfs_result.unwrap_err())

    easysfs_cmd = easysfs_result.unwrap()

    # Always invoke through the active Python interpreter so that the correct
    # conda/virtualenv environment (including easySFS dependencies) is used,
    # regardless of the shebang line in the script.
    if easysfs_cmd.endswith(".py"):
        cmd = [sys.executable, easysfs_cmd, "-i", str(vcf_path), "-p", str(pops_file)]
    else:
        cmd = [easysfs_cmd, "-i", str(vcf_path), "-p", str(pops_file)]

    if preview_only:
        cmd.append("--preview")
    else:
        cmd.extend(["-o", str(output_dir)])
        if projection is not None:
            cmd.extend(["--proj", projection])
        cmd.extend(["-a", "-f", "-y"])  # -a: all SNPs, -f: force overwrite, -y: no prompt

    logger.info(f"Running easySFS: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            cmd_str = ' '.join(cmd)
            stderr = result.stderr.strip() or "(empty)"
            stdout = result.stdout.strip() or "(empty)"
            return Err(
                f"easySFS failed (exit {result.returncode}):\n"
                f"  cmd:    {cmd_str}\n"
                f"  stderr: {stderr}\n"
                f"  stdout: {stdout}"
            )
        if preview_only:
            # Return preview output as text
            logger.info(f"easySFS preview:\n{result.stdout}")
            return Ok(output_dir)  # Return dir path even for preview
        
        return Ok(output_dir)
        
    except Exception as e:
        return Err(f"Failed to run easySFS: {e}")


def parse_easysfs_preview(
    vcf_path: Path,
    pops_file: Path,
) -> Result[Dict[str, List[Tuple[int, int]]], str]:
    """
    Run easySFS preview to get optimal projection values.
    
    Args:
        vcf_path: Path to VCF file
        pops_file: Population assignment file
        
    Returns:
        Result containing dict of pop_id -> [(projection, n_segregating_sites), ...]
    """
    import subprocess
    import sys

    easysfs_result = check_easysfs_available()
    if easysfs_result.is_err():
        return Err(easysfs_result.unwrap_err())

    easysfs_cmd = easysfs_result.unwrap()

    if easysfs_cmd.endswith(".py"):
        cmd = [sys.executable, easysfs_cmd, "-i", str(vcf_path), "-p", str(pops_file), "--preview", "-a"]
    else:
        cmd = [easysfs_cmd, "-i", str(vcf_path), "-p", str(pops_file), "--preview", "-a"]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            stderr = result.stderr.strip() or "(empty)"
            stdout = result.stdout.strip() or "(empty)"
            return Err(
                f"easySFS preview failed (exit {result.returncode}):\n"
                f"  stderr: {stderr}\n"
                f"  stdout: {stdout}"
            )
        # Parse preview output
        # Format: pop_name\n(proj, n_sites) (proj, n_sites) ...
        preview_dict = {}
        current_pop = None
        
        for line in result.stdout.split('\n'):
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('('):
                # This is projection data
                if current_pop:
                    projections = []
                    for match in re.findall(r'\((\d+),\s*([\d.]+)\)', line):
                        proj, n_sites = int(match[0]), int(float(match[1]))
                        projections.append((proj, n_sites))
                    preview_dict[current_pop] = projections
            else:
                # This is a population name
                current_pop = line
        
        return Ok(preview_dict)
        
    except Exception as e:
        return Err(f"Failed to parse easySFS preview: {e}")


def find_optimal_projection(
    preview_dict: Dict[str, List[Tuple[int, int]]],
    max_projection: Optional[int] = None,
) -> Dict[str, int]:
    """
    Find optimal projection values that maximize segregating sites.
    
    Args:
        preview_dict: Dict from parse_easysfs_preview
        max_projection: Optional maximum projection value to consider
        
    Returns:
        Dict mapping pop_id -> optimal projection value
    """
    optimal = {}
    
    for pop_id, projections in preview_dict.items():
        if not projections:
            continue
        
        # Filter by max_projection if specified
        if max_projection is not None:
            projections = [(p, n) for p, n in projections if p <= max_projection]
        
        if projections:
            # Find projection that maximizes segregating sites
            best = max(projections, key=lambda x: x[1])
            optimal[pop_id] = best[0]
    
    return optimal


def subset_vcf_by_positions(
    vcf_path: Path,
    positions: Set[Tuple[str, int]],
    output_path: Path,
    threads: int = 1,
) -> Result[Path, str]:
    """
    Create a subset VCF containing only specified positions using bcftools.
    
    Args:
        vcf_path: Input VCF path
        positions: Set of (chrom, pos) tuples to include
        output_path: Output VCF path
        
    Returns:
        Result containing output path
    """
    import subprocess
    import shutil
    import tempfile
    
    if shutil.which("bcftools") is None:
        return Err("bcftools not found in PATH")
    
    try:
        # Create regions file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            regions_file = f.name
            for chrom, pos in sorted(positions):
                f.write(f"{chrom}\t{pos}\n")
        
        # Run bcftools view with regions
        cmd = [
            "bcftools", "view",
            "-T", regions_file,
            "--threads", str(threads),
            "-O", "z",
            "-o", str(output_path),
            str(vcf_path)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        
        # Clean up
        Path(regions_file).unlink()
        
        if result.returncode != 0:
            return Err(f"bcftools failed: {result.stderr}")
        
        # Index the output
        subprocess.run(["bcftools", "index", "--threads", str(threads), "-f", str(output_path)], check=False)

        return Ok(output_path)

    except Exception as e:
        return Err(f"Failed to subset VCF: {e}")


def subset_vcf_by_samples(
    vcf_path: Path,
    sample_ids: List[str],
    output_path: Path,
    threads: int = 1,
) -> Result[Path, str]:
    """
    Create a subset VCF containing only specified samples using bcftools.

    Args:
        vcf_path: Input VCF path
        sample_ids: List of sample IDs to keep
        output_path: Output VCF path

    Returns:
        Result containing output path
    """
    import subprocess
    import shutil
    import tempfile

    if shutil.which("bcftools") is None:
        return Err("bcftools not found in PATH")

    try:
        # Write sample list to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            samples_file = f.name
            for s in sample_ids:
                f.write(f"{s}\n")

        cmd = [
            "bcftools", "view",
            "-S", samples_file,
            "--force-samples",
            "--threads", str(threads),
            "-O", "z",
            "-o", str(output_path),
            str(vcf_path)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        Path(samples_file).unlink(missing_ok=True)

        if result.returncode != 0:
            return Err(f"bcftools view -S failed (exit {result.returncode}): {result.stderr.strip()}")

        subprocess.run(["bcftools", "index", "--threads", str(threads), "-f", str(output_path)], check=False)
        return Ok(output_path)

    except Exception as e:
        return Err(f"Failed to subset VCF by samples: {e}")


def subset_vcf_by_samples_and_positions(
    vcf_path: Path,
    sample_ids: List[str],
    positions: Set[Tuple[str, int]],
    output_path: Path,
    threads: int = 1,
) -> Result[Path, str]:
    """
    Subset VCF by samples AND positions in a single bcftools call.

    Equivalent to running bcftools view -S samples.txt -T positions.txt
    in one pass, which avoids writing/reading an intermediate VCF.
    """
    import subprocess
    import shutil
    import tempfile

    if shutil.which("bcftools") is None:
        return Err("bcftools not found in PATH")

    samples_file = None
    regions_file = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            samples_file = f.name
            for s in sample_ids:
                f.write(f"{s}\n")

        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            regions_file = f.name
            for chrom, pos in sorted(positions):
                f.write(f"{chrom}\t{pos}\n")

        cmd = [
            "bcftools", "view",
            "-S", samples_file,
            "--force-samples",
            "-T", regions_file,
            "--threads", str(threads),
            "-O", "z",
            "-o", str(output_path),
            str(vcf_path)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            return Err(f"bcftools view -S -T failed (exit {result.returncode}): {result.stderr.strip()}")

        subprocess.run(["bcftools", "index", "--threads", str(threads), "-f", str(output_path)], check=False)
        return Ok(output_path)

    except Exception as e:
        return Err(f"Failed to subset VCF by samples+positions: {e}")
    finally:
        if samples_file:
            Path(samples_file).unlink(missing_ok=True)
        if regions_file:
            Path(regions_file).unlink(missing_ok=True)


def parse_dadi_1d_sfs(sfs_dir: Path, pop_id: str) -> Result[np.ndarray, str]:
    """
    Parse 1D SFS from easySFS dadi output.
    
    Args:
        sfs_dir: easySFS output directory
        pop_id: Population ID
        
    Returns:
        Result containing 1D SFS array (counts)
    """
    dadi_dir = sfs_dir / "dadi"
    
    # Find the 1D SFS file for this population
    # Format: {pop_id}-{projection}.sfs  (e.g., Pop_A-6.sfs)
    # Exclude joint SFS files which have format {pop1}-{pop2}.sfs
    sfs_files = []
    if dadi_dir.exists():
        for f in dadi_dir.glob(f"{pop_id}-*.sfs"):
            # Check if it's a 1D SFS (ends with number) not a 2D SFS (ends with pop name)
            name = f.stem  # e.g., "Pop_A-6" or "Pop_A-Pop_B"
            suffix = name.split('-')[-1]
            if suffix.isdigit():  # 1D SFS has projection number
                sfs_files.append(f)
    
    if not sfs_files:
        return Err(f"No SFS file found for population {pop_id}")
    
    sfs_path = sfs_files[0]
    
    try:
        with open(sfs_path, 'r') as f:
            lines = f.readlines()
        
        # Skip header line, parse counts
        data_lines = [l.strip() for l in lines if l.strip() and not l.startswith('#')]
        
        if len(data_lines) < 2:
            return Err(f"Invalid SFS format in {sfs_path}")
        
        # First data line (line 1) = header with size
        # Second line (line 2) = SFS counts  
        # Third line (line 3) = mask (if present)
        header = data_lines[0].split()
        n_bins = int(header[0])
        
        counts_str = data_lines[1]
        counts = np.array([float(x) for x in counts_str.split()])
        
        # Validate length
        if len(counts) != n_bins:
            return Err(f"SFS length mismatch: expected {n_bins}, got {len(counts)}")
        
        return Ok(counts)
        
    except Exception as e:
        return Err(f"Failed to parse SFS: {e}")


def parse_dadi_2d_sfs(sfs_dir: Path, pop1: str, pop2: str) -> Result[np.ndarray, str]:
    """
    Parse 2D joint SFS from easySFS dadi output.
    
    Args:
        sfs_dir: easySFS output directory
        pop1: First population ID
        pop2: Second population ID
        
    Returns:
        Result containing 2D SFS array (counts)
    """
    dadi_dir = sfs_dir / "dadi"
    
    # Try both orderings
    patterns = [f"{pop1}-{pop2}.sfs", f"{pop2}-{pop1}.sfs"]
    sfs_path = None
    swapped = False
    
    for i, pattern in enumerate(patterns):
        candidate = dadi_dir / pattern
        if candidate.exists():
            sfs_path = candidate
            swapped = (i == 1)
            break
    
    if sfs_path is None:
        return Err(f"No joint SFS file found for {pop1}-{pop2}")
    
    try:
        with open(sfs_path, 'r') as f:
            content = f.read()
        
        # Parse dadi 2D SFS format
        # First line format: "n1 n2 folded|unfolded \"Pop1\" \"Pop2\""
        # Then a flattened array of n1*n2 values
        lines = content.strip().split('\n')
        if not lines:
            return Err(f"Empty SFS file: {sfs_path}")
        
        # Parse header - extract dimensions (first two integers)
        header = lines[0].split()
        dim1, dim2 = int(header[0]), int(header[1])
        
        # Collect all values from remaining content
        values_str = ' '.join(lines[1:])
        all_values = [float(x) for x in values_str.split()]
        
        # dadi format contains counts + mask (2 * n1 * n2 values)
        expected = dim1 * dim2
        if len(all_values) == 2 * expected:
            # Data + mask format: take only counts (first half)
            counts_flat = all_values[:expected]
        elif len(all_values) == expected:
            # Just counts
            counts_flat = all_values
        else:
            return Err(f"SFS shape mismatch: expected {expected} or {2*expected} values, got {len(all_values)}")
        
        counts = np.array(counts_flat).reshape(dim1, dim2)
        
        if swapped:
            counts = counts.T
        
        return Ok(counts)
        
    except Exception as e:
        return Err(f"Failed to parse 2D SFS: {e}")


def generate_multipop_sfs_heatmap(
    sfs_dir: Path,
    pop_ids: List[str],
    output_path: Path,
    title: str = "Site Frequency Spectrum",
    projection_values: Optional[Dict[str, int]] = None,
) -> Result[Path, str]:
    """
    Generate multi-population SFS heatmap like the reference image.
    
    Diagonal blocks show 1D SFS for each population.
    Off-diagonal blocks show pairwise joint 2D SFS.
    
    Args:
        sfs_dir: easySFS output directory
        pop_ids: List of population IDs
        output_path: Output image path
        title: Plot title
        projection_values: Dict of pop_id -> projection value for axis labels
        
    Returns:
        Result containing output path
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.gridspec as gridspec
    except ImportError:
        return Err("Plotting requires matplotlib")
    
    # Validate inputs
    if not pop_ids:
        return Err("No populations provided")
    
    n_pops = len(pop_ids)
    
    # Load all 1D and 2D SFS
    sfs_1d = {}
    sfs_2d = {}
    
    for pop in pop_ids:
        result = parse_dadi_1d_sfs(sfs_dir, pop)
        if result.is_ok():
            sfs_1d[pop] = result.unwrap()
        else:
            logger.warning(f"Could not load 1D SFS for {pop}: {result.unwrap_err()}")
    
    for i, pop1 in enumerate(pop_ids):
        for j, pop2 in enumerate(pop_ids):
            if i < j:
                result = parse_dadi_2d_sfs(sfs_dir, pop1, pop2)
                if result.is_ok():
                    sfs_2d[(pop1, pop2)] = result.unwrap()
                else:
                    logger.warning(f"Could not load 2D SFS for {pop1}-{pop2}")
    
    # Check if we have any data
    if not sfs_1d and not sfs_2d:
        return Err("No SFS data could be loaded. Check easySFS output.")
    
    # Custom colormap: dark red -> red -> orange -> yellow -> white
    colors = ['#8B0000', '#FF0000', '#FF4500', '#FFA500', '#FFD700', '#FFFF00', '#FFFACD', '#FFFFFF']
    cmap = LinearSegmentedColormap.from_list('sfs_cmap', colors[::-1])  # Reverse for high=red
    cmap.set_bad('lightgray')  # Color for masked (monomorphic) cells

    # Log-transform for better dynamic range (standard SFS visualization)
    plot_sfs_1d = {}
    for pop, arr in sfs_1d.items():
        log_arr = np.log10(1.0 + arr)
        masked = np.ma.array(log_arr)
        masked[0] = np.ma.masked   # mask monomorphic (all ancestral)
        masked[-1] = np.ma.masked  # mask fixed (all derived)
        plot_sfs_1d[pop] = masked

    plot_sfs_2d = {}
    for key, arr in sfs_2d.items():
        log_arr = np.log10(1.0 + arr)
        masked = np.ma.array(log_arr)
        masked[0, 0] = np.ma.masked    # mask all-ancestral corner
        masked[-1, -1] = np.ma.masked  # mask all-derived corner
        plot_sfs_2d[key] = masked

    # Global vmax from unmasked log-transformed values
    all_values = []
    for arr in plot_sfs_1d.values():
        all_values.extend(arr.compressed())
    for arr in plot_sfs_2d.values():
        all_values.extend(arr.compressed())
    
    vmax = max(all_values) if all_values else 1.0
    
    # Create figure with custom layout
    # Main grid for heatmaps + space for colorbar
    fig_width = 2.5 * n_pops + 1.0  # Extra space for colorbar
    fig_height = 2.5 * n_pops + 0.8  # Extra space for title
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    # Create gridspec: n_pops rows x (n_pops + 1) cols (last col for colorbar)
    gs = gridspec.GridSpec(
        n_pops, n_pops + 1,
        figure=fig,
        width_ratios=[1] * n_pops + [0.1],  # Small width for colorbar column
        wspace=0,  # No horizontal gap
        hspace=0,  # No vertical gap
        left=0.12, right=0.88, top=0.92, bottom=0.08
    )
    
    axes = {}
    
    for i, pop_row in enumerate(pop_ids):
        for j, pop_col in enumerate(pop_ids):
            ax = fig.add_subplot(gs[i, j])
            axes[(i, j)] = ax
            
            if i == j:
                # Diagonal: 1D SFS as column heatmap
                if pop_row in plot_sfs_1d:
                    sfs_data = plot_sfs_1d[pop_row]
                    im = ax.imshow(
                        sfs_data.reshape(-1, 1),
                        aspect='auto',
                        cmap=cmap,
                        vmin=0,
                        vmax=vmax,
                        origin='upper'
                    )
                    ax.set_xticks([])
                else:
                    ax.text(0.5, 0.5, 'N/A', ha='center', va='center', fontsize=10)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    
            elif i > j:
                # Lower triangle: 2D joint SFS
                key = (pop_col, pop_row) if (pop_col, pop_row) in plot_sfs_2d else (pop_row, pop_col)
                if key in plot_sfs_2d:
                    sfs_data = plot_sfs_2d[key]
                    if key == (pop_row, pop_col):
                        sfs_data = sfs_data.T
                    
                    im = ax.imshow(
                        sfs_data,
                        aspect='auto',
                        cmap=cmap,
                        vmin=0,
                        vmax=vmax,
                        origin='upper'
                    )
                else:
                    ax.text(0.5, 0.5, 'N/A', ha='center', va='center', fontsize=10)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    
            else:
                # Upper triangle: empty (white), no border
                ax.axis('off')
            
            # Y-axis: only show ticks on leftmost column
            if j == 0 and i >= j:
                if i == j and pop_row in plot_sfs_1d:
                    n_bins = len(plot_sfs_1d[pop_row])
                    ax.set_yticks(range(n_bins))
                    ax.set_yticklabels(range(n_bins), fontsize=7)
                elif i > j:
                    key = (pop_col, pop_row) if (pop_col, pop_row) in plot_sfs_2d else (pop_row, pop_col)
                    if key in plot_sfs_2d:
                        sfs_data = plot_sfs_2d[key]
                        if key == (pop_row, pop_col):
                            sfs_data = sfs_data.T
                        ax.set_yticks(range(sfs_data.shape[0]))
                        ax.set_yticklabels(range(sfs_data.shape[0]), fontsize=7)
            elif i >= j:
                ax.set_yticks([])
                # Also hide tick marks on non-leftmost columns
                ax.tick_params(left=False)
            
            # X-axis: hide all ticks
            if i >= j:
                ax.set_xticks([])
                ax.tick_params(bottom=False)
            
            # Black border only for diagonal and lower triangle
            if i >= j:
                for spine in ax.spines.values():
                    spine.set_visible(True)
                    spine.set_color('black')
                    spine.set_linewidth(1.5)
    
    # Add population labels on top (X-axis titles)
    for j, pop_col in enumerate(pop_ids):
        ax = axes[(0, j)]
        ax.set_title(pop_col, fontsize=11, fontweight='bold', pad=8)
    
    # Add population labels on left (Y-axis titles)
    for i, pop_row in enumerate(pop_ids):
        ax = axes[(i, 0)]
        ax.set_ylabel(pop_row, fontsize=11, fontweight='bold', labelpad=8)
    
    # Add colorbar on the right side (shortened to half height)
    cbar_ax = fig.add_axes([0.91, 0.65, 0.025, 0.175])  # Height reduced from 0.35 to 0.175
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, vmax)),
        cax=cbar_ax
    )
    cbar.set_label('log\u2081\u2080(count + 1)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    
    # Add title
    fig.suptitle(title, fontsize=13, fontweight='bold', y=0.98)
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    return Ok(output_path)


def run_sfs_assessment(
    snp_tsv_path: Path,
    vcf_path: Path,
    output_prefix: Path,
    pop_file: Path,
    n_samples_per_pop: int = 5,
    specified_samples: Optional[List[str]] = None,
    projection: Optional[str] = None,
    pops: Optional[List[str]] = None,
    seed: int = 42,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Compare Site Frequency Spectrum between full VCF and probe-covered SNPs.
    
    Uses easySFS to calculate multi-population joint SFS with proper
    down-projection to handle missing data.
    
    Args:
        snp_tsv_path: Path to selected SNPs TSV file
        vcf_path: Path to original VCF file
        output_prefix: Output file prefix
        pop_file: Population assignment file (sample_id<tab>pop_id)
        n_samples_per_pop: Number of samples to randomly select per population
        specified_samples: Optional list of specific samples to use
        projection: Projection string (e.g., "10,10,10"), auto-calculated if None
        seed: Random seed for reproducibility
        threads: Number of threads for bcftools
        verbose: Enable verbose logging
        
    Returns:
        Result containing SFS assessment statistics
    """
    import tempfile
    import shutil
    
    logger.info("Starting SFS assessment using easySFS")
    logger.info(f"SNP TSV: {snp_tsv_path}")
    logger.info(f"VCF: {vcf_path}")
    logger.info(f"Pop file: {pop_file}")
    logger.info(f"Using {threads} thread(s) for bcftools")
    
    # Check easySFS and bcftools
    easysfs_check = check_easysfs_available()
    if easysfs_check.is_err():
        return Err(easysfs_check.unwrap_err())
    
    if shutil.which("bcftools") is None:
        return Err("bcftools not found in PATH")
    
    # Validate input files exist
    if not snp_tsv_path.exists():
        return Err(f"SNP TSV file not found: {snp_tsv_path}")
    if not vcf_path.exists():
        return Err(f"VCF file not found: {vcf_path}")
    if not pop_file.exists():
        return Err(f"Population file not found: {pop_file}")
    
    # Create output directory
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    # Load population file
    pop_result = load_pop_file(pop_file)
    if pop_result.is_err():
        return Err(pop_result.unwrap_err())
    
    pop_dict = pop_result.unwrap()

    # ── Validate pop file samples against actual VCF header ──────────
    # Pop files often list samples that are absent from the VCF (different
    # cohort, renamed, etc.).  Subsampling from phantom samples causes
    # bcftools --force-samples to silently drop them, and easySFS then
    # crashes with a KeyError because the population has zero genotypes.
    import subprocess as _sp
    try:
        _hdr = _sp.run(
            ["bcftools", "query", "-l", str(vcf_path)],
            capture_output=True, text=True, check=False,
        )
        vcf_samples = set(s for s in _hdr.stdout.strip().splitlines() if s)
    except Exception as e:
        return Err(f"Failed to read sample list from VCF: {e}")

    if not vcf_samples:
        return Err("No samples found in VCF header")

    filtered_pop_dict: Dict[str, List[str]] = {}
    for pop_id, members in pop_dict.items():
        valid = [m for m in members if m in vcf_samples]
        if valid:
            filtered_pop_dict[pop_id] = valid
        else:
            logger.warning(f"Population '{pop_id}': none of {len(members)} samples found in VCF — skipped")

    if not filtered_pop_dict:
        return Err(
            "No population samples match VCF header. "
            "Check that sample IDs in the pop file match those in the VCF."
        )

    n_dropped = sum(len(pop_dict[p]) for p in pop_dict) - sum(len(filtered_pop_dict[p]) for p in filtered_pop_dict)
    if n_dropped > 0:
        logger.info(f"Dropped {n_dropped} pop-file samples not present in VCF")
    pop_dict = filtered_pop_dict

    # Select populations (exclude 'unknown', limit to 3)
    pop_select_result = select_sfs_populations(pop_dict, pops)
    if pop_select_result.is_err():
        return Err(pop_select_result.unwrap_err())
    pop_dict, pop_ids = pop_select_result.unwrap()

    logger.info(f"Populations: {pop_ids}")
    logger.info(f"Samples per pop: {[len(pop_dict[p]) for p in pop_ids]}")
    
    # Subsample populations
    subsample_result = subsample_populations(
        pop_dict, n_samples_per_pop, specified_samples, seed
    )
    if subsample_result.is_err():
        return Err(subsample_result.unwrap_err())
    
    subsampled_pops = subsample_result.unwrap()
    logger.info(f"Subsampled to {n_samples_per_pop} samples per pop")
    
    # Load probe positions
    snp_df_result = SNPDataFrame.from_tsv(snp_tsv_path)
    if snp_df_result.is_err():
        return Err(f"Failed to load SNP TSV: {snp_df_result.unwrap_err()}")
    
    snp_df = snp_df_result.unwrap()
    probe_positions = set()
    for _, row in snp_df.df.iterrows():
        probe_positions.add((str(row['chr']), int(row['pos'])))
    
    logger.info(f"Loaded {len(probe_positions)} probe positions")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Create subsampled pop file
        sub_pop_file = tmpdir / "subsampled_pops.txt"
        pop_file_result = create_subsampled_pop_file(subsampled_pops, sub_pop_file)
        if pop_file_result.is_err():
            return Err(pop_file_result.unwrap_err())
        
        # Determine projection values
        if projection is None:
            # Use 2 * n_samples_per_pop as default projection (diploid)
            proj_val = 2 * n_samples_per_pop
            projection = ",".join([str(proj_val)] * len(pop_ids))
            logger.info(f"Auto projection: {projection}")

        # Build flat list of all subsampled sample IDs (for bcftools -S)
        all_subsampled = [s for samples in subsampled_pops.values() for s in samples]

        # === Subset full VCF to subsampled individuals before running easySFS ===
        # This avoids OOM: original VCF has all samples; we only need the subsampled ones
        logger.info(f"Subsetting full VCF to {len(all_subsampled)} subsampled individuals...")
        full_vcf_sub = tmpdir / "full_samples_subset.vcf.gz"
        full_sample_result = subset_vcf_by_samples(vcf_path, all_subsampled, full_vcf_sub, threads=threads)
        if full_sample_result.is_err():
            return Err(f"Failed to subset full VCF by samples: {full_sample_result.unwrap_err()}")

        # === Run easySFS on sample-subsetted full VCF ===
        logger.info("Running easySFS on full VCF...")
        full_sfs_dir = tmpdir / "full_sfs"

        full_result = run_easysfs(
            full_vcf_sub, sub_pop_file, full_sfs_dir,
            projection=projection, preview_only=False
        )
        
        if full_result.is_err():
            return Err(f"easySFS failed on full VCF: {full_result.unwrap_err()}")
        
        # Count sites in full VCF
        n_sites_full = _count_vcf_sites(vcf_path)

        # === Create probe subset VCF (samples + positions in one bcftools call) ===
        # Read directly from original VCF instead of re-reading the intermediate
        # full_vcf_sub: this avoids one full I/O pass over that large intermediate file.
        logger.info("Creating probe-subset VCF (samples + positions in one pass)...")
        probe_vcf = tmpdir / "probe_subset.vcf.gz"
        subset_result = subset_vcf_by_samples_and_positions(
            vcf_path, all_subsampled, probe_positions, probe_vcf, threads=threads
        )

        if subset_result.is_err():
            return Err(f"VCF subsetting failed: {subset_result.unwrap_err()}")

        # === Run easySFS on probe VCF ===
        logger.info("Running easySFS on probe-subset VCF...")
        probe_sfs_dir = tmpdir / "probe_sfs"

        probe_result = run_easysfs(
            probe_vcf, sub_pop_file, probe_sfs_dir,
            projection=projection, preview_only=False
        )
        
        if probe_result.is_err():
            return Err(f"easySFS failed on probe VCF: {probe_result.unwrap_err()}")
        
        n_sites_probe = len(probe_positions)
        
        # === Generate comparison plots ===
        logger.info("Generating SFS heatmaps...")
        
        # Full VCF heatmap
        full_plot_path = Path(str(output_prefix) + ".sfs_full.png")
        full_plot_result = generate_multipop_sfs_heatmap(
            full_sfs_dir, pop_ids, full_plot_path,
            title=f"Full VCF SFS (n={n_sites_full:,} sites)"
        )
        if full_plot_result.is_err():
            logger.warning(f"Full SFS plot failed: {full_plot_result.unwrap_err()}")
        
        # Probe VCF heatmap
        probe_plot_path = Path(str(output_prefix) + ".sfs_probe.png")
        probe_plot_result = generate_multipop_sfs_heatmap(
            probe_sfs_dir, pop_ids, probe_plot_path,
            title=f"Probe Set SFS (n={n_sites_probe:,} sites)"
        )
        if probe_plot_result.is_err():
            logger.warning(f"Probe SFS plot failed: {probe_plot_result.unwrap_err()}")
        
        # Calculate SFS correlation (using 1D SFS)
        sfs_correlations = {}
        for pop in pop_ids:
            full_1d = parse_dadi_1d_sfs(full_sfs_dir, pop)
            probe_1d = parse_dadi_1d_sfs(probe_sfs_dir, pop)
            
            if full_1d.is_ok() and probe_1d.is_ok():
                f_sfs = full_1d.unwrap()
                p_sfs = probe_1d.unwrap()
                
                # Normalize
                f_sfs = f_sfs / np.sum(f_sfs) if np.sum(f_sfs) > 0 else f_sfs
                p_sfs = p_sfs / np.sum(p_sfs) if np.sum(p_sfs) > 0 else p_sfs
                
                # Ensure same length
                min_len = min(len(f_sfs), len(p_sfs))
                f_sfs = f_sfs[:min_len]
                p_sfs = p_sfs[:min_len]
                
                if not np.all(f_sfs == 0) and not np.all(p_sfs == 0):
                    corr, _ = stats.pearsonr(f_sfs, p_sfs)
                    sfs_correlations[pop] = corr
        
        avg_correlation = np.mean(list(sfs_correlations.values())) if sfs_correlations else np.nan
        
        # Copy dadi output files to output directory
        dadi_output = Path(str(output_prefix) + "_dadi")
        if (full_sfs_dir / "dadi").exists():
            shutil.copytree(full_sfs_dir / "dadi", dadi_output / "full")
        if (probe_sfs_dir / "dadi").exists():
            shutil.copytree(probe_sfs_dir / "dadi", dadi_output / "probe")
    
    # Save summary
    summary_path = Path(str(output_prefix) + ".sfs_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("Site Frequency Spectrum Assessment Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Full VCF: {vcf_path}\n")
        f.write(f"Probe SNPs: {snp_tsv_path}\n")
        f.write(f"Pop file: {pop_file}\n\n")
        f.write(f"Populations: {', '.join(pop_ids)}\n")
        f.write(f"Samples per pop: {n_samples_per_pop}\n")
        f.write(f"Projection: {projection}\n\n")
        f.write(f"Full VCF sites: {n_sites_full:,}\n")
        f.write(f"Probe sites: {n_sites_probe:,}\n")
        f.write(f"Probe coverage: {n_sites_probe / n_sites_full * 100:.2f}%\n\n")
        f.write("SFS Correlations per Population:\n")
        f.write("-" * 30 + "\n")
        for pop, corr in sfs_correlations.items():
            f.write(f"  {pop}: {corr:.4f}\n")
        f.write(f"\nAverage correlation: {avg_correlation:.4f}\n\n")
        f.write("Interpretation:\n")
        f.write("-" * 30 + "\n")
        if np.isnan(avg_correlation):
            f.write("Unable to calculate correlation.\n")
        elif avg_correlation > 0.95:
            f.write("Excellent: SFS highly similar across populations.\n")
        elif avg_correlation > 0.85:
            f.write("Good: SFS similar, probe set captures most signal.\n")
        elif avg_correlation > 0.7:
            f.write("Moderate: Some deviation, consider more probes.\n")
        else:
            f.write("Low: SFS differs substantially.\n")
    
    logger.info(f"Saved summary to {summary_path}")
    
    stats_dict = {
        "populations": pop_ids,
        "n_samples_per_pop": n_samples_per_pop,
        "projection": projection,
        "n_sites_full": n_sites_full,
        "n_sites_probe": n_sites_probe,
        "coverage_percent": n_sites_probe / n_sites_full * 100,
        "sfs_correlations": sfs_correlations,
        "avg_correlation": float(avg_correlation) if not np.isnan(avg_correlation) else None,
        "full_plot": str(full_plot_path) if full_plot_result.is_ok() else None,
        "probe_plot": str(probe_plot_path) if probe_plot_result.is_ok() else None,
        "dadi_output": str(dadi_output),
        "summary_file": str(summary_path),
    }
    
    return Ok(stats_dict)


def _count_vcf_sites(vcf_path: Path) -> int:
    """Count number of variant sites in VCF using bcftools index (fast, no decompression)."""
    import subprocess
    import shutil

    if shutil.which("bcftools") is None:
        return 0

    try:
        # --nrecords reads the index file directly — O(1), no decompression needed
        result = subprocess.run(
            ["bcftools", "index", "--nrecords", str(vcf_path)],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode == 0 and result.stdout.strip().isdigit():
            return int(result.stdout.strip())

        # Fallback: bcftools stats (reads full file but doesn't buffer into Python)
        result2 = subprocess.run(
            ["bcftools", "stats", str(vcf_path)],
            capture_output=True,
            text=True,
            check=False,
        )
        for line in result2.stdout.splitlines():
            if line.startswith("SN") and "number of records" in line:
                return int(line.split()[-1])
    except Exception:
        pass

    return 0


def run_distance_assessment(
    snp_tsv_path: Path,
    vcf_path: Path,
    output_prefix: Path,
    max_sites: int = 100000,
    pop_file: Optional[Path] = None,
    n_samples: int = 100,
    seed: int = 42,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Run IBS distance comparison between full VCF and probe-covered SNPs using PLINK.

    Compares pairwise genetic distances calculated from:
    1. Full VCF (all SNPs)
    2. Probe-covered SNPs only

    Uses PLINK's 1-IBS distance calculation for consistency with standard
    population genetics workflows.

    Args:
        snp_tsv_path: Path to selected SNPs TSV file
        vcf_path: Path to original VCF file
        output_prefix: Output file prefix
        max_sites: Not used (kept for backwards compatibility)
        pop_file: Optional population file. If provided, stratified subsampling
            is used to ensure each population is represented (n_samples per pop).
        n_samples: If pop_file is None: max total individuals (random subsample).
            If pop_file is provided: max individuals per population.
        seed: Random seed for subsampling
        threads: Number of threads for PLINK and bcftools
        verbose: Enable verbose logging

    Returns:
        Result containing distance assessment statistics
    """
    import tempfile
    import shutil
    
    logger.info("Starting IBS distance assessment using PLINK")
    logger.info(f"SNP TSV: {snp_tsv_path}")
    logger.info(f"VCF: {vcf_path}")
    logger.info(f"Using {threads} thread(s) for PLINK and bcftools")
    logger.info(f"Using {threads} thread(s) for PLINK and bcftools")
    
    # Check plink and bcftools
    plink_check = check_plink_available()
    if plink_check.is_err():
        return Err(plink_check.unwrap_err())
    
    if shutil.which("bcftools") is None:
        return Err("bcftools not found in PATH")
    
    # Validate input files
    if not snp_tsv_path.exists():
        return Err(f"SNP TSV file not found: {snp_tsv_path}")
    if not vcf_path.exists():
        return Err(f"VCF file not found: {vcf_path}")
    
    # Load selected SNP positions from TSV (only requires chr and pos columns)
    try:
        snp_df_pandas = pd.read_csv(snp_tsv_path, sep='\t')
        if 'chr' not in snp_df_pandas.columns or 'pos' not in snp_df_pandas.columns:
            return Err(f"SNP TSV must have 'chr' and 'pos' columns. Found: {list(snp_df_pandas.columns)}")
        
        probe_positions = set()
        for _, row in snp_df_pandas.iterrows():
            probe_positions.add((str(row['chr']), int(row['pos'])))
    except Exception as e:
        return Err(f"Failed to load SNP TSV: {e}")
    
    logger.info(f"Loaded {len(probe_positions)} probe positions")

    # Create output directory
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # === Subsample individuals ===
        # If pop_file provided: stratified sampling (n_samples per population).
        # Otherwise: pure random subsample up to n_samples total.
        # Same individuals are used for both full-VCF and probe analyses.
        import subprocess, random as _random
        _random.seed(seed)
        try:
            header_result = subprocess.run(
                ["bcftools", "query", "-l", str(vcf_path)],
                capture_output=True, text=True, check=False
            )
            all_samples = [s for s in header_result.stdout.strip().splitlines() if s]
        except Exception as e:
            return Err(f"Failed to read sample IDs from VCF: {e}")

        if not all_samples:
            return Err("No samples found in VCF header")

        vcf_sample_set = set(all_samples)

        if pop_file is not None:
            # Stratified: sample n_samples per population so every pop is represented
            # Only consider members that actually exist in the VCF header
            pop_result = load_pop_file(pop_file)
            if pop_result.is_err():
                return Err(f"Failed to load pop file for stratified sampling: {pop_result.unwrap_err()}")
            pop_dict_strat = pop_result.unwrap()
            subsampled = []
            skipped_pops = []
            for pop_id, members in pop_dict_strat.items():
                if pop_id.lower() == 'unknown':
                    continue
                vcf_members = [m for m in members if m in vcf_sample_set]
                if not vcf_members:
                    skipped_pops.append(pop_id)
                    continue
                chosen = sorted(_random.sample(vcf_members, min(n_samples, len(vcf_members))))
                subsampled.extend(chosen)
            if skipped_pops:
                logger.warning(f"Populations with no matching VCF samples (skipped): {skipped_pops}")
            subsampled = sorted(set(subsampled))
            logger.info(
                f"Stratified subsampling: {len(subsampled)} individuals across "
                f"{sum(1 for p in pop_dict_strat if p.lower() != 'unknown')} populations "
                f"(≤{n_samples}/pop, seed={seed})"
            )
        elif len(all_samples) > n_samples:
            subsampled = sorted(_random.sample(all_samples, n_samples))
            logger.info(f"Subsampled {n_samples}/{len(all_samples)} individuals (seed={seed})")
        else:
            subsampled = all_samples
            logger.info(f"Using all {len(all_samples)} individuals")

        full_vcf_sub = tmpdir / "full_samples_subset.vcf.gz"
        sample_sub_result = subset_vcf_by_samples(vcf_path, subsampled, full_vcf_sub, threads=threads)
        if sample_sub_result.is_err():
            return Err(f"Failed to subset VCF by samples: {sample_sub_result.unwrap_err()}")

        # === Calculate IBS distance for sample-subsetted full VCF ===
        logger.info("Calculating IBS distances for full VCF using PLINK...")
        full_ibs_result = run_plink_ibs_distance(full_vcf_sub, tmpdir / "full", threads=threads)

        if full_ibs_result.is_err():
            return Err(f"Full VCF IBS failed: {full_ibs_result.unwrap_err()}")

        dist_full, samples_full = full_ibs_result.unwrap()
        n_sites_full = _count_vcf_sites(vcf_path)
        logger.info(f"Full VCF: {len(samples_full)} samples, ~{n_sites_full} sites")

        # === Create probe subset VCF (positions, from sample-subsetted VCF) ===
        logger.info("Creating probe-subset VCF...")
        probe_vcf = tmpdir / "probe_subset.vcf.gz"
        subset_result = subset_vcf_by_positions(full_vcf_sub, probe_positions, probe_vcf, threads=threads)

        if subset_result.is_err():
            err_msg = f"VCF subsetting failed: {subset_result.unwrap_err()}"
            logger.error(err_msg)
            return Err(err_msg)

        logger.info(f"Probe-subset VCF created: {probe_vcf}")

        # === Calculate IBS distance for probe VCF ===
        logger.info("Calculating IBS distances for probe set using PLINK...")
        probe_ibs_result = run_plink_ibs_distance(probe_vcf, tmpdir / "probe", threads=threads)

        if probe_ibs_result.is_err():
            err_msg = f"Probe VCF IBS failed: {probe_ibs_result.unwrap_err()}"
            logger.error(err_msg)
            return Err(err_msg)

        dist_probe, samples_probe = probe_ibs_result.unwrap()
        n_sites_probe = len(probe_positions)
        logger.info(f"Probe VCF: {len(samples_probe)} samples, {n_sites_probe} sites")
    
    # Verify samples match
    if samples_full != samples_probe:
        return Err("Sample mismatch between full and probe VCF parsing")
    
    # Calculate correlation statistics
    logger.info("Calculating correlation statistics...")
    corr_stats = calculate_distance_correlation(dist_full, dist_probe)
    
    logger.info(f"Pearson r: {corr_stats['pearson_r']:.4f}")
    logger.info(f"Spearman ρ: {corr_stats['spearman_r']:.4f}")
    logger.info(f"Manhattan distance: {corr_stats['manhattan_distance']:.4f}")
    
    # Generate plots
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    # Combined heatmap
    heatmap_path = Path(str(output_prefix) + ".ibs_heatmap.png")
    logger.info(f"Generating combined heatmap: {heatmap_path}")
    heatmap_result = generate_combined_heatmap(
        dist_full, dist_probe, samples_full, heatmap_path, corr_stats,
        pop_file=pop_file,
    )
    if heatmap_result.is_err():
        logger.warning(f"Failed to generate heatmap: {heatmap_result.unwrap_err()}")
    
    # Save distance matrices
    dist_full_path = Path(str(output_prefix) + ".ibs_full.tsv")
    dist_probe_path = Path(str(output_prefix) + ".ibs_probe.tsv")
    
    df_full = pd.DataFrame(dist_full, index=samples_full, columns=samples_full)
    df_probe = pd.DataFrame(dist_probe, index=samples_probe, columns=samples_probe)
    
    df_full.to_csv(dist_full_path, sep='\t')
    df_probe.to_csv(dist_probe_path, sep='\t')
    logger.info(f"Saved distance matrices to {dist_full_path} and {dist_probe_path}")
    
    # Save summary statistics
    summary_path = Path(str(output_prefix) + ".distance_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("1-IBS Distance Assessment Summary (PLINK)\n")
        f.write("=" * 50 + "\n\n")
        f.write("Distance metric: 1-IBS (one minus Identity By State)\n")
        f.write("  0 = identical genotypes, 1 = completely different\n\n")
        f.write(f"Full VCF: {vcf_path}\n")
        f.write(f"Probe SNPs: {snp_tsv_path}\n\n")
        f.write(f"Samples: {len(samples_full)}\n")
        f.write(f"Full VCF sites: {n_sites_full}\n")
        f.write(f"Probe sites: {n_sites_probe}\n")
        coverage_pct = (n_sites_probe / n_sites_full * 100) if n_sites_full > 0 else 0
        f.write(f"Probe coverage: {coverage_pct:.2f}%\n\n")
        f.write("Correlation Statistics (Full VCF vs Probe Set):\n")
        f.write("-" * 30 + "\n")
        f.write(f"Pearson r:         {corr_stats['pearson_r']:.6f}\n")
        f.write(f"Pearson p-value:   {corr_stats['pearson_p']:.2e}\n")
        f.write(f"Spearman ρ:        {corr_stats['spearman_r']:.6f}\n")
        f.write(f"Spearman p-value:  {corr_stats['spearman_p']:.2e}\n")
        f.write(f"Manhattan dist:    {corr_stats['manhattan_distance']:.6f}\n")
        f.write(f"Sample pairs:      {corr_stats['n_pairs']}\n")
    
    logger.info(f"Saved summary to {summary_path}")
    
    stats = {
        "n_samples": len(samples_full),
        "n_sites_full": n_sites_full,
        "n_sites_probe": n_sites_probe,
        "coverage_percent": coverage_pct,
        "correlation": corr_stats,
        "heatmap_file": str(heatmap_path) if heatmap_result.is_ok() else None,
        "dist_full_file": str(dist_full_path),
        "dist_probe_file": str(dist_probe_path),
        "summary_file": str(summary_path),
    }
    
    return Ok(stats)


def _smart_plot_limits(
    values: np.ndarray,
    n_bins: int = 50,
) -> Tuple[Tuple[float, float], np.ndarray]:
    """Auto-detect discrete vs continuous data and build appropriate bins.

    When unique values <= 60, builds bins aligned to value boundaries so
    every bar sits on an actual value (fixes GC choppy-bar issue).
    Otherwise, falls back to P1-P99-clipped limits.
    """
    if len(values) == 0:
        return (0.0, 1.0), np.linspace(0.0, 1.0, n_bins + 1)

    unique_vals = np.sort(np.unique(values))

    # Discrete-aware binning when few unique values
    if len(unique_vals) <= 60:
        if len(unique_vals) == 1:
            hw = max(abs(unique_vals[0]) * 0.1, 0.5)
            edges = np.array([unique_vals[0] - hw, unique_vals[0] + hw])
        else:
            midpoints = (unique_vals[:-1] + unique_vals[1:]) / 2
            first_gap = midpoints[0] - unique_vals[0]
            last_gap = unique_vals[-1] - midpoints[-1]
            edges = np.concatenate([
                [unique_vals[0] - first_gap],
                midpoints,
                [unique_vals[-1] + last_gap],
            ])
        return (float(edges[0]), float(edges[-1])), edges

    # Standard continuous binning for high-cardinality data
    xmin = float(np.percentile(values, 1))
    xmax = float(np.percentile(values, 99))
    if xmax <= xmin:
        xmax = float(np.max(values))
    if xmax <= xmin:
        xmax = xmin + 1.0
    span = xmax - xmin
    xmin = xmin - span * 0.05
    xmax = xmax + span * 0.05
    return (xmin, xmax), np.linspace(xmin, xmax, n_bins + 1)


def generate_assessment_plots(
    results: List[AssessmentResult],
    tags: List[str],
    output_prefix: Path,
    compare_results: Optional[List[AssessmentResult]] = None,
    xlim_overrides: Optional[Dict[str, Tuple[float, float]]] = None,
    ylim_overrides: Optional[Dict[str, Tuple[float, float]]] = None,
    bins_overrides: Optional[Dict[str, int]] = None,
) -> Result[List[Path], str]:
    """
    Generate distribution plots for assessment metrics.

    Args:
        results: List of assessment results
        tags: List of metrics to plot
        output_prefix: Output file prefix
        compare_results: Optional second result set to overlay (e.g. post-filter)
        xlim_overrides: Per-tag x-axis limits, e.g. {"gc": (0.0, 1.0)}
        ylim_overrides: Per-tag y-axis limits, e.g. {"gc": (0.0, 50.0)}
        bins_overrides: Per-tag bin count, e.g. {"gc": 100}

    Returns:
        Result containing list of generated plot paths
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        return Err("Plotting requires matplotlib and seaborn. Install with: pip install eprobe[plot]")

    plot_paths = []

    # Legacy color palette
    colors = ['#83639F', '#EA7827', '#C22f2F', '#449945', '#1F70A9']

    for idx, tag in enumerate(tags):
        values = np.array([getattr(r, tag) for r in results if getattr(r, tag) is not None], dtype=float)

        if len(values) == 0:
            continue

        tag_bins = (bins_overrides or {}).get(tag, 50)
        xlim_fixed = (xlim_overrides or {}).get(tag)

        if xlim_fixed is not None:
            xlim = xlim_fixed
            bin_edges = np.linspace(xlim[0], xlim[1], tag_bins + 1)
        else:
            xlim, bin_edges = _smart_plot_limits(values, n_bins=tag_bins)

        # When overlaying, expand axis to cover both datasets (unless x-axis is fixed)
        cmp_values = np.array([], dtype=float)
        if compare_results is not None:
            cmp_values = np.array(
                [getattr(r, tag) for r in compare_results if getattr(r, tag) is not None], dtype=float
            )
            if len(cmp_values) > 0 and xlim_fixed is None:
                combined = np.concatenate([values, cmp_values])
                xlim, bin_edges = _smart_plot_limits(combined, n_bins=tag_bins)

        # --- Hairpin-specific: equal-spaced stem-bp transform ---
        tick_positions_custom = None
        tick_labels_custom = None
        plot_vals = values
        plot_cmp = cmp_values
        if tag == 'hairpin' and xlim_fixed is None:
            unique_raw = np.sort(
                np.unique(np.concatenate([values, cmp_values]))
                if len(cmp_values) > 0 else np.unique(values)
            )
            val_to_idx = {v: i for i, v in enumerate(unique_raw)}
            plot_vals = np.array([val_to_idx[v] for v in values], dtype=float)
            plot_cmp = (np.array([val_to_idx[v] for v in cmp_values], dtype=float)
                        if len(cmp_values) > 0 else cmp_values)
            n_cat = len(unique_raw)
            bin_edges = np.arange(-0.5, n_cat, 1.0)
            xlim = (-0.5, n_cat - 0.5)
            tick_positions_custom = np.arange(n_cat, dtype=float)
            # Convert scores to stem bp labels
            _log4 = np.log(4)
            nonzero = unique_raw[unique_raw > 0]
            if len(nonzero) > 0:
                min_nz = float(nonzero.min())
                labels = []
                for v in unique_raw:
                    if v == 0:
                        labels.append('<4bp')
                    else:
                        n = round(np.log(v / min_nz) / _log4) + 1 if v > min_nz * 0.5 else 1
                        labels.append(f'{n + 3}bp')
                tick_labels_custom = labels
            else:
                tick_labels_custom = [f'{v:.2f}' for v in unique_raw]

        color = colors[idx % len(colors)]
        sns.set(style='darkgrid')
        fig, ax = plt.subplots(figsize=(8, 6))

        overlay = compare_results is not None and len(cmp_values) > 0
        pre_label = f'Input (N={len(values):,})' if overlay else None
        sns.histplot(plot_vals, alpha=0.65, stat='percent',
                     bins=bin_edges, edgecolor=color, color=color, element='step',
                     line_kws={'linewidth': 2}, ax=ax, label=pre_label)

        if overlay:
            sns.histplot(plot_cmp, alpha=0.35, stat='percent',
                         bins=bin_edges, edgecolor=color, color=color, element='step',
                         line_kws={'linewidth': 1.5, 'linestyle': '--'}, ax=ax,
                         label=f'Compare (N={len(cmp_values):,})')
            ax.legend(fontsize=14)

        ax.set_xlim(xlim)
        ylim_fixed = (ylim_overrides or {}).get(tag)
        if ylim_fixed is not None:
            ax.set_ylim(ylim_fixed)
        ax.set_title('Distribution Plot', fontsize=22, fontname='Arial', pad=20)
        ax.set_xlabel(
            'Estimated stem length (bp)' if tick_positions_custom is not None else 'Values',
            fontsize=20, fontname='Arial', labelpad=10,
        )
        ax.set_ylabel('Percent (%)', fontsize=18, fontname='Arial', labelpad=10)
        ax.tick_params(axis='both', labelsize=18)
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        if tick_positions_custom is not None:
            ax.set_xticks(tick_positions_custom)
            ax.set_xticklabels(tick_labels_custom, fontsize=15, rotation=30, ha='right',
                               fontname='Arial')
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontname('Arial')

        plt.tight_layout()

        plot_path = Path(str(output_prefix) + f".{tag}_dist.jpg")
        fig.savefig(plot_path, dpi=600)
        plt.close(fig)

        plot_paths.append(plot_path)
        logger.info(f"Saved plot: {plot_path}")

    return Ok(plot_paths)


def run_assess(
    input_path: Path,
    output_prefix: Path,
    mode: str = "tags",
    vcf_path: Optional[Path] = None,
    reference_path: Optional[Path] = None,
    tags: Optional[List[str]] = None,
    generate_plots: bool = True,
    probe_length: int = 81,
    sample_dimer: int = 1000,
    max_vcf_sites: int = 100000,
    pop_file: Optional[Path] = None,
    n_samples_per_pop: int = 5,
    specified_samples: Optional[List[str]] = None,
    projection: Optional[str] = None,
    pops: Optional[List[str]] = None,
    seed: int = 42,
    threads: int = 1,
    verbose: bool = False,
    compare_path: Optional[Path] = None,
    plot_xlim: Optional[Dict[str, Tuple[float, float]]] = None,
    plot_ylim: Optional[Dict[str, Tuple[float, float]]] = None,
    plot_bins: Optional[Dict[str, int]] = None,
) -> Result[Dict[str, Any], str]:
    """
    Assess quality of probe set.
    
    Main entry point for the assess command. Supports multiple assessment modes:
    - tags: Calculate biophysical metrics and distributions
    - distance: Compare IBS distances between full VCF and probe-covered SNPs (uses PLINK)
    - pca: Compare Principal Component Analysis between full VCF and probe-covered SNPs (uses PLINK)
    - sfs: Compare multi-population joint Site Frequency Spectrum using easySFS
    - all: Run all assessments
    
    Args:
        input_path: Input file (FASTA or TSV for tags, TSV for distance/pca/sfs)
        output_prefix: Output file prefix
        mode: Assessment mode (tags, distance, pca, sfs, all)
        vcf_path: Original VCF file (required for distance/pca/sfs mode)
        reference_path: Reference FASTA (optional, for generating probe sequences)
        tags: List of metrics to calculate (default: gc, tm, complexity, hairpin)
        generate_plots: Generate distribution plots
        probe_length: Probe length when generating sequences from TSV (default: 81)
        sample_dimer: Number of pairs to sample for dimer calculation
        max_vcf_sites: Maximum sites to load from full VCF
        pop_file: Population assignment file for SFS/PCA mode (optional for PCA coloring)
        n_samples_per_pop: Number of samples per population for SFS
        specified_samples: Specific samples to use (overrides random subsampling)
        projection: easySFS projection values (e.g., "10,10,10")
        seed: Random seed for sample subsampling
        threads: Number of threads
        verbose: Enable verbose logging
        compare_path: Second TSV to overlay against input in tags mode
        plot_xlim: Per-tag x-axis limits, e.g. {"gc": (0.0, 1.0)}
        plot_ylim: Per-tag y-axis limits, e.g. {"gc": (0.0, 50.0)}
        plot_bins: Per-tag bin count, e.g. {"gc": 100}
        
    Returns:
        Result containing assessment statistics
    """
    # Configure logging
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        logger.addHandler(handler)
    
    if verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    
    mode = mode.lower()
    all_stats = {}
    
    # === DISTANCE MODE ===
    if mode in ["distance", "all"]:
        if vcf_path is None:
            return Err("Distance mode requires --vcf parameter")
        
        logger.info("=" * 50)
        logger.info("Running IBS Distance Assessment")
        logger.info("=" * 50)
        
        dist_result = run_distance_assessment(
            snp_tsv_path=input_path,
            vcf_path=vcf_path,
            output_prefix=output_prefix,
            max_sites=max_vcf_sites,
            pop_file=pop_file,
            n_samples=n_samples_per_pop,
            seed=seed,
            threads=threads,
            verbose=verbose,
        )
        
        if dist_result.is_err():
            return Err(f"Distance assessment failed: {dist_result.unwrap_err()}")
        
        all_stats["distance"] = dist_result.unwrap()
        logger.info("Distance assessment complete")
    
    # === SFS MODE ===
    if mode in ["sfs", "all"]:
        if vcf_path is None:
            return Err("SFS mode requires --vcf parameter")
        if pop_file is None:
            return Err("SFS mode requires --pop_file parameter")
        
        logger.info("=" * 50)
        logger.info("Running SFS Assessment (easySFS)")
        logger.info("=" * 50)
        
        sfs_result = run_sfs_assessment(
            snp_tsv_path=input_path,
            vcf_path=vcf_path,
            output_prefix=output_prefix,
            pop_file=pop_file,
            n_samples_per_pop=n_samples_per_pop,
            specified_samples=specified_samples,
            projection=projection,
            pops=pops,
            seed=seed,
            threads=threads,
            verbose=verbose,
        )
        
        if sfs_result.is_err():
            return Err(f"SFS assessment failed: {sfs_result.unwrap_err()}")
        
        all_stats["sfs"] = sfs_result.unwrap()
        logger.info("SFS assessment complete")
    
    # === PCA MODE ===
    if mode in ["pca", "all"]:
        if vcf_path is None:
            return Err("PCA mode requires --vcf parameter")
        
        logger.info("=" * 50)
        logger.info("Running PCA Assessment (PLINK)")
        logger.info("=" * 50)
        
        pca_result = run_pca_assessment(
            snp_tsv_path=input_path,
            vcf_path=vcf_path,
            output_prefix=output_prefix,
            pop_file=pop_file,
            n_samples=n_samples_per_pop,
            seed=seed,
            threads=threads,
            verbose=verbose,
        )
        
        if pca_result.is_err():
            return Err(f"PCA assessment failed: {pca_result.unwrap_err()}")
        
        all_stats["pca"] = pca_result.unwrap()
        logger.info("PCA assessment complete")
    
    # === TAGS MODE ===
    if mode in ["tags", "all"]:
        logger.info("=" * 50)
        logger.info("Running Biophysical Tags Assessment")
        logger.info("=" * 50)
        
        tags_result = run_tags_assessment(
            input_path=input_path,
            output_prefix=output_prefix,
            reference_path=reference_path,
            tags=tags,
            generate_plots=generate_plots,
            sample_dimer=sample_dimer,
            compare_path=compare_path,
            plot_xlim=plot_xlim,
            plot_ylim=plot_ylim,
            plot_bins=plot_bins,
            probe_length=probe_length,
            verbose=verbose,
        )
        
        if tags_result.is_err():
            return Err(f"Tags assessment failed: {tags_result.unwrap_err()}")
        
        all_stats["tags"] = tags_result.unwrap()
        logger.info("Tags assessment complete")
    
    if not all_stats:
        return Err(f"Unknown mode: {mode}. Use: tags, distance, sfs, all")
    
    return Ok(all_stats)


def run_tags_assessment(
    input_path: Path,
    output_prefix: Path,
    reference_path: Optional[Path] = None,
    tags: Optional[List[str]] = None,
    generate_plots: bool = True,
    sample_dimer: int = 1000,
    compare_path: Optional[Path] = None,
    plot_xlim: Optional[Dict[str, Tuple[float, float]]] = None,
    plot_ylim: Optional[Dict[str, Tuple[float, float]]] = None,
    plot_bins: Optional[Dict[str, int]] = None,
    probe_length: int = 81,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Run biophysical tags assessment.
    
    Calculates GC, Tm, complexity, hairpin, dimer scores for probes.
    
    Args:
        input_path: Input file (FASTA or TSV)
        output_prefix: Output file prefix
        reference_path: Reference FASTA (for generating sequences from TSV)
        tags: List of metrics to calculate
        generate_plots: Generate distribution plots
        sample_dimer: Number of pairs for dimer sampling
        compare_path: Optional second TSV to overlay distributions against
        plot_xlim: Per-tag x-axis limits override. Preset defaults: gc=(20,80), tm=(60,80), complexity=(0,2)
        plot_ylim: Per-tag y-axis limits override, e.g. {"gc": (0.0, 50.0)}
        plot_bins: Per-tag bin count override, e.g. {"gc": 100}
        verbose: Enable verbose logging
        
    Returns:
        Result containing assessment statistics
    """
    # Default to all tags except dimer (expensive)
    if tags is None:
        tags = ["gc", "tm", "complexity", "hairpin"]
    
    # Normalize tag names
    tags = [t.lower() for t in tags]
    
    # Validate tags
    invalid_tags = [t for t in tags if t not in AVAILABLE_TAGS]
    if invalid_tags:
        return Err(f"Invalid tags: {invalid_tags}. Available: {list(AVAILABLE_TAGS.keys())}")
    
    logger.info(f"Starting tags assessment from {input_path}")
    logger.info(f"Tags: {tags}")
    
    # Load input
    suffix = input_path.suffix.lower()
    
    if suffix in ['.fa', '.fasta', '.fna']:
        # Load from FASTA
        fasta_result = read_fasta(input_path)
        if fasta_result.is_err():
            return Err(f"Failed to read FASTA: {fasta_result.unwrap_err()}")
        sequences = fasta_result.unwrap()
    elif suffix == '.tsv':
        # Always generate probe sequences fresh from reference + coordinates.
        # Never reuse biophysical columns from filter step — assess is an
        # independent quality check that must recalculate from scratch.
        snp_df_result = SNPDataFrame.from_tsv(input_path)
        if snp_df_result.is_err():
            return Err(f"Failed to read TSV: {snp_df_result.unwrap_err()}")

        snp_df = snp_df_result.unwrap()
        df = snp_df.df

        if reference_path is None:
            return Err("TSV input requires --reference to generate probe sequences")

        ref_result = read_fasta(reference_path)
        if ref_result.is_err():
            return Err(f"Failed to read reference: {ref_result.unwrap_err()}")

        reference = ref_result.unwrap()

        # Generate probe sequences matching filter.py logic with asymmetric flanking for even lengths:
        # For ODD probe_length (e.g., 81): centered flanking
        #   flank_size=40 → 40 + 1(ref) + 40 = 81 ✓ (perfect centering)
        # For EVEN probe_length (e.g., 80): SNP shifted left by 1 bp
        #   flank_size=40 → 39 + 1(ref@left) + 40 = 80 ✓ (precise, SNP not at center)
        # This ensures all probes meet exact length requirements without truncation.
        is_odd = probe_length % 2 == 1
        if is_odd:
            left_flank_len = (probe_length - 1) // 2
            right_flank_len = (probe_length - 1) // 2
        else:
            left_flank_len = probe_length // 2 - 1
            right_flank_len = probe_length // 2
        
        has_ref_col = "ref" in df.columns

        sequences = {}
        skipped = 0
        skipped_by_chrom = 0
        skipped_by_boundary = 0
        observed_chroms = set()
        boundary_examples = []  # Store examples of boundary-clipped probes
        
        # Calculate expected flank sizes for diagnostics
        is_even = probe_length % 2 == 0
        expected_pattern = f"{'asymmetric' if is_even else 'centered'} ({left_flank_len}+1+{right_flank_len})"
        
        for _, row in df.iterrows():
            chrom = str(row['chr'])
            pos = int(row['pos'])
            observed_chroms.add(chrom)

            if chrom not in reference:
                skipped += 1
                skipped_by_chrom += 1
                continue

            chrom_seq = reference[chrom]
            chrom_len = len(chrom_seq)
            left_end = pos - 1       # 0-based position of SNP
            right_start = pos        # 0-based position after SNP
            
            start_pos = max(0, left_end - left_flank_len)
            end_pos = min(chrom_len, right_start + right_flank_len)

            left_flank = chrom_seq[start_pos:left_end] if left_end > start_pos else ""
            right_flank = chrom_seq[right_start:end_pos] if end_pos > right_start else ""

            ref_base = str(row['ref']).upper() if has_ref_col else chrom_seq[pos - 1].upper()
            probe_seq = left_flank + ref_base + right_flank

            if len(probe_seq) != probe_length:
                skipped += 1
                skipped_by_boundary += 1
                
                # Collect examples of boundary clipping
                if len(boundary_examples) < 3:
                    boundary_examples.append({
                        'chrom': chrom,
                        'pos': pos,
                        'chrom_len': chrom_len,
                        'expected_len': probe_length,
                        'actual_len': len(probe_seq),
                        'left_flank_available': len(left_flank),
                        'right_flank_available': len(right_flank),
                    })
                continue

            probe_id = f"{chrom}_{pos}"
            sequences[probe_id] = probe_seq.upper()

        # Verify all sequences have the expected length (for data quality)
        if sequences:
            seq_lengths = {len(seq) for seq in sequences.values()}
            if len(seq_lengths) > 1:
                logger.warning(f"UNEXPECTED: Generated {len(seq_lengths)} different probe lengths!")
                for length in sorted(seq_lengths):
                    count = sum(1 for seq in sequences.values() if len(seq) == length)
                    logger.warning(f"  - {length}bp: {count} probes")
            elif list(seq_lengths)[0] != probe_length:
                logger.warning(f"UNEXPECTED: All {len(sequences)} probes are {list(seq_lengths)[0]}bp, not {probe_length}bp!")

        # Diagnose if chromosome mismatch caused all probes to be skipped
        if len(sequences) == 0 and len(df) > 0:
            ref_chroms = set(reference.keys())
            missing_chroms = observed_chroms - ref_chroms
            if missing_chroms:
                logger.error(f"CHROMOSOME MISMATCH: {len(missing_chroms)} missing in reference")
                logger.error(f"  TSV chromosomes: {sorted(observed_chroms)}")
                logger.error(f"  Reference chromosomes: {sorted(ref_chroms)}")
                logger.error(f"  Missing in reference: {sorted(missing_chroms)}")
            
            # If chromosome names match but all probes clipped, diagnose flanking issue
            if not missing_chroms and skipped_by_boundary > 0:
                logger.error(f"ALL PROBES CLIPPED: All {skipped_by_boundary} SNPs skipped due to boundaries")
                logger.error(f"Probe configuration: length={probe_length} ({expected_pattern})")
                logger.error(f"  Required flank: left {left_flank_len}bp, right {right_flank_len}bp")
                logger.error(f"Examples of clipped SNPs:")
                for ex in boundary_examples:
                    logger.error(f"  • {ex['chrom']}:{ex['pos']} (chromosome {ex['chrom_len']}bp)")
                    logger.error(f"    Expected: {ex['expected_len']}bp, got: {ex['actual_len']}bp")
                    logger.error(f"    Left flank available: {ex['left_flank_available']}bp, Right: {ex['right_flank_available']}bp")
                
                logger.error(f"\nTo fix this:")
                logger.error(f"  1. Use a shorter probe_length (e.g., -l 50 instead of 60)")
                logger.error(f"  2. Pre-filter TSV to exclude SNPs near boundaries")
                logger.error(f"  3. Check if positions are 0-based (should be 1-based)")
                logger.error(f"\nFor reference, with probe_length={probe_length}:")
                logger.error(f"  - SNPs in first {left_flank_len}bp of chromosome will be clipped on left")
                logger.error(f"  - SNPs in last {right_flank_len}bp of chromosome will be clipped on right")
        
        if skipped > 0:
            logger.warning(f"Skipped {skipped} probes total:")
            logger.warning(f"  - {skipped_by_chrom} missing chromosome in reference")
            logger.warning(f"  - {skipped_by_boundary} boundary-clipped (insufficient flanking)")
        logger.info(f"Generated {len(sequences)} probe sequences (length={probe_length}) from reference")
    else:
        return Err(f"Unsupported input format: {suffix}")
    
    logger.info(f"Loaded {len(sequences)} sequences")
    
    # Calculate individual probe metrics
    logger.info("Calculating probe metrics...")
    tags_without_dimer = [t for t in tags if t != "dimer"]
    results = assess_probes(sequences, tags_without_dimer)
    
    # Calculate dimer scores separately (sampling-based)
    if "dimer" in tags:
        logger.info(f"Calculating dimer scores (sampling {sample_dimer} pairs)...")
        dimer_scores = calculate_dimer_scores(sequences, sample_dimer)
        
        for result in results:
            result.dimer = dimer_scores.get(result.probe_id, 0.0)
    
    # Calculate summary statistics
    summary = calculate_summary_statistics(results, tags)
    
    # Save detailed results to TSV
    output_path = Path(str(output_prefix) + ".tags_stats.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Convert results to DataFrame
    data = []
    for r in results:
        row = {"probe_id": r.probe_id}
        for tag in tags:
            row[tag] = getattr(r, tag)
        data.append(row)
    
    df = pd.DataFrame(data)
    df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Saved detailed stats to {output_path}")
    
    # Save summary
    summary_path = Path(str(output_prefix) + ".tags_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("eProbe Biophysical Tags Assessment Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Input: {input_path}\n")
        f.write(f"Total probes: {len(sequences)}\n\n")
        
        for tag, stats in summary.items():
            f.write(f"{AVAILABLE_TAGS.get(tag, tag)}:\n")
            f.write(f"  Mean:   {stats['mean']:.3f}\n")
            f.write(f"  Std:    {stats['std']:.3f}\n")
            f.write(f"  Min:    {stats['min']:.3f}\n")
            f.write(f"  Max:    {stats['max']:.3f}\n")
            f.write(f"  Median: {stats['median']:.3f}\n\n")
    
    logger.info(f"Saved summary to {summary_path}")
    
    # Apply default plot preset ranges (user overrides take precedence)
    default_xlim = {
        "gc": (20.0, 80.0),
        "tm": (60.0, 80.0),
        "complexity": (0.0, 2.0),
    }
    xlim_final = {**default_xlim, **(plot_xlim or {})}
    ylim_final = plot_ylim or {}
    bins_final = plot_bins or {}
    
    # Generate plots
    plot_paths = []
    if generate_plots:
        plot_result = generate_assessment_plots(
            results, tags, output_prefix,
            xlim_overrides=xlim_final,
            ylim_overrides=ylim_final,
            bins_overrides=bins_final,
        )
        if plot_result.is_ok():
            plot_paths = plot_result.unwrap()
    
    stats = {
        "probe_count": len(sequences),
        "tags": tags,
        "summary": summary,
        "stats_file": str(output_path),
        "summary_file": str(summary_path),
        "plot_files": [str(p) for p in plot_paths],
    }
    
    logger.info(f"Tags assessment complete for {len(sequences)} probes")
    
    return Ok(stats)


def run_tags_from_dataframe(
    df: pd.DataFrame,
    tags: List[str],
    output_prefix: Path,
    generate_plots: bool = True,
    compare_df: Optional[pd.DataFrame] = None,
    xlim_overrides: Optional[Dict[str, Tuple[float, float]]] = None,
    ylim_overrides: Optional[Dict[str, Tuple[float, float]]] = None,
    bins_overrides: Optional[Dict[str, int]] = None,
) -> Result[Dict[str, Any], str]:
    """
    Run tags assessment using existing biophysical columns in DataFrame.

    This is used when the TSV from filter step already has gc, tm, etc. columns.

    Args:
        df: DataFrame with biophysical columns
        tags: List of tags to analyze
        output_prefix: Output prefix
        generate_plots: Whether to generate plots
        compare_df: Optional second DataFrame to overlay in plots (e.g. post-filter)
        xlim_overrides: Per-tag x-axis limits override
        ylim_overrides: Per-tag y-axis limits override
        bins_overrides: Per-tag bin count override

    Returns:
        Result with assessment statistics
    """
    logger.info("Using existing biophysical data from DataFrame")
    
    summary = {}
    for tag in tags:
        if tag not in df.columns:
            continue
        
        values = df[tag].dropna().values
        if len(values) == 0:
            continue
        
        summary[tag] = {
            "mean": float(np.mean(values)),
            "std": float(np.std(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "median": float(np.median(values)),
            "count": len(values),
        }
    
    # Save summary
    summary_path = Path(str(output_prefix) + ".tags_summary.txt")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(summary_path, 'w') as f:
        f.write("eProbe Biophysical Tags Assessment Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total SNPs: {len(df)}\n")
        f.write(f"Tags analyzed: {tags}\n\n")
        
        for tag, stats in summary.items():
            f.write(f"{AVAILABLE_TAGS.get(tag, tag)}:\n")
            f.write(f"  Mean:   {stats['mean']:.3f}\n")
            f.write(f"  Std:    {stats['std']:.3f}\n")
            f.write(f"  Min:    {stats['min']:.3f}\n")
            f.write(f"  Max:    {stats['max']:.3f}\n")
            f.write(f"  Median: {stats['median']:.3f}\n\n")
    
    logger.info(f"Saved summary to {summary_path}")
    
    # Generate plots
    plot_paths = []
    if generate_plots:
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns

            # Legacy color palette
            colors = ['#83639F', '#EA7827', '#C22f2F', '#449945', '#1F70A9']

            for idx, tag in enumerate(tags):
                if tag not in df.columns:
                    continue

                values = df[tag].dropna().values
                if len(values) == 0:
                    continue

                tag_bins = (bins_overrides or {}).get(tag, 50)
                xlim_fixed = (xlim_overrides or {}).get(tag)

                if xlim_fixed is not None:
                    xlim = xlim_fixed
                    bin_edges = np.linspace(xlim[0], xlim[1], tag_bins + 1)
                else:
                    xlim, bin_edges = _smart_plot_limits(values, n_bins=tag_bins)

                # When overlaying, expand axis to cover both datasets (unless x-axis is fixed)
                cmp_values = np.array([], dtype=float)
                if compare_df is not None and tag in compare_df.columns:
                    cmp_values = compare_df[tag].dropna().values
                    if len(cmp_values) > 0 and xlim_fixed is None:
                        combined = np.concatenate([values, cmp_values])
                        xlim, bin_edges = _smart_plot_limits(combined, n_bins=tag_bins)

                # --- Hairpin-specific: equal-spaced stem-bp transform ---
                tick_positions_custom = None
                tick_labels_custom = None
                plot_vals = values
                plot_cmp = cmp_values
                if tag == 'hairpin' and xlim_fixed is None:
                    unique_raw = np.sort(
                        np.unique(np.concatenate([values, cmp_values]))
                        if len(cmp_values) > 0 else np.unique(values)
                    )
                    val_to_idx = {v: i for i, v in enumerate(unique_raw)}
                    plot_vals = np.array([val_to_idx[v] for v in values], dtype=float)
                    plot_cmp = (np.array([val_to_idx[v] for v in cmp_values], dtype=float)
                                if len(cmp_values) > 0 else cmp_values)
                    n_cat = len(unique_raw)
                    bin_edges = np.arange(-0.5, n_cat, 1.0)
                    xlim = (-0.5, n_cat - 0.5)
                    tick_positions_custom = np.arange(n_cat, dtype=float)
                    _log4 = np.log(4)
                    nonzero = unique_raw[unique_raw > 0]
                    if len(nonzero) > 0:
                        min_nz = float(nonzero.min())
                        labels = []
                        for v in unique_raw:
                            if v == 0:
                                labels.append('<4bp')
                            else:
                                n = round(np.log(v / min_nz) / _log4) + 1 if v > min_nz * 0.5 else 1
                                labels.append(f'{n + 3}bp')
                        tick_labels_custom = labels
                    else:
                        tick_labels_custom = [f'{v:.2f}' for v in unique_raw]

                color = colors[idx % len(colors)]
                sns.set(style='darkgrid')
                fig, ax = plt.subplots(figsize=(8, 6))

                overlay = len(cmp_values) > 0
                pre_label = f'Input (N={len(values):,})' if overlay else None
                sns.histplot(plot_vals, alpha=0.65, stat='percent',
                             bins=bin_edges, edgecolor=color, color=color, element='step',
                             line_kws={'linewidth': 2}, ax=ax, label=pre_label)

                if overlay:
                    sns.histplot(plot_cmp, alpha=0.35, stat='percent',
                                 bins=bin_edges, edgecolor=color, color=color, element='step',
                                 line_kws={'linewidth': 1.5, 'linestyle': '--'}, ax=ax,
                                 label=f'Compare (N={len(cmp_values):,})')
                    ax.legend(fontsize=14)

                ax.set_xlim(xlim)
                ylim_fixed = (ylim_overrides or {}).get(tag)
                if ylim_fixed is not None:
                    ax.set_ylim(ylim_fixed)
                ax.set_title('Distribution Plot', fontsize=22, fontname='Arial', pad=20)
                ax.set_xlabel(
                    'Estimated stem length (bp)' if tick_positions_custom is not None else 'Values',
                    fontsize=20, fontname='Arial', labelpad=10,
                )
                ax.set_ylabel('Percent (%)', fontsize=18, fontname='Arial', labelpad=10)
                ax.tick_params(axis='both', labelsize=18)
                ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
                if tick_positions_custom is not None:
                    ax.set_xticks(tick_positions_custom)
                    ax.set_xticklabels(tick_labels_custom, fontsize=15, rotation=30, ha='right',
                                       fontname='Arial')
                for label in ax.get_xticklabels() + ax.get_yticklabels():
                    label.set_fontname('Arial')

                plt.tight_layout()

                plot_path = Path(str(output_prefix) + f".{tag}_dist.jpg")
                fig.savefig(plot_path, dpi=600)
                plt.close(fig)

                plot_paths.append(str(plot_path))
                logger.info(f"Saved plot: {plot_path}")

        except ImportError:
            logger.warning("Plotting requires matplotlib and seaborn")
    
    return Ok({
        "probe_count": len(df),
        "tags": tags,
        "summary": summary,
        "summary_file": str(summary_path),
        "plot_files": plot_paths,
    })
