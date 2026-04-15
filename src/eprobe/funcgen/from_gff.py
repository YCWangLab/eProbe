"""
Probe generation from GFF/GTF annotations.

Workflow:
  1. Parse GFF/GTF, filter by feature type (CDS, exon, gene, mRNA, etc.)
  2. Optionally filter by gene IDs
  3. Merge overlapping features by gene into BedRegion objects
  4. Delegate to from_bed's process_regions() for probe generation
     (including optional VCF haplotyping)

Legacy equivalent: no direct legacy counterpart (new functionality)
"""

import logging
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from eprobe.core.result import Err, Ok, Result
from eprobe.funcgen.from_bed import BedRegion, process_regions

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GFF data model
# ---------------------------------------------------------------------------

@dataclass
class GffFeature:
    """Representation of a GFF/GTF feature line."""
    seqid: str
    source: str
    feature_type: str
    start: int          # 1-based inclusive
    end: int            # 1-based inclusive
    score: Optional[float]
    strand: str
    phase: Optional[int]
    attributes: Dict[str, str]

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def get_id(self) -> Optional[str]:
        """Get feature ID from attributes (tries multiple conventions)."""
        return (
            self.attributes.get("ID")
            or self.attributes.get("gene_id")
            or self.attributes.get("Name")
            or self.attributes.get("transcript_id")
        )

    def get_parent(self) -> Optional[str]:
        return self.attributes.get("Parent")

    def get_gene_name(self) -> Optional[str]:
        return (
            self.attributes.get("gene_name")
            or self.attributes.get("gene")
            or self.attributes.get("Name")
        )


# ---------------------------------------------------------------------------
# GFF parsing
# ---------------------------------------------------------------------------

def _parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    """
    Parse GFF attribute string into dictionary.

    Handles both GFF3 (key=value;) and GTF (key "value";) formats.
    """
    attributes: Dict[str, str] = {}

    if not attr_string or attr_string == ".":
        return attributes

    for part in attr_string.split(";"):
        part = part.strip()
        if not part:
            continue

        if "=" in part:
            # GFF3: key=value
            key, value = part.split("=", 1)
            attributes[key.strip()] = value.strip()
        elif " " in part:
            # GTF: key "value"
            match = re.match(r'(\S+)\s+"?([^"]+)"?', part)
            if match:
                attributes[match.group(1)] = match.group(2)

    return attributes


def parse_gff_file(
    gff_path: Path,
    feature_types: Optional[Set[str]] = None,
) -> Result[List[GffFeature], str]:
    """
    Parse GFF/GFF3/GTF file.

    Args:
        gff_path: Path to GFF file
        feature_types: Set of feature types to keep (None = all)

    Returns:
        Ok(list of GffFeature), Err(message) on failure
    """
    features: List[GffFeature] = []

    try:
        with open(gff_path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) < 9:
                    continue

                try:
                    feature_type = parts[2]

                    if feature_types and feature_type not in feature_types:
                        continue

                    features.append(GffFeature(
                        seqid=parts[0],
                        source=parts[1],
                        feature_type=feature_type,
                        start=int(parts[3]),
                        end=int(parts[4]),
                        score=float(parts[5]) if parts[5] != "." else None,
                        strand=parts[6],
                        phase=int(parts[7]) if parts[7] != "." else None,
                        attributes=_parse_gff_attributes(parts[8]),
                    ))
                except (ValueError, IndexError) as e:
                    logger.debug(f"GFF line {line_num}: parse error ({e})")

        return Ok(features)

    except Exception as e:
        return Err(f"Failed to parse GFF file: {e}")


def filter_features_by_gene_ids(
    features: List[GffFeature],
    gene_ids: Set[str],
) -> List[GffFeature]:
    """Filter features matching specified gene IDs."""
    filtered = []
    for f in features:
        fid = f.get_id()
        gname = f.get_gene_name()
        parent = f.get_parent()
        if fid in gene_ids or gname in gene_ids or parent in gene_ids:
            filtered.append(f)
    return filtered


def features_to_bed_regions(
    features: List[GffFeature],
    merge_by_gene: bool = True,
) -> List[BedRegion]:
    """
    Convert GFF features to BedRegion objects.

    When merge_by_gene=True, overlapping/adjacent features of the same gene
    are merged into continuous regions.

    Args:
        features: List of GFF features
        merge_by_gene: Group and merge by gene identifier

    Returns:
        List of BedRegion objects (0-based half-open coordinates)
    """
    if not merge_by_gene:
        return [
            BedRegion(
                chrom=f.seqid,
                start=f.start - 1,  # GFF 1-based → BED 0-based
                end=f.end,
                name=f.get_id() or f"{f.seqid}:{f.start}",
                strand=f.strand,
            )
            for f in features
        ]

    # Group by gene — prioritize Parent for sub-features (CDS, exon, mRNA)
    gene_features: Dict[str, List[GffFeature]] = defaultdict(list)
    for f in features:
        # For sub-features (CDS, exon, mRNA), use Parent gene ID
        # For top-level features (gene), use own ID
        gene_id = (
            f.get_parent() or f.get_gene_name() or f.get_id()
            or f"{f.seqid}:{f.start}"
        )
        gene_features[gene_id].append(f)

    regions: List[BedRegion] = []

    for gene_id, feats in gene_features.items():
        feats.sort(key=lambda x: x.start)

        chrom = feats[0].seqid
        strand = feats[0].strand

        # Merge overlapping/adjacent intervals
        merged_start = feats[0].start - 1  # 0-based
        merged_end = feats[0].end

        for f in feats[1:]:
            f_start = f.start - 1
            if f_start <= merged_end:
                # Overlapping or adjacent
                merged_end = max(merged_end, f.end)
            else:
                # Gap — save current, start new
                regions.append(BedRegion(
                    chrom=chrom, start=merged_start, end=merged_end,
                    name=gene_id, strand=strand,
                ))
                merged_start = f_start
                merged_end = f.end

        # Save last interval
        regions.append(BedRegion(
            chrom=chrom, start=merged_start, end=merged_end,
            name=gene_id, strand=strand,
        ))

    return regions


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_from_gff(
    gff_path: Path,
    reference_path: Path,
    output_prefix: Path,
    probe_length: int = 81,
    step_size: int = 30,
    feature_type: str = "CDS",
    gene_ids: Optional[List[str]] = None,
    vcf_path: Optional[Path] = None,
    phase: bool = False,
    min_freq: float = 0.05,
    variant_only: bool = False,
    aligner: str = "muscle",
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Generate probes from GFF annotations.

    Workflow:
      1. Parse GFF, filter by feature_type (and optionally gene_ids)
      2. Merge overlapping features by gene → BedRegion list
      3. Delegate to process_regions() (shared with from_bed)
         - Simple: extract ref sequences, tile
         - Haplotype: phase VCF, extract haplotypes, variant-aware probes

    Args:
        gff_path: Input GFF/GFF3/GTF file
        reference_path: Reference genome FASTA (with .fai)
        output_prefix: Output file prefix
        probe_length: Probe length (bp)
        step_size: Sliding window step (bp)
        feature_type: GFF feature type to extract (e.g., "CDS", "exon")
        gene_ids: Optional list of gene IDs to restrict to
        vcf_path: VCF for haplotype inference (optional)
        phase: Phase VCF with shapeit5 first
        min_freq: Minimum haplotype frequency
        variant_only: Only probes at variant positions
        aligner: MSA aligner
        threads: Number of threads
        verbose: Enable verbose logging

    Returns:
        Ok(stats dict) on success, Err(message) on failure

    Output files:
        {output_prefix}.probes.fasta - Probe sequences
    """
    if verbose:
        logger.setLevel(logging.DEBUG)

    logger.info(f"Starting probe generation from {gff_path}")
    logger.info(f"Feature type: {feature_type}")
    logger.info(f"Probe length: {probe_length}, Step: {step_size}")

    # Parse GFF
    # Support comma-separated feature types (e.g., "CDS,exon")
    feature_type_set = {ft.strip() for ft in feature_type.split(",")}
    gff_result = parse_gff_file(gff_path, feature_types=feature_type_set)
    if gff_result.is_err():
        return Err(gff_result.unwrap_err())

    features = gff_result.unwrap()
    logger.info(f"Loaded {len(features)} {feature_type} features from GFF")

    if not features:
        return Err(f"No {feature_type} features found in GFF file")

    # Filter by gene IDs
    if gene_ids:
        gene_id_set = set(gene_ids)
        features = filter_features_by_gene_ids(features, gene_id_set)
        logger.info(f"Filtered to {len(features)} features for {len(gene_ids)} genes")
        if not features:
            return Err("No features found for specified gene IDs")

    # Convert to BED regions (merge overlapping by gene)
    regions = features_to_bed_regions(features, merge_by_gene=True)
    logger.info(f"Merged into {len(regions)} regions")

    # Delegate to shared region processor
    result = process_regions(
        regions=regions,
        reference_path=reference_path,
        output_prefix=output_prefix,
        probe_length=probe_length,
        step_size=step_size,
        vcf_path=vcf_path,
        phase=phase,
        min_freq=min_freq,
        variant_only=variant_only,
        aligner=aligner,
        threads=threads,
        verbose=verbose,
    )

    # Enrich stats with GFF-specific info
    if result.is_ok():
        stats = result.unwrap()
        stats["feature_type"] = feature_type
        stats["feature_count"] = len(features)
        stats["gene_count"] = len(regions)
        return Ok(stats)

    return result
