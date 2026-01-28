"""
Probe generation from GFF annotations.

Extracts features from GFF files, retrieves sequences from reference,
and generates tiled probes.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, Set
from dataclasses import dataclass
import re

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta
from eprobe.core.models import Probe
from eprobe.funcgen.from_fasta import tile_sequence, TilingConfig
from eprobe.funcgen.from_bed import BedRegion, extract_region_sequence

logger = logging.getLogger(__name__)


@dataclass
class GffFeature:
    """Representation of a GFF feature."""
    seqid: str
    source: str
    feature_type: str
    start: int  # 1-based
    end: int    # 1-based, inclusive
    score: Optional[float]
    strand: str
    phase: Optional[int]
    attributes: Dict[str, str]
    
    @property
    def length(self) -> int:
        return self.end - self.start + 1
    
    def get_id(self) -> Optional[str]:
        """Get feature ID from attributes."""
        return (
            self.attributes.get('ID') or 
            self.attributes.get('gene_id') or
            self.attributes.get('Name') or
            self.attributes.get('transcript_id')
        )
    
    def get_parent(self) -> Optional[str]:
        """Get parent ID from attributes."""
        return self.attributes.get('Parent')
    
    def get_gene_name(self) -> Optional[str]:
        """Get gene name from attributes."""
        return (
            self.attributes.get('gene_name') or
            self.attributes.get('gene') or
            self.attributes.get('Name')
        )


def parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    """
    Parse GFF attribute string into dictionary.
    
    Handles both GFF3 (key=value;) and GTF (key "value";) formats.
    """
    attributes = {}
    
    if not attr_string or attr_string == '.':
        return attributes
    
    # Try GFF3 format first (key=value)
    for part in attr_string.split(';'):
        part = part.strip()
        if not part:
            continue
        
        if '=' in part:
            # GFF3 format
            key, value = part.split('=', 1)
            attributes[key.strip()] = value.strip()
        elif ' ' in part:
            # GTF format
            match = re.match(r'(\S+)\s+"?([^"]+)"?', part)
            if match:
                attributes[match.group(1)] = match.group(2)
    
    return attributes


def parse_gff_file(
    gff_path: Path,
    feature_types: Optional[Set[str]] = None,
) -> Result[List[GffFeature], str]:
    """
    Parse GFF/GFF3/GTF file into list of GffFeature objects.
    
    Args:
        gff_path: Path to GFF file
        feature_types: Set of feature types to extract (None = all)
        
    Returns:
        Result containing list of GffFeature objects
    """
    features = []
    
    try:
        with open(gff_path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                
                if len(parts) < 9:
                    logger.debug(f"Line {line_num}: insufficient columns, skipping")
                    continue
                
                try:
                    feature_type = parts[2]
                    
                    # Filter by feature type
                    if feature_types and feature_type not in feature_types:
                        continue
                    
                    feature = GffFeature(
                        seqid=parts[0],
                        source=parts[1],
                        feature_type=feature_type,
                        start=int(parts[3]),
                        end=int(parts[4]),
                        score=float(parts[5]) if parts[5] != '.' else None,
                        strand=parts[6],
                        phase=int(parts[7]) if parts[7] != '.' else None,
                        attributes=parse_gff_attributes(parts[8]),
                    )
                    features.append(feature)
                    
                except (ValueError, IndexError) as e:
                    logger.debug(f"Line {line_num}: parse error ({e}), skipping")
                    continue
        
        return Ok(features)
        
    except Exception as e:
        return Err(f"Failed to parse GFF file: {e}")


def filter_features_by_gene_ids(
    features: List[GffFeature],
    gene_ids: Set[str],
) -> List[GffFeature]:
    """
    Filter features to only those matching specified gene IDs.
    
    Args:
        features: List of GFF features
        gene_ids: Set of gene IDs to keep
        
    Returns:
        Filtered list of features
    """
    filtered = []
    
    for feature in features:
        feature_id = feature.get_id()
        gene_name = feature.get_gene_name()
        parent = feature.get_parent()
        
        # Check if any identifier matches
        if (feature_id in gene_ids or 
            gene_name in gene_ids or
            parent in gene_ids):
            filtered.append(feature)
    
    return filtered


def merge_overlapping_features(
    features: List[GffFeature],
    by_gene: bool = True,
) -> List[BedRegion]:
    """
    Merge overlapping features into regions.
    
    For CDS features of the same gene, merges into continuous regions.
    
    Args:
        features: List of GFF features
        by_gene: Group and merge by gene ID
        
    Returns:
        List of merged BED regions
    """
    from collections import defaultdict
    
    if by_gene:
        # Group by gene
        gene_features: Dict[str, List[GffFeature]] = defaultdict(list)
        
        for f in features:
            gene_id = f.get_id() or f.get_gene_name() or f.get_parent() or f"{f.seqid}:{f.start}"
            gene_features[gene_id].append(f)
        
        regions = []
        for gene_id, gene_feats in gene_features.items():
            # Sort by start position
            gene_feats.sort(key=lambda x: x.start)
            
            # Get chromosome and strand from first feature
            chrom = gene_feats[0].seqid
            strand = gene_feats[0].strand
            
            # Merge overlapping ranges
            merged_start = gene_feats[0].start - 1  # Convert to 0-based
            merged_end = gene_feats[0].end
            
            for f in gene_feats[1:]:
                if f.start - 1 <= merged_end:
                    # Overlapping or adjacent - extend
                    merged_end = max(merged_end, f.end)
                else:
                    # Non-overlapping - save current and start new
                    regions.append(BedRegion(
                        chrom=chrom,
                        start=merged_start,
                        end=merged_end,
                        name=gene_id,
                        strand=strand,
                    ))
                    merged_start = f.start - 1
                    merged_end = f.end
            
            # Save last region
            regions.append(BedRegion(
                chrom=chrom,
                start=merged_start,
                end=merged_end,
                name=gene_id,
                strand=strand,
            ))
        
        return regions
    else:
        # No merging - convert each feature to region
        return [
            BedRegion(
                chrom=f.seqid,
                start=f.start - 1,  # Convert to 0-based
                end=f.end,
                name=f.get_id() or f"{f.seqid}:{f.start}",
                strand=f.strand,
            )
            for f in features
        ]


def run_from_gff(
    gff_path: Path,
    reference_path: Path,
    output_prefix: Path,
    probe_length: int = 81,
    step_size: int = 30,
    feature_type: str = "CDS",
    gene_ids: Optional[List[str]] = None,
    threads: int = 1,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Generate probes from GFF annotations.
    
    Main entry point for the from_gff command. Extracts features
    from GFF, retrieves sequences, and tiles probes.
    
    Args:
        gff_path: Input GFF/GFF3/GTF file
        reference_path: Reference genome FASTA
        output_prefix: Output file prefix
        probe_length: Probe length
        step_size: Sliding window step
        feature_type: GFF feature type to extract
        gene_ids: Optional list of gene IDs to extract
        threads: Number of threads
        verbose: Enable verbose logging
        
    Returns:
        Result containing generation statistics
    """
    if verbose:
        logger.setLevel(logging.DEBUG)
    
    logger.info(f"Starting probe generation from {gff_path}")
    logger.info(f"Feature type: {feature_type}")
    logger.info(f"Probe length: {probe_length}, Step: {step_size}")
    
    # Parse GFF file
    gff_result = parse_gff_file(gff_path, feature_types={feature_type})
    if gff_result.is_err():
        return Err(gff_result.unwrap_err())
    
    features = gff_result.unwrap()
    logger.info(f"Loaded {len(features)} {feature_type} features")
    
    if not features:
        return Err(f"No {feature_type} features found in GFF file")
    
    # Filter by gene IDs if specified
    if gene_ids:
        gene_id_set = set(gene_ids)
        features = filter_features_by_gene_ids(features, gene_id_set)
        logger.info(f"Filtered to {len(features)} features for {len(gene_ids)} genes")
        
        if not features:
            return Err(f"No features found for specified gene IDs")
    
    # Merge features into regions
    regions = merge_overlapping_features(features, by_gene=True)
    logger.info(f"Merged into {len(regions)} regions")
    
    # Load reference genome
    logger.info("Loading reference genome...")
    ref_result = read_fasta(reference_path)
    if ref_result.is_err():
        return Err(f"Failed to read reference: {ref_result.unwrap_err()}")
    
    reference_sequences = ref_result.unwrap()
    logger.info(f"Loaded {len(reference_sequences)} chromosomes")
    
    # Configure tiling
    tiling_config = TilingConfig(
        probe_length=probe_length,
        step_size=step_size,
    )
    
    # Extract sequences and generate probes
    all_probes: List[Probe] = []
    extraction_errors = 0
    
    for region in regions:
        # Extract sequence
        seq_result = extract_region_sequence(region, reference_sequences)
        
        if seq_result.is_err():
            logger.warning(f"Failed to extract {region.name}: {seq_result.unwrap_err()}")
            extraction_errors += 1
            continue
        
        sequence = seq_result.unwrap()
        
        # Tile probes
        probes = tile_sequence(region.name, sequence, tiling_config)
        
        # Update coordinates
        for probe in probes:
            if region.strand != '-':
                probe.start = region.start + probe.start
                probe.end = region.start + probe.end
            probe.chrom = region.chrom
        
        all_probes.extend(probes)
    
    logger.info(f"Generated {len(all_probes)} probes from {len(regions) - extraction_errors} genes")
    
    # Save outputs
    fasta_output = Path(str(output_prefix) + ".probes.fa")
    fasta_output.parent.mkdir(parents=True, exist_ok=True)
    
    probe_sequences = {p.probe_id: p.sequence for p in all_probes}
    write_result = write_fasta(probe_sequences, fasta_output)
    if write_result.is_err():
        return Err(f"Failed to save FASTA: {write_result.unwrap_err()}")
    
    # Save TSV
    tsv_output = Path(str(output_prefix) + ".probes.tsv")
    with open(tsv_output, 'w') as f:
        f.write("probe_id\tsequence\tchrom\tstart\tend\tgene\tlength\n")
        for p in all_probes:
            f.write(f"{p.probe_id}\t{p.sequence}\t{p.chrom}\t{p.start}\t{p.end}\t{p.snp_id}\t{len(p.sequence)}\n")
    
    stats = {
        "gene_count": len(regions) - extraction_errors,
        "feature_count": len(features),
        "probe_count": len(all_probes),
        "extraction_errors": extraction_errors,
        "feature_type": feature_type,
        "probe_length": probe_length,
        "step_size": step_size,
        "fasta_file": str(fasta_output),
        "tsv_file": str(tsv_output),
    }
    
    return Ok(stats)
