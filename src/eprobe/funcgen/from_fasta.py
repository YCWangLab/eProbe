"""
Probe generation from FASTA sequences.

Entry points:
  - Simple mode: tile each input sequence with sliding window
  - Haplotype mode: group sequences by gene (name separator),
    align alleles, generate haplotype-aware probes

Legacy equivalent: Get_probe_from_seq.py + Get_probe_from_allele.py (seq mode)
"""

import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from eprobe.core.fasta import read_fasta, write_fasta
from eprobe.core.models import Probe
from eprobe.core.result import Err, Ok, Result
from eprobe.funcgen.haplotype import (
    generate_haplotype_probes,
    tile_sequence,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Haplotype grouping from FASTA sequence IDs
# ---------------------------------------------------------------------------

def parse_haplotype_id(
    seq_id: str,
    separator: str = "_",
) -> Tuple[str, Optional[str]]:
    """
    Parse sequence ID to extract gene and allele components.

    Convention: GeneID{sep}AlleleNum  (e.g., "BRCA1_1", "BRCA1_2")
    If no allele suffix found, returns (seq_id, None).

    Args:
        seq_id: Sequence identifier
        separator: Character separating gene from allele number

    Returns:
        (gene_id, allele_id or None)
    """
    if separator in seq_id:
        parts = seq_id.rsplit(separator, 1)
        if len(parts) == 2 and parts[1].isdigit():
            return parts[0], parts[1]
    return seq_id, None


def group_by_gene(
    sequences: Dict[str, str],
    separator: str = "_",
) -> Dict[str, Dict[str, str]]:
    """
    Group FASTA sequences by gene ID.

    Args:
        sequences: seq_id -> sequence
        separator: Separator for gene_allele naming

    Returns:
        gene_id -> {allele_label: sequence}
    """
    groups: Dict[str, Dict[str, str]] = defaultdict(dict)

    for seq_id, sequence in sequences.items():
        gene_id, allele_id = parse_haplotype_id(seq_id, separator)

        if allele_id is None:
            # No allele suffix — treat entire seq_id as gene, allele "H1"
            allele_label = "H1"
        else:
            allele_label = f"H{allele_id}"

        groups[gene_id][allele_label] = sequence

    return dict(groups)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_from_fasta(
    fasta_path: Path,
    output_prefix: Path,
    probe_length: int = 81,
    step_size: int = 30,
    haplotyping: bool = False,
    haplotype_sep: str = "_",
    variant_only: bool = False,
    aligner: str = "muscle",
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Generate probes from FASTA sequences.

    Simple mode:
      Tiles probes across each input sequence using a sliding window.

    Haplotype mode (--haplotyping):
      Groups sequences by gene ID (using separator convention),
      aligns alleles of the same gene, and generates probes that:
        - At non-variant positions: 1 shared probe
        - At variant positions: 1 probe per unique haplotype
      With --variant_only: skips non-variant positions entirely.

    Args:
        fasta_path: Input FASTA file
        output_prefix: Output file prefix
        probe_length: Probe length (bp)
        step_size: Sliding window step (bp)
        haplotyping: Enable haplotype-aware processing
        haplotype_sep: Separator for gene_allele naming (e.g., "_")
        variant_only: Only generate probes at variant positions
        aligner: MSA tool ('muscle', 'clustalo', 'mafft')
        verbose: Enable verbose logging

    Returns:
        Ok(stats dict) on success, Err(message) on failure

    Output files:
        {output_prefix}.probes.fa   - Probe sequences in FASTA format
        {output_prefix}.probes.tsv  - Probe metadata table
    """
    if verbose:
        logger.setLevel(logging.DEBUG)

    logger.info(f"Starting probe generation from {fasta_path}")
    logger.info(f"Probe length: {probe_length}, Step: {step_size}")

    # Load sequences
    fasta_result = read_fasta(fasta_path)
    if fasta_result.is_err():
        return Err(f"Failed to read FASTA: {fasta_result.unwrap_err()}")

    sequences = fasta_result.unwrap()
    logger.info(f"Loaded {len(sequences)} sequences")

    if not sequences:
        return Err("No sequences found in input FASTA")

    # Generate probes
    all_probes: List[Probe] = []
    n_genes = 0
    n_multi_allele = 0

    if haplotyping:
        # Haplotype-aware mode
        logger.info(f"Haplotype mode: separator='{haplotype_sep}', variant_only={variant_only}")
        groups = group_by_gene(sequences, haplotype_sep)
        n_genes = len(groups)

        for gene_id, alleles in groups.items():
            if len(alleles) > 1:
                n_multi_allele += 1
                logger.info(
                    f"  {gene_id}: {len(alleles)} alleles "
                    f"({', '.join(f'{k}={len(v)}bp' for k, v in alleles.items())})"
                )

            result = generate_haplotype_probes(
                region_id=gene_id,
                alleles=alleles,
                probe_length=probe_length,
                step_size=step_size,
                variant_only=variant_only,
                aligner=aligner,
            )

            if result.is_err():
                logger.warning(f"  {gene_id}: probe generation failed: {result.unwrap_err()}")
                continue

            probes = result.unwrap()
            all_probes.extend(probes)
            logger.debug(f"  {gene_id}: {len(probes)} probes")

        logger.info(
            f"Processed {n_genes} genes ({n_multi_allele} multi-allele)"
        )
    else:
        # Simple tiling mode
        n_genes = len(sequences)
        for seq_id, sequence in sequences.items():
            probes = tile_sequence(seq_id, sequence, probe_length, step_size)
            all_probes.extend(probes)

    logger.info(f"Generated {len(all_probes)} probes total")

    if not all_probes:
        return Err("No probes generated (sequences may be too short or all filtered)")

    # --- Save FASTA output ---
    fasta_output = Path(str(output_prefix) + ".probes.fa")
    fasta_output.parent.mkdir(parents=True, exist_ok=True)

    probe_seqs = {p.id: p.sequence for p in all_probes}
    wr = write_fasta(probe_seqs, fasta_output)
    if wr.is_err():
        return Err(f"Failed to write FASTA: {wr.unwrap_err()}")

    logger.info(f"Saved probes to {fasta_output}")

    # --- Save TSV metadata ---
    tsv_output = Path(str(output_prefix) + ".probes.tsv")
    with open(tsv_output, "w") as f:
        f.write("probe_id\tsequence\tregion\tchrom\tstart\tend\tlength\n")
        for p in all_probes:
            f.write(
                f"{p.id}\t{p.sequence}\t{p.source_chrom}\t"
                f"{p.source_chrom}\t{p.source_start}\t{p.source_end}\t"
                f"{len(p.sequence)}\n"
            )

    stats = {
        "gene_count": n_genes,
        "multi_allele_count": n_multi_allele,
        "probe_count": len(all_probes),
        "probe_length": probe_length,
        "step_size": step_size,
        "haplotyping": haplotyping,
        "variant_only": variant_only,
        "fasta_file": str(fasta_output),
        "tsv_file": str(tsv_output),
    }

    return Ok(stats)
