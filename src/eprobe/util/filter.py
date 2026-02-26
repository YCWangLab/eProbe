"""
Standalone FASTA-level filtering using popgen filter stages.

Applies any popgen filter step (BG, AC, TX, biophysical) directly to
a FASTA file — no SNP TSV or reference genome required.

Design: wraps existing popgen/filter.py functions via a lightweight
FastaEntry adapter class whose .id matches FASTA headers, so all
ID-based filtering logic works unchanged.
"""

import logging
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple

from eprobe.core.result import Result, Ok, Err
from eprobe.core.fasta import read_fasta, write_fasta

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# SNP-compatible adapter for raw FASTA entries
# ---------------------------------------------------------------------------

class FastaEntry:
    """
    Lightweight SNP-compatible wrapper around a FASTA sequence.

    Provides ``id``, ``chrom``, ``pos``, ``ref``, ``alt``, ``mutation_type``
    and ``to_dict()`` — just enough interface for the popgen filter functions
    which only use ``.id`` for set-membership lookups.
    """

    __slots__ = ("_id", "chrom", "pos", "ref", "alt", "mutation_type")

    def __init__(self, seq_id: str) -> None:
        self._id = seq_id
        self.chrom = seq_id
        self.pos = 1
        self.ref = "A"
        self.alt = "T"
        self.mutation_type = "tv"

    @property
    def id(self) -> str:  # noqa: A003
        return self._id

    def to_dict(self) -> dict:
        return {
            "chr": self.chrom,
            "pos": self.pos,
            "type": self.mutation_type,
            "ref": self.ref,
            "alt": self.alt,
        }


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_util_filter(
    input_fasta: Path,
    output_prefix: Path,
    steps: List[str],
    threads: int = 1,
    # BG params
    bg_db: Optional[str] = None,
    # AC params
    ac_db: Optional[str] = None,
    ac_mode: str = "strict",
    ac_score_diff: int = 10,
    ac_min_genomes: Optional[int] = None,
    # TX params
    tx_db: Optional[str] = None,
    tx_ids: Optional[List[int]] = None,
    tx_mode: str = "lca",
    names_dmp: Optional[str] = None,
    nodes_dmp: Optional[str] = None,
    acc2tax: Optional[str] = None,
    tx_min_edit: int = 0,
    tx_max_edit: int = 2,
    tx_keep_hits: int = 100,
    tx_outgroup_ids: Optional[List[int]] = None,
    tx_exclude_mode: str = "diff",
    tx_exclude_val: float = 2,
    tx_target_min_sim: float = 90.0,
    tx_target_max_nm: int = 5,
    # Biophysical params
    gc_range: Tuple[float, float] = (35.0, 65.0),
    tm_range: Tuple[float, float] = (55.0, 75.0),
    max_complexity: float = 2.0,
    max_hairpin: Optional[float] = None,
    max_dimer: Optional[float] = None,
    nn_table: str = "DNA_NN4",
    na_conc: float = 50.0,
    # General
    keep_temp: bool = False,
    verbose: bool = False,
) -> Result[Dict[str, Any], str]:
    """
    Filter FASTA sequences through one or more popgen filter stages.

    Parameters
    ----------
    input_fasta : Path
        Input FASTA file.
    output_prefix : Path
        Output prefix.  Produces ``{prefix}.filtered.fa``,
        ``{prefix}.rejected.fa``, and ``{prefix}.filter.log``.
    steps : list of str
        Filter stages to apply **in order**.
        Allowed: ``"bg"``, ``"ac"``, ``"tx"``, ``"biophysical"``.
    threads : int
        Parallel threads for external tools.

    Returns
    -------
    Result containing stats dict on success.
    """
    # Lazy imports — only pull popgen/filter when actually needed
    from eprobe.popgen.filter import (
        filter_background_noise,
        filter_accessibility,
        filter_biophysical,
        filter_taxonomy,
        filter_taxonomy_besthit,
        BiophysicalThresholds,
    )

    # ---- setup logging ----
    log_path = Path(f"{output_prefix}.filter.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(str(log_path), mode="w")
    fh.setLevel(logging.DEBUG if verbose else logging.INFO)
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)
    if verbose:
        logger.setLevel(logging.DEBUG)

    try:
        return _run_impl(
            input_fasta=input_fasta,
            output_prefix=output_prefix,
            steps=steps,
            threads=threads,
            bg_db=bg_db,
            ac_db=ac_db,
            ac_mode=ac_mode,
            ac_score_diff=ac_score_diff,
            ac_min_genomes=ac_min_genomes,
            tx_db=tx_db,
            tx_ids=tx_ids,
            tx_mode=tx_mode,
            names_dmp=names_dmp,
            nodes_dmp=nodes_dmp,
            acc2tax=acc2tax,
            tx_min_edit=tx_min_edit,
            tx_max_edit=tx_max_edit,
            tx_keep_hits=tx_keep_hits,
            tx_outgroup_ids=tx_outgroup_ids,
            tx_exclude_mode=tx_exclude_mode,
            tx_exclude_val=tx_exclude_val,
            tx_target_min_sim=tx_target_min_sim,
            tx_target_max_nm=tx_target_max_nm,
            gc_range=gc_range,
            tm_range=tm_range,
            max_complexity=max_complexity,
            max_hairpin=max_hairpin,
            max_dimer=max_dimer,
            nn_table=nn_table,
            na_conc=na_conc,
            keep_temp=keep_temp,
            verbose=verbose,
            # popgen imports
            filter_background_noise=filter_background_noise,
            filter_accessibility=filter_accessibility,
            filter_biophysical=filter_biophysical,
            filter_taxonomy=filter_taxonomy,
            filter_taxonomy_besthit=filter_taxonomy_besthit,
            BiophysicalThresholds=BiophysicalThresholds,
        )
    finally:
        logger.removeHandler(fh)
        fh.close()


# ---------------------------------------------------------------------------
# Internal implementation
# ---------------------------------------------------------------------------

def _run_impl(  # noqa: C901 — complex by nature
    *,
    input_fasta: Path,
    output_prefix: Path,
    steps: List[str],
    threads: int,
    bg_db, ac_db, ac_mode, ac_score_diff, ac_min_genomes,
    tx_db, tx_ids, tx_mode, names_dmp, nodes_dmp, acc2tax,
    tx_min_edit, tx_max_edit, tx_keep_hits,
    tx_outgroup_ids, tx_exclude_mode, tx_exclude_val,
    tx_target_min_sim, tx_target_max_nm,
    gc_range, tm_range, max_complexity, max_hairpin, max_dimer,
    nn_table, na_conc, keep_temp, verbose,
    # injected functions
    filter_background_noise,
    filter_accessibility,
    filter_biophysical,
    filter_taxonomy,
    filter_taxonomy_besthit,
    BiophysicalThresholds,
) -> Result[Dict[str, Any], str]:
    """Core implementation — separated for testability."""

    # ---- read input FASTA ----
    fasta_result = read_fasta(input_fasta)
    if fasta_result.is_err():
        return Err(f"Cannot read input FASTA: {fasta_result.unwrap_err()}")
    sequences: OrderedDict = fasta_result.unwrap()

    if not sequences:
        return Err("Input FASTA is empty")

    initial_count = len(sequences)
    logger.info(f"Loaded {initial_count} sequences from {input_fasta}")

    # Build adapter objects + probe-sequences dict
    entries = [FastaEntry(sid) for sid in sequences]
    probe_seqs: Dict[str, str] = dict(sequences)

    # Prepare a temporary working FASTA that gets rewritten after each stage
    work_dir = Path(tempfile.mkdtemp(prefix="eprobe_util_filter_"))
    work_fasta = work_dir / "probes.fa"
    write_fasta(sequences, work_fasta)

    stats: Dict[str, Any] = {
        "initial_count": initial_count,
        "steps": {},
    }

    # ---- apply requested filter stages in order ----
    valid_steps = {"bg", "ac", "tx", "biophysical"}
    step_counter: Dict[str, int] = {}
    for step in steps:
        step_lower = step.strip().lower()
        if step_lower not in valid_steps:
            return Err(f"Unknown filter step '{step}'. Choose from: {sorted(valid_steps)}")

        # Generate unique label for repeated steps
        step_counter[step_lower] = step_counter.get(step_lower, 0) + 1
        step_n = step_counter[step_lower]

        def _label(base: str) -> str:
            return base if step_n == 1 else f"{base}({step_n})"
        before = len(entries)

        # ---- Background (Kraken2) ----
        if step_lower == "bg":
            if not bg_db:
                return Err("BG filter requires --bg_db")
            db_paths = [Path(d.strip()) for d in bg_db.split(",")]
            for db_path in db_paths:
                result = filter_background_noise(
                    snps=entries,
                    db_path=db_path,
                    fasta_path=work_fasta,
                    threads=threads,
                    work_dir=work_dir,
                    keep_temp=keep_temp,
                )
                if result.is_err():
                    return Err(f"BG filter failed ({db_path}): {result.unwrap_err()}")
                entries = result.unwrap()
                # Rewrite working FASTA with survivors
                _rewrite_fasta(entries, sequences, work_fasta)
            stats["steps"][_label("BG")] = {"before": before, "after": len(entries),
                                    "removed": before - len(entries)}
            logger.info(f"BG: {before} → {len(entries)} ({before - len(entries)} removed)")

        # ---- Accessibility (Bowtie2) ----
        elif step_lower == "ac":
            if not ac_db:
                return Err("AC filter requires --ac_db")
            index_paths = [Path(d.strip()) for d in ac_db.split(",")]
            result = filter_accessibility(
                snps=entries,
                index_paths=index_paths,
                threads=threads,
                mode=ac_mode,
                score_diff_threshold=ac_score_diff,
                min_genomes=ac_min_genomes,
                fasta_path=work_fasta,
                work_dir=work_dir,
                keep_temp=keep_temp,
            )
            if result.is_err():
                return Err(f"AC filter failed: {result.unwrap_err()}")
            entries = result.unwrap()
            _rewrite_fasta(entries, sequences, work_fasta)
            stats["steps"][_label("AC")] = {"before": before, "after": len(entries),
                                    "removed": before - len(entries)}
            logger.info(f"AC: {before} → {len(entries)} ({before - len(entries)} removed)")

        # ---- Taxonomic (Bowtie2 + LCA / best-hit) ----
        elif step_lower == "tx":
            if not tx_db:
                return Err("TX filter requires --tx_db")
            if not tx_ids:
                return Err("TX filter requires --tx_taxid or --tx_target_names")
            index_paths = [Path(d.strip()) for d in tx_db.split(",")]

            if tx_mode == "besthit":
                if not acc2tax or not names_dmp or not nodes_dmp:
                    return Err("TX best-hit mode requires --tx_acc2tax, --tx_names, --tx_nodes")
                result = filter_taxonomy_besthit(
                    snps=entries,
                    index_paths=index_paths,
                    names_dmp=Path(names_dmp),
                    nodes_dmp=Path(nodes_dmp),
                    acc2tax=Path(acc2tax),
                    target_taxids=tx_ids,
                    fasta_path=work_fasta,
                    outgroup_taxids=tx_outgroup_ids,
                    exclude_mode=tx_exclude_mode,
                    exclude_val=tx_exclude_val,
                    target_min_sim=tx_target_min_sim,
                    target_max_nm=tx_target_max_nm,
                    keep_hits=tx_keep_hits,
                    threads=threads,
                    work_dir=work_dir,
                    keep_temp=keep_temp,
                )
                if result.is_err():
                    return Err(f"TX(besthit) filter failed: {result.unwrap_err()}")
                entries, tx_stats = result.unwrap()
                stats["steps"][_label("TX(besthit)")] = {
                    "before": before, "after": len(entries),
                    "removed": before - len(entries), "details": tx_stats,
                }
            else:
                # LCA mode
                if not names_dmp or not nodes_dmp or not acc2tax:
                    return Err("TX LCA mode requires --tx_names, --tx_nodes, --tx_acc2tax")
                result = filter_taxonomy(
                    snps=entries,
                    index_paths=index_paths,
                    names_dmp=Path(names_dmp),
                    nodes_dmp=Path(nodes_dmp),
                    acc2tax=Path(acc2tax),
                    target_taxids=tx_ids,
                    fasta_path=work_fasta,
                    min_edit=tx_min_edit,
                    max_edit=tx_max_edit,
                    keep_hits=tx_keep_hits,
                    threads=threads,
                    work_dir=work_dir,
                    keep_temp=keep_temp,
                )
                if result.is_err():
                    return Err(f"TX(lca) filter failed: {result.unwrap_err()}")
                entries, taxid_counts = result.unwrap()
                stats["steps"][_label("TX(lca)")] = {
                    "before": before, "after": len(entries),
                    "removed": before - len(entries), "taxid_counts": taxid_counts,
                }

            _rewrite_fasta(entries, sequences, work_fasta)
            mode_label = "TX(besthit)" if tx_mode == "besthit" else "TX(lca)"
            logger.info(f"{mode_label}: {before} → {len(entries)} ({before - len(entries)} removed)")

        # ---- Biophysical ----
        elif step_lower == "biophysical":
            # Rebuild probe_seqs from current survivors
            surviving_ids = {e.id for e in entries}
            current_seqs = {k: v for k, v in sequences.items() if k in surviving_ids}

            thresholds = BiophysicalThresholds(
                gc_min=gc_range[0], gc_max=gc_range[1],
                tm_min=tm_range[0], tm_max=tm_range[1],
                complexity_max=max_complexity,
                hairpin_max=max_hairpin,
                dimer_max=max_dimer,
                nn_table=nn_table,
                na_conc=na_conc,
            )
            result = filter_biophysical(
                snps=entries,
                thresholds=thresholds,
                probe_sequences=current_seqs,
                threads=threads,
            )
            if result.is_err():
                return Err(f"Biophysical filter failed: {result.unwrap_err()}")
            entries, biophys_stats = result.unwrap()
            _rewrite_fasta(entries, sequences, work_fasta)
            stats["steps"][_label("Biophysical")] = {
                "before": before, "after": len(entries),
                "removed": before - len(entries), "details": biophys_stats,
            }
            logger.info(f"Biophysical: {before} → {len(entries)} ({before - len(entries)} removed)")

    # ---- write outputs ----
    surviving_ids = {e.id for e in entries}
    final_seqs = OrderedDict((k, v) for k, v in sequences.items() if k in surviving_ids)
    rejected_seqs = OrderedDict((k, v) for k, v in sequences.items() if k not in surviving_ids)

    out_filtered = Path(f"{output_prefix}.filtered.fa")
    out_rejected = Path(f"{output_prefix}.rejected.fa")

    write_fasta(final_seqs, out_filtered)
    logger.info(f"Wrote {len(final_seqs)} filtered sequences to {out_filtered}")

    if rejected_seqs:
        write_fasta(rejected_seqs, out_rejected)
        logger.info(f"Wrote {len(rejected_seqs)} rejected sequences to {out_rejected}")

    stats["final_count"] = len(final_seqs)
    stats["total_removed"] = initial_count - len(final_seqs)

    # Cleanup temp dir
    if not keep_temp:
        import shutil
        shutil.rmtree(work_dir, ignore_errors=True)

    logger.info(f"Done — {initial_count} → {len(final_seqs)} "
                f"({stats['total_removed']} removed total)")
    return Ok(stats)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rewrite_fasta(
    entries: list,
    original_sequences: Dict[str, str],
    fasta_path: Path,
) -> None:
    """Rewrite working FASTA with only surviving entries."""
    surviving = OrderedDict(
        (e.id, original_sequences[e.id])
        for e in entries
        if e.id in original_sequences
    )
    write_fasta(surviving, fasta_path)
