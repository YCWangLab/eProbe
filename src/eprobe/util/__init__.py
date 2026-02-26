"""
Utility implementation modules.

Post-design probe manipulation toolkit:
  merge    - Merge probe sets with sequence-level dedup
  tile     - Simple FASTA sliding-window tiling
  adapter  - Add 5'/3' adapter sequences
  rename   - Batch rename with prefix or ID map
  assess   - Biophysical tags assessment & filtering
  sample   - Random subsampling & CD-HIT clustering
  target   - Bowtie2 mapping → target BED generation
  filter   - FASTA-level popgen filter stages (BG/AC/TX/biophysical)

Legacy tools (still available):
  tiling   - Multi-position SNP-based probe tiling
  dedup    - Exact/CD-HIT deduplication
  subset   - Filter/sample probes
  convert  - Format conversion
  validate - File validation
"""

# --- New modules ---
from eprobe.util.merge import run_merge
from eprobe.util.tile import run_tile
from eprobe.util.adapter import run_adapter
from eprobe.util.rename import run_rename
from eprobe.util.assess import run_assess as run_util_assess
from eprobe.util.sample import run_sample
from eprobe.util.target import run_target
from eprobe.util.filter import run_util_filter

# --- Legacy modules ---
from eprobe.util.tiling import run_tiling
from eprobe.util.dedup import run_dedup
from eprobe.util.subset import run_subset
from eprobe.util.convert import run_convert
from eprobe.util.validate import run_validate

__all__ = [
    # New
    "run_merge",
    "run_tile",
    "run_adapter",
    "run_rename",
    "run_util_assess",
    "run_sample",
    "run_target",
    "run_util_filter",
    # Legacy
    "run_tiling",
    "run_dedup",
    "run_subset",
    "run_convert",
    "run_validate",
]
