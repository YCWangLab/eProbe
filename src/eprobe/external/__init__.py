"""
External tools bundled with eProbe.

This module contains third-party tools that have been integrated
into eProbe for convenience.

- easySFS: Site Frequency Spectrum calculation tool
  Original: https://github.com/isaacovercast/easySFS
  License: See easySFS.py header
"""

from pathlib import Path

# Get the path to this directory for finding bundled scripts
EXTERNAL_DIR = Path(__file__).parent

def get_easysfs_path() -> Path:
    """Get path to bundled easySFS.py script."""
    return EXTERNAL_DIR / "easySFS.py"
