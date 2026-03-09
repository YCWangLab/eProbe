#!/usr/bin/env python3
"""Quick test for new features."""
import sys
sys.path.insert(0, "/Users/zh384/Desktop/scripts_dev/vs_code/developed_package/eProbe/src")

from eprobe.popgen.filter import BiophysicalThresholds
from eprobe.biophysics.biophysics import SmartDimerFilter

# Test 1: New defaults
t = BiophysicalThresholds()
print("=== New Defaults ===")
print(f"  GC:         {t.gc_min}-{t.gc_max}%")
print(f"  Tm:         {t.tm_min}-{t.tm_max}C")
print(f"  Complexity: {t.complexity_max}")
print(f"  Hairpin:    {t.hairpin}")
print(f"  Dimer:      {t.dimer} (smart: k-mer sharing fraction)")
assert t.dimer == 0.15, f"Expected dimer=0.15, got {t.dimer}"

# Test 2: SmartDimerFilter
sdf = SmartDimerFilter(k=5, min_shared_fraction=0.3)
seqs = [
    "ACGTACGTACGTACGT",       # probe 0
    "ACGTACGTACGTACGT",       # probe 1 (identical -> group)
    "TTTTTGGGGGTTTTTG",       # probe 2 (unrelated)
    "ACGTACGTACGAAAAA",       # probe 3 (partially similar)
]
keep_idx, stats = sdf.filter(seqs)
print(f"\n=== SmartDimerFilter Test (k=5, shared>=30%) ===")
print(f"  Input: {len(seqs)} probes")
print(f"  Kept: {keep_idx}")
print(f"  Groups: {stats['dimer_groups']}, Removed: {stats['dimer_failed']}")
print(f"  Edges: {stats['dimer_edges']}, Mode: {stats['mode']}")
assert stats["dimer_groups"] >= 1, "Should detect at least 1 dimer group"
assert len(keep_idx) < len(seqs), "Should remove some probes"

# Test 3: -1 disable
t2 = BiophysicalThresholds(gc_min=-1, gc_max=-1, tm_min=-1, tm_max=-1,
                            complexity_max=-1, hairpin=-1, dimer=-1)
assert t2.gc_min < 0
assert t2.hairpin <= 0
assert t2.dimer <= 0
print(f"\n=== Disable All (-1) ===")
print(f"  GC skip: {t2.gc_min < 0}")
print(f"  Hairpin skip: {t2.hairpin <= 0}")
print(f"  Dimer skip: {t2.dimer <= 0}")

print("\n=== ALL TESTS PASSED ===")
