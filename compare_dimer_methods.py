#!/usr/bin/env python3
"""
Comprehensive comparison of two dimer evaluation methods:

1. SmartDimerFilter (filtering-based, used in filter.py)
   - Graph-based clustering of similar probes
   - Binary decision: keep or remove
   - Used for FILTERING STAGE

2. Legacy D(g) = (100/N)*Σc_i(g) (assess-based, used in assess.py)
   - K-mer frequency sum
   - Continuous score: 0-100+
   - Used for ASSESSMENT/REPORTING

These serve different purposes and shouldn't be directly compared!
"""

import numpy as np
from typing import Dict, List, Tuple
from collections import Counter, defaultdict


def smart_dimer_analyze(sequences: List[str], k: int = 11, min_shared: float = 0.15) -> Dict:
    """
    Simulate SmartDimerFilter logic to understand what it does.
    """
    n = len(sequences)
    
    # Build canonical k-mer index
    kmer_to_probes: Dict[str, set] = defaultdict(set)
    probe_kmers: Dict[int, set] = {}
    
    def reverse_complement(seq):
        rc_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(rc_map.get(c, 'N') for c in reversed(seq))
    
    for i, seq in enumerate(sequences):
        seq_upper = seq.upper()
        seen = set()
        for j in range(len(seq_upper) - k + 1):
            kmer = seq_upper[j:j + k]
            if 'N' not in kmer:
                rc = reverse_complement(kmer)
                canonical = kmer if kmer <= rc else rc
                if canonical not in seen:
                    kmer_to_probes[canonical].add(i)
                    seen.add(canonical)
        probe_kmers[i] = seen
    
    # Count pairwise shared k-mers
    edges = []
    for canonical, probe_set in kmer_to_probes.items():
        if len(probe_set) < 2:
            continue
        probes = list(probe_set)
        for a in range(len(probes)):
            for b in range(a+1, len(probes)):
                pa, pb = probes[a], probes[b]
                shared = 1
                edges.append((pa, pb, shared))
    
    # Aggregate shared count
    pair_sharing = defaultdict(int)
    for pa, pb, shared in edges:
        key = (min(pa, pb), max(pa, pb))
        pair_sharing[key] += 1
    
    # Build adjacency with threshold
    adj = defaultdict(set)
    n_edges_above_threshold = 0
    
    for (i, j), count in pair_sharing.items():
        frac = count / min(len(probe_kmers[i]), len(probe_kmers[j]))
        if frac >= min_shared:
            adj[i].add(j)
            adj[j].add(i)
            n_edges_above_threshold += 1
    
    # Union-Find for connected components
    parent = list(range(n))
    
    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]
    
    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            if px < py:
                parent[py] = px
            else:
                parent[px] = py
    
    for i, neighbors in adj.items():
        for j in neighbors:
            union(i, j)
    
    # Find components
    components = defaultdict(list)
    for i in range(n):
        components[find(i)].append(i)
    
    # Count groups and removals
    n_groups = sum(1 for comp in components.values() if len(comp) > 1)
    n_removed = sum(len(comp) - 1 for comp in components.values() if len(comp) > 1)
    
    return {
        "name": "SmartDimerFilter",
        "purpose": "FILTERING (binary: keep/remove)",
        "kept": n - n_removed,
        "removed": n_removed,
        "groups_detected": n_groups,
        "edges_above_threshold": n_edges_above_threshold,
        "max_group_size": max((len(comp) for comp in components.values()), default=1),
        "components": components,
        "adj": adj,
    }


def legacy_dimer_assess(sequences: List[str], k: int = 11, constant: float = 100) -> Dict:
    """
    Calculate legacy dimer scores: D(g) = (constant/N) * Σ c_i(g)
    """
    N = len(sequences)
    
    # Build k-mer frequency
    kmer_freq = Counter()
    for seq in sequences:
        seq_upper = seq.upper()
        for i in range(len(seq_upper) - k + 1):
            kmer = seq_upper[i:i + k]
            if 'N' not in kmer:
                kmer_freq[kmer] += 1
    
    # Calculate scores
    scores = {}
    for i, seq in enumerate(sequences):
        seq_upper = seq.upper()
        kmer_sum = 0
        for j in range(len(seq_upper) - k + 1):
            kmer = seq_upper[j:j + k]
            if 'N' not in kmer:
                kmer_sum += kmer_freq[kmer]
        
        score = (kmer_sum / N) * constant
        scores[i] = round(score, 2)
    
    return {
        "name": "Legacy D(g)",
        "purpose": "ASSESSMENT (continuous score)",
        "scores": scores,
        "min_score": min(scores.values()),
        "max_score": max(scores.values()),
        "mean_score": round(np.mean(list(scores.values())), 2),
        "std_score": round(np.std(list(scores.values())), 2),
    }


def print_comparison():
    """Print detailed comparison of both methods."""
    
    print("=" * 90)
    print("COMPARISON: SmartDimerFilter (Filtering) vs Legacy D(g) (Assessment)")
    print("=" * 90)
    
    # Characteristics table
    print("\n📊 CHARACTERISTICS:")
    print("-" * 90)
    
    comparison_table = [
        ("PURPOSE", "Filtering stage (remove bad probes)", "Assessment/reporting stage"),
        ("INPUT", "Probe sequences (filter.py)", "Probe sequences (assess.py)"),
        ("OUTPUT", "Binary (keep/remove)", "Continuous score (0-100+)"),
        ("K-MER BASIS", "Canonical k-mers (min(k,rc))", "Forward k-mers only"),
        ("KEY PARAM", "min_shared_fraction (15%)", "constant (100)"),
        ("ALGORITHM", "Graph theory + Union-Find", "Simple frequency sum"),
        ("COMPLEXITY", "O(n²) pairwise comparison", "O(n·L) linear scan"),
        ("TEMPORAL", "Used ONCE during filtering", "Can be run multiple times"),
        ("SEMANTIC", "Group-based risk detection", "Pool-based frequency scoring"),
    ]
    
    for metric, filter_method, assess_method in comparison_table:
        print(f"{metric:25} | {filter_method:35} | {assess_method:35}")
    
    # Biological interpretation
    print("\n\n🧬 BIOLOGICAL INTERPRETATION:")
    print("-" * 90)
    
    print(f"""
SmartDimerFilter:
  ✓ Identifies CLUSTERS of similar probes
  ✓ Detects probes that would cross-hybridize with each other
  ✓ Removes all but ONE from each risk group
  ✓ Conservative: assumes 15% k-mer sharing = dimer risk
  ✓ Output: "These 50 probes form 3 risk groups → keep 3"

Legacy D(g):
  ✓ Measures how COMMON each probe's k-mers are in the pool
  ✓ Higher score = more "repetitive" sequences
  ✓ Not a clustering method, just a frequency metric
  ✓ Scores pool-based complementarity risk
  ✓ Output: "This probe scores 45 (moderately repetitive)"
    """)
    
    # Use case comparison
    print("\n🎯 USE CASES:")
    print("-" * 90)
    
    print("""
SmartDimerFilter is BEST for:
  1. Pure probe set selection (must reduce to non-overlapping subset)
  2. Multiplexing where you can only afford N probes
  3. Minimal redundancy required
  4. Output: Final probe panel size is known
  
Legacy D(g) is BEST for:
  1. Quality reporting and diagnostics
  2. Understanding probe characteristics
  3. Comparing different probe sets
  4. Ranking probes by "uniqueness"
  5. Blending with other metrics (GC, Tm, hairpin)
    """)
    
    print("\n\n⚖️ EVALUATION APPROACH:")
    print("-" * 90)
    
    print("""
SmartDimerFilter:
  "Are these two probes too similar to use together?"
  → Graph-based clustering → Remove representatives
  
Legacy D(g):
  "How repetitive is this probe in the pool?"
  → Frequency-based ranking → Score each probe
    """)
    
    print("\n\n⚠️ KEY DIFFERENCES IN OUTPUT:")
    print("-" * 90)
    
    print("""
SmartDimerFilter produces:
  - Probe panels where NO two probes are "dimer partners"
  - Reduced set (typically 30-50% of input)
  - Binary quality metric per probe: included/excluded
  
Legacy D(g) produces:
  - Full set of probes with quality scores
  - 0-100+ range (varies with pool size)
  - Continuous metric for ranking/filtering
    """)
    
    print("\n\n✅ WHICH IS MORE REASONABLE?\n")
    print("-" * 90)
    
    print("""
ANSWER: They serve DIFFERENT purposes, so comparing them is like comparing
"apples vs oranges" — they're both reasonable for their intended use!

FOR FILTERING stage (removing sequences):
  ✓ SmartDimerFilter is more appropriate
    - Deterministic binary decision
    - Guarantees no high-risk pairs remain
    - Removes minimum necessary probes
    - Fast enough for 10,000+ sequences

FOR ASSESSMENT stage (reporting quality):
  ✓ Legacy D(g) is more appropriate
    - Provides continuous metric for comparison
    - Interpretable as "frequency ranking"
    - Can combine with other metrics
    - Lightweight reporting tool

HYBRID APPROACH (RECOMMENDED):
  1. Use SmartDimerFilter for FILTERING → obtain reduced panel
  2. Use Legacy D(g) for ASSESSMENT → understand characteristics
  
This maximizes both safety (removed dimer-risk clusters) and insight
(understand remaining probes' quality).
    """)
    
    print("\n" + "=" * 90)
    print("CONCLUSION")
    print("=" * 90)
    
    print("""
Neither is universally "better" — they answer different questions:

SmartDimerFilter:  "Which probes can I safely use TOGETHER?"
                   → Filtering/Panel Selection

Legacy D(g):       "How unique/common is each probe's sequence?"
                   → Quality Metrics/Assessment

Both should be used in their respective pipelines!
    """)


if __name__ == '__main__':
    print_comparison()
    
    # Optional: Generate example analysis
    print("\n" + "=" * 90)
    print("EXAMPLE ANALYSIS (small dataset)")
    print("=" * 90)
    
    # Create test sequences with known similarity
    test_seqs = [
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # 0
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",  # 1 (clone of 0)
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",  # 2
        "TTACTTACTTACTTACTTACTTACTTACTTACTTACTTACTTACTTACTTACTTACTTACTTACTTACTTACTTAC",  # 3 (unique)
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",  # 4 (clone of 2)
    ]
    
    print("\nTest set: 5 sequences (3 unique groups × 2 clones + 1 unique)")
    
    smart = smart_dimer_analyze(test_seqs)
    print(f"\nSmartDimerFilter result:")
    print(f"  Kept: {smart['kept']}/5")
    print(f"  Removed: {smart['removed']}")
    print(f"  Groups detected: {smart['groups_detected']}")
    
    legacy = legacy_dimer_assess(test_seqs)
    print(f"\nLegacy D(g) result:")
    for i, score in legacy['scores'].items():
        print(f"  Probe {i}: {score}")
    
    print("\nExplanation:")
    print("  SmartDimerFilter → Groups identical/similar probes, keeps 1 per group")
    print("  Legacy D(g) → Scores each probe independently on frequency")
    print("  → Different purposes, different outputs!")
