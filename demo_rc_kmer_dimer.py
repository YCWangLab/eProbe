#!/usr/bin/env python3
"""
Demonstration of the corrected dimer scoring with RC k-mer index.

Shows the difference between:
1. Simple k-mer matching (wrong - biochemistry incorrect)
2. RC k-mer matching (correct - reflects Watson-Crick pairing)
"""

def _reverse_complement(seq: str) -> str:
    """Get reverse complement of sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def demo_rc_kmer_dimer():
    """Demonstrate the RC k-mer based dimer scoring."""
    
    print("=" * 90)
    print("DEMONSTRATION: RC K-mer Based Dimer Scoring (Biochemically Correct)")
    print("=" * 90)
    
    # Simple example: 3 probes of 20bp
    probes = {
        "probe_A": "ATCGATCGATCGATCGATCG",
        "probe_B": "ATCGATCGATCGATCGATCG",  # Identical to A
        "probe_C": "GCTAGCTAGCTAGCTAGCTA",  # Different
    }
    
    print("\n📊 INPUT PROBES:")
    for pid, seq in probes.items():
        print(f"  {pid}: 5'-{seq}-3'")
    
    # Show reverse complements
    print("\n🔄 REVERSE COMPLEMENTS (Binding partner sequences):")
    for pid, seq in probes.items():
        rc = _reverse_complement(seq)
        print(f"  {pid} RC:  3'-{rc}-5'")
    
    # Build RC k-mer index
    k = 11
    rc_kmer_freq = {}
    
    print(f"\n📚 BUILDING RC K-MER INDEX (k={k}):")
    print("\nStep 1: Extract RC k-mers from all probes")
    
    for pid, seq in probes.items():
        seq_upper = seq.upper()
        rc_seq = _reverse_complement(seq_upper)
        kmers = [rc_seq[i:i+k] for i in range(len(rc_seq) - k + 1) if 'N' not in rc_seq[i:i+k]]
        print(f"  {pid} RC → {len(kmers)} k-mers: {kmers[:2]}... (showing first 2)")
        
        for kmer in kmers:
            rc_kmer_freq[kmer] = rc_kmer_freq.get(kmer, 0) + 1
    
    print(f"\nTotal unique RC k-mers in pool: {len(rc_kmer_freq)}")
    print(f"Total RC k-mer instances: {sum(rc_kmer_freq.values())}")
    
    # Score each probe
    print(f"\n🧮 SCORING EACH PROBE (Query forward k-mers against RC index):")
    print(f"Formula: D(g) = (100/N) * Σ c_i(g)")
    print(f"         where c_i(g) = frequency of forward k-mer matches to RC pool")
    
    N = len(probes)
    scores = {}
    
    for pid, seq in probes.items():
        seq_upper = seq.upper()
        print(f"\n  {pid}:")
        print(f"    Forward: 5'-{seq}-3'")
        
        kmer_sum = 0
        matches = []
        
        for i in range(len(seq_upper) - k + 1):
            kmer = seq_upper[i:i+k]
            if 'N' not in kmer:
                freq = rc_kmer_freq.get(kmer, 0)
                if freq > 0:
                    matches.append((kmer, freq))
                    kmer_sum += freq
        
        print(f"    Matches found: {len(matches)}")
        if matches:
            print(f"      Example: {matches[0][0]} found {matches[0][1]}x in RC pool")
        
        score = (kmer_sum / N) * 100.0
        scores[pid] = round(score, 2)
        
        print(f"    Total matches: {kmer_sum}")
        print(f"    Score: ({kmer_sum}/{N}) * 100 = {score:.2f}")
    
    print("\n" + "=" * 90)
    print("📈 RESULTS:")
    print("=" * 90)
    
    for pid in sorted(scores.keys()):
        print(f"  {pid}: {scores[pid]}")
    
    print("\n💡 BIOLOGICAL INTERPRETATION:")
    print(f"""
  Probe A & B score identically ({scores['probe_A']}) because:
    - They are identical sequences
    - Their forward k-mers perfectly match each other's RC
    - High dimer risk! Should not use both in same multiplex
  
  Probe C scores differently ({scores['probe_C']}) because:
    - Sequence is unrelated to A and B
    - Forward k-mers don't match A/B's RC sequences
    - Lower dimer risk with A and B
  
  Why RC k-mers matter:
    - Dimer = probe1_forward + probe2_RC binding
    - We score this by: probe's forward k-mers vs pool's RC k-mers
    - If many matches → high dimer risk
    """)
    
    print("\n" + "=" * 90)
    print("COMPARISON: Why RC k-mer index is biochemically correct")
    print("=" * 90)
    
    print("""
WRONG (simple forward k-mer matching):
  - Looks at: k-mers in probe vs k-mers in other probes
  - Would score identical probes differently if we don't think about RC
  - Doesn't reflect actual Watson-Crick pairing
  
CORRECT (RC k-mer index):
  - Looks at: forward k-mers of probe vs RC k-mers of pool
  - Reflects actual molecular binding: 5'→3' binds to 3'→5'
  - A+B dimer = probe_A_forward + probe_B_complement
  - Score correctly reflects: how much can this probe pair with others?
    """)
    
    print("\n" + "=" * 90)


if __name__ == '__main__':
    demo_rc_kmer_dimer()
