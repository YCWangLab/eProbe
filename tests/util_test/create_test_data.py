"""Create test data for util module testing."""
import os
import random

d = os.path.dirname(os.path.abspath(__file__))

# set_a.fa: 5 probes, 81bp each
seqs_a = {
    "probe_A1": "ATCGATCG" * 10 + "A",
    "probe_A2": "GCTAGCTA" * 10 + "G",
    "probe_A3": "TTTTAAAA" * 10 + "T",
    "probe_A4": "AACCGGTT" * 10 + "A",
    "probe_A5": "GCGCATAT" * 10 + "G",
}
with open(os.path.join(d, "set_a.fa"), "w") as f:
    for k, v in seqs_a.items():
        f.write(f">{k}\n{v}\n")

# set_b.fa: 4 probes, probe_B1 = dup of A1, probe_B4 = dup of A2
seqs_b = {
    "probe_B1": "ATCGATCG" * 10 + "A",
    "probe_B2": "CCCAATTG" * 10 + "C",
    "probe_B3": "TATAGATC" * 10 + "T",
    "probe_B4": "GCTAGCTA" * 10 + "G",
}
with open(os.path.join(d, "set_b.fa"), "w") as f:
    for k, v in seqs_b.items():
        f.write(f">{k}\n{v}\n")

# gene.fa: 2 longer sequences for tiling
with open(os.path.join(d, "gene.fa"), "w") as f:
    f.write(">gene1\n" + "ATCGATCG" * 34 + "AT\n")
    f.write(">gene2\n" + "GCTAGCTA" * 33 + "GCTA\n")

# probes50.fa: 50 random probes for assess/sample testing
random.seed(42)
with open(os.path.join(d, "probes50.fa"), "w") as f:
    for i in range(50):
        seq = "".join(random.choices("ATCG", k=81))
        f.write(f">probe_{i + 1:04d}\n{seq}\n")

print("Test data created:")
for fn in sorted(os.listdir(d)):
    if fn.endswith(".fa"):
        path = os.path.join(d, fn)
        with open(path) as fh:
            lines = fh.readlines()
        n_seqs = sum(1 for l in lines if l.startswith(">"))
        print(f"  {fn}: {n_seqs} sequences")
