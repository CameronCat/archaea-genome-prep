"""
generate_example_data.py
------------------------
Generate a synthetic archaeal genome FASTA for testing.
The sequence is random but with GC content typical of Methanosarcina (~42%).

Run with:
    python generate_example_data.py
"""

import random
import os

random.seed(42)
OUTPUT_DIR = "data"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Methanosarcina-like GC content
GC_FRACTION = 0.42
AT_FRACTION = 1 - GC_FRACTION

BASES_GC = ["G", "C"]
BASES_AT = ["A", "T"]


def random_seq(length: int) -> str:
    bases = []
    for _ in range(length):
        if random.random() < GC_FRACTION:
            bases.append(random.choice(BASES_GC))
        else:
            bases.append(random.choice(BASES_AT))
    return "".join(bases)


# Generate a synthetic genome with 3 contigs
# (mimics a draft assembly)
contigs = [
    ("contig_1", random_seq(2_500_000)),  # main chromosome fragment
    ("contig_2", random_seq(1_800_000)),  # second large contig
    ("contig_3", random_seq(450_000)),    # smaller contig
    ("contig_4", random_seq(120)),        # short contig - will be filtered
]

output_path = f"{OUTPUT_DIR}/example_genome.fasta"
with open(output_path, "w") as f:
    for seq_id, seq in contigs:
        f.write(f">{seq_id} synthetic Methanosarcina-like genome\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

total = sum(len(s) for _, s in contigs)
print(f"Generated {output_path}")
print(f"  Contigs: {len(contigs)}")
print(f"  Total length: {total:,} bp")
print(f"  (contig_4 is intentionally short and will be filtered)")
print(f"\nRun: python -c \"from archaea_genome_prep import prepare; r = prepare('data/example_genome.fasta'); print(r.summary())\"")
