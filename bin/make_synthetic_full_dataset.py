#!/usr/bin/env python3
"""
Build a synthetic "full" dataset for testing, consistent with the committed
test outputs of Phase 1 in test_data/ (vjdb1_merged_reps.csv / vjdb1_merged_reps.fna.gz).

Only the cluster representatives have real sequences in the repo. For every
other member listed in the CSV this script creates a sequence:
  - exact copies (copy_of != vjdb_id) get the sequence of their copy_of entry
  - all other members get their cluster rep's sequence with ~1% substitutions,
    trimmed or padded to the seq_len recorded in the CSV

It then writes vjdb1_test_full_dataset.fna.gz and adds/refreshes the
seq_hash (md5) column in the CSV. Test data only, run once from the repo root.
"""

import gzip
import hashlib
import random

import pandas as pd
from Bio import SeqIO

random.seed(1)

reps_file = "test_data/vjdb1_merged_reps.fna.gz"
csv_file = "test_data/vjdb1_merged_reps.csv"
out_fasta = "test_data/vjdb1_test_full_dataset.fna.gz"
MUTATION_RATE = 0.01

with gzip.open(reps_file, "rt") as handle:
    rep_seqs = {r.id: str(r.seq) for r in SeqIO.parse(handle, "fasta")}

df = pd.read_csv(csv_file, dtype=str)
df = df.drop(columns=["seq_hash"], errors="ignore")


def mutate(seq, length):
    s = list(seq[:length])
    while len(s) < length:
        s.append(random.choice("ACGT"))
    for i in range(len(s)):
        if random.random() < MUTATION_RATE:
            s[i] = random.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


seqs = {}
# Unique sequences first (reps and non-copy members), then the exact copies
for row in df.itertuples():
    if row.copy_of != row.vjdb_id:
        continue
    if row.vjdb_id in rep_seqs:
        seqs[row.vjdb_id] = rep_seqs[row.vjdb_id]
    else:
        seqs[row.vjdb_id] = mutate(rep_seqs[row.cluster_rep], int(row.seq_len))
for row in df.itertuples():
    if row.copy_of != row.vjdb_id:
        seqs[row.vjdb_id] = seqs[row.copy_of]

hashes = [hashlib.md5(seqs[i].encode()).hexdigest() for i in df["vjdb_id"]]
cols = list(df.columns)
df.insert(cols.index("cluster_rep"), "seq_hash", hashes)
df.to_csv(csv_file, index=False)

with gzip.open(out_fasta, "wt", 2) as handle:
    for vjdb_id in df["vjdb_id"]:
        handle.write(f">{vjdb_id}\n{seqs[vjdb_id]}\n")

print(f"Wrote {len(df)} sequences to {out_fasta} and added seq_hash to {csv_file}")
