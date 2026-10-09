#!/usr/bin/env python3
"""
Generate test inputs for the update step (run from the repo root).

vjdb1_test_new_sequences.fna.gz covers every case 5_prepare_update.py and
6_finalize_update.py handle:
  - 20 random sequences                          -> new clusters
  - 4 exact copies of random ones (new IDs)      -> new_copy
  - 1 random sequence repeated with the same ID  -> duplicate_in_batch
  - 10 exact copies of existing reps (new IDs)   -> existing_copy
  - 1 exact copy of an existing non-rep member   -> existing_copy
  - 1 existing member with its own ID            -> resubmitted
  - 3 reps with ~1% substitutions                -> grow existing clusters
  - 1 chimera of two singleton reps              -> merges two existing clusters

vjdb1_test_clash.fna.gz holds one existing ID with a different sequence; the
update must refuse it.
"""

import gzip
import random

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

random.seed(42)

reps_file = "vjdb1_merged_reps.fna.gz"
full_file = "vjdb1_test_full_dataset.fna.gz"
csv_file = "vjdb1_merged_reps.csv"
out_file = "vjdb1_test_new_sequences.fna.gz"
clash_file = "vjdb1_test_clash.fna.gz"


def rec(seq, seq_id):
    return SeqRecord(Seq(seq), id=seq_id, description="")


def mutate(seq, rate=0.01):
    return "".join(random.choice([b for b in "ACGT" if b != c]) if random.random() < rate else c for c in seq)


with gzip.open(reps_file, "rt") as handle:
    reps = [(r.id, str(r.seq)) for r in SeqIO.parse(handle, "fasta")]
with gzip.open(full_file, "rt") as handle:
    full = {r.id: str(r.seq) for r in SeqIO.parse(handle, "fasta")}
df = pd.read_csv(csv_file, dtype=str)
size = df["cluster_rep"].value_counts()
rep_ids = set(size.index)
members = df[~df["vjdb_id"].isin(rep_ids) & (df["vjdb_id"] == df["copy_of"])]["vjdb_id"].tolist()

records = []

# New random sequences, some copied under new IDs, one repeated with the same ID
randoms = [rec("".join(random.choices("ACGT", k=random.randint(500, 50000))), f"random_seq_{i}") for i in range(20)]
records += randoms
for i in [3, 7, 7, 12]:
    records.append(rec(str(randoms[i].seq), f"dup_of_{randoms[i].id}_{len(records)}"))
records.append(rec(str(randoms[5].seq), randoms[5].id))

# Exact copies of existing data under new IDs
for rep_id, seq in reps[:10]:
    records.append(rec(seq, f"existing_{rep_id}"))
records.append(rec(full[members[0]], f"copy_of_{members[0]}"))

# Resubmission of an existing member
records.append(rec(full[members[1]], members[1]))

# Near-copies of reps of larger clusters -> grown
big = [r for r, n in size.items() if n > 1][:3]
seqs = dict(reps)
for rep_id in big:
    records.append(rec(mutate(seqs[rep_id]), f"variant_of_{rep_id}"))

# Chimera of two mid-length singleton reps -> merges their clusters
singles = sorted((r for r, n in size.items() if n == 1 and 2000 < len(seqs[r]) < 5000), key=lambda r: len(seqs[r]))
a, b = singles[0], singles[1]
records.append(rec(seqs[a] + seqs[b], f"chimera_{a}_{b}"))

random.shuffle(records)
with gzip.open(out_file, "wt", 2) as handle:
    SeqIO.write(records, handle, "fasta")

with gzip.open(clash_file, "wt", 2) as handle:
    SeqIO.write([rec("".join(random.choices("ACGT", k=1000)), members[2])], handle, "fasta")

print(f"Wrote {len(records)} sequences to {out_file} (expected merge: {a} + {b})")
print(f"Wrote 1 clashing sequence ({members[2]}) to {clash_file}")
