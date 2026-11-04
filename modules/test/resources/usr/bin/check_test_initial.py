#!/usr/bin/env python3
"""Check an initial clustering run on test_data/vjdb1_test_full_dataset.fna.gz.

Usage: check_test_initial.py <merged_reps.csv> <merged_reps.fna.gz> [expected.csv]

The synthetic dataset is built so that clustering it recovers the committed
test clustering. Leiden may pick other members as representatives, so the
check compares cluster membership, not representative IDs.
"""

import gzip
import sys

import pandas as pd

csv_file, fasta_file = sys.argv[1], sys.argv[2]
expected_file = sys.argv[3] if len(sys.argv) > 3 else "test_data/vjdb1_merged_reps.csv"

got = pd.read_csv(csv_file, dtype=str)
expected = pd.read_csv(expected_file, dtype=str)

assert got["vjdb_id"].is_unique, "duplicated vjdb_id in output CSV"
assert set(got["vjdb_id"]) == set(expected["vjdb_id"]), "output CSV has different sequences than expected"

got = got.set_index("vjdb_id").loc[expected["vjdb_id"]]
expected = expected.set_index("vjdb_id")
for col in ("copy_of", "seq_hash"):
    diff = (got[col] != expected[col]).sum()
    assert diff == 0, f"{diff} sequences have a different {col} than expected"


def partition(reps):
    return {frozenset(members.index) for _, members in reps.groupby(reps)}


assert partition(got["cluster_rep"]) == partition(expected["cluster_rep"]), "clusters differ from the expected clustering"

with gzip.open(fasta_file, "rt") as handle:
    fasta_ids = {line[1:].split()[0] for line in handle if line.startswith(">")}
reps = set(got["cluster_rep"])
assert fasta_ids == reps, (f"rep FASTA and CSV disagree: {len(reps - fasta_ids)} reps missing from FASTA, "
                           f"{len(fasta_ids - reps)} extra sequences in FASTA")

print(f"OK: all checks passed ({len(got)} sequences in {len(reps)} clusters, same clusters as expected)")
