#!/usr/bin/env python3
"""Check the outputs of an update run on the data from make_test_data.py.

Usage: check_test_update.py <output-prefix> [existing-csv]
"""

import sys

import pandas as pd

prefix = sys.argv[1]
existing_csv = sys.argv[2] if len(sys.argv) > 2 else "vjdb1_merged_reps.csv"

existing = pd.read_csv(existing_csv, dtype=str)
merged = pd.read_csv(f"{prefix}_merged_reps.csv", dtype=str)
new = pd.read_csv(f"{prefix}_new_clusters.csv", dtype=str)
events = pd.read_csv(f"{prefix}_cluster_events.tsv", sep="\t", dtype=str)

expected = {"new_unique": 24, "new_copy": 4, "existing_copy": 11, "resubmitted": 1, "duplicate_in_batch": 1}
got = new["status"].value_counts().to_dict()
assert got == expected, f"status counts {got} != {expected}"

assert merged["vjdb_id"].is_unique, "duplicated vjdb_id in updated CSV"
assert len(merged) == len(existing) + 39, f"expected {len(existing) + 39} rows, got {len(merged)}"

merges = events[events["event"] == "merged"]
chimera = new.loc[new["vjdb_id"].str.startswith("chimera_"), "vjdb_id"].iat[0]
pair = set(chimera.split("_")[1:])
assert len(merges) == 1 and set(merges["old_cluster_reps"].iat[0].split(";")) == pair, \
    f"expected one merge of {pair}, got {merges.to_dict('records')}"

assert (events["event"] == "new").sum() == 20, "expected 20 new clusters"

variants = new[new["vjdb_id"].str.startswith("variant_of_")]
assert (variants["cluster_rep"] == variants["vjdb_id"].str.replace("variant_of_", "")).all(), \
    "variants should keep their rep's cluster ID"

# Stable IDs: apart from the merge, existing sequences keep their cluster ID
absorbed = set(merges["old_cluster_reps"].iat[0].split(";")) - {merges["cluster_rep"].iat[0]}
old = existing.set_index("vjdb_id")["cluster_rep"]
now = merged.set_index("vjdb_id").loc[old.index, "cluster_rep"]
changed = old[old != now]
assert set(changed) <= absorbed, f"unexpected cluster ID changes: {changed.head().to_dict()}"

print(f"OK: all checks passed ({len(merged)} sequences, events: "
      + ", ".join(f"{k}={v}" for k, v in events["event"].value_counts().items()) + ")")
