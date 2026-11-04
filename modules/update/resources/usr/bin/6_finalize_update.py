#!/usr/bin/env python3
"""
Finalize a clustering update.

Inputs come from 5_prepare_update.py and the vclust run on the combined FASTA.

Cluster IDs stay stable: a cluster that contains exactly one old representative
keeps it as its ID, even if Leiden picked another sequence. When several old
clusters merge, the merged cluster keeps the old representative of the largest
old cluster. Clusters made only of new sequences use the Leiden representative.

Outputs (<prefix>_...):
  merged_reps.csv     updated membership table (vjdb_id is unique)
  new_clusters.csv    every record of the new batch with status and cluster_rep
  cluster_events.tsv  one row per changed cluster: merged, grown or new
  merged_reps.fna.gz  representative sequences of all clusters

Splits can't happen: only representatives are re-clustered, and all members of
an old cluster follow their representative into exactly one new cluster.
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Finalize clustering update: remap clusters, report events, extract reps.")
parser.add_argument("--leiden-tsv", required=True, help="vclust Leiden TSV from the combined run")
parser.add_argument("--new-hashed-csv", required=True, help="new_hashed_reps.csv from 5_prepare_update.py")
parser.add_argument("--existing-csv", required=True, help="Existing membership CSV (e.g. vjdb1_merged_reps.csv)")
parser.add_argument("--combined-fasta", required=True, help="combined_for_vclust.fna.gz from 5_prepare_update.py")
parser.add_argument("--output-prefix", required=True, help="Prefix for output files")
args = parser.parse_args()

out_csv = f"{args.output_prefix}_merged_reps.csv"
out_new = f"{args.output_prefix}_new_clusters.csv"
out_events = f"{args.output_prefix}_cluster_events.tsv"
out_fasta = f"{args.output_prefix}_merged_reps.fna.gz"

inputs = {os.path.realpath(p) for p in (args.leiden_tsv, args.new_hashed_csv, args.existing_csv, args.combined_fasta)}
if any(os.path.realpath(p) in inputs for p in (out_csv, out_new, out_events, out_fasta)):
    sys.exit("Error: an output file would overwrite an input file; choose a different --output-prefix")

leiden_df = pd.read_csv(args.leiden_tsv, sep="\t", dtype=str)
new_df = pd.read_csv(args.new_hashed_csv, dtype=str)
existing_df = pd.read_csv(args.existing_csv, dtype=str)

with gzip.open(args.combined_fasta, "rt") as handle:
    combined_ids = [r.id for r in SeqIO.parse(handle, "fasta")]

# --- Every sequence given to vclust must come back in the Leiden output --------
leiden = dict(zip(leiden_df["object"], leiden_df["cluster"]))
unmapped = [i for i in combined_ids if i not in leiden]
if unmapped:
    sys.exit(f"Error: {len(unmapped)} sequence(s) given to vclust are missing from {args.leiden_tsv}, "
             f"e.g. {unmapped[:10]}. Check for ID mismatches or sequences vclust skipped.")

# --- Stable cluster IDs ------------------------------------------------------
old_size = existing_df["cluster_rep"].value_counts()
members = defaultdict(list)
for obj, cl in leiden.items():
    members[cl].append(obj)

final_id = {}  # sequence in vclust input -> final cluster ID
old_reps_of = {}  # final cluster ID -> old reps it contains
for cl, objs in members.items():
    old_reps = sorted((o for o in objs if o in old_size.index), key=lambda o: (-old_size[o], o))
    cid = old_reps[0] if old_reps else cl
    old_reps_of[cid] = old_reps
    for o in objs:
        final_id[o] = cid

# --- Existing sequences follow their old representative -------------------------
existing_df["cluster_rep"] = existing_df["cluster_rep"].map(final_id)
cluster_of = dict(zip(existing_df["vjdb_id"], existing_df["cluster_rep"]))

# --- New sequences -----------------------------------------------------------
# Rows are in batch order, so a copy always comes after the sequence it copies
added_statuses = ("new_unique", "new_copy", "existing_copy")
new_clusters = []
for r in new_df.itertuples():
    if r.status == "new_unique":
        c = final_id[r.vjdb_id]
    elif r.status in ("new_copy", "existing_copy"):
        c = cluster_of[r.copy_of]
    else:  # resubmitted, duplicate_in_batch: same ID seen before
        c = cluster_of[r.vjdb_id]
    if r.status in added_statuses:
        cluster_of[r.vjdb_id] = c
    new_clusters.append(c)
new_df["cluster_rep"] = new_clusters
new_df.to_csv(out_new, index=False)
print(f"Wrote {len(new_df)} new batch assignments to {out_new}")

added = new_df[new_df["status"].isin(added_statuses)].drop(columns="status")
updated_df = pd.concat([existing_df, added], ignore_index=True)
if not updated_df["vjdb_id"].is_unique:
    dups = updated_df.loc[updated_df["vjdb_id"].duplicated(), "vjdb_id"].unique()
    sys.exit(f"Error: duplicated vjdb_id(s) in updated table, e.g. {list(dups[:10])}")
updated_df.to_csv(out_csv, index=False)
print(f"Wrote {len(updated_df)} total sequences to {out_csv}")

# --- Cluster events ----------------------------------------------------------
n_existing = existing_df["cluster_rep"].value_counts()
n_new = added["cluster_rep"].value_counts()
events = []
for cid, old_reps in old_reps_of.items():
    if len(old_reps) >= 2:
        event = "merged"
    elif not old_reps:
        event = "new"
    elif n_new.get(cid, 0) > 0:
        event = "grown"
    else:
        event = "unchanged"
    events.append({"cluster_rep": cid, "event": event, "old_cluster_reps": ";".join(old_reps),
                   "n_existing": int(n_existing.get(cid, 0)), "n_new": int(n_new.get(cid, 0))})
events_df = pd.DataFrame(events, columns=["cluster_rep", "event", "old_cluster_reps", "n_existing", "n_new"])
event_order = {"merged": 0, "grown": 1, "new": 2, "unchanged": 3}
events_df = events_df.sort_values(["event", "cluster_rep"], key=lambda c: c.map(event_order) if c.name == "event" else c)
events_df[events_df["event"] != "unchanged"].to_csv(out_events, sep="\t", index=False)
counts = events_df["event"].value_counts()
print("Cluster events: " + ", ".join(f"{e}={counts.get(e, 0)}" for e in ("merged", "grown", "new", "unchanged")))
if counts.get("merged", 0):
    n_absorbed = sum(len(r) - 1 for r in old_reps_of.values() if len(r) > 1)
    print(f"  {n_absorbed} old cluster(s) were merged into others; see {out_events}")

# --- Representative FASTA ----------------------------------------------------
final_reps = set(old_reps_of)
extracted = 0
with gzip.open(out_fasta, "wt", 2) as out_handle, gzip.open(args.combined_fasta, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.id in final_reps:
            out_handle.write(record.format("fasta"))
            extracted += 1
if extracted != len(final_reps):
    sys.exit(f"Error: extracted {extracted} representatives but expected {len(final_reps)}")
print(f"Extracted {extracted} cluster representative sequences to {out_fasta}")
