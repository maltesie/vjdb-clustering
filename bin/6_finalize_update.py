#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import os, gzip
import pandas as pd

parser = argparse.ArgumentParser(description="Finalize clustering update: remap clusters, update CSVs, extract reps.")
parser.add_argument("--outdir", required=True, help="Directory with intermediate files (vclust output, new_hashed_reps.csv)")
parser.add_argument("--existing-csv", required=True, help="Path to existing vjdb1_merged_reps.csv")
parser.add_argument("--output-prefix", required=True, help="Prefix for output files")
args = parser.parse_args()

# Read vclust leiden results from the combined run
leiden_df = pd.read_csv(os.path.join(args.outdir, "clusterreps_leiden.tsv"), sep="\t", dtype=str)
leiden_dict = {row.object: row.cluster for row in leiden_df.itertuples()}
cluster_set = set(leiden_df["cluster"])

print(f"Leiden clustering: {len(leiden_df)} objects in {len(cluster_set)} clusters")

# Read new sequence hash-dedup info
new_hashed_df = pd.read_csv(os.path.join(args.outdir, "new_hashed_reps.csv"), dtype={"seq_id": str, "copy_of": str})

# Read existing full CSV
existing_df = pd.read_csv(args.existing_csv, dtype={"vjdb_id": str, "copy_of": str, "cluster_rep": str})

# Remap existing sequences: old cluster_rep -> new leiden cluster
# Old cluster_rep IDs were in the combined vclust input, so they appear in leiden_dict
cant_remap = 0
updated_existing_clusters = []
for row in existing_df.itertuples():
    old_rep = str(row.cluster_rep)
    if old_rep in leiden_dict:
        updated_existing_clusters.append(leiden_dict[old_rep])
    else:
        updated_existing_clusters.append(old_rep)
        cant_remap += 1

existing_df["cluster_rep"] = updated_existing_clusters
if cant_remap:
    print(f"Warning: {cant_remap} existing sequences could not be remapped (kept old cluster_rep)")

# Map new sequences: copy_of -> leiden cluster
cant_map = 0
new_clusters = []
for row in new_hashed_df.itertuples():
    rep = str(row.copy_of)
    if rep in leiden_dict:
        new_clusters.append(leiden_dict[rep])
    else:
        new_clusters.append(rep)
        cant_map += 1

new_hashed_df["cluster_rep"] = new_clusters
if cant_map:
    print(f"Warning: {cant_map} new sequences could not be mapped to a cluster (using copy_of as fallback)")

# Write small CSV for just the new sequences
new_out = f"{args.output_prefix}_new_clusters.csv"
new_hashed_df.to_csv(new_out, index=False)
print(f"Wrote {len(new_hashed_df)} new sequence assignments to {new_out}")

# Append new sequences to existing CSV (rename seq_id -> vjdb_id to match schema)
new_for_merge = new_hashed_df.rename(columns={"seq_id": "vjdb_id"})
updated_df = pd.concat([existing_df, new_for_merge], ignore_index=True)
updated_csv = f"{args.output_prefix}_merged_reps.csv"
updated_df.to_csv(updated_csv, index=False)
print(f"Wrote {len(updated_df)} total sequences to {updated_csv}")

# Extract updated cluster representative sequences from the combined FASTA
combined_file = os.path.join(args.outdir, "combined_for_vclust.fna.gz")
out_fasta = f"{args.output_prefix}_merged_reps.fna.gz"
extracted = 0
with gzip.open(out_fasta, "wt", 2) as out_handle:
    with gzip.open(combined_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in cluster_set:
                out_handle.write(record.format("fasta"))
                extracted += 1

print(f"Extracted {extracted} cluster representative sequences to {out_fasta}")
