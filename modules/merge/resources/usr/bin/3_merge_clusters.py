#!/usr/bin/env python3

import argparse
import glob
import os
import pandas as pd
from collections import Counter

parser = argparse.ArgumentParser(description="Merge per-chunk clusters into final cluster assignments.")
parser.add_argument("--hashed-csv", default=None, help="Hash-dedup CSV (default: chunks/vjdb1_hashed_reps.csv)")
parser.add_argument("--merged-leiden", default=None, help="Merged-rep Leiden TSV (default: chunks/vjdb1_merged/clusterreps_leiden.tsv)")
parser.add_argument("--chunk-leidens", nargs="+", default=None, help="Per-chunk Leiden TSV files (default: glob chunks/)")
parser.add_argument("--chunk-linclusts", nargs="+", default=None, help="Per-chunk linclust cluster TSV files (default: glob chunks/)")
parser.add_argument("-o", "--output", default="vjdb1_merged_reps.csv", help="Output CSV")
args = parser.parse_args()

# Fall back to hardcoded paths if no explicit files given (backward compat)
if args.hashed_csv is None:
    args.hashed_csv = "chunks/vjdb1_hashed_reps.csv"
if args.merged_leiden is None:
    args.merged_leiden = "chunks/vjdb1_merged/clusterreps_leiden.tsv"
if args.chunk_leidens is None:
    args.chunk_leidens = sorted(glob.glob("chunks/vjdb1_[0-9]*/clusterreps_leiden.tsv"))
if args.chunk_linclusts is None:
    args.chunk_linclusts = sorted(glob.glob("chunks/vjdb1_[0-9]*/linclust_cluster.tsv"))

vjdb_hashed_df = pd.read_csv(args.hashed_csv, sep=",")
merged_df = pd.read_csv(args.merged_leiden, sep="\t")

chunk_leidens = sorted(args.chunk_leidens)
chunk_linclusts = sorted(args.chunk_linclusts)

if len(chunk_leidens) != len(chunk_linclusts):
    raise ValueError(f"Mismatch: {len(chunk_leidens)} leiden files vs {len(chunk_linclusts)} linclust files")

print(f"Collecting {len(chunk_leidens)} chunks")

linclust_count = 0
cluster_dicts = []
for leiden_file, linclust_file in zip(chunk_leidens, chunk_linclusts):
    linclust_df = pd.read_csv(linclust_file, sep="\t", names=["cluster", "object"])
    linclust_count += len(Counter(linclust_df["cluster"]))
    vclust_df = pd.read_csv(leiden_file, sep="\t")
    vclust_dict = {row.object: row.cluster for row in vclust_df.itertuples()}
    cluster_dicts.append({row.object: vclust_dict[row.cluster] for row in linclust_df.itertuples()})

merged_dict = {row.object: row.cluster for row in merged_df.itertuples()}
vclust_dict = {}
for d in cluster_dicts:
    vclust_dict |= d

print("Mapping to final clusters")
final_clusters = [""] * len(vjdb_hashed_df)
cant_map = []
cant_map_merged = []
for i, row in enumerate(vjdb_hashed_df.itertuples()):
    if row.copy_of in vclust_dict:
        a = vclust_dict[row.copy_of]
        if a in merged_dict:
            final_clusters[i] = merged_dict[a]
        else:
            final_clusters[i] = a
            cant_map_merged.append(a)
    else:
        cant_map.append(row.copy_of)
        final_clusters[i] = row.copy_of

print("cant map merged:", Counter(cant_map_merged).most_common(10))
print("cant map:", Counter(cant_map).most_common(10))

all_count = len(vjdb_hashed_df)
hash_count = len(Counter(vjdb_hashed_df["copy_of"]))
vclust_count = len(merged_df)
merged_count = len(Counter(merged_df["cluster"]))

print(f"all: {all_count}, hashed: {hash_count}, linclust: {linclust_count}, vclust: {vclust_count}, merged: {merged_count}")

vjdb_hashed_df["cluster_rep"] = final_clusters
vjdb_hashed_df.to_csv(args.output, index=False)
