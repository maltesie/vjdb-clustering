#!/usr/bin/env python3

import argparse
import hashlib
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os, gzip
import pandas as pd

parser = argparse.ArgumentParser(description="Prepare new sequences for clustering with existing representatives.")
parser.add_argument("--new-seqs", required=True, help="Path to new sequences FASTA.gz")
parser.add_argument("--existing-reps", required=True, help="Path to existing merged rep FASTA.gz (e.g. vjdb1_merged_reps.fna.gz)")
parser.add_argument("--outdir", required=True, help="Output directory for intermediate files")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

hash_to_id = dict()
seq_ids = []
rep_ids = []
lens = []
gcs = []

combined_file = os.path.join(args.outdir, "combined_for_vclust.fna.gz")
new_hashed_file = os.path.join(args.outdir, "new_hashed_reps.csv")

existing_count = 0
with gzip.open(combined_file, "wt", 2) as out_handle:
    # Write all existing cluster reps first
    with gzip.open(args.existing_reps, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            out_handle.write(record.format("fasta"))
            existing_count += 1

    print(f"Wrote {existing_count} existing cluster representatives")

    # Hash-deduplicate and append unique new sequences
    with gzip.open(args.new_seqs, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_hash = hashlib.md5(str(record.seq).encode()).hexdigest()
            record.id = record.id.split(",")[0]

            seq_ids.append(record.id)
            lens.append(len(record))
            gcs.append(gc_fraction(str(record.seq)))

            if seq_hash not in hash_to_id:
                hash_to_id[seq_hash] = record.id
                out_handle.write(record.format("fasta"))

            rep_ids.append(hash_to_id[seq_hash])

print(f"{len(hash_to_id)} unique new sequences out of {len(seq_ids)} total new sequences")

out_df = pd.DataFrame({"seq_id": seq_ids, "copy_of": rep_ids, "seq_gc": gcs, "seq_len": lens})
out_df.to_csv(new_hashed_file, index=False)
