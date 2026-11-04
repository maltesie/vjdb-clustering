#!/usr/bin/env python3
"""
Prepare new sequences for clustering with the existing representatives.

Every new sequence is hashed and compared against the full existing dataset
(the seq_hash column of the existing membership CSV) and against the rest of
the new batch. Each new record gets a status:

  new_unique          new sequence not seen before -> goes to vclust
  new_copy            identical to an earlier sequence in this batch
  existing_copy       identical to a sequence already in the dataset
  resubmitted         ID and sequence already in the dataset -> not added again
  duplicate_in_batch  ID and sequence repeated within this batch -> not added again

An ID that is already used for a different sequence (in the dataset or in the
batch) is a clash: the script lists all clashes and stops before vclust runs.
"""

import argparse
import gzip
import hashlib
import os
import sys

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

parser = argparse.ArgumentParser(description="Prepare new sequences for clustering with existing representatives.")
parser.add_argument("--new-seqs", required=True, help="Path to new sequences FASTA.gz")
parser.add_argument("--existing-reps", required=True, help="Path to existing merged rep FASTA.gz (e.g. vjdb1_merged_reps.fna.gz)")
parser.add_argument("--existing-csv", required=True, help="Path to existing membership CSV with a seq_hash column (e.g. vjdb1_merged_reps.csv)")
parser.add_argument("--outdir", default=".", help="Output directory (default: current directory)")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)
combined_file = os.path.join(args.outdir, "combined_for_vclust.fna.gz")
new_hashed_file = os.path.join(args.outdir, "new_hashed_reps.csv")

# --- Existing dataset -------------------------------------------------------
existing_df = pd.read_csv(args.existing_csv, dtype=str)
missing_cols = {"vjdb_id", "copy_of", "seq_hash", "cluster_rep"} - set(existing_df.columns)
if missing_cols:
    sys.exit(f"Error: {args.existing_csv} is missing column(s) {sorted(missing_cols)}. "
             "Re-run Phase 1 so the CSV contains seq_hash.")
if not existing_df["vjdb_id"].is_unique:
    dups = existing_df.loc[existing_df["vjdb_id"].duplicated(), "vjdb_id"].unique()
    sys.exit(f"Error: existing CSV has {len(dups)} duplicated vjdb_id(s), e.g. {list(dups[:10])}")

existing_hash_by_id = dict(zip(existing_df["vjdb_id"], existing_df["seq_hash"]))
existing_copy_of_by_id = dict(zip(existing_df["vjdb_id"], existing_df["copy_of"]))
# hash -> canonical ID (copy_of) of the first existing sequence with that hash
existing_id_by_hash = dict(zip(existing_df["seq_hash"][::-1], existing_df["copy_of"][::-1]))
existing_reps = set(existing_df["cluster_rep"])

# --- New sequences ----------------------------------------------------------
rows = []
new_hash_by_id = {}
new_id_by_hash = {}
new_copy_of_by_id = {}
clashes = []
fasta_ids = set()

with gzip.open(combined_file, "wt", 2) as out_handle:
    # Existing cluster reps first; they must match the reps in the CSV
    with gzip.open(args.existing_reps, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            out_handle.write(record.format("fasta"))
            fasta_ids.add(record.id)
    if fasta_ids != existing_reps:
        sys.exit(f"Error: rep FASTA and CSV disagree: {len(fasta_ids - existing_reps)} FASTA IDs are not a "
                 f"cluster_rep in the CSV, {len(existing_reps - fasta_ids)} cluster_reps are missing from the FASTA")
    print(f"Wrote {len(fasta_ids)} existing cluster representatives")

    with gzip.open(args.new_seqs, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq)
            seq_hash = hashlib.md5(seq.encode()).hexdigest()
            seq_id = record.id.split(",")[0]
            record.id, record.description = seq_id, ""

            if seq_id in existing_hash_by_id:
                if existing_hash_by_id[seq_id] == seq_hash:
                    status, copy_of = "resubmitted", existing_copy_of_by_id[seq_id]
                else:
                    clashes.append((seq_id, "already in dataset with a different sequence"))
                    continue
            elif seq_id in new_hash_by_id:
                if new_hash_by_id[seq_id] == seq_hash:
                    status, copy_of = "duplicate_in_batch", new_copy_of_by_id[seq_id]
                else:
                    clashes.append((seq_id, "used twice in new batch for different sequences"))
                    continue
            elif seq_hash in existing_id_by_hash:
                status, copy_of = "existing_copy", existing_id_by_hash[seq_hash]
            elif seq_hash in new_id_by_hash:
                status, copy_of = "new_copy", new_id_by_hash[seq_hash]
            else:
                status, copy_of = "new_unique", seq_id
                new_id_by_hash[seq_hash] = seq_id
                out_handle.write(record.format("fasta"))

            new_hash_by_id[seq_id] = seq_hash
            new_copy_of_by_id.setdefault(seq_id, copy_of)
            rows.append({"vjdb_id": seq_id, "copy_of": copy_of, "seq_gc": gc_fraction(seq),
                         "seq_len": len(seq), "seq_hash": seq_hash, "status": status})

if clashes:
    os.remove(combined_file)
    print(f"Error: {len(clashes)} ID clash(es) in {args.new_seqs}:", file=sys.stderr)
    for seq_id, reason in clashes:
        print(f"  {seq_id}: {reason}", file=sys.stderr)
    sys.exit(1)

new_df = pd.DataFrame(rows, columns=["vjdb_id", "copy_of", "seq_gc", "seq_len", "seq_hash", "status"])
new_df.to_csv(new_hashed_file, index=False)
print(f"{len(new_df)} new records: " + ", ".join(f"{k}={v}" for k, v in new_df["status"].value_counts().items()))
