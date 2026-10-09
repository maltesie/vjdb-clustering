#!/usr/bin/env python3

import argparse
import hashlib
import os, gzip
import pandas as pd
from collections import Counter
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

parser = argparse.ArgumentParser(description="Hash-deduplicate and chunk VJDB sequences.")
parser.add_argument("--input", required=True, help="Path to input FASTA.gz")
parser.add_argument("--outdir", default="chunks",
                    help="Output directory for chunks (default: ./chunks)")
parser.add_argument("--block-size", type=int, default=1000000, help="Sequences per chunk (default: 1000000)")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

hash_to_id = dict()
seq_ids = []
rep_ids = []
lens = []
gcs = []
out_handle = None

with gzip.open(args.input, "rt") as handle:
    for i, record in enumerate(SeqIO.parse(handle, "fasta")):

        if i % args.block_size == 0:
            if out_handle is not None:
                out_handle.close()
            fname = f"vjdb1_{i // args.block_size}.fna.gz"
            out_handle = gzip.open(os.path.join(args.outdir, fname), "wt", 2)

        seq_hash = hashlib.md5(str(record.seq).encode()).hexdigest()

        record.id = record.id.split(',')[0]

        seq_ids.append(record.id)
        lens.append(len(record))
        gcs.append(gc_fraction(str(record.seq)))

        if seq_hash not in hash_to_id:
            hash_to_id[seq_hash] = record.id
            out_handle.write(record.format("fasta"))

        rep_ids.append(hash_to_id[seq_hash])

    if out_handle is not None:
        out_handle.close()

print(Counter(seq_ids).most_common(10))
print(f"{len(hash_to_id)} unique sequences out of {i + 1} sequences")
out_df = pd.DataFrame({"vjdb_id": seq_ids, "copy_of": rep_ids, "seq_gc": gcs, "seq_len": lens})
out_df.to_csv(os.path.join(args.outdir, "vjdb1_hashed_reps.csv"), index=False)
