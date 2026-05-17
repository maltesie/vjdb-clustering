#!/usr/bin/env python3

import argparse
import glob
import gzip
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Extract and merge vclust cluster representatives across chunks.")
parser.add_argument("--leiden-tsvs", nargs="+", default=None, help="Per-chunk Leiden cluster TSV files (default: glob chunks/)")
parser.add_argument("--rep-fastas", nargs="+", default=None, help="Per-chunk representative FASTA.gz files (default: glob chunks/)")
parser.add_argument("-o", "--output", default="vjdb1_vclust_reps_merged.fna.gz", help="Output merged FASTA.gz")
args = parser.parse_args()

# Fall back to glob if no explicit files given (backward compat with bash scripts)
if args.leiden_tsvs is None:
    args.leiden_tsvs = sorted(glob.glob("chunks/vjdb1_[0-9]*/clusterreps_leiden.tsv"))
if args.rep_fastas is None:
    args.rep_fastas = sorted(glob.glob("chunks/vjdb1_[0-9]*/linclust_rep_seq.fasta.gz"))

leiden_files = sorted(args.leiden_tsvs)
fasta_files = sorted(args.rep_fastas)

if len(leiden_files) != len(fasta_files):
    raise ValueError(f"Mismatch: {len(leiden_files)} leiden files vs {len(fasta_files)} fasta files")

print(f"Processing {len(leiden_files)} chunks")

with gzip.open(args.output, "wt", 2) as out_handle:
    for leiden_file, fasta_file in zip(leiden_files, fasta_files):
        clusters = set(pd.read_csv(leiden_file, sep="\t", dtype=str)["cluster"])
        print(f"{leiden_file}: {len(clusters)} clusters")
        with gzip.open(fasta_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in clusters:
                    out_handle.write(record.format("fasta"))
