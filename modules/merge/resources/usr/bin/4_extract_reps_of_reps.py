#!/usr/bin/env python3

import argparse
import gzip
import os
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Extract final representative sequences.")
parser.add_argument("--input-fasta", default=None, help="Merged representative FASTA.gz (default: chunks/vjdb1_vclust_reps_merged.fna.gz)")
parser.add_argument("--leiden-tsv", default=None, help="Leiden TSV from merged-rep vclust run (default: chunks/vjdb1_merged/clusterreps_leiden.tsv)")
parser.add_argument("-o", "--output", default="vjdb1_merged_reps.fna.gz", help="Output FASTA.gz")
args = parser.parse_args()

# Fall back to hardcoded paths if no explicit files given (backward compat)
if args.input_fasta is None:
    args.input_fasta = os.path.join("chunks", "vjdb1_vclust_reps_merged.fna.gz")
if args.leiden_tsv is None:
    args.leiden_tsv = os.path.join("chunks", "vjdb1_merged", "clusterreps_leiden.tsv")

clusters = set(pd.read_csv(args.leiden_tsv, sep="\t", dtype=str)["cluster"])

with gzip.open(args.output, "wt", 2) as out_handle:
    with gzip.open(args.input_fasta, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in clusters:
                out_handle.write(record.format("fasta"))
