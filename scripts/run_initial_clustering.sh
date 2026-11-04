#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scripts/run_initial_clustering.sh <input.fasta.gz>
# Requires: VCLUST, MMSEQS, and (optionally) THREADS as environment variables.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(dirname "$SCRIPT_DIR")"  # repo root

INPUT_FASTA="${1:?Usage: bash scripts/run_initial_clustering.sh <input.fasta.gz>}"
THREADS="${THREADS:-24}"
VCLUST="${VCLUST:?Set VCLUST to the vclust command (e.g. 'python3 /path/to/vclust.py')}"
MMSEQS="${MMSEQS:?Set MMSEQS to the mmseqs binary (e.g. '/path/to/mmseqs')}"

# --- Activate Python environment ---
if [ ! -d "py3env" ]; then
    echo "Error: py3env not found. Run scripts/setup_environment.sh first." >&2
    exit 1
fi
source py3env/bin/activate

mkdir -p chunks

# --- Step 1: Hash-dedup and split into chunks ---
python3 modules/hash_dedup/resources/usr/bin/1_prepare_clustering.py --input "$INPUT_FASTA"

# --- Step 2: Linclust per chunk ---
NUM_CHUNKS=$(ls chunks/vjdb1_*.fna.gz 2>/dev/null | wc -l)
echo "Detected $NUM_CHUNKS chunks"
for ((i=0; i<NUM_CHUNKS; i++)); do
    zippedfile="chunks/vjdb1_${i}.fna.gz"
    outdir="chunks/vjdb1_${i}/"
    mkdir -p "$outdir"

    [[ -f "${outdir}linclust_rep_seq.fasta.gz" ]] && continue

    tmpdir="lintmp"
    mkdir -p "$tmpdir"

    $MMSEQS easy-linclust --threads "$THREADS" --min-seq-id 0.95 \
        "$zippedfile" "${outdir}linclust" "$tmpdir"

    rm -rf "$tmpdir"
    rm -f "${outdir}linclust_all_seqs.fasta"
    gzip -2 "${outdir}linclust_rep_seq.fasta"
done

# --- Step 3: Vclust per chunk ---
for ((i=0; i<NUM_CHUNKS; i++)); do
    infile="chunks/vjdb1_${i}/linclust_rep_seq.fasta.gz"
    outdir="chunks/vjdb1_${i}"
    mkdir -p "$outdir"

    [[ -f "$outdir/clusterreps_leiden.tsv" ]] && continue

    $VCLUST prefilter -i "$infile" -o "$outdir/fltr.txt" -t "$THREADS" \
        --min-ident 0.95 --batch-size 100000 --kmers-fraction 0.2 --max-seqs 2000
    $VCLUST align -i "$infile" -o "$outdir/ani.tsv" -t "$THREADS" \
        --filter "$outdir/fltr.txt" --filter-threshold 0.95 \
        --out-ani 0.95 --out-qcov 0.85
    $VCLUST cluster -i "$outdir/ani.tsv" -o "$outdir/clusterreps_leiden.tsv" \
        --ids "$outdir/ani.ids.tsv" --algorithm leiden --metric ani \
        --ani 0.95 --qcov 0.85 --out-repr

    rm -f "$outdir/ani.tsv" "$outdir/fltr.txt" "$outdir/ani.ids.tsv"
done

# --- Step 4: Merge chunk-level reps ---
python3 modules/merge/resources/usr/bin/2_extract_and_merge_vclust_reps.py

# --- Step 5: Vclust on merged reps ---
infile="chunks/vjdb1_vclust_reps_merged.fna.gz"
outdir="chunks/vjdb1_merged"
mkdir -p "$outdir"

$VCLUST prefilter -i "$infile" -o "$outdir/fltr.txt" -t "$THREADS" \
    --min-ident 0.95 --batch-size 100000 --kmers-fraction 0.2 --max-seqs 2000
$VCLUST align -i "$infile" -o "$outdir/ani.tsv" -t "$THREADS" \
    --filter "$outdir/fltr.txt" --filter-threshold 0.95 \
    --out-ani 0.95 --out-qcov 0.85
$VCLUST cluster -i "$outdir/ani.tsv" -o "$outdir/clusterreps_leiden.tsv" \
    --ids "$outdir/ani.ids.tsv" --algorithm leiden --metric ani \
    --ani 0.95 --qcov 0.85 --out-repr

rm -f "$outdir/ani.tsv" "$outdir/fltr.txt" "$outdir/ani.ids.tsv"

# --- Step 6: Final merge and extract reps ---
python3 modules/merge/resources/usr/bin/3_merge_clusters.py
python3 modules/merge/resources/usr/bin/4_extract_reps_of_reps.py

echo "Done. Outputs: vjdb1_merged_reps.csv, vjdb1_merged_reps.fna.gz"
