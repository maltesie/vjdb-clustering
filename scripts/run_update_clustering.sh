#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scripts/run_update_clustering.sh <new_sequences.fasta.gz> [workdir]
# Requires: VCLUST and (optionally) THREADS as environment variables.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(dirname "$SCRIPT_DIR")"  # repo root

NEW_SEQS="${1:?Usage: bash scripts/run_update_clustering.sh <new_sequences.fasta.gz> [workdir]}"
WORKDIR="${2:-.}"

THREADS="${THREADS:-24}"
VCLUST="${VCLUST:?Set VCLUST to the vclust command (e.g. 'python3 /path/to/vclust.py')}"

EXISTING_REPS="${WORKDIR}/vjdb1_merged_reps.fna.gz"
EXISTING_CSV="${WORKDIR}/vjdb1_merged_reps.csv"
TMPDIR="${WORKDIR}/update_tmp"
OUTPUT_PREFIX="${WORKDIR}/vjdb1_new"

# --- Activate Python environment ---
if [ ! -d "py3env" ]; then
    echo "Error: py3env not found. Run scripts/setup_environment.sh first." >&2
    exit 1
fi
source py3env/bin/activate

# --- Step 1: Hash-dedup new sequences, combine with existing reps ---
python3 bin/5_prepare_update.py \
    --new-seqs "$NEW_SEQS" \
    --existing-reps "$EXISTING_REPS" \
    --outdir "$TMPDIR"

# --- Step 2: Run vclust on the combined set ---
COMBINED="${TMPDIR}/combined_for_vclust.fna.gz"

$VCLUST prefilter -i "$COMBINED" -o "$TMPDIR/fltr.txt" -t "$THREADS" \
    --min-ident 0.95 --batch-size 100000 --kmers-fraction 0.2 --max-seqs 2000

$VCLUST align -i "$COMBINED" -o "$TMPDIR/ani.tsv" -t "$THREADS" \
    --filter "$TMPDIR/fltr.txt" --filter-threshold 0.95 \
    --out-ani 0.95 --out-qcov 0.85

$VCLUST cluster -i "$TMPDIR/ani.tsv" -o "$TMPDIR/clusterreps_leiden.tsv" \
    --ids "$TMPDIR/ani.ids.tsv" --algorithm leiden --metric ani \
    --ani 0.95 --qcov 0.85 --out-repr

rm -f "$TMPDIR/ani.tsv" "$TMPDIR/fltr.txt" "$TMPDIR/ani.ids.tsv"

# --- Step 3: Remap clusters, update CSVs, extract new reps ---
python3 bin/6_finalize_update.py \
    --outdir "$TMPDIR" \
    --existing-csv "$EXISTING_CSV" \
    --output-prefix "$OUTPUT_PREFIX"

echo "Done. Outputs:"
echo "  Updated full CSV:    ${OUTPUT_PREFIX}_merged_reps.csv"
echo "  New sequences CSV:   ${OUTPUT_PREFIX}_new_clusters.csv"
echo "  Updated rep FASTA:   ${OUTPUT_PREFIX}_merged_reps.fna.gz"

# --- Step 4: Validate consistency ---
echo ""
echo "=== Validating clustering consistency ==="
python3 bin/validate_update.py \
    --before "$EXISTING_CSV" \
    --after "${OUTPUT_PREFIX}_merged_reps.csv" \
    --new "${OUTPUT_PREFIX}_new_clusters.csv"
