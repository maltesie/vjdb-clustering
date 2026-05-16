#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scripts/run_test.sh
# Requires: VCLUST and (optionally) THREADS as environment variables.
# Expects vjdb1_merged_reps.fna.gz and vjdb1_merged_reps.csv in the repo root.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$(dirname "$SCRIPT_DIR")"  # repo root

THREADS="${THREADS:-24}"
VCLUST="${VCLUST:?Set VCLUST to the vclust command (e.g. 'python3 /path/to/vclust.py')}"
TMPDIR="update_tmp"

# --- Activate Python environment ---
if [ ! -d "py3env" ]; then
    echo "Error: py3env not found. Run scripts/setup_environment.sh first." >&2
    exit 1
fi
source py3env/bin/activate

# --- Generate test data ---
python3 bin/make_test_data.py

echo "Test data generated. Running update pipeline..."

# --- Step 1: Prepare combined FASTA ---
python3 bin/5_prepare_update.py \
    --new-seqs vjdb1_test_new_sequences.fna.gz \
    --existing-reps vjdb1_merged_reps.fna.gz \
    --outdir "$TMPDIR"

# --- Step 2: Run vclust on combined set ---
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

# --- Step 3: Finalize ---
python3 bin/6_finalize_update.py \
    --outdir "$TMPDIR" \
    --existing-csv vjdb1_merged_reps.csv \
    --output-prefix vjdb1_test

echo ""
echo "=== Validating clustering consistency ==="
python3 bin/validate_update.py \
    --before vjdb1_merged_reps.csv \
    --after vjdb1_test_merged_reps.csv \
    --new vjdb1_test_new_clusters.csv
