#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scripts/run_test.sh
# Requires: VCLUST and (optionally) THREADS as environment variables.
# Expects vjdb1_merged_reps.fna.gz, vjdb1_merged_reps.csv and
# vjdb1_test_full_dataset.fna.gz in the repo root.

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

# --- Clashing IDs must be rejected ---
if python3 bin/5_prepare_update.py \
    --new-seqs vjdb1_test_clash.fna.gz \
    --existing-reps vjdb1_merged_reps.fna.gz \
    --existing-csv vjdb1_merged_reps.csv \
    --outdir "${TMPDIR}_clash" 2>/dev/null; then
    echo "FAIL: clashing ID was not rejected" >&2
    exit 1
fi
rm -rf "${TMPDIR}_clash"
echo "OK: clashing ID rejected"

# --- Step 1: Prepare combined FASTA ---
python3 bin/5_prepare_update.py \
    --new-seqs vjdb1_test_new_sequences.fna.gz \
    --existing-reps vjdb1_merged_reps.fna.gz \
    --existing-csv vjdb1_merged_reps.csv \
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
    --leiden-tsv "$TMPDIR/clusterreps_leiden.tsv" \
    --new-hashed-csv "$TMPDIR/new_hashed_reps.csv" \
    --combined-fasta "$COMBINED" \
    --existing-csv vjdb1_merged_reps.csv \
    --output-prefix vjdb1_test

# --- Check the expected outcome ---
python3 bin/check_test_update.py vjdb1_test

echo "Done. See vjdb1_test_new_clusters.csv and vjdb1_test_cluster_events.tsv."
