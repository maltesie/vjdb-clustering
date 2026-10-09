# VJDB Sequence Clustering

Hierarchical dereplication and clustering pipeline for the VirJenDB viral genome database, using hash deduplication, MMseqs2 linclust, and vclust (Leiden community detection).

## Dependencies

- Python 3.9+ with `biopython` and `pandas`
- [MMseqs2](https://github.com/soedinglab/MMseqs2) (v14+, only needed for Phase 1)
- [vclust](https://github.com/refresh-bio/vclust) (v1.3+)

## Quick Start

```bash
bash scripts/setup_environment.sh     # creates py3env/

export VCLUST="python3 /path/to/vclust.py"
export MMSEQS="/path/to/mmseqs"       # only needed for Phase 1
export THREADS=24                     # optional, defaults to 24
```

## Repo Structure

```
bin/           Python scripts (core logic)
scripts/       Bash reference implementation
nextflow/      Nextflow pipeline skeleton
test_data/     Test datasets
```

## Hackathon

See [CONTRIBUTING.md](CONTRIBUTING.md) for task descriptions, branch conventions, and getting started. The hackathon focuses on Phase 2 (update clustering) and the Nextflow implementation — MMseqs2 is not required.

---

## Phase 1 — Initial Clustering

Clusters the full VJDB dataset from scratch. Run once to produce the baseline cluster representatives and membership CSV.

```
VJDB FASTA
  │
  ├─ 1_prepare_clustering.py ──► chunked FASTAs + hash-dedup CSV
  │
  ├─ mmseqs linclust (per chunk) ──► linclust representatives
  │
  ├─ vclust prefilter/align/cluster (per chunk)
  │        ──► Leiden clusters per chunk
  │
  ├─ 2_extract_and_merge_vclust_reps.py
  │        ──► merged representative FASTA across chunks
  │
  ├─ vclust prefilter/align/cluster (merged reps)
  │        ──► global Leiden clusters
  │
  ├─ 3_merge_clusters.py ──► vjdb1_merged_reps.csv
  │
  └─ 4_extract_reps_of_reps.py ──► vjdb1_merged_reps.fna.gz
```

```bash
bash scripts/run_initial_clustering.sh <input.fasta.gz>
```

**Outputs:** `vjdb1_merged_reps.csv` (full membership table) and `vjdb1_merged_reps.fna.gz` (cluster representative sequences).

---

## Phase 2 — Update Clustering (new sequences)

Integrates new sequences into the existing clusters without re-running from scratch. Only the existing representatives and the genuinely new sequences are re-clustered with vclust.

```
New FASTA + existing reps + existing CSV
  │
  ├─ 5_prepare_update.py ──► combined FASTA + new_hashed_reps.csv
  │     hashes new sequences against the full dataset (seq_hash column),
  │     stops on ID clashes
  │
  ├─ vclust prefilter/align/cluster (combined)
  │        ──► updated Leiden clusters
  │
  └─ 6_finalize_update.py ──► updated CSVs, cluster events, updated rep FASTA
```

```bash
bash scripts/run_update_clustering.sh <new_sequences.fasta.gz> [workdir]
# or
nextflow run nextflow/main.nf --mode update --new_seqs <new.fasta.gz> \
    --existing_reps vjdb1_merged_reps.fna.gz --existing_csv vjdb1_merged_reps.csv
```

**Outputs:**

| File | Content |
|------|---------|
| `<prefix>_merged_reps.csv` | Updated membership table; `vjdb_id` is unique |
| `<prefix>_merged_reps.fna.gz` | Representative sequences of all clusters |
| `<prefix>_new_clusters.csv` | Every record of the new batch with its `status` and `cluster_rep` |
| `<prefix>_cluster_events.tsv` | One row per changed cluster: `merged`, `grown` or `new` |

**Status of new records:** `new_unique` (clustered with vclust), `new_copy` (identical to another new sequence), `existing_copy` (identical to a sequence already in the dataset), `resubmitted` (same ID and sequence already in the dataset, not added again), `duplicate_in_batch` (same ID and sequence repeated in the batch, added once). An ID that is already used for a different sequence stops the update with a list of the clashing IDs.

**Stable cluster IDs:** a cluster keeps its representative ID across updates. When clusters merge, the merged cluster keeps the representative of the largest old cluster; the others are listed in `old_cluster_reps` of the events file. Splits can't happen, because members of an old cluster always follow their representative.

**Sanity checks:** the update stops if the rep FASTA and CSV disagree, if a sequence given to vclust is missing from its output (vclust lists singletons as their own cluster, so this indicates an ID mismatch), or if the final table would contain a duplicated `vjdb_id`.

---

## Testing

```bash
bash scripts/run_test.sh
```

Generates a batch covering every update case (new sequences, copies, a resubmission, cluster growth and a merge caused by a chimera), checks that a clashing ID is rejected, runs Phase 2 end-to-end and verifies the outputs with `bin/check_test_update.py`.

`vjdb1_test_full_dataset.fna.gz` is a synthetic full dataset matching the test CSV (members are ~1% mutated copies of their representative), built by `bin/make_synthetic_full_dataset.py`. It is only used to create test inputs.

---

## Clustering Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| linclust `--min-seq-id` | 0.95 | Pre-clustering identity threshold |
| vclust `--ani` | 0.95 | ANI threshold for Leiden clustering |
| vclust `--qcov` | 0.85 | Query coverage threshold |
| Leiden algorithm | community detection | Resolves transitive clusters |
