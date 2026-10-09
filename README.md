# VJDB Sequence Clustering

Nextflow pipeline for hierarchical dereplication and clustering of the VirJenDB viral genome database, using hash deduplication, MMseqs2 linclust and vclust (Leiden community detection). It has two modes:

- **`initial`** clusters a full dataset from scratch.
- **`update`** adds new sequences to an existing clustering and reports which clusters merged, grew or are new.

## Requirements

- [Nextflow](https://www.nextflow.io/) 23.04 or newer
- [micromamba](https://mamba.readthedocs.io/), which the pipeline uses to install vclust, MMseqs2, biopython and pandas from `environment.yml`. If those tools are already on your `PATH`, add `-profile no_conda` instead.

## Quick start

Run the update mode on the bundled test data. It checks the results automatically:

```bash
nextflow run maltesie/vjdb-clustering -profile test
```

or, from a clone of this repository:

```bash
nextflow run main.nf -profile test
```

The run ends with `OK: all checks passed` and writes its outputs to `results_test/`.

---

## Phase 1 — Initial clustering

Clusters the full dataset from scratch. Run it once to produce the baseline cluster representatives and membership table.

```bash
nextflow run maltesie/vjdb-clustering --mode initial --input vjdb.fasta.gz --outdir results
```

```
input FASTA.gz
  │
  ├─ HASH_DEDUP             (1_prepare_clustering.py)       hash-dedup, split into chunks
  ├─ LINCLUST               (mmseqs easy-linclust)          per chunk, in parallel
  ├─ VCLUST_CLUSTER         (vclust)                        Leiden clusters per chunk
  ├─ EXTRACT_MERGE_CHUNK_REPS (2_extract_and_merge_vclust_reps.py)
  ├─ VCLUST_CLUSTER_MERGED  (vclust)                        global Leiden clusters
  ├─ MERGE_FINAL_CLUSTERS   (3_merge_clusters.py)        ──► vjdb1_merged_reps.csv
  └─ EXTRACT_REPS_OF_REPS   (4_extract_reps_of_reps.py)  ──► vjdb1_merged_reps.fna.gz
```

**Outputs** in `--outdir`:

| File | Content |
|------|---------|
| `vjdb1_merged_reps.csv` | Membership table: `vjdb_id`, `copy_of`, `seq_gc`, `seq_len`, `seq_hash`, `cluster_rep` |
| `vjdb1_merged_reps.fna.gz` | Representative sequence of each cluster |

These two files are the input for Phase 2.

---

## Phase 2 — Update clustering

Adds new sequences to an existing clustering without starting from scratch. Only the existing representatives and the genuinely new sequences are re-clustered.

```bash
nextflow run maltesie/vjdb-clustering --mode update \
    --new_seqs new_sequences.fasta.gz \
    --existing_reps results/vjdb1_merged_reps.fna.gz \
    --existing_csv results/vjdb1_merged_reps.csv \
    --prefix_update vjdb2 --outdir results_update
```

```
new FASTA.gz + existing reps + existing CSV
  │
  ├─ PREPARE_UPDATE         (5_prepare_update.py)   hash against the full dataset, check IDs
  ├─ VCLUST_CLUSTER_UPDATE  (vclust)                Leiden clusters of reps + new sequences
  └─ FINALIZE_UPDATE        (6_finalize_update.py)  stable IDs, events, updated tables
```

**Outputs** in `--outdir`:

| File | Content |
|------|---------|
| `<prefix>_merged_reps.csv` | Updated membership table; `vjdb_id` is unique |
| `<prefix>_merged_reps.fna.gz` | Representative sequences of all clusters |
| `<prefix>_new_clusters.csv` | Every record of the new batch with its `status` and `cluster_rep` |
| `<prefix>_cluster_events.tsv` | One row per changed cluster: `merged`, `grown` or `new` |

The updated CSV and FASTA can be used as `--existing_csv` and `--existing_reps` for the next update.

**Status of new records:** `new_unique` (clustered with vclust), `new_copy` (identical to another new sequence), `existing_copy` (identical to a sequence already in the dataset), `resubmitted` (same ID and sequence already in the dataset, not added again), `duplicate_in_batch` (same ID and sequence repeated in the batch, added once). An ID that is already used for a different sequence stops the update with a list of the clashing IDs.

**Stable cluster IDs:** a cluster keeps its representative ID across updates. When clusters merge, the merged cluster keeps the representative of the largest old cluster; the others are listed in `old_cluster_reps` of the events file. Splits can't happen, because members of an old cluster always follow their representative.

**Sanity checks:** the update stops if the representative FASTA and CSV disagree, if a sequence given to vclust is missing from its output (vclust lists singletons as their own cluster, so this indicates an ID mismatch), or if the final table would contain a duplicated `vjdb_id`.

---

## Parameters

All defaults are set in `nextflow.config`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mode` | – | `initial` or `update` |
| `--input` | – | Phase 1: FASTA.gz with all sequences |
| `--block_size` | 1000000 | Phase 1: sequences per chunk |
| `--new_seqs` | – | Phase 2: FASTA.gz with new sequences |
| `--existing_reps` | – | Phase 2: representative FASTA.gz from Phase 1 or a previous update |
| `--existing_csv` | – | Phase 2: membership CSV from Phase 1 or a previous update |
| `--prefix_update` | `vjdb1_new` | Phase 2: prefix of the output files |
| `--outdir` | `results` | Output directory |
| `--threads` | 10 | CPUs for MMseqs2 and vclust |
| `--linclust_id` | 0.95 | MMseqs2 linclust `--min-seq-id` (Phase 1) |
| `--ani` | 0.95 | vclust ANI threshold for Leiden clustering |
| `--qcov` | 0.85 | vclust query coverage threshold |
| `--mmseqs` | `mmseqs` | MMseqs2 command |
| `--vclust` | `vclust` | vclust command |

## Running on your system

- **Profiles:** `-profile test` (bundled test data), `no_conda` (tools from your `PATH`), `slurm` (example for a cluster; adjust the queue in `nextflow.config`). Profiles can be combined, e.g. `-profile test,no_conda`.
- **Resources:** each process asks for 16 GB of memory and `--threads` CPUs by default. Nextflow refuses to start a process that asks for more than your machine has, so on a small machine lower them in your own config and pass it with `-c my.config`:
  ```groovy
  process.memory = '4 GB'
  ```
- **Adjusting memory:** the requirements are set in [`nextflow.config`](nextflow.config): the default for all processes in the `process` block, and the MMseqs2 and vclust values in the `slurm` profile. Edit them there, or override a single step in your own config, e.g. `process { withName: 'VCLUST_CLUSTER.*' { memory = '300 GB' } }`.
- **Resuming:** add `-resume` to continue an interrupted run; finished steps are taken from the `work/` directory.
- **Files:** results are copied to `--outdir`. Intermediate files stay in `work/`, which you can delete after a successful run.

## Test data

`test_data/` holds a small clustering (100 clusters, 667 sequences) and inputs for an update that covers every case: new sequences, copies, a resubmission, cluster growth and a merge caused by a chimera, plus a file with a clashing ID that the update must reject.

- `vjdb1_test_full_dataset.fna.gz` is a synthetic full dataset matching the test CSV (members are ~1% mutated copies of their representative), built by `bin/make_synthetic_full_dataset.py`. It can also be used as `--input` for Phase 1.
- `vjdb1_test_new_sequences.fna.gz` and `vjdb1_test_clash.fna.gz` are built by `bin/make_test_data.py`; `bin/check_test_update.py` checks the outputs.

## Repository structure

```
main.nf            Pipeline entry point
nextflow.config    Parameters, resources and profiles
modules/           Nextflow processes
environment.yml    Conda environment for all processes
bin/               Python scripts used by the processes
test_data/         Test datasets
scripts/           Bash reference implementation
```

## Bash reference implementation

The original bash scripts in `scripts/` still run the same Python steps without Nextflow. They need a Python environment (`bash scripts/setup_environment.sh`) and the tools set as environment variables:

```bash
export VCLUST=vclust MMSEQS=mmseqs THREADS=24
bash scripts/run_initial_clustering.sh <input.fasta.gz>
bash scripts/run_update_clustering.sh <new_sequences.fasta.gz> [workdir]
bash scripts/run_test.sh
```

## Background

The bash scripts were turned into this Nextflow pipeline at the VirJenDB hackathon at ViBiom 2026.
