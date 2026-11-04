#!/usr/bin/env nextflow
nextflow.enable.types = true
// Scripts in modules/*/resources/usr/bin are added to the task PATH. Set here,
// not in nextflow.config, so that it also applies when the pipeline is included.
nextflow.enable.moduleBinaries = true

/*
 * VJDB Sequence Clustering Pipeline
 *
 *   nextflow run main.nf --mode initial --input <fasta.gz>
 *   nextflow run main.nf --mode update  --new_seqs <fasta.gz> --existing_reps <fasta.gz> --existing_csv <csv>
 *   nextflow run main.nf -profile test           (update on test data, with checks)
 *   nextflow run main.nf -profile test_initial   (initial clustering on test data, with checks)
 *
 * The pipeline can also be included in another pipeline (Nextflow >= 26.10):
 *
 *   include { params as VjdbParams ; workflow as VJDB_CLUSTERING } from './pipelines/vjdb-clustering/main.nf'
 *
 *   vjdb = VJDB_CLUSTERING(record(mode: 'initial', input: viral_fasta))
 *   vjdb.clusters_csv   // Value<Path>
 *   vjdb.reps_fasta     // Value<Path>
 *
 * Inputs are Value<Path> so that they can be given as a file on the command
 * line or as a process output from the including pipeline. Outputs are
 * declared in the output block; when the pipeline is included they become
 * the fields of the returned record and are not published.
 */

include { HASH_DEDUP }                              from './modules/hash_dedup'
include { LINCLUST }                                from './modules/linclust'
include { VCLUST_CLUSTER }                          from './modules/vclust'
include { VCLUST_CLUSTER as VCLUST_CLUSTER_MERGED } from './modules/vclust'
include { VCLUST_CLUSTER as VCLUST_CLUSTER_UPDATE } from './modules/vclust'
include { EXTRACT_MERGE_CHUNK_REPS }                from './modules/merge'
include { MERGE_FINAL_CLUSTERS }                    from './modules/merge'
include { EXTRACT_REPS_OF_REPS }                    from './modules/merge'
include { PREPARE_UPDATE }                          from './modules/update'
include { FINALIZE_UPDATE }                         from './modules/update'
include { CHECK_TEST_INITIAL }                      from './modules/test'
include { CHECK_TEST_UPDATE }                       from './modules/test'

params {
    // Mode: 'initial' (Phase 1) or 'update' (Phase 2)
    mode: String

    // Phase 1 inputs
    input: Value<Path>?            // FASTA.gz with all sequences
    block_size: Integer = 1000000  // sequences per chunk

    // Phase 2 inputs
    new_seqs: Value<Path>?         // FASTA.gz with new sequences
    existing_reps: Value<Path>?    // representative FASTA.gz from Phase 1 or a previous update
    existing_csv: Value<Path>?     // membership CSV from Phase 1 or a previous update
    prefix_update: String = 'vjdb1_new'

    // Clustering thresholds
    linclust_id: Float = 0.95      // MMseqs2 linclust --min-seq-id (Phase 1)
    ani: Float = 0.95              // vclust ANI threshold
    qcov: Float = 0.85             // vclust query coverage threshold

    // Tool commands (provided by the conda environment by default)
    mmseqs: String = 'mmseqs'
    vclust: String = 'vclust'

    // Check the outputs against the expected test results (set by the test profiles)
    check_test: Boolean = false
}

workflow INITIAL_CLUSTERING {
    take:
    input_fasta: Value<Path>
    block_size: Integer
    linclust_id: Float
    ani: Float
    qcov: Float
    mmseqs: String
    vclust: String

    main:
    // Step 1: hash-dedup and chunk
    chunks = HASH_DEDUP(input_fasta, block_size)

    // Step 2: linclust per chunk
    linclust_out = LINCLUST(chunks.flatMap { r -> r.fasta_chunks }, linclust_id, mmseqs)

    // Step 3: vclust per chunk
    chunk_leiden_tsvs = VCLUST_CLUSTER(linclust_out.map { r -> r.rep_fasta }, ani, qcov, vclust).collect()

    // Step 4: merge chunk-level reps
    merged_reps_fasta = EXTRACT_MERGE_CHUNK_REPS(
        chunk_leiden_tsvs,
        linclust_out.map { r -> r.rep_fasta }.collect()
    )

    // Step 5: vclust on merged reps
    merged_leiden_tsv = VCLUST_CLUSTER_MERGED(merged_reps_fasta, ani, qcov, vclust)

    // Step 6: final merge
    final_csv = MERGE_FINAL_CLUSTERS(
        chunks.map { r -> r.hashed_csv },
        merged_leiden_tsv,
        chunk_leiden_tsvs,
        linclust_out.map { r -> r.cluster_tsv }.collect()
    )
    final_fasta = EXTRACT_REPS_OF_REPS(merged_reps_fasta, merged_leiden_tsv)

    emit:
    csv: Value<Path> = final_csv
    fasta: Value<Path> = final_fasta
}

workflow UPDATE_CLUSTERING {
    take:
    new_seqs: Value<Path>
    existing_reps: Value<Path>
    existing_csv: Value<Path>
    prefix: String
    ani: Float
    qcov: Float
    vclust: String

    main:
    // Step 1: hash new seqs against the full dataset, check IDs, combine with existing reps
    prepared = PREPARE_UPDATE(new_seqs, existing_reps, existing_csv)
    combined_fasta = prepared.map { r -> r.combined_fasta }

    // Step 2: vclust on combined set
    leiden_tsv = VCLUST_CLUSTER_UPDATE(combined_fasta, ani, qcov, vclust)

    // Step 3: remap with stable IDs, report cluster events
    finalized = FINALIZE_UPDATE(
        leiden_tsv,
        prepared.map { r -> r.hashed_csv },
        existing_csv,
        combined_fasta,
        prefix
    )

    emit:
    csv: Value<Path> = finalized.map { r -> r.merged_csv }
    fasta: Value<Path> = finalized.map { r -> r.rep_fasta }
    new_clusters: Value<Path> = finalized.map { r -> r.new_clusters }
    events: Value<Path> = finalized.map { r -> r.events }
}

workflow {
    main:
    if (params.mode == 'initial') {
        if (!params.input) error "--mode initial needs --input"
    } else if (params.mode == 'update') {
        def missing = [
            params.new_seqs ? null : 'new_seqs',
            params.existing_reps ? null : 'existing_reps',
            params.existing_csv ? null : 'existing_csv'
        ].findAll { name -> name != null }
        if (missing) error "--mode update needs " + missing.collect { name -> "--${name}" }.join(', ')
    } else {
        error "Set --mode to 'initial' or 'update', or use -profile test / test_initial"
    }

    initial = params.mode == 'initial'
        ? INITIAL_CLUSTERING(
            params.input,
            params.block_size,
            params.linclust_id,
            params.ani,
            params.qcov,
            params.mmseqs,
            params.vclust
        )
        : null

    update = params.mode == 'update'
        ? UPDATE_CLUSTERING(
            params.new_seqs,
            params.existing_reps,
            params.existing_csv,
            params.prefix_update,
            params.ani,
            params.qcov,
            params.vclust
        )
        : null

    if (params.check_test && initial) {
        CHECK_TEST_INITIAL(
            initial.csv,
            initial.fasta,
            file("${moduleDir}/test_data/vjdb1_merged_reps.csv", checkIfExists: true)
        )
    }
    if (params.check_test && update) {
        CHECK_TEST_UPDATE(
            update.csv,
            update.new_clusters,
            update.events,
            params.existing_csv,
            params.prefix_update
        )
    }

    publish:
    clusters_csv   = initial?.csv ?: update.csv
    reps_fasta     = initial?.fasta ?: update.fasta
    new_clusters   = update?.new_clusters
    cluster_events = update?.events
}

output {
    // Membership CSV: <prefix>_merged_reps.csv
    clusters_csv: Value<Path> {
        path '.'
    }
    // Representative sequences: <prefix>_merged_reps.fna.gz
    reps_fasta: Value<Path> {
        path '.'
    }
    // Update mode only: clusters created by the update
    new_clusters: Value<Path?> {
        path '.'
    }
    // Update mode only: merges, splits and new clusters
    cluster_events: Value<Path?> {
        path '.'
    }
}
