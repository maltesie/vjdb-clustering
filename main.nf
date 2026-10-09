#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * VJDB Sequence Clustering Pipeline
 *
 *   nextflow run main.nf --mode initial --input <fasta.gz>
 *   nextflow run main.nf --mode update  --new_seqs <fasta.gz> --existing_reps <fasta.gz> --existing_csv <csv>
 *   nextflow run main.nf -profile test           (update on test data, with checks)
 *   nextflow run main.nf -profile test_initial   (initial clustering on test data, with checks)
 *
 * All parameter defaults are defined in nextflow.config.
 */

include { HASH_DEDUP }                   from './modules/hash_dedup'
include { LINCLUST }                     from './modules/linclust'
include { VCLUST_CLUSTER }                              from './modules/vclust'
include { VCLUST_CLUSTER as VCLUST_CLUSTER_MERGED }      from './modules/vclust'
include { VCLUST_CLUSTER as VCLUST_CLUSTER_UPDATE }      from './modules/vclust'
include { EXTRACT_MERGE_CHUNK_REPS }     from './modules/merge'
include { MERGE_FINAL_CLUSTERS }         from './modules/merge'
include { EXTRACT_REPS_OF_REPS }         from './modules/merge'
include { PREPARE_UPDATE }               from './modules/update'
include { FINALIZE_UPDATE }              from './modules/update'
include { CHECK_TEST_INITIAL }           from './modules/test'
include { CHECK_TEST_UPDATE }            from './modules/test'

workflow INITIAL_CLUSTERING {
    take:
        input_fasta

    main:
        // Step 1: hash-dedup and chunk
        chunks = HASH_DEDUP(input_fasta, params.block_size)

        // Step 2: linclust per chunk
        linclust_out = LINCLUST(chunks.fasta_chunks.flatten())

        // Step 3: vclust per chunk
        vclust_chunk_out = VCLUST_CLUSTER(linclust_out.rep_fasta, params.ani, params.qcov)

        // Step 4: merge chunk-level reps
        merged_reps_fasta = EXTRACT_MERGE_CHUNK_REPS(
            vclust_chunk_out.leiden_tsv.collect(),
            linclust_out.rep_fasta.collect()
        )

        // Step 5: vclust on merged reps
        vclust_merged_out = VCLUST_CLUSTER_MERGED(merged_reps_fasta, params.ani, params.qcov)

        // Step 6: final merge
        final_csv = MERGE_FINAL_CLUSTERS(
            chunks.hashed_csv,
            vclust_merged_out.leiden_tsv,
            vclust_chunk_out.leiden_tsv.collect(),
            linclust_out.cluster_tsv.collect()
        )
        final_fasta = EXTRACT_REPS_OF_REPS(
            merged_reps_fasta,
            vclust_merged_out.leiden_tsv
        )

    emit:
        csv   = final_csv
        fasta = final_fasta
}

workflow UPDATE_CLUSTERING {
    take:
        new_seqs
        existing_reps
        existing_csv

    main:
        // Step 1: hash new seqs against the full dataset, check IDs, combine with existing reps
        prepared = PREPARE_UPDATE(new_seqs, existing_reps, existing_csv)

        // Step 2: vclust on combined set
        vclust_out = VCLUST_CLUSTER_UPDATE(prepared.combined_fasta, params.ani, params.qcov)

        // Step 3: remap with stable IDs, report cluster events
        finalized = FINALIZE_UPDATE(
            vclust_out.leiden_tsv,
            prepared.hashed_csv,
            existing_csv,
            prepared.combined_fasta
        )

    emit:
        csv          = finalized.merged_csv
        fasta        = finalized.rep_fasta
        new_clusters = finalized.new_clusters
        events       = finalized.events
}

workflow {
    if (params.mode == 'initial') {
        if (!params.input) error "--mode initial needs --input"
        INITIAL_CLUSTERING(Channel.fromPath(params.input, checkIfExists: true))
        if (params.check_test) {
            CHECK_TEST_INITIAL(
                INITIAL_CLUSTERING.out.csv,
                INITIAL_CLUSTERING.out.fasta,
                Channel.fromPath("${projectDir}/test_data/vjdb1_merged_reps.csv", checkIfExists: true)
            )
        }
    } else if (params.mode == 'update') {
        def missing = ['new_seqs', 'existing_reps', 'existing_csv'].findAll { !params[it] }
        if (missing) error "--mode update needs " + missing.collect { "--${it}" }.join(', ')
        def existing_csv = Channel.fromPath(params.existing_csv, checkIfExists: true)
        UPDATE_CLUSTERING(
            Channel.fromPath(params.new_seqs, checkIfExists: true),
            Channel.fromPath(params.existing_reps, checkIfExists: true),
            existing_csv
        )
        if (params.check_test) {
            CHECK_TEST_UPDATE(
                UPDATE_CLUSTERING.out.csv,
                UPDATE_CLUSTERING.out.new_clusters,
                UPDATE_CLUSTERING.out.events,
                existing_csv
            )
        }
    } else {
        error "Set --mode to 'initial' or 'update', or use -profile test / test_initial"
    }
}
