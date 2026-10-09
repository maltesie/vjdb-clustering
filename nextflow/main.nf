#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * VJDB Sequence Clustering Pipeline
 *
 * Two entry points:
 *   nextflow run main.nf --mode initial --input <fasta.gz>
 *   nextflow run main.nf --mode update  --new_seqs <fasta.gz> --existing_reps <fasta.gz> --existing_csv <csv>
 */

params.mode           = null   // 'initial' or 'update'
params.input          = null   // Phase 1: input FASTA.gz
params.block_size     = 1000000
params.new_seqs       = null   // Phase 2: new sequences FASTA.gz
params.existing_reps  = null   // Phase 2: existing representative FASTA.gz
params.existing_csv   = null   // Phase 2: existing merged CSV
params.threads        = 24
params.outdir         = 'results'
params.prefix_update = params.prefix_update ?: 'vjdb1_new'

// Clustering thresholds (single source of truth)
params.ani            = 0.95
params.qcov           = 0.85

// Tool paths — override in nextflow.config or on the command line
params.mmseqs         = 'mmseqs'
params.vclust         = 'vclust.py'

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
        // Step 1: hash-dedup new, combine with existing reps
        prepared = PREPARE_UPDATE(new_seqs, existing_reps)

        // Step 2: vclust on combined set
        vclust_out = VCLUST_CLUSTER_UPDATE(prepared.combined_fasta, params.ani, params.qcov)

        // Step 3: remap and finalize
        finalized = FINALIZE_UPDATE(
            prepared.update_tmp_dir,
            vclust_out.leiden_tsv,
            prepared.hashed_csv,
            existing_csv,
            prepared.combined_fasta
        )

    emit:
        csv   = finalized.merged_csv
        fasta = finalized.rep_fasta
}

workflow {
    if (params.mode == 'initial') {
        INITIAL_CLUSTERING(Channel.fromPath(params.input))
    } else if (params.mode == 'update') {
        UPDATE_CLUSTERING(
            Channel.fromPath(params.new_seqs),
            Channel.fromPath(params.existing_reps),
            Channel.fromPath(params.existing_csv)
        )
    } else {
        error "Set --mode to 'initial' or 'update'"
    }
}
