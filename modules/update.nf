process PREPARE_UPDATE {
    /*
     * Hash new sequences against the full existing dataset, check for ID
     * clashes and combine unique new sequences with the existing reps.
     * Wraps: bin/5_prepare_update.py
     */

    input:
    path new_seqs
    path existing_reps
    path existing_csv

    output:
    path "combined_for_vclust.fna.gz", emit: combined_fasta
    path "new_hashed_reps.csv",        emit: hashed_csv

    script:
    """
    5_prepare_update.py \
        --new-seqs ${new_seqs} \
        --existing-reps ${existing_reps} \
        --existing-csv ${existing_csv} \
        --outdir .
    """
}

process FINALIZE_UPDATE {
    /*
     * Remap clusters with stable IDs, report cluster events, write the
     * updated membership CSV and representative FASTA.
     * Wraps: bin/6_finalize_update.py
     */

    publishDir params.outdir, mode: 'copy'

    input:
    path leiden_tsv
    path new_hashed_csv
    path existing_csv
    path combined_fasta

    output:
    path "${params.prefix_update}_merged_reps.csv",    emit: merged_csv
    path "${params.prefix_update}_merged_reps.fna.gz", emit: rep_fasta
    path "${params.prefix_update}_new_clusters.csv",   emit: new_clusters
    path "${params.prefix_update}_cluster_events.tsv", emit: events

    script:
    """
    6_finalize_update.py \
        --leiden-tsv ${leiden_tsv} \
        --new-hashed-csv ${new_hashed_csv} \
        --existing-csv ${existing_csv} \
        --combined-fasta ${combined_fasta} \
        --output-prefix ${params.prefix_update}
    """
}
