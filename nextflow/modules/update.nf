process PREPARE_UPDATE {
    /*
     * Hash-dedup new sequences and combine with existing reps.
     * Wraps: bin/5_prepare_update.py
     *
     * TODO: Implement this process.
     */

    input:
    path new_seqs
    path existing_reps

    output:
    path "combined_for_vclust.fna.gz", emit: combined_fasta
    path "new_hashed_reps.csv",        emit: hashed_csv

    script:
    """
    # TODO
    """
}

process FINALIZE_UPDATE {
    /*
     * Remap clusters after update, produce final CSVs and rep FASTA.
     * Wraps: bin/6_finalize_update.py
     *
     * TODO: Implement this process.
     */

    input:
    path leiden_tsv
    path new_hashed_csv
    path existing_csv
    path combined_fasta

    output:
    path "*_merged_reps.csv",    emit: merged_csv
    path "*_merged_reps.fna.gz", emit: rep_fasta

    script:
    """
    # TODO
    """
}
