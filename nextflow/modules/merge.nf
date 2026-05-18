process EXTRACT_MERGE_CHUNK_REPS {
    /*
     * Collect Leiden cluster reps from all chunks into a single FASTA.
     * Wraps: bin/2_extract_and_merge_vclust_reps.py
     *
     * TODO: Implement this process.
     */

    input:
    path leiden_tsvs   // collected from all chunks
    path rep_fastas    // collected from all chunks

    output:
    path "vjdb1_vclust_reps_merged.fna.gz"

    script:
    """
    # TODO
    """
}

process MERGE_FINAL_CLUSTERS {
    /*
     * Map all sequences through the clustering hierarchy to final cluster IDs.
     * Wraps: bin/3_merge_clusters.py
     *
     * TODO: Implement this process.
     */

    input:
    path hashed_csv
    path merged_leiden_tsv
    path chunk_leiden_tsvs   // collected
    path chunk_linclust_tsvs // collected

    output:
    path "vjdb1_merged_reps.csv"

    script:
    """
    # TODO
    """
}

process EXTRACT_REPS_OF_REPS {
    /*
     * Extract final representative sequences from the merged FASTA.
     * Wraps: bin/4_extract_reps_of_reps.py
     *
     * TODO: Implement this process.
     */

    input:
    path merged_fasta
    path leiden_tsv

    output:
    path "vjdb1_merged_reps.fna.gz"

    script:
    """
    # TODO
    """
}
