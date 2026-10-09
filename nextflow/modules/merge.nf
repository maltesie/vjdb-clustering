process EXTRACT_MERGE_CHUNK_REPS {
    /*
     * Collect Leiden cluster reps from all chunks into a single FASTA.
     * Wraps: bin/2_extract_and_merge_vclust_reps.py
     */
     
    input:
    path leiden_tsvs   // collected from all chunks
    path rep_fastas    // collected from all chunks

    output:
    path "vjdb1_vclust_reps_merged.fna.gz"

    script:
    """
    2_extract_and_merge_vclust_reps.py \
        --leiden-tsvs ${leiden_tsvs} \
        --rep-fastas ${rep_fastas} \
        -o vjdb1_vclust_reps_merged.fna.gz
    """
}

process MERGE_FINAL_CLUSTERS {
    /*
     * Map all sequences through the clustering hierarchy to final cluster IDs.
     * Wraps: bin/3_merge_clusters.py
     */
    
    publishDir params.outdir, mode: 'move'
    
    input:
    path hashed_csv
    path merged_leiden_tsv
    path chunk_leiden_tsvs   // collected
    path chunk_linclust_tsvs // collected

    output:
    path "vjdb1_merged_reps.csv"

    script:
    """
    3_merge_clusters.py \
        --hashed-csv ${hashed_csv} \
        --merged-leiden ${merged_leiden_tsv} \
        --chunk-leidens ${chunk_leiden_tsvs} \
        --chunk-linclusts ${chunk_linclust_tsvs} \
        -o vjdb1_merged_reps.csv
    """
}

process EXTRACT_REPS_OF_REPS {
    /*
     * Extract final representative sequences from the merged FASTA.
     * Wraps: bin/4_extract_reps_of_reps.py
     */

    publishDir params.outdir, mode: 'move'
    
    input:
    path merged_fasta
    path leiden_tsv

    output:
    path "vjdb1_merged_reps.fna.gz"

    script:
    """
    4_extract_reps_of_reps.py \
        --input-fasta ${merged_fasta} \
        --leiden-tsv ${leiden_tsv} \
        -o vjdb1_merged_reps.fna.gz
    """
}
