nextflow.enable.types = true

process EXTRACT_MERGE_CHUNK_REPS {
    /*
     * Collect Leiden cluster reps from all chunks into a single FASTA.
     * Wraps: resources/usr/bin/2_extract_and_merge_vclust_reps.py
     */

    container 'community.wave.seqera.io/library/biopython_pandas:a2cf2d914972cc52'
    conda "${moduleDir}/../../environment.yml"

    input:
    leiden_tsvs: Bag<Path>   // collected from all chunks
    rep_fastas: Bag<Path>    // collected from all chunks

    output:
    file("vjdb1_vclust_reps_merged.fna.gz")

    script:
    """
    2_extract_and_merge_vclust_reps.py \
        --leiden-tsvs ${leiden_tsvs.join(' ')} \
        --rep-fastas ${rep_fastas.join(' ')} \
        -o vjdb1_vclust_reps_merged.fna.gz
    """

    stub:
    """
    echo "" | gzip > vjdb1_vclust_reps_merged.fna.gz
    """
}

process MERGE_FINAL_CLUSTERS {
    /*
     * Map all sequences through the clustering hierarchy to final cluster IDs.
     * Wraps: resources/usr/bin/3_merge_clusters.py
     */

    container 'community.wave.seqera.io/library/biopython_pandas:a2cf2d914972cc52'
    conda "${moduleDir}/../../environment.yml"

    input:
    hashed_csv: Path
    merged_leiden_tsv: Path
    chunk_leiden_tsvs: Bag<Path>    // collected
    chunk_linclust_tsvs: Bag<Path>  // collected

    output:
    file("vjdb1_merged_reps.csv")

    script:
    """
    3_merge_clusters.py \
        --hashed-csv ${hashed_csv} \
        --merged-leiden ${merged_leiden_tsv} \
        --chunk-leidens ${chunk_leiden_tsvs.join(' ')} \
        --chunk-linclusts ${chunk_linclust_tsvs.join(' ')} \
        -o vjdb1_merged_reps.csv
    """

    stub:
    """
    echo "vjdb_id,copy_of,seq_gc,seq_len,seq_hash,cluster_rep" > vjdb1_merged_reps.csv
    """
}

process EXTRACT_REPS_OF_REPS {
    /*
     * Extract final representative sequences from the merged FASTA.
     * Wraps: resources/usr/bin/4_extract_reps_of_reps.py
     */

    container 'community.wave.seqera.io/library/biopython_pandas:a2cf2d914972cc52'
    conda "${moduleDir}/../../environment.yml"

    input:
    merged_fasta: Path
    leiden_tsv: Path

    output:
    file("vjdb1_merged_reps.fna.gz")

    script:
    """
    4_extract_reps_of_reps.py \
        --input-fasta ${merged_fasta} \
        --leiden-tsv ${leiden_tsv} \
        -o vjdb1_merged_reps.fna.gz
    """

    stub:
    """
    echo "" | gzip > vjdb1_merged_reps.fna.gz
    """
}
