process HASH_DEDUP {
    /*
     * Hash-deduplicate input sequences and split into chunks.
     * Wraps: bin/1_prepare_clustering.py
     */

    input:
    path fasta
    val  block_size

    output:
    path "chunks/vjdb1_*.fna.gz", emit: fasta_chunks
    path "chunks/vjdb1_hashed_reps.csv", emit: hashed_csv

    script:
    """
    1_prepare_clustering.py --input ${fasta} --outdir chunks --block-size ${block_size}
    """
}
