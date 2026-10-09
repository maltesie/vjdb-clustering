process LINCLUST {
    /*
     * Run MMseqs2 linclust on a single chunk.
     * This process is called once per chunk and can run in parallel.
     * Outputs are prefixed with the chunk name to avoid collisions when collected.
     */

    tag "${chunk_fasta.simpleName}"

    input:
    path chunk_fasta

    output:
    path "${chunk_fasta.simpleName}_linclust_rep_seq.fasta.gz", emit: rep_fasta
    path "${chunk_fasta.simpleName}_linclust_cluster.tsv",      emit: cluster_tsv

    script:
    def prefix = chunk_fasta.simpleName
    """
    mkdir -p tmpdir
    ${params.mmseqs} easy-linclust --threads ${task.cpus} --min-seq-id 0.95 \
        ${chunk_fasta} ${prefix}_linclust tmpdir
    gzip -2 ${prefix}_linclust_rep_seq.fasta
    rm -rf tmpdir ${prefix}_linclust_all_seqs.fasta
    """
}
