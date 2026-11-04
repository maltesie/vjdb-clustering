nextflow.enable.types = true

process LINCLUST {
    /*
     * Run MMseqs2 linclust on a single chunk.
     * This process is called once per chunk and can run in parallel.
     * Outputs are prefixed with the chunk name to avoid collisions when collected.
     */

    tag "${chunk_fasta.simpleName}"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe49c17754753d6cd9a31e5894117edaf1c81e3d6053a12bf6dc8f3af1dffe23/data'
        : 'community.wave.seqera.io/library/mmseqs2:18.8cc5c--af05c9a98d9f6139'}"
    conda "${moduleDir}/../../environment.yml"

    input:
    chunk_fasta: Path
    min_seq_id: Float
    mmseqs: String

    output:
    record(
        rep_fasta:   file("${chunk_fasta.simpleName}_linclust_rep_seq.fasta.gz"),
        cluster_tsv: file("${chunk_fasta.simpleName}_linclust_cluster.tsv")
    )

    script:
    def prefix = chunk_fasta.simpleName
    """
    mkdir -p tmpdir
    ${mmseqs} easy-linclust --threads ${task.cpus} --min-seq-id ${min_seq_id} \
        ${chunk_fasta} ${prefix}_linclust tmpdir
    gzip -2 ${prefix}_linclust_rep_seq.fasta
    rm -rf tmpdir ${prefix}_linclust_all_seqs.fasta
    """

    stub:
    def prefix = chunk_fasta.simpleName
    """
    echo "" | gzip > ${prefix}_linclust_rep_seq.fasta.gz
    touch ${prefix}_linclust_cluster.tsv
    """
}
