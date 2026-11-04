nextflow.enable.types = true

process VCLUST_CLUSTER {
    /*
     * Run vclust prefilter → align → cluster (Leiden) on a FASTA.
     * Used both per-chunk and on the merged representatives.
     *
     * All three steps use consistent ANI/qcov thresholds passed as inputs.
     * Output is renamed with a prefix to avoid collisions when collected.
     */

    tag "${rep_fasta.simpleName}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vclust:1.3.1--py313h9ee0642_0':
        'biocontainers/vclust:1.3.1--py313h9ee0642_0' }"

    conda "${moduleDir}/../../environment.yml"

    input:
    rep_fasta: Path
    ani: Float
    qcov: Float
    vclust: String

    output:
    file("${rep_fasta.simpleName}_leiden.tsv")

    script:
    def prefix = rep_fasta.simpleName
    """
    ${vclust} prefilter \
        -i ${rep_fasta} \
        -o fltr.txt \
        -t ${task.cpus} \
        --min-ident ${ani} \
        --batch-size 100000 \
        --kmers-fraction 0.2 \
        --max-seqs 2000

    ${vclust} align \
        -i ${rep_fasta} \
        -o ani.tsv \
        -t ${task.cpus} \
        --filter fltr.txt \
        --filter-threshold ${ani} \
        --out-ani ${ani} \
        --out-qcov ${qcov}

    ${vclust} cluster \
        -i ani.tsv \
        -o clusterreps_leiden.tsv \
        --ids ani.ids.tsv \
        --algorithm leiden \
        --metric ani \
        --ani ${ani} \
        --qcov ${qcov} \
        --out-repr

    mv clusterreps_leiden.tsv ${prefix}_leiden.tsv
    rm -f fltr.txt ani.tsv ani.ids.tsv
    """

    stub:
    """
    touch ${rep_fasta.simpleName}_leiden.tsv
    """
}
