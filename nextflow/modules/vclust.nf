process VCLUST_CLUSTER {
    /*
     * Run vclust prefilter → align → cluster (Leiden) on a FASTA.
     * Used both per-chunk and on the merged representatives.
     *
     * All three steps use consistent ANI/qcov thresholds passed as parameters.
     */

    tag "${rep_fasta.simpleName}"

    input:
    path rep_fasta
    val  ani
    val  qcov

    output:
    path "clusterreps_leiden.tsv", emit: leiden_tsv

    script:
    """
    ${params.vclust} prefilter \
        -i ${rep_fasta} \
        -o fltr.txt \
        -t ${task.cpus} \
        --min-ident ${ani} \
        --batch-size 100000 \
        --kmers-fraction 0.2 \
        --max-seqs 2000

    ${params.vclust} align \
        -i ${rep_fasta} \
        -o ani.tsv \
        -t ${task.cpus} \
        --filter fltr.txt \
        --filter-threshold ${ani} \
        --out-ani ${ani} \
        --out-qcov ${qcov}

    ${params.vclust} cluster \
        -i ani.tsv \
        -o clusterreps_leiden.tsv \
        --ids ani.ids.tsv \
        --algorithm leiden \
        --metric ani \
        --ani ${ani} \
        --qcov ${qcov} \
        --out-repr

    rm -f fltr.txt ani.tsv ani.ids.tsv
    """
}
