process VCLUST_CLUSTER {
    /*
     * Run vclust prefilter → align → cluster (Leiden) on a FASTA.
     * Used both per-chunk and on the merged representatives.
     *
     * All three steps use consistent ANI/qcov thresholds passed as parameters.
     * Output is renamed with a prefix to avoid collisions when collected.
     *
     * If vclust align produces no pairs, falls back to treating each sequence
     * as its own cluster (vclust cluster crashes on empty ANI input).
     */

    tag "${rep_fasta.simpleName}"

    input:
    path rep_fasta
    val  ani
    val  qcov

    output:
    path "${rep_fasta.simpleName}_leiden.tsv", emit: leiden_tsv

    script:
    def prefix = rep_fasta.simpleName
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

    ANI_LINES=\$(wc -l < ani.tsv)
    if [ "\${ANI_LINES}" -le 1 ]; then
        echo "No ANI pairs above threshold — each sequence is its own cluster"
        echo -e "object\\tcluster" > ${prefix}_leiden.tsv
        zgrep '^>' ${rep_fasta} | sed 's/^>//' | cut -d' ' -f1 | \
            awk '{print \$1 "\\t" \$1}' >> ${prefix}_leiden.tsv
    else
        ${params.vclust} cluster \
            -i ani.tsv \
            -o clusterreps_leiden.tsv \
            --ids ani.ids.tsv \
            --algorithm leiden \
            --metric ani \
            --ani ${ani} \
            --qcov ${qcov} \
            --out-repr
        mv clusterreps_leiden.tsv ${prefix}_leiden.tsv
    fi

    rm -f fltr.txt ani.tsv ani.ids.tsv
    """
}
