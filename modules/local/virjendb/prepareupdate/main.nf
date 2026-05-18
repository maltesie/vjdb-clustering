process VIRJENDB_PREPAREUPDATE {
    label 'process_single'
    tag "${fna_gz}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'quay.io/biocontainers/biopython:1.81' }"

    input:
    path fna_gz
    path reps_fna_gz

    output:
    path "combined_for_vclust.fna.gz", emit: combined_reps_fna
    path "*.csv", emit: hashed_reps
    tuple val("${task.process}"), val('biopython'), eval('python -c "import Bio; print(Bio.__version__)"'), emit: versions_biopython, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    5_prepare_update.py \\
        --new-seqs $fna_gz \\
        --existing-reps $reps_fna_gz
    """

    stub:
    """
    echo | gzip > combined_for_vclust.fna.gz
    touch new_hashed_reps.csv
    """
}
