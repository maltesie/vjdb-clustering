/*
 * Checks against the expected results on the bundled test data.
 * Only used by -profile test and -profile test_initial (params.check_test).
 */

nextflow.enable.types = true

process CHECK_TEST_INITIAL {
    /*
     * Wraps: resources/usr/bin/check_test_initial.py
     */

    debug true
    container 'community.wave.seqera.io/library/biopython_pandas:a2cf2d914972cc52'
    conda "${moduleDir}/../../environment.yml"

    input:
    merged_csv: Path
    rep_fasta: Path
    expected_csv: Path

    stage:
    stageAs expected_csv, 'expected/vjdb1_merged_reps.csv'

    script:
    """
    check_test_initial.py ${merged_csv} ${rep_fasta} expected/vjdb1_merged_reps.csv
    """
}

process CHECK_TEST_UPDATE {
    /*
     * Wraps: resources/usr/bin/check_test_update.py
     */

    debug true
    container 'community.wave.seqera.io/library/biopython_pandas:a2cf2d914972cc52'
    conda "${moduleDir}/../../environment.yml"

    input:
    merged_csv: Path
    new_clusters: Path
    events: Path
    existing_csv: Path
    prefix: String

    script:
    """
    check_test_update.py ${prefix} ${existing_csv}
    """
}
