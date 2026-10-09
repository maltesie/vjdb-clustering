/*
 * Checks against the expected results on the bundled test data.
 * Only used by -profile test and -profile test_initial (params.check_test).
 */

process CHECK_TEST_INITIAL {
    /*
     * Wraps: bin/check_test_initial.py
     */

    debug true

    input:
    path merged_csv
    path rep_fasta
    path expected_csv, stageAs: 'expected/*'

    script:
    """
    check_test_initial.py ${merged_csv} ${rep_fasta} ${expected_csv}
    """
}

process CHECK_TEST_UPDATE {
    /*
     * Wraps: bin/check_test_update.py
     */

    debug true

    input:
    path merged_csv
    path new_clusters
    path events
    path existing_csv

    script:
    """
    check_test_update.py ${params.prefix_update} ${existing_csv}
    """
}
