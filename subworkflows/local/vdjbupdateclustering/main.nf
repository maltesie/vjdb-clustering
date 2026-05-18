// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include {PREPARE_UPDATE}    from 

workflow VDJBUPDATECLUSTERING {

    take:
    // TODO nf-core: edit input (take) channels
    new_seqs // channel or file - is it always going to be one FASTA? If not, should be channels
    existing_reps
    existing_csv
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:
    // TODO nf-core: substitute modules here for the modules of your subworkflow

    // Step 1: hash-dedup new, combine with existing reps
        prepared = PREPARE_UPDATE(new_seqs, existing_reps)

        // Step 2: vclust on combined set
        vclust_out = VCLUST_CLUSTER(prepared.combined_fasta)

    // copied the take and emit of the VCLUST subworkflow for reference
    FASTA_VCLUST_PREFILTER_ALIGN_CLUSTER         
    take:
    ch_fasta       // channel: [ val(meta), [ fasta ] ]
    save_alignment // boolean
    metric         // string
    tani           // float
    gani           // float
    ani            // float
    emit:
    fasta     = ch_out.fasta           // channel: [ val(meta), [ fasta ] ]
    clusters  = ch_out.tsv             // channel: [ val(meta), [ tsv ] ]
    alignment = VCLUST_ALIGN.out.aln   // channel: [ val(meta), [ aln ] ]
    ids       = VCLUST_ALIGN.out.ids   // channel: [ val(meta), [ ids ] ]
    versions  = ch_versions            // channel: [ versions.yml ]
    // end VCLUST subworkflow reference
    

        // Step 3: remap and finalize
        finalized = FINALIZE_UPDATE(
            vclust_out.clusters,
            prepared.hashed_csv,
            existing_csv,
            prepared.combined_fasta
        )

    emit:
        csv   = finalized.merged_csv
        fasta = finalized.rep_fasta
}
