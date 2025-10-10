/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_leoniechrissi_pipeline'
include { TRIMGALORE                } from '../modules/nf-core/trimgalore/main'
include { HISAT2_ALIGN              } from '../modules/nf-core/hisat2/align/main'
include { HISAT2_EXTRACTSPLICESITES } from '../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD              } from '../modules/nf-core/hisat2/build/main'
include { SUBREAD_FEATURECOUNTS     } from '../modules/nf-core/subread/featurecounts/main'

include { FEATURECOUNTS_MERGE       } from '../modules/local/merge/main'

include { getGenomeAttribute      } from '../subworkflows/local/utils_nfcore_leoniechrissi_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')
params.gtf = getGenomeAttribute('gtf')


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LEONIECHRISSI {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:


    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // unique ids
    ch_samplesheet = ch_samplesheet.unique { meta -> meta.id }


    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    //MODULE: TRIMGALORE
    //
    TRIMGALORE(
        ch_samplesheet
    )

    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]})
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())

    //
    // MODULE: HISAT2_EXTRACTSPLICESITES
    //

    path_fasta = getGenomeAttribute("fasta")
    ch_fasta = Channel
        .fromPath(path_fasta)
        .map { fasta_file -> tuple(id:fasta_file.getSimpleName(), fasta_file) }

    path_gtf = getGenomeAttribute("gtf")
    ch_gtf = Channel
        .fromPath(path_gtf)
        .map { gtf_file -> tuple(id:gtf_file.getSimpleName(), gtf_file) }

    HISAT2_EXTRACTSPLICESITES(ch_gtf)
    ch_splicesites = HISAT2_EXTRACTSPLICESITES.out.txt
    ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions.first())

    //
    // MODULE: HISAT2_BUILT
    //

    HISAT2_BUILD(
        fasta = ch_fasta,
        gtf = ch_gtf,
        splicesites = ch_splicesites
    )
    ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions.first())

    //
    // MODULE: HISAT2_ALIGN
    //

    ch_trimmed_reads = TRIMGALORE.out.reads
    ch_ht2_file = HISAT2_BUILD.out.index
    

    HISAT2_ALIGN(
        ch_trimmed_reads,
        ch_ht2_file.collect(),
        ch_splicesites.collect()
    )
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.collect{it[1]})

    //
    // MODULE: FeatureCounts
    //

    ch_featurecount  = HISAT2_ALIGN.out.bam
    ch_featurecount_in = ch_featurecount.merge(ch_gtf.collect{it[1]})


    SUBREAD_FEATURECOUNTS(
        ch_featurecount_in
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]})

    ch_files = SUBREAD_FEATURECOUNTS.out.counts.collect{it[1]}


    // merge process
    FEATURECOUNTS_MERGE(
        ch_files.collect()
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'leoniechrissi_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
