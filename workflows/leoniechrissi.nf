/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { PICARD_MARKDUPLICATES     } from '../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_FAIDX            } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_SORT             } from '../modules/nf-core/samtools/sort/main'
include { TRIMGALORE                } from '../modules/nf-core/trimgalore/main'
include { HISAT2_ALIGN              } from '../modules/nf-core/hisat2/align/main'
include { HISAT2_EXTRACTSPLICESITES } from '../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD              } from '../modules/nf-core/hisat2/build/main'
include { SUBREAD_FEATURECOUNTS     } from '../modules/nf-core/subread/featurecounts/main'

include { FEATURECOUNTS_MERGE } from '../modules/local/merge/main'

include { paramsSummaryMap } from 'plugin/nf-schema'

include { paramsSummaryMultiqc             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore'
include { methodsDescriptionText           } from '../subworkflows/local/utils_nfcore_leoniechrissi_pipeline'
include { getGenomeAttribute               } from '../subworkflows/local/utils_nfcore_leoniechrissi_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta = getGenomeAttribute('fasta')
params.gtf = getGenomeAttribute('gtf')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LEONIECHRISSI {

    take:
    // channel: samplesheet read in from --input
    ch_samplesheet 

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // unique ids
    ch_samplesheet = ch_samplesheet.unique { meta -> meta.id }

    //
    // SUBWORKFLOW: fastqc, umnitools, trimgalore
    //

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
        ch_samplesheet,         // reads
        false,                  // skip_fastqc
        params.umi,             // with_umi
        false,                  // skip_umi_extract
        false,                  // skip_trimming
        2,                      // umi_discard_read
        3                       // min_trimmed_reads
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)

    //
    // MODULE: HISAT2_EXTRACTSPLICESITES
    // extracts splice sites
    //

    // preparare for extract splice sites
    path_gtf = getGenomeAttribute("gtf")
    ch_gtf = Channel
        .fromPath(path_gtf)
        .map { gtf_file -> tuple(id:gtf_file.getSimpleName(), gtf_file) }

    HISAT2_EXTRACTSPLICESITES(ch_gtf)
    ch_splicesites = HISAT2_EXTRACTSPLICESITES.out.txt
    ch_versions = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions.first())

    //
    // MODULE: HISAT2_BUILT
    // build index for reference genome
    //

    // prepare for hisat2 build
    path_fasta = getGenomeAttribute("fasta")
    ch_fasta = Channel
        .fromPath(path_fasta)
        .map { fasta_file -> tuple(id:fasta_file.getSimpleName(), fasta_file) }


    HISAT2_BUILD(
        ch_fasta,           // fasta
        ch_gtf,             // gtf
        ch_splicesites      // splicesites
    )
    ch_versions = ch_versions.mix(HISAT2_BUILD.out.versions.first())

    //
    // MODULE: HISAT2_ALIGN
    // align reads to reference genome
    //

    //prepare for hisat2_align
    ch_trimmed_reads = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
    ch_ht2_file = HISAT2_BUILD.out.index

    HISAT2_ALIGN(
        ch_trimmed_reads,           // reads
        ch_ht2_file.collect(),      // index 
        ch_splicesites.collect()    // splice sites
    )
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.collect{it[1]})

    // 
    // MODULE: samtools/faidx
    // get correct file format for the reference genome
    // 

    SAMTOOLS_FAIDX (
        ch_fasta,   // fasta
        [[],[]],    // fai
        false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    // 
    // MODULE: samtools sort
    // sort in preperation for Markduplicates
    //

    // prepare for samtools sort
    bam_files = HISAT2_ALIGN.out.bam
    bam_files_sorted = bam_files.map { meta, bam ->
        //meta_new = meta.clone()
        meta.id = "${meta.id}.sorted"
        tuple(meta, bam)
    }

    SAMTOOLS_SORT(
        bam_files_sorted,       // bam files
        ch_fasta.collect(),     // fasta
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    //  
    // MODULE: Markduplicates
    //
    
    // prepare for Markduplicates
    sort_bam = SAMTOOLS_SORT.out.bam
    bam_files_with_prefix = sort_bam.map { meta, bam ->
    meta.id = "${meta.id}.markdup"
    tuple(meta, bam)
    }
    fai_genome = SAMTOOLS_FAIDX.out.fai

    PICARD_MARKDUPLICATES(
        bam_files_with_prefix,  // reads
        ch_fasta.collect(),     // fasta
        fai_genome.collect()    // fai
    )
   
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collect{it[1]})

    //
    // MODULE: FeatureCounts
    // get read counts
    //

    // prepare for FeatureCounts
    ch_featurecount  = PICARD_MARKDUPLICATES.out.bam
    ch_featurecount_in = ch_featurecount.merge(ch_gtf.collect{it[1]})


    SUBREAD_FEATURECOUNTS(
        ch_featurecount_in
    )

    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]})


    //
    // MODULE: featurecounts merge
    // merge process
    //
    
    ch_files = SUBREAD_FEATURECOUNTS.out.counts.collect{it[1]}
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
