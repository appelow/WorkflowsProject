<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-leoniechrissi_logo_dark.png">
    <img alt="nf-core/leoniechrissi" src="docs/images/nf-core-leoniechrissi_logo_light.png">
  </picture>
</h1>

<h1>
  <picture>
    <img alt="nf-core/leoniechrissi" src="docs/images/docs/images/nf-core-leoniechrissi.drawio.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/leoniechrissi/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/leoniechrissi/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/leoniechrissi/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/leoniechrissi/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/leoniechrissi/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/leoniechrissi)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23leoniechrissi-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/leoniechrissi)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/leoniechrissi** is a bioinformatics pipeline designed for the end-to-end processing and analysis of RNA sequencing (RNA-seq) data. It takes a samplesheet with FASTQ files as input and performs read quality control, adapter and quality trimming, splice-aware alignment to a reference genome, duplicate marking, and gene-level quantification.  
The pipeline outputs a merged count table and a comprehensive MultiQC report summarizing all processing and quality control steps.

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->

### Main workflow steps
1. **Read QC** for initial quality assessment ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. **Adapter trimming** ([`TrimGalore`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)) 
3. **Genome index building** (`HISAT2 build`)
4. **Splice site extraction** (`HISAT2 extract_splice_sites.py`)
5. **Alignment** with splice-aware mapping([`HISAT2`](https://daehwankimlab.github.io/hisat2/manual/))  
6. **BAM processing** for sorting and indexing ([`SAMtools`](https://www.htslib.org/)) 
7. **Duplicate marking** ([`Picard MarkDuplicates`](https://broadinstitute.github.io/picard/))
8. **Feature counting**for gene-level quantification  ([`Subread FeatureCounts`] (https://pubmed.ncbi.nlm.nih.gov/24227677/))
9. **Result aggregation** for summarizing all reports into one interactive HTML output ([`MultiQC`](https://multiqc.info/))


## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows: <br>
**samplesheet.csv**

```
sample,fastq_1,fastq_2,strandedness,group,replicate
WT_REP1,SRR6357070_1.fastq.gz,SRR6357070_2.fastq.gz,reverse,control,1
WT_REP2,SRR6357072_1.fastq.gz,SRR6357072_2.fastq.gz,reverse,starvation,2
RAP1_UNINDUCED_REP1,SRR6357073_1.fastq.gz,,reverse,,
RAP1_UNINDUCED_REP2,SRR6357074_1.fastq.gz,,reverse,,

```
Each row represents a fastq file (single-end) or a pair of fastq files (paired end). Each row needs a unique sample identifier. The strandedness refers to the library preparation and can be set to "forward", "reverse" or "unstranded". Group and replicate are optional.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

Now, you can run the pipeline using:

<!-- TODO: run only with test profile? -->

```bash
nextflow run ./main.nf \
   -profile test,docker 
```

## Pipeline output

<!-- TODO: merged_counts.tsv is the final table we need, correct?-->
All output files are stored in the results directory.
The final quality control summary, multiqc_report.html, is located in the results/multiqc folder.The merged read count table, merged_counts.tsv, can be found in the results/featurecounts folder.

## Credits

nf-core/leoniechrissi was originally written by Leonie Wehnert, Christina Parpoulas.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
