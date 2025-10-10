process FEATURECOUNTS_MERGE {
    tag "$input_folder"

    input:
    path input_folder

    output:
    path "merged_counts.tsv", emit: counts

    script:
    """
    python3 ${projectDir}/modules/local/merge/merge_counts.py ${input_folder} merged_counts.tsv
    """
}