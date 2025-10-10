process FEATURECOUNTS_MERGE {
    tag "$count_files"

    input:
    path count_files

    output:
    path "merged_counts.tsv", emit: counts

    script:
    """
    python3 ${projectDir}/modules/local/merge/merge_counts.py ${count_files.join(' ')} merged_counts.tsv
    """
    // python3 ${projectDir}/modules/local/merge/merge_counts.py ${count_files} merged_counts.tsv
}