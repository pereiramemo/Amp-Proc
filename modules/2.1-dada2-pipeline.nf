// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2.1: DADA2 ASV inference
// Input:  all trimmed paired-end reads (staged together)
// Output: ASV table + filtered reads + plots/tables
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_2_1_DADA2_PIPELINE {

    container "ghcr.io/pereiramemo/amp-proc/2.1-dada2-pipeline:${params.container_tag}"
    publishDir "${params.output_dir}/",
           mode: "copy",
           enabled: params.full_output.toBoolean()

    tag "ASVs"

    input:
    path all_trimmed

    output:
    path "2.1-dada2-pipeline-out",                        emit: dir
    path "2.1-dada2-pipeline-out/output/tables/asv_table.csv",   emit: asv_table

    script:
    """
    2.1-dada2-pipeline.R \
        --input_dir       . \
        --output_dir      2.1-dada2-pipeline-out \
        --pattern_r1      _R1_trimmed.fastq.gz \
        --pattern_r2      _R2_trimmed.fastq.gz \
        --trunc_r1        ${params.trunc_r1} \
        --trunc_r2        ${params.trunc_r2} \
        --min_overlap     ${params.dada2_min_overlap} \
        --bimeras_method  ${params.bimeras_method} \
        --nslots          ${params.nslots} \
        --no_qual_plot \
        --no_err_plot \
        --no_save_workspace \
        --overwrite
    """

}
