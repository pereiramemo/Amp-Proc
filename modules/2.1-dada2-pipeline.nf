// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2.1: DADA2 ASV inference
// Input:  all trimmed paired-end reads (staged together)
// Output: ASV table + filtered reads + plots/tables
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_2_1 {

    container "ghcr.io/epereira/amp-proc/module-2.1:latest"
    publishDir "${params.output_dir}/2.1-dada2-pipeline",
           mode: "copy",
           enabled: params.full_output.toBoolean()

    input:
    path reads

    output:
    path "dada2",                 emit: dir
    path "dada2/asv_table.csv",   emit: asv_table

    script:
    """
    2.1-dada2-pipeline.R \
        --input_dir       . \
        --output_dir      dada2 \
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
