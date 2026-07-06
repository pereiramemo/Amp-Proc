// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1.2: comparative quality-check across all samples (R)
// Input:  all raw paired-end reads (staged together)
// Output: cross-sample QC plots (quality vs read count, histograms, PhiX) + stats
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_1_2_QUALITY_CHECK {

    container "ghcr.io/pereiramemo/amp-proc/1.2-quality-check:${params.container_tag}"
    publishDir "${params.output_dir}/",
           mode: "copy",
           enabled: params.full_output.toBoolean()

    tag "all samples"

    input:
    path reads

    output:
    path "1.2-quality-check-out"

    script:
    """
    1.2-quality-check.R \
        --input_dir      . \
        --output_dir     1.2-quality-check-out \
        --reads_pattern  "${params.reads_pattern}" \
        --nslots         ${task.cpus} \
        --overwrite      t
    """

}
