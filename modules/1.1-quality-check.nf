// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1.1: quality-check report (fastp)
// Input:  paired-end reads (R1, R2)
// Output: per-sample fastp HTML/JSON reports + summary
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_1_1_QUALITY_CHECK {

    container "ghcr.io/pereiramemo/amp-proc/1.1-quality-check:latest"
    publishDir "${params.output_dir}/1.1-quality-check-out",
           mode: "copy",
           enabled: params.full_output.toBoolean()

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)

    output:
    path "${sample_name}"

    script:
    """
    1.1-quality-check.py \
        --reads1                     ${reads[0]} \
        --reads2                     ${reads[1]} \
        --sample_name                 ${sample_name} \
        --output_dir                 ${sample_name} \
        --nslots                     ${params.nslots} \
        --overwrite                  t
    """

}
