// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1.1: quality-check report (fastp)
// Input:  paired-end reads (R1, R2)
// Output: per-sample fastp HTML/JSON reports + summary
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_1_1 {

    container "ghcr.io/epereira/amp-proc/module-1.1:latest"
    publishDir "${params.output_dir}/1.1-quality-check",
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
        --output_dir                 ${sample_name} \
        --nslots                     ${params.nslots} \
        --min_length                 ${params.qc_min_length} \
        --qualified_quality_phred    ${params.qualified_quality_phred} \
        --unqualified_percent_limit  ${params.unqualified_percent_limit} \
        --disable_adapter_trimming   ${params.disable_adapter_trimming} \
        --overwrite                  t
    """

}
