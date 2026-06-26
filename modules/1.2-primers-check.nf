// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1.2: IUPAC-aware primer-orientation check
// Input:  paired-end reads (R1, R2) + publish sub-directory label
// Output: per-sample primer-hit tables + summary
// Called twice (before / after primer removal) via aliased imports.
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_1_2 {

    container "ghcr.io/epereira/amp-proc/module-1.2:latest"
    publishDir { "${params.output_dir}/${publish_subdir}" },
           mode: "copy",
           enabled: params.full_output.toBoolean()

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)
    val publish_subdir

    output:
    path "${sample_name}"

    script:
    """
    1.2-check-primers.py \
        --reads1          ${reads[0]} \
        --reads2          ${reads[1]} \
        --output_dir      ${sample_name} \
        --primer_fwd      ${params.primer_fwd} \
        --primer_rev      ${params.primer_rev} \
        --subsample_size  ${params.subsample_size}
    """

}
