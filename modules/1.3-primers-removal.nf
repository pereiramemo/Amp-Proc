// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1.3: primer removal (cutadapt)
// Input:  paired-end reads (R1, R2)
// Output: trimmed R1/R2 (emitted downstream) + per-sample logs/stats
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_1_3 {

    container "ghcr.io/epereira/amp-proc/module-1.3:latest"
    publishDir "${params.output_dir}/1.3-remove-primers",
           mode: "copy",
           enabled: params.full_output.toBoolean()

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name),
          path("${sample_name}/trimmed/*_R1_trimmed.fastq.gz"),
          path("${sample_name}/trimmed/*_R2_trimmed.fastq.gz")

    script:
    """
    1.3-remove-primers.py \
        --reads1            ${reads[0]} \
        --reads2            ${reads[1]} \
        --output_dir        ${sample_name} \
        --trimmed_dir       ${sample_name}/trimmed \
        --primer_fwd        ${params.primer_fwd} \
        --primer_rev        ${params.primer_rev} \
        --nslots            ${params.nslots} \
        --error_rate        ${params.error_rate} \
        --min_overlap       ${params.min_overlap} \
        --min_length        ${params.min_length} \
        --discard_untrimmed ${params.discard_untrimmed} \
        --compress          t \
        --overwrite         t
    """

}
