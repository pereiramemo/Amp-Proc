// ─────────────────────────────────────────────────────────────────────────────
// MODULE 1.3: primer removal (cutadapt)
// Input:  paired-end reads (R1, R2)
// Output: trimmed R1/R2 (emitted downstream) + per-sample logs/stats
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_1_3_PRIMERS_REMOVAL {

    container "ghcr.io/pereiramemo/amp-proc/1.3-primers-removal:${params.container_tag}"
    publishDir "${params.output_dir}/1.3-primers-removal-out",
           mode: "copy",
           enabled: params.full_output.toBoolean()

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name),
          path("${sample_name}/output/${sample_name}_R1_trimmed.fastq.gz"),
          path("${sample_name}/output/${sample_name}_R2_trimmed.fastq.gz"),
          path("${sample_name}/stats/1.3-primers-removal-${sample_name}-stats.tsv"),
          path("${sample_name}/logs/1.3-primers-removal-${sample_name}.log")

    script:
    """
    1.3-primers-removal.py \
        --reads1            ${reads[0]} \
        --reads2            ${reads[1]} \
        --sample_name       ${sample_name} \
        --output_dir        ${sample_name} \
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
