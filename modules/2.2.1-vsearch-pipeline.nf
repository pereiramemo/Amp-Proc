// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2.2.1: VSEARCH per-sample processing
// Input:  trimmed paired-end reads (R1, R2)
// Output: per-sample dir (merged -> filtered -> derep -> chimera-checked)
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_2_2_1_VSEARCH_PIPELINE {

    container "ghcr.io/epereira/amp-proc/2.2.1-vsearch-pipeline:latest"
    publishDir "${params.output_dir}/2.2.1-vsearch-pipeline-out",
           mode: "copy",
           enabled: params.full_output.toBoolean()

    tag "${sample_name}"

    input:
    tuple val(sample_name), path(reads)

    output:
    path "${sample_name}"

    script:
    """
    2.2.1-vsearch-pipeline.py \
        --reads1          ${reads[0]} \
        --reads2          ${reads[1]} \
        --sample_name     ${sample_name} \
        --output_dir      ${sample_name} \
        --nslots          ${params.nslots} \
        --min_length      ${params.vsearch_min_length} \
        --fastq_minovlen  ${params.fastq_minovlen} \
        --fastq_maxdiffs  ${params.fastq_maxdiffs} \
        --fastq_maxee     ${params.fastq_maxee} \
        --min_size        ${params.min_size} \
        --abskew          ${params.abskew} \
        --overwrite       t
    """

}
