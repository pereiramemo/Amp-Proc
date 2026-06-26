// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2.2.2: VSEARCH OTU construction
// Input:  all per-sample 2.2.1 output directories (staged together)
// Output: pooled FASTA, OTU centroids, OTU table
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_2_2_2 {

    container "ghcr.io/epereira/amp-proc/module-2.2.2:latest"
    publishDir "${params.output_dir}/2.2.2-vsearch-pipeline",
           mode: "copy"

    input:
    path samples

    output:
    path "otu"

    script:
    """
    2.2.2-vsearch-pipeline.py \
        --samples_dir  . \
        --output_dir   otu \
        --nslots       ${params.nslots} \
        --identity     ${params.identity} \
        --overwrite    t
    """

}
