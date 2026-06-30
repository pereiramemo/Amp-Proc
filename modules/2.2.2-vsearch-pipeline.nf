// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2.2.2: VSEARCH OTU construction
// Input:  all per-sample 2.2.1 output directories (staged together)
// Output: pooled FASTA, OTU centroids, OTU table
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_2_2_2_VSEARCH_PIPELINE {

    container "ghcr.io/epereira/amp-proc/2.2.2-vsearch-pipeline:latest"
    publishDir "${params.output_dir}/",
           mode: "copy"
           
    tag "OTUs"       

    input:
    path samples

    output:               
    path "2.2.2-vsearch-pipeline-out",                       emit: dir
    path "2.2.2-vsearch-pipeline-out/output/otu_table.tsv",  emit: otu_table

    script:
    """
    2.2.2-vsearch-pipeline.py \
        --samples_dir  . \
        --output_dir   2.2.2-vsearch-pipeline-out \
        --nslots       ${params.nslots} \
        --identity     ${params.identity} \
        --overwrite    t
    """

}
