// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2.2.2: VSEARCH OTU construction
// Input:  all per-sample 2.2.1 output directories (staged together)
// Output: pooled FASTA, OTU centroids, OTU table
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_2_2_2_VSEARCH_PIPELINE {

    container "ghcr.io/pereiramemo/amp-proc/2.2.2-vsearch-pipeline:${params.container_tag}"
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
        --nslots       ${task.cpus} \
        --identity     ${params.identity} \
        --overwrite    t
    """

}
