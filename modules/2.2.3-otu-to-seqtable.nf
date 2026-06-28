// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2.2.3: OTU -> sequence-keyed count table
// Input:  the "otu" directory produced by MODULE_2_2_2_VSEARCH_PIPELINE
// Output: DADA2-style sequence-keyed CSV (centroid sequence + per-sample counts)
//         so OTU centroids can be annotated by MODULE_3_TAXA_ANNOT (3-taxa-annot.R).
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_2_2_3_OTU_TO_SEQTABLE {

    container "ghcr.io/epereira/amp-proc/2.2.3-otu-to-seqtable:latest"
    publishDir "${params.output_dir}/2.2.3-otu-to-seqtable",
           mode: "copy",
           enabled: params.full_output.toBoolean()

    input:
    path otu_dir

    output:
    path "otu_seqtable.csv"

    script:
    """
    2.2.3-otu-to-seqtable.py \
        --otus_fasta  ${otu_dir}/output/otus/otus.fasta.gz \
        --otu_table   ${otu_dir}/output/otus/otu_table.tsv \
        --output      otu_seqtable.csv
    """

}
