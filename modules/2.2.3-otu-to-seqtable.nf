// ─────────────────────────────────────────────────────────────────────────────
// MODULE 2.2.3: OTU -> sequence-keyed count table
// Input:  the "otu" directory produced by MODULE_2_2_2
// Output: DADA2-style sequence-keyed CSV (centroid sequence + per-sample counts)
//         so OTU centroids can be annotated by MODULE_3 (3-taxa_annot.R).
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_2_2_3 {

    container "ghcr.io/epereira/amp-proc/module-2.2.3:latest"
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
        --otus_fasta  ${otu_dir}/otus/otus.fasta.gz \
        --otu_table   ${otu_dir}/otus/otu_table.tsv \
        --output      otu_seqtable.csv
    """

}
