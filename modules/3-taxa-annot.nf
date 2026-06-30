// ─────────────────────────────────────────────────────────────────────────────
// MODULE 3: taxonomic annotation (SILVA via DADA2 / BLAST)
// Input:  tuple(label, sequence-keyed count table)
//           - label "asv": ASV table from MODULE_2_1_DADA2_PIPELINE
//           - label "otu": OTU centroid table from MODULE_2_2_3_OTU_TO_SEQTABLE
// Output: count table with taxonomic annotation
// Reference databases are mounted via containerOptions (see nextflow.config).
// Called once per branch via aliased imports (MODULE_3_TAXA_ANNOT_ASV / MODULE_3_TAXA_ANNOT_OTU).
// ─────────────────────────────────────────────────────────────────────────────

process MODULE_3_TAXA_ANNOT {

    container "ghcr.io/epereira/amp-proc/3-taxa-annot:latest"
    publishDir { "${params.output_dir}/3-taxa-annot/${label}" },
           mode: "copy"

    tag "${label.toUpperCase()}s"

    input:
    tuple val(table_delim), val(label), path(seq_table)

    output:
    path "3-taxa-annot-${label}-out",                                        emit: dir

    script:
    """
    3-taxa-annot.R \
        --input_asv_table   ${seq_table} \
        --table_delim       ${table_delim} \
        --output_dir        3-taxa-annot-${label}-out \
        --method            ${params.taxa_method} \
        --train_db          ${params.train_db} \
        --ref_db            ${params.ref_db} \
        --nslots            ${params.nslots} \
        --no_save_workspace \
        --overwrite
    """

}
