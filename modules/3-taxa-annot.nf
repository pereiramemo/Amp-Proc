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

    tag "${label}"

    input:
    tuple val(label), path(seq_table)

    output:
    path "${label}_tax_annot.csv"

    script:
    """
    # toolbox.R is sourced from the script's directory; stage it next to the run
    cp \$(dirname \$(command -v 3-taxa-annot.R))/toolbox.R .

    3-taxa-annot.R \
        --input_asv_table   ${seq_table} \
        --output_asv_table  ${label}_tax_annot.csv \
        --method            ${params.taxa_method} \
        --evalue            ${params.evalue} \
        --min_identity      ${params.min_identity} \
        --train_db          ${params.train_db} \
        --ref_db            ${params.ref_db} \
        --nslots            ${params.nslots} \
        --no_save_workspace \
        --overwrite
    """

}
