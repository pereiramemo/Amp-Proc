#!/usr/bin/env -S nextflow run

// Include modules
include { MODULE_1_1_QUALITY_CHECK }                                from './modules/1.1-quality-check.nf'
include { MODULE_1_2_PRIMERS_CHECK as MODULE_1_2_PRIMERS_CHECK_BEFORE }           from './modules/1.2-primers-check.nf'
include { MODULE_1_2_PRIMERS_CHECK as MODULE_1_2_PRIMERS_CHECK_AFTER }            from './modules/1.2-primers-check.nf'
include { MODULE_1_3_PRIMERS_REMOVAL }                                from './modules/1.3-primers-removal.nf'
include { MODULE_2_1_DADA2_PIPELINE }                                from './modules/2.1-dada2-pipeline.nf'
include { MODULE_2_2_1_VSEARCH_PIPELINE }                              from './modules/2.2.1-vsearch-pipeline.nf'
include { MODULE_2_2_2_VSEARCH_PIPELINE }                              from './modules/2.2.2-vsearch-pipeline.nf'
include { MODULE_2_2_3_OTU_TO_SEQTABLE }                              from './modules/2.2.3-otu-to-seqtable.nf'
include { MODULE_3_TAXA_ANNOT as MODULE_3_TAXA_ANNOT_ASV }                  from './modules/3-taxa-annot.nf'
include { MODULE_3_TAXA_ANNOT as MODULE_3_TAXA_ANNOT_OTU }                  from './modules/3-taxa-annot.nf'

workflow {

    main:
    if (params.help) {
        log.info """
        Amp-Proc: amplicon processing from paired-end reads to ASV/OTU tables

        Usage: nextflow run main.nf [options]

        General:
          --input_dir       DIR   Input directory with paired-end FASTQ files (default: ${params.input_dir})
          --reads_pattern   STR   Glob pattern for fromFilePairs (default: ${params.reads_pattern})
          --output_dir      DIR   Output directory (default: ${params.output_dir})
          --nslots          INT   CPU threads per tool (default: ${params.nslots})
          --method          STR   Denoising branch: dada2 | vsearch | both (default: ${params.method})
          --full_output     BOOL  Publish all intermediate outputs (default: ${params.full_output})
          --skip_tax_annot  BOOL  Skip MODULE_3_TAXA_ANNOT taxonomic annotation (default: ${params.skip_tax_annot})
          --maxForks        INT   Max parallel process instances (default: ${params.maxForks})

        Primers (MODULE_1_2_PRIMERS_CHECK, MODULE_1_3_PRIMERS_REMOVAL):
          --primer_fwd      STR   Forward primer 5'->3' (default: ${params.primer_fwd})
          --primer_rev      STR   Reverse primer 5'->3' (default: ${params.primer_rev})

        MODULE_1_2_PRIMERS_CHECK — primer check:
          --subsample_size  INT   Reads to subsample per file (default: ${params.subsample_size})

        MODULE_1_3_PRIMERS_REMOVAL — cutadapt primer removal:
          --error_rate        NUM  Max allowed error rate (default: ${params.error_rate})
          --min_overlap       INT  Min primer-read overlap (default: ${params.min_overlap})
          --min_length        INT  Discard reads shorter than this (default: ${params.min_length})
          --discard_untrimmed STR  Discard reads with no primer, t/f (default: ${params.discard_untrimmed})

        MODULE_2_1_DADA2_PIPELINE — DADA2 (ASV):
          --trunc_r1          INT  Truncate R1 from 3' end (default: ${params.trunc_r1})
          --trunc_r2          INT  Truncate R2 from 3' end (default: ${params.trunc_r2})
          --dada2_min_overlap INT  Min overlap when merging (default: ${params.dada2_min_overlap})
          --bimeras_method    STR  pooled | consensus | per-sample (default: ${params.bimeras_method})

        MODULE_2_2_1_VSEARCH_PIPELINE — VSEARCH per-sample:
          --fastq_minovlen     INT  Min overlap for PE merging (default: ${params.fastq_minovlen})
          --fastq_maxdiffs     INT  Max mismatches in overlap (default: ${params.fastq_maxdiffs})
          --fastq_maxee        NUM  Max expected errors per read (default: ${params.fastq_maxee})
          --min_size           INT  Min abundance after derep (default: ${params.min_size})
          --abskew             NUM  Min parent/child ratio, chimeras (default: ${params.abskew})
          --vsearch_min_length INT  Min merged-read length (default: ${params.vsearch_min_length})

        MODULE_2_2_2_VSEARCH_PIPELINE — VSEARCH OTU construction:
          --identity        NUM   OTU clustering identity 0-1 (default: ${params.identity})

        MODULE_3_TAXA_ANNOT — taxonomic annotation (SILVA):
          --taxa_method     STR   NBC | NBCandEM | BLAST (default: ${params.taxa_method})
          --evalue          NUM   BLAST e-value (default: ${params.evalue})
          --min_identity    NUM   BLAST min identity (default: ${params.min_identity})
          --train_db        PATH  NBC training database (default: ${params.train_db})
          --ref_db          PATH  EM reference database (default: ${params.ref_db})
        """.stripIndent()
        exit 0
    }

    // Read paired-end samples from input_dir using the pattern defined in nextflow.config
    reads_ch = channel.fromFilePairs(
        "${params.input_dir}/${params.reads_pattern}",
        checkIfExists: true
    )

    method = params.method.toString().toLowerCase()
    log.info "Denoising branch(es): ${method}"
    if ("${params.skip_tax_annot}" == "true") {
        log.info "Taxonomic annotation will be skipped"
    }

    // MODULE_1_1_QUALITY_CHECK: fastp quality-check report (diagnostic, always runs)
    MODULE_1_1_QUALITY_CHECK(reads_ch)

    // MODULE_1_2_PRIMERS_CHECK: primer check before trimming (diagnostic, always runs)
    MODULE_1_2_PRIMERS_CHECK_BEFORE(reads_ch, "1.2-primers-check-before-out")

    // MODULE_1_3_PRIMERS_REMOVAL: primer removal with cutadapt
    trimmed_out = MODULE_1_3_PRIMERS_REMOVAL(reads_ch)
    trimmed_reads = trimmed_out.map { sample_name, r1, r2, _stats, _log -> tuple(sample_name, [r1, r2]) }

    // MODULE_1_2_PRIMERS_CHECK: primer check after trimming (diagnostic)
    MODULE_1_2_PRIMERS_CHECK_AFTER(trimmed_reads, "1.2-primers-check-after-out")

    // DADA2 (ASV) branch
    if (method in ['dada2', 'both']) {
        all_trimmed = trimmed_out.flatMap { _sample_name, r1, r2, _stats, _log -> [r1, r2] }.collect()
        dada2_out   = MODULE_2_1_DADA2_PIPELINE(all_trimmed)
    }    

  // MODULE_3_TAXA_ANNOT: taxonomic annotation of the ASV table
  //      if (!params.skip_tax_annot.toBoolean()) {
  //          MODULE_3_TAXA_ANNOT_ASV(dada2_out.asv_table.map { tbl -> tuple("asv", tbl) })
  //      }
  //  }

    // VSEARCH (OTU) branch
    if (method in ['vsearch', 'both']) {
        vsearch_sample = MODULE_2_2_1_VSEARCH_PIPELINE(trimmed_reads)
        all_samples    = vsearch_sample.collect()
        otu_out        = MODULE_2_2_2_VSEARCH_PIPELINE(all_samples)
    }    

        // MODULE_3_TAXA_ANNOT: taxonomic annotation of the OTU centroids
  //      if (!params.skip_tax_annot.toBoolean()) {
  //          otu_seqtable = MODULE_2_2_3_OTU_TO_SEQTABLE(otu_out)
  //          MODULE_3_TAXA_ANNOT_OTU(otu_seqtable.map { tbl -> tuple("otu", tbl) })
  //      }
 //   }

}