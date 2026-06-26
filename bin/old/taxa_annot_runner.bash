#!/bin/bash

###############################################################################
### 1. Load env
###############################################################################

set -o pipefail 

source "/home/epereira/workspace/repositories/amplicon_pipelines/scripts_v2/conf.sh"

###############################################################################
# 2. Define help
###############################################################################

show_usage(){
  cat <<EOF
Usage: ./dada2_pipelin_runnar.bash <options>
--help                          print this help
--input_asv_table CHAR          asv table generated with DADA2 (required)
--input_fasta CHAR              fasta file with sequences to be annotated (if not given it will obtain from the ASVs)
--output_asv_table CHAR         output asv table with taxonomic annotation (required)
--method CHAR                   method used to annotate the amplicon sequences. One of NBC, NBCandEM, BLAST. NBC: Naive Bayes Classifier; EM: Exact Matching. (Default NBC).
--evalue NUM                    evalue used in BLAST search (default 1e-10)
--min_identity NUM              minimum identity used in BLAST search (default 97)
--train_db CHAR                 training database to run NBC (default silva_nr_v138_train_seq.fa.gz)
--ref_db CHAR                   reference database to run EM (default silva_species_assignment_v138.fa.gz)
--blast_db CHAR                 blast formatted database to run BLAST (default  SILVA_138_SSURef_NR99_tax_silva.fasta)
--taxa_map CHAR                 tsv file mapping silva acc with taxonomy (used when running BLAST; default taxmap_slv_ssu_ref_nr_138.txt)
--blout CHAR                    tabular blast output file
--nslots NUM                    number of threads used (default 12)
--save_workspace T|F            save R workspace image (default T)
--overwrite T|F                 overwrite previous output (default F)
EOF
}

###############################################################################
# 3. Parse input parameters
###############################################################################

while :; do
  case "${1}" in
    --help) # Call a "show_help" function to display a synopsis, then exit.
    show_usage
    exit 1;
    ;;
#############
  --input_asv_table)
  if [[ -n "${2}" ]]; then
    INPUT_ASV_TABLE="${2}"
    shift
  fi
  ;;
  --input_asv_table=?*)
  INPUT_ASV_TABLE="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --input_asv_table=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --input_fasta)
  if [[ -n "${2}" ]]; then
    INPUT_FASTA="${2}"
    shift
  fi
  ;;
  --input_fasta=?*)
  INPUT_FASTA="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --input_fasta=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;  
#############
  --output_asv_table)
  if [[ -n "${2}" ]]; then
    OUTPUT_ASV_TABLE="${2}"
    shift
  fi
  ;;
  --output_asv_table=?*)
  OUTPUT_ASV_TABLE="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --output_asv_table=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --method)
  if [[ -n "${2}" ]]; then
    METHOD="${2}"
    shift
  fi
  ;;
  --method=?*)
  METHOD="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --method=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --evalue)
  if [[ -n "${2}" ]]; then
    EVALUE="${2}"
    shift
  fi
  ;;
  --evalue=?*)
  EVALUE="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --evalue=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --min_identity)
  if [[ -n "${2}" ]]; then
    MIN_IDENTITY="${2}"
    shift
  fi
  ;;
  --min_identity=?*)
  MIN_IDENTITY="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --min_identity=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;; 
#############
  --train_db)
  if [[ -n "${2}" ]]; then
    TRAIN_DB="${2}"
    shift
  fi
  ;;
  --train_db=?*)
  TRAIN_DB="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --train_db=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;  
#############
  --ref_db)
  if [[ -n "${2}" ]]; then
    REF_DB="${2}"
    shift
  fi
  ;;
  --ref_db=?*)
  REF_DB="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --ref_db=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --blast_db)
  if [[ -n "${2}" ]]; then
    BLAST_DB="${2}"
    shift
  fi
  ;;
  --blast_db=?*)
  BLAST_DB="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --blast_db=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --taxa_map)
  if [[ -n "${2}" ]]; then
    TAXA_MAP="${2}"
    shift
  fi
  ;;
  --taxa_map=?*)
  TAXA_MAP="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --taxa_map=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --blout)
  if [[ -n "${2}" ]]; then
    BLOUT="${2}"
    shift
  fi
  ;;
  --blout=?*)
  BLOUT="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --blout=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --nslots)
  if [[ -n "${2}" ]]; then
    NSLOTS="${2}"
    shift
  fi
  ;;
  --nslots=?*)
  NSLOTS="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --nslots=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --save_workspace)
  if [[ -n "${2}" ]]; then
    SAVE_WORKSPACE="${2}"
    shift
  fi
  ;;
  --save_workspace=?*)
  SAVE_WORKSPACE="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --save_workspace=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --overwrite)
  if [[ -n "${2}" ]]; then
    OVERWRITE="${2}"
    shift
  fi
  ;;
  --overwrite=?*)
  OVERWRITE="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --overwrite=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;; 
############ End of all options.
  --)       
  shift
  break
  ;;
  -?*)
  printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
  ;;
  *) # Default case: If no more options, then break out of the loop.
  break
  esac
  shift
done  

###############################################################################
# 4. Define defaults
###############################################################################

if [[ -z "${METHOD}" ]]; then
  METHOD="NBC"
fi

if [[ -z "${EVALUE}" ]]; then
  EVALUE="1e-10"
fi

if [[ -z "${MIN_IDENTITY}" ]]; then
  MIN_IDENTITY="97"
fi

if [[ -z "${TRAIN_DB}" ]]; then
  TRAIN_DB="${SILVA_DIR}/silva_nr_v138_train_seq.fa.gz"
fi

if [[ -z "${REF_DB}" ]]; then
  REF_DB="${SILVA_DIR}/silva_species_assignment_v138.fa.gz"
fi

if [[ -z "${BLAST_DB}" ]]; then
  BLAST_DB="${BLAST_DB_DIR}/SILVA_138_SSURef_NR99_tax_silva.fasta"
fi

if [[ -z "${TAXA_MAP}" ]]; then
  TAXA_MAP="${BLAST_DB_DIR}/taxmap_slv_ssu_ref_nr_138.txt"
fi

if [[ -z "${NSLOTS}" ]]; then
  NSLOTS="12"
fi

if [[ -z "${SAVE_WORKSPACE}" ]]; then
  SAVE_WORKSPACE="T"
fi

if [[ -z "${OVERWRITE}" ]]; then
  OVERWRITE="F"
fi

if [[ -z "${BLOUT}" ]]; then
  BLOUT=$(mktemp)
  TMP_BLOUT=1
fi


###############################################################################
# 5. Check input asv file
###############################################################################

if [[ ! -f "${INPUT_ASV_TABLE}" ]]; then
  echo "no input asv table"
  exit 1
fi

###############################################################################
# 6. Check output 
###############################################################################

if [[ -f "${OUTPUT_ASV_TABLE}" && "${OVERWRITE}" == "T" ]]; then
  rm "${OUTPUT_ASV_TABLE}"
fi

if [[ -f "${OUTPUT_ASV_TABLE}" && "${OVERWRITE}" == "F" ]]; then
  echo "${OUTPUT_ASV_TABLE} already exists. Use --overwrite T to overwrite"
  exit
fi

###############################################################################
# 5. Create fasta and seq map files
###############################################################################

if [[ ! -f "${INPUT_FASTA}" ]]; then

  INPUT_FASTA=$(mktemp)
  cut -f1 -d"," "${INPUT_ASV_TABLE}" | \
  awk 'NR>1 {gsub(/\"/,"",$0); 
             n++; printf "%s\n%s\n", ">asv_"n,$0 }' > "${INPUT_FASTA}"
             
  if [[ $? -ne "0" ]]; then
    echo "create fasta file failed"
    exit 1
  fi  

  TMP_FASTA=1  
fi 

SEQ_MAP=$(mktemp)      
awk 'BEGIN {printf "%s\t%s\n", "qseqid", "asv"} 
     $0 ~ />/ {header = $0; sub(">","",header) }; 
     $0 !~ />/ { seq = $0; printf "%s\t%s\n", header,seq}' \
     "${INPUT_FASTA}" > "${SEQ_MAP}"
       
if [[ $? -ne "0" ]]; then
  echo "create seq map file failed"
  exit 1
fi  

###############################################################################
# 6. Run rscript
###############################################################################

Rscript --vanilla \
"${taxa_annot}" \
"${INPUT_ASV_TABLE}" \
"${INPUT_FASTA}" \
"${OUTPUT_ASV_TABLE}" \
"${METHOD}" \
"${TRAIN_DB}" \
"${REF_DB}" \
"${BLAST_DB}" \
"${TAXA_MAP}" \
"${SEQ_MAP}" \
"${EVALUE}" \
"${MIN_IDENTITY}" \
"${BLOUT}" \
"${NSLOTS}" \
"${SAVE_WORKSPACE}"

if [[ $? != "0" ]]; then
  echo "taxa_annot.R failed"
  rm -r "${OUTPUT_ASV_TABLE}"
  exit 1
fi

###############################################################################
# 7. Clean
###############################################################################

if [[ -f "${INPUT_FASTA}" && "${TMP_FASTA}" -eq "1" ]]; then
  rm "${INPUT_FASTA}"
fi

if [[ -f "${BLOUT}" && "${TMP_BLOUT}" -eq "1" ]]; then
  rm "${BLOUT}"
fi

rm "${SEQ_MAP}"
