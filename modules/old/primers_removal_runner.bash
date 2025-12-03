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
--input_dir CHAR                directory with the input raw fastq files
--output_dir CHAR               directory to output generated data (i.e., preprocessed data, plots, tables)
--pattern_r1 CHAR               patter of R1 reads to load fastq files (default _L001_R1_001.fastq.gz)
--pattern_r2 CHAR               patter of R2 reads to load fastq files (default _L001_R2_001.fastq.gz)
--primer_fwd                    forward primer sequence used in PCR
--primer_rev                    reverse primer sequence used in PCR
--nslots NUM                    number of threads used (default 12)
--save_workspace T|F            save R workspace image
--overwrite T|F                 overwrite previous directory
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
  --input_dir)
  if [[ -n "${2}" ]]; then
    INPUT_DIR="${2}"
    shift
  fi
  ;;
  --input_dir=?*)
  INPUT_DIR="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --input_dir=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --output_dir)
  if [[ -n "${2}" ]]; then
    OUTPUT_DIR="${2}"
    shift
  fi
  ;;
  --output_dir=?*)
  OUTPUT_DIR="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --output_dir=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --pattern_r1)
  if [[ -n "${2}" ]]; then
    PATTERN_R1="${2}"
    shift
  fi
  ;;
  --pattern_r1=?*)
  PATTERN_R1="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --pattern_r1=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --pattern_r2)
  if [[ -n "${2}" ]]; then
    PATTERN_R2="${2}"
    shift
  fi
  ;;
  --pattern_r2=?*)
  PATTERN_R2="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --pattern_r2=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --primer_fwd)
  if [[ -n "${2}" ]]; then
    PRIMER_FWD="${2}"
    shift
  fi
  ;;
  --primer_fwd=?*)
  PRIMER_FWD="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --primer_fwd=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;; 
#############
  --primer_rev)
  if [[ -n "${2}" ]]; then
    PRIMER_REV="${2}"
    shift
  fi
  ;;
  --primer_rev=?*)
  PRIMER_REV="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --primer_rev=) # Handle the empty case
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

if [[ -z "${PATTERN_R1}" ]]; then
  PATTERN_R1="_L001_R1_001.fastq.gz"
fi

if [[ -z "${PATTERN_R2}" ]]; then
  PATTERN_R2="_L001_R2_001.fastq.gz"
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

###############################################################################
# 5. Check input dir
###############################################################################

if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "invalid input directory"
  exit 1
fi

###############################################################################
# 6. Check output dir
###############################################################################

if [[ -d "${OUTPUT_DIR}" && "${OVERWRITE}" == "T" ]]; then
  rm -r "${OUTPUT_DIR}"
fi

if [[ -d "${OUTPUT_DIR}" && "${OVERWRITE}" == "F" ]]; then
  echo "${OUTPUT_DIR} already exists. Use --overwrite T to overwrite"
  exit
fi

###############################################################################
# 5. Run rscript
###############################################################################

Rscript --vanilla \
"${primers_removal}" \
"${INPUT_DIR}" \
"${OUTPUT_DIR}" \
"${PATTERN_R1}" \
"${PATTERN_R2}" \
"${PRIMER_FWD}" \
"${PRIMER_REV}" \
"${NSLOTS}" \
"${SAVE_WORKSPACE}"

if [[ $? != "0" ]]; then
  echo "primers_removal.R failed"
  rm -r "${OUTPUT_DIR}"
  exit 1
fi
