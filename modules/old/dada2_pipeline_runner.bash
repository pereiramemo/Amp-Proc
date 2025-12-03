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
--nslots NUM                    number of threads used (default 12)
--trunc_r1 NUM                  number of nuc to remove in R1 from the 3' end (default 250)
--trunc_r2 NUM                  number of nuc to remove in R2 from the 3' end (default 200)
--pattern_r1 CHAR               patter of R1 reads to load fastq files (default _L001_R1_001.fastq.gz)
--pattern_r2 CHAR               patter of R2 reads to load fastq files (default _L001_R2_001.fastq.gz)
--bimeras_method CHAR           method to check bimeras: pooled, consensus, per-sample (default consensus)
--min_overlap NUM               minimum number of nucletides to overlap in merging (default 12)
--pooled T|F                    use pooled option when running dada2 (default T)
--qual_plot T|F                 create quality plots
--err_plot T|F                  create error plots
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
  --trunc_r1)
  if [[ -n "${2}" ]]; then
    TRUNC_R1="${2}"
    shift
  fi
  ;;
  --trunc_r1=?*)
  TRUNC_R1="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --trunc_r1=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --trunc_r2)
  if [[ -n "${2}" ]]; then
    TRUNC_R2="${2}"
    shift
  fi
  ;;
  --trunc_r2=?*)
  TRUNC_R2="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --trunc_r2=) # Handle the empty case
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
  --bimeras_method)
  if [[ -n "${2}" ]]; then
    BIMERAS_METHOD="${2}"
    shift
  fi
  ;;
  --bimeras_method=?*)
  BIMERAS_METHOD="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --bimeras_method=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --min_overlap)
  if [[ -n "${2}" ]]; then
    MIN_OVERLAP="${2}"
    shift
  fi
  ;;
  --min_overlap=?*)
  MIN_OVERLAP="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --min_overlap=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --pooled)
  if [[ -n "${2}" ]]; then
    POOLED="${2}"
    shift
  fi
  ;;
  --pooled=?*)
  POOLED="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --pooled=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --qual_plot)
  if [[ -n "${2}" ]]; then
    QUAL_PLOT="${2}"
    shift
  fi
  ;;
  --qual_plot=?*)
  QUAL_PLOT="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --qual_plot=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --err_plot)
  if [[ -n "${2}" ]]; then
    ERR_PLOT="${2}"
    shift
  fi
  ;;
  --err_plot=?*)
  ERR_PLOT="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --err_plot=) # Handle the empty case
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

if [[ -z "${NSLOTS}" ]]; then
  NSLOTS="12"
fi

if [[ -z "${TRUNC_R1}" ]]; then
  TRUNC_R1="250"
fi

if [[ -z "${TRUNC_R2}" ]]; then
  TRUNC_R2="200"
fi

if [[ -z "${PATTERN_R1}" ]]; then
  PATTERN_R1="_L001_R1_001.fastq.gz"
fi

if [[ -z "${PATTERN_R2}" ]]; then
  PATTERN_R2="_L001_R2_001.fastq.gz"
fi

if [[ -z "${MIN_OVERLAP}" ]]; then
  MIN_OVERLAP="12"
fi

if [[ -z "${BIMERAS_METHOD}" ]]; then
  BIMERAS_METHOD="consensus"
fi

if [[ -z "${POOLED}" ]]; then
  POOLED="T"
fi

if [[ -z "${QUAL_PLOT}" ]]; then
  QUAL_PLOT="T"
fi

if [[ -z "${ERR_PLOT}" ]]; then
  ERR_PLOT="T"
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
"${dada2_pipeline}" \
"${INPUT_DIR}" \
"${OUTPUT_DIR}" \
"${PATTERN_R1}" \
"${PATTERN_R2}" \
"${NSLOTS}" \
"${TRUNC_R1}" \
"${TRUNC_R2}" \
"${MIN_OVERLAP}" \
"${BIMERAS_METHOD}" \
"${POOLED}" \
"${QUAL_PLOT}" \
"${ERR_PLOT}" \
"${SAVE_WORKSPACE}"

if [[ $? != "0" ]]; then
  echo "dada2_pipeline.R failed"
  rm -r "${OUTPUT_DIR}"
  exit 1
fi
