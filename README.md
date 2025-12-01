# Amplicon Pipelines

This repo includes pipelines for the preprocessing, clustering and annotation of amplicon sequence data. 

## Repository Structure

```
.
├── LICENSE                  # License file
├── README.md               # This file
├── environment.yml         # Conda/Mamba environment specification
├── documentation/          # Additional documentation
│   └── documentation.md   # Detailed pipeline documentation
└── scripts/               # Main analysis scripts
    ├── dada2_pipeline.R   # DADA2 ASV generation pipeline
    ├── taxa_annot.R       # Taxonomic annotation pipeline
    ├── toolbox.R          # Utility functions
    ├── check_primers.py   # Primer checking utility
    └── conf.sh           # Configuration file
```

## Installation

All dependencies can be installed using [Mamba](https://mamba.readthedocs.io/) (or Conda) with the provided environment file:

```bash
# Install the environment
mamba env create -f environment.yml

# Activate the environment
conda activate metabarcoding-processing-obm
```



