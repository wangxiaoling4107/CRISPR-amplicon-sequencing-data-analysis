# CRISPR Amplicon Sequencing Data Analysis

## Overview

This repository contains scripts used to analyse amplicon sequencing data generated to evaluate CRISPR-based genome editing, including base editing (BE) and prime editing (PE). 

FASTQ files are initially processed using CRISPResso2. Custom Python and shell scripts are then used to process CRISPResso2 outputs, classify editing outcomes, quantify editing efficiency and product purity, and generate visualizations.

## Workflow

FASTQ
  ↓
CRISPResso2
  ↓
CRISPResso2 output
  ↓
Custom analysis scripts
  ↓
Editing outcome classification
  ↓
Quantification
  ↓
Visualization

## Repository structure

data/      Input/example data (fastq)
meta/      Metadata
script/    Analysis scripts
result/    Example results

## Example Output
See result folder

## Requirements
- Python 3.9+
- CRISPResso2
- pandas
- numpy
- matplotlib
- seaborn

## Instructions
1. Prepare metadata (format see meta/)
2. Process metadata, generate scripts for CRISPResso2 analysis and run them
   e.g. crispresso2_BE_src.py -f info_BE.xlsx or crispresso2_PE_src.py -f info_PE.xlsx
3. Generate scripts to process CRISPResso2 output for BE analysis and run to merge results
   e.g. BE_category_src.sh
4. Visualize BE results
   e.g. histogram_BE.py -f info_BE.xlsx
   e.g. heatmap_BE.py -f info_BE.xlsx
5. Merge PE results
   e.g. merge_PE_result.py -f info_PE.xlsx

## Reproducibility

The version of this repository corresponding to the
manuscript is archived at Zenodo:

DOI: [to be added]

## Manuscript

The scripts in this repository were used to generate
the analyses presented in:

[manuscript citation / title]

## Author
Xiaoling Wang
