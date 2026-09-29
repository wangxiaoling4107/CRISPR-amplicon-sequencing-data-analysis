# CRISPR Amplicon Sequencing Data Analysis

## Overview

This repository contains scripts used to analyse amplicon sequencing data generated to evaluate gene editing outcomes of base editors (BE) and prime editors (PE). 

FASTQ files are initially processed using CRISPResso2. Custom Python and shell scripts are then used to process CRISPResso2 outputs, classify editing outcomes, quantify editing efficiency and product purity, and generate visualizations.

## Workflow

```text
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
```

## Repository structure

```text
data/      Input example (fastq)
meta/      Metadata example
script/    Analysis scripts
result/    Example results
```

## Example Output
See result/

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
```text
  crispresso2_BE_src.py -f info_BE.xlsx
  crispresso2_PE_src.py -f info_PE.xlsx
```
3. Generate scripts to process CRISPResso2 output for BE analysis and run to merge results
```text
  BE_category_src.sh
```
4. Visualize BE results
```text
  histogram_BE.py -f info_BE.xlsx
  heatmap_BE.py -f info_BE.xlsx
```
5. Merge PE results
```text
  merge_PE_result.py -f info_PE.xlsx
```

## Reproducibility

The version of this repository corresponding to the
manuscript is archived at Zenodo:

DOI: [to be added]

## Manuscript

The scripts in this repository were used to generate
the analyses presented in:

*The chromatin context differently impacts prime editors and base editors and controls the fidelity and purity of base editing*

## Author
Xiaoling Wang
