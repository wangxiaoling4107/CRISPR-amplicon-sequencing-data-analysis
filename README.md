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
3. Generate scripts to process CRISPResso2 output for BE analysis, run and merge results
```text
  BE_category_src.sh sample.info.csv
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

## Relationship to manuscript figures

The scripts in this repository were used to process and quantify
amplicon sequencing data underlying the manuscript figures.
Some quantitative results were subsequently visualized using GraphPad Prism.

### Base editing

| Manuscript figure(s) | Script | Description |
|---|---|---|
| Figure 6B, 6C, Figure S14, S17, S18C | `script/BE_category.py` | Classification and quantification of base-editing outcomes. The resulting data were subsequently visualized using GraphPad Prism. |
| Figure 7 | `script/heatmap_BE.py` | Generation of the base-editing sequencing result heatmap. |

### Prime editing

| Manuscript figure(s) | Script | Description |
|---|---|---|
| Figure 4F, Figure S7B, S8B–C, S9B–C, S10C | `script/crispresso2_PE_src.py` and `script/merge_PE_result.py` | Processing and merging of prime-editing amplicon sequencing results. The resulting data were subsequently visualized using GraphPad Prism. |


