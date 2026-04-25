## Overview
This repository provides scripts for the analysis of base editing (BE) and prime editing (PE) sequencing data. FASTQ files are initially processed using CRISPResso2. For base editing applications, custom in-house scripts further process CRISPResso2 outputs, categorizing editing events to facilitate quantitative evaluation of editing efficiency and product purity.
## Workflow
1. Prepare metadata and run CRISPResso2
2. Process CRISPResso2 output file (Alleles_frequency_table.zip)
3. Run BE analysis scripts
4. Merge result and generate visualizations

## Instructions
1. Prepare metadata (see meta). 
2. Generate scripts for CRISPResso2 analysis and run. e.g. crispresso2_BE_src.py -f info_BE.xlsx or crispresso2_PE_src.py -f info_PE.xlsx
3. Generate scripts to process CRISPResso2 output for BE analysis and run to merge result. e.g. BE_category_src.sh
4. Visualize BE result, histogram_BE.py -f info_BE.xlsx; heatmap_BE.py -f info_BE.xlsx
5. Merge PE result, merge_PE_result.py -f info_PE.xlsx

## Example Output
See result folder

## Requirements
- Python 3.9+
- CRISPResso2

## Author
Xiaoling Wang
