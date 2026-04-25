import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description="Print BE scripts")
parser.add_argument("-f","--info",type=str,help="Sample information, xlsx file",required=True)
args = parser.parse_args()

df_info = pd.read_excel(args.info, sheet_name="info",header=0)
df_gRNA = pd.read_excel(args.info, sheet_name="gRNA",header=0)
df_amp = pd.read_excel(args.info, sheet_name="amplicon",header=0)

dict_gRNA = {}
for index,row in df_gRNA.iterrows():
    dict_gRNA[row['Name']] = row['Sequence'].upper()

dict_amplicon = {}
for index,row in df_amp.iterrows():
    dict_amplicon[row['Primer']] = row['Sequence'].upper()

script_dir = os.path.join(os.getcwd(),'scripts')
if not os.path.exists(script_dir):
    os.makedirs(script_dir)

result_dir = os.path.join(os.getcwd(),'result')
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

for index,row in df_info.iterrows():
    sampleID = row['GE_Tool'] + "_" + row['gRNA'] + "_" + row['Condition'] + "_" + str(row['NGS_Code'])
    r1 = os.path.join(os.getcwd(),'data',row['Read1'])
    r2 = os.path.join(os.getcwd(),'data',row['Read2'])
    gRNA_name = row['gRNA'].replace('.','-')
    f = open(os.path.join(script_dir, sampleID + ".sh"), "w")
    f.write('#!/bin/sh\n')
    f.write('CRISPResso --fastq_r1 %s \\\n' % (r1))
    f.write('--fastq_r2 %s \\\n' % (r2))
    f.write('--amplicon_seq %s \\\n' % (dict_amplicon[row['Primer']]))
    f.write('--amplicon_name %s \\\n' % (row['Primer']))
    f.write('--guide_seq %s \\\n' % (dict_gRNA[row['gRNA']]))
    f.write('--guide_name %s \\\n' % (gRNA_name))
    f.write('--trim_sequences \\\n')
    f.write('--trimmomatic_options_string ILLUMINACLIP:/home/xwang3/project/bin/adapter.fa:0:90:10:0:true \\\n')
    f.write('--min_average_read_quality 20 \\\n')
    f.write('--min_paired_end_reads_overlap 50 \\\n')
    f.write('--max_paired_end_reads_overlap 300 \\\n')
    f.write('--default_min_aln_score 80 \\\n')
    f.write('--quantification_window_coordinates %s \\\n' % (row['Win40_gCenter']))
    f.write('--write_detailed_allele_table \\\n')
    f.write('--base_editor_output \\\n')
    f.write('--suppress_report \\\n')
    f.write('--exclude_bp_from_left 0 \\\n')
    f.write('--exclude_bp_from_right 0 \\\n')
    f.write('--quantification_window_center -10 \\\n')
    f.write('-n %s \\\n' % (sampleID))
    f.write('-o %s \n' % (result_dir))
