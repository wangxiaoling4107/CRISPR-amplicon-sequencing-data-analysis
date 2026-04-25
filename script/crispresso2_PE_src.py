import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description="Print PE2 scripts")
parser.add_argument("-f","--info",type=str,help="Sample information, xlsx file",required=True)
args = parser.parse_args()

df_info = pd.read_excel(args.info, sheet_name="info",header=0)
df_gRNA = pd.read_excel(args.info, sheet_name="gRNA",header=0)
df_amp = pd.read_excel(args.info, sheet_name="amplicon",header=0)
df_ext = pd.read_excel(args.info, sheet_name="extention",header=0)
df_sca = pd.read_excel(args.info, sheet_name="scaffold",header=0)

dict_gRNA = {}
for index,row in df_gRNA.iterrows():
    dict_gRNA[row['Name']] = row['Sequence'].upper()

dict_amplicon = {}
for index,row in df_amp.iterrows():
    dict_amplicon[row['Primer']] = row['Sequence'].upper()

dict_extension = {}
for index,row in df_ext.iterrows():
    dict_extension[row['Name']] = row['Sequence'].upper()

dict_scaffold = {}
for index,row in df_sca.iterrows():
    dict_scaffold[row['Name']] = row['Sequence'].upper()

script_dir = os.path.join(os.getcwd(),'script')
if not os.path.exists(script_dir):
    os.makedirs(script_dir)

result_dir = os.path.join(os.getcwd(),'result')
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

for index,row in df_info.iterrows():
    sampleID = row['GE_Tool'] + "_" + row['gRNA'] + "_" + row['Condition'] + "_" + str(row['NGS_Code'])
    r1 = os.path.join(os.getcwd(),'data',row['Read1'])
    r2 = os.path.join(os.getcwd(),'data',row['Read2'])
    f = open(os.path.join(script_dir, sampleID + ".sh"), "w")
    f.write('#!/bin/sh\n')
    f.write('CRISPResso --fastq_r1 %s \\\n' % (r1))
    f.write('--fastq_r2 %s \\\n' % (r2))
    f.write('--amplicon_seq %s \\\n' % (dict_amplicon[row['Primer']]))
    f.write('--amplicon_name %s \\\n' % (row['Primer']))
    f.write('--prime_editing_pegRNA_spacer_seq %s \\\n' % (dict_gRNA[row['gRNA']]))
    f.write('--prime_editing_pegRNA_extension_seq %s \\\n' % (dict_extension[row['Extention']]))
    f.write('--prime_editing_pegRNA_scaffold_seq %s \\\n' % (dict_scaffold[row['Scaffold']]))
    f.write('--prime_editing_pegRNA_extension_quantification_window_size 10 \\\n')
    f.write('--prime_editing_pegRNA_scaffold_min_match_length 3 \\\n')
    f.write('--trim_sequences \\\n')
    f.write('--trimmomatic_options_string ILLUMINACLIP:/home/xwang3/project/bin/adapter.fa:0:90:10:0:true \\\n')
    f.write('--min_average_read_quality 20 \\\n')
    f.write('--amplicon_min_alignment_score 30 \\\n')
    f.write('--min_paired_end_reads_overlap 50 \\\n')
    f.write('--max_paired_end_reads_overlap 300 \\\n')
    f.write('--write_detailed_allele_table \\\n')
    f.write('--write_cleaned_report \\\n')
    f.write('--place_report_in_output_folder \\\n')
    f.write('--suppress_report \\\n')
    f.write('-n %s \\\n' % (sampleID))
    f.write('-o %s \n' % (result_dir))
