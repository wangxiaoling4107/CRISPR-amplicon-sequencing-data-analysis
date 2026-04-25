import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description="Get PE results")
parser.add_argument("-f","--info",type=str,help="Sample information, xlsx file",required=True)
args = parser.parse_args()

def read_edit_frequency(res):
    df = pd.read_csv(res,sep="\t")
    total_input = df['Reads_in_input'][0]
    total_aligned = df['Reads_aligned_all_amplicons'][0]
    Ref_aligned = df['Reads_aligned'][0]
    PE_aligned = df['Reads_aligned'][1]
    Scaffold_incorporated = df['Reads_aligned'][2]
    Ref_unmodified = df['Unmodified'][0]
    PE_unmodified = df['Unmodified'][1]
    Ref_modified = df['Modified'][0]
    PE_modified = df['Modified'][1]
    ambiguous = df['Reads_aligned_all_amplicons'][0] - sum(df['Reads_aligned'])
    if total_aligned != 0:
        Ref_unmodified_ratio = (100*Ref_unmodified/total_aligned).round(1)
        Ref_modified_ratio = (100*Ref_modified/total_aligned).round(1) #revised 04-Jun-2024
        PE_unmodified_ratio = (100*PE_unmodified/total_aligned).round(1)
        PE_modified_ratio = (100*PE_modified/total_aligned).round(1)
        Scaffold_incorporated_ratio = (100*Scaffold_incorporated/total_aligned).round(1)
        ambiguous_ratio = (100*ambiguous/total_aligned).round(1)
    print(total_input,"\t",total_aligned,"\t",Ref_aligned,"\t",Ref_unmodified,"\t",Ref_unmodified_ratio,"\t",Ref_modified,"\t",Ref_modified_ratio,"\t",PE_aligned,"\t",PE_unmodified,"\t",PE_unmodified_ratio,"\t",PE_modified,"\t",PE_modified_ratio,"\t",Scaffold_incorporated,"\t",Scaffold_incorporated_ratio,"\t",ambiguous,"\t",ambiguous_ratio)

df_info = pd.read_excel(args.info, sheet_name="info",header=0)
result_dir = os.path.join(os.getcwd(),'result')
if not os.path.exists(result_dir):
    print('Result path not exists!')

##print title
print("Total_input\tTotal_aligned\tRef_aligned\tRef_unmodified\tRef_unmodified%\tRef_modified\tRef_modified%\tPE_aligned\tPE_unmodified\tPE_unmodified%\tPE_modified\tPE_modified%\tScaffold_incorporated\tScaffold_incorporated%\tAmbigous\tAmbigous%")

for index,row in df_info.iterrows():
    sampleID = row['GE_Tool'] + "_" + row['gRNA'] + "_" + row['Condition'] + "_" + str(row['NGS_Code'])
    sampleID = sampleID.replace(".","_")
    res = os.path.join(result_dir,"CRISPResso_on_" + sampleID, "CRISPResso_quantification_of_editing_frequency.txt")
    if not os.path.exists(res):
        print("{} file not exists!".format(sampleID))
    else:
        print(sampleID,end="\t")
        read_edit_frequency(res)
