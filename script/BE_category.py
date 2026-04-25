###this script focuses on the 40bp quantification window
###classify reads into wildtype, indel, inside_desired, inside_desired_undesired, inside_undesired, outside, inside and outside categories (7 categories in total), and calculate its ratio
import argparse
import os
import pandas as pd
import zipfile

parser = argparse.ArgumentParser(description="Get statistics from Alleles_frequency_table.zip")
parser.add_argument("-f","--CRISPResso2_folder",type=str,help="CRISPResso output folder containing finished analysis",required=True)
parser.add_argument("-w1","--edit_window",type=str,help="Editing window coordinate, eg, 20-25",required=True)
parser.add_argument("-w2","--wide_window",type=str,help="Wider window coordinate, eg, 10-50",required=True)
parser.add_argument("-d","--desired_edit",type=str,help="Desired edit, CT for CBE, AG for ABE",required=True)
args = parser.parse_args()

def snv_type(in_desire,in_undesire,out_desire,out_undesire):
    in_total = in_desire + in_undesire
    out_total = out_desire + out_undesire
    if in_total==0 and out_total==0:                       #=0,=0,=0,=0
        flag="wildtype"
    elif in_total==0 and out_undesire>0:                   #=0,=0,=0,>0  =0,=0,>0,>0
        flag="outside_other"
    elif in_total==0 and out_desire>0:                     #=0,=0,>0,=0
        flag="outside_desire"
    elif in_desire==0 and in_undesire>0 and out_total==0:  #=0,>0,=0,=0
        flag="inside_undesire"
    elif in_desire>0 and in_undesire==0 and out_total==0:  #>0,=0,=0,=0
        flag="inside_desire"
    elif in_undesire>0 and out_total==0:                   #>0,>0,=0,=0
        flag="inside_desire_undesire"
    elif in_undesire>0 and out_total>0:                    #>0,>0,>0,>0  >0,>0,=0,>0  >0,>0,>0,=0  =0,>0,>0,>0  =0,>0,=0,>0  =0,>0,>0,=0  >0,>0,=0,>0  
        flag="inside_outside_other"
    elif in_total>0 and out_undesire>0:                    #>0,0,0,>0
        flag="inside_outside_other"
    elif in_desire>0 and out_desire>0:                     #>0,=0,>0,=0
        flag="inside_outside_desire"
    return flag

def key_value(df,key):
    value=df[key] if key in df.keys() else 0
    return value

edit_win=args.edit_window
wide_win=args.wide_window
desired_edit=args.desired_edit
sample=os.path.basename(args.CRISPResso2_folder)
sample=sample.replace("CRISPResso_on_","")

df_map = pd.read_csv(os.path.join(args.CRISPResso2_folder,"CRISPResso_mapping_statistics.txt"),sep="\t")
n_input=df_map['READS AFTER PREPROCESSING'] #should be READS IN INPUTS
n_aligned=df_map['READS ALIGNED']

z = zipfile.ZipFile(os.path.join(args.CRISPResso2_folder,"Alleles_frequency_table.zip"))
zf=z.open("Alleles_frequency_table.txt")
df_alleles = pd.read_csv(zf,sep="\t")

edit_win_start=int(edit_win.split("-")[0])
edit_win_end=int(edit_win.split("-")[1])
wide_win_start=int(wide_win.split("-")[0])
wide_win_end=int(wide_win.split("-")[1])

df_alleles['flag']=''
flag_idx=df_alleles.columns.get_loc("flag")

for idx in df_alleles.index:
    snv_pos=df_alleles['all_substitution_positions'][idx][1:-1]
    ins_pos=df_alleles['all_insertion_positions'][idx][1:-1]
    del_pos=df_alleles['all_deletion_positions'][idx][1:-1]
    if ins_pos or del_pos:  #combine insertion and deletion positions across whole amplicon
        if ins_pos and del_pos:
            all_indel=ins_pos+','+del_pos
        else:
            all_indel=ins_pos+del_pos
        n_indel_win=0
        for pos in all_indel.split(','):
            pos=int(float(pos))
            if pos in range(wide_win_start,wide_win_end+1):  
                n_indel_win+=1
        if n_indel_win>0:   #indel within 40bp win
            df_alleles.iloc[idx,flag_idx]="indel"
        else:     #no indel within 40bp win
            if snv_pos:   #snv exists
                n_in_desire=0
                n_in_undesire=0
                n_out_desire=0
                n_out_undesire=0
                #print(df_alleles['Aligned_Sequence'][idx][edit_win_start,edit_win_end+1],end="\t")
                #print(df_alleles['Reference_Sequence'][idx][edit_win_start,edit_win_end+1])
                for pos in snv_pos.split(','):
                    pos=int(float(pos))
                    amplicon_base=df_alleles['Aligned_Sequence'][idx][pos]
                    ref_base=df_alleles['Reference_Sequence'][idx][pos]
                    if pos in range(edit_win_start,edit_win_end+1):  #within canonical win
                        if amplicon_base == desired_edit[1] and ref_base == desired_edit[0]:  #desired edit
                            n_in_desire+=1
                        else:
                            n_in_undesire+=1
                    elif pos in range(wide_win_start,edit_win_start) or pos in range(edit_win_end+1,wide_win_end+1): #without canonical win but within 40bp win
                        if amplicon_base == desired_edit[1] and ref_base == desired_edit[0]:  #desired edit
                            n_out_desire+=1
                        else:
                            n_out_undesire+=1
                df_alleles.iloc[idx,flag_idx]=snv_type(n_in_desire,n_in_undesire,n_out_desire,n_out_undesire)
            else:          #no snv
                df_alleles.iloc[idx,flag_idx]="wildtype"
    else:         #no indel across whole amplicom
        if snv_pos:
            n_in_desire=0
            n_in_undesire=0
            n_out_desire=0
            n_out_undesire=0
            for pos in snv_pos.split(','):
                pos=int(float(pos))
                amplicon_base=df_alleles['Aligned_Sequence'][idx][pos]
                ref_base=df_alleles['Reference_Sequence'][idx][pos]
                if pos in range(edit_win_start,edit_win_end+1):
                    if amplicon_base == desired_edit[1] and ref_base == desired_edit[0]:
                        n_in_desire+=1
                    else:
                        n_in_undesire+=1
                elif pos in range(wide_win_start,edit_win_start) or pos in range(edit_win_end+1,wide_win_end+1):
                    if amplicon_base == desired_edit[1] and ref_base == desired_edit[0]:  #desired edit
                        n_out_desire+=1
                    else:
                        n_out_undesire+=1
            df_alleles.iloc[idx,flag_idx]=snv_type(n_in_desire,n_in_undesire,n_out_desire,n_out_undesire)
        else:
            df_alleles.iloc[idx,flag_idx]="wildtype"

df=df_alleles[['#Reads','flag']].groupby(['flag']).sum()['#Reads']
n_total=df.sum()
total_indel=key_value(df,'indel')
total_wildtype=key_value(df,'wildtype')
total_inside_desire=key_value(df,'inside_desire')
total_inside_undesire=key_value(df,'inside_undesire')
total_inside_desire_undesire=key_value(df,'inside_desire_undesire')
total_inside_outside_desire=key_value(df,'inside_outside_desire')
total_inside_outside_other=key_value(df,'inside_outside_other')
total_outside_desire=key_value(df,'outside_desire')
total_outside_other=key_value(df,'outside_other')

r_indel=100*total_indel/n_total
r_wt=100*total_wildtype/n_total
r_in_des=100*total_inside_desire/n_total
r_in_undes=100*total_inside_undesire/n_total
r_in_des_undes=100*total_inside_desire_undesire/n_total
r_in_out_des=100*total_inside_outside_desire/n_total
r_in_out_other=100*total_inside_outside_other/n_total
r_out_des=100*total_outside_desire/n_total
r_out_other=100*total_outside_other/n_total

print("%s\t%d\t%d" %(sample,n_input,n_aligned),end='\t')  #mapping reads
print("%d\t%d\t%4.2f\t%d\t%4.2f\t%d\t%4.2f\t%d\t%4.2f" %(n_total,total_indel,r_indel,total_wildtype,r_wt,total_inside_desire,r_in_des,total_inside_undesire,r_in_undes),end='\t') 
print("%d\t%4.2f\t%d\t%4.2f\t%d\t%4.2f\t%d\t%4.2f\t%d\t%4.2f" %(total_inside_desire_undesire,r_in_des_undes,total_inside_outside_desire,r_in_out_des,total_inside_outside_other,r_in_out_other,total_outside_desire,r_out_des,total_outside_other,r_out_other))
