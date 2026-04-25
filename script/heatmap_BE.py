import os
import glob
import zipfile
import argparse
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

parser = argparse.ArgumentParser(description="Plot BE heatmap...")
parser.add_argument("-f","--info",type=str,help="Sample information for plotting, xlsx file",required=True)
args = parser.parse_args()

def arrStr_to_arr(val):
    return [int(x) for x in val[1:-1].split(",")]

def get_complementary_seq(seq):
    complement_seq=''
    pair={'A':'T','G':'C','C':'G','T':'A'}
    for base in seq:
        complement_seq += pair[base]
    return complement_seq

#Calculate base percentages of each position in wide window, reads with InDels or Ns are excluded first.
def read_single_df(zip_file,pos):
   # print('Reading {}'.format(zip_file))
    z = zipfile.ZipFile(zip_file)
    zf=z.open("Alleles_frequency_table.txt")
    df=pd.read_csv(zf,sep="\t")
    w_start=int(pos[2]) #wide_window start for reference sequence
    w_end=int(pos[3])+1 #wide_window end for reference sequence
    ref=df['Reference_Sequence'][0][w_start:w_end]
    df['ref_positions']=df['ref_positions'].apply(arrStr_to_arr)
    frequency_table=[]
    for i in range(len(ref)):
        frequency_table.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})
    for idx, row in df.iterrows():
        # get wide window start for each aligned sequence; the wide window might be shifted if there are insertions before it
        row_aligned_w_start=row['ref_positions'].index(w_start) 
        row_aligned_w_end=row_aligned_w_start+40
        row_align_seq=row['Aligned_Sequence'][row_aligned_w_start:row_aligned_w_end]
        row_ref_seq=row['Reference_Sequence'][row_aligned_w_start:row_aligned_w_end]
        #print(row_align_seq,end="\t")
        #print(row_ref_seq)
        if "-" in row_align_seq or "-" in row_ref_seq: #exclude reads with InDels
            continue
        if "N" in row_align_seq: #exclude reads with Ns
            continue
        for position,base in enumerate(row_align_seq):
            frequency_table[position][base] += row['#Reads']
    nLines = float(df['#Reads'].sum()) #total reads
    df_frequency=pd.DataFrame(frequency_table)
    df_row_sum=df_frequency.sum(axis=1)
    df_percent=df_frequency.div(df_row_sum,axis='rows').round(4).T
    #df_percent[df_percent>0.50]=0  
    #arbitrary; set 0.5 as shreshold, because all editing percentages are below 40%; in this way, bases which are the same with reference will not show on the picture. 
    df_percent=df_percent*100
    df_percent.columns=list(ref) #column names
    for idx, row in df_percent.iterrows():
        df_percent.loc[idx,idx]=0 #set 0 if no mutations
    return df_percent

def figure(df,ax,color):
#    plt.figure(figsize = (4,9)) #figure size
    y_axis_labels=[word[0] for word in df.index.array] #keep the first string in label
    ax=sns.heatmap(data=df,cmap="Blues",annot=True,yticklabels=y_axis_labels,linewidths=0.5,linecolor='gray',vmin=0, vmax=25) #set x,y limits
    ax.set_title(k)
    plt.yticks(rotation=0)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 10)
    ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize = 10)
    for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(),row_colors):
        ticklabel.set_color(tickcolor)
    l_width=2 #set wideth of heatmap frame line
    ax.axhline(y=40, color='gray',linewidth=l_width)
    ax.axvline(x=4, color='gray',linewidth=l_width)
    return ax
#    plt.savefig(filename)
#    plt.cla()

path=os.path.dirname(os.getcwd())
df_info = pd.read_excel(args.info, sheet_name="info",header=0)
dict_exp={}
for index,row in df_info.iterrows():
    #exp=sample+"_"+dox+"_"+code #
    exp=row['GE_Tool']+"_"+row['gRNA']+"_"+row['Condition']+"_"+str(row['NGS_Code'])
    infos=row['WinEdit']+"-"+row['Win40_gCenter']+"-"+row['gRNA_win']+"-"+row['gRNA_Strand']
    dict_exp[exp]={}
    #infos=items[4]+"-"+items[5]+"-"+items[8]+"-"+strand
    dict_exp[exp]=infos #editWin-wideWin-gRNAwin-strand

plt.figure(figsize=(4*2, 9*len(dict_exp)/2))
count_fig = 0
for k,v in dict_exp.items():
    # get one zip file for each replicate
    res=glob.glob('{}/result/CRISPResso_on_{}/Alleles_frequency_table.zip'.format(path,k.replace(".","_"))) 
    pos=v.split("-")[0:-1]
    pos=[ int(x) for x in pos ]
    df=read_single_df(res[0],pos)
    df=(df).round(1)  #average
    strand=v.split("-")[-1].strip()
    if 'reverse' in strand: #use corresponding complementary sequence if it is a reverse strand, so it is always A>G for ABE, C>T for CBE
        df.columns=list(get_complementary_seq(df.columns))
        df.index=list(get_complementary_seq(df.index))
    df=df.sort_index()
    new_edit_start=pos[0]-pos[2] #editint window position relative to the first base of wide window
    new_edit_end=pos[1]-pos[2]
    new_gRNA_start=pos[4]-pos[2]
    new_gRNA_end=pos[5]-pos[2]
    row_colors=['black']*40 #row color
    row_colors[new_gRNA_start:(new_gRNA_end+1)]=['green']*(new_gRNA_end-new_gRNA_start+1)  #gRNA color
    row_colors[new_edit_start:(new_edit_end+1)]=['red']*(new_edit_end-new_edit_start+1)  #edit window color
    filename=k+'.pdf'
    count_fig+=1
    ax = plt.subplot(int(len(dict_exp)/2), 2, count_fig)
    figure(df.T,ax,row_colors)

plt.savefig('heatmap.png')

