import os
import re
import glob
import zipfile
import argparse
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import defaultdict

parser = argparse.ArgumentParser(description="Plot BE histogram...")
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

#Calculate substitution frequencies for each position, reads with InDels or Ns are excluded first.
def get_all_substitution_vectors(zip_file,pos):
   # print('Reading {}'.format(zip_file))
    z = zipfile.ZipFile(zip_file)
    zf=z.open("Alleles_frequency_table.txt")
    df=pd.read_csv(zf,sep="\t")
    w_start=int(pos[2]) #wide_window start for reference sequence
    w_end=int(pos[3]) #wide_window end for reference sequence
    ref=df['Reference_Sequence'][0][w_start:(w_end+1)] #wild type sequence, the first line
    df['ref_positions']=df['ref_positions'].apply(arrStr_to_arr)
    frequency_table=[]
    for i in range(len(ref)):
        frequency_table.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})
    for idx, row in df.iterrows():
        row_aligned_w_start=row['ref_positions'].index(w_start) #get wide window start for each aligned sequence; the wide window might be shifted if there are insertions before it
        row_aligned_w_end=row_aligned_w_start+39   #40bp window on aligned sequence
        row_align_seq=row['Aligned_Sequence'][row_aligned_w_start:(row_aligned_w_end+1)]
        row_ref_seq=row['Reference_Sequence'][row_aligned_w_start:(row_aligned_w_end+1)]
        if "-" in row_align_seq or "-" in row_ref_seq: #exclude reads with InDels
            continue
        if "N" in row_align_seq: #exclude reads with Ns
            continue
        for position,base in enumerate(row_align_seq):
            frequency_table[position][base] += row['#Reads']
    substitution_vectors={}
    nucs = ['A', 'T', 'C', 'G']
    for nuc in nucs:
        substitution_vectors[nuc]=[]
    for idx in range(len(ref)):
        for nuc in nucs:
            if nuc != ref[idx]:
                substitution_vectors[nuc].append(frequency_table[idx][nuc])
            else:
                substitution_vectors[nuc].append(0)
    for nuc in nucs:
            substitution_vectors[nuc] = np.array(substitution_vectors[nuc])
                ###Add substitutions of all line together!
    nLines = df['#Reads'].sum()
    return substitution_vectors, nLines,ref ##ref sequence 

def get_nuc_color(nuc, alpha):
    get_color=lambda x, y, z: (x/255.0, y/255.0, z/255.0, alpha)
    if nuc == "A":
        return get_color(127, 201, 127)
    elif nuc == "T":
        return get_color(190, 174, 212)
    elif nuc == "C":
        return get_color(253, 192, 134)
    elif nuc == "G":
        return get_color(255, 255, 153)

def get_color_lookup(nucs, alpha):
    colorLookup = {}
    for nuc in nucs:
        colorLookup[nuc] = get_nuc_color(nuc, alpha)
    return colorLookup

def plot_subs_across_wide_window(ref_len, ref_count, all_substitution_base_vectors, fig_filename, edit_win_start, edit_win_end, max_percent, ref_seq):
    '''
    plots substitutions across the reference sequece - each position on the x axis reprsents a nucleotide in the reference
    bars at each x posion show the number of times the reference nucleotide was substituted for another reference
    '''
    #if max_percent<0.1:
     #   fig, ax = plt.subplots(figsize=(10, 3))
    #else:
    #    fig, ax = plt.subplots(figsize=(10, 6))
    #fig, ax = plt.subplots(figsize=(10, 3)) #individual figure size
    ind = np.arange(ref_len+1)
    ind = np.delete(ind, 0)
    alph = ['A', 'C', 'G', 'T']
    bar_width = 0.3
    color_lookup = get_color_lookup(alph, alpha=1)
    pA = ax.bar(ind, all_substitution_base_vectors['A'], color=color_lookup['A'], width = bar_width)
    pC = ax.bar(ind, all_substitution_base_vectors['C'], color=color_lookup['C'], bottom=all_substitution_base_vectors['A'], width = bar_width)
    pG = ax.bar(ind, all_substitution_base_vectors["G"], color=color_lookup['G'], bottom=all_substitution_base_vectors["A"]+all_substitution_base_vectors["C"], width = bar_width)
    pT = ax.bar(ind, all_substitution_base_vectors["T"], color=color_lookup['T'], bottom=all_substitution_base_vectors["A"]+all_substitution_base_vectors["C"]+all_substitution_base_vectors["G"], width = bar_width)
    tots = all_substitution_base_vectors["A"]+all_substitution_base_vectors["C"]+all_substitution_base_vectors["G"]+all_substitution_base_vectors["T"]
    #y_max = max(15, (max(max(tots), 1))*1.1)
    y_max = np.round(ref_count*max_percent*1.1) #set y_max
    ax.set_ylim(0, y_max)
    ax.set_xlim([0, ref_len+1]) ##keep barplot inside the frame

    ref_seq_list = [*ref_seq] #ref_seq to list
    ax.set_xticks(ind)
    ax.set_xticklabels(ref_seq_list)
    
    legend_patches = [pA[0], pC[0], pG[0], pT[0]]
    legend_labels = ['A', 'C', 'G', 'T']
    
    p = matplotlib.patches.Rectangle((edit_win_start+0.5, 0), 1+(edit_win_end-edit_win_start), y_max, facecolor=(0, 0, 0, 0.05), edgecolor=(0, 0, 0, 0.25), linestyle=(0, (5, 2)), linewidth=1, label='Editing window')
    ax.add_patch(p)
    edit_win_patch = patches.Patch(fill=None, facecolor=(0, 0, 0, 0.05), edgecolor=(0, 0, 0, 0.25), linestyle=(0, (5, 2)), linewidth=1, label='Editing window')
    legend_patches.append(edit_win_patch)
    legend_labels.append('Editing window')
    ax.set_title(re.sub("_histogram.pdf","",fig_filename), fontsize=12)
    #ax.set_ylabel('% of total bases (Number of substitutions)', fontsize=18)
    #ax.set_xlabel('Quantification window position', fontsize=18)
    ####lgd = ax.legend(handles=legend_patches, labels=legend_labels, loc='upper right', ncol=2, fancybox=True, shadow=True) #remove legend
    #lgd = ax.legend(handles=legend_patches, labels=legend_labels, loc='upper right', bbox_to_anchor=(0.5, 0.5), ncol=2, fancybox=True, shadow=True)
    y_label_values= np.round(np.linspace(0, max(y_max, min(max(tots), max(ax.get_yticks()))), 6))
    ax.set_yticks(y_label_values)
    ax.set_yticklabels(
        [
            '%.1f%% (%d)' % (n_reads / ref_count * 100, n_reads)
            for n_reads in y_label_values
        ],
    )
    ax.tick_params(left=True, bottom=True)
    return ax
    #fig.savefig(fig_filename, bbox_extra_artists=(lgd,), bbox_inches='tight') #save separately
    #plt.close(fig) ##save separately

path=os.getcwd()
df_info = pd.read_excel(args.info, sheet_name="info",header=0)
dict_exp={}
for index,row in df_info.iterrows():
    #exp=row['GE_Tool']+"_"+row['gRNA']+"_"+row['Condition']+"_"+str(row['NGS_Code']) ## single replicate
    exp=row['GE_Tool']+"_"+row['gRNA']+"_"+row['Condition']  ## average of replicates
    infos=row['WinEdit']+"-"+row['Win40_gCenter']+"-"+row['gRNA_win']+"-"+row['gRNA_Strand']
    dict_exp[exp]={}
    dict_exp[exp]=infos #editWin-wideWin-gRNAwin-strand

dict_fig_info={}
for k,v in dict_exp.items():
    res=glob.glob('{}/result/CRISPResso_on_{}_*/Alleles_frequency_table.zip'.format(path,k.replace(".","_"))) #get one zip file for each replicate
    pos=v.split("-")[0:-1]
    pos=[ int(x) for x in pos ]
    #all_substitution_vectors,total_reads=get_all_substitution_vectors(res[0],pos) # single replicate
    all_substitution_vectors, total_reads, ref_seq=get_all_substitution_vectors(res[0],pos) #ref_seq
    for i in range(len(res)):
        if i==0:
            continue
        else:
            for nuc in ['A', 'C', 'G', 'T']:
                all_substitution_vectors[nuc] += (get_all_substitution_vectors(res[i],pos))[0][nuc]
            total_reads = total_reads + (get_all_substitution_vectors(res[i],pos))[1]
    for nuc in ['A', 'C', 'G', 'T']:
        all_substitution_vectors[nuc] = (all_substitution_vectors[nuc]/len(res)).round(0)
    total_reads = (total_reads/len(res)).round(0)
    strand=v.split("-")[-1].strip()
    complement_substition_vectors={}
    for nuc in ['A', 'C', 'G', 'T']:
        complement_substition_vectors[nuc]=[]

    if 'reverse' in strand: #use corresponding complementary sequence if it is a reverse strand, so it is always A>G for ABE, C>T for CBE
        ref_seq = get_complementary_seq(ref_seq) #use complement sequence
        for key in all_substitution_vectors.keys():
            key_complement=get_complementary_seq(key)
     #       print('complement key %s, key %s' % (key_complement,key))
     #       print(all_substitution_vectors[key])
            complement_substition_vectors[key_complement] = all_substitution_vectors[key]
    else:
        complement_substition_vectors = all_substitution_vectors 
    
    new_edit_start=pos[0]-pos[2] #editint window position relative to the first base of wide window
    new_edit_end=pos[1]-pos[2]
    ref_len = pos[3] - pos[2] + 1
    filename=k+'_histogram.pdf'
    dict_fig_info[k]={}
    dict_fig_info[k]=[ref_len, total_reads, complement_substition_vectors, filename, new_edit_start, new_edit_end, ref_seq]

plt.figure(figsize=(10*2, 4*len(dict_fig_info)/2))
plt.tight_layout(h_pad = 2)
count_fig = 0
for k,v in sorted(dict_fig_info.items()): #sort keys
    sub_vectors = v[2] #get substitution vector for this condition
    max_sub_percent = np.round(max(sub_vectors['A'] + sub_vectors['C'] + sub_vectors['G'] + sub_vectors['T'])/v[1],3)
    if k.endswith('N'):
        k_Y = re.sub("N$","Y",k)
        sub_vectors_Y = dict_fig_info[k_Y][2]
        total_reads_Y = dict_fig_info[k_Y][1]
        max_sub_percent_Y = np.round(max(sub_vectors_Y['A'] + sub_vectors_Y['C'] + sub_vectors_Y['G'] + sub_vectors_Y['T'])/total_reads_Y,3)
        max_percent = max(max_sub_percent, max_sub_percent_Y)
    else:
        k_N = re.sub("Y$","N",k)
        sub_vectors_N = dict_fig_info[k_N][2]
        total_reads_N = dict_fig_info[k_N][1]
        max_sub_percent_N = np.round(max(sub_vectors_N['A'] + sub_vectors_N['C'] + sub_vectors_N['G'] + sub_vectors_N['T'])/total_reads_N, 3)
        max_percent = max(max_sub_percent, max_sub_percent_N)
    count_fig+=1
    ax = plt.subplot(int(len(dict_fig_info)/2), 2, count_fig)
    #ax.tight_layout(h_pad = 2)
    plot_subs_across_wide_window(v[0], v[1], v[2], v[3], v[4], v[5], max_percent, v[6])

plt.subplots_adjust(top = 0.85)
plt.savefig('all_histograms.pdf')

