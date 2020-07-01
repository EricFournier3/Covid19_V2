import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import matplotlib
import argparse
"""
Eric Fournier 2020-06-23

"""

parser = argparse.ArgumentParser(
description="Modify tree based on exposure",
formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument('--input', type=str, required=True, help="input alignment")
parser.add_argument('--out-sarbeco', type=str, metavar="TSV", required=True, help="output primer sarbeco profile")
parser.add_argument('--out-lspq', type=str, metavar="TSV", required=True, help="output primer lspq profile")
args = parser.parse_args()


fasta_file = args.input

wuhan_ref = str(SeqIO.read('/data/Users/Eric/Covid19/reference.gb','genbank').seq)


forward_primer_sarbeco = "ACAGGTACGTTAATAGTTAATAGCGT"
rev_primer_sarbeco = str(Seq('ATATTGCAGCAGTACGCACACA', generic_dna).reverse_complement())
probe_primer_sarbeco = "ACACTAGCCATCCTTACTGCGCTTCG"
primer_sarbeco_fig = args.out_sarbeco
forwad_primer_sarbeco_name = 'E_Sarbeco_F1'
rev_primer_sarbeco_name = 'E_Sarbeco_R2'
probe_primer_sarbeco_name = 'E_Sarbeco_P1'

forward_primer_lspq = "AACCAGAATGGAGAACGCAGTG"
rev_primer_lspq = str(Seq('CGGTGAACCAAGACGCAGTATTAT', generic_dna).reverse_complement())
probe_primer_lspq = "CGATCAAAACAACGTCGGCCCCAAGGTTTAC"
primer_lspq_fig = args.out_lspq
forwad_primer_lspq_name = 'WuhanCoVNf'
rev_primer_lspq_name = 'WuhanCoVNr'
probe_primer_lspq_name = 'CoVNp'

nuc = {'A':0,'C':1,'G':2,'T':3,'N':4}
nuc_colors = {'A':'#0000FF','C':'#FF0000','G':'#008000','T':'#FFFF00','N':'#FFC0CB'}


def PrepareMinorFreqNuc(nuc_freq):

    nuc_freq_arr = np.array(nuc_freq)

    max_by_col = np.amax(nuc_freq_arr,axis=0)
    min_freq_nuc = np.zeros(nuc_freq_arr.shape)

    for i in range(0,nuc_freq_arr.shape[0]):
        for j in range(0,nuc_freq_arr.shape[1]):
            val = nuc_freq_arr[i][j]
            if (val != max_by_col[j]) and (val != 0):
                min_freq_nuc[i][j] = 1
            else:
                min_freq_nuc[i][j] = 0

    min_freq_nuc = pd.DataFrame(data=min_freq_nuc)
    return(min_freq_nuc)
    
                
def FindPrimerBindingRange(primer):
    start = wuhan_ref.find(primer,0,len(wuhan_ref))
    end = start + len(primer)
    return((start,end))

def UpdateNucFrequency(_nuc,_index,nuc_frequency):
    if _nuc.upper() not in ['A','C','G','T']:
        _nuc = 'N'

    line = nuc[_nuc.upper()]
    col = _index

    nuc_frequency[line][col] += 1


def BuildNucFrequency():
    for rec in SeqIO.parse(fasta_file,'fasta'):
        site_index = 0

        for pos_site in range(*forward_primer_sarbeco_range):
            UpdateNucFrequency(rec.seq[pos_site],site_index,nuc_freq_by_site_forward_primer_sarbeco)
            site_index += 1

        site_index = 0

        for pos_site in range(*rev_primer_sarbeco_range):
            UpdateNucFrequency(rec.seq[pos_site],site_index,nuc_freq_by_site_rev_primer_sarbeco)
            site_index += 1

        site_index = 0

        for pos_site in range(*probe_primer_sarbeco_range):
            UpdateNucFrequency(rec.seq[pos_site],site_index,nuc_freq_by_site_probe_primer_sarbeco)
            site_index += 1

        site_index = 0

        for pos_site in range(*forward_primer_lspq_range):
            UpdateNucFrequency(rec.seq[pos_site],site_index,nuc_freq_by_site_forward_primer_lspq)
            site_index += 1

        site_index = 0

        for pos_site in range(*rev_primer_lspq_range):
            UpdateNucFrequency(rec.seq[pos_site],site_index,nuc_freq_by_site_rev_primer_lspq)
            site_index += 1

        site_index = 0

        for pos_site in range(*probe_primer_lspq_range):
            UpdateNucFrequency(rec.seq[pos_site],site_index,nuc_freq_by_site_probe_primer_lspq)
            site_index += 1

matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Helvetica Neue')
matplotlib.rc('text', usetex='false')
matplotlib.rcParams.update({'font.size': 10})



def PlotNucFrequency(primer,ax,nuc_freq,primer_name):

    base_title = "Nucleotids distribution for primer "

    table_row_name = [y[0] for y in sorted(nuc.items(), key=lambda x:x[1])]
    table_col_name = [primer[i] for i in range(0,len(primer))]

    nuc_freq_text = []
    y_offset = np.zeros(len(table_col_name))
    colors = [nuc_colors['A'],nuc_colors['C'],nuc_colors['G'],nuc_colors['T'],nuc_colors['N']]
    
    bar_index = np.arange(len(table_col_name)) + 0.25
    bar_width = 0.3
    spare_width = (1 - bar_width * 2) / 2

    ax.set_xlim(-spare_width, len(table_col_name) - spare_width)

    for row in range(len(table_row_name)):

        dat = nuc_freq[row]
        
        plt.bar(bar_index,dat,bar_width,bottom=y_offset,color=colors[row])
        y_offset = y_offset + dat
        nuc_freq_text.append([x for x in dat])

    minor_freq_nuc = PrepareMinorFreqNuc(nuc_freq_text)

    the_table = plt.table(cellText=nuc_freq_text,rowLabels=table_row_name,colLabels=table_col_name,rowColours=colors,cellLoc='center',loc='bottom',bbox=[0,-0.95,1,0.95],cellColours=plt.cm.YlOrRd(minor_freq_nuc))

    
    #voir https://stackoverflow.com/questions/27972524/matplotlib-row-heights-table-property et https://stackoverflow.com/questions/13921338/matplotlib-table-formatting-change-width-of-row-label-cells
    cellDict = the_table.get_celld()
    #header line
    for column in range(0,len(table_col_name)):
        cellDict[(0,column)].set_height(.8)

    #lines A, C, G, T and N
    for column in range(-1,len(table_col_name)):
        for line in range(1,len(table_row_name)+1):
            cellDict[(line,column)].set_height(.75)


    the_table.scale(1,3.5)
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(10)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    
    plt.title(base_title + primer_name)
            
forward_primer_sarbeco_range = FindPrimerBindingRange(forward_primer_sarbeco)
rev_primer_sarbeco_range = FindPrimerBindingRange(rev_primer_sarbeco)
probe_primer_sarbeco_range = FindPrimerBindingRange(probe_primer_sarbeco)
forward_primer_lspq_range = FindPrimerBindingRange(forward_primer_lspq)
rev_primer_lspq_range = FindPrimerBindingRange(rev_primer_lspq)
probe_primer_lspq_range = FindPrimerBindingRange(probe_primer_lspq)

'''
print("forward_primer_sarbeco_range: ",forward_primer_sarbeco_range)
print("rev_primer_sarbeco_range: ",rev_primer_sarbeco_range)
print("probe_primer_sarbeco_range: ",probe_primer_sarbeco_range)
print("forward_primer_lspq_range: ",forward_primer_lspq_range)
print("rev_primer_lspq_range: ",rev_primer_lspq_range)
print("probe_primer_lspq_range: ",probe_primer_lspq_range)
'''

forward_primer_sarbeco_range_len = len(range(*forward_primer_sarbeco_range))
rev_primer_sarbeco_range_len = len(range(*rev_primer_sarbeco_range))
probe_primer_sarbeco_range_len = len(range(*probe_primer_sarbeco_range))
forward_primer_lspq_range_len = len(range(*forward_primer_lspq_range))
rev_primer_lspq_range_len = len(range(*rev_primer_lspq_range))
probe_primer_lspq_range_len = len(range(*probe_primer_lspq_range))

nuc_freq_by_site_forward_primer_sarbeco = np.zeros(len(nuc) * forward_primer_sarbeco_range_len).reshape(len(nuc),forward_primer_sarbeco_range_len)
nuc_freq_by_site_rev_primer_sarbeco = np.zeros(len(nuc) * rev_primer_sarbeco_range_len).reshape(len(nuc),rev_primer_sarbeco_range_len)
nuc_freq_by_site_probe_primer_sarbeco = np.zeros(len(nuc) * probe_primer_sarbeco_range_len).reshape(len(nuc),probe_primer_sarbeco_range_len)
nuc_freq_by_site_forward_primer_lspq = np.zeros(len(nuc) * forward_primer_lspq_range_len).reshape(len(nuc),forward_primer_lspq_range_len)
nuc_freq_by_site_rev_primer_lspq = np.zeros(len(nuc) * rev_primer_lspq_range_len).reshape(len(nuc),rev_primer_lspq_range_len)
nuc_freq_by_site_probe_primer_lspq = np.zeros(len(nuc) * probe_primer_lspq_range_len).reshape(len(nuc),probe_primer_lspq_range_len)

BuildNucFrequency()


fig_sarbeco = plt.figure()
#fig_sarbeco.set_size_inches(20.0,7.5)
fig_sarbeco.set_size_inches(20.0,12.5)

forward_sarbeco_ax = fig_sarbeco.add_subplot(311)
PlotNucFrequency(forward_primer_sarbeco,forward_sarbeco_ax,nuc_freq_by_site_forward_primer_sarbeco,forwad_primer_sarbeco_name)

rev_sarbeco_ax = fig_sarbeco.add_subplot(312)
PlotNucFrequency(rev_primer_sarbeco,rev_sarbeco_ax,nuc_freq_by_site_rev_primer_sarbeco,rev_primer_sarbeco_name)

probe_sarbeco_ax = fig_sarbeco.add_subplot(313)
PlotNucFrequency(probe_primer_sarbeco,probe_sarbeco_ax,nuc_freq_by_site_probe_primer_sarbeco,probe_primer_sarbeco_name)

fig_sarbeco.tight_layout(pad=0.2)
plt.subplots_adjust(left=0.05, bottom=0.2)
plt.savefig(primer_sarbeco_fig)

#plt.show()
plt.close()

fig_lspq = plt.figure()
fig_lspq.set_size_inches(20.0,12.5)

forward_lspq_ax = fig_lspq.add_subplot(311)
PlotNucFrequency(forward_primer_lspq,forward_lspq_ax,nuc_freq_by_site_forward_primer_lspq,forwad_primer_lspq_name)

rev_lspq_ax = fig_lspq.add_subplot(312)
PlotNucFrequency(rev_primer_lspq,rev_lspq_ax,nuc_freq_by_site_rev_primer_lspq,rev_primer_lspq_name)

probe_lspq_ax = fig_lspq.add_subplot(313)
PlotNucFrequency(probe_primer_lspq,probe_lspq_ax,nuc_freq_by_site_probe_primer_lspq,probe_primer_lspq_name)

fig_lspq.tight_layout(pad=0.2)
plt.subplots_adjust(left=0.05, bottom=0.2)
plt.savefig(primer_lspq_fig)

#plt.show()
plt.close()
