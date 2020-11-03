
'''
Eric Fournier 2020-11-03
Note: need to run this script in nextstrain conda env

'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import matplotlib
import argparse
from Bio import pairwise2
import subprocess
import os
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
import subprocess



gisaid_seq = "sequences_small.fasta"

wuhan_ref = str(SeqIO.read('/data/Users/Eric/Covid19/reference.gb','genbank').seq)
wuhan_rec = SeqIO.read('/data/Users/Eric/Covid19/reference.gb','genbank')

wuhan_rec_sub_fasta = "WuhanPrimerRegion.fasta"
wuhan_rec_sub_list = []
wuhan_rec_sub_sarbeco = wuhan_rec[26200:26400]
wuhan_rec_sub_sarbeco.id = "wuhan_MN908947_sarbeco_region"
wuhan_rec_sub_sarbeco.description = ""
wuhan_rec_sub_lspq = wuhan_rec[28300:28500]
wuhan_rec_sub_lspq.id = "wuhan_MN908947_lspq_region"
wuhan_rec_sub_lspq.description = ""
wuhan_rec_sub_list.extend([wuhan_rec_sub_sarbeco,wuhan_rec_sub_lspq])


SeqIO.write(wuhan_rec_sub_list,wuhan_rec_sub_fasta,'fasta')


def FindPrimerBindingRange(primer):
    start = wuhan_ref.find(primer,0,len(wuhan_ref))
    end = start + len(primer)
    return([start,end])


forward_primer_sarbeco = "ACAGGTACGTTAATAGTTAATAGCGT"
rev_primer_sarbeco = str(Seq('ATATTGCAGCAGTACGCACACA', generic_dna).reverse_complement())
probe_primer_sarbeco = "ACACTAGCCATCCTTACTGCGCTTCG"
forward_primer_sarbeco_name = 'E_Sarbeco_F1'
rev_primer_sarbeco_name = 'E_Sarbeco_R2'
probe_primer_sarbeco_name = 'E_Sarbeco_P1'

forward_primer_lspq = "AACCAGAATGGAGAACGCAGTG"
rev_primer_lspq = str(Seq('CGGTGAACCAAGACGCAGTATTAT', generic_dna).reverse_complement())
probe_primer_lspq = "CGATCAAAACAACGTCGGCCCCAAGGTTTAC"
forward_primer_lspq_name = 'WuhanCoVNf'
rev_primer_lspq_name = 'WuhanCoVNr'
probe_primer_lspq_name = 'CoVNp'

forward_primer_sarbeco_range = FindPrimerBindingRange(forward_primer_sarbeco)
rev_primer_sarbeco_range = FindPrimerBindingRange(rev_primer_sarbeco)
probe_primer_sarbeco_range = FindPrimerBindingRange(probe_primer_sarbeco)
forward_primer_lspq_range = FindPrimerBindingRange(forward_primer_lspq)
rev_primer_lspq_range = FindPrimerBindingRange(rev_primer_lspq)
probe_primer_lspq_range = FindPrimerBindingRange(probe_primer_lspq)


def RunCmd(cmd):
    env = os.environ.copy()
    shellexec = ['env','bash']
    shargs = ['-c', "set -euo pipefail; " + cmd]

    try:
        subprocess.check_output(shellexec + shargs,shell = False, stderr = subprocess.STDOUT,env = env)
    except subprocess.CalledProcessError as error:
        print_error(
            "{out}\nshell exited {rc} when running: {cmd}{extra}",
            out = error.output,
            rc  = error.returncode,
            cmd = cmd,
            extra = "\nAre you sure this program is installed?" if error.returncode==127 else "",
        )
        if raise_errors:
            raise
        else:
            return False
    else:
        return True

seq_primer_dict = dict()


for rec in SeqIO.parse(gisaid_seq,'fasta'):
    print("Work on ",rec.id)
    rec_list = []
    rec_list.append(wuhan_rec)
    rec_list.append(rec)
    SeqIO.write(rec_list,"nonalign.fasta","fasta")
    align_cmd="mafft --reorder --anysymbol --nomemsave --adjustdirection --thread 40 nonalign.fasta > align.fasta"
    success = RunCmd(align_cmd)
    my_align = AlignIO.read("align.fasta",'fasta')
    #print("Length %i" % my_align.get_alignment_length())

    for align_rec in my_align:
        
        if align_rec.id != wuhan_rec.id:
            seq_primer_dict[align_rec.id] = dict()

            forward_primer_sarbeco_rec = SeqRecord(align_rec.seq[forward_primer_sarbeco_range[0]:forward_primer_sarbeco_range[1]])
            forward_primer_sarbeco_rec.id = align_rec.id + "_" + forward_primer_sarbeco_name
            forward_primer_sarbeco_rec.description =  ""
            seq_primer_dict[align_rec.id][forward_primer_sarbeco_name] = forward_primer_sarbeco_rec

            rev_primer_sarbeco_rec = SeqRecord(align_rec.seq[rev_primer_sarbeco_range[0]:rev_primer_sarbeco_range[1]])
            rev_primer_sarbeco_rec.id = align_rec.id + "_" + rev_primer_sarbeco_name
            rev_primer_sarbeco_rec.description =  ""
            seq_primer_dict[align_rec.id][rev_primer_sarbeco_name] = rev_primer_sarbeco_rec

            probe_primer_sarbeco_rec = SeqRecord(align_rec.seq[probe_primer_sarbeco_range[0]:probe_primer_sarbeco_range[1]])
            probe_primer_sarbeco_rec.id = align_rec.id + "_" + probe_primer_sarbeco_name
            probe_primer_sarbeco_rec.description =  ""
            seq_primer_dict[align_rec.id][probe_primer_sarbeco_name] = probe_primer_sarbeco_rec

            forward_primer_lspq_rec = SeqRecord(align_rec.seq[forward_primer_lspq_range[0]:forward_primer_lspq_range[1]])
            forward_primer_lspq_rec.id = align_rec.id + "_" + forward_primer_lspq_name
            forward_primer_lspq_rec.description =  ""
            seq_primer_dict[align_rec.id][forward_primer_lspq_name] = forward_primer_lspq_rec

            rev_primer_lspq_rec = SeqRecord(align_rec.seq[rev_primer_lspq_range[0]:rev_primer_lspq_range[1]])
            rev_primer_lspq_rec.id = align_rec.id + "_" + rev_primer_lspq_name
            rev_primer_lspq_rec.description =  ""
            seq_primer_dict[align_rec.id][rev_primer_lspq_name] = rev_primer_lspq_rec

            probe_primer_lspq_rec = SeqRecord(align_rec.seq[probe_primer_lspq_range[0]:probe_primer_lspq_range[1]])
            probe_primer_lspq_rec.id = align_rec.id + "_" + probe_primer_lspq_name
            probe_primer_lspq_rec.description =  ""
            seq_primer_dict[align_rec.id][probe_primer_lspq_name] = probe_primer_lspq_rec

primer_rec_list = []
primer_fasta = "primer.fasta"

def CreateFastaPrimer():
    for seq_id,primer_dict in seq_primer_dict.items():
        #print(" seq_id ", seq_id)
        #print(" primer ",primer_dict)
        for primer_name, primer_rec in primer_dict.items():
            primer_rec_list.append(primer_rec)

CreateFastaPrimer()

SeqIO.write(primer_rec_list,primer_fasta,'fasta')

rc = subprocess.call("/data/PROJETS/Covid19_NextStrainBuilds/TestPrimerBind_20201102/MapPrimerToGISAID.sh")
#print(rc)


exit(0)




