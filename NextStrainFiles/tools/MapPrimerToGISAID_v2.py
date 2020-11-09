
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
import sys
import re

          

align_in = "/data/PROJETS/Covid19_NextStrainBuilds/TestPrimerBind_20201102/gisaid_qc_aligned.fasta"

wuhan_ref = str(SeqIO.read('/data/Users/Eric/Covid19/reference.gb','genbank').seq)
wuhan_rec = SeqIO.read('/data/Users/Eric/Covid19/reference.gb','genbank')
wuhan_rec_sub_fasta = "/data/PROJETS/Covid19_NextStrainBuilds/TestPrimerBind_20201102/WuhanPrimerRegion.fasta"
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

seq_primer_dict = dict()

nb_rec_treated = 0


base_dir = "/data/PROJETS/Covid19_NextStrainBuilds/TestPrimerBind_20201102/"


for rec in AlignIO.read(align_in,'fasta'):
    #print("Work on ",rec.id)
    #print("Length is ",len(rec.seq))
   
    nb_rec_treated += 1
 
    sys.stdout.write("Work on record > %d\r"%nb_rec_treated)
    sys.stdout.flush()

    seq_primer_dict[rec.id] = dict()
    
    forward_primer_sarbeco_rec = SeqRecord(rec.seq[forward_primer_sarbeco_range[0]:forward_primer_sarbeco_range[1]])
    forward_primer_sarbeco_rec.id = rec.id + "_" + forward_primer_sarbeco_name
    forward_primer_sarbeco_rec.description =  ""
    seq_primer_dict[rec.id][forward_primer_sarbeco_name] = forward_primer_sarbeco_rec

    rev_primer_sarbeco_rec = SeqRecord(rec.seq[rev_primer_sarbeco_range[0]:rev_primer_sarbeco_range[1]])
    rev_primer_sarbeco_rec.id = rec.id + "_" + rev_primer_sarbeco_name
    rev_primer_sarbeco_rec.description =  ""
    seq_primer_dict[rec.id][rev_primer_sarbeco_name] = rev_primer_sarbeco_rec 

    probe_primer_sarbeco_rec = SeqRecord(rec.seq[probe_primer_sarbeco_range[0]:probe_primer_sarbeco_range[1]])
    probe_primer_sarbeco_rec.id = rec.id + "_" + probe_primer_sarbeco_name
    probe_primer_sarbeco_rec.description =  ""
    seq_primer_dict[rec.id][probe_primer_sarbeco_name] = probe_primer_sarbeco_rec

    forward_primer_lspq_rec = SeqRecord(rec.seq[forward_primer_lspq_range[0]:forward_primer_lspq_range[1]])
    forward_primer_lspq_rec.id = rec.id + "_" + forward_primer_lspq_name
    forward_primer_lspq_rec.description =  ""
    seq_primer_dict[rec.id][forward_primer_lspq_name] = forward_primer_lspq_rec

    rev_primer_lspq_rec = SeqRecord(rec.seq[rev_primer_lspq_range[0]:rev_primer_lspq_range[1]])
    rev_primer_lspq_rec.id = rec.id + "_" + rev_primer_lspq_name
    rev_primer_lspq_rec.description =  ""
    seq_primer_dict[rec.id][rev_primer_lspq_name] = rev_primer_lspq_rec

    probe_primer_lspq_rec = SeqRecord(rec.seq[probe_primer_lspq_range[0]:probe_primer_lspq_range[1]])
    probe_primer_lspq_rec.id = rec.id + "_" + probe_primer_lspq_name
    probe_primer_lspq_rec.description =  ""
    seq_primer_dict[rec.id][probe_primer_lspq_name] = probe_primer_lspq_rec


#print(seq_primer_dict)

primer_rec_list = []
primer_fasta = "/data/PROJETS/Covid19_NextStrainBuilds/TestPrimerBind_20201102/primer.fasta"

good_character = ['a','c','g','t','n','A','C','G','T','N']
bad_character = '-'

def CheckPrimerSeq(seq):
    test = bad_character in seq
    if test:
        return False
    else:
        return True

out_check_mismatch_sarbeco_r = os.path.join(base_dir,"Mismatches_SarbecoR.fasta")
out_check_mismatch_sarbeco_f = os.path.join(base_dir,"Mismatches_SarbecoF.fasta")
out_check_mismatch_lspq_r = os.path.join(base_dir,"Mismatches_LspqR.fasta")
out_check_mismatch_lspq_f = os.path.join(base_dir,"Mismatches_LspqF.fasta")

mismatch_sarbeco_r = []
mismatch_sarbeco_f = []
mismatch_lspq_r = []
mismatch_lspq_f = []

def CheckMismatchSarbecoR(primer_rec):
    #print("Sarbeco R ",primer_rec.id)
    #rev_primer_sarbeco
    first_nuc = str(primer_rec.seq)[0]    

    if (rev_primer_sarbeco[0] != first_nuc) and (first_nuc in ['A','C','G','T','a','c','g','t']):
        mismatch_sarbeco_r.append(primer_rec)

def CheckMismatchSarbecoF(primer_rec):
    #print("Sarbeco F ",primer_rec.id)
    last_nuc = str(primer_rec.seq)[-1]   

    if (forward_primer_sarbeco[-1] != last_nuc) and (last_nuc  in ['A','C','G','T','a','c','g','t']):
        mismatch_sarbeco_f.append(primer_rec)

def CheckMismatchLspqR(primer_rec):
    #print("Lspq R ",primer_rec.id)
    first_nuc = str(primer_rec.seq)[0]  

    if (rev_primer_lspq[0] != first_nuc) and (first_nuc in ['A','C','G','T','a','c','g','t']):
        mismatch_lspq_r.append(primer_rec)

def CheckMismatchLspqF(primer_rec):
    #print("Lspq F ",primer_rec.id)
    last_nuc = str(primer_rec.seq)[-1]   

    if (forward_primer_lspq[-1] != last_nuc) and (last_nuc  in ['A','C','G','T','a','c','g','t']):
        mismatch_lspq_f.append(primer_rec)

def CreateFastaPrimer():
    for seq_id,primer_dict in seq_primer_dict.items():
        for primer_name, primer_rec in primer_dict.items():
            if CheckPrimerSeq(str(primer_rec.seq)):
                if re.search(r'E_Sarbeco_R2$',primer_rec.id):
                    CheckMismatchSarbecoR(primer_rec)
                elif re.search(r'E_Sarbeco_F1$',primer_rec.id):
                    CheckMismatchSarbecoF(primer_rec)
                elif re.search(r'WuhanCoVNf$',primer_rec.id):
                    CheckMismatchLspqF(primer_rec)
                elif re.search(r'WuhanCoVNr$',primer_rec.id):
                    #print(primer_rec.id)
                    CheckMismatchLspqR(primer_rec)

                primer_rec_list.append(primer_rec)


CreateFastaPrimer()

SeqIO.write(mismatch_sarbeco_r,out_check_mismatch_sarbeco_r,'fasta')
#print("UN ",mismatch_sarbeco_r)
SeqIO.write(mismatch_sarbeco_f,out_check_mismatch_sarbeco_f,'fasta')
#print("deux ",mismatch_sarbeco_f)
SeqIO.write(mismatch_lspq_r,out_check_mismatch_lspq_r,'fasta')
#print("trois ",mismatch_lspq_r)
SeqIO.write(mismatch_lspq_f,out_check_mismatch_lspq_f,'fasta')
#print("quatre ",mismatch_lspq_f)


SeqIO.write(primer_rec_list,primer_fasta,'fasta')


step = "1"
rc = subprocess.check_call(["/data/PROJETS/Covid19_NextStrainBuilds/TestPrimerBind_20201102/MapPrimerToGISAID.sh", step])

sam_file_in = "/data/PROJETS/Covid19_NextStrainBuilds/TestPrimerBind_20201102/primer_map_temp.sam"
sam_file_out = "/data/PROJETS/Covid19_NextStrainBuilds/TestPrimerBind_20201102/primer_map.sam"


with open(sam_file_out,'w') as wf:
    with open(sam_file_in) as rf:
        for line in rf:
            l = line.split("\t")
            try:
                if re.search(r'E_Sarbeco_F1$',l[0]):
                    l[3] = "69"
                    l[5] = "26M"
                elif re.search(r'E_Sarbeco_P1$',l[0]):
                    l[3] = "132"
                    l[5] = "26M"
                elif re.search(r'E_Sarbeco_R2$',l[0]):
                    l[3] = "160"
                    l[5] = "22M"
                elif re.search(r'WuhanCoVNf$',l[0]):
                    l[3] = "52"
                    l[5] = "22M"
                elif re.search(r'CoVNp$',l[0]):
                    l[3] = "79"
                    l[5] = "31M"
                elif re.search(r'WuhanCoVNr$',l[0]):
                    l[3] = "113"
                    l[5] = "24M"
            except:
                pass

            l = "\t".join(l)
            wf.write(l)

step = "2"

rc = subprocess.check_call(["/data/PROJETS/Covid19_NextStrainBuilds/TestPrimerBind_20201102/MapPrimerToGISAID.sh", step])



exit(0)



