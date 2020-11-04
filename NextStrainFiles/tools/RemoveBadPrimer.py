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


fasta_in = "/data/PROJETS/Covid19_NextStrainBuilds/TestPrimerBind_20201102/primer.fasta"

out_primer = "/data/PROJETS/Covid19_NextStrainBuilds/TestPrimerBind_20201102/good_primer.fasta"


good_primer_rec_list = []

good_character = ['a','c','g','t','n','A','C','G','T','N']
bad_character = '-'

for rec in SeqIO.parse(fasta_in,'fasta'):
    seq = str(rec.seq)
    #print(seq)
    '''
    test = [c in good_character for c in seq]
    if not all(test):
        print("Bad primer seq ", seq)
    '''
    test = bad_character in seq
    if test:
        #print("Bad primer seq ", seq)
        pass
    else:
        good_primer_rec_list.append(rec)

SeqIO.write(good_primer_rec_list,out_primer,'fasta')

