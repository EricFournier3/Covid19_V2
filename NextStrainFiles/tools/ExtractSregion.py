from Bio import AlignIO
from Bio import SeqIO
import argparse
import math

parser = argparse.ArgumentParser(description="")
parser.add_argument("--alignment",required=True,help="FASTA alignment")
parser.add_argument("--output-s-region",required=True,help="output s region aligned")
args = parser.parse_args()

rec_s_list = []

for rec in AlignIO.read(args.alignment,'fasta'):
    rec_s = rec[21562:25384]
    #print(rec_s)
    rec_s_list.append(rec_s)
    #print("S Length ",len(rec_s.seq))

SeqIO.write(rec_s_list,args.output_s_region,'fasta')

nb_total = len(rec_s_list)
print("TOTAL ",nb_total)

nb_var = 0

for rec in AlignIO.read(args.output_s_region,'fasta'):
    del_region = rec.seq[202:208]
    if 'N' in del_region:
        nb_var += 1

print("NB VAR ",nb_var)

frac = format((float(nb_var) / float(nb_total)) * 100,'.2f')
print("FRAC ",frac)

