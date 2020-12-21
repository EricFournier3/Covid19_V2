from Bio import AlignIO
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("--alignment",required=True,help="FASTA alignment")
parser.add_argument("--output-s-region",required=True,help="output s region aligned")
args = parser.parse_args()

rec_s_list = []

for rec in AlignIO.read(args.alignment,'fasta'):
    rec_s = rec[21563:25385]
    #print(rec_s)
    rec_s_list.append(rec_s)

SeqIO.write(rec_s_list,args.output_s_region,'fasta')
