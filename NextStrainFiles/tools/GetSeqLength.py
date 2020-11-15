from Bio import AlignIO
import argparse


parser = argparse.ArgumentParser(description="")
parser.add_argument("--alignment",required=True,help="FASTA alignment")
args = parser.parse_args()


out_file = "seq_length.txt"

out_file_handle = open(out_file,'w')


for rec in AlignIO.read(args.alignment,'fasta'):
    out_file_handle.write(rec.id + " : " + str(len(rec.seq)) + " b\n")
    
out_file_handle.close()
