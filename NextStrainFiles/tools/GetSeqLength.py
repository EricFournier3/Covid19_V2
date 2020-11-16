from Bio import AlignIO
import argparse


parser = argparse.ArgumentParser(description="")
parser.add_argument("--alignment",required=True,help="FASTA alignment")
parser.add_argument("--output-seq-length",required=True,help="output seq length file")
args = parser.parse_args()

out_file_handle = open(args.output_seq_length,'w')


for rec in AlignIO.read(args.alignment,'fasta'):
    out_file_handle.write(rec.id + " : " + str(len(rec.seq)) + " bases\n")
    
out_file_handle.close()

