"""
Mask initial bases from alignment FASTA
"""
import argparse
import Bio
import Bio.SeqIO
from Bio.Seq import Seq
from Bio import SeqIO
import os

def check_cobas_mutation(fasta_seq):
    check_dir = {}
    cobas_mutant = {}
    cobas_undefined = {}

    res_file = os.path.join(os.path.dirname(fasta_seq),"Cobas_mutant.txt")
    res_file_undefined = os.path.join(os.path.dirname(fasta_seq),"Cobas_undefined.txt")

    res_file_handle = open(res_file,'w')
    res_file_undefined_handle = open(res_file_undefined,'w')

    for record in SeqIO.parse(fasta_seq,'fasta'):
        target_nuc_at_26340 = str(record.seq[26339])
        check_dir[record.id] = target_nuc_at_26340

        if target_nuc_at_26340.upper() not in ['C','N']:
            cobas_mutant[record.id] = target_nuc_at_26340

        if target_nuc_at_26340.upper() == 'N':
            cobas_undefined[record.id] = target_nuc_at_26340
        
    for spec_id,nuc in cobas_mutant.items():
        res_file_handle.write(spec_id + "\t" + nuc + "\n")

    for spec_id,nuc in cobas_undefined.items():
        res_file_undefined_handle.write(spec_id + "\t" + nuc + "\n")

    res_file_handle.close()
    res_file_undefined_handle.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Mask initial bases from alignment FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", required=True, help="FASTA file of alignment")
    parser.add_argument("--mask-from-beginning", type = int, required=True, help="number of bases to mask from start")
    parser.add_argument("--mask-from-end", type = int, help="number of bases to mask from end")
    parser.add_argument("--mask-sites", nargs='+', type = int,  help="list of sites to mask")
    parser.add_argument("--output", required=True, help="FASTA file of output alignment")
    args = parser.parse_args()

    check_cobas_mutation(args.alignment)
    



    being_length = 0
    if args.mask_from_beginning:
        begin_length = args.mask_from_beginning
    end_length = 0
    if args.mask_from_end:
        end_length = args.mask_from_end

    with open(args.output, 'w') as outfile:
        for record in Bio.SeqIO.parse(args.alignment, 'fasta'):
            seq = str(record.seq)
            start = "N" * begin_length
            middle = seq[begin_length:-end_length]
            end = "N" * end_length
            seq_list = list(start + middle + end)
            if args.mask_sites:
                for site in args.mask_sites:
                    seq_list[site-1] = "N"
            record.seq = Seq("".join(seq_list))
            Bio.SeqIO.write(record, outfile, 'fasta')
