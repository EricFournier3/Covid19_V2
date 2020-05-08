"""
Eric Fournier 20200415
"""

import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
import logging

basedir = "/data/PROJETS/COVID_19/TestMetadata"

ref_ids = ["Wuhan-Hu-1/2019","Wuhan/WH01/2019"]

metadata = os.path.join(basedir,"metadata.tsv")
metadata_handle = open(metadata)
metadata_handle.readline()

fake_seq_file = os.path.join(basedir,"sequences_fakeid.fasta")
new_seq_file = os.path.join(basedir,"sequences.fasta")

new_rec_list = []
new_rec_list_id =  []

logging.basicConfig(level = logging.DEBUG)

metadata_id = []

for line in metadata_handle:
    id = line.split("\t")[0]
    metadata_id.append(id)

metadata_handle.close()

inc = 0

for rec in SeqIO.parse(fake_seq_file,'fasta'):

    if(inc < len(metadata_id)):
        new_rec = SeqRecord(seq = rec.seq, id = metadata_id[inc], description = "" )
        new_rec_list.append(new_rec)
        inc += 1

SeqIO.write(new_rec_list,new_seq_file,"fasta")

logging.info("End")

