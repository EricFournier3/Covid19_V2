# coding=utf-8

"""
Eric Fournier 2020-12-08

"""

import datetime
import pandas as pd
import os
import numpy as np
import re
import sys
import logging
import gc
import argparse
import time
import glob
from Bio import SeqIO

beluga_server = "fournie1@beluga.computecanada.ca:/home/fournie1"
mnt_beluga_server = "/mnt/BelugaEric/"

def MountBelugaServer():
    logging.info("Try to mount Beluga")
    os.system("sudo umount " + mnt_beluga_server)
    os.system("sudo sshfs -o allow_other -o follow_symlinks {0} {1}".format(beluga_server,mnt_beluga_server))
    logging.info("Beluga mounted")

#MountBelugaServer()


compare_outdir = "/data/PROJETS/COVID-19_Beluga/Consensus/CompareFreeze1_OLDvsNEW_pipeline/"
new_consensus_outdir = os.path.join(compare_outdir,"NewConsensus")

illumina_reprocess_path_file = os.path.join(compare_outdir,"illumina_reprocess_path.txt")
mgi_reprocess_path_file = os.path.join(compare_outdir,"mgi_reprocess_path.txt")

illumina_reprocess_path_df = pd.read_table(illumina_reprocess_path_file)
illumina_reprocess_path_df.columns = ['PATH']

mgi_reprocess_path_df = pd.read_table(mgi_reprocess_path_file)
mgi_reprocess_path_df.columns = ['PATH']

base_dir_gisaid_sub = "/data/PROJETS/COVID-19_Beluga/Gisaid/FinalPublished/release1/"

base_dir_gisaid_sub_1 = os.path.join(base_dir_gisaid_sub,"20200520")
base_dir_gisaid_sub_2 = os.path.join(base_dir_gisaid_sub,"20200610")
base_dir_gisaid_sub_3 = os.path.join(base_dir_gisaid_sub,"20200914")

excel_gisaid_metadata_sub_1 = os.path.join(base_dir_gisaid_sub_1,"20200520_ncov19_metadata.xls")
excel_gisaid_metadata_sub_2 = os.path.join(base_dir_gisaid_sub_2,"20200610_ncov19_metadata.xls")
excel_gisaid_metadata_sub_3 = os.path.join(base_dir_gisaid_sub_3,"20200914_ncov19_metadata.xls")

excel_gisaid_metadata_sub_1_df = pd.read_excel(excel_gisaid_metadata_sub_1,sheet_name=0)
excel_gisaid_metadata_sub_2_df = pd.read_excel(excel_gisaid_metadata_sub_2,sheet_name=0)
excel_gisaid_metadata_sub_3_df = pd.read_excel(excel_gisaid_metadata_sub_3,sheet_name=0)


excel_gisaid_metadata_all_df_all = pd.concat([excel_gisaid_metadata_sub_1_df,excel_gisaid_metadata_sub_2_df,excel_gisaid_metadata_sub_3_df])
excel_gisaid_metadata_all_df_all =  excel_gisaid_metadata_all_df_all.drop(excel_gisaid_metadata_all_df_all.index[0])
#excel_gisaid_metadata_all_df_all covv_seq_technology
excel_gisaid_metadata_all_df_all = excel_gisaid_metadata_all_df_all[['covv_virus_name','covv_seq_technology']]
excel_gisaid_metadata_all_df_all =    excel_gisaid_metadata_all_df_all.reset_index(drop=True)
#print(excel_gisaid_metadata_all_df_all)

fasta_gisaid_sequence_sub_1 = os.path.join(base_dir_gisaid_sub_1,"all_sequences.fasta")
fasta_gisaid_sequence_sub_2 = os.path.join(base_dir_gisaid_sub_2,"all_sequences.fasta")
fasta_gisaid_sequence_sub_3 = os.path.join(base_dir_gisaid_sub_3,"all_sequences.fasta")


gisaid_rec_dict = {}


def GetTechno(seq_id):
    techno = excel_gisaid_metadata_all_df_all.loc[excel_gisaid_metadata_all_df_all['covv_virus_name'] == seq_id,'covv_seq_technology']
    #print("TECHNO IS ",(list(techno)[0]))
    return(list(techno)[0])


for gisaid_fasta in [fasta_gisaid_sequence_sub_1,fasta_gisaid_sequence_sub_2,fasta_gisaid_sequence_sub_3]:
    for rec in SeqIO.parse(gisaid_fasta,'fasta'):
        techno = GetTechno(rec.id)
        gisaid_rec_dict[rec.id] = techno

print("len(gisaid_rec_dict) : ",len(gisaid_rec_dict))

new_rej_consensus = []
new_multiple_consensus = []
new_keeped_consensus_list = []

gisaid_consensus_not_found_in_new = []
nb_gisaid_nanapore_consensus = 0

for seq_id in  gisaid_rec_dict.keys():
    beluga_seq_id = re.search(r'^hCoV-19/Canada/Qc-(\S+)/\d{4}',seq_id).group(1)
    techno = gisaid_rec_dict[seq_id]
    if techno == "Illumina_NexteraFlex":
        path = illumina_reprocess_path_df.loc[illumina_reprocess_path_df['PATH'].str.contains(beluga_seq_id),'PATH']
    elif techno == "MGI_CleanPlex":
        path = mgi_reprocess_path_df.loc[mgi_reprocess_path_df['PATH'].str.contains(beluga_seq_id),'PATH']
    else:
        nb_gisaid_nanapore_consensus += 1
        continue

    path_list = list(path)

   
    if len(path_list) == 0:
        gisaid_consensus_not_found_in_new += 1
    elif len(path_list) > 1:
        new_multiple_consensus.append(path_list)
        for i in path_list:
            if re.search(r'pass',os.path.basename(i)):
                keeped_path = i
                break
            elif re.search(r'flag',os.path.basename(i)):
                keeped_path = i
                break
    else:
        keeped_path = path_list[0]

    if re.search(r'rej',os.path.basename(keeped_path)):
        print("REJ FOR ",keeped_path)
        new_rej_consensus.append(keeped_path)
        continue

    keeped_path = re.sub(r'/home/fournie1/',r'/mnt/BelugaEric/',keeped_path)
    #print(path)
    new_keeped_consensus_list.append(keeped_path)    
        
print("len(new_keeped_consensus_list) : ", len(new_keeped_consensus_list))
print("nb_gisaid_nanapore_consensus  : ",nb_gisaid_nanapore_consensus)
print("gisaid_consensus_not_found_in_new : ", gisaid_consensus_not_found_in_new)
print("new_rej_consensus.append : ",new_rej_consensus.append)
