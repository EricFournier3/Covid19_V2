# coding=utf-8
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import re
import os

base_dir = "/data/PROJETS/Covid19_NextStrainBuilds/data_20201030_PASS_FLAG_minmaxSampleDate_2020-02-01_2020-06-01_WithContext_v2/temp/"

metadata = os.path.join(base_dir,"metadata_2020-10-12_07-15_small.tsv")
seq_file = os.path.join(base_dir,"subsampled_alignment_quebec_small.fasta")
lineage_out = os.path.join(base_dir,"lineage_res.tsv")


pd_lineage = pd.DataFrame(columns=['strain','pangolin_lineage','gisaid_clade'])
nb_rec_pd_lineage = 0

gisaid_id_list = []

for rec in SeqIO.parse(seq_file,'fasta'):
    rec_id = rec.id
    #print("REC ID ",rec_id)
    gisaid_id_list.append(rec_id)

#print("LIST ",gisaid_id_list)

pd_metadata = pd.read_csv(metadata,sep="\t",index_col=False)
pd_metadata = pd_metadata.loc[~pd_metadata['strain'].str.contains('Canada/Qc-'),:]

#print(pd_metadata)

for index,row in pd_metadata.loc[:,].iterrows():
    #print(row['strain']) 
    strain = row['strain']
    pangolin_lineage = row['pangolin_lineage']
    gisaid_clade = row['GISAID_clade']


    if strain in gisaid_id_list:
        pd_lineage.loc[nb_rec_pd_lineage] = {'strain': strain,'pangolin_lineage':pangolin_lineage,'gisaid_clade':gisaid_clade}
        nb_rec_pd_lineage += 1


print(pd_lineage)
pd_lineage.to_csv(lineage_out,index=False,sep="\t")




