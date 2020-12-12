# coding=utf-8

"""
Eric Fournier 2020-11-21
"""


import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import re
import os

_debug = False


base_dir = "/data/PROJETS/Covid19_NextStrainBuilds/data_20201030_PASS_FLAG_minmaxSampleDate_2020-02-01_2020-06-01_WithContext_v2/LineageAndClade/"

outdir = os.path.join(base_dir,"LineageAndClade")


metadata = os.path.join("/data/Applications/GitScript/Covid19_V2/NextStrainFiles/data/gisaid/all/","metadata_2020-10-12_07-15.tsv")
seq_file = os.path.join(base_dir,"results","subsampled_alignment_quebec.fasta")
lineage_out = os.path.join(outdir,"lineage_res.tsv")


pd_lineage = pd.DataFrame(columns=['strain','pangolin_lineage','gisaid_clade'])
nb_rec_pd_lineage = 0

gisaid_id_list = []

for rec in SeqIO.parse(seq_file,'fasta'):
    rec_id = rec.id
    gisaid_id_list.append(rec_id)


pd_metadata = pd.read_csv(metadata,sep="\t",index_col=False)
pd_metadata = pd_metadata.loc[~pd_metadata['strain'].str.contains('Canada/Qc-'),:]


check = 0

for index,row in pd_metadata.loc[:,].iterrows():
    check += 1
    sys.stdout.write("Check >>  %d\r"%check)
    strain = row['strain']
    pangolin_lineage = row['pangolin_lineage']
    gisaid_clade = row['GISAID_clade']


    if strain in gisaid_id_list:
        pd_lineage.loc[nb_rec_pd_lineage] = {'strain': strain,'pangolin_lineage':pangolin_lineage,'gisaid_clade':gisaid_clade}
        nb_rec_pd_lineage += 1

pd_lineage.to_csv(lineage_out,index=False,sep="\t")




