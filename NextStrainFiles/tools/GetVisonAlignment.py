# coding=utf-8

'''
Eric Fournier 2020-11-12
Note: need to run this script in nextstrain conda env

'''

import logging
import datetime
import pandas as pd
import os
import numpy as np
import re
import sys
import logging
import time
import glob
from Bio import SeqIO
import subprocess
logging.basicConfig(level=logging.DEBUG)


base_dir = "/data/Users/Eric/Covid19/Vison/"

metadata_all = os.path.join(base_dir,"metadata_2020-11-12_10-27.tsv")
pd_metadata_all = pd.read_csv(metadata_all,sep="\t",index_col=False)
pd_metadata_mink = pd_metadata_all.loc[pd_metadata_all['host'] == 'Neovison vison',:]
#print(pd_metadata_mink.shape[0]) #226

subprocess.call(["augur","-h"],stdout=open("/data/Users/Eric/Covid19/Vison/augur_out.txt","w"))

exit(0)
for index,row in pd_metadata_mink.iterrows():
    strain = str(row['strain'])
    #print(strain)
    subprocess.call(["seqkit","grep","-r","-p","^" + str(strain),"/data/Users/Eric/Covid19/Vison/sequences_2020-11-12_07-33_small.fasta"],stdout=open("/data/Users/Eric/Covid19/Vison/myres.fasta","a"))


