# coding=utf-8

"""
Eric Fournier 2020-09-24

"""

import mysql.connector
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
from Covid19DB import MySQLcovid19,MySQLcovid19Selector


"""
TODO
- Ajuster Usage exemple selon options ajoute
- ajouter option max sample date

"""


"""
Usage exemple:
python PrepareMetadataFromDSPdb.py --dest LSPQ --debug (pour mode debug et metadata forat lspqo)

python PrepareMetadataFromDSPdb.py --dest GC  (pour mode production et metadata format Genome Center)

python PrepareMetadataFromDSPdb.py --dest LSPQ --debug --all (en mode debug, format lspq, extraire tous les samples de la base de donnees)
"""

global base_dir

pd.options.display.max_columns = 100
logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser(description="Create metadata file")
parser.add_argument('--dest', '-d', required=True,help="choose lspq (lspq) or genome center (gc)",choices=['LSPQ','GC'])
parser.add_argument('--all',help='extract all samples in database',action='store_true')
parser.add_argument('--debug',help="run in debug mode",action='store_true')
args = parser.parse_args()

_debug_ = args.debug
metadata_destination = args.dest
extract_all_samples = args.all

base_dir = "/home/foueri01@inspq.qc.ca/temp/20200924/" if  _debug_ else "/data/PROJETS/COVID-19_Beluga/Metadata"


def BuildSeqList(pd_seq_list):

    start = time.time()

    seq_list = np.array([])

    for index, row in pd_seq_list.loc[:,].iterrows():
        pass
        seq_list = np.append(seq_list,row[2]) 
    
    end = time.time()
    print("In BuildSeqList  for ",end - start, " seconds")
    
    return(list(seq_list))


def CheckMissingSpec(pd_found_spec,all_spec):
    
    all_spec = set(all_spec)
    found_spec = set(pd_found_spec['sample'])

    missing_spec = list(all_spec - found_spec)
    if len(missing_spec) > 0:
        missing_spec = pd.Series(missing_spec)
    else:
        missing_spec = pd.Series([],dtype='object')

    return(missing_spec)


def CreateMetadata(pd_seq_list,metadata_destination,extract_all_samples):
    
    MySQLcovid19.SetConnection()
    seq_list = BuildSeqList(pd_seq_list)
    pd_metadata = MySQLcovid19Selector.GetMetadataAsPdDataFrame(MySQLcovid19.GetConnection(),seq_list,metadata_destination,extract_all_samples)
    
    pd_missing_spec = CheckMissingSpec(pd_metadata,seq_list)
    return([pd_metadata,pd_missing_spec])


class PdWriter:

    global today
    today = datetime.datetime.now().strftime("%Y-%m-%d")

    def __init__(self,destination):
        self.destination = destination
        self.message_no_missing_samples = "************** No missing samples ***************"

    def WritePdToFile(self,pd_metadata):

        self.metadata_out = os.path.join(base_dir,"OUT",self.destination,"metadata_out_{0}.tsv".format(today))
        pd_metadata.to_csv(self.metadata_out,sep="\t",index=False)


    def WritePdMissingSamplesToFile(self,pd_missing_sample):
        self.missing_spec_out = os.path.join(base_dir,"OUT",self.destination,"missing_samples_{0}.tsv".format(today))

        if(pd_missing_sample.size == 0):
            pd_missing_sample[1] = self.message_no_missing_samples
            pd_missing_sample.to_csv(self.missing_spec_out,index=False,header=False)
        else:
            pd_missing_sample.to_csv(self.missing_spec_out,header=["Missing_samples"],index=False)


class PdBuilderFromFile:
    def __init__(self):
        pass
    def ReadSeqFileList(self):
        
        self.seq_file_list = os.path.join(base_dir,"IN","AllSample.20200912_small.list")
        return pd.read_csv(self.seq_file_list,sep="\t",index_col=False)


def Main():
    logging.info("In Main()")
    pd_builder_from_file = PdBuilderFromFile()
    pd_seq_list = pd_builder_from_file.ReadSeqFileList()

    pd_writer =  PdWriter(metadata_destination)
    pd_metadata,pd_missing_spec =  CreateMetadata(pd_seq_list,metadata_destination,extract_all_samples)

    pd_writer.WritePdToFile(pd_metadata)
    pd_writer.WritePdMissingSamplesToFile(pd_missing_spec)
    logging.info("Terminé")

if __name__ == '__main__':
    Main()
