# coding=utf-8

"""
Eric Fournier 2021-01-12


TODO
- check args
"""

import shutil
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
import glob
from Bio import SeqIO

pd.options.display.max_columns = 100
logging.basicConfig(level=logging.DEBUG)

class Tools:
    @staticmethod
    def CheckDateFormat(mydate):
        try:
            return(datetime.datetime.strptime(mydate,"%Y-%m-%d"))
        except ValueError:
            msg = 'Format de date incorrect pour {0} : should be YYYY-MM-DD'.format(mydate)
            raise argparse.ArgumentTypeError(msg)
    @staticmethod
    def CheckUserArgs(gisaid_metadata,mode):
            if mode == 'submission':
                if not gisaid_metadata:
                    logging.error('Gisaid metadata argument missing')
                    exit(1)

parser = argparse.ArgumentParser(description="Download Genome Center data for Nextstrain, Pangolin, variant analysis and GISAID/NML submission")
parser.add_argument('--debug',help="run in debug mode",action='store_true')
parser.add_argument('--mode',help='execution mode',choices=['nextstrain','pangolin','variant','submission'],required=True)
parser.add_argument('--input','-i',help="genome center report input file name",required=True)
parser.add_argument('--gisaid-metadata','-g',help="gisaid metadata file name")
parser.add_argument('--user','-u',help='user name',required=True,choices=['foueri01','morsan01'])
parser.add_argument('--minsampledate',help="Minimum sampled date YYYY-MM-DD",type=Tools.CheckDateFormat)
parser.add_argument('--maxsampledate',help="Maximum sampled date YYYY-MM-DD",type=Tools.CheckDateFormat)
parser.add_argument('--keep',help="Minimum Qc status to keep",choices=['ALL','PASS','FLAG'])
parser.add_argument('--beluga-run',help="beluga run name")
parser.add_argument('--techno',help="sequencing technologie",choices=['illumina','nanopore','mgi'])

args = parser.parse_args()

_debug_ = args.debug

input_file_name = args.input
gisaid_metadata_filename = args.gisaid_metadata

global seq_techno
seq_techno = args.techno

global beluga_run
beluga_run = args.beluga_run

global user_name
user_name = args.user

global gisaid_publication_basedir

global gisaid_metadata_basedir


global mode
mode = args.mode

global min_sample_date
min_sample_date = args.minsampledate

global max_sample_date
max_sample_date = args.maxsampledate


global base_dir
global input_file
global gisaid_metadata


if _debug_:
    input_file = '/data/PROJETS/ScriptDebug/GetGenomeCenterData/IN/2021-01-08_LSPQReport_test.tsv'
    gisaid_metadata =  '/data/PROJETS/ScriptDebug/GetGenomeCenterData/IN/metadata_2020-12-20_12-24_test.tsv'
else:
    base_dir = "/data/PROJETS/COVID-19_Beluga/"


#print('INPUT ',input_file)


Tools.CheckUserArgs(gisaid_metadata_filename,mode)


class CovBankDB:
    def __init__(self):
        self.self.yaml_conn_param = open('/data/Databases/CovBanQ_Epi/CovBankParam.yaml')
        self.ReadConnParam()
        self.connection = self.SetConnection()

    def SetConnection(self):
        return mysql.connector.connect(host=self.host,user=self.user,password=self.password,database=self.database)

    def CloseConnection(self):
        self.GetConnection().close()

    def GetConnection(self):
        return self.connection

    def ReadConnParam(self):
        param = yaml.load(self.yaml_conn_param,Loader=yaml.FullLoader)
        self.host = param['host']
        self.user = param['user']
        self.password = param['password']
        self.database = param['database']


class NextstrainDataManager:
    def __init__(self):
        pass


class PangolinDataManager:
    def __init__(self):
        pass

class VariantDataManager: 
    def __init__(self):
        pass


class DataSubmissionManager:
    def __init__(self,input_file,gisaid_metadata,run_name,techno):
        self.input_file = input_file
        self.gisaid_metadata = gisaid_metadata
       
        self.techno = techno
        self.run_name = run_name
        self.today = datetime.datetime.now().strftime("%Y%m%d")

        self.sample_to_submit_dict = {}

        self.SetDataSubmissionDf()

        self.CreateSubmissionDirectory()

    def CreateSubmissionDirectory(self):
        if _debug_:
            self.gisaid_publication_basedir = '/data/PROJETS/ScriptDebug/GetGenomeCenterData/OUT/Gisaid/'
        else:
            self.gisaid_publication_basedir = '/data/PROJETS/ScriptDebug/GetGenomeCenterData/OUT/Gisaid/'

        self.submission_dir = os.path.join(self.gisaid_publication_basedir,self.techno,self.run_name,self.today)

        
        try:
            os.makedirs(self.submission_dir,exist_ok=True)
        except:
            logging.error("Impossible de creer " + self.submission_dir)
        
    def SetDataSubmissionDf(self):
        self.input_file_df = pd.read_csv(self.input_file,sep="\t",index_col=False)
        self.input_file_df = self.input_file_df.loc[self.input_file_df['run_name'] == beluga_run,: ]
        #print(self.input_file_df)

        self.input_file_df['ncov_tools.pass'] = self.input_file_df['ncov_tools.pass'].astype(str)

        self.gisaid_metadata_df = pd.read_csv(self.gisaid_metadata,sep="\t",index_col=False)
        self.gisaid_metadata_df = self.gisaid_metadata_df.loc[self.gisaid_metadata_df['strain'].str.contains('^Canada/Qc-\S+/\S+',regex=True),:]
        self.gisaid_metadata_qc_list = list(self.gisaid_metadata_df['strain'].str.replace(r'^Canada/Qc-(\S+)/\S+',r'\1',regex=True))

        self.input_file_df = self.input_file_df.loc[(~self.input_file_df['Sample Name'].isin(self.gisaid_metadata_qc_list)) & (self.input_file_df['PASS/FLAG/REJ'].str.upper() == 'PASS') & (self.input_file_df['ncov_tools.pass'].str.upper() ==  'TRUE'),['Sample Name','PASS/FLAG/REJ','run_name','ncov_tools.pass','platform']]

        #print(self.input_file_df)

    def GetConsensus(self):
        
        for index,row in self.input_file_df.iterrows():
            print("ROW ",row)
            sample_name = row['Sample Name']
            run_name = row['run_name']
            
            #TODO if not none
            consensus = GenomeCenterConnector.GetConsensusPath('illumina','L00232955','20200609_illumina_LSPQPlate05_HM2CTDRXX')
            rec = SeqIO.read(consensus,'fasta')
            parsed_header = re.search(r'(Canada/Qc-)(\S+)/(\d{4}) seq_method:(\S+)\|assemb_method:\S+\|snv_call_method:\S+ ',rec.description)
            #print("rec id ",rec.id)
            method = parsed_header.group(4)
            self.sample_to_submit_dict[rec.id] = {}
            self.sample_to_submit_dict[rec.id]['method'] = method
            self.sample_to_submit_dict[rec.id]['gisaid_id'] = 'hCoV-19/' + rec.id



class GenomeCenterConnector:
    beluga_user_dict = {'foueri01':['fournie1','BelugaEric'],'morsan01':['moreiras','BelugaSam']}
    beluga_server =  "{0}@beluga.computecanada.ca:/home/{0}".format(beluga_user_dict[user_name][0])
    mnt_beluga_server = "/mnt/{}/".format(beluga_user_dict[user_name][1])

    full_processing_path = os.path.join(mnt_beluga_server,'COVID_full_processing')
    

    @staticmethod
    def MountBelugaServer():
        logging.info("Try to mount Beluga")
        os.system("sudo umount " + GenomeCenterConnector.mnt_beluga_server)
        os.system("sudo sshfs -o allow_other -o follow_symlinks {0} {1}".format(GenomeCenterConnector.beluga_server,GenomeCenterConnector.mnt_beluga_server))
        logging.info("Beluga mounted")

    @staticmethod
    def GetConsensusPath(techno,sample,run_name):
        sample = sample.split('_')[0]
        if techno == 'nanopore':
            consensus_path = os.path.join(GenomeCenterConnector.full_processing_path,techno,run_name,'analysis','*_nanopolish_800x',sample + "*",sample + ".consensus." + techno + ".*.fasta")
        elif techno in ['mgi','illumina']:
            consensus_path = os.path.join(GenomeCenterConnector.full_processing_path,techno + "_reprocess", run_name,'consensus',sample,sample + ".consensus." + str(techno.upper() if techno == 'mgi' else techno) + ".*.fasta")
        try:
            return(glob.glob(consensus_path)[0])
        except:
            return(None)
            

def Main():
    logging.info('Begin')
    GenomeCenterConnector.MountBelugaServer()
    data_submission_manager = DataSubmissionManager(input_file,gisaid_metadata,beluga_run,seq_techno)
    #print(GenomeCenterConnector.GetConsensusPath(seq_techno,'L00214634',beluga_run))

    data_submission_manager.GetConsensus()

if __name__ == '__main__':
    Main()
    logging.info("Terminé")

