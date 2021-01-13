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
args = parser.parse_args()

_debug_ = args.debug

input_file_name = args.input
gisaid_metadata_filename = args.gisaid_metadata


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
    gisaid_publication_basedir = '/data/PROJETS/ScriptDebug/GetGenomeCenterData/OUT/Gisaid/'
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
    def __init__(self):
        pass

class GenomeCenterConnector:
    beluga_user_dict = {'foueri01':['fournie1','BelugaEric'],'morsan01':['moreiras','BelugaSam']}
    beluga_server =  "{0}@beluga.computecanada.ca:/home/{0}".format(beluga_user_dict[user_name][0])
    mnt_beluga_server = "/mnt/{}/".format(beluga_user_dict[user_name][1])

    @staticmethod
    def MountBelugaServer():
        logging.info("Try to mount Beluga")
        os.system("sudo umount " + GenomeCenterConnector.mnt_beluga_server)
        os.system("sudo sshfs -o allow_other -o follow_symlinks {0} {1}".format(GenomeCenterConnector.beluga_server,GenomeCenterConnector.mnt_beluga_server))
        logging.info("Beluga mounted")

def Main():
    logging.info('Begin')
    #GenomeCenterConnector.MountBelugaServer()

if __name__ == '__main__':
    Main()
    logging.info("Terminé")

