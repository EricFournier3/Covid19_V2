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
import yaml
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
    input_file = '/data/PROJETS/ScriptDebug/GetGenomeCenterData/IN/2021-01-08_LSPQReport_test2.tsv'
    gisaid_metadata =  '/data/PROJETS/ScriptDebug/GetGenomeCenterData/IN/metadata_2020-12-20_12-24_test.tsv'
else:
    base_dir = "/data/PROJETS/COVID-19_Beluga/"


#print('INPUT ',input_file)


Tools.CheckUserArgs(gisaid_metadata_filename,mode)


class CovBankDB:
    def __init__(self):
        self.yaml_conn_param = open('/data/Databases/CovBanQ_Epi/CovBankParam.yaml')
        self.ReadConnParam()
        self.connection = self.SetConnection()

    def GetCursor(self):
        return self.GetConnection().cursor()

    def Commit(self):
        self.connection.commit()

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

    def GetGisaidMetadataAsPdDataFrame(self,sample_list):
        sample_list = '|'.join(sample_list)
        
        Prelevements_alias = "pr"
        Patients_alias = "p"

        join_column = 'ID_PATIENT'

        lspq_name = "Laboratoire de santé publique du Québec"
        lspq_addr = "20045, chemin Sainte-Marie, Sainte-Anne-de-Bellevue, QC, Canada"

        authors = "Sandrine Moreira, Ioannis Ragoussis, Guillaume Bourque, Jesse Shapiro, Mark Lathrop and Michel Roger on behalf of the CoVSeQ research group (http://covseq.ca/researchgroup)"

        SUBMITTER = "\'SandrineMoreiraLSPQ\' as submitter"
        FN = "\'all_sequences.fasta\' as fn"
        COVV_VIRUS_NAME = "{0}.GENOME_QUEBEC_REQUETE as covv_virus_name".format(Prelevements_alias) #decoration ajouter dans CovGisaidSubmissionV2.py
        COVV_TYPE = "\'betacoronavirus\' as covv_type"
        COVV_PASSAGE = "\'Original\' as covv_passage"
        COVV_COLLECTION_DATE = "{0}.DATE_PRELEV as covv_collection_date".format(Prelevements_alias)
        COVV_LOCATION = "\'North-America / Canada / Quebec\' as covv_location"
        COVV_ADD_LOCATION = "\' \' as covv_add_location"
        COVV_HOST = "\'Human\' as covv_host"
        COVV_ADD_HOST_INFO = "\' \' as covv_add_host_info"
        COVV_GENDER = "{0}.SEXE as covv_gender".format(Patients_alias)

        DTNAISSINFO = "{0}.DTNAISS".format(Patients_alias) # va permettre de determiner l age

        COVV_PATIENT_STATUS = "\'unknown\' as covv_patient_status"
        COVV_SPECIMEN = "\' \' as covv_specimen"
        COVV_OUTBREAK = "\' \' as covv_outbreak"
        COVV_LAST_VACCINATED = "\' \' as covv_last_vaccinated"
        COVV_TREATMENT = "\' \' as covv_treatment"

        COVV_ASSEMBLY_METHOD = "\' \' as covv_assembly_method"
        COVV_COVERAGE = "\' \' as covv_coverage"
        COVV_ORIG_LAB = "{0}.NOM_HOPITAL as covv_orig_lab".format(Prelevements_alias)
        COVV_ORIG_LAB_ADDR = "{0}.ADRESSE_HOPITAL as covv_orig_lab_addr".format(Prelevements_alias)
        COVV_PROVIDER_SAMPLE_ID = "\' \' as covv_provider_sample_id"
        COVV_SUBM_LAB = "'" + lspq_name  + "' as covv_subm_lab"
        COVV_SUBM_LAB_ADDR = "'" + lspq_addr + "' as covv_subm_lab_addr"
        COVV_SUBM_SAMPLE_ID = "{0}.GENOME_QUEBEC_REQUETE as covv_subm_sample_id".format(Prelevements_alias)
        COVV_AUTHORS = "'" + authors + "' as covv_authors"

        sql = "SELECT {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25} FROM Prelevements {26} inner join Patients {27} on {27}.{28} = {26}.{28} WHERE {26}.GENOME_QUEBEC_REQUETE REGEXP '{29}'".format(SUBMITTER,FN,COVV_VIRUS_NAME,COVV_TYPE,COVV_PASSAGE,COVV_COLLECTION_DATE,COVV_LOCATION,COVV_ADD_LOCATION,COVV_HOST,COVV_ADD_HOST_INFO,COVV_GENDER,DTNAISSINFO,COVV_PATIENT_STATUS,COVV_SPECIMEN,COVV_OUTBREAK,COVV_LAST_VACCINATED,COVV_TREATMENT,COVV_ASSEMBLY_METHOD,COVV_COVERAGE,COVV_ORIG_LAB,COVV_ORIG_LAB_ADDR,COVV_PROVIDER_SAMPLE_ID,COVV_SUBM_LAB,COVV_SUBM_LAB_ADDR,COVV_SUBM_SAMPLE_ID,COVV_AUTHORS,Prelevements_alias,Patients_alias,join_column,sample_list)

        #print("SQL ",sql)
        df = pd.read_sql(sql,con=self.GetConnection())
        return(df)
        


class NextstrainDataManager:
    def __init__(self):
        pass


class PangolinDataManager:
    def __init__(self):
        pass

class VariantDataManager: 
    def __init__(self):
        pass


class GisaidDataSubmissionManager:
    def __init__(self,input_file,gisaid_metadata,run_name,techno,cov_bank_db_obj):
        self.input_file = input_file
        self.gisaid_metadata = gisaid_metadata
       
        self.cov_bank_db_obj = cov_bank_db_obj
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
        consensus_out = os.path.join(self.submission_dir,"all_sequences.fasta")
        rec_list =  [] 
        for index,row in self.input_file_df.iterrows():
            #print("ROW ",row)
            sample_name = row['Sample Name']
            run_name = row['run_name']
            
            consensus = GenomeCenterConnector.GetConsensusPath(self.techno,sample_name,run_name)
            if consensus is None:
                logging.warning("No consensus for " + sample_name)
                return()

            rec = SeqIO.read(consensus,'fasta')
            try:
                print(rec.description)
                parsed_header = re.search(r'(Canada/Qc-)(\S+)/(\d{4}) seq_method:(\S+)\|assemb_method:\S+\|snv_call_method:\S+',rec.description)
                #print("rec id ",rec.id)
                method = parsed_header.group(4)
                self.sample_to_submit_dict[rec.id] = {}
                self.sample_to_submit_dict[rec.id]['method'] = method
                self.sample_to_submit_dict[rec.id]['gisaid_id'] = 'hCoV-19/' + rec.id
                rec.description = ""
                rec.id = self.sample_to_submit_dict[rec.id]['gisaid_id']
                rec_list.append(rec)
                 
            except:
                logging.error("Bug parse " + rec.description)
         
        #print(self.sample_to_submit_dict)
        SeqIO.write(rec_list,consensus_out,'fasta')

    def BuildSampleList(self):
        self.sample_list = []
        for sample in self.sample_to_submit_dict:
            print("SAMPLE ,", sample)
            short_sample_name = re.search(r'Canada/Qc-(\S+)/\d+',sample).group(1)
            self.sample_list.append(short_sample_name)

        print(self.sample_list)

    def CreateMetadata(self):
        metadata_out = os.path.join(self.submission_dir,"{0}_ncov19_metadata.xls".format(self.today))
        print("METADATA OUT ",metadata_out)
        metadata_df = self.cov_bank_db_obj.GetGisaidMetadataAsPdDataFrame(self.sample_list) 
        #TODO enrergistrer les sample non trouve dans db

        self.CheckMissingMetadata(metadata_df['covv_subm_sample_id'])

    def CheckMissingMetadata(self,sample_in_metadata):
        #print(set(sample_in_metadata))
        #print(set(self.sample_list))
        missing = set(self.sample_list) - set(sample_in_metadata) 
        test_set = set(['UN','DEUX'])
        df = pd.DataFrame({'MissingSample':list(test_set)})
        df.to_csv(os.path.join(self.submission_dir,"MissingMetadata.tsv"),sep="\t",index=False)

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
            #TODO save consensus non trouve
            return(None)
            

def Main():
    logging.info('Begin')

    cov_bank_db_obj = CovBankDB()

    GenomeCenterConnector.MountBelugaServer()
    gisaid_data_submission_manager = GisaidDataSubmissionManager(input_file,gisaid_metadata,beluga_run,seq_techno,cov_bank_db_obj)

    gisaid_data_submission_manager.GetConsensus()
    gisaid_data_submission_manager.BuildSampleList()
    gisaid_data_submission_manager.CreateMetadata()

    cov_bank_db_obj.CloseConnection()

if __name__ == '__main__':
    Main()
    logging.info("Terminé")

