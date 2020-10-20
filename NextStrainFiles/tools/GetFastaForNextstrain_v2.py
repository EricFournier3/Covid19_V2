# coding=utf-8

"""
Eric Fournier 2020-10-19

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
import glob
from Bio import SeqIO
from Covid19DB import MySQLcovid19,MySQLcovid19Selector


"""
TODO

"""


"""
Usage exemple:

python GetFastaForNextstrain_v2.py -i /home/foueri01@inspq.qc.ca/temp/TEST_GET_FASTA_FOR_NEXTSTRAIN/IN/updateListOfRuns_v2_Dev_2020-10-17_consensusList_small2.list --debug --keepflag


"""

base_dir_dsp_db = "/data/Databases/CovBanQ_Epi/"

global base_dir


pd.options.display.max_columns = 100
logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser(description="Download Beluga consensus and Create metadata file")
parser.add_argument('--debug',help="run in debug mode",action='store_true')
parser.add_argument('--input','-i',help="Beluga fasta list file",required=True)
parser.add_argument('--onlymetadata',help="Only produce metadata",action='store_true')
parser.add_argument('--maxsampledate',help="Maximum sampled date YYYY-MM-DD",required=True)
parser.add_argument('--keep',help="Minimum Qc status to keep",choices=['ALL','PASS','FLAG'],required=True)
args = parser.parse_args()

_debug_ = args.debug
beluga_fasta_file  = args.input
only_metadata = args.onlymetadata
max_sample_date = args.maxsampledate

global qc_keep

qc_keep = args.keep


try:
    max_sample_date = datetime.datetime.strptime(max_sample_date,"%Y-%m-%d")
except:
    logging.error("format de date incorrect : YYYY-MM-DD")
    exit(1)

global fasta_qual_status_to_keep


if qc_keep == 'ALL':
    fasta_qual_status_to_keep = ['PASS','FLAG','REJ','NA','MISSING_METRICS_HEADER','MISSING_CONS_PERC_N']
elif qc_keep == 'FLAG':
    fasta_qual_status_to_keep = ['PASS','FLAG']
else:
    fasta_qual_status_to_keep = ['PASS']

global beluga_server
beluga_server = "fournie1@beluga.computecanada.ca:/home/fournie1"

global mnt_beluga_server
mnt_beluga_server = "/mnt/BelugaEric/"

global liste_envoi_genome_quebec
global sgil_extract
global fasta_outdir
global metadata_ourdir

if _debug_:
    sgil_extract = "/home/foueri01@inspq.qc.ca/temp/20200924/extract_with_Covid19_extraction_v2_20200923_CovidPos_small.txt"
    #liste_envoi_genome_quebec = os.path.join(base_dir_dsp_db,"/home/foueri01@inspq.qc.ca/temp/20200924/ListeEnvoisGenomeQuebec_small6.xlsx")
    liste_envoi_genome_quebec = os.path.join(base_dir_dsp_db,"LISTE_ENVOIS_GENOME_QUEBEC","EnvoiSmall.xlsx")
else:
    sgil_extract = os.path.join(base_dir_dsp_db,"SGIL_EXTRACT","extract_with_Covid19_extraction_v2_20200923_CovidPos.txt")
    liste_envoi_genome_quebec = os.path.join(base_dir_dsp_db,"LISTE_ENVOIS_GENOME_QUEBEC","ListeEnvoisGenomeQuebec_2020-08-28CORR.xlsx")

base_dir = "/home/foueri01@inspq.qc.ca/temp/TEST_GET_FASTA_FOR_NEXTSTRAIN/" if  _debug_ else "/data/PROJETS/COVID-19_Beluga/Metadata"
in_dir = os.path.join(base_dir,"IN")
out_dir = os.path.join(base_dir,"OUT")

fasta_outdir = os.path.join(out_dir,'FASTA')
metadata_outdir = os.path.join(out_dir,'METADATA')

if not os.path.isfile(os.path.join(in_dir,beluga_fasta_file)):
    logging.error("File missing : " + os.path.join(in_dir,beluga_fasta_file))
    exit(1)


##################################################################################
def MountBelugaServer():
    os.system("sudo umount " + mnt_beluga_server)
    os.system("sudo sshfs -o allow_other -o follow_symlinks {0} {1}".format(beluga_server,mnt_beluga_server))

class MetadataManager():
    def __init__(self,samples_list,pd_fasta_list):
        self.samples_list = samples_list
        self.pd_fasta_list = pd_fasta_list

        self.pd_envoi_qenome_quebec = None
        self.pd_sgil_extract = None
        self.pd_metadata = None

        self.SetPdEnvoisGenomeQuebec()
        self.SetPdSgilExtract()

    def SetPdEnvoisGenomeQuebec(self):
        self.pd_envoi_genome_quebec = pd.read_excel(liste_envoi_genome_quebec,sheet_name=0)

    def SetPdSgilExtract(self):
        self.pd_sgil_extract = pd.read_csv(sgil_extract,sep="\t",index_col=False)

    def CheckMissingSpec(self,pd_found_spec,all_spec):
        all_spec = set(all_spec)
        found_spec = set(pd_found_spec['sample'])
        missing_spec = list(all_spec - found_spec)
    
        if len(missing_spec) > 0:
            missing_spec = pd.Series(missing_spec)
        else:
            missing_spec = pd.Series([],dtype='object')

        return(missing_spec)

    def AddMissingFromSgilExtract(self,pd_missing_samples,pd_sgil_extract,metadata_columns):
        res_df = self.pd_sgil_extract.loc[self.pd_sgil_extract['NUMERO_SGIL'].isin(list(pd_missing_samples)),['NUMERO_SGIL','SAMPLED_DATE','TRAVEL_HISTORY','CT','SEX','RSS_PATIENT','DATE_NAISS']].copy()
        res_df = res_df.rename(columns={'NUMERO_SGIL':'sample','SAMPLED_DATE':'sample_date','TRAVEL_HISTORY':'country_exposure','CT':'ct','SEX':'sex','RSS_PATIENT':'rss','DATE_NAISS':'date_naiss'})
        res_df['country'] = 'Canada'
        res_df['division'] = 'Quebec'

        return(res_df)

    def AddMissingFromEnvoisGenomeQuebec(self,pd_missing_samples,pd_envoi_genome_quebec,metadata_columns):
        res_df = self.pd_envoi_genome_quebec.loc[self.pd_envoi_genome_quebec['# Requête'].isin(list(pd_missing_samples)),['# Requête','Date de prélèvement','Date de naissance']].copy()
        res_df = res_df.rename(columns={'# Requête':'sample','Date de prélèvement':'sample_date','Date de naissance':'date_naiss'})
        res_df['country'] = 'Canada'
        res_df['country_exposure'] = 'INDETERMINE'
        res_df['ct'] = '0'
        res_df['rss'] = 'INDETERMINE'
        res_df['sex'] = 'INDETERMINE'
        res_df['division'] = 'Quebec'

        return(res_df)

    def CreateMetadata(self,max_sample_date):
        MySQLcovid19.SetConnection()
        self.pd_metadata = MySQLcovid19Selector.GetMetadataAsPdDataFrame(MySQLcovid19.GetConnection(),self.samples_list,'LSPQ',False)
        #print(self.pd_metadata)
        self.pd_metadata['sample'] = self.pd_metadata['sample'].str.replace('LSPQ-','')
        self.pd_metadata['sample'] = self.pd_metadata['sample'].str.strip(' ')
        pd_missing_samples = self.CheckMissingSpec(self.pd_metadata,self.samples_list)

        pd_missing_get_from_sgil_extract = self.AddMissingFromSgilExtract(pd_missing_samples,self.pd_sgil_extract,self.pd_metadata.columns)
        self.pd_metadata = pd.concat([self.pd_metadata,pd_missing_get_from_sgil_extract])

        self.pd_metadata['sample'] = self.pd_metadata['sample'].str.replace('LSPQ-','')
        self.pd_metadata['sample'] = self.pd_metadata['sample'].str.strip(' ')
        pd_missing_samples = self.CheckMissingSpec(self.pd_metadata,self.samples_list)

        pd_missing_get_from_EnvoisGenomeQuebec = self.AddMissingFromEnvoisGenomeQuebec(pd_missing_samples,self.pd_envoi_qenome_quebec,self.pd_metadata.columns)
        self.pd_metadata = pd.concat([self.pd_metadata,pd_missing_get_from_EnvoisGenomeQuebec])

        self.pd_metadata['sample'] = self.pd_metadata['sample'].str.replace('LSPQ-','')
        self.pd_metadata['sample'] = self.pd_metadata['sample'].str.strip(' ')
        self.pd_missing_samples = self.CheckMissingSpec(self.pd_metadata,self.samples_list)

        self.pd_metadata['sample_date'] = pd.to_datetime(self.pd_metadata.sample_date)
        self.pd_metadata['sample_date'] = self.pd_metadata['sample_date'].dt.strftime('%Y-%m-%d')

        self.pd_metadata = self.pd_metadata.drop_duplicates(subset='sample',keep='first')

        self.pd_samples_missing_rss = self.pd_metadata.loc[self.pd_metadata['rss'] == 'INDETERMINE',['sample']]
        self.pd_metadata = self.pd_metadata.loc[self.pd_metadata['rss'] != 'INDETERMINE',:]
        self.pd_metadata.loc[self.pd_metadata['sample'].str.contains('HGA-'),'sample']= self.pd_metadata['sample'] + '2D'
        self.pd_metadata = self.pd_metadata.sort_values(by=['sample'])

        self.pd_metadata['sample_date']= self.pd_metadata['sample_date'].astype('datetime64[ns]')
        #print(self.pd_metadata.dtypes)

        self.pd_metadata = self.pd_metadata.loc[self.pd_metadata['sample_date'] <= max_sample_date,:]
        #print(self.pd_metadata)
        #print(self.pd_samples_missing_rss)

    def GetBelugaRunsWithTargetSamples(self):
        self.pd_samples_with_run_name = pd.merge(self.pd_metadata,self.pd_fasta_list,left_on='sample',right_on='SAMPLE',how='left',indicator=True)
        self.pd_samples_with_run_name = self.pd_samples_with_run_name[['sample','sample_date','STATUS','TECHNO','RUN_NAME']]
        #print(self.pd_with_run_name)
        #self.pd_run_list = pd.DataFrame({'RUN_NAME':self.pd_samples_with_run_name['RUN_NAME'].unique()})
        self.pd_run_list = self.pd_samples_with_run_name.drop_duplicates(['RUN_NAME','TECHNO'])
        self.pd_run_list = self.pd_run_list[['RUN_NAME','TECHNO']]
        #print(self.pd_run_list)

    def GetPdMetadata(self):
        return self.pd_metadata

    def GetPdMissingSamples(self):
        return self.pd_missing_samples

    def WriteMetadata(self,beluga_fasta_file,max_sample_date):

        if qc_keep  == 'FLAG':
            qc_status_suffix = "PASS_FLAG"
        else:
            qc_status_suffix = qc_keep


        today = datetime.datetime.now().strftime("%Y-%m-%d")

        self.metadata_out = os.path.join(metadata_outdir,"metadata_{0}_from_{1}.tsv".format(today,os.path.basename(beluga_fasta_file)))
        self.missing_samples_out = os.path.join(metadata_outdir,"missing_samples_{0}_from_{1}.tsv".format(today,os.path.basename(beluga_fasta_file))) 
        self.samples_missing_rss_out = os.path.join(metadata_outdir,"samples_missing_rss_{0}_from_{1}.tsv".format(today,os.path.basename(beluga_fasta_file)))
        self.sample_with_run_name_out = os.path.join(metadata_outdir,"Samples_{2}_maxSampleDate_{0}_from_{1}".format(str(max_sample_date),os.path.basename(beluga_fasta_file),qc_status_suffix))
        self.runs_with_samples_max_sample_date_out = os.path.join(metadata_outdir,"Runs_with_samples_{2}_maxSampleDate_{0}_from_{1}".format(max_sample_date,os.path.basename(beluga_fasta_file),qc_status_suffix))

        self.pd_metadata.to_csv(self.metadata_out,sep="\t",index=False)
        self.pd_missing_samples.to_csv(self.missing_samples_out,sep="\t",index=False)
        self.pd_samples_missing_rss.to_csv(self.samples_missing_rss_out,sep="\t",index=False)
        self.pd_samples_with_run_name.to_csv(self.sample_with_run_name_out,sep="\t",index=False)
        self.pd_run_list.to_csv(self.runs_with_samples_max_sample_date_out,sep="\t",index=False)


class FastaGetter():
    def __init__(self,pd_metadata,pd_fasta_list):
        self.pd_metadata = pd_metadata
        self.pd_fasta_list = pd_fasta_list

        self.fasta_rec_list = []

        self.samples_to_keep = list(self.pd_metadata['sample'])
        self.pd_selected_fasta = None

    def SelectFasta(self):
        self.pd_selected_fasta = self.pd_fasta_list.loc[self.pd_fasta_list['SAMPLE'].isin(self.samples_to_keep),:]
        
        #on groupe par sample/statut et on garde celui avec le moins de N
        idx  = self.pd_selected_fasta.groupby(['SAMPLE','STATUS'])['PERC_N'].transform(min) == self.pd_selected_fasta['PERC_N']
        self.pd_selected_fasta = self.pd_selected_fasta[idx]

        #on groupe ensuite par sample et on conserve les PASS si pour un sample on a un PASS et un FLAG 
        idx  = self.pd_selected_fasta.groupby(['SAMPLE'])['STATUS'].transform(max) == self.pd_selected_fasta['STATUS']
        self.pd_selected_fasta = self.pd_selected_fasta[idx]

        #si on a une egalite de PERC_N, on garde illumina
        idx  = self.pd_selected_fasta.groupby(['SAMPLE'])['TECHNO'].transform(min) == self.pd_selected_fasta['TECHNO']
        self.pd_selected_fasta = self.pd_selected_fasta[idx]


    def GetFastaFromBeluga(self,beluga_fasta_file):
        for index, row in self.pd_selected_fasta.iterrows():
            fasta_path = str(row['PATH'])
            fasta_path = re.sub(r'/genfs/projects/',mnt_beluga_server,fasta_path)
            logging.info("Get " + fasta_path)
            rec = SeqIO.read(fasta_path,'fasta')
            self.fasta_rec_list.append(rec)

        fasta_out = os.path.join(fasta_outdir,"temp.fasta")
        SeqIO.write(self.fasta_rec_list,fasta_out,'fasta')
        self.AddPathToHeader(fasta_out,os.path.basename(beluga_fasta_file),fasta_path)
        os.remove(os.path.join(fasta_outdir,"temp.fasta"))

    def AddPathToHeader(self,fasta_in,beluga_fasta_file,fasta_path):
        rec_list = []
        for rec in SeqIO.parse(fasta_in,'fasta'):
            rec.description = rec.description + " " +  fasta_path
            rec_list.append(rec)

        SeqIO.write(rec_list,os.path.join(fasta_outdir,"sequences_from_" + beluga_fasta_file),'fasta')

        

class FastaListManager():
    def __init__(self,fasta_list_file):
        self.fasta_list_file = fasta_list_file
        self.pd_fasta_list = self.GetPdFastaList()
        self.BuildSamplesList()

    def GetPdFastaList(self):
        pd_df = pd.read_csv(self.fasta_list_file,sep="\t",index_col=False)
        return(pd_df.loc[pd_df['STATUS'].isin(fasta_qual_status_to_keep),:])

    def GetSamplesList(self):
        return self.samples_list

    def BuildSamplesList(self):
        start = time.time()
        self.samples_list = np.array([])
        self.samples_list = self.pd_fasta_list.loc[self.pd_fasta_list['STATUS'].isin(fasta_qual_status_to_keep),'SAMPLE'].unique()

        i = 0
        for sample in self.samples_list:
            if sample.startswith('HGA-'):
                self.samples_list[i] = re.sub(r'2D$','',sample)
            i+=1

        end = time.time()

        print("In BuildSamplesList  for ",end - start, " seconds")


def Main():
    logging.info("In Main()")

    #MountBelugaServer()
    fasta_list_manager = FastaListManager(os.path.join(in_dir,beluga_fasta_file))
    samples_list = fasta_list_manager.GetSamplesList()
    metadata_manager = MetadataManager(samples_list,fasta_list_manager.GetPdFastaList())
    metadata_manager.CreateMetadata(max_sample_date)
    metadata_manager.GetBelugaRunsWithTargetSamples()
    metadata_manager.WriteMetadata(beluga_fasta_file,max_sample_date.strftime("%Y-%m-%d"))

    if not only_metadata and qc_keep in ['PASS','FLAG']:
        fasta_getter = FastaGetter(metadata_manager.GetPdMetadata(),fasta_list_manager.GetPdFastaList())
        fasta_getter.SelectFasta()
        fasta_getter.GetFastaFromBeluga(beluga_fasta_file)


if __name__ == '__main__':
    Main()
    logging.info("Terminé")
