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
import glob
from Bio import SeqIO
from Covid19DB import MySQLcovid19,MySQLcovid19Selector


"""
TODO
- Ajuster Usage exemple selon options ajoute
- ajouter option max sample date

"""


"""
Usage exemple:
python PrepareMetadataFromDSPdb.py --dest LSPQ --debug  -i AllSample.20200912.list (pour mode debug et metadata forat lspqo) 

python PrepareMetadataFromDSPdb.py --dest GC -i AllSample.20200912.list  (pour mode production et metadata format Genome Center)

python PrepareMetadataFromDSPdb.py --dest LSPQ --debug --all   (en mode debug, format lspq, extraire tous les samples de la base de donnees)


Les fichiers AllSample.YYYYMMDD.list sont localise sur Compute Canada/Beluga dans /home/fournie1/COVID_LSPQ/TRACE et sont importe dans /data/PROJETS/COVID-19_Beluga/Metadata/IN

"""

base_dir_dsp_db = "/data/Databases/COVID19_DSP/"

global base_dir


pd.options.display.max_columns = 100
logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser(description="Create metadata file")
parser.add_argument('--dest', '-d', required=True,help="choose lspq (lspq) or genome center (gc)",choices=['LSPQ','GC'])
parser.add_argument('--all',help='extract all samples in database',action='store_true')
parser.add_argument('--debug',help="run in debug mode",action='store_true')
parser.add_argument('--input','-i',help="input file")
parser.add_argument('--keepflag',help="keep flag assemblies in fasta file",action='store_true')
args = parser.parse_args()

_debug_ = args.debug
metadata_destination = args.dest
extract_all_samples = args.all
genome_center_file  = args.input
keep_flag_assemblies = args.keepflag


if not genome_center_file and not extract_all_samples:
   logging.error("Genome Center file missing as argument")
   exit(0)

global fasta_qual_status_to_keep

if keep_flag_assemblies:
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

if _debug_:
    sgil_extract = "/home/foueri01@inspq.qc.ca/temp/20200924/extract_with_Covid19_extraction_v2_20200923_CovidPos.txt"
    #liste_envoi_genome_quebec = os.path.join(base_dir_dsp_db,"/home/foueri01@inspq.qc.ca/temp/20200924/ListeEnvoisGenomeQuebec_small6.xlsx")
    liste_envoi_genome_quebec = os.path.join(base_dir_dsp_db,"LISTE_ENVOIS_GENOME_QUEBEC","ListeEnvoisGenomeQuebec_2020-08-28CORR.xlsx")
else:
    sgil_extract = os.path.join(base_dir_dsp_db,"SGIL_EXTRACT","extract_with_Covid19_extraction_v2_20200923_CovidPos.txt")
    liste_envoi_genome_quebec = os.path.join(base_dir_dsp_db,"LISTE_ENVOIS_GENOME_QUEBEC","ListeEnvoisGenomeQuebec_2020-08-28CORR.xlsx")

base_dir = "/home/foueri01@inspq.qc.ca/temp/TEST_GET_FASTA_FOR_NEXTSTRAIN/" if  _debug_ else "/data/PROJETS/COVID-19_Beluga/Metadata"

#TODO A MODIFIER
fasta_outdir = os.path.join(base_dir,'FASTA') if _debug_ else os.path.join(base_dir,'FASTA') 


if genome_center_file and (not os.path.isfile(os.path.join(base_dir,"IN",genome_center_file))):
    logging.error("File missing : " + os.path.join(base_dir,"IN",genome_center_file))
    exit(1)

def BuildSeqList(pd_seq_list):

    start = time.time()
    seq_list = np.array([])

    for index, row in pd_seq_list.loc[:,].iterrows():
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

def AddMissingFromEnvoisGenomeQuebec(pd_missing_spec,pd_envoi_genome_quebec,metadata_columns):
    res_df = pd_envoi_genome_quebec.loc[pd_envoi_genome_quebec['# Requête'].isin(list(pd_missing_spec)),['# Requête','Date de prélèvement']].copy()
    res_df = res_df.rename(columns={'# Requête':'sample','Date de prélèvement':'sample_date'})
    res_df['country'] = 'Canada'
    res_df['country_exposure'] = 'INDETERMINE'
    res_df['ct'] = '0'
    res_df['rss'] = 'INDETERMINE'
    res_df['sex'] = 'INDETERMINE'
    res_df['division'] = 'Quebec'
 
    return(res_df)


def AddMissingFromSgilExtract(pd_missing_spec,pd_sgil_extract,metadata_columns):

    res_df = pd_sgil_extract.loc[pd_sgil_extract['NUMERO_SGIL'].isin(list(pd_missing_spec)),['NUMERO_SGIL','SAMPLED_DATE','TRAVEL_HISTORY','CT','SEX','RSS_PATIENT']].copy()
    res_df = res_df.rename(columns={'NUMERO_SGIL':'sample','SAMPLED_DATE':'sample_date','TRAVEL_HISTORY':'country_exposure','CT':'ct','SEX':'sex','RSS_PATIENT':'rss'})
    res_df['country'] = 'Canada'
    res_df['division'] = 'Quebec'

    return(res_df)


def CreateMetadataForAllDspDbSamples(metadata_destination):
    MySQLcovid19.SetConnection()
    pd_metadata = MySQLcovid19Selector.GetMetadataAsPdDataFrame(MySQLcovid19.GetConnection(),[""],metadata_destination,True)

    return(pd_metadata)

def CreateMetadata(pd_seq_list,metadata_destination,extract_all_samples,pd_envoi_qenome_quebec,pd_sgil_extract):
    
    MySQLcovid19.SetConnection()
    seq_list = BuildSeqList(pd_seq_list)
    pd_metadata = MySQLcovid19Selector.GetMetadataAsPdDataFrame(MySQLcovid19.GetConnection(),seq_list,metadata_destination,extract_all_samples)
    
    pd_missing_spec = CheckMissingSpec(pd_metadata,seq_list)

    if (extract_all_samples,pd_envoi_qenome_quebec is not None) and (pd_sgil_extract is not None):
        pd_missing_get_from_sgil_extract = AddMissingFromSgilExtract(pd_missing_spec,pd_sgil_extract,pd_metadata.columns)
        pd_metadata = pd.concat([pd_metadata,pd_missing_get_from_sgil_extract])

        pd_missing_spec = CheckMissingSpec(pd_metadata,seq_list)
        pd_missing_get_from_EnvoisGenomeQuebec = AddMissingFromEnvoisGenomeQuebec(pd_missing_spec,pd_envoi_qenome_quebec,pd_metadata.columns)

        pd_metadata = pd.concat([pd_metadata,pd_missing_get_from_EnvoisGenomeQuebec])
        pd_missing_spec = CheckMissingSpec(pd_metadata,seq_list)

    pd_metadata['sample_date'] = pd.to_datetime(pd_metadata.sample_date)
    pd_metadata['sample_date'] = pd_metadata['sample_date'].dt.strftime('%Y-%m-%d')
    return([pd_metadata,pd_missing_spec])


class PdWriter:

    global today
    today = datetime.datetime.now().strftime("%Y-%m-%d")

    def __init__(self,destination,in_file):
        self.in_file = in_file
        self.destination = destination
        self.message_no_missing_samples = "************** No missing samples ***************"

    def WritePdToFile(self,pd_metadata):

        self.metadata_out = os.path.join(base_dir,"OUT",self.destination,"metadata_out_{0}_from_{1}.tsv".format(today,self.in_file))
        pd_metadata = pd_metadata.sort_values(by=['sample'])
        pd_metadata.to_csv(self.metadata_out,sep="\t",index=False)

    def WritePdMissingSamplesToFile(self,pd_missing_sample):
        self.missing_spec_out = os.path.join(base_dir,"OUT",self.destination,"missing_samples_{0}_from_{1}.tsv".format(today,self.in_file))

        if(pd_missing_sample.size == 0):
            pd_missing_sample[1] = self.message_no_missing_samples
            pd_missing_sample.to_csv(self.missing_spec_out,index=False,header=False)
        else:
            pd_missing_sample = pd_missing_sample.sort_values(ascending=True)
            pd_missing_sample.to_csv(self.missing_spec_out,header=["Missing_samples"],index=False)

class FastaGetter():
    def __init__(self,pd_seq_list,pd_metadata):
        self.pd_seq_list = pd_seq_list
        self.pd_metadata = pd_metadata
        self.sample_with_complete_metadata = self.BuildSampleWithCompleteMetadata()
        self.pd_final_target_fasta_sample = self.GetFinalTargetFastaSample()
        self.fasta_rec = []


    def GetFinalTargetFastaSample(self):
        pd_temp = self.pd_seq_list.loc[(self.pd_seq_list['Sample[2]'].isin(self.sample_with_complete_metadata)) & (self.pd_seq_list['QCStatus{PASS,FLAG,REJ}[5]'].isin(fasta_qual_status_to_keep)),['Sample[2]','QCStatus{PASS,FLAG,REJ}[5]','REPOPATH[8]','RunDate[4]']]

        idx = pd_temp.groupby('Sample[2]')['RunDate[4]'].transform(max) == pd_temp['RunDate[4]']
        pd_fasta_target = pd_temp[idx]
        return(pd_fasta_target)

    def BuildSampleWithCompleteMetadata(self):
        #todo mettre les indetermine dans un fichier
        return list(self.pd_metadata.loc[self.pd_metadata['rss'] != 'INDETERMINE','sample'])

    def GetFastaFromBeluga(self):
        pass
        for index, row in self.pd_final_target_fasta_sample.iterrows():
            #print(index, " ____________ ", row[['Sample[2]','QCStatus{PASS,FLAG,REJ}[5]','REPOPATH[8]']])
             print("SAMPLE : ", str(row['Sample[2]']), " ---- PATH ",str(row['REPOPATH[8]']))
             fasta_path = re.sub(r'/genfs/projects/',mnt_beluga_server,str(row['REPOPATH[8]']))
             #print(fasta_path)
             fasta_list = glob.glob(fasta_path + "/*.fasta")
             #print(fasta_list)
             if len(fasta_list) > 0:
                 self.fasta_rec.append(SeqIO.read(fasta_list[0],'fasta'))

        #print(self.fasta_rec)
        SeqIO.write(self.fasta_rec,os.path.join(fasta_outdir,"sequences.fasta"),'fasta')


class PdBuilderFromFile:
    def __init__(self):
        pass

    def ReadSeqFileList(self,in_file):        
        self.seq_file_list = os.path.join(base_dir,"IN",in_file)
        return pd.read_csv(self.seq_file_list,sep="\t",index_col=False)

    def ReadEnvoiGenomeQuebecFile(self,in_file):
    
        return pd.read_excel(in_file,sheet_name=0)

    def ReadSgilExtract(self,in_file):
        
        return pd.read_csv(in_file,sep="\t",index_col=False)

def MountBelugaServer():
    os.system("sudo umount " + mnt_beluga_server)
    os.system("sudo sshfs -o allow_other -o follow_symlinks {0} {1}".format(beluga_server,mnt_beluga_server))

def Main():
    logging.info("In Main()")

    MountBelugaServer()

    if extract_all_samples:
        pd_writer =  PdWriter(metadata_destination,None)
        pd_metadata = CreateMetadataForAllDspDbSamples(metadata_destination)
        pd_writer.WritePdToFile(pd_metadata)
        return

    pd_writer =  PdWriter(metadata_destination,os.path.basename(genome_center_file))
    pd_builder_from_file = PdBuilderFromFile()
    
    pd_seq_list = pd_builder_from_file.ReadSeqFileList(genome_center_file)
    #print(pd_seq_list.loc[pd_seq_list['RunDate[4]'].astype(str) == '200619',['Sample[2]','RunDate[4]']])
    #il y a un typo dans la run /genfs/projects/COVID_full_processing/illumina/200619_M03555_0529_SeqWell_testing
    pd_seq_list.loc[pd_seq_list['RunDate[4]'].astype(str) == '200619',['RunDate[4]']] = '20200619'
    pd_seq_list['RunDate[4]'] = pd.to_datetime(pd_seq_list['RunDate[4]'],format='%Y%m%d')

    pd_envoi_qenome_quebec = None
    pd_sgil_extract = None

    if metadata_destination == 'GC':
        pd_envoi_qenome_quebec =  pd_builder_from_file.ReadEnvoiGenomeQuebecFile(liste_envoi_genome_quebec)
        pd_sgil_extract = pd_builder_from_file.ReadSgilExtract(sgil_extract)
     
    pd_metadata,pd_missing_spec =  CreateMetadata(pd_seq_list,metadata_destination,extract_all_samples,pd_envoi_qenome_quebec,pd_sgil_extract)

    pd_writer.WritePdToFile(pd_metadata)
    pd_writer.WritePdMissingSamplesToFile(pd_missing_spec)

    fasta_getter = FastaGetter(pd_seq_list,pd_metadata)
    fasta_getter.GetFastaFromBeluga()

if __name__ == '__main__':
    Main()
    logging.info("Terminé")
