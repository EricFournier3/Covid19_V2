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

python GetFastaForNextstrain_v2.py -i /home/foueri01@inspq.qc.ca/temp/TEST_GET_FASTA_FOR_NEXTSTRAIN/IN/updateListOfRuns_v2_Dev_2020-10-20_consensusList_small.list --debug   --maxsampledate 2020-10-20  --minsampledate 2020-05-15 --keep 'FLAG' --pangolin

"""

base_dir_dsp_db = "/data/Databases/CovBanQ_Epi/"

global base_dir


pd.options.display.max_columns = 100
logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser(description="Download Beluga consensus and Create metadata file")
parser.add_argument('--debug',help="run in debug mode",action='store_true')
parser.add_argument('--input','-i',help="Beluga fasta list file",required=True)
parser.add_argument('--onlymetadata',help="Only produce metadata",action='store_true')
parser.add_argument('--pangolin',help="Produce pangolin mapping file for LNM",action='store_true')
parser.add_argument('--maxsampledate',help="Maximum sampled date YYYY-MM-DD",required=True)
parser.add_argument('--minsampledate',help="Minimum sampled date YYYY-MM-DD",required=True)
parser.add_argument('--keep',help="Minimum Qc status to keep",choices=['ALL','PASS','FLAG'],required=True)
args = parser.parse_args()

_debug_ = args.debug
beluga_fasta_file  = args.input
only_metadata = args.onlymetadata
max_sample_date = args.maxsampledate
min_sample_date = args.minsampledate
pangolin = args.pangolin


global qc_keep

qc_keep = args.keep


try:
    max_sample_date = datetime.datetime.strptime(max_sample_date,"%Y-%m-%d")
except:
    logging.error("format de date incorrect pour maxsampledate: YYYY-MM-DD")
    exit(1)

try:
    min_sample_date = datetime.datetime.strptime(min_sample_date,"%Y-%m-%d")
except:
    logging.error("format de date incorrect pour minsampledate: YYYY-MM-DD")
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
global pangolin_outdir


if _debug_:
    #sgil_extract = "/home/foueri01@inspq.qc.ca/temp/20200924/extract_with_Covid19_extraction_v2_20200923_CovidPos_small.txt"
    #liste_envoi_genome_quebec = os.path.join(base_dir_dsp_db,"/home/foueri01@inspq.qc.ca/temp/20200924/ListeEnvoisGenomeQuebec_small6.xlsx")
    #liste_envoi_genome_quebec = os.path.join(base_dir_dsp_db,"LISTE_ENVOIS_GENOME_QUEBEC","EnvoiSmall.xlsx")
    sgil_extract = os.path.join(base_dir_dsp_db,"SGIL_EXTRACT","extract_with_Covid19_extraction_v2_20201110_CovidPos.txt")
    liste_envoi_genome_quebec = os.path.join(base_dir_dsp_db,"LISTE_ENVOIS_GENOME_QUEBEC","ListeEnvoisGenomeQuebec_2020-11-06.xlsx")
else:
    sgil_extract = os.path.join(base_dir_dsp_db,"SGIL_EXTRACT","extract_with_Covid19_extraction_v2_20201110_CovidPos.txt")
    liste_envoi_genome_quebec = os.path.join(base_dir_dsp_db,"LISTE_ENVOIS_GENOME_QUEBEC","ListeEnvoisGenomeQuebec_2020-11-06.xlsx")

#TODO CHANGER base_dir pour production
base_dir = "/home/foueri01@inspq.qc.ca/temp/TEST_GET_FASTA_FOR_NEXTSTRAIN/" if  _debug_ else "/data/PROJETS/COVID-19_Beluga/Metadata"
in_dir = os.path.join(base_dir,"IN")
out_dir = os.path.join(base_dir,"OUT")

fasta_outdir = os.path.join(out_dir,'FASTA')
metadata_outdir = os.path.join(out_dir,'METADATA')
pangolin_outdir = os.path.join(out_dir,'PANGOLIN')


if not os.path.isfile(os.path.join(in_dir,beluga_fasta_file)):
    logging.error("File missing : " + os.path.join(in_dir,beluga_fasta_file))
    exit(1)


##################################################################################
def MountBelugaServer():
    logging.info("Try to mount Beluga")
    os.system("sudo umount " + mnt_beluga_server)
    os.system("sudo sshfs -o allow_other -o follow_symlinks {0} {1}".format(beluga_server,mnt_beluga_server))
    logging.info("Beluga mounted")

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
        #print(list(pd_missing_samples))
        res_df = self.pd_sgil_extract.loc[self.pd_sgil_extract['NUMERO_SGIL'].isin(list(pd_missing_samples)),['NUMERO_SGIL','SAMPLED_DATE','TRAVEL_HISTORY','CT','SEX','RSS_PATIENT','DATE_NAISS','POSTAL_CODE']].copy()
        res_df = res_df.rename(columns={'NUMERO_SGIL':'sample','SAMPLED_DATE':'sample_date','TRAVEL_HISTORY':'country_exposure','CT':'ct','SEX':'sex','RSS_PATIENT':'rss','DATE_NAISS':'date_naiss','POSTAL_CODE':'rta'})
        res_df['country'] = 'Canada'
        res_df['division'] = 'Quebec'
        res_df['rta'] = res_df['rta'].str.slice(0,3)
        #print(res_df)
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
        res_df['rta'] = 'G8P'

        return(res_df)

    def CreateMetadata(self,max_sample_date):
        MySQLcovid19.SetConnection()

        year_2020 = datetime.datetime.strptime("2020","%Y")

        self.pd_metadata = MySQLcovid19Selector.GetMetadataAsPdDataFrame(MySQLcovid19.GetConnection(),self.samples_list,'LSPQ',False)
        #print(self.pd_metadata)
        #self.pd_metadata['sample'] = self.pd_metadata['sample'].str.replace('LSPQ-','')

        self.pd_metadata['sample'] = self.pd_metadata['sample'].str.strip(' ')

        pd_missing_samples = self.CheckMissingSpec(self.pd_metadata,self.samples_list)

        pd_missing_get_from_sgil_extract = self.AddMissingFromSgilExtract(pd_missing_samples,self.pd_sgil_extract,self.pd_metadata.columns)
        self.pd_metadata = pd.concat([self.pd_metadata,pd_missing_get_from_sgil_extract])
        
        #self.pd_metadata.to_csv("/home/foueri01@inspq.qc.ca/temp/20201111/test.tsv",sep="\t",index=False)
        #self.pd_metadata['sample'] = self.pd_metadata['sample'].str.replace('LSPQ-','')
        self.pd_metadata['sample'] = self.pd_metadata['sample'].str.strip(' ')
        pd_missing_samples = self.CheckMissingSpec(self.pd_metadata,self.samples_list)

        pd_missing_get_from_EnvoisGenomeQuebec = self.AddMissingFromEnvoisGenomeQuebec(pd_missing_samples,self.pd_envoi_qenome_quebec,self.pd_metadata.columns)
        self.pd_metadata = pd.concat([self.pd_metadata,pd_missing_get_from_EnvoisGenomeQuebec])

        #self.pd_metadata['sample'] = self.pd_metadata['sample'].str.replace('LSPQ-','')
        self.pd_metadata['sample'] = self.pd_metadata['sample'].str.strip(' ')
        self.pd_missing_samples = self.CheckMissingSpec(self.pd_metadata,self.samples_list)

        self.pd_metadata['sample_date'] = pd.to_datetime(self.pd_metadata.sample_date)
        self.pd_metadata['sample_date'] = self.pd_metadata['sample_date'].dt.strftime('%Y-%m-%d')

        self.pd_metadata = self.pd_metadata.drop_duplicates(subset='sample',keep='first')

        self.pd_metadata.reset_index(drop=True,inplace=True) 

        self.pd_samples_missing_rss = self.pd_metadata.loc[self.pd_metadata['rss'] == 'INDETERMINE',['sample']]
        self.pd_metadata = self.pd_metadata.loc[self.pd_metadata['rss'] != 'INDETERMINE',:]
        #print(self.pd_metadata.index)
        #print(self.pd_metadata[self.pd_metadata.index.duplicated()])
        self.pd_metadata.loc[self.pd_metadata['sample'].str.contains('HGA-'),'sample']= self.pd_metadata['sample'] + '2D'
        self.pd_metadata = self.pd_metadata.sort_values(by=['sample'])

        self.pd_metadata['sample_date']= self.pd_metadata['sample_date'].astype('datetime64[ns]')
        #print(self.pd_metadata.dtypes)

        #self.pd_metadata = self.pd_metadata.loc[(self.pd_metadata['sample_date'] <= max_sample_date) & (self.pd_metadata['sample_date'] >= year_2020 ),:]
        self.pd_metadata = self.pd_metadata.loc[(self.pd_metadata['sample_date'] >= min_sample_date) & (self.pd_metadata['sample_date'] <= max_sample_date) & (self.pd_metadata['sample_date'] >= year_2020 ),:]
        #print(self.pd_metadata)
        #print(self.pd_samples_missing_rss)

    def GetBelugaRunsWithTargetSamples(self):
        self.pd_samples_with_run_name = pd.merge(self.pd_metadata,self.pd_fasta_list,left_on='sample',right_on='SAMPLE',how='left',indicator=True)
        self.pd_samples_with_run_name = self.pd_samples_with_run_name[['sample','sample_date','STATUS','TECHNO','RUN_NAME']]
        self.pd_samples_with_run_name = self.pd_samples_with_run_name.sort_values(by=['sample_date'],ascending=True)
        #self.pd_samples_with_run_name = self.pd_samples_with_run_name.sort_values(by=['RUN_NAME'],ascending=True)
        #print(self.pd_with_run_name)
        #self.pd_run_list = pd.DataFrame({'RUN_NAME':self.pd_samples_with_run_name['RUN_NAME'].unique()})
        self.pd_run_list = self.pd_samples_with_run_name.drop_duplicates(['RUN_NAME','TECHNO'])
        self.pd_run_list = self.pd_run_list[['RUN_NAME','TECHNO']]
        #print(self.pd_run_list)

    def GetPdMetadata(self):
        return self.pd_metadata

    def GetPdMissingSamples(self):
        return self.pd_missing_samples

    def WriteMetadata(self,beluga_fasta_file,max_sample_date,min_sample_date):

        if qc_keep  == 'FLAG':
            qc_status_suffix = "PASS_FLAG"
        else:
            qc_status_suffix = qc_keep


        today = datetime.datetime.now().strftime("%Y-%m-%d")

        file_out_suffix = "{0}_{1}_minmaxSampleDate_{2}_{3}".format(today,qc_status_suffix,min_sample_date,max_sample_date)

        self.metadata_out = os.path.join(metadata_outdir,"metadata_{0}.tsv".format(file_out_suffix))
        self.missing_samples_out = os.path.join(metadata_outdir,"missing_samples_{0}.txt".format(file_out_suffix))
        self.samples_missing_rss_out = os.path.join(metadata_outdir,"samples_missing_rss_{0}.txt".format(file_out_suffix))
        self.samples_with_run_name_out = os.path.join(metadata_outdir,"samples_runs_{0}.txt".format(file_out_suffix))
        self.runs_out = os.path.join(metadata_outdir,"runs_{0}.txt".format(file_out_suffix))

        self.pd_metadata.to_csv(self.metadata_out,sep="\t",index=False)
        self.pd_missing_samples.to_csv(self.missing_samples_out,sep="\t",index=False)
        self.pd_samples_missing_rss.to_csv(self.samples_missing_rss_out,sep="\t",index=False)
        self.pd_samples_with_run_name.to_csv(self.samples_with_run_name_out,sep="\t",index=False)
        self.pd_run_list.to_csv(self.runs_out,sep="\t",index=False)


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


    def GetFastaFromBeluga(self,beluga_fasta_file,max_sample_date,min_sample_date):

        id_pattern = r'(^Canada/Qc-)(\S+)(/\d{4})'

        for index, row in self.pd_selected_fasta.iterrows():
            fasta_path = str(row['PATH'])
            sample = str(row['SAMPLE'])
            sample_date = self.pd_metadata.loc[self.pd_metadata['sample'] == sample,'sample_date']
            sample_date = sample_date.to_numpy()[0]
            sample_date = np.datetime_as_string(sample_date,unit='D')
            month = datetime.datetime.strptime(sample_date,'%Y-%m-%d').strftime("%B")
            fasta_path = re.sub(r'/genfs/projects/',mnt_beluga_server,fasta_path)
            qc_status = str(row['STATUS'])
            logging.info("Get " + fasta_path)
            rec = SeqIO.read(fasta_path,'fasta')
            id_short = re.search(id_pattern, rec.description).group(2)
            rec.id = re.sub(id_pattern,r"\1" + str(id_short).upper() + r"\3",rec.id)
            rec.description = re.sub(id_pattern,r"\1" + str(id_short).upper() + r"\3",rec.description)
            rec.description = rec.description + " " +  fasta_path + " " + qc_status + " " + month
            self.fasta_rec_list.append(rec)

        if qc_keep  == 'FLAG':
            qc_status_suffix = "PASS_FLAG"
        else:
            qc_status_suffix = qc_keep

        today = datetime.datetime.now().strftime("%Y-%m-%d")

        file_out_suffix = "{0}_{1}_minmaxSampleDate_{2}_{3}".format(today,qc_status_suffix,min_sample_date,max_sample_date)
        
        out_fasta = os.path.join(fasta_outdir,"sequences_{0}_temp.fasta".format(file_out_suffix))
        SeqIO.write(self.fasta_rec_list,out_fasta,'fasta')
        
        self.Dedup(out_fasta)

    def Dedup(self,fasta_file):
        
        out_fasta = re.sub(r'_temp.fasta','.fasta',fasta_file)

        final_rec_list= []
        rec_dict_dedup = {}

        for rec in SeqIO.parse(fasta_file,"fasta"):
            rec_id = rec.id
            rec_desc = rec.description
            desc_pattern = r'^Canada/Qc-\S+/\d{4} seq_method:(\S+)\|assemb_method:\S+\|snv_call_method:\S+ \S+ (PASS|FLAG) \S+$'
            seq_method = re.search(desc_pattern,rec_desc).group(1)
            qc_status = re.search(desc_pattern,rec_desc).group(2)

            if not rec_id  in rec_dict_dedup:
                rec_dict_dedup[rec_id] = [qc_status,seq_method,rec]
            else:
                val_in = rec_dict_dedup[rec_id] 
                current_val = [qc_status,seq_method,rec_desc]

                if val_in[0] == "FLAG" and current_val[0] == "PASS":
                    rec_dict_dedup[rec_id] = [qc_status,seq_method,rec]                 

        for val in rec_dict_dedup.values():
            final_rec_list.append(val[2])

        SeqIO.write(final_rec_list,out_fasta,'fasta')
        os.remove(fasta_file)



class FastaListManager():
    def __init__(self,fasta_list_file):
        self.fasta_list_file = fasta_list_file
        self.pd_fasta_list = self.GetPdFastaList()
        self.BuildSamplesList()

    def GetPdFastaList(self):
        pd_df = pd.read_csv(self.fasta_list_file,sep="\t",index_col=False)
        pd_df['SAMPLE'] = pd_df['SAMPLE'].str.upper()
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


class PangolinManager():
    def __init__(self,min_sample_date,max_sample_date):
        self.min_sample_date = min_sample_date
        self.max_sample_date = max_sample_date
        self.out_map = os.path.join(pangolin_outdir,"pangolin_map_minmaxSampleDate_" + self.min_sample_date + "_" + self.max_sample_date + ".tsv")

    def CreatePangolinMapForLNM(self,pd_metadata):
        self.pd_pangolin_map = pd.DataFrame()
        self.pd_pangolin_map['taxon'] = "Canada/Qc-" +  pd_metadata['sample'] + "/" + pd_metadata['sample_date'].astype(str).str.slice(0,4)
        self.pd_pangolin_map['SampleDate'] = pd_metadata['sample_date']
        self.pd_pangolin_map.to_csv(self.out_map,sep="\t",index=False)


def Main():
    logging.info("In Main()")

    fasta_list_manager = FastaListManager(os.path.join(in_dir,beluga_fasta_file))
    samples_list = fasta_list_manager.GetSamplesList()
    metadata_manager = MetadataManager(samples_list,fasta_list_manager.GetPdFastaList())
    metadata_manager.CreateMetadata(max_sample_date)
    metadata_manager.GetBelugaRunsWithTargetSamples()
    metadata_manager.WriteMetadata(beluga_fasta_file,max_sample_date.strftime("%Y-%m-%d"),min_sample_date.strftime("%Y-%m-%d"))

    if not only_metadata and qc_keep in ['PASS','FLAG']:
        MountBelugaServer()
        fasta_getter = FastaGetter(metadata_manager.GetPdMetadata(),fasta_list_manager.GetPdFastaList())
        fasta_getter.SelectFasta()
        fasta_getter.GetFastaFromBeluga(beluga_fasta_file,max_sample_date.strftime("%Y-%m-%d"),min_sample_date.strftime("%Y-%m-%d"))

    if pangolin:
        pangolin_manager = PangolinManager(min_sample_date.strftime("%Y-%m-%d"),max_sample_date.strftime("%Y-%m-%d"))
        pangolin_manager.CreatePangolinMapForLNM(metadata_manager.GetPdMetadata())


if __name__ == '__main__':
    Main()
    logging.info("Terminé")
