# coding=utf-8

"""
Eric Fournier 2021-01-12


Command example:
python GetGenomeCenterData.py  --lspq-report /data/PROJETS/COVID-19_Beluga/LSPQ_REPORT/2021-01-08_LSPQReport_test2.tsv --mode submission --user foueri01  --gisaid-qc-metadata /data/Applications/GitScript/Covid19_V2/NextStrainFiles/data/gisaid/all/QC_20210113/gisaid_hcov-19_2021_01_13_22.tsv --beluga-run 20200814_LSPQ_GQ0001-0008_CTL  --techno illumina

TODO
- check args
- enregistrer commande dans repertoire de soumission
- ajouter option freeze1 resubmission
- PLus besoin de faire la correspondance pour numero HCLM ET AUTRE CAR BIOBANK_ID INTEGRE DANS COVBANK
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
parser.add_argument('--lspq-report','-l',help="genome center LSPQ report input file path")
parser.add_argument('--gisaid-qc-metadata','-g',help="gisaid quebec metadata file path")
parser.add_argument('--user','-u',help='user name',required=True,choices=['foueri01','morsan01'])
parser.add_argument('--minsampledate',help="Minimum sampled date YYYY-MM-DD",type=Tools.CheckDateFormat)
parser.add_argument('--maxsampledate',help="Maximum sampled date YYYY-MM-DD",type=Tools.CheckDateFormat)
parser.add_argument('--keep',help="Minimum Qc status to keep",choices=['PASS','FLAG'])
parser.add_argument('--beluga-run',help="beluga run name")
parser.add_argument('--techno',help="sequencing technologie",choices=['illumina','nanopore','mgi'])
parser.add_argument('--freeze1',help="for mgi and illumina freeze1 resubmission",action='store_true')
parser.add_argument('--cons-file',help="path to updateListOfRuns_v3_YYYY-MM-DD_consensusList.list")
parser.add_argument('--vcf-file',help="path to updateListOfRuns_v3_YYYY-MM-DD_vcfList.list")
parser.add_argument('--max-perc-n',type=int,help='Maximun percentage of N tolerated in reject consensus for nextstrain')


args = parser.parse_args()

_debug_ = args.debug

global fasta_qual_status_to_keep
if args.keep == 'FLAG':
    fasta_qual_status_to_keep = ['PASS','FLAG']
else:
    fasta_qual_status_to_keep = ['PASS']

global max_perc_n
max_perc_n = args.max_perc_n

global consensus_list_file
consensus_list_file = args.cons_file

global freeze1
freeze1 = args.freeze1

global lspq_report
lspq_report = args.lspq_report

global gisaid_qc_metadata
gisaid_qc_metadata = args.gisaid_qc_metadata

global seq_techno
seq_techno = args.techno

global beluga_run
beluga_run = args.beluga_run

global user_name
user_name = args.user

global gisaid_publication_basedir

global mode
mode = args.mode

global min_sample_date
min_sample_date = args.minsampledate

global max_sample_date
max_sample_date = args.maxsampledate


if _debug_:
    #lspq_report = '/data/PROJETS/ScriptDebug/GetGenomeCenterData/IN/2021-01-08_LSPQReport_test2.tsv' # avec 20200814_LSPQ_GQ0001-0008_CTL
    lspq_report = "/data/PROJETS/ScriptDebug/GetGenomeCenterData/IN/2021-01-08_LSPQReport.tsv"
    #lspq_report = '/data/PROJETS/ScriptDebug/GetGenomeCenterData/IN/test_HCLM.tsv' # avec 20201113_LSPQ_GQ0046-0049_CTL
    gisaid_qc_metadata =  '/data/PROJETS/ScriptDebug/GetGenomeCenterData/IN/gisaid_hcov-19_2021_01_13_22.tsv'
    consensus_list_file = "/data/PROJETS/ScriptDebug/GetGenomeCenterData/IN/updateListOfRuns_v3_2021-01-18_consensusList.list"   
    vcf_list_file = "/data/PROJETS/ScriptDebug/GetGenomeCenterData/IN/updateListOfRuns_v3_2021-01-25_vcfList.list"

 
#TODO A REACTIVER
#Tools.CheckUserArgs(gisaid_qc_metadata,mode)

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

    def GetVcfMetadataAsPdDataFrame(self):
        columns_renamed = {'GENOME_QUEBEC_REQUETE':'# Requête','NOM':'Nom','PRENOM':'Prénom','DTNAISS':'Date de naissance','DATE_PRELEV':'Date de prélèvement'}

        sql = "SELECT pr.GENOME_QUEBEC_REQUETE, pa.NOM, pa.PRENOM, pa.DTNAISS, pr.DATE_PRELEV, pa.NAM FROM Prelevements pr inner join Patients pa on pa.ID_PATIENT = pr.ID_PATIENT WHERE pr.DATE_PRELEV BETWEEN date('{0}') AND date('{1}')".format(min_sample_date,max_sample_date)

        df = pd.read_sql(sql,con=self.GetConnection())
        df = df.rename(columns=columns_renamed)
        df['# Requête'] = df['# Requête'].str.strip(' ')
        return(df)

    def GetNextstrainMetadataAsPdDataFrame(self):
        start = time.time()

        columns_renamed = {'RTA':'rta','GENOME_QUEBEC_REQUETE':'sample','DATE_PRELEV':'sample_date','TRAVEL_HISTORY':'country_exposure','CT':'ct','RSS':'rss','SEXE':'sex','COUNTRY':'country','DIVISION':'division','DTNAISS':'date_naiss','OUTBREAK':'OUTBREAK'}

        sql = "SELECT pr.GENOME_QUEBEC_REQUETE, pr.DATE_PRELEV, pr.TRAVEL_HISTORY, pr.CT, pr.OUTBREAK, pa.RSS, pa.DTNAISS,pa.SEXE,{2} as COUNTRY,{3} as DIVISION,{4} as RTA FROM Prelevements pr inner join Patients pa on pa.ID_PATIENT = pr.ID_PATIENT WHERE pr.DATE_PRELEV BETWEEN date('{0}') AND date('{1}')".format(min_sample_date,max_sample_date,"'Canada'","'Quebec'","'G8P'")
        df = pd.read_sql(sql,con=self.GetConnection())
        df = df.rename(columns=columns_renamed)
        df['sample'] = df['sample'].str.strip(' ')

        end = time.time()

        return(df)

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
        #COVV_VIRUS_NAME = "{0}.GENOME_QUEBEC_REQUETE as covv_virus_name".format(Prelevements_alias) #decoration ajouter dans CovGisaidSubmissionV2.py
        COVV_VIRUS_NAME = "{0}.COVBANK_ID as covv_virus_name".format(Prelevements_alias) #decoration ajouter dans CovGisaidSubmissionV2.py
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
        #COVV_ORIG_LAB = "{0}.NOM_HOPITAL as covv_orig_lab".format(Prelevements_alias)
        #COVV_ORIG_LAB_ADDR = "{0}.ADRESSE_HOPITAL as covv_orig_lab_addr".format(Prelevements_alias)
        COVV_ORIG_LAB = "'" + lspq_name  + "' as covv_orig_lab"
        COVV_ORIG_LAB_ADDR = "'" + lspq_addr + "' as covv_orig_lab_addr"
        COVV_PROVIDER_SAMPLE_ID = "\' \' as covv_provider_sample_id"
        COVV_SUBM_LAB = "'" + lspq_name  + "' as covv_subm_lab"
        COVV_SUBM_LAB_ADDR = "'" + lspq_addr + "' as covv_subm_lab_addr"
        #COVV_SUBM_SAMPLE_ID = "{0}.GENOME_QUEBEC_REQUETE as covv_subm_sample_id".format(Prelevements_alias)
        COVV_SUBM_SAMPLE_ID = "{0}.COVBANK_ID as covv_subm_sample_id".format(Prelevements_alias)
        COVV_AUTHORS = "'" + authors + "' as covv_authors"

        #pour case sensitive on a besoin de BINARY
        #sql = "SELECT {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25} FROM Prelevements {26} inner join Patients {27} on {27}.{28} = {26}.{28} WHERE BINARY {26}.GENOME_QUEBEC_REQUETE REGEXP '{29}'".format(SUBMITTER,FN,COVV_VIRUS_NAME,COVV_TYPE,COVV_PASSAGE,COVV_COLLECTION_DATE,COVV_LOCATION,COVV_ADD_LOCATION,COVV_HOST,COVV_ADD_HOST_INFO,COVV_GENDER,DTNAISSINFO,COVV_PATIENT_STATUS,COVV_SPECIMEN,COVV_OUTBREAK,COVV_LAST_VACCINATED,COVV_TREATMENT,COVV_ASSEMBLY_METHOD,COVV_COVERAGE,COVV_ORIG_LAB,COVV_ORIG_LAB_ADDR,COVV_PROVIDER_SAMPLE_ID,COVV_SUBM_LAB,COVV_SUBM_LAB_ADDR,COVV_SUBM_SAMPLE_ID,COVV_AUTHORS,Prelevements_alias,Patients_alias,join_column,sample_list)


        #TODO UTILISER  IN au lieu de REGEXP
        sql = "SELECT {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25} FROM Prelevements {26} inner join Patients {27} on {27}.{28} = {26}.{28} WHERE  {26}.GENOME_QUEBEC_REQUETE REGEXP '{29}'".format(SUBMITTER,FN,COVV_VIRUS_NAME,COVV_TYPE,COVV_PASSAGE,COVV_COLLECTION_DATE,COVV_LOCATION,COVV_ADD_LOCATION,COVV_HOST,COVV_ADD_HOST_INFO,COVV_GENDER,DTNAISSINFO,COVV_PATIENT_STATUS,COVV_SPECIMEN,COVV_OUTBREAK,COVV_LAST_VACCINATED,COVV_TREATMENT,COVV_ASSEMBLY_METHOD,COVV_COVERAGE,COVV_ORIG_LAB,COVV_ORIG_LAB_ADDR,COVV_PROVIDER_SAMPLE_ID,COVV_SUBM_LAB,COVV_SUBM_LAB_ADDR,COVV_SUBM_SAMPLE_ID,COVV_AUTHORS,Prelevements_alias,Patients_alias,join_column,sample_list)
        #print("SQL ",sql)
        df = pd.read_sql(sql,con=self.GetConnection())
        return(df)

class VariantDataManager:
    def __init__(self,cov_bank_db_obj):
        self.vcf_list_df = pd.read_csv(vcf_list_file,sep="\t",index_col=False)
        self.cov_bank_db_obj = cov_bank_db_obj
        
        self.spec_mut_dict = {}
        self.mut_set = set([])

        if _debug_:
            self.base_dir_out = "/data/PROJETS/ScriptDebug/GetGenomeCenterData/OUT/Variant/"
        else:
            self.base_dir_out = "/data/PROJETS/ScriptDebug/GetGenomeCenterData/OUT/Variant/"

        self.parsed_outdir = os.path.join(self.base_dir_out,"PARSED/")
       
        self.SelectVcfToDownload()
 
    def GetNumericStatus(self,status):
        if status == 'PASS':
            return(3)
        elif status == 'FLAG':
            return(2)
        else:
            return(1)

    def SelectVcfToDownload(self):

        self.vcf_to_download_df = self.vcf_list_df.loc[self.vcf_list_df['STATUS'].isin(fasta_qual_status_to_keep),:]
        self.vcf_to_download_df = self.vcf_to_download_df.copy()

        self.vcf_to_download_df['STATUS_VAL'] = self.vcf_to_download_df['STATUS'].apply(self.GetNumericStatus)
        idx  = self.vcf_to_download_df.groupby(['SAMPLE','STATUS'])['PERC_N'].transform(min) == self.vcf_to_download_df['PERC_N']
        self.vcf_to_download_df = self.vcf_to_download_df[idx]
        idx  = self.vcf_to_download_df.groupby(['SAMPLE'])['STATUS_VAL'].transform(max) == self.vcf_to_download_df['STATUS_VAL']
        self.vcf_to_download_df = self.vcf_to_download_df[idx]
        idx  = self.vcf_to_download_df.groupby(['SAMPLE'])['TECHNO'].transform(min) == self.vcf_to_download_df['TECHNO']
        self.vcf_to_download_df = self.vcf_to_download_df[idx]

    def PrepareMetadataBasedOnSampledDate(self):
        self.metadata = self.cov_bank_db_obj.GetVcfMetadataAsPdDataFrame()

    def MergeVcfToDownloadWithMetadata(self):
        self.vcf_to_download_df =  pd.merge(self.vcf_to_download_df,self.metadata,left_on='SAMPLE',right_on='# Requête',how='inner')

    def GetVcfFromBeluga(self):
        for index,row in self.vcf_to_download_df.iterrows():
            sample = str(row['SAMPLE'])
            path = str(row['PATH'])
            nom = str(row['Nom'])
            prenom = str(row['Prénom'])
            dt_naiss = str(row['Date de naissance'])
            nam = str(row['NAM'])
            dt_prelev = str(row['Date de prélèvement'])

            vcf = GenomeCenterConnector.GetVcfPath(path)
            logging.info("Get " + vcf)
            shutil.copy(vcf,os.path.join(self.base_dir_out,os.path.basename(vcf)))

            self.UpdateSpecMutDict(os.path.join(self.base_dir_out,os.path.basename(vcf)))
            self.spec_mut_dict[sample]['# Requête'] = sample
            self.spec_mut_dict[sample]['Nom'] = nom
            self.spec_mut_dict[sample]['Prénom'] = prenom
            self.spec_mut_dict[sample]['Date de naissance'] = dt_naiss
            self.spec_mut_dict[sample]['Date de prélèvement'] = dt_prelev
            self.spec_mut_dict[sample]['NAM'] = nam
            
    def CompileMutation(self):
        for spec,spec_info in self.spec_mut_dict.items():
            for mut in spec_info['mutlist']:
                self.mut_set.add(mut)

    def CreateSummaryFile(self):
        
        mut_df_columns = ['# Requête','Nom','Prénom','Date de naissance','Date de prélèvement','NAM',"AA_MUTATIONS"]
        mut_df_columns.extend(list(self.mut_set))
        
        mut_df = pd.DataFrame(columns=mut_df_columns)

        for spec,spec_info in self.spec_mut_dict.items():
            my_mut_dict = {}
            spec_mut_list = spec_info['mutlist']
            for uniq_mut in list(self.mut_set):
                if uniq_mut in spec_mut_list:
                    my_mut_dict[uniq_mut] = "1"
                else:
                    my_mut_dict[uniq_mut] = "0"
            my_mut_dict.update({'# Requête':spec,'Nom':spec_info['Nom'],'Prénom':spec_info['Prénom'],'Date de naissance':spec_info['Date de naissance'],'Date de prélèvement':spec_info['Date de prélèvement'],'NAM':spec_info['NAM'],"AA_MUTATIONS":str(spec_mut_list)})
        
            mut_df = mut_df.append(my_mut_dict,ignore_index=True)

        if fasta_qual_status_to_keep  == ['PASS','FLAG']:
            qc_status_suffix = "PASS_FLAG"
        else:
            qc_status_suffix = "PASS"

        today = datetime.datetime.now().strftime("%Y-%m-%d")

        min_sample_date_s = min_sample_date.strftime("%Y-%m-%d")
        max_sample_date_s = max_sample_date.strftime("%Y-%m-%d")

        file_out_suffix = "{0}_{1}_minmaxSampleDate_{2}_{3}".format(today,qc_status_suffix,min_sample_date_s,max_sample_date_s)

        out_summary = os.path.join(self.parsed_outdir,"variant_{0}.xlsx".format(file_out_suffix))

        mut_df.to_excel(out_summary,sheet_name='Sheet1',index=False)

    def UpdateSpecMutDict(self,vcf_file):
        spec = os.path.basename(vcf_file).split('.')[0].split('_')[0]
        self.spec_mut_dict[spec] = {}
        self.spec_mut_dict[spec]['mutlist'] = []
        header_read = False
        is_header = {}
        is_nanopore = False
        line_nb = 0

        with open(vcf_file) as vcff:
            for line in vcff:
                line_nb += 1
                is_header = False

                if (line == 2) and (not re.search(r'^##source=ivar',line)):
                    is_nanopore = True 

                if re.search(r'^#CHROM\t',line):
                    header_read = True
                    is_header = True

                if header_read and not is_header:
                    aa_change_set = set([])
                    split_line = line.split('\t')
                    nuc_change_qual = split_line[6]

                    if nuc_change_qual != "PASS":
                        continue

                    aa_change_annotations = split_line[7]

                    aa_change_annotations = re.search(r'\S+ANN=(\S+)',aa_change_annotations).group(1)
                    aa_change_first_annotation = aa_change_annotations.split(',')[0]
                    aa_change_first_annotation = aa_change_first_annotation.split('|')

                    if aa_change_first_annotation[1] == "missense_variant":
                        gene_name = aa_change_first_annotation[3]
                        aa_mut = aa_change_first_annotation[10][2:]
                        self.spec_mut_dict[spec]['mutlist'].append(gene_name + ":" + aa_mut)

class NextstrainDataManager:
    def __init__(self,cov_bank_db_obj):
        self.consensus_list_df = pd.read_csv(consensus_list_file,sep="\t",index_col=False)
        self.cov_bank_db_obj = cov_bank_db_obj

        if _debug_:
            self.base_dir_out = "/data/PROJETS/ScriptDebug/GetGenomeCenterData/OUT/Nextstrain/"
            self.base_dir_out_metadata = os.path.join(self.base_dir_out,"metadata")
            self.base_dir_out_consensus = os.path.join(self.base_dir_out,"consensus")
            self.base_dir_out_pangolin = os.path.join(self.base_dir_out,"pangolin")
        else:
            self.base_dir_out = "/data/PROJETS/COVID-19_Beluga/"
            self.base_dir_out_metadata = os.path.join(self.base_dir_out,"Metadata/")
            self.base_dir_out_consensus = os.path.join(self.base_dir_out,"Consensus/Fasta/")
            self.base_dir_out_pangolin = os.path.join(self.base_dir_out,"Pangolin/")


        self.SelectConsensusToDownload()

    def GetNumericStatus(self,status):
        if status == 'PASS':
            return(3)
        elif status == 'FLAG':
            return(2)
        else:
            return(1)

    def SelectConsensusToDownload(self):

        self.consensus_to_download_df = self.consensus_list_df.loc[self.consensus_list_df['STATUS'].isin(fasta_qual_status_to_keep),:]
        self.consensus_to_download_df = self.consensus_to_download_df.copy()

        if max_perc_n:
            self.pd_df_rej_samples_tolerated = self.consensus_list_df.loc[(self.consensus_list_df['STATUS'] == 'REJ' ) & (self.consensus_list_df['PERC_N'] <= max_perc_n)]
            self.consensus_to_download_df = pd.concat([self.consensus_to_download_df,self.pd_df_rej_samples_tolerated])
           
        self.consensus_to_download_df['STATUS_VAL'] = self.consensus_to_download_df['STATUS'].apply(self.GetNumericStatus) 
 
        idx  = self.consensus_to_download_df.groupby(['SAMPLE','STATUS'])['PERC_N'].transform(min) == self.consensus_to_download_df['PERC_N']
        self.consensus_to_download_df = self.consensus_to_download_df[idx]

        idx  = self.consensus_to_download_df.groupby(['SAMPLE'])['STATUS_VAL'].transform(max) == self.consensus_to_download_df['STATUS_VAL']
        self.consensus_to_download_df = self.consensus_to_download_df[idx]       

        idx  = self.consensus_to_download_df.groupby(['SAMPLE'])['TECHNO'].transform(min) == self.consensus_to_download_df['TECHNO']
        self.consensus_to_download_df = self.consensus_to_download_df[idx]

    def PrepareMetadataBasedOnSampledDate(self):
        self.metadata = self.cov_bank_db_obj.GetNextstrainMetadataAsPdDataFrame()


    def MergeConsensusToDownloadWithMetadata(self):
        consensus_columns = list(self.consensus_to_download_df.columns).copy()
        consensus_columns.append('sample_date')

        metadata_columns = list(self.metadata.columns).copy()

        self.consensus_to_download_df = pd.merge(self.consensus_to_download_df,self.metadata,left_on='SAMPLE',right_on='sample',how='inner')
        self.final_metadata = self.consensus_to_download_df[metadata_columns]
        self.consensus_to_download_df = self.consensus_to_download_df[consensus_columns]

        self.final_metadata['fasta_id'] = "temp"

    def GetConsensusFromBeluga(self):
        id_pattern = r'(^Canada/Qc-)(\S+)(/\d{4})'
        self.rec_list = []

        for index,row in self.consensus_to_download_df.iterrows():
            sample = str(row['SAMPLE'])
            path = str(row['PATH'])
            techno = str(row['TECHNO'])
            run_name = str(row['RUN_NAME'])
            qc_status = str(row['STATUS'])
            sample_date = row['sample_date']
            month = sample_date.strftime("%B")
            consensus = GenomeCenterConnector.GetConsensusPath(techno,sample,run_name)
            if consensus is None:
                msg = "No consensus found for " + sample
                logging.warning(msg)
                continue
            
            logging.info("Get " + consensus)            
            rec =  SeqIO.read(consensus,'fasta')
            id_short = re.search(id_pattern, rec.description).group(2)
            rec.id = re.sub(id_pattern,r"\1" + str(id_short) + r"\3",rec.id)
            rec.description = re.sub(id_pattern,r"\1" + str(id_short) + r"\3",rec.description)
            rec.description = rec.description + " " +  consensus + " " + qc_status + " " + month
            self.rec_list.append(rec)

            self.final_metadata.loc[self.final_metadata['sample'] == sample,['fasta_id']] = rec.id

        if fasta_qual_status_to_keep  == ['PASS','FLAG']:
            qc_status_suffix = "PASS_FLAG"
        else:
            qc_status_suffix = "PASS"

        today = datetime.datetime.now().strftime("%Y-%m-%d")
           
        min_sample_date_s = min_sample_date.strftime("%Y-%m-%d") 
        max_sample_date_s = max_sample_date.strftime("%Y-%m-%d") 
        file_out_suffix = "{0}_{1}_minmaxSampleDate_{2}_{3}".format(today,qc_status_suffix,min_sample_date_s,max_sample_date_s)
        out_consensus =  os.path.join(self.base_dir_out_consensus,"sequences_{0}.fasta".format(file_out_suffix))
        SeqIO.write(self.rec_list,out_consensus,'fasta')

        self.final_metadata.loc[(self.final_metadata['OUTBREAK'].isnull()) | (self.final_metadata['OUTBREAK'] == 'NA'),['OUTBREAK']] = 'NoOutbreakRelated'
        self.final_metadata.to_csv(os.path.join(self.base_dir_out_metadata,"metadata_{0}.tsv".format(file_out_suffix)),sep="\t",index=False)

        self.pangolin_map_df = self.final_metadata[['fasta_id','sample_date']] 
        self.pangolin_map_df = self.pangolin_map_df.rename(columns={'fasta_id':'taxon','sample_date':'SampleDate'}) 
        out_pangolin = os.path.join(self.base_dir_out_pangolin,"pangolin_map_{0}".format(file_out_suffix))
        self.pangolin_map_df.to_csv(out_pangolin,sep="\t",index=False)
 
class GisaidSubmissionTraceLogger:
    def __init__(self,output_dir):
        self.output_dir = output_dir

        self.log_out = os.path.join(self.output_dir,'log.txt')
        self.missing_consensus_out = os.path.join(self.output_dir,'MissingConsensus.txt') 

        self.SetLogHandler()
        self.SetMissingConsensusHandler()

    def SetLogHandler(self):
        self.log_handler = open(self.log_out,'w')

    def SetMissingConsensusHandler(self):
        self.missing_consensus_handler = open(self.missing_consensus_out,'w')

    def GetLogHandler(self):
        return(self.log_handler)

    def GetMissingConsensusHandler(self):
        return(self.missing_consensus_handler)

    def CloseLogHandler(self):
         self.log_handler.close()

    def CloseMissingConsensusHandler(self):
         self.missing_consensus_handler.close()


class GisaidFreeze1DataSubmissionManager:
    def __init__(self,lspq_report,techno,cov_bank_db_obj):
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

        self.logger = GisaidSubmissionTraceLogger(self.submission_dir)
        
    def CloseLogger(self):
        self.logger.CloseLogHandler()
        self.logger.CloseMissingConsensusHandler()

    def CreateSubmissionDirectory(self):
        if _debug_:
            self.gisaid_publication_basedir = '/data/PROJETS/ScriptDebug/GetGenomeCenterData/OUT/Gisaid/'
        else:
            self.gisaid_publication_basedir = '/data/PROJETS/COVID-19_Beluga/Gisaid/'

        self.submission_dir = os.path.join(self.gisaid_publication_basedir,self.techno,self.run_name,self.today)
        
        try:
            os.makedirs(self.submission_dir,exist_ok=True)
        except:
            logging.error("Impossible de creer " + self.submission_dir)
        
    def SetDataSubmissionDf(self):
        self.input_file_df = pd.read_csv(self.input_file,sep="\t",index_col=False)
        self.input_file_df = self.input_file_df.loc[self.input_file_df['run_name'] == beluga_run,: ]
        self.input_file_df['ncov_tools.pass'] = self.input_file_df['ncov_tools.pass'].astype(str)

        self.gisaid_metadata_df = pd.read_csv(self.gisaid_metadata,sep="\t",index_col=False)
        self.gisaid_metadata_df = self.gisaid_metadata_df.loc[self.gisaid_metadata_df['Virus name'].str.contains('^hCoV-19/Canada/Qc-\S+/\S+',flags=re.IGNORECASE,regex=True),:]
        self.gisaid_metadata_df['Virus name'] = self.gisaid_metadata_df['Virus name'].str.replace(r'^hCoV-19/Canada/Qc-(\S+)/\S+',r'\1',flags=re.IGNORECASE,regex=True)
        self.gisaid_metadata_qc_list = list(self.gisaid_metadata_df['Virus name'].str.replace(r'LSPQ-',r'',flags=re.IGNORECASE,regex=True))

        self.input_file_df = self.input_file_df.loc[(~self.input_file_df['Sample Name'].isin(self.gisaid_metadata_qc_list)) & (self.input_file_df['PASS/FLAG/REJ'].str.upper() == 'PASS') & (self.input_file_df['ncov_tools.pass'].str.upper() ==  'TRUE'),['Sample Name','PASS/FLAG/REJ','run_name','ncov_tools.pass','platform']]

    def GetConsensus(self):
        self.rec_list =  [] 
        for index,row in self.input_file_df.iterrows():
            sample_name = row['Sample Name']
            run_name = row['run_name']
            
            consensus = GenomeCenterConnector.GetConsensusPath(self.techno,sample_name,run_name)
            if consensus is None:
                msg = "No consensus found for " + sample_name
                logging.warning(msg)
                self.logger.GetMissingConsensusHandler().write(msg + "\n")
                continue

            rec = SeqIO.read(consensus,'fasta')
            try:
                msg = "Try to get " + rec.id
                logging.info(msg)
                self.logger.GetLogHandler().write(msg + "\n")
                parsed_header = re.search(r'(Canada/Qc-)(\S+)/(\d{4}) seq_method:(\S+)\|assemb_method:\S+\|snv_call_method:\S+',rec.description)
                method = parsed_header.group(4)
                self.sample_to_submit_dict[rec.id] = {}
                self.sample_to_submit_dict[rec.id]['method'] = method
                self.sample_to_submit_dict[rec.id]['gisaid_id'] = 'hCoV-19/' + rec.id
                rec.description = ""
                rec.id = self.sample_to_submit_dict[rec.id]['gisaid_id']
                self.rec_list.append(rec)
            except:
                msg = "Unable to get " + rec.description
                logging.warning(msg)
                self.logger.GetMissingConsensusHandler().write(msg + "\n")

    def SaveConsensus(self):
        consensus_out = os.path.join(self.submission_dir,"all_sequences.fasta")
        final_rec_list = []
        for index, row in self.metadata_df.iterrows():
            gisaid_id = row['covv_virus_name']
            for rec in self.rec_list:
                if rec.id == gisaid_id:
                    final_rec_list.append(rec)

        SeqIO.write(final_rec_list,consensus_out,'fasta')

    def GetCovBankId(self,biobank_id):
        covbank_id = biobank_id.upper()

        if (re.search(r'^HCLM-',covbank_id)) and (len(covbank_id) == 15) :
            covbank_id = re.sub(r'(^HCLM-)\S{3}(\S+)',r'\1\2',biobank_id)

        if (re.search(r'^JUS-',covbank_id)) and (len(covbank_id) == 14):
            covbank_id = biobank_id[:-2]

        return(covbank_id)

    def BuildBioBankToCovBankDict(self):
       self.biobankId_2_covbankId_dict = {}
       for biobank_id in self.sample_list:
           self.biobankId_2_covbankId_dict[biobank_id] = self.GetCovBankId(biobank_id)

    def BuildSampleList(self):
        self.sample_list = []
        for sample in self.sample_to_submit_dict:
            short_sample_name = re.search(r'Canada/Qc-(\S+)/\d+',sample).group(1)
            self.sample_list.append(short_sample_name)

    def GetBiobankIdFromCovbankID(self,covbank_id_from_db):
        covbank_id_from_db = covbank_id_from_db.upper()
        for biobank_id,covbank_id in self.biobankId_2_covbankId_dict.items():
            if covbank_id_from_db == covbank_id:
                return(biobank_id)

    def CreateMetadata(self):
        metadata_out = os.path.join(self.submission_dir,"{0}_ncov19_metadata.xls".format(self.today))
        self.metadata_df = self.cov_bank_db_obj.GetGisaidMetadataAsPdDataFrame(self.biobankId_2_covbankId_dict.values()) 
        self.metadata_df['covv_subm_sample_id'] = self.metadata_df['covv_subm_sample_id'].apply(self.GetBiobankIdFromCovbankID)
        self.CheckMissingMetadata(self.metadata_df['covv_subm_sample_id'])
        self.metadata_df['covv_virus_name'] = self.metadata_df['covv_virus_name'].apply(self.GetBiobankIdFromCovbankID)
        self.metadata_df['covv_virus_name'] = self.metadata_df['covv_virus_name'].apply(self.GetVirusName)
        self.metadata_df.insert(loc=11,column='covv_patient_age',value=self.metadata_df['DTNAISS'].apply(lambda x: self.from_dob_to_age(x)))
        self.metadata_df.insert(loc=18,column='covv_seq_technology',value=self.metadata_df['covv_virus_name'].apply(self.GetSequencingMethod))
        self.metadata_df['covv_collection_date'] = self.metadata_df['covv_collection_date'].astype(str) 
        del self.metadata_df['DTNAISS']

        added_header = pd.DataFrame({'submitter':['Submitter'],'fn':['FASTA filename'],'covv_virus_name':['Virus name'],'covv_type':['Type'],'covv_passage':['Passage details/history'],'covv_collection_date':['Collection date'],'covv_location':['Location'],'covv_add_location':['Additionnal location information'],'covv_host':['Host'],'covv_add_host_info':['Additional host info'], 'covv_gender':['Gender'],'covv_patient_age':['Patient age'],'covv_patient_status':['Patient status'],'covv_specimen':['Specimen source'],'covv_outbreak':['Outbreak'],'covv_last_vaccinated':['Last vaccinated'],'covv_treatment':['Treatment'],'covv_seq_technology':['Sequencing technology'],'covv_assembly_method':['Assembly method'],'covv_coverage':['Coverage'],'covv_orig_lab':['Originating lab'],'covv_orig_lab_addr':['Address'],'covv_provider_sample_id':['Sample ID given by the sample provider'],'covv_subm_lab':['Submitting lab'],'covv_subm_lab_addr':['Address'],'covv_subm_sample_id':['Sample ID given by the submitting laboratory'],'covv_authors':['Authors']})

        self.metadata_df = pd.concat([added_header,self.metadata_df])
        self.metadata_df.to_excel(metadata_out,index=False,sheet_name='Submission')

    def GetSequencingMethod(self,sample_name):
        try:
            sample_name_ori = sample_name
            sample_name = re.search(r'hCoV-19/(Canada/Qc-(\S+)/\d+)',sample_name).group(1)
            return(self.sample_to_submit_dict[sample_name]['method'])
        except:
            logging.error("Problem get sequencing method for " + sample_name_ori)
            

    def from_dob_to_age(self,born):
        today = datetime.date.today()
        return today.year - born.year - ((today.month, today.day) < (born.month, born.day))

    def GetVirusName(self,req_number):
        for rec_id,d in self.sample_to_submit_dict.items():
            sample_name = re.search(r'^Canada/Qc-(\S+)/\d+',rec_id).group(1) 
            if sample_name == req_number:
                return(d['gisaid_id']) 

    def CheckMissingMetadata(self,sample_in_metadata):
        missing = set(self.sample_list) - set(sample_in_metadata) 
        df = pd.DataFrame({'MissingSample':list(missing)})
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
            return(None)

    @staticmethod
    def GetVcfPath(raw_beluga_path):
        mnt_beluga_path = re.sub(r'/genfs/projects/',GenomeCenterConnector.mnt_beluga_server,raw_beluga_path)
        return(mnt_beluga_path)

def Main():
    logging.info('Begin')

    cov_bank_db_obj = CovBankDB()

    if mode == 'nextstrain':
        nextstrain_data_manager = NextstrainDataManager(cov_bank_db_obj)
        nextstrain_data_manager.PrepareMetadataBasedOnSampledDate()
        nextstrain_data_manager.MergeConsensusToDownloadWithMetadata()

        GenomeCenterConnector.MountBelugaServer()
        nextstrain_data_manager.GetConsensusFromBeluga()
        cov_bank_db_obj.CloseConnection()

    elif mode == "variant":
        variant_data_manager = VariantDataManager(cov_bank_db_obj)
        variant_data_manager.PrepareMetadataBasedOnSampledDate()
        variant_data_manager.MergeVcfToDownloadWithMetadata()

        GenomeCenterConnector.MountBelugaServer()
        variant_data_manager.GetVcfFromBeluga()
        variant_data_manager.CompileMutation()
        variant_data_manager.CreateSummaryFile()
        cov_bank_db_obj.CloseConnection()


    elif mode == 'submission':
        if not freeze1:
            GenomeCenterConnector.MountBelugaServer()
            gisaid_data_submission_manager = GisaidDataSubmissionManager(lspq_report,gisaid_qc_metadata,beluga_run,seq_techno,cov_bank_db_obj)

            gisaid_data_submission_manager.GetConsensus()
            gisaid_data_submission_manager.BuildSampleList()
            gisaid_data_submission_manager.BuildBioBankToCovBankDict()
            gisaid_data_submission_manager.CreateMetadata()
            gisaid_data_submission_manager.SaveConsensus()

            gisaid_data_submission_manager.CloseLogger()
            cov_bank_db_obj.CloseConnection()

        else:
            logging.info("Freeze1 {0} resubmission".format(seq_techno))

if __name__ == '__main__':
    Main()
    logging.info("Terminé")

