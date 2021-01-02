# coding=utf-8

import pandas as pd
import logging
import numpy as np
import glob
import argparse
import os
import shutil
import re
import time
import sys
from datetime import datetime

'''
Created on 29 December 2020
@author: Eric Fournier

Note:
    - pour executer ce script, il faut au prealable executer <module load scipy-stack>


TODO:
    - combiner updateSampleList_v2.py et updateListOfRuns_v2.py
    - modifier le README
    - ne pas oublier le NOTFOUND dans REPOSITORY
    - corriger le logging stdout
'''

pd.options.display.max_columns = 100
logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser(description="Parse sequencing runs")

parser.add_argument('--debug',help='run in debug mode',action='store_true')
parser.add_argument('--input',help='listPlateSampleYYYYMMDD.tsv file name',required=True)

args = parser.parse_args()

global _debug
_debug = args.debug

global input_file_name
input_file_name = args.input

global full_processing_path
global repository_path
global trace_path

if _debug:
    debug_base_dir = "/genfs/projects/COVID_LSPQ/DEBUG/"
    full_processing_path = os.path.join(debug_base_dir,"COVID_full_processing/")
    repository_path = os.path.join(debug_base_dir,"REPOSITORY/")
    trace_path = os.path.join(debug_base_dir,"TRACE/")
else:
    full_processing_path = "/genfs/projects/COVID_full_processing/"
    repository_path = "/genfs/projects/COVID_LSPQ/REPOSITORY/"
    trace_path = "/genfs/projects/COVID_LSPQ/TRACE/"

global illumina_base_dir
illumina_base_dir = os.path.join(full_processing_path,"illumina_reprocess")

global mgi_base_dir
mgi_base_dir = os.path.join(full_processing_path,"mgi_reprocess")

global nanopore_base_dir
nanopore_base_dir = os.path.join(full_processing_path,"nanopore")

global input_file_path
input_file_path = os.path.join(trace_path,input_file_name)

class Run():
    def __init__(self,run_name,run_path,techno):
        self.run_name = run_name
        self.techno = techno
        self.run_path = run_path

        self.SetPath()
        self.SetSamplesObjList()

    def GetRunPath(self):
        return(self.run_path)

    def GetRunName(self):
        return(self.run_name)

    def GetRunDate(self):
            return(self.run_name.split('_')[-1] if self.techno == 'nanopore' else self.run_name.split('_')[0])

    def GetRunTechno(self):
        return(self.techno.lower())

    def GetSamplesWithoutPlate(self):
        list_samples_without_plate = []

        for sample in self.samples_obj_list:
            if not sample.GetIsInPlate():
                list_samples_without_plate.append(sample.GetSampleName())

        return(list_samples_without_plate)

    def SetSamplesObjList(self):
        self.samples_obj_list = []

        consensus_suffix = ".consensus.{0}.*.fasta".format(self.techno)
        clean_raw_reads_suffix = "*trim.pair*.fastq.gz"
        host_removal_suffix = "*host_removed.pair*.fastq.gz"
        metrics_suffix = "*metrics.csv"

        if hasattr(self,"analysis_dir"): #nanopore
            search_dir = self.analysis_dir
            variant_snpeff_suffix = "*.pass.SnpEff.vcf"
        else:
            search_dir = self.consensus_dir
            variant_snpeff_suffix = "*sorted.filtered.primerTrim.annotate.vcf" 

        for sample in [x for x in os.listdir(search_dir) if os.path.isdir(os.path.join(search_dir,x))]:

            sample_consensus = glob.glob(self.analysis_dir + "/" + sample + "/" + re.sub(r'_\d$','',sample)  + consensus_suffix) if hasattr(self,"analysis_dir") else glob.glob(self.consensus_dir + "/" + sample + "/" + re.sub(r'_\d$','',sample) + consensus_suffix)

            try:
                sample_consensus = sample_consensus[0]
            except:
                sample_consensus = ""

            cleaned_raw_reads = glob.glob(self.cleaned_raw_reads_dir  + "/" + sample + "/" + re.sub(r'_\d$','',sample) + clean_raw_reads_suffix) if not hasattr(self,"analysis_dir") else []

            host_removal = glob.glob(self.host_removal_dir  + "/" + sample + "/" + re.sub(r'_\d$','',sample) + host_removal_suffix) if not hasattr(self,"analysis_dir") else []

            metrics = glob.glob(self.analysis_dir  + "/" + sample + "/" + re.sub(r'_\d$','',sample) + metrics_suffix) if  hasattr(self,"analysis_dir") else glob.glob(self.metrics_dir + "/" + metrics_suffix)

            try:
                metrics = metrics[0]
            except:
                metrics = ""

            variant_snpeff = glob.glob(self.analysis_dir  + "/" + sample + "/" + re.sub(r'_\d$','',sample) + variant_snpeff_suffix) if  hasattr(self,"analysis_dir") else glob.glob(self.variant_dir + "/" + sample  +  "/" + re.sub(r'_\d$','',sample) + variant_snpeff_suffix)

            try:
                variant_snpeff = variant_snpeff[0]
            except:
                variant_snpeff = ""

            qc_status, perc_n = self.GetSampleQcStatus(metrics,sample)

            #print("QC STATUS ",qc_status)
            self.samples_obj_list.append(Sample(re.sub(r'_\d$','',sample),sample_consensus,cleaned_raw_reads,host_removal,metrics,variant_snpeff,qc_status,perc_n))

    def GetSampleQcStatus(self,metrics,sample_name):
        
        if len(metrics) == 0:
            return("missing","")
        #print(sample_name," ",self.run_name) 
        metrics_df = pd.read_csv(metrics,sep=",",index_col=False)

        if self.techno in ['illumina','MGI']:
            metrics_df = metrics_df[['sample','cons.per.N','bam.perc.50x']]
            metrics_df['cons.per.N'] = pd.to_numeric(metrics_df['cons.per.N'],errors='coerce')
            metrics_df['bam.perc.50x'] =  pd.to_numeric(metrics_df['bam.perc.50x'],errors='coerce')
            metrics_df = metrics_df.loc[metrics_df['sample'] == sample_name,['cons.per.N','bam.perc.50x']]
            if metrics_df.shape[0] == 0:
                return("missing","")
            cons_perc_n = list(metrics_df['cons.per.N'])[0]
            bam_perc_n = list(metrics_df['bam.perc.50x'])[0]
            #print("CONS PERC N ", cons_perc_n, " BAM ",bam_perc_n)
            #print("TYPE ",type(cons_perc_n), " TYPE ",type(bam_perc_n))
            #print("METRICS DF ",metrics_df)
        else: # nanopore
            metrics_df = metrics_df[['sample','cons.perc.N','bam.perc.50x']]
            metrics_df['cons.perc.N'] = pd.to_numeric(metrics_df['cons.perc.N'],errors='coerce')
            metrics_df['bam.perc.50x'] = pd.to_numeric(metrics_df['bam.perc.50x'],errors='coerce')
            metrics_df = metrics_df.loc[metrics_df['sample'] == sample_name,['cons.perc.N','bam.perc.50x']]
            if metrics_df.shape[0] == 0:
                return("missing","")
            cons_perc_n = list(metrics_df['cons.perc.N'])[0]
            bam_perc_n = list(metrics_df['bam.perc.50x'])[0]
            #print("CONS PERC N ", cons_perc_n, " BAM ",bam_perc_n)
            #print("TYPE ",type(cons_perc_n), " TYPE ",type(bam_perc_n))
            #print("METRICS DF ",metrics_df)

        if (str(cons_perc_n) != "nan") and (str(bam_perc_n) != "nan"):
            if cons_perc_n > 5:
                return("REJ",cons_perc_n)
            elif cons_perc_n > 1:
                return("FLAG",cons_perc_n)
            elif bam_perc_n < 90:
                return("FLAG",cons_perc_n)
            else:
                return("PASS",cons_perc_n)
        else:
            return("missing",cons_perc_n)

    def SetPath(self):
        if self.techno in ['illumina','MGI']:
            self.consensus_dir = os.path.join(self.run_path,"consensus") 
            self.alignment_dir = os.path.join(self.run_path,"alignment") 
            self.cleaned_raw_reads_dir = os.path.join(self.run_path,"cleaned_raw_reads") 
            self.host_removal_dir = os.path.join(self.run_path,"host_removal") 
            self.metrics_dir = os.path.join(self.run_path,"metrics") 
            self.variant_dir = os.path.join(self.run_path,"variant") 
            self.data_dir = os.path.join(self.run_path,"data") 
        else: #nanopore
            self.analysis_dir = glob.glob(os.path.join(self.run_path,"analysis" + "/*nanopolish_800x"))[0]

    def FindThisSample(self,sample_name):
        sample_obj_found = None

        for sample_obj in self.samples_obj_list:
            if sample_obj.GetSampleName() == sample_name:
                sample_obj_found = sample_obj
                sample_obj_found.SetIsInPlate(True)
                return(sample_obj_found)
        return(sample_obj_found)

class Plate():
    def __init__(self,plate_name,df):
        self.plate_name = plate_name
        self.df = df
        self.samples_list = list(self.df['Name'])

    def GetPlateName(self):
        return(self.plate_name)

    def GetSamplesList(self):
        return(self.samples_list)

class Sample():
    def __init__(self,sample_name,consensus,cleaned_raw_reads,host_removal,metrics,variant_snpeff,qc_status,cons_perc_n):
        self.sample_name = sample_name
        self.consensus_path = consensus
        self.cleaned_raw_reads = cleaned_raw_reads
        self.host_removal = host_removal
        self.metrics = metrics
        self.variant_snpeff = variant_snpeff
        self.qc_status = qc_status
        self.cons_perc_n = cons_perc_n
        self.is_in_plate = False

    def GetConsPercN(self):
        return(self.cons_perc_n)

    def GetQcStatus(self):
        return(self.qc_status)

    def SetIsInPlate(self,is_in_plate):
        self.is_in_plate = is_in_plate

    def GetIsInPlate(self):
        return(self.is_in_plate)

    def GetSampleName(self):
        return(self.sample_name)

    def GetConsensusPath(self):
        return(self.consensus_path)

    def GetCleanedRawReads(self):
        return(self.cleaned_raw_reads)

    def GetHostRemoval(self):
        return(self.host_removal)

    def GetMetrics(self):
        return(self.metrics)

    def GetVariantSnpeff(self):
        return(self.variant_snpeff)

class ListPlateSampleManager():
    def __init__(self,listPlateSample_name):
        self.listPlateSample_name = listPlateSample_name
        self.SetPdDf()

    def SetPdDf(self):
        self.pd_df = pd.read_csv(self.listPlateSample_name,sep="\t",index_col=False)

    def GetPdDf(self):
        return(self.pd_df)

class FileOutputManager():
    def __init__(self):
        self.missing_samples_filename = "MissingSamples.tsv"
        self.missing_qc_status_filename = "MissingQcStatus.tsv"
        self.missing_files_filename = "MissingFiles.tsv"
        self.consensus_list_filename = os.path.basename(__file__)[:-3] + "_" + datetime.now().strftime('%Y-%m-%d') + "_consensusList" +".list"
        self.vcf_list_filename = os.path.basename(__file__)[:-3] + "_" + datetime.now().strftime('%Y-%m-%d') + "_vcfList" +".list"

        self.missing_samples_handler = None
        self.missing_qc_status_handler = None
        self.missing_files_handler = None
        self.consensus_list_handler = None
        self.vcf_list_handler = None

        self.SetFilesHandler()

    def WriteMissingSample(self,sample_name,plate_name):
        self.missing_samples_handler.write(sample_name + "\t" + plate_name + "\n")

    def WriteMissingQcStatus(self,sample_name,plate_name,run_name):
        self.missing_qc_status_handler.write(sample_name + "\t" + plate_name + "\t" + run_name + "\n")

    def WriteMissingFile(self,sample_name,plate_name,run_name,filetype):
        self.missing_files_handler.write(sample_name + "\t" + plate_name + "\t" + run_name + "\t" + filetype + "\n")

    def WriteConsensusList(self,sample_name,qc_status,consensus_path,techno,perc_n,run_name,plate_name):
        self.consensus_list_handler.write(sample_name + "\t" + qc_status + "\t" + consensus_path + "\t" + techno + "\t" + str(perc_n) + "\t" + run_name + "\t" + plate_name + "\n")

    def WriteVcfList(self,sample_name,qc_status,vcf_path,techno,perc_n,run_name,plate_name):
        self.vcf_list_handler.write(sample_name + "\t" + qc_status + "\t" + vcf_path + "\t" + techno + "\t" + str(perc_n) + "\t" + run_name + "\t" + plate_name + "\n")

    def SetFilesHandler(self):
        self.missing_samples_handler = open(os.path.join(trace_path,self.missing_samples_filename),'w')
        self.missing_samples_handler.write("SAMPLE\tPLATE\n")

        self.missing_qc_status_handler = open(os.path.join(trace_path,self.missing_qc_status_filename),'w')
        self.missing_qc_status_handler.write("SAMPLE\tPLATE\tRUN\n")

        self.missing_files_handler = open(os.path.join(trace_path,self.missing_files_filename),'w')
        self.missing_files_handler.write("SAMPLE\tPLATE\tRUN\tFILETYPE\n")

        self.consensus_list_handler = open(os.path.join(trace_path,self.consensus_list_filename),'w')
        self.consensus_list_handler.write("SAMPLE\tSTATUS\tPATH\tTECHNO\tPERC_N\tRUN_NAME\tPATE_NAME\n")

        self.vcf_list_handler = open(os.path.join(trace_path,self.vcf_list_filename),'w')
        self.vcf_list_handler.write("SAMPLE\tSTATUS\tPATH\tTECHNO\tPERC_N\tRUN_NAME\tPATE_NAME\n")


    def CloseFilesHandler(self):
        self.missing_samples_handler.close()
        self.missing_files_handler.close()
        self.missing_qc_status_handler.close()
        self.consensus_list_handler.close()

class Logger():
    def __init__(self,logger_name,output):
        self.log_level = logging.INFO
        self.logger_name = logger_name
        self.logger = None
        self.formatter = None
        self.file_handler = None
        self.stdout_handler = None
        self.output = output
        self.Configure()

    def LogMessage(self,message):
        self.logger.info(message)

    def SetLogger(self):
        self.logger = logging.getLogger(self.logger_name)
        self.logger.setLevel(self.log_level)

    def SetStdOutHandler(self):
        self.stdout_handler = logging.StreamHandler()
        self.stdout_handler.setLevel(self.log_level)
        self.stdout_handler.setFormatter(self.formatter)

    def SetFormatter(self):
        self.formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', '%m/%d/%Y %I:%M:%S')

    def AddHandlerToLogger(self):
        self.logger.addHandler(self.file_handler)
        self.logger.addHandler(self.stdout_handler)

    def SetFileHandler(self):
        self.file_handler = logging.FileHandler(self.output)
        self.file_handler.setLevel(self.log_level)
        self.file_handler.setFormatter(self.formatter)

    def Configure(self):
        self.SetLogger()
        self.SetFormatter()
        self.SetStdOutHandler()
        self.SetFileHandler()
        self.AddHandlerToLogger()

def BuildNotFoundPlateObj(list_samples_without_plate):
    plate_name = 'NOTFOUND'
    plate_df = pd.DataFrame({'PlateOrSet':[plate_name]*len(list_samples_without_plate),'Name':list_samples_without_plate})
    #print("PLATE_DF \n",plate_df)
    return(Plate(plate_name,plate_df))

def BuildPlateObjList(listplate_sample_manager_obj):

    plate_obj_list = []

    df = listplate_sample_manager_obj.GetPdDf()
    uniq_plate_name_df = pd.DataFrame(df['PlateOrSet'].unique(),columns=['PlateOrSet'])

    for index, row in uniq_plate_name_df.loc[:,].iterrows():
        plate_name = row['PlateOrSet']
        plate_df = row.to_frame()
        plate_df = plate_df.T
        plate_df = pd.merge(plate_df,df,how='inner')[['PlateOrSet','Name']]
        plate_obj_list.append(Plate(plate_name,plate_df))

    return(plate_obj_list)

def FindThisSampleInRuns(sample,run_obj_list):
    
    found_runs_dict = {}

    for run in run_obj_list:
        sample_obj_found = run.FindThisSample(sample)
        if sample_obj_found:
            found_runs_dict[run] = sample_obj_found 

    return(found_runs_dict)

def BuildRepository(plate_obj_list,run_obj_list,logger,output_manager):
    for plate in plate_obj_list:
        logger.LogMessage("In plate " + plate.GetPlateName())

        plate_name = plate.GetPlateName()
        samples_list = plate.GetSamplesList()

        plate_path = os.path.join(repository_path,plate_name)

        try:
            os.mkdir(plate_path)
        except OSError as error:
            #logging.warning("Impossible de crée " + plate_name)
            pass

        for sample in samples_list:
            #TODO enregistrer les sample trouve null part
            found_runs_dict = FindThisSampleInRuns(sample,run_obj_list)

            if len(found_runs_dict) == 0:
                output_manager.WriteMissingSample(sample,plate_name)

            for found_run,samples_obj in found_runs_dict.items():
                sample_dir_name = ""
                run_techno = found_run.GetRunTechno()
                run_name = found_run.GetRunName()
                run_date = found_run.GetRunDate()
                run_path = found_run.GetRunPath()

                sample_dir_name = plate_name + "." + sample + "." + run_techno + "." + run_date 
                sample_dir_path = os.path.join(plate_path,sample_dir_name)

                try:
                    os.mkdir(sample_dir_path)
                except OSError as error:
                    #logging.warning("Impossible de crée " + sample_dir_path)
                    pass

                ################# Symlink consensus ################
                try:
                    src_consensus = samples_obj.GetConsensusPath()

                    if len(src_consensus) == 0:
                        output_manager.WriteMissingFile(sample,plate_name,run_name,"consensus")
                    else:
                        output_manager.WriteConsensusList(sample,samples_obj.GetQcStatus(),src_consensus,run_techno,samples_obj.GetConsPercN(),run_name,plate_name)
                        symlink_consensus = os.path.join(sample_dir_path,os.path.basename(src_consensus))
                        os.symlink(src_consensus,symlink_consensus)
                except OSError as error:
                    #logging.warning("Impossible de créer le symlink consensus pour " + sample_dir_path)
                    pass

                ################# Symlink cleaned raw reads ################
                try:
                    src_cleaned_raw_reads = samples_obj.GetCleanedRawReads()

                    if (len(src_cleaned_raw_reads) != 2) and (run_techno != "nanopore"):
                        output_manager.WriteMissingFile(sample,plate_name,run_name,"cleaned_raw_reads")
                    else:
                        for src in src_cleaned_raw_reads:
                            symlink = os.path.join(sample_dir_path,os.path.basename(src))
                            os.symlink(src,symlink)
                except OSError as error:
                    #logging.warning("Impossible de créer le symlink cleaned raw reds pour " + sample_dir_path)
                    pass


                ################# Symlink host removal ################
                try:
                    src_host_removal = samples_obj.GetHostRemoval()

                    if (len(src_host_removal) != 2) and (run_techno != "nanopore"):
                        output_manager.WriteMissingFile(sample,plate_name,run_name,"host_removal")
                    else:
                        for src in src_host_removal:
                            symlink = os.path.join(sample_dir_path,os.path.basename(src))
                            os.symlink(src,symlink)
                except OSError as error:
                    #logging.warning("Impossible de créer le symlink host removal pour " + sample_dir_path)
                    pass
                 
                ################# Symlink metrics ################
                try:
                    src_metrics = samples_obj.GetMetrics()

                    if len(src_metrics) == 0:
                        output_manager.WriteMissingFile(sample,plate_name,run_name,"metrics")
                        output_manager.WriteMissingQcStatus(sample,plate_name,run_name)
                    else:
                        if samples_obj.GetQcStatus() == "missing":
                            output_manager.WriteMissingQcStatus(sample,plate_name,run_name)

                        symlink_metrics = os.path.join(sample_dir_path,re.sub(r'_\d','',os.path.basename(src_metrics)))
                        os.symlink(src_metrics,symlink_metrics)
                except OSError as error:
                    #logging.warning("Impossible de créer le symlink metrics pour " + sample_dir_path)
                    pass

                ################# Symlink variant snpeff ################
                try:
                    src_variant_snpeff = samples_obj.GetVariantSnpeff()

                    if len(src_variant_snpeff) == 0:
                        output_manager.WriteMissingFile(sample,plate_name,run_name,"variant_snpeff")
                    else:
                        output_manager.WriteVcfList(sample,samples_obj.GetQcStatus(),src_variant_snpeff,run_techno,samples_obj.GetConsPercN(),run_name,plate_name)
                        symlink_variant_snpeff = os.path.join(sample_dir_path,re.sub(r'_\d','',os.path.basename(src_variant_snpeff)))
                        os.symlink(src_variant_snpeff,symlink_variant_snpeff)
                except OSError as error:
                    #logging.warning("Impossible de créer le symlink variant snpeff pour " + sample_dir_path)
                    pass

def MakeListOfSamplesWithoutPlate(run_obj_list):
    list_samples_without_plate = []

    for run in run_obj_list:
        l = run.GetSamplesWithoutPlate()
        list_samples_without_plate.extend(l)

    return(list_samples_without_plate)

def Main():
    today = datetime.today().strftime('%Y%m%d')

    process_logger = Logger("ProcessLogger",os.path.join(trace_path,"Process_{0}.log".format(today)))
    file_output_manager = FileOutputManager()

    listplate_sample_manager_obj = ListPlateSampleManager(input_file_path)
   
    process_logger.LogMessage("Bluild plates")
    plate_obj_list = BuildPlateObjList(listplate_sample_manager_obj)

    run_obj_list = []

    process_logger.LogMessage("Bluild runs")

    for key, value in {illumina_base_dir:'illumina',mgi_base_dir:'MGI',nanopore_base_dir:'nanopore'}.items():
        techno = value

        for run_name in os.listdir(key):
            run_path = os.path.join(key,run_name)
            run_obj = Run(run_name,run_path,techno)
            run_obj_list.append(run_obj)

    process_logger.LogMessage("Bluild repository")

    BuildRepository(plate_obj_list,run_obj_list,process_logger,file_output_manager)

    list_samples_without_plate = list(set(MakeListOfSamplesWithoutPlate(run_obj_list)))
    not_found_plate_obj = BuildNotFoundPlateObj(list_samples_without_plate)

    BuildRepository([not_found_plate_obj],run_obj_list,process_logger,file_output_manager)

    file_output_manager.CloseFilesHandler()
    process_logger.LogMessage("Terminé")

if __name__ == '__main__':
    Main()

