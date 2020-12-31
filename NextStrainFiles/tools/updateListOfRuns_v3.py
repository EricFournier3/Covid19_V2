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


'''
Created on 29 December 2020
@author: Eric Fournier

Note:
    - pour executer ce script, il faut au prealable executer <module load scipy-stack>


TODO:
    - combiner updateSampleList_v2.py et updateListOfRuns_v2.py
    - modifier le README
    - ne pas oublier le NOTFOUND dans REPOSITORY
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
        #print("self.run_path >> ",self.run_path)

        self.SetPath()
        self.SetSamplesObjList()

        #print(self.samples_obj_list)

    def GetRunPath(self):
        return(self.run_path)

    def GetRunName(self):
        return(self.run_name)

    def GetRunDate(self):
            return(self.run_name.split('_')[-1] if self.techno == 'nanopore' else self.run_name.split('_')[0])

    def GetRunTechno(self):
        return(self.techno.lower())

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

            #print("APPEND ",sample)
            self.samples_obj_list.append(Sample(re.sub(r'_\d$','',sample),sample_consensus,cleaned_raw_reads,host_removal,metrics,variant_snpeff))

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
        #print("TRY TO FIND ",sample_name)
        sample_obj_found = None
        for sample_obj in self.samples_obj_list:
            if sample_obj.GetSampleName() == sample_name:
                #print("YES FOUND ",sample_name)
                sample_obj_found = sample_obj
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
    def __init__(self,sample_name,consensus,cleaned_raw_reads,host_removal,metrics,variant_snpeff):
        self.sample_name = sample_name
        self.consensus_path = consensus
        self.cleaned_raw_reads = cleaned_raw_reads
        self.host_removal = host_removal
        self.metrics = metrics
        self.variant_snpeff = variant_snpeff

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
        pass

def BuildPlateObjList(listplate_sample_manager_obj):

    plate_obj_list = []

    df = listplate_sample_manager_obj.GetPdDf()
    uniq_plate_name_df = pd.DataFrame(df['PlateOrSet'].unique(),columns=['PlateOrSet'])
    #print("uniq_plate_name \n",uniq_plate_name_df)

    for index, row in uniq_plate_name_df.loc[:,].iterrows():
        plate_name = row['PlateOrSet']
        plate_df = row.to_frame()
        plate_df = plate_df.T
        plate_df = pd.merge(plate_df,df,how='inner')[['PlateOrSet','Name']]
        #print(plate_df)
        plate_obj_list.append(Plate(plate_name,plate_df))

    return(plate_obj_list)


def FindThisSampleInRuns(sample,run_obj_list):
    
    found_runs_dict = {}

    for run in run_obj_list:
        #print("RUN ",run.GetRunName())
        sample_obj_found = run.FindThisSample(sample)
        if sample_obj_found:
            found_runs_dict[run] = sample_obj_found 

    return(found_runs_dict)

def BuildRepository(plate_obj_list,run_obj_list):
    for plate in plate_obj_list:
        print("PLATE ",plate.GetPlateName())
        print("    SAMPLE LIST ",plate.GetSamplesList())
        plate_name = plate.GetPlateName()
        samples_list = plate.GetSamplesList()

        plate_path = os.path.join(repository_path,plate_name)
        try:
            os.mkdir(plate_path)
        except OSError as error:
            logging.warning("Impossible de crée " + plate_name)

        for sample in samples_list:
            #TODO enregistrer les sample trouve null part
            found_runs_dict = FindThisSampleInRuns(sample,run_obj_list)
            #print("For ", sample , " => FOUND RUNS ",[x.GetRunName() for x in found_runs])
            for found_run,samples_obj in found_runs_dict.items():
                sample_dir_name = ""
                run_techno = found_run.GetRunTechno()
                run_name = found_run.GetRunName()
                run_date = found_run.GetRunDate()
                run_path = found_run.GetRunPath()

                sample_dir_name = plate_name + "." + sample + "." + run_techno + "." + run_date 
                sample_dir_path = os.path.join(plate_path,sample_dir_name)
                #print("SAMPLE DIR PATH ",sample_dir_path)

                try:
                    os.mkdir(sample_dir_path)
                except OSError as error:
                    logging.warning("Impossible de crée " + sample_dir_path)

                ################# Symlink consensus ################
                try:
                    src_consensus = samples_obj.GetConsensusPath()
                    #print("SRC CONSENSUS ", src_consensus)
                    symlink_consensus = os.path.join(sample_dir_path,os.path.basename(src_consensus))
                    #print("SYMLINK CONSENSUS ",symlink_consensus)
                    os.symlink(src_consensus,symlink_consensus)
                except OSError as error:
                    logging.warning("Impossible de créer le symlink consensus pour " + sample_dir_path)


                ################# Symlink cleaned raw reads ################
                try:
                    src_cleaned_raw_reads = samples_obj.GetCleanedRawReads()
                    for src in src_cleaned_raw_reads:
                        print("SRC CLEANED RAW READS ",src)
                        symlink = os.path.join(sample_dir_path,os.path.basename(src))
                        print("SYMLINK CLEANED RAW READS ",symlink)
                        os.symlink(src,symlink)
                except OSError as error:
                    logging.warning("Impossible de créer le symlink cleaned raw reds pour " + sample_dir_path)
                 


        print("--------------------------------------------------------")

def Main():
    listplate_sample_manager_obj = ListPlateSampleManager(input_file_path)
    
    plate_obj_list = BuildPlateObjList(listplate_sample_manager_obj)

    run_obj_list = []

    for key, value in {illumina_base_dir:'illumina',mgi_base_dir:'MGI',nanopore_base_dir:'nanopore'}.items():
        techno = value
        #print("techno >> ",techno)
        for run_name in os.listdir(key):
            #print("run_name >> ",run_name)
            run_path = os.path.join(key,run_name)
            #print("run_path >> ",run_path)
            run_obj = Run(run_name,run_path,techno)
            #print("HAS ATTR >> ",hasattr(run_obj,"analysis_dir"))
            run_obj_list.append(run_obj)

    #print(run_obj_list)
    BuildRepository(plate_obj_list,run_obj_list)

if __name__ == '__main__':
    Main()
    logging.info('Terminé')

