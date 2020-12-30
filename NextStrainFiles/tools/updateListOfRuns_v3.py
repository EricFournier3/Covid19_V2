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
        self.SetSamplesList()

        #print(self.samples_list)

    def SetSamplesList(self):
        self.samples_list = []

        if hasattr(self,"analysis_dir"): #nanopore
            search_dir = self.analysis_dir
        else:
            search_dir = self.consensus_dir

        for sample in [x for x in os.listdir(search_dir) if os.path.isdir(os.path.join(search_dir,x))]:
            self.samples_list.append(sample)
            #TODO creer des object Sample ici

    def SetPath(self):
        if self.techno in ['illumina','mgi']:
            self.consensus_dir = os.path.join(self.run_path,"consensus") 
            self.alignment_dir = os.path.join(self.run_path,"alignment") 
            self.cleaned_raw_reads_dir = os.path.join(self.run_path,"cleaned_raw_reads") 
            self.host_removal_dir = os.path.join(self.run_path,"host_removal") 
            self.metrics_dir = os.path.join(self.run_path,"metrics") 
            self.variant_dir = os.path.join(self.run_path,"variant") 
            self.data_dir = os.path.join(self.run_path,"data") 
        else: #nanopore
            self.analysis_dir = glob.glob(os.path.join(self.run_path,"analysis" + "/*nanopolish_800x"))[0]


class Plate():
    def __init__(self,plate_name,df):
        self.plate_name = plate_name
        self.df = df

class Sample():
    def __init__(self,sample_name):
        self.sample_name = sample_name

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

    print(plate_obj_list)

def Main():
    listplate_sample_manager_obj = ListPlateSampleManager(input_file_path)
    
    plate_obj_list = BuildPlateObjList(listplate_sample_manager_obj)

    run_obj_list = []

    for key, value in {illumina_base_dir:'illumina',mgi_base_dir:'mgi',nanopore_base_dir:'nanopore'}.items():
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


if __name__ == '__main__':
    Main()
    logging.info('Terminé')

