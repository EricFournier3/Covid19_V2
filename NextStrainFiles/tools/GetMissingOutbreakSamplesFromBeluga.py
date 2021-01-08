# -*- coding: utf-8 -*-

"""
Eric Fournier 2020-12-07
"""

import datetime
import pandas as pd
import os
import numpy as np
import re
import sys
import logging
import argparse
import glob

outbreak_dir = "/data/Databases/CovBanQ_Epi/OUTBREAK_NEW/"
metadata = "/data/PROJETS/COVID-19_Beluga/Metadata/metadata_2021-01-07_PASS_FLAG_minmaxSampleDate_2020-02-01_2021-01-06.tsv"
beluga_consensus = "/data/PROJETS/COVID-19_Beluga/UpdateListOfRunsConsensusAndVcfList/updateListOfRuns_v3_2021-01-07_consensusList.list" 
missing_run_samples = "/data/PROJETS/COVID-19_Beluga/UpdateListOfRunsConsensusAndVcfList/MissingSamples_20210107.tsv"
missing_qc_status = "/data/PROJETS/COVID-19_Beluga/UpdateListOfRunsConsensusAndVcfList/MissingQcStatus_20210107.tsv"

metadata_df = pd.read_csv(metadata,sep="\t",index_col=False)
beluga_consensus_df = pd.read_csv(beluga_consensus,sep="\t",index_col=False)
missing_run_samples_df = pd.read_csv(missing_run_samples,sep="\t",index_col=False)
missing_qc_status_df = pd.read_csv(missing_qc_status,sep="\t",index_col=False)

temp_dir = "/home/foueri01@inspq.qc.ca/temp/20210108/"

stat_file = os.path.join(temp_dir,"StatMissingOutbreakSamples.stat")
stat_file_handler = open(stat_file,'w')

outbreak_total_file = os.path.join(temp_dir,"OutbreakTotal.tsv")
missing_from_metadata_file = os.path.join(temp_dir,"MissingFromMetadata.tsv")
missing_from_consensus_file = os.path.join(temp_dir,"MissingFromMConsensus.tsv")
missing_from_metadata_in_consensus_file = os.path.join(temp_dir,"MissingFromMetadataInConsensus.tsv")

missing_from_missing_run_samples = os.path.join(temp_dir,"MissingFromRunSamples.tsv")
missing_from_consensus_in_missing_run_samples = os.path.join(temp_dir,"MissingFromConsensusInMissingRunSamples.tsv")

missing_from_missing_qc_status = os.path.join(temp_dir,"MissingFromMissingQcStatus.tsv")
missing_from_missing_run_samples_in_missing_qc_status = os.path.join(temp_dir,"MissingFromMissingRunSamplesInMissinQcStatus.tsv")



outbreak_file_list = []

for outbreak_file in glob.glob(outbreak_dir + "*.list"):
    outbreak_file_list.append(outbreak_file)

outbreak_df =  pd.DataFrame(columns=['COVBANK','BIOBANK','Outbreak'])

for outbreak_file in outbreak_file_list:
    file_name = os.path.basename(outbreak_file)
    outbreak_name = re.sub(r'\.list','',file_name)
    temp_df = pd.read_csv(outbreak_file,sep="\t",index_col=False)
    temp_df['Outbreak'] = outbreak_name
    outbreak_df = pd.concat([outbreak_df,temp_df])

outbreak_df.to_csv(outbreak_total_file,sep="\t",index=False)


################### Outbreak total ###################
stat_file_handler.write("Number of outbreak samples : " + str(outbreak_df.shape[0]) + "\n")
biobank_set_total = set(outbreak_df['BIOBANK'])


################### Search metadata ###################
merge_df_1 = pd.merge(outbreak_df,metadata_df,how='inner',left_on='BIOBANK',right_on='sample')
biobank_set_1  = set(merge_df_1['BIOBANK'])
biobank_missing_set_1 = biobank_set_total - biobank_set_1
missing_df_1 = outbreak_df.loc[outbreak_df['BIOBANK'].isin(biobank_missing_set_1),:] #enregister ce df
missing_df_1.to_csv(missing_from_metadata_file,sep="\t",index=False)
stat_file_handler.write("Number of missing from metadata : " + str(missing_df_1.shape[0]) + "\n")


################### Search consensus beluga ###################
merge_df_2 = pd.merge(missing_df_1,beluga_consensus_df,how='inner',left_on='BIOBANK',right_on='SAMPLE')
merge_df_2 = merge_df_2.sort_values(by=['BIOBANK']) #enregistrer ce df
biobank_set_2 = set(merge_df_2['BIOBANK'])
biobank_missing_set_2 = biobank_missing_set_1 - biobank_set_2
missing_df_2 = outbreak_df.loc[outbreak_df['BIOBANK'].isin(biobank_missing_set_2),:]#enregistrer ce df
missing_df_2.to_csv(missing_from_consensus_file,sep="\t",index=False)
merge_df_2.to_csv(missing_from_metadata_in_consensus_file,sep="\t",index=False)
stat_file_handler.write("Number of missing from consensus : " + str(missing_df_2.shape[0]) + "\n")
stat_file_handler.write("Number of missing from metadata in consensus : " + str(merge_df_2.shape[0]) + "\n")


################### Search in MissingSamples from beluga runs ###################
merge_df_3 = pd.merge(missing_df_2,missing_run_samples_df,how='inner',left_on='BIOBANK',right_on='SAMPLE')
merge_df_3 = merge_df_3.sort_values(by=['BIOBANK'])
biobank_set_3 = set(merge_df_3['BIOBANK'])
biobank_missing_set_3 = biobank_missing_set_2 - biobank_set_3
missing_df_3 = outbreak_df.loc[outbreak_df['BIOBANK'].isin(biobank_missing_set_3),:]#enregistrer ce df
missing_df_3.to_csv(missing_from_missing_run_samples,sep="\t",index=False)
merge_df_3.to_csv(missing_from_consensus_in_missing_run_samples,sep="\t",index=False)
stat_file_handler.write("Number of missing from missing run samples : " + str(missing_df_3.shape[0]) + "\n")
stat_file_handler.write("Number of missing from consensus in missing run samples : " + str(merge_df_3.shape[0]) + "\n")

################### Search in MissingQcStatus from beluga runs ###################
merge_df_4 = pd.merge(missing_df_2,missing_qc_status_df,how='inner',left_on='BIOBANK',right_on='SAMPLE')
merge_df_4 = merge_df_4.sort_values(by=['BIOBANK'])
biobank_set_4 = set(merge_df_4['BIOBANK'])
biobank_missing_set_4 = biobank_missing_set_3 - biobank_set_4
missing_df_4 = outbreak_df.loc[outbreak_df['BIOBANK'].isin(biobank_missing_set_4),:]#enregistrer ce df
missing_df_4.to_csv(missing_from_missing_qc_status,sep="\t",index=False)
merge_df_4.to_csv(missing_from_missing_run_samples_in_missing_qc_status,sep="\t",index=False)
stat_file_handler.write("Number of missing from missing qc status : " + str(missing_df_4.shape[0]) + "\n")
stat_file_handler.write("Number of missing from run samples in missing qc status : " + str(merge_df_4.shape[0]) + "\n")

stat_file_handler.close()
