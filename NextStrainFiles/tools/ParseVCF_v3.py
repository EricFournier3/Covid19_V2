# coding=utf-8

"""
Eric Fournier 2020-12-02

"""

import os
import re
import sys
import pandas as pd
import argparse
import glob

#vcf_basedir = "/home/foueri01@inspq.qc.ca/temp/TEST_GET_FASTA_FOR_NEXTSTRAIN/OUT/VCF/"
vcf_basedir = "/home/foueri01@inspq.qc.ca/temp/TEST_GET_FASTA_FOR_NEXTSTRAIN/OUT/VCF/small/"
#parsed_outdir = os.path.join(vcf_basedir,"PARSED")

parsed_outdir = os.path.join("/home/foueri01@inspq.qc.ca/temp/TEST_GET_FASTA_FOR_NEXTSTRAIN/OUT/VCF","PARSED","small")
parsed_outfile = os.path.join(parsed_outdir,"Variant_20200201_20201218_small.xlsx")
req_not_found_in_envois_gq = os.path.join(parsed_outdir,"ReqNotFound_20200201_20201218_small.xlsx")

vcf_infiles = glob.glob(vcf_basedir + "/*.vcf")

#print(vcf_infiles)

global mut_set
global spec_mut_dict

mut_set = set([])
spec_mut_dict = {}

def CompileMut():
    for spec,muts in spec_mut_dict.items():
       #print(muts)
       for mut in muts:
           #print(mut)
           mut_set.add(mut)
            
       


for vcf_file in vcf_infiles:
    spec = os.path.basename(vcf_file).split('.')[0].split('_')[0]
    #print(vcf_file, " => ",spec)
    spec_mut_dict[spec] = []
    header_read = False
    is_header = {}
    is_nanopore = False
    line_nb = 0
    #print(vcf_file)

    with open(vcf_file) as vcff:
        
        for line in vcff:
            line_nb += 1 
            is_header = False

            if (line == 2) and (not re.search(r'^##source=ivar',line)):
                is_nanopore = True
                #print(re.search(r'^##source=ivar',line))
                #print(is_nanopore)

            if re.search(r'^#CHROM\t',line):
                header_read = True
                is_header = True
            if header_read and not is_header:
                #nuc_change = ""
                aa_change_set = set([])

                split_line = line.split('\t')
                #nuc_change_pos = split_line[1]
                #nuc_ref = split_line[3]
                #nuc_mut = split_line[4]
                #nuc_change = nuc_ref + nuc_change_pos + nuc_mut
                nuc_change_qual = split_line[6]

                if nuc_change_qual != "PASS":
                    continue

                aa_change_annotations = split_line[7]
                
                aa_change_annotations = re.search(r'\S+ANN=(\S+)',aa_change_annotations).group(1)
                aa_change_first_annotation = aa_change_annotations.split(',')[0]
                aa_change_first_annotation = aa_change_first_annotation.split('|')
                if aa_change_first_annotation[1] == "missense_variant":
                    #print(">>>> ",aa_change_first_annotation)
                    gene_name = aa_change_first_annotation[3]
                    aa_mut = aa_change_first_annotation[10][2:]
                    spec_mut_dict[spec].append(gene_name + ":" + aa_mut) 

        #print(spec_mut_dict[spec])     
        #print("********************************************")

#print(spec_mut_dict) 
CompileMut()
#print(mut_seg)
#exit(1)
mut_df_columns = ['# Requête',"AA_MUTATIONS"]
mut_df_columns.extend(list(mut_set)) 
#print(mut_df_columns)
#mut_df = pd.DataFrame(columns=['# Requête',"AA_MUTATIONS"].extend(list(mut_set)))
mut_df = pd.DataFrame(columns=mut_df_columns)
#print(mut_df)


#TODO 20201221 REACTIVE
envois_gq_df = pd.read_excel("/data/Databases/CovBanQ_Epi/LISTE_ENVOIS_GENOME_QUEBEC/ListeEnvoisGenomeQuebec_2020-12-02.xlsx",sheet_name=0) 
#envois_gq_df = pd.read_excel("/data/Databases/CovBanQ_Epi/LISTE_ENVOIS_GENOME_QUEBEC/ListeEnvoisGenomeQuebec_test_variant.xlsx",sheet_name=0) 

#TODO 20201221 REACTIVE
envois_gq_df['# Requête'] = envois_gq_df['# Requête'].str.upper().str.strip(' ')

#TODO 20201221 REACTIVE
sgil_df = pd.read_table("/data/Databases/CovBanQ_Epi/SGIL_EXTRACT/extract_with_Covid19_extraction_v2_20201203_CovidPos.txt")


for spec,mut in spec_mut_dict.items():
    #print(spec, " : ",str(mut))
    my_mut_dict = {}
    for uniq_mut in list(mut_set):
        if uniq_mut in mut:
            my_mut_dict[uniq_mut] = "1"
        else:
            my_mut_dict[uniq_mut] = "0"
    #print(my_mut_dict) 
    my_mut_dict.update({'# Requête':spec,'AA_MUTATIONS':str(mut)})
    #mut_df = mut_df.append({'# Requête':spec,'AA_MUTATIONS':str(mut)},ignore_index=True)
    mut_df = mut_df.append(my_mut_dict,ignore_index=True)
    mut_df['# Requête'] = mut_df['# Requête'].str.upper()

#print(mut_df)


mut_df['# Requête'] = mut_df['# Requête'].str.replace(r'(^HGA-\S+)2D$',r'\1',regex=True)
mut_df['# Requête'] = mut_df['# Requête'].str.strip(' ')

temp_df = mut_df.loc[mut_df['# Requête'].str.contains('HGA'),'# Requête']
#print(temp_df)

all_cols = ['# Requête','Nom','Prénom','Date de naissance','NAM','Date de prélèvement','AA_MUTATIONS']
all_cols.extend(list(mut_set))
final_mut_df_1 = pd.merge(mut_df,envois_gq_df,how='inner',on = ['# Requête'])
final_mut_df_1 = final_mut_df_1[all_cols]

print("SHAPE 1",final_mut_df_1.shape[0])

final_mut_df_2 = pd.merge(mut_df,sgil_df,how='inner',left_on = '# Requête',right_on='NUMERO_SGIL')
final_mut_df_2 = final_mut_df_2.rename(columns={'DATE_NAISS':'Date de naissance','NOM':'Nom','PRENOM':'Prénom','SAMPLED_DATE':'Date de prélèvement'})
final_mut_df_2 = final_mut_df_2[all_cols]

print("SHAPE 2",final_mut_df_2.shape[0])
final_mut_df = pd.concat([final_mut_df_1,final_mut_df_2])

req_found = list(final_mut_df['# Requête'])
#print(final_mut_df)

req_not_found_in_envois_gq_df = mut_df.loc[~mut_df['# Requête'].isin(req_found)]
req_not_found_in_envois_gq_df = req_not_found_in_envois_gq_df['# Requête']

final_mut_df.to_excel(parsed_outfile,sheet_name='Sheet1',index=False)
req_not_found_in_envois_gq_df.to_excel(req_not_found_in_envois_gq,sheet_name='Sheet1')
