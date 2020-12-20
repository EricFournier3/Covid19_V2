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

vcf_basedir = "/home/foueri01@inspq.qc.ca/temp/TEST_GET_FASTA_FOR_NEXTSTRAIN/OUT/VCF/"
parsed_outdir = os.path.join(vcf_basedir,"PARSED")
parsed_outfile = os.path.join(parsed_outdir,"Variant_20200201_20201218.xlsx")

vcf_infiles = glob.glob(vcf_basedir + "/*.vcf")

#print(vcf_infiles)

spec_mut_dict = {}

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
     
mut_df = pd.DataFrame(columns=['# Requête',"AA_MUTATIONS"])
envois_gq_df = pd.read_excel("/data/Databases/CovBanQ_Epi/LISTE_ENVOIS_GENOME_QUEBEC/ListeEnvoisGenomeQuebec_test_variant.xlsx",sheet_name=0) 

for spec,mut in spec_mut_dict.items():
    #print(spec, " : ",str(mut))
    mut_df = mut_df.append({'# Requête':spec,'AA_MUTATIONS':str(mut)},ignore_index=True)
    mut_df['# Requête'] = mut_df['# Requête'].str.upper()

final_mut_df = pd.merge(mut_df,envois_gq_df,how='inner',on = ['# Requête'])
final_mut_df = final_mut_df[['# Requête','Nom','Prénom','Date de naissance','NAM','AA_MUTATIONS']]
#print(final_mut_df)
final_mut_df.to_excel(parsed_outfile,sheet_name='Sheet1')

