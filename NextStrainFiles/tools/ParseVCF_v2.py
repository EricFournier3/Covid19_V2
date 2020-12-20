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
vcf_outdir = os.path.join(vcf_basedir,"PARSED")


vcf_infiles = glob.glob(vcf_basedir + "/*.vcf")

#print(vcf_infiles)

for vcf_file in vcf_infiles:
    header_read = False
    is_header = {}
    is_nanopore = False
    line_nb = 0
    print(vcf_file)

    with open(vcf_file) as vcff:
        aa_mut_list =  []
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
                    aa_mut_list.append(gene_name + ":" + aa_mut) 

        print(aa_mut_list)     
        print("********************************************")
                


exit(0)

vcf_file = "/data/PROJETS/COVID-19_Beluga/VCF/JUS-V5290618.sorted.filtered.primerTrim.annotateTEST.vcf"
out_test = "/data/PROJETS/COVID-19_Beluga/VCF/variant_test.xlsx"


header_read = False
is_header = False

mut = {}

with open(vcf_file) as vcff:
    for line in vcff:

        is_header = False

        if re.search(r'^#CHROM\t',line):
            header_read = True
            is_header = True
        if header_read and not is_header:
            nuc_change = ""
            aa_change_set = set([])

            split_line = line.split('\t')
            nuc_change_pos = split_line[1]
            nuc_ref = split_line[3]
            nuc_mut = split_line[4]
            nuc_change = nuc_ref + nuc_change_pos + nuc_mut
            nuc_change_qual = split_line[6]

            if nuc_change_qual != "PASS":
                continue

            mut[nuc_change] = aa_change_set

            aa_change_annotations = split_line[7]
            aa_change_annotations = aa_change_annotations.split(',')

            for aa_change_annotation in aa_change_annotations:
                aa_change_annotation_split = aa_change_annotation.split('|')
                nuc_mut_type = aa_change_annotation_split[1]
                if (len(aa_change_annotation_split[15]) == 0):
                    if (nuc_mut_type == 'missense_variant') :
                        aa_change = aa_change_annotation_split[3] + ":" + aa_change_annotation_split[10][2:]
                    else:
                        aa_change = nuc_mut_type
                    aa_change_set.add(aa_change)
    
nuc_mut_string = ""
aa_mut_string = ""
        
for nuc_mut in mut.keys():
    nuc_mut_string = nuc_mut_string + nuc_mut + "|"
    for aa_mut in mut[nuc_mut]:
        aa_mut_string = aa_mut_string + aa_mut + ","
    
    aa_mut_string = aa_mut_string.strip(',')
    aa_mut_string = aa_mut_string + "|"

pd_df = pd.DataFrame({"NUC_VARIANT":[nuc_mut_string],"AA_VARIANT":[aa_mut_string]})
pd_df.to_excel(out_test,sheet_name='Sheet1',index=False)

