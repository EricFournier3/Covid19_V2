# coding=utf-8

"""
Eric Fournier 2020-12-02

"""

import os
import re
import sys

vcf_file = "/data/PROJETS/COVID-19_Beluga/VCF/JUS-V5290618.sorted.filtered.primerTrim.annotateTEST.vcf"

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
                        aa_change = aa_change_annotation_split[3] + " => " + aa_change_annotation_split[10][2:]
                    else:
                        aa_change = nuc_mut_type
                    aa_change_set.add(aa_change)
            
print(mut)
