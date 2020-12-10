# coding=utf-8

"""
Eric Fournier 2020-12-08
Note: need to run this script in nextstrain conda env

"""

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
import subprocess
import shutil
from Bio.SeqIO.FastaIO import SimpleFastaParser
from scipy import sparse


logging.basicConfig(level=logging.INFO)

beluga_server = "fournie1@beluga.computecanada.ca:/home/fournie1"
mnt_beluga_server = "/mnt/BelugaEric/"

def MountBelugaServer():
    logging.info("Try to mount Beluga")
    os.system("sudo umount " + mnt_beluga_server)
    os.system("sudo sshfs -o allow_other -o follow_symlinks {0} {1}".format(beluga_server,mnt_beluga_server))
    logging.info("Beluga mounted")

#MountBelugaServer()


compare_outdir = "/data/PROJETS/COVID-19_Beluga/Consensus/CompareFreeze1_OLDvsNEW_pipeline/"
new_consensus_outdir = os.path.join(compare_outdir,"NewConsensus") #TODO CHANGE TO NewConsensus
temp_dir = os.path.join(compare_outdir,"temp")
align_dir = os.path.join(compare_outdir,"Aligned")
nuc_diff_df_out = os.path.join(align_dir,"Freeze1_ConsensusDiff.xlsx")
nuc_diff_df_out_no_N = os.path.join(align_dir,"Freeze1_ConsensusDiff_no_N.xlsx")

illumina_reprocess_path_file = os.path.join(compare_outdir,"illumina_reprocess_path.txt")
mgi_reprocess_path_file = os.path.join(compare_outdir,"mgi_reprocess_path.txt")

illumina_reprocess_path_df = pd.read_table(illumina_reprocess_path_file)
illumina_reprocess_path_df.columns = ['PATH']

mgi_reprocess_path_df = pd.read_table(mgi_reprocess_path_file)
mgi_reprocess_path_df.columns = ['PATH']

base_dir_gisaid_sub = "/data/PROJETS/COVID-19_Beluga/Gisaid/FinalPublished/release1/"

base_dir_gisaid_sub_1 = os.path.join(base_dir_gisaid_sub,"20200520")
base_dir_gisaid_sub_2 = os.path.join(base_dir_gisaid_sub,"20200610")
base_dir_gisaid_sub_3 = os.path.join(base_dir_gisaid_sub,"20200914")

excel_gisaid_metadata_sub_1 = os.path.join(base_dir_gisaid_sub_1,"20200520_ncov19_metadata.xls")
excel_gisaid_metadata_sub_2 = os.path.join(base_dir_gisaid_sub_2,"20200610_ncov19_metadata.xls")
excel_gisaid_metadata_sub_3 = os.path.join(base_dir_gisaid_sub_3,"20200914_ncov19_metadata.xls")

excel_gisaid_metadata_sub_1_df = pd.read_excel(excel_gisaid_metadata_sub_1,sheet_name=0)
excel_gisaid_metadata_sub_2_df = pd.read_excel(excel_gisaid_metadata_sub_2,sheet_name=0)
excel_gisaid_metadata_sub_3_df = pd.read_excel(excel_gisaid_metadata_sub_3,sheet_name=0)

excel_gisaid_metadata_all_df_all = pd.concat([excel_gisaid_metadata_sub_1_df,excel_gisaid_metadata_sub_2_df,excel_gisaid_metadata_sub_3_df])
excel_gisaid_metadata_all_df_all =  excel_gisaid_metadata_all_df_all.drop(excel_gisaid_metadata_all_df_all.index[0])
excel_gisaid_metadata_all_df_all = excel_gisaid_metadata_all_df_all[['covv_virus_name','covv_seq_technology']]
excel_gisaid_metadata_all_df_all =    excel_gisaid_metadata_all_df_all.reset_index(drop=True)

fasta_gisaid_sequence_sub_1 = os.path.join(base_dir_gisaid_sub_1,"all_sequences.fasta")
fasta_gisaid_sequence_sub_2 = os.path.join(base_dir_gisaid_sub_2,"all_sequences.fasta")
fasta_gisaid_sequence_sub_3 = os.path.join(base_dir_gisaid_sub_3,"all_sequences.fasta")


gisaid_rec_dict = {}
gisaid_rec_dict_2 = {}


def GetTechno(seq_id):
    techno = excel_gisaid_metadata_all_df_all.loc[excel_gisaid_metadata_all_df_all['covv_virus_name'] == seq_id,'covv_seq_technology']
    return(list(techno)[0])

for gisaid_fasta in [fasta_gisaid_sequence_sub_1,fasta_gisaid_sequence_sub_2,fasta_gisaid_sequence_sub_3]:
    for rec in SeqIO.parse(gisaid_fasta,'fasta'):
        techno = GetTechno(rec.id)
        gisaid_rec_dict[rec.id] = techno
        gisaid_rec_dict_2[re.sub(r'hCoV-19/','',rec.id)] = rec

print("len(gisaid_rec_dict) : ",len(gisaid_rec_dict)) # len(gisaid_rec_dict) :  734

new_rej_consensus = []
new_multiple_consensus = []
new_keeped_consensus_list = []

gisaid_consensus_not_found_in_new = []
nb_gisaid_nanapore_consensus = 0

check = 0

check2 = 0
check3 = 0
check4 = 0
check5 = 0
check5a = 0
check6 = 0

for seq_id in  gisaid_rec_dict.keys():
    check += 1
    keeped_path = ""
    path = ""
    path_list = [] 
    beluga_seq_id = ""
    techno = ""   

    beluga_seq_id = re.search(r'^hCoV-19/Canada/Qc-(\S+)/\d{4}',seq_id).group(1)
    techno = gisaid_rec_dict[seq_id]
    if techno == "Illumina_NexteraFlex":
        path = illumina_reprocess_path_df.loc[illumina_reprocess_path_df['PATH'].str.contains(beluga_seq_id),'PATH']
    elif techno == "MGI_CleanPlex":
        path = mgi_reprocess_path_df.loc[mgi_reprocess_path_df['PATH'].str.contains(beluga_seq_id),'PATH']
    else:
        nb_gisaid_nanapore_consensus += 1
        continue

    path_list = list(path)

    if len(path_list) == 0:
        check2 += 1
        gisaid_consensus_not_found_in_new += 1
        continue
    elif len(path_list) > 1:
        check3 += 1
        new_multiple_consensus.append(path_list)
        #print("PATH LIST ",path_list)
        pass_list = []
        flag_list = []
        rej_list = []

        for i in path_list:
            if re.search(r'pass',os.path.basename(i)):
                pass_list.append(i)                
            elif re.search(r'flag',os.path.basename(i)):
                flag_list.append(i)
            elif re.search(r'rej',os.path.basename(i)):
                rej_list.append(i)

        if len(pass_list) > 0:
            keeped_path = pass_list[0]
        elif len(flag_list) > 0:
            keeped_path = flag_list[0]       
        elif len(rej_list) > 0:
            #print("REJ LIST ",rej_list)
            new_rej_consensus.append(rej_list)
            continue 
    else:
        check3 += 1
        keeped_path = path_list[0]

    if re.search(r'rej',os.path.basename(keeped_path)):
        check5 += 1
        #print("REJ FOR ",keeped_path)
        new_rej_consensus.append(keeped_path)
        continue

    check5a += 1
    keeped_path = re.sub(r'/home/fournie1/',r'/mnt/BelugaEric/',keeped_path)
    if len(keeped_path) > 0:
        check6 += 1
        new_keeped_consensus_list.append(keeped_path)  
    else:
        print("STRANGE PATH for ", beluga_seq_id, " ", techno ) 


duplicated_path = []

def FindDuplicatedPath():
    temp = []
    for p in new_keeped_consensus_list:
        if p in temp:
            duplicated_path.append(p)
        temp.append(p)
   
FindDuplicatedPath()
 
new_keeped_consensus_set = set(new_keeped_consensus_list)
#print("len(new_keeped_consensus_set)", len(new_keeped_consensus_set)) # len(new_keeped_consensus_set) 517
#print("duplicated_path ", duplicated_path) # duplicated_path  []
#print("len(duplicated_path)",len(duplicated_path)) # len(duplicated_path) 0

#print("Check ",check, " check2 ",check2," check3 ",check3, " check4 ",check4," check5 ",check5, " check5a ",check5a, " check6 ",check6)

for fasta in new_keeped_consensus_list:
    pass
    #logging.info("Get " + fasta)
    #shutil.copy2(fasta,new_consensus_outdir) 
           
#print("len(new_multiple_consensus) ", len(new_multiple_consensus)) # len(new_multiple_consensus)  67
#print("len(new_keeped_consensus_list) : ", len(new_keeped_consensus_list)) #len(new_keeped_consensus_list) :  517
#print("nb_gisaid_nanapore_consensus  : ",nb_gisaid_nanapore_consensus) # nb_gisaid_nanapore_consensus  :  210
#print("gisaid_consensus_not_found_in_new : ", gisaid_consensus_not_found_in_new) # gisaid_consensus_not_found_in_new :  []
#print("new_rej_consensus : ",new_rej_consensus) # [['/home/fournie1/COVID_full_processing/mgi_reprocess/20200605_mgi_LSPQPlate03_V300063030/consensus/L00227094/L00227094.consensus.MGI.rej.fasta', '/home/fournie1/COVID_full_processing/mgi_reprocess/20200608/consensus/L00227094/L00227094.consensus.MGI.rej.fasta'], ['/home/fournie1/COVID_full_processing/mgi_reprocess/20200605_mgi_LSPQPlate03_V300063030/consensus/L00227481/L00227481.consensus.MGI.rej.fasta', '/home/fournie1/COVID_full_processing/mgi_reprocess/20200608/consensus/L00227481/L00227481.consensus.MGI.rej.fasta'], ['/home/fournie1/COVID_full_processing/mgi_reprocess/20200605_mgi_LSPQPlate05_V300063030/consensus/L00232614/L00232614.consensus.MGI.rej.fasta', '/home/fournie1/COVID_full_processing/mgi_reprocess/20200608/consensus/L00232614/L00232614.consensus.MGI.rej.fasta'], ['/home/fournie1/COVID_full_processing/mgi_reprocess/20200605_mgi_LSPQPlate05_V300063030/consensus/L00232813/L00232813.consensus.MGI.rej.fasta', '/home/fournie1/COVID_full_processing/mgi_reprocess/20200608/consensus/L00232813/L00232813.consensus.MGI.rej.fasta'], '/home/fournie1/COVID_full_processing/mgi_reprocess/20200424_mgi_LSPQPlate09_V300035341/consensus/L00241242/L00241242.consensus.MGI.rej.fasta', '/home/fournie1/COVID_full_processing/illumina_reprocess/20200609_illumina_LSPQPlate06_HM2CTDRXX/consensus/L00235085/L00235085.consensus.illumina.rej.fasta', '/home/fournie1/COVID_full_processing/illumina_reprocess/20200616_illumina_LSPQPlate11_HM275DRXX/consensus/L00243973/L00243973.consensus.illumina.rej.fasta']
#print("len(new_rej_consensus) ",len(new_rej_consensus)) # 7


def RunCmd(cmd):
    env = os.environ.copy()
    shellexec = ['env','bash']
    shargs = ['-c', "set -euo pipefail; " + cmd]

    try:
        subprocess.check_output(shellexec + shargs,shell = False, stderr = subprocess.STDOUT,env = env)
    except subprocess.CalledProcessError as error:
        print("BUG")
        print(cmd)
        return False
    else:
        return True

#new_fasta_consensus = glob.glob(new_consensus_outdir + "/*.fasta")[0:3]
new_fasta_consensus = glob.glob(new_consensus_outdir + "/*.fasta")

nuc_diff_df_list = []
nuc_diff_df_no_N_list = []


def CheckSnp(fasta_align,spec_id):
    logging.info("Work on " + spec_id) # voir /data/Applications/GitScript/Covid19_V2/NextStrainFiles/scripts/priorities.py
    
    fh = open(fasta_align, 'rt')
    
    generator = SimpleFastaParser(fh)

    old = next(generator)
    old_s = old[1]
    old_np_s = np.array(list(old_s))

    new = next(generator)
    new_s = new[1]
    new_np_s = np.array(list(new_s))

    old_np  = np.frombuffer(old_s.lower().encode('utf-8'), dtype=np.int8).copy()
    new_np  = np.frombuffer(new_s.lower().encode('utf-8'), dtype=np.int8).copy()

    snps_np_bool = new_np!=old_np

    new_snps_np = new_np[snps_np_bool]

    old_snps_np = old_np[snps_np_bool]

    new_snps_np_s = new_np_s[snps_np_bool]

    old_snps_np_s = old_np_s[snps_np_bool]

    snps_np_pos = np.nonzero(snps_np_bool)[0]

    id_np = np.full((len(snps_np_pos),1),spec_id)

    final_np = np.array(list(zip(snps_np_pos,new_snps_np_s,old_snps_np_s)))
    final_np = np.append(final_np,id_np,axis=1)

    df = pd.DataFrame(final_np,columns=['POS','NEW_NUC','OLD_NUC','ID'])
    df_no_N = df.loc[(df['NEW_NUC'] != 'N') & (df['OLD_NUC'] != 'N'),:]
    nuc_diff_df_list.append(df)
    nuc_diff_df_no_N_list.append(df_no_N)

    fh.close()



for fasta in new_fasta_consensus:
    #print(fasta)
    rec_new = SeqIO.read(fasta,'fasta')
    #print(rec.id) 
    spec_id = re.search(r'Canada/Qc-(\S+)/\d{4}',rec_new.id).group(1)
    #print(gisaid_rec_dict_2[rec.id])
    rec_list = []
    rec_list.extend([rec_new,gisaid_rec_dict_2[rec_new.id]])
    
    not_align = os.path.join(temp_dir,'not_align.fasta')
    align = os.path.join(align_dir,spec_id + '_align.fasta')
    SeqIO.write(rec_list,not_align,'fasta')

    align_cmd="mafft --reorder --anysymbol --nomemsave --adjustdirection --thread 40 {0} > {1}".format(not_align,align)
    success = RunCmd(align_cmd)
    
    CheckSnp(align,spec_id)

nuc_diff_df = pd.concat(nuc_diff_df_list)
nuc_diff_df.to_excel(nuc_diff_df_out,sheet_name='Sheet1',index=False)

nuc_diff_df_no_N = pd.concat(nuc_diff_df_no_N_list)
nuc_diff_df_no_N.to_excel(nuc_diff_df_out_no_N,sheet_name='Sheet1',index=False)
exit(0)

# *******************************************   TEST   **************************************
'''
global nuc_diff_df
nuc_diff_df = pd.DataFrame(columns=['POS','NEW','OLD','ID'])


def CheckSnpTest(fastafile,spec_id,final_df):
    logging.info("Work on " + spec_id)
    #voir /data/Applications/GitScript/Covid19_V2/NextStrainFiles/scripts/priorities.py
    fh = open(fastafile, 'rt')
   
    print(nuc_diff_df)
 
    gen = SimpleFastaParser(fh)
    
    old = next(gen)
    old_s = old[1]
    print("OLD ",old_s) # aagaatttatgaacgt
    old_np_s = np.array(list(old_s))        

    new = next(gen)
    new_s = new[1]
    print("NEW ",new_s) # aaaaatttgtgaacgt
    new_np_s = np.array(list(new_s))

    print("old_np_s ",old_np_s) # ['a' 'a' 'g' 'a' 'a' 't' 't' 't' 'a' 't' 'g' 'a' 'a' 'c' 'g' 't']
    print("new_np_s ",new_np_s) # ['a' 'a' 'a' 'a' 'a' 't' 't' 't' 'g' 't' 'g' 'a' 'a' 'c' 'g' 't']

    old_np  = np.frombuffer(old_s.lower().encode('utf-8'), dtype=np.int8).copy()
    #old_np[(old_np!=97) & (old_np!=99) & (old_np!=103) & (old_np!=116)] = 97 PAS BESOIN DE CA
    print(old_np) # [ 97  97 103  97  97 116 116 116  97 116 103  97  97  99 103 116]

    new_np  = np.frombuffer(new_s.lower().encode('utf-8'), dtype=np.int8).copy()
    #new_np[(new_np!=97) & (new_np!=99) & (new_np!=103) & (new_np!=116)] = 97  PAS BESOIN DE CA
    print(new_np) # [ 97  97  97  97  97 116 116 116 103 116 103  97  97  99 103 116]


    snps_np_bool = new_np!=old_np
    print("snps_np_bool", snps_np_bool) # [False False  True False False False False False  True False False False False False False False]
    
    new_snps_np = new_np[snps_np_bool]
    print("new_snps_np ",new_snps_np) # [ 97 103]

    old_snps_np = old_np[snps_np_bool]
    print("old_snps_np ",old_snps_np) #  [103  97] 

    print(" ******* ")    

    new_snps_np_s = new_np_s[snps_np_bool]
    print("new_snp_np_s ",new_snps_np_s) # ['a' 'g']

    old_snps_np_s = old_np_s[snps_np_bool]
    print("old_snp_np_s ",old_snps_np_s) # ['g' 'a']

    snps_np_pos = np.nonzero(snps_np_bool)[0]
    print("snp_np_pos ",snps_np_pos) # [2 8]

    id_np = np.full((len(snps_np_pos),1),spec_id)
    print("id_np ",id_np)

    final_np = np.array(list(zip(snps_np_pos,new_snps_np_s,old_snps_np_s)))
    final_np = np.append(final_np,id_np,axis=1)
    print("Final np \n",final_np) # [['2' 'a' 'g' 'L00241115Old'] ['8' 'g' 'a' 'L00241115Old']] 
    

    df = pd.DataFrame(final_np,columns=['POS','NEW','OLD','ID'])
    
    final_df = pd.concat([final_df,df])
    
    #print("final_df \n ",final_df)
    # POS NEW OLD ID
    # 2    a  g  L00241115
    # 8    g  a  L00241115

    #mettre dans un pd.dataframe et ensuite dans excel voir tableau dans cahier
    #TODO sauvegarde des alignements
    fh.close()


rec_list = []
rec1 = SeqIO.read(os.path.join(temp_dir,'fasta1Old.fasta'),'fasta')
print("ID ",rec1.id)
spec_id = re.search(r'Canada/Qc-(\S+)/\d{4}',rec1.id).group(1)
rec2 = SeqIO.read(os.path.join(temp_dir,'fasta1New.fasta'),'fasta')
rec_list.extend([rec1,rec2])
not_align = os.path.join(temp_dir,'not_align.fasta')
align = os.path.join(temp_dir,spec_id + '_aligned.fasta')
SeqIO.write(rec_list,not_align,'fasta')
align_cmd="mafft --reorder --anysymbol --nomemsave --adjustdirection --thread 40 {0} > {1}".format(not_align,align)
success = RunCmd(align_cmd)

CheckSnpTest(align,spec_id,nuc_diff_df)

rec_list = []
rec1 = SeqIO.read(os.path.join(temp_dir,'fasta2Old.fasta'),'fasta')
print("ID ",rec1.id)
spec_id = re.search(r'Canada/Qc-(\S+)/\d{4}',rec1.id).group(1)
rec2 = SeqIO.read(os.path.join(temp_dir,'fasta2New.fasta'),'fasta')
rec_list.extend([rec1,rec2])
not_align = os.path.join(temp_dir,'not_align.fasta')
align = os.path.join(temp_dir,spec_id + '_aligned.fasta')
SeqIO.write(rec_list,not_align,'fasta')
align_cmd="mafft --reorder --anysymbol --nomemsave --adjustdirection --thread 40 {0} > {1}".format(not_align,align)
success = RunCmd(align_cmd)

CheckSnpTest(align,spec_id,nuc_diff_df)
print("nuc_diff_df ",nuc_diff_df)
exit(0)
'''
