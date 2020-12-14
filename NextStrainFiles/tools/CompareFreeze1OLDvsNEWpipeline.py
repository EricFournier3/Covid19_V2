# coding=utf-8

"""
Eric Fournier 2020-12-08
Note: need to run this script in nextstrain conda env

"""
import matplotlib.pyplot as plt
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
new_consensus_outdir = os.path.join(compare_outdir,"NewConsensus") 
temp_dir = os.path.join(compare_outdir,"temp")
align_dir = os.path.join(compare_outdir,"Aligned")
nuc_diff_df_out = os.path.join(align_dir,"Freeze1_ConsensusDiff.xlsx")
nuc_diff_df_out_no_N = os.path.join(align_dir,"Freeze1_ConsensusDiff_no_N.xlsx")
no_diff_out = os.path.join(align_dir,"NoDiffSpecList.txt")
no_diff_out_handle = open(no_diff_out,'a')

illumina_reprocess_path_file = os.path.join(compare_outdir,"illumina_reprocess_path.txt")
mgi_reprocess_path_file = os.path.join(compare_outdir,"mgi_reprocess_path.txt")

illumina_reprocess_path_df = pd.read_table(illumina_reprocess_path_file)
illumina_reprocess_path_df.columns = ['PATH']

mgi_reprocess_path_df = pd.read_table(mgi_reprocess_path_file)
mgi_reprocess_path_df.columns = ['PATH']

base_dir_gisaid_sub = "/data/PROJETS/COVID-19_Beluga/Gisaid/FinalPublished/release1/"

old_beluga_fasta_path_file = "/data/PROJETS/COVID-19_Beluga/Gisaid/FinalPublished/release1/Freeze1FastaPath.txt"
old_beluga_fasta_path_df = pd.read_table(old_beluga_fasta_path_file,header=None)
old_beluga_fasta_path_df.columns = ['PATH']


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
rec_dict = {}

def GetGisaidTechnoFromFastaHeader(seq_id):
    techno = excel_gisaid_metadata_all_df_all.loc[excel_gisaid_metadata_all_df_all['covv_virus_name'] == seq_id,'covv_seq_technology']
    return(list(techno)[0])

def GetTechoAndQcStatusFromBelugaFastaPath(seq_id):
    seq_id = re.search(r'Canada/Qc-(\S+)/\d{4}',seq_id).group(1)
    path = old_beluga_fasta_path_df.loc[old_beluga_fasta_path_df['PATH'].str.contains(seq_id),'PATH']
    path = list(path)[0]
    fasta_name = os.path.basename(path)
    techno = ""
    qc_status = ""
    try:
        search_obj = re.search(r'L\S+\.consensus\.(\S+)\.(\S+)\.fasta',fasta_name)
        techno = search_obj.group(1)
        qc_status = search_obj.group(2)
    except:
        print("No parse for fasta_name ",fasta_name)
    return((techno,qc_status))


def CompareOldTechnoFromFastaHeader2FastaName():
    for my_id in rec_dict.keys():
        techno_from_header = rec_dict[my_id]['old_techno_from_fasta_header']
        if techno_from_header == 'ONT_ARTIC':
            techno_from_header = 'nanopore'

        techno_from_name = rec_dict[my_id]['old_techno_from_beluga_fasta_path']
        if not re.search(techno_from_name,techno_from_header,re.IGNORECASE):
            print("Mismatch old techno for ", my_id, " header: ",techno_from_header, " name : ",techno_from_name)

for gisaid_fasta in [fasta_gisaid_sequence_sub_1,fasta_gisaid_sequence_sub_2,fasta_gisaid_sequence_sub_3]:
    for rec in SeqIO.parse(gisaid_fasta,'fasta'):
        techno = GetGisaidTechnoFromFastaHeader(rec.id)
        gisaid_rec_dict[rec.id] = techno
        id_minimal = re.sub(r'hCoV-19/','',rec.id) 
        old_techno_from_beluga_fasta_path, old_qc_status_from_beluga_fasta_path = GetTechoAndQcStatusFromBelugaFastaPath(id_minimal)
        gisaid_rec_dict_2[id_minimal] = rec
        rec_dict[id_minimal] = {'old_techno_from_beluga_fasta_path':old_techno_from_beluga_fasta_path,'new_techno_from_beluga_fasta_path':"",'old_techno_from_fasta_header':techno,'old_qcstatus_from_beluga_fasta_path':old_qc_status_from_beluga_fasta_path,'new_qcstatus_from_beluga_fasta_path':""}
        
CompareOldTechnoFromFastaHeader2FastaName()


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


#Get new consensus from Beluga
for seq_id in  gisaid_rec_dict.keys():
    check += 1
    keeped_path = ""
    path = ""
    path_list = [] 
    beluga_seq_id = ""
    techno = ""   

    beluga_seq_id = re.search(r'^hCoV-19/Canada/Qc-(\S+)/\d{4}',seq_id).group(1)
    minimal_id = re.sub(r'hCoV-19/','',seq_id)
    techno = gisaid_rec_dict[seq_id]
    if techno == "Illumina_NexteraFlex":
        path = illumina_reprocess_path_df.loc[illumina_reprocess_path_df['PATH'].str.contains(beluga_seq_id),'PATH']
    elif techno in ["MGI_CleanPlex","MGI"]:
        path = mgi_reprocess_path_df.loc[mgi_reprocess_path_df['PATH'].str.contains(beluga_seq_id),'PATH']
    else:
        #on ne tient pas compte de nanopore
        nb_gisaid_nanapore_consensus += 1
        continue

    #print(beluga_seq_id)

    path_list = list(path)

    if len(path_list) == 0:
        check2 += 1
        gisaid_consensus_not_found_in_new += 1
        continue
    elif len(path_list) > 1:
        check3 += 1
        new_multiple_consensus.append(path_list)
        pass_list = []
        flag_list = []
        rej_list = []

        for i in path_list:
            if re.search(r'L\S{8}\.consensus\.\S+\.\S+.fasta',i):
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
            keeped_path = rej_list[0]       
        else:
            print("No new fasta found for ",seq_id)
            continue
    else:
        check3 += 1
        keeped_path = path_list[0]


    #print("KEEPED ",keeped_path)


    check5a += 1
    keeped_path = re.sub(r'/home/fournie1/',r'/mnt/BelugaEric/',keeped_path)
    if len(keeped_path) > 0:
        check6 += 1
        search_obj = re.search(r'L\S{8}\.consensus\.(\S+)\.(\S+).fasta',os.path.basename(keeped_path))
        new_techno = search_obj.group(1)
        new_qc_status = search_obj.group(2)
        rec_dict[minimal_id]['new_techno_from_beluga_fasta_path'] = new_techno
        rec_dict[minimal_id]['new_qcstatus_from_beluga_fasta_path'] = new_qc_status
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
print("len(new_keeped_consensus_set)", len(new_keeped_consensus_set)) #  len(new_keeped_consensus_set) 541
print("duplicated_path ", duplicated_path) #  duplicated_path  []
print("len(duplicated_path)",len(duplicated_path)) # len(duplicated_path) 0

#print("Check ",check, " check2 ",check2," check3 ",check3, " check4 ",check4," check5 ",check5, " check5a ",check5a, " check6 ",check6)

for fasta in new_keeped_consensus_list:
    #logging.info("Get " + fasta)
    #shutil.copy2(fasta,new_consensus_outdir) 
    pass
           
def RunCmd(cmd):
    env = os.environ.copy()
    shellexec = ['env','bash']
    shargs = ['-c', "set -euo pipefail; " + cmd]

    try:
        subprocess.check_output(shellexec + shargs,shell = False, stderr = subprocess.STDOUT,env = env)
    except subprocess.CalledProcessError as error:
        logging.error("Bug in RunCmd")
        print(cmd)
        return False
    else:
        return True

#new_fasta_consensus = glob.glob(new_consensus_outdir + "/*.fasta")[0:3]
new_fasta_consensus = glob.glob(new_consensus_outdir + "/*.fasta")

nuc_diff_df_list = []
nuc_diff_df_no_N_list = []

#Check differences between old and new consensus
def CheckSnp(fasta_align,spec_id):
    #logging.info("Work on " + spec_id) 
    
    fh = open(fasta_align, 'rt')
    
    generator = SimpleFastaParser(fh)

    new = next(generator)
    my_seq_id = re.search(r'(^\S+) ',new[0]).group(1)
    new_s = new[1]
    new_np_s = np.array(list(new_s))

    old = next(generator)
    old_s = old[1]
    old_np_s = np.array(list(old_s))


    old_np  = np.frombuffer(old_s.lower().encode('utf-8'), dtype=np.int8).copy()
    new_np  = np.frombuffer(new_s.lower().encode('utf-8'), dtype=np.int8).copy()

    snps_np_bool = new_np!=old_np

    new_snps_np = new_np[snps_np_bool]

    old_snps_np = old_np[snps_np_bool]

    new_snps_np_s = new_np_s[snps_np_bool]

    old_snps_np_s = old_np_s[snps_np_bool]

    snps_np_pos = np.nonzero(snps_np_bool)[0]

    id_np = np.full((len(snps_np_pos),1),spec_id)

    old_qc_status = rec_dict[my_seq_id]['old_qcstatus_from_beluga_fasta_path']
    old_qc_status_np = np.full((len(snps_np_pos),1),old_qc_status) 
    #print(old_qc_status_np)

    old_techno_from_fasta_header = rec_dict[my_seq_id]['old_techno_from_fasta_header']
    old_techno_from_fasta_header_np = np.full((len(snps_np_pos),1),old_techno_from_fasta_header)
    #print(old_techno_from_fasta_header_np)

    old_techno_from_fasta_name = rec_dict[my_seq_id]['old_techno_from_beluga_fasta_path']
    old_techno_from_fasta_name_np = np.full((len(snps_np_pos),1),old_techno_from_fasta_name)
    #print(old_techno_from_fasta_name_np)    

    new_qc_status = rec_dict[my_seq_id]['new_qcstatus_from_beluga_fasta_path']
    new_qc_status_np = np.full((len(snps_np_pos),1),new_qc_status)

    new_techno = rec_dict[my_seq_id]['new_techno_from_beluga_fasta_path']
    new_techno_np = np.full((len(snps_np_pos),1),new_techno)

    final_np = np.array(list(zip(snps_np_pos,new_snps_np_s,old_snps_np_s)))
    #print(final_np)

    try:
        final_np = np.append(final_np,id_np,axis=1)
        final_np = np.append(final_np,old_qc_status_np,axis=1)
        final_np = np.append(final_np,new_qc_status_np,axis=1)
        final_np = np.append(final_np,old_techno_from_fasta_header_np,axis=1)
        final_np = np.append(final_np,old_techno_from_fasta_name_np,axis=1)
        final_np = np.append(final_np,new_techno_np,axis=1)

        df = pd.DataFrame(final_np,columns=['POS','NEW_NUC','OLD_NUC','ID','OLD_QC_STATUS','NEW_QC_STATUS','OLD_TECHNO_FROM_HEADER','OLD_TECHNO_FROM_NAME','NEW_TECHNO'])
        df_no_N = df.loc[(df['NEW_NUC'] != 'N') & (df['OLD_NUC'] != 'N'),:]
        nuc_diff_df_list.append(df)
        nuc_diff_df_no_N_list.append(df_no_N)
    except:
        print(">>>>>>>>>>>>>>>>>>>>> No difference for  ", spec_id)
        #print(id_np)
        #print("---------")
        #print(final_np)
        no_diff_out_handle.write(spec_id + "\n")
    fh.close()

for fasta in new_fasta_consensus:
    rec_new = SeqIO.read(fasta,'fasta')
    spec_id = re.search(r'Canada/Qc-(\S+)/\d{4}',rec_new.id).group(1)
    rec_list = []
    rec_list.extend([rec_new,gisaid_rec_dict_2[rec_new.id]])
    
    not_align = os.path.join(temp_dir,'not_align.fasta')
    align = os.path.join(align_dir,spec_id + '_align.fasta')
    SeqIO.write(rec_list,not_align,'fasta')

    align_cmd="mafft --reorder --anysymbol --nomemsave --adjustdirection --thread 40 {0} > {1}".format(not_align,align)
    #success = RunCmd(align_cmd)
    
    CheckSnp(align,spec_id)


nuc_diff_df = pd.concat(nuc_diff_df_list)

nuc_diff_df.to_excel(nuc_diff_df_out,sheet_name='Sheet1',index=False)
nuc_diff_df_no_N = pd.concat(nuc_diff_df_no_N_list)
nuc_diff_df_no_N.to_excel(nuc_diff_df_out_no_N,sheet_name='Sheet1',index=False)
no_diff_out_handle.close()


qc_status_df_excel = os.path.join(align_dir,"QcStatus.xlsx") 
old_qc_status_df = nuc_diff_df[['ID','OLD_QC_STATUS']].drop_duplicates()
old_qc_status_df = old_qc_status_df.rename(columns={'OLD_QC_STATUS':'QC_STATUS'})
old_qc_status_df = old_qc_status_df.reset_index(drop=True)
old_qc_status_df['VERSION'] = 'OLD'

new_qc_status_df = nuc_diff_df[['ID','NEW_QC_STATUS']].drop_duplicates()
new_qc_status_df = new_qc_status_df.rename(columns={'NEW_QC_STATUS':'QC_STATUS'})
new_qc_status_df = new_qc_status_df.reset_index(drop=True)
new_qc_status_df['VERSION'] = 'NEW'

qc_status_df = pd.concat([old_qc_status_df,new_qc_status_df])
qc_status_df = qc_status_df.reset_index(drop=True)
qc_status_df.to_excel(qc_status_df_excel,index=False,sheet_name='Sheet1')

out_qc_status_hist = os.path.join(align_dir,"QcStatus.png")
ax = qc_status_df['QC_STATUS'].hist(by=qc_status_df['VERSION'])
for i,my_ax in enumerate(ax):
    my_ax.set_xlabel("Qc status")
    rects = my_ax.patches
    for rect in rects:
        #print("RECT ",rect)
        #x_val = rect.get_width()
        x_val = rect.get_x() + rect.get_width() / 2
        #print("X VAL ",x_val)
        #y_val = rect.get_y() + rect.get_height() / 2
        y_val = rect.get_height()
        #print("Y VAL ",y_val)
        space = 2
        ha = 'left'
        if x_val < 0:
            space *= -1
            ha = 'right'

        #label = "{:.1f}".format(y_val)
        label = str(int(y_val))
        #print("label ",label)
        if y_val > 0:
            my_ax.annotate(label,(x_val,y_val),xytext=(-1,10),textcoords="offset points",va='center',ha='center',rotation=60) # voir https://matplotlib.org/3.3.3/api/_as_gen/matplotlib.axes.Axes.annotate.html
    
    if i == 0:
        my_ax.set_ylabel("Samples Frequency")
plt.suptitle("QC Status Distribution")
plt.savefig(out_qc_status_hist)
plt.show()

out_nb_sites_diff = os.path.join(align_dir,"NbSitesDiff.png")
group_by_id_df = nuc_diff_df.groupby(['ID']).agg('count')
ax = group_by_id_df['NEW_NUC'].hist(bins=20)
ax.set_xlabel("Number of different positions")
ax.set_ylabel("Samples Frequency")
plt.suptitle("New versus old consensus")
plt.savefig(out_nb_sites_diff)
plt.show()


out_nb_sites_diff_no_N = os.path.join(align_dir,"NbSitesDiffNoN.png")
group_by_id_df_no_N = nuc_diff_df_no_N.groupby(['ID']).agg('count')
ax = group_by_id_df_no_N['NEW_NUC'].hist(bins=20)
ax.set_xlabel("Number of different positions")
ax.set_ylabel("Samples Frequency")
plt.suptitle("New versus old consensus (excluding Ns)")
plt.savefig(out_nb_sites_diff_no_N)
plt.show()

def CheckIfTechnoFromHeaderSameAsTechnoFromFastaName(row):
    if row.OLD_TECHNO_FROM_NAME.upper() in row.OLD_TECHNO_FROM_HEADER.upper():
        return True
    else:
        return False

techno_qc_status_df_out = os.path.join(align_dir,"TechnoAndQcStatus.xlsx")

techno_qc_status_df = nuc_diff_df[['ID','OLD_QC_STATUS','NEW_QC_STATUS','OLD_TECHNO_FROM_HEADER','OLD_TECHNO_FROM_NAME','NEW_TECHNO']].drop_duplicates()
#print(techno_qc_status_df)

techno_qc_status_df['TECHO_HEADER_SAMEAS_TECHNO_FASTA_NAME'] = techno_qc_status_df.apply(CheckIfTechnoFromHeaderSameAsTechnoFromFastaName,axis=1)

#print(techno_qc_status_df)
techno_qc_status_df.to_excel(techno_qc_status_df_out,index=False,sheet_name='Sheet1')

exit(0)

