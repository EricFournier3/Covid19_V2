# coding=utf-8

"""
Eric Fournier 2020-10-06

TODO LOGGIN TO FILE

"""

import glob
import datetime
import re
import os
import sys
import argparse
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(description="Prepare and send Covid-19 data to NML")
parser.add_argument('--task','-t',required=True,help=" <1> Prepare sample list for Beluga server <2> Transfer fastq from Beluga server and prepare data <3> Send data to LNM",choices=['1','2','3'])
parser.add_argument('-d','--debug',help='run in debug mode',action='store_true')
parser.add_argument('--release','-r',help='data release',required=True,type=int)

args = parser.parse_args()


global beluga_server
beluga_server = "fournie1@beluga.computecanada.ca:/home/fournie1"

global mnt_beluga_server
mnt_beluga_server = "/mnt/BelugaEric"

global __debug__
_debug_ = args.debug

release = "release" + str(args.release)

global raw_data_base_dir
raw_data_basedir = os.path.join("/home/foueri01@inspq.qc.ca/temp/TEST_TRANSFER_LNM/Rawdata",release)


task = args.task

def MountBelugaServer():
    os.system("sudo umount " + mnt_beluga_server)
    os.system("sudo sshfs -o allow_other -o follow_symlinks {0} {1}".format(beluga_server,mnt_beluga_server))

#MountBelugaServer()


class Logger:
    pass

class ComputeCanadaToLSPQManager:
    def __init__(self,release):
        self.release = release

        if __debug__:
            self.gisaid_sample_list_file = os.path.join("/home/foueri01@inspq.qc.ca/temp/TEST_TRANSFER_LNM/Gisaid/FinalPublished/",self.release,"SampleListPublished_" + self.release + ".txt")
            self.beluga_fastq_basedir = os.path.join(mnt_beluga_server,"test2",self.release)
            self.raw_data_base_dir = os.path.join("/home/foueri01@inspq.qc.ca/temp/TEST_TRANSFER_LNM/Rawdata",self.release)
        else:
            self.gisaid_sample_list_file = os.path.join("/home/foueri01@inspq.qc.ca/temp/TEST_TRANSFER_LNM/Gisaid/FinalPublished/",self.release,"SampleListPublished_" + self.release + ".txt")
            self.beluga_fastq_basedir = os.path.join(mnt_beluga_server,"test2",self.release)
            self.raw_data_base_dir = os.path.join("/home/foueri01@inspq.qc.ca/temp/TEST_TRANSFER_LNM/Rawdata",self.release)

        self.gisaid_mgi_id_list = []
        self.gisaid_illumina_id_list = []

        self.missing_beluga_mgi = []
        self.missing_beluga_illumina = []


    def GetMissingBelugaSamples(self):
        for sample in self.gisaid_mgi_id_list:
            if not re.search(sample,self.beluga_mgi_fastq_basename_list):
                self.missing_beluga_mgi.append(sample)

        for sample in self.gisaid_illumina_id_list:
            if not re.search(sample,self.beluga_illumina_fastq_basename_list):
                self.missing_beluga_illumina.append(sample)

    def SaveMissingBelugaSamples(self):
        missing_beluga_sample_file = os.path.join(self.raw_data_base_dir,"MissingBelugaSamples.txt")

        with open(missing_beluga_sample_file,'w') as out_file:
            for spec in self.missing_beluga_mgi:
                out_file.write(spec + "\t" + "MGI\n")
            for spec in self.missing_beluga_illumina:
                out_file.write(spec + "\t" + "Illumina\n")


    def BuildGisaidIdList(self):
        line_pattern = re.compile("^Canada/Qc-(?P<id>\S+)/\d{4} seq_method:(?P<techno>\S+)\|assemb_method")

        with open(self.gisaid_sample_list_file) as gisaid_id_file:
            for line in gisaid_id_file:
                match = line_pattern.search(line)
                spec_id = match.group('id')
                techno = match.group('techno')

                if techno == "MGI" or techno == "MGI_CleanPlex":
                    self.gisaid_mgi_id_list.append(spec_id)
                elif techno == "Illumina_NexteraFlex":
                    self.gisaid_illumina_id_list.append(spec_id)

    def GetComputeCanadaFastq(self):
        #todo log transfer file et verbose
        slbio_mgi_fastq_dir = os.path.join(self.raw_data_base_dir,"mgi")
        slbio_illumina_fastq_dir = os.path.join(self.raw_data_base_dir,"illumina")
        
        #todo renommer les fastq
   
        print("***************  MGI  ********************")
        for fastq in self.beluga_mgi_fastq_list:
            print(fastq)
            print(os.path.dirname(fastq), " ____ ",os.path.basename(fastq)) 
            parsed_list = str(os.path.basename(fastq)).split('_')
            if len(parsed_list) == 2:
                spec_name = re.sub('\.host','',parsed_list[0])
            elif len(parsed_list) > 2:
                spec_name = parsed_list[0]

            print("SAMPLE NAME ",spec_name)


        print("***************  ILLUMINA  ********************")
        for fastq in self.beluga_illumina_fastq_list:
            print(fastq) 


    def GetComputeCanadaFastqList(self):
        self.beluga_mgi_fastq_list = glob.glob(os.path.join(self.beluga_fastq_basedir,"mgi") + "/*.fastq.gz")
        self.beluga_mgi_fastq_basename_list = "-".join([os.path.basename(x) for x in self.beluga_mgi_fastq_list if re.search('pair1',x)])
        self.beluga_illumina_fastq_list = glob.glob(os.path.join(self.beluga_fastq_basedir,"illumina") + "/*.fastq.gz")
        self.beluga_illumina_fastq_basename_list = "-".join([os.path.basename(x) for x in self.beluga_illumina_fastq_list if re.search('pair1',x)])
        self.GetMissingBelugaSamples()
        self.SaveMissingBelugaSamples()


class GisaidFastaManager:
    def __init__(self,release):
        self.release = release

        if __debug__:
            self.beluga_sample_list_out = os.path.join(mnt_beluga_server,"test2",self.release,"SampleListPublished_" + self.release + ".txt")
            self.basedir_gisaid = os.path.join("/home/foueri01@inspq.qc.ca/temp/TEST_TRANSFER_LNM/Gisaid/FinalPublished/",release)
        else:
            self.beluga_sample_list_out = os.path.join(mnt_beluga_server,self.release + "_dehosted_raw_reads","SampleListPublished_" + self.release + ".txt")
            self.basedir_gisaid = os.path.join("/home/foueri01@inspq.qc.ca/temp/TEST_TRANSFER_LNM/Gisaid/FinalPublished/",release)

        self.slbio_sample_list_out = os.path.join(self.basedir_gisaid,"SampleListPublished_" + self.release + ".txt")
        
        self.rec_list = []

    def WriteList(self):
        with open(self.slbio_sample_list_out,'w') as sample_file_out:
            for rec in self.rec_list:
                sample_file_out.write(rec.description + "\n")
    
        os.system("cp " +  self.slbio_sample_list_out + " " +  self.beluga_sample_list_out)

            
    def BuildList(self):
        for fasta_dir in os.listdir(self.basedir_gisaid):
            for fasta in glob.glob(os.path.join(self.basedir_gisaid,fasta_dir) + "/*.fasta"):
                if os.path.basename(fasta) != "all_sequences.fasta":
                    for rec in SeqIO.parse(fasta,'fasta'):
                        self.rec_list.append(rec)

    def CopyFastaToRawData(self):
        if os.path.isdir(raw_data_basedir):
            pass
        else:
            os.makedirs(os.path.join(raw_data_basedir,"fasta_temp"))
            os.makedirs(os.path.join(raw_data_basedir,"illumina"))
            os.makedirs(os.path.join(raw_data_basedir,"mgi"))
            os.makedirs(os.path.join(raw_data_basedir,"nanopore"))

        for rec in self.rec_list:
            fasta_name = re.search(r'Canada/Qc-(\S+)/\d{4}',rec.id).group(1) + ".fasta"
            fasta_path = os.path.join(raw_data_basedir,"fasta_temp",fasta_name)
            SeqIO.write(rec,fasta_path,'fasta')


if str(task) == "1":
    logging.info("Prepare Sample list for Beluga server")
    gisaid_fasta_manager = GisaidFastaManager(release)    
    gisaid_fasta_manager.BuildList()
    gisaid_fasta_manager.WriteList()
    #gisaid_fasta_manager.CopyFastaToRawData()
elif str(task) == "2":
    logging.info("Get fastq from Beluga")
    cc2lspq_manager = ComputeCanadaToLSPQManager(release)    
    cc2lspq_manager.BuildGisaidIdList()
    cc2lspq_manager.GetComputeCanadaFastqList()
    cc2lspq_manager.GetComputeCanadaFastq()
