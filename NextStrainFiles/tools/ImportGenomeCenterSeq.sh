#!/bin/bash

<<COM
Eric Fournier 2020-05-11
Import Covid-19 assemblies from ComputeCanada:Beluga

COM


green_message="\e[32m"
white_message="\e[39m"

lspq_data_nextstrain_dir="/data/Applications/GitScript/Covid19_V2/NextStrainFiles/data/lspq/"

#TODO Change fasta name
lspq_seq_for_nextstrain=${lspq_data_nextstrain_dir}"sequences_temp.fasta"

base_dir_slbio="/data/Runs/SARS-CoV-2/GenomeCenterSeq/"

out_seq_slbio=${base_dir_slbio}"FinalRelease/"
log_slbio="${out_seq_slbio}BelugaImport.log"

beluga_server="fournie1@beluga.computecanada.ca"
#TODO Change path
beluga_seq_dir="/home/fournie1/test"

slbio_user=$(whoami)
beluga_pass_file="/home/${slbio_user}/compute_canada_pass.txt"
read beluga_pw < ${beluga_pass_file}

scp_cmd="sshpass -p ${beluga_pw}  scp -r ${beluga_server}:${beluga_seq_dir} ${out_seq_slbio}"

echo -e "Begin import from ${beluga_server}:${beluga_seq_dir}\t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}

#echo -e "${green_message}INFO: " "Begin import from ${beluga_server}:${beluga_seq_dir}"
#TODO remove comment
#eval ${scp_cmd}
#echo -e "${green_message}INFO: " "End import from ${beluga_server}:${beluga_seq_dir}"

echo -e "End import from ${beluga_server}:${beluga_seq_dir}\t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}

#echo -e "${white_message}"

echo -e "Begin concat fasta to ${lspq_seq_for_nextstrain} \t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}
#TODO remove comment
#cat "${out_seq_slbio}$(basename ${beluga_seq_dir})/"*".fasta"  > ${lspq_seq_for_nextstrain}
echo -e "End concat fasta to ${lspq_seq_for_nextstrain} \t$(date "+%Y-%m-%d @ %H:%M$S")\n" >> ${log_slbio}


