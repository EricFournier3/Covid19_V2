#!/bin/bash

<<COM
Eric Fournier 2020-05-11
Import Covid-19 assemblies from ComputeCanada:Beluga

COM

green_message="\e[32m"
white_message="\e[39m"

lspq_data_nextstrain_dir="/data/Applications/GitScript/Covid19_V2/NextStrainFiles/data/lspq/"

lspq_seq_for_nextstrain=${lspq_data_nextstrain_dir}"sequences.fasta"

base_dir_slbio="/data/Runs/SARS-CoV-2/GenomeCenterSeq/"

#TODO enlever TEST
out_seq_slbio=${base_dir_slbio}"FinalRelease/TEST/"
log_slbio="${out_seq_slbio}BelugaImport.log"

beluga_server="fournie1@beluga.computecanada.ca"
#TODO Change path
beluga_seq_dir="/home/fournie1/COVID_consensus/"

slbio_user=$(whoami)
beluga_pass_file="/home/${slbio_user}/compute_canada_pass.txt"
read beluga_pw < ${beluga_pass_file}

final_unpublished="FinalUnpublished"
final_published="FinalPublished"
transfered_2_lspq_dir=${final_unpublished}"/Transfered2Lspq/"

ImportSeqFromBeluga(){
  #echo -e "${green_message}INFO: " "Begin import from ${beluga_server}:${beluga_seq_dir}"
  echo -e "Begin import from ${beluga_server}:${beluga_seq_dir}\t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}

  #scp_cmd_obsolete="sshpass -p ${beluga_pw}  scp -r ${beluga_server}:\"${beluga_seq_dir}{${final_unpublished}/*.fasta,${final_published}/*.fasta}\" ${out_seq_slbio}"
  #scp_cmd_obsolete_2="sshpass -p ${beluga_pw}  scp -r ${beluga_server}:\"${beluga_seq_dir}{${final_published}/,${final_unpublished}/}\" ${out_seq_slbio}"

  
  scp_cmd="sshpass -p ${beluga_pw}  scp -r ${beluga_server}:\"${beluga_seq_dir}${final_unpublished}/*.fasta\" ${out_seq_slbio}${final_unpublished}"

  #TODO remove comment
  eval ${scp_cmd}

  #echo -e "${green_message}INFO: " "End import from ${beluga_server}:${beluga_seq_dir}"
  echo -e "End import from ${beluga_server}:${beluga_seq_dir}\t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}

  #echo -e "${white_message}"
}


ExportSeqForNextstrain(){
  
 echo "In ExportSeqForNextstrain"

  #echo -e "Begin concat fasta to ${lspq_seq_for_nextstrain} \t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}

  #TODO remove comment
  #cat ${out_seq_slbio}{${final_unpublished}/*.fasta,${final_published}/*.fasta} >${lspq_seq_for_nextstrain}

  #echo -e "End concat fasta to ${lspq_seq_for_nextstrain} \t$(date "+%Y-%m-%d @ %H:%M$S")\n" >> ${log_slbio}
}


FinaliseImportOnBeluga(){
  echo -e "Begin move sequence from ${beluga_server}:${beluga_seq_dir}${final_unpublished} to ${beluga_server}:${transfered_2_lspq_dir} \t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}  

  transfer_cmd="sshpass -p ${beluga_pw} ssh ${beluga_server} 'mv  ${beluga_seq_dir}${final_unpublished}/*.fasta ${beluga_seq_dir}${transfered_2_lspq_dir}' "

  eval ${transfer_cmd}
 
  echo -e "End move sequence from ${beluga_server}:${beluga_seq_dir}${final_unpublished} to ${beluga_server}:${transfered_2_lspq_dir} \t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}  
  
}


ImportSeqFromBeluga
#ExportSeqForNextstrain(){
FinaliseImportOnBeluga
exit 0





