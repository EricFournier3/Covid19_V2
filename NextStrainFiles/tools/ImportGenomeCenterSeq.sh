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

out_seq_slbio=${base_dir_slbio}"FinalRelease/"
log_slbio="${out_seq_slbio}BelugaImport.log"

beluga_server="fournie1@beluga.computecanada.ca"
beluga_seq_dir="/home/fournie1/COVID_consensus/"

slbio_user=$(whoami)
beluga_pass_file="/home/${slbio_user}/compute_canada_pass.txt"
read beluga_pw < ${beluga_pass_file}

final_unpublished="FinalUnpublished"
final_published="FinalPublished"
final_submitted="FinalSubmitted"
final_do_not_publish="DoNotPublish"
transfered_2_lspq_dir=${final_unpublished}"/Transfered2Lspq/"

ImportSeqFromBeluga(){
  echo -e "Begin import from ${beluga_server}:${beluga_seq_dir}\t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}

  #scp_cmd_obsolete="sshpass -p ${beluga_pw}  scp -r ${beluga_server}:\"${beluga_seq_dir}{${final_unpublished}/*.fasta,${final_published}/*.fasta}\" ${out_seq_slbio}"
  #scp_cmd_obsolete_2="sshpass -p ${beluga_pw}  scp -r ${beluga_server}:\"${beluga_seq_dir}{${final_published}/,${final_unpublished}/}\" ${out_seq_slbio}"

  
  scp_cmd="sshpass -p ${beluga_pw}  scp -r ${beluga_server}:\"${beluga_seq_dir}${final_unpublished}/*.fasta\" ${out_seq_slbio}${final_unpublished}"

  eval ${scp_cmd}

  echo -e "End import from ${beluga_server}:${beluga_seq_dir}\t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}

}


ExportSeqForNextstrain(){
  echo "In ExportSeqForNextstrain"
  echo -e "Begin concat fasta to ${lspq_seq_for_nextstrain} \t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}

  cat ${out_seq_slbio}{${final_unpublished}/L*.fasta,${final_published}/L*.fasta,${final_submitted}/L*.fasta} >${lspq_seq_for_nextstrain}

  echo -e "End concat fasta to ${lspq_seq_for_nextstrain} \t$(date "+%Y-%m-%d @ %H:%M$S")\n" >> ${log_slbio}
}


SelectSamplesOnBeluga(){
  echo -e "Begin select sequences on ${beluga_server}  \t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}  

  transfer_cmd="sshpass -p ${beluga_pw} ssh ${beluga_server} /home/fournie1/SelectSamplesForLSPQimport.sh"

  eval ${transfer_cmd}
  
  echo -e "End select sequences on ${beluga_server}  \t$(date "+%Y-%m-%d @ %H:%M$S")" >> ${log_slbio}  
}

MakeAlreadyTransferList(){
  seq_list=()
  
  already_transfered_file=${out_seq_slbio}"already_transfered_file.txt"
 
  if [ -f $already_transfered_file ]
      then
      rm $already_transfered_file
      touch $already_transfered_file
  else
     touch $already_transfered_file
  fi


  for seq in $(ls "${out_seq_slbio}${final_unpublished}/L"*".fasta" 2>/dev/null)
    do
    seq=$(basename $seq | cut -d '.' -f1)
    #seq_list+=($seq)
    echo $seq >> ${already_transfered_file}
  done

  for seq in $(ls "${out_seq_slbio}${final_published}/L"*".fasta" 2>/dev/null)
    do
    seq=$(basename $seq | cut -d '.' -f1)
    #seq_list+=($seq)
    echo $seq >> ${already_transfered_file}
  done

  for seq in $(ls "${out_seq_slbio}${final_submitted}/"*"/L"*".fasta" 2>/dev/null)
    do
    seq=$(basename $seq | cut -d '.' -f1)
    #seq_list+=($seq)
    echo $seq >> ${already_transfered_file}
  done

  for seq in $(ls "${out_seq_slbio}${final_do_not_publish}/L"*".fasta" 2>/dev/null)
    do
    seq=$(basename $seq | cut -d '.' -f1)
    #seq_list+=($seq)
    echo $seq >> ${already_transfered_file}
  done

  scp_cmd="sshpass -p ${beluga_pw}  scp   ${already_transfered_file}    ${beluga_server}:${beluga_seq_dir}"
  eval $scp_cmd
}

MakeAlreadyTransferList
SelectSamplesOnBeluga
ImportSeqFromBeluga
ExportSeqForNextstrain

echo "Finish"




