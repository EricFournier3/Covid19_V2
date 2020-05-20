#!/bin/bash

<<COM
Eric Fournier 20200519
Script to run on fournie1@beluga.computecanada.ca
COM

my_home_basedir="/home/fournie1/COVID_consensus/"
all_final_unpublished_dir="/genfs/projects/COVID_consensus/FinalUnpublished/"

already_import_to_lspq_files=${my_home_basedir}"already_transfered_file.txt"
already_import_lspq_seq=()

final_ready_to_import_dir=${my_home_basedir}"FinalUnpublished/"

BuildAlreadyImportSeqList(){
  while read seq
   do
   already_import_lspq_seq+=($seq)
  done <  ${already_import_to_lspq_files}	   
 
  echo ${already_import_lspq_seq[@]}

}


BuildAlreadyImportSeqList


for fasta in $(ls ${all_final_unpublished_dir})
  do 
	  seq_name=$(echo $fasta | cut -d '.' -f1)

	  if ! [[ "${already_import_lspq_seq[@]}" =~  "${seq_name}" ]]
		 then
	         cp ${all_final_unpublished_dir}$fasta ${final_ready_to_import_dir}
	  fi
  done
