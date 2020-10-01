#!/bin/bash

######################################################################
# Eric Fournier 2020-10-01                                           #
# Tansfert de fichiers fastq a partir d une liste de samples id      # 
# Script to run on fournie1@beluga.computecanada.ca                                                                    #
######################################################################


out_dir="/home/fournie1/test2/"
out_duplicate_sample=${out_dir}"DuplicateSamples.txt"

beluga_server="fournie1@beluga.computecanada.ca"



while read sample techno
  do
  sample=$(echo $sample | cut -d '-' -f2 )
  sample=$(echo $sample | cut -d '/' -f1)
  
  echo "SAMPLE IS ${sample}"
 
  if [[ "$techno" =~ (Illumina) ]]
    then
    techno="illumina"
  elif  [[ "$techno" =~ (ONT_ARTIC) ]]
    then
    techno="nanopore"
  elif [[ "$techno" =~ (MGI) ]]
    then
    techno="mgi"
  else
    :
    echo "****************** Bug get techno for $sample *********************"
    continue 
  fi

  n=0


  for i  in $(find "./COVID_LSPQ/REPOSITORY/" -name "$sample"*fastq.gz | grep  "$sample"'.'"$techno")
    do
    n=$((n+1))

    echo "copy Symlink for $i"
    
    cp  -d $i $out_dir
  done


  if [ "$techno" == "illumina" ] || [ "$techno" == "mgi" ]
    then
    if [ $n -gt 2 ]
      then
      echo -e $(find "./COVID_LSPQ/REPOSITORY/" -name "$sample"*fastq.gz | grep  "$sample"'.'"$techno") "\n" >> $out_duplicate_sample
    fi
  elif [ "$techno" == "nanopore" ]
    then
    if [ $n -gt 1 ]
      then
      echo -e $(find "./COVID_LSPQ/REPOSITORY/" -name "$sample"*fastq.gz | grep  "$sample"'.'"$techno") "\n"   >> $out_duplicate_sample
    fi
  fi


done < SamplesListTEST2.txt

echo "Terminé"
