#!/bin/bash

#############################################################################################
# Eric Fournier 2020-10-07                                                                  #
# Tansfert de fichiers fastq de beluga vers slbio00d a partir d une liste de samples id     # 
#                                                                                           #
#############################################################################################

today=$(echo $(date "+%Y%m%d"))

beluga_server="fournie1@beluga.computecanada.ca"
mnt_beluga_repository="/mnt/BelugaEric/"

spec_list_file="/data/PROJETS/COVID-19_Beluga/Gisaid/SamplesListPublishedTEST_20200520_20200610_20200914.txt"
missing_spec_list=()
missing_spec_file="/data/PROJETS/COVID-19_Beluga/RawData/missing_spec_from_beluga_${today}.txt"
log_file="/data/PROJETS/COVID-19_Beluga/RawData/fastq_transfer_beluga_to_slbio00d_${today}.log"

echo "" > $missing_spec_file
echo "" > $log_file

sudo umount $mnt_beluga_repository
sudo sshfs  -o allow_other -o follow_symlinks ${beluga_server}:/home/fournie1 $mnt_beluga_repository

while read sample techno
  do
  sample=$(echo $sample | cut -d '-' -f2 )
  sample=$(echo $sample | cut -d '/' -f1)
  found="False"

  if  [[ "$techno" =~ (Illumina) ]]
    then
    for fastq in $(ls /mnt/BelugaEric/release1_dehosted_raw_reads/illumina/*.fastq.gz):
      do
      if [[ "$fastq" =~ "$sample" ]]
        then
        found="True"
        echo "********** Get Illumina fastq for ${sample} ${techno} ***************"
        echo -e "********** Get Illumina fastq for ${sample} ${techno} ***************\n" >> ${log_file}
        #cp $fastq /data/PROJETS/COVID-19_Beluga/RawData/illumina/fastq_dehosted
      fi
    done

    if [ "${found}" = "False" ]
     then
      missing_spec_list+=("${sample}_Illumina")
    fi

  elif [[ "$techno" =~ (MGI) ]]
    then
    for fastq in $(ls /mnt/BelugaEric/release1_dehosted_raw_reads/mgi/*.fastq.gz):
      do
      if [[ "$fastq" =~ "$sample" ]]
        then
        found="True"
         echo "********** Get MGI fastq for ${sample} ${techno} ***************"
         echo -e "********** Get MGI fastq for ${sample} ${techno} ***************\n" >> ${log_file}
         #cp $fastq /data/PROJETS/COVID-19_Beluga/RawData/mgi/fastq_dehosted
      fi
    done

    if [ "${found}" = "False" ]
     then
      missing_spec_list+=("${sample}_MGI")
    fi
  else
      missing_spec_list+=("${sample}_Nanopore")
  fi
done < ${spec_list_file}

for spec in ${missing_spec_list[@]}
  do
  echo -e "${spec}\n" >> $missing_spec_file
done

echo "Terminé"
