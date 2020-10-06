#!/bin/bash

#############################################################################################
# Eric Fournier 2020-10-07                                                                  #
# Tansfert de fichiers fastq de beluga vers slbio00d a partir d une liste de samples id     # 
#                                                                                           #
#############################################################################################

today=$(echo $(date "+%Y%m%d"))
hour=$(echo $(date "+%Y-%m-%d @ %H:%M$S"))

base_dir_raw_data="/data/PROJETS/COVID-19_Beluga/RawData/"

beluga_server="fournie1@beluga.computecanada.ca"
mnt_beluga_repository="/mnt/BelugaEric/"

subset="subset_2" #TODO to adjust accordingly
gisaid_dirs="20200520_20200610_20200914" #TODO to adjust accordingly
suffix="_${gisaid_dirs}_${subset}"

spec_list_file="/data/PROJETS/COVID-19_Beluga/Gisaid/SamplesListPublished${suffix}.txt"
missing_spec_list=()
missing_spec_file="${base_dir_raw_data}missing_spec_from_beluga_${today}${suffix}.txt"
log_file="${base_dir_raw_data}fastq_transfer_beluga_to_slbio00d_${today}${suffix}.log"
space_file="${base_dir_raw_data}amount_transfered_${today}${suffix}.txt"

illumina_out="${base_dir_raw_data}illumina/fastq_dehosted_${subset}/"
mgi_out="${base_dir_raw_data}mgi/fastq_dehosted_${subset}/"

echo "" > $missing_spec_file
echo -e "Start : ${hour}\n" > $log_file

echo -e "Start : ${hour}\n"

sudo umount $mnt_beluga_repository
sudo sshfs  -o allow_other -o follow_symlinks ${beluga_server}:/home/fournie1 $mnt_beluga_repository

while read sample techno
  do
  sample=$(echo $sample | cut -d '-' -f2 )
  sample=$(echo $sample | cut -d '/' -f1)
  found="False"

  if  [[ "$techno" =~ (Illumina) ]]
    then
    for fastq in $(ls ${mnt_beluga_repository}release1_dehosted_raw_reads/illumina/*.fastq.gz):
      do
      if [[ "$fastq" =~ "$sample" ]]
        then
        found="True"
        echo "********** Get Illumina fastq for ${sample} ${techno} ***************"
        echo -e "********** Get Illumina fastq for ${sample} ${techno} ***************\n" >> ${log_file}
        cp $fastq $illumina_out
      fi
    done

    if [ "${found}" = "False" ]
     then
      missing_spec_list+=("${sample}_Illumina")
    fi

  elif [[ "$techno" =~ (MGI) ]]
    then
    for fastq in $(ls ${mnt_beluga_repository}release1_dehosted_raw_reads/mgi/*.fastq.gz):
      do
      if [[ "$fastq" =~ "$sample" ]]
        then
        found="True"
         echo "********** Get MGI fastq for ${sample} ${techno} ***************"
         echo -e "********** Get MGI fastq for ${sample} ${techno} ***************\n" >> ${log_file}
         cp $fastq $mgi_out
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
  echo -e "${spec}" >> $missing_spec_file
done

hour=$(echo $(date "+%Y-%m-%d @ %H:%M$S"))
echo -e "End : ${hour}\n" >> $log_file
echo -e "End : ${hour}\n" 

du -ch  -L ${illumina_out}*.fastq.gz ${mgi_out}*.fastq.gz > ${space_file}
echo "Terminé"
