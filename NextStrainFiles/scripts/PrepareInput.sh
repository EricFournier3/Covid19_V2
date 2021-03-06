#!/bin/bash

#TODO wget au lieu de copy


_debug="False"


work_dir=$1"/"
temp_dir=${work_dir}"temp/"
config_dir=${work_dir}"config"/
data_dir=${work_dir}"data/"
script_dir=${work_dir}"scripts/"
gisaid_dir=${work_dir}"gisaid/"
lspq_dir=${work_dir}"lspq/"
init_file=${work_dir}"init.txt"
create_new_config_script=${script_dir}"CreateNextstrainConfigV3.py"
prepare_metadata_script=${script_dir}"PrepareMetadataV2.py"

lat_long_out=${config_dir}"lat_longs.tsv"
ordering_out=${config_dir}"ordering.tsv"

metadata_out=${data_dir}"metadata.tsv"
sequences_out=${data_dir}"sequences.fasta"

nextstrain_files_base_dir="${2}"

nb_lspq_seq_to_keep=$3
echo ">>>>>>>>>>>>>>>> $nb_lspq_seq_to_keep "
nb_gisaid_seq_to_keep=$4

nb_canadian_gisaid_seq_to_keep=$5

root=$6
max_sample_date=$7
minimum_qc_status=$8
original_data_dir=$9
in_data_version=${10}
min_sample_date=${11}


if [ ${minimum_qc_status} == "FLAG" ]
  then
  minimum_qc_status="PASS_FLAG"
fi
    

BuildFramework(){
   echo "In BuildFrameWork"

   touch ${init_file}
   mkdir ${temp_dir}
   mkdir ${data_dir}
   mkdir ${gisaid_dir}
   mkdir ${lspq_dir}  
}

TransferConfigFiles(){
  #TODO wget
  ext=(".txt" ".tsv" ".json" ".gb" ".md" ".yaml")

  for _ext in ${ext[@]}
    do
    cp "${nextstrain_files_base_dir}config/"*"${_ext}" ${config_dir}
  done 
  
  cp ${temp_dir}"lat_longs.tsv" ${config_dir}
}

TransferScripts(){
  ext=(".py" ".sh")

  for _ext in ${ext[@]}
    do
    cp "${nextstrain_files_base_dir}scripts/"*"${_ext}" ${script_dir}
  done
}


ImportLatLong(){
  #wget https://raw.githubusercontent.com/nextstrain/ncov/master/config/lat_longs.tsv -O ${temp_dir}"lat_longs.tsv" 1>/dev/null  2>&1
   cp "${nextstrain_files_base_dir}config/lat_longs.tsv" ${temp_dir}
}

CreateNewConfigFiles(){

  lat_long_ori="${config_dir}lat_longs.tsv"
  missing_sgil_countries=${config_dir}"MissingSgilCountries.tsv"
  country_lat_long="${config_dir}country_lat_long.tsv"
  division_lat_long="${config_dir}division_lat_long.tsv"

  awk 'BEGIN{FS="\t"}/country/{print $2"\t"$3"\t"$4}' ${lat_long_ori} > ${country_lat_long}
  cat ${missing_sgil_countries} >> ${country_lat_long}
  awk 'BEGIN{FS="\t"}/division/{print $2"\t"$3"\t"$4}' ${lat_long_ori} > ${division_lat_long}

  python ${create_new_config_script} --lat-long-out ${lat_long_out} --ordering-out ${ordering_out} --base-dir ${work_dir} 
}

TransferGisaidFiles(){
   _dir="${nextstrain_files_base_dir}data/gisaid/"
   fasta_path=$(ls ${_dir}all/*.fasta)
   fasta_name=$(basename ${fasta_path})
   metadata_path=$(ls ${_dir}all/*.tsv)
   metadata_name=$(basename ${metadata_path})
   echo -e "Fasta version : ${fasta_name}\n" > ${gisaid_dir}"data_version.txt"
   echo -e "Metadata version : ${metadata_name}\n" >> ${gisaid_dir}"data_version.txt"
   #cp "${_dir}"{"all/gisaid_all.fasta","all/metadata.tsv","gisaid_wuhan_ref_20200425.fasta"} ${gisaid_dir}

   cp ${fasta_path} ${gisaid_dir}"gisaid_all.fasta"
   cp ${metadata_path} ${gisaid_dir}"metadata.tsv"
   cp ${_dir}"gisaid_wuhan_ref_20200425.fasta" ${gisaid_dir}
}

TransferLspqFiles(){

    if [ "$_debug" = "False" ]
     then
      nextstrain_lspq_data_dir="${nextstrain_files_base_dir}data/lspq/"
    else
      nextstrain_lspq_data_dir="${nextstrain_files_base_dir}data/lspq/debug/"
    fi

    original_metadata=${original_data_dir}"METADATA/metadata_${in_data_version}_${minimum_qc_status}_minmaxSampleDate_${min_sample_date}_${max_sample_date}.tsv"
    original_sequences=${original_data_dir}"FASTA/sequences_${in_data_version}_${minimum_qc_status}_minmaxSampleDate_${min_sample_date}_${max_sample_date}.fasta"
 
    if ! [ -f ${original_metadata} ]
        then
        echo "${original_metadata} est absent"
        exit 1
    fi

    if ! [ -f ${original_sequences} ]
        then
        echo "${original_sequences} est absent"
        exit 1
    fi

    lspq_data_version_file=${lspq_dir}"data_version.txt"

    original_metadata_basename=$(echo $(basename ${original_metadata}))
    original_sequences_basename=$(echo $(basename ${original_sequences}))

    echo -e "Metadata version : ${original_metadata_basename}" > ${lspq_data_version_file} 
    echo -e "Sequences version : ${original_sequences_basename}" >> ${lspq_data_version_file} 
 
    cp ${original_metadata} ${lspq_dir}"metadata.tsv"
    cp ${original_sequences} ${lspq_dir}"sequences.fasta"
    cp "${nextstrain_lspq_data_dir}""sequences_force_include.fasta" ${lspq_dir}

    #ligne ci-dessous obsolete
    #cp "${nextstrain_lspq_data_dir}"{"sequences.fasta","sequences_force_include.fasta","sgil_extract.tsv"} ${lspq_dir}
}

PrepareMetadata(){

  if [ "${root}" = "Wuhan/WH01/2019" ]
    then
    python ${prepare_metadata_script} --metadata-out ${metadata_out}   --sequences-out ${sequences_out}  --base-dir ${work_dir} --wuhan-root  --metadata-lspq "${lspq_dir}metadata.tsv"
  else
    python ${prepare_metadata_script} --metadata-out ${metadata_out}   --sequences-out ${sequences_out}  --base-dir ${work_dir} --metadata-lspq "${lspq_dir}metadata.tsv"
    sed -i "s/Wuhan\/WH01\/2019/${root//\//\\/}/g" ${work_dir}Snakefile
  fi
}

ExcludeSamples(){

 mv ${gisaid_dir}"gisaid_all.fasta"   ${gisaid_dir}"gisaid_all_temp.fasta"
 mv ${lspq_dir}"sequences.fasta" ${lspq_dir}"sequences_temp.fasta"

 non_canadian_seq=${gisaid_dir}"non_canadian_seq.fasta"
 canadian_seq=${gisaid_dir}"canadian_seq.fasta"

 seqkit grep -r -p 'Canada' ${gisaid_dir}"gisaid_all_temp.fasta" | seqkit grep -r -p 'Canada/Qc-' -v | seqkit grep -r -p 'Canada/QC_' -v > ${canadian_seq} 
 seqkit grep -r -p 'Canada' ${gisaid_dir}"gisaid_all_temp.fasta" -v > ${non_canadian_seq} 

 if ! [ "${nb_canadian_gisaid_seq_to_keep}" == "all" ]
   then
   seqtk sample -s $RANDOM "${canadian_seq}" ${nb_canadian_gisaid_seq_to_keep} > ${gisaid_dir}"gisaid_all.fasta"
 else
   cat "${canadian_seq}" > ${gisaid_dir}"gisaid_all.fasta"
 fi

 seqtk sample -s $RANDOM "${non_canadian_seq}" ${nb_gisaid_seq_to_keep} >> ${gisaid_dir}"gisaid_all.fasta"
 seqtk sample -s $RANDOM ${lspq_dir}"sequences_temp.fasta" ${nb_lspq_seq_to_keep} > ${lspq_dir}"sequences_temp2.fasta"

 cat ${lspq_dir}"sequences_temp2.fasta" ${lspq_dir}"sequences_force_include.fasta" > ${lspq_dir}"sequences_temp3.fasta"
 seqkit rmdup -n ${lspq_dir}"sequences_temp3.fasta" > ${lspq_dir}"sequences.fasta"

 n_qc_lspq=$(grep ">" ${lspq_dir}"sequences.fasta" | wc -l)

 echo "Number of Qc sequences from LSPQ is ${n_qc_lspq}"

 rm ${canadian_seq}
 rm ${non_canadian_seq}
 rm ${gisaid_dir}"gisaid_all_temp.fasta"
 rm ${lspq_dir}"sequences_temp.fasta" ${lspq_dir}"sequences_temp2.fasta"  ${lspq_dir}"sequences_temp3.fasta"


}


if ! [ -f ${init_file} ]
  then
  BuildFramework
  ImportLatLong
  TransferGisaidFiles
  TransferLspqFiles
  ExcludeSamples
fi

TransferConfigFiles
TransferScripts
CreateNewConfigFiles
PrepareMetadata


