#!/bin/bash

green_message="\e[32m"
white_message="\e[39m"

base_dir_slbio="/data/Runs/SARS-CoV-2/GenomeCenterSeq/FinalRelease"

beluga_server="fournie1@beluga.computecanada.ca"
beluga_seq_dir="/home/fournie1/test"

slbio_user=$(whoami)
beluga_pass_file="/home/${slbio_user}/compute_canada_pass.txt"
read beluga_pw < ${beluga_pass_file}

scp_cmd="sshpass -p ${beluga_pw}  scp -r ${beluga_server}:${beluga_seq_dir} ${base_dir_slbio}"

echo -e "${green_message}INFO: " "Begin import from ${beluga_server}:${beluga_seq_dir}"
eval ${scp_cmd}
echo -e "${green_message}INFO: " "End import from ${beluga_server}:${beluga_seq_dir}"

echo -e "${white_message}"
