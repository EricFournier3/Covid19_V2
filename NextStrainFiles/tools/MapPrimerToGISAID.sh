#!/bin/bash

echo "Start primer mapping"

index_ref_cmd="smalt index -k 13 -s 2 WuhanPrimerRegion WuhanPrimerRegion.fasta"
eval $index_ref_cmd

map_cmd="smalt map -f sam -o primer_map.sam WuhanPrimerRegion primer.fasta"
eval $map_cmd

sam2bam_cmd="samtools view -bS primer_map.sam  > primer_map.bam"
eval $sam2bam_cmd

sort_cmd="samtools sort primer_map.bam   primer_map_sort"
eval $sort_cmd

index_cmd="samtools index primer_map_sort.bam"
eval $index_cmd

echo "End"
