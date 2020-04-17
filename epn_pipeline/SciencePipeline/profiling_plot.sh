#!/bin/bash

read -p "control bam/merged bam file(xxx.bam):" input_file_1
read -p "treated bam/merged bam file(xxx.bam):" input_file_2
read -p "before region start length(b):" upstream_kb
read -p "region body length(b):" region_length
read -p "after region start length(b):" downstream_kb
read -p "bin size(b):" bin_size
read -p "output file prefix:" output_file_prefix

bamCoverage -b $input_file_1 -o ${input_file_1%.bam*}.bw
bamCoverage -b $input_file_2 -o ${input_file_2%.bam*}.bw

computeMatrix scale-regions -S ${input_file_1%.bam*}.bw ${input_file_2%.bam*}.bw -R /database/hg19/hg19_refseq -b $upstream_kb -m $region_length -a $downstream_kb --skipZeros -bs $bin_size -p 10 -o $output_file_prefix".mat.gz"
plotProfile -m $output_file_prefix".mat.gz" -o $output_file_prefix"_TSS_TES.pdf" --perGroup
#plotProfile -m $output_file_prefix".mat.gz" -o $output_file_prefix"_TSS_TES.pdf" --perGroup  --samplesLabel control treated