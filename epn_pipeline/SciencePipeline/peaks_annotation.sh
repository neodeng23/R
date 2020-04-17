#!/bin/bash

read -p "bam/merged bam file(xxx.bam):" input_bam_file
read -p "output file prefix:" output_file_prefix

macs2 callpeak -t $input_bam_file -f BAMPE -n $output_file_prefix --nomodel --nolambda --broad
annotatePeaks.pl $output_file_prefix"_peaks.broadPeak" hg19 > $output_file_prefix"_peaks_anno.bed"
python /database/scripts/peak_annotation_pie_chart.py -i $output_file_prefix"_peaks_anno.bed" -o $output_file_prefix