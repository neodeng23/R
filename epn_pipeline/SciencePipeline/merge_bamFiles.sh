#i/bin/bash

read -p "sample list file(seq_num):" seq_file
read -p "output file name(xxx.bam):" output_file
read -p "work directory:" work_dir
all_files=""
cd $work_dir
for line in `cat $seq_file`
do
  cp /data/rawdata/bam/SEQ"$line"/SEQ"$line".bam $work_dir
  samtools sort -@ 30 SEQ${line}.bam SEQ${line}_sort
  this_file="SEQ"$line"_sort.bam"
  #this_file="/data/rawdata/bam/"$line"/"$line".bam"
  all_files="$all_files $this_file"
done
samtools merge -@ 30 $output_file $all_files
samtools index $output_file

