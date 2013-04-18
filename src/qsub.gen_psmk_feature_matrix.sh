#!/bin/bash
#$ -S /bin/bash

src_file=$1
job_name=$2
sequence_file=$3
feature_file=$4
output_file=$5.$SGE_TASK_ID.pkl
K=$6
M=$7
chunk_file=$8

echo $src_file
echo $job_name
echo $sequence_file
echo $feature_file
echo $output_file
echo $K
echo $M
echo $chunk_file

# read range for chunks
foo=`sed -n "$SGE_TASK_ID"p $chunk_file`
IFS=' ' read -a bar <<< $foo
seq_start=${bar[0]}
seq_end=${bar[1]}
feat_start=${bar[2]}
feat_end=${bar[3]}

/ifs/home/c2b2/cw_lab/kje2109/apps/python27/bin/python $src_file -J $job_name -IF $sequence_file -FF $feature_file -K $K -M $M -OF $output_file -SR $seq_start $seq_end -FR $feat_start $feat_end

############################################
############################################
