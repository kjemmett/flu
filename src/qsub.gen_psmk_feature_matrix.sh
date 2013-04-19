#!/bin/bash
#$ -S /bin/bash

# parse input
job_name=$1
src_file=$2
sequence_file=$3
feature_file=$4
output_file=$5.$SGE_TASK_ID.pkl
K=$6
M=$7
seq_chunks=$8
feature_chunks=$9

num_chunks=$SGE_TASK_LAST
num_sequences=`grep ">" $sequence_file | wc -l`
num_features=`wc -l < $feature_file`

seq_chunk_size=`python -c "from math import ceil; print int(ceil($num_sequences/$seq_chunks.))"`
feature_chunk_size=`python -c "from math import ceil; print int(ceil($num_features/$feature_chunks.))"`

# compute chunk ranges
seq_start=`expr '(' $SGE_TASK_ID - 1 ')' % $seq_chunks \* $seq_chunk_size`
seq_end=`expr '(' $SGE_TASK_ID - 1 ')' % $seq_chunks \* $seq_chunk_size + $seq_chunk_size`
feat_start=`expr $SGE_TASK_ID / $feature_chunks \* $feature_chunk_size`
feat_end=`expr $SGE_TASK_ID / $feature_chunks \* $feature_chunk_size + $feature_chunk_size`

/ifs/home/c2b2/cw_lab/kje2109/apps/python27/bin/python $src_file -J $job_name -IF $sequence_file -FF $feature_file -K $K -M $M -OF $output_file -SR $seq_start $seq_end -FR $feat_start $feat_end

############################################
############################################
