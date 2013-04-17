#!/bin/bash

# parse input
job_name=ncbi.allsubtypes.PB1.aa
prj_path=/ifs/scratch/c2b2/cw_lab/kje2109/proj/flu
src_path={$prj_path}/src
log_path={$prj_path}/log

# log files
ef=$log_path/e
of=$log_path/o

sequence_file={$prj_path}/data/{$job_name}.sequences.fasta
feature_file={$prj_path}/data/{$job_name}.K{$K}.feature_list.csv
matrix_file={$prj_path}/data/${job_name}.K{$K}.M{$M}.feature_matrix.pkl

PYBIN=/ifs/home/c2b2/cw_lab/kje2109/apps/python27/bin/python


# determine ranges for array jobs
seq_chunk_size=100
feature_chunk_size=100

K=10

# step1: run gen_psmk_feature_list.py
m1=1G;t1=:10:;
echo "$PYBIN $src_path/gen_psmk_feature_list.py -J $job_name -I $sequence_file -K $K -O $feature_file" | qsub -l mem=${m1},time=${t1} -R y -V -N job01 -e $ef -o $of 

# step2: run gen_psmk_feature_matrix.py, hold on job1
#echo "$PYBIN $src_path/gen_psmk_feature_matrix.py" | qsub -l mem=${m2},time=${t2} -R y -V -N job02 -e $ef -o $of -hold_jid job01

# step3: run gen_psmk_features_assemble.py
#echo "$PYBIN " | qsub -l mem=${m3},time=${t3} -R y -V -N job03 -e $ef -o $of -hold_jid 
