#!/bin/sh

PYBIN=/ifs/home/c2b2/cw_lab/kje2109/apps/python27/bin/python

# job parameters
job_id=ncbi.allsubtypes.PB1.aa
K=10
M=4

# parse input
job_name=ncbi.allsubtypes.PB1.aa.K$K.M$M
prj_path=/ifs/scratch/c2b2/cw_lab/kje2109/proj/flu
src_path=$prj_path/src
log_path=$prj_path/log
tmp_path=$prj_path/tmp
dat_path=$prj_path/data

sequence_file=$dat_path/$job_id.sequences.fasta
feature_file=$dat_path/$job_id.K$K.feature_list.csv
matrix_file=$dat_path/$job_id.K$K.M$M.feature_matrix.pkl

#####################################
# step1: run gen_psmk_feature_list.py
#####################################

e1=$log_path/e.$job_name.step01
o1=$log_path/o.$job_name.step01
m1=1G;t1=:10:;
#echo "$PYBIN $src_path/gen_psmk_feature_list.py -J $job_name -I $sequence_file -K $K -O $feature_file" | qsub -l mem=${m1},time=${t1} -R y -V -N step01 -e $e1 -o $o1 


#######################################################
# step02: run gen_psmk_feature_matrix.py, hold on job01
#######################################################

# determine ranges for array jobs
seq_chunk_size=1000
feature_chunk_size=500

num_sequences=`grep ">" $sequence_file | wc -l`
num_features=`wc -l < $feature_file`

seq_range=`python -c "from math import ceil; print int(ceil($num_sequences/$seq_chunk_size.))"`
feature_range=`python -c "from math import ceil; print int(ceil($num_features/$feature_chunk_size.))"`
num_chunks=`expr $seq_range \* $feature_range`

# write chunk ranges to $chunk_file
chunk_file=$tmp_path/$job_name.chunks
if [[ -f $chunk_file ]]; then rm $chunk_file; fi
touch $chunk_file
for (( i = 1; i <= $num_chunks; i++ ))
do
    echo $i
    seq_start=`expr '(' $i - 1 ')' % $seq_range \* $seq_chunk_size`
    seq_end=`expr '(' $i - 1 ')' % $seq_range \* $seq_chunk_size + $seq_chunk_size`
    feat_start=`expr $i / $feature_range \* $feature_chunk_size`
    feat_end=`expr $i / $feature_range \* $feature_chunk_size + $feature_chunk_size`
    echo $seq_start $seq_end $feat_start $feat_end >> $chunk_file
done

e2=$log_path/e.$job_name.step02
o2=$log_path/o.$job_name.step02
m2=2G;t2=1::;
src_file=$src_path/gen_psmk_feature_matrix.py
chunk_prefix=$tmp_path/$job_name.chunked_feature_matrix
echo $num_chunks
#qsub -l mem=${m2},time=${t2} -R y -V -N step02 -e $e2 -o $o2 -hold_jid step01 -t 1-$num_chunks qsub.gen_psmk_feature_matrix.sh $src_file $job_name $sequence_file $feature_file $chunk_prefix $K $M $chunk_file


##########################################
# step3: run gen_psmk_features_assemble.py
##########################################

e3=$log_path/e.$job_name.step03
o3=$log_path/o.$job_name.step03
m3=2G;t3=1::;
echo "$PYBIN $src_path/gen_psmk_feature_assemble.py -J $job_name -SF $sequence_file -FF $feature_file -NC $num_chunks -CP $chunk_prefix -OF $matrix_file" | qsub -l mem=${m3},time=${t3} -R y -V -N step03 -e $e3 -o $o3 -hold_jid step02
