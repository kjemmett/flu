#!/bin/sh

PYBIN=/ifs/home/c2b2/cw_lab/kje2109/apps/python27/bin/python

# job parameters
job_id=ncbi.allsubtypes.PB2.aa
#K=4
#M=2

for K in 4; do
    M_MAX=`expr $K / 2`
    for (( M = 0; M <= $M_MAX; M++ )); do

        # parse input
        job_name=$job_id.K$K.M$M
        prj_path=/ifs/scratch/c2b2/cw_lab/kje2109/proj/flu
        src_path=$prj_path/src
        log_path=$prj_path/log
        tmp_path=$prj_path/tmp
        dat_path=$prj_path/data/$job_id
        chunk_path=$tmp_path/$job_name
        mkdir -p $chunk_path

        # relevant filepaths
        sequence_file=$dat_path/$job_id.sequences.fasta
        feature_file=$dat_path/$job_id.K$K.feature_list.csv
        matrix_file=$dat_path/$job_id.K$K.M$M.feature_matrix.pkl

        #####################################
        # step1: run gen_psmk_feature_list.py
        #####################################


        e1=$log_path/e.$job_name.step01
        o1=$log_path/o.$job_name.step01
        m1=1G;t1=:10:;
        echo "$PYBIN $src_path/gen_psmk_feature_list.py -J $job_name -I $sequence_file -K $K -O $feature_file" | qsub -l mem=${m1},time=${t1} -R y -V -N $job_name.step01 -e $e1 -o $o1 


        #######################################################
        # step02: run gen_psmk_feature_matrix.py, hold on job01
        #######################################################

        # set chunk sizes
        # feature_chunk_size is set based on num_chunks

        seq_chunks=10
        feature_chunks=100
        num_chunks=`expr $seq_chunks \* $feature_chunks`

        e2=$log_path/e.$job_name.step02
        o2=$log_path/o.$job_name.step02
        m2=2G;t2=1::;
        src_file=$src_path/gen_psmk_feature_matrix.py
        
        chunk_prefix=$chunk_path/$job_name.chunked_feature_matrix
        export MPLCONFIGDIR=/ifs/scratch/c2b2/cw_lab/kje2109/.matplotlib/
        qsub -l mem=${m2},time=${t2} -R y -V -N $job_name.step02 -e $e2 -o $o2 -hold_jid $job_name.step01 -t 1-$num_chunks qsub.gen_psmk_feature_matrix.sh $job_name $src_file $sequence_file $feature_file $chunk_prefix $K $M $seq_chunks $feature_chunks


        ##########################################
        # step3: run gen_psmk_features_assemble.py
        ##########################################

        e3=$log_path/e.$job_name.step03
        o3=$log_path/o.$job_name.step03
        m3=2G;t3=1::;
        echo "$PYBIN $src_path/gen_psmk_feature_assemble.py -J $job_name -SF $sequence_file -FF $feature_file -NC $num_chunks -CP $chunk_prefix -OF $matrix_file" | qsub -l mem=${m3},time=${t3} -R y -V -N $job_name.step03 -e $e3 -o $o3 -hold_jid $job_name.step02
    done
done
