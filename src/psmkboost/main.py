import pdb
from Bio import SeqIO
import pandas as pd
import numpy as np
import cPickle
import splitdata
import boost
import sys
import os
import csv
import time
from argparse import ArgumentParser

def gen_label_matrix(lf, N, C):
    '''
    generates the label matrix of appropriate dimension

    inputs:
    lf (str) : path to label csv
    C (int) : number of classes

    outputs:
    Y (NxC) : label matrix
    '''

    Y = -1 * np.ones((N,C),dtype='int')

    f = open(lf,'r')
    i = 0

    for row in csv.reader(open(lf,'r'),delimiter=','):
        for l in row[1:]:
            Y[i,int(l)-1] = 1
        i += 1

    return Y

def main():

    # default proj_path
    def_pp = os.sep.join(os.getcwd().split(os.sep)[:-1])

    # parse arguments
    parser = ArgumentParser()
    parser.add_argument("-J","--jobid", dest="jobid",
            type=str, help="dataset id")
    parser.add_argument("-K", "--kmerlength", dest="K",
            type=int, help="kmer feature length")
    parser.add_argument("-M", "--mismatch", dest="M",
            type=int, help="max allowed mismatch")
    parser.add_argument("-P", "--projpath", dest="proj_path",
            type=str, default=def_pp, help="project path (defaults to cwd)")
    parser.add_argument("-T","--rounds", dest="T",
            type=int, help="number of boosting rounds")
    parser.add_argument("-N","--folds", dest="Nfolds",
            type=int, default=5, help="number of cv folds to execute")
    parser.add_argument("-R","--runid", dest="runid",
            type=str, help="run id")
    parser.add_argument("-O","--outdir",dest="outdir",
            type=str, help="output directory")
    parser.add_argument("-Z","--model", dest="model",
            type=str, default="tree", help="model type (tree or stump)")
    args = parser.parse_args()

    jobid = args.jobid
    runid = args.runid
    K = args.K
    M = args.M
    T = args.T
    Nfolds = args.Nfolds
    model = args.model
    outdir = args.outdir

    # set up paths
    proj_path = args.proj_path
    src_path = proj_path + '/src/psmkboost'
    C_path = '%s/get_new_function.c' % (src_path)
    feat_path = '%s/data' % (proj_path)
    # run_path = '%s/cache/runs/%s' % (proj_path, runid)
    run_path = outdir
    if not os.path.exists(run_path):
        os.makedirs(run_path)

    print 'Running adaboost on %s using position specific mismatch kmer feature space' % (jobid)

    # load feature matrix
    XDF = pd.load('%s/data/%s.K%d.M%d.feature_matrix.pkl' % (proj_path, jobid, K, M))
    (N, P) = XDF.shape

    # load feature list
    features = XDF.columns
    feat_dict = dict(zip(range(P), features))

    # load label dict
    ldf = '%s/data/%s.labeldict.csv' % (proj_path, jobid)
    label_dict = {}
    for row in csv.reader(open(ldf,'r'), delimiter=','):
        label_dict[int(row[0])] = row[1]
    C = len(label_dict)

    # load label matrix
    lf = '%s/data/%s.labels.csv' % (proj_path, jobid)
    Y = gen_label_matrix(lf, N, C)
    Yt = Y.T

    # in this case we will loop over binary thresholds
    threshold_list = range(2)

    #holds predicted label at each round
    predicted_labels = np.zeros((N,T),dtype='int')

    # split the data indices into `Nfold` random disjoint sets
    Fidx = splitdata.cv_multiclass_fold(Yt,Nfolds)

    for fold in range(Nfolds):

        print 'executing fold %d'%(fold+1)
        # split the data and labels into train and test sets
        # skip this and return only indices?
        print 'splitting data...'
        tt = time.time()
        train_data, train_labels, test_data, test_labels \
            = splitdata.cv_split(XDF.as_matrix().T, Yt, Fidx[fold])
        print 'split data time=%.2f seconds'%(time.time()-tt)

        # specify output file names
        filetag = '%s_K%d_M%d_fold%d'%(model,K,M,fold)
        output_file = '%s/%s.%s.outputsummary_%s.txt' % (run_path, jobid, runid, filetag)
        handle = open(output_file,'w')
        to_write = ['round', 'kmer', 'threshold', 'train_auc', 
                    'train_acc', 'test_auc', 'test_acc', 'runtime']
        handle.write('\t'.join(to_write)+'\n')
        handle.close()

        # parse the C code from get_new_function.c
        f = open(C_path,'r')
        C_code = '\n'.join([line for line in f if '//' not in line])
        f.close()

        pdb.set_trace()
        # run Adaboost
        print 'entering adaboost'
        adt, adt_outputs, performance, predicted_labels = boost.adaboost(C_code, \
            train_data, train_labels, test_data, test_labels, T, \
            output_file=output_file, kmer_dict=feat_dict, model=model, \
            predicted_labels=predicted_labels, test_indices=Fidx[fold])

        # save the learned model
        model_file = '%s/%s.%s.adt_%s.pkl' % (run_path, jobid, runid, filetag)
        handle = open(model_file,'w')
        cPickle.dump(adt,handle)
        handle.close()

        # save algorithm performance (errors, runtime, etc)
        results_file = '%s/%s.%s.performance_%s.pkl' % (run_path, jobid, runid, filetag)
        handle = open(results_file,'w')
        cPickle.Pickler(handle,protocol=2).dump(adt_outputs)
        cPickle.Pickler(handle,protocol=2).dump(performance)
        handle.close()

    # output predicted labels on test data for each CV fold
    output_file = '%s/%s.%s.testsetpredictions_%d.pkl' \
        % (run_path, jobid, runid, K)
    handle = open(output_file,'w')
    cPickle.Pickler(handle,protocol=2).dump(Fidx)
    cPickle.Pickler(handle,protocol=2).dump(predicted_labels)
    handle.close()
            
if __name__=="__main__":
    main()
