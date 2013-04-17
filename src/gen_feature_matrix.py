import cPickle
import numpy as np
import pandas as pd
from Bio import SeqIO
from argparse import ArgumentParser

def calc_mismatch(str1, str2):
    '''
    calculates mismatch between two strings

    input:
    str1 : string one
    str2 : string two

    output:
    m : mismatch
    '''

    return sum(ch1 != ch2 for ch1, ch2 in zip(str1, str2))


def load_sequences(ff):
    '''
    loads sequences from fastas file

    input:
    ff (str) : path to fasta file containing sequences

    output (obj) : SeqRecord object containing sequence information
    '''

    sh = open(ff, 'r')
    seq = list(SeqIO.parse(sh, "fasta"))
    sh.close()
    return seq


def generate_feature_matrix(seqrecords, K, M, features):
    '''
    generates psmk feature matrix from set of
    seq records

    input:
    seqrecords (obj) : SeqRecord Object

    output:
    X (array) : numpy array (size a priori unknown
    '''

    # initialize feature matrix
    print 'initializing feature matrix'
    seqids = [seq.id for seq in seqrecords]
    X = pd.DataFrame(index=seqids, columns=features, dtype='bool')

    print 'computing mismatches'
    for idx, feature in enumerate(features):
        # parse feature
        start, kmer, stop = feature.split('_')
        i = int(start) - 1
        # compute mismatches
        mval = [calc_mismatch(kmer, str(seq.seq[i:i + K])) for \
                seq in seqrecords]
        X[feature] = np.array(map(lambda x: x <= M, mval))

    return X


def main():

    # parse arguments
    parser = ArgumentParser()
    parser.add_argument("-J", "--jobname", dest="job_name",
        type=str, help="job name")
    parser.add_argument("-IF", "--sequencefile", dest="sequence_file",
        type=str, help="input sequence file")
    parser.add_argument("-FF", "--featurefile", dest="feature_file",
        type=str, help="intput feature file")
    parser.add_argument("-K", "--kmerlength", dest="K",
        type=int, help="kmer length")
    parser.add_argument("-M", "--mismatch", dest="M",
        type=int, help="mismatch value")
    parser.add_argument("-OF", "--output", dest="output_file",
        type=str, help="output path")
    parser.add_argument("-SR", "--seqrange", nargs=2, dest="seqrange",
        type=int, help="range of sequences (default=all)")
    parser.add_argument("-FR", "--featurerange", nargs=2, dest="featurerange",
        type=int, help="range of eatures (default=all)")
    args = parser.parse_args()

    # load sequences
    seqrecords = load_sequences(args.sequence_file)

    # slice requested seqrange (if necessary)
    if args.seqrange is not None:
        if args.seqrange[1] > len(seqrecords):
            args.seqrange[1] = len(seqrecords)
        start = args.seqrange[0]
        stop = args.seqrange[1]
        seqrecords = seqrecords[start:stop]

    # load features
    ff = open(args.feature_file, 'r').readlines()
    features = [fi.strip() for fi in ff]

    if args.featurerange is not None:
        if args.featurerange[1] > len(features):
            args.featurerange[1] = len(features)
        features = features[args.featurerange[0]:args.featurerange[1]]

    # compute feature matrix
    X = generate_feature_matrix(seqrecords, args.K, args.M, features)

    # save data
    cPickle.Pickler(open(args.output_file, 'w'), protocol=2).dump(X)


if __name__=="__main__":
    main()
