import os
import pdb
from argparse import ArgumentParser
from Bio import SeqIO
import pandas as pd
import numpy as np

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


def main():

    # parse arguments
    parser = ArgumentParser()
    parser.add_argument("-J", "--jobname", dest="job_name",
        type=str, help="job name")
    parser.add_argument("-SF", "--sequencefile", dest="sequence_file",
        type=str, help="sequence file")
    parser.add_argument("-FF", "--featurefile", dest="feature_file",
        type=str, help="feature file")
    parser.add_argument("-NC", "--numchunks", dest="num_chunks",
        type=int, help="number of chunks")
    parser.add_argument("-CP", "--chunkprefix", dest="chunk_prefix",
        type=str, help="input chunk prefix")
    parser.add_argument("-OF", "--output", dest="output_file",
        type=str, help="output file")
    parser.add_argument("-DC", "--delete_chunks", dest="delete_chunks",
        type=int, help="delete chunks? (0 = no, 1 = yes (default)",
        default=1)
    args = parser.parse_args()

    # load sequences
    seqrecords = load_sequences(args.sequence_file)
    seqids = [seq.id for seq in seqrecords]

    # load features
    ff = open(args.feature_file, 'r').readlines()
    ff = [fi.strip() for fi in ff]

    # initialize dataframe
    X = pd.DataFrame(columns=ff, index=seqids, dtype='bool')
    
    # load chunks and place in appropriate spot
    for i in range(1, args.num_chunks+1):
        chunk_file = '%s.%d.pkl' % (args.chunk_prefix, i)
        chunk = pd.load(chunk_file) 
        X.ix[chunk.index, chunk.columns] = chunk.as_matrix()
        if args.delete_chunks:
            os.remove(chunk_file)

    # save
    pd.save(X, args.output_file)

if __name__=="__main__":
    main()
