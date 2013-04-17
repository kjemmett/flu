import os
from Bio import SeqIO
from argparse import ArgumentParser

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


def generate_feature_list(seqrecords, K):
    '''
    generates a list of features from input sequence data.
    in this case, the set of unique kmers at each position

    input:
    seqrecords (obj) : SeqRecord object of sequences
    K (int) : kmer size

    output:
    feat_columns (list) : list of features
    '''

    features = []

    L = len(seqrecords[0].seq)

    for i in range(0, L-K):
        print '%s / %s' % (str(i+1), str(L-K+1))
        featlist = list(set([str(seqrecord.seq[i:i + K]) for seqrecord in seqrecords]))
        features.extend(['_'.join([str(i + 1), feature, str(i + 1 + K)])\
                for feature in featlist])

    return features


def main():

    # parse arguments
    parser = ArgumentParser()
    parser.add_argument("-J", "--jobname", dest="job_name",
        type=str, help="job name")
    parser.add_argument("-I", "--sequencefile", dest="sequence_file",
        type=str, help="input sequence file")
    parser.add_argument("-K", "--kmerlength", dest="K",
        type=int, help="kmer length")
    parser.add_argument("-O", "--output", dest="output_path",
        type=str, help="output path")
    args = parser.parse_args()

    # load sequences
    seqrecords = load_sequences(args.sequence_file)

    # generate features
    feature_list = generate_feature_list(seqrecords, args.K)

    # save
    output_file = os.path.join(args.output_path, '%s.K%02s.features.csv' \
            % (args.job_name, args.K))
    f = open(output_file, 'w')
    for feature in feature_list:
        f.write('%s\n' % feature)


if __name__=="__main__":
    main()
