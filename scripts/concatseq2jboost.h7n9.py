from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import csv
import pdb


if __name__=="__main__":

    # segments to concatenate
    projpath = '/Users/kje/work/proj/flu/data/jboost'
    segments = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']

    segfiles = dict()
    # load data
    for seg in segments:
        segfile = os.path.join(projpath, seg+'.h7n9.test')
        f = open(segfile, 'r').readlines()
        segfiles[seg] = [fi.strip() for fi in f] 
        
    num_sequences = len(segfiles[segments[0]])

    # concat and write
    outfile = os.path.join(projpath, 'concat.' + '.'.join([seg for seg in segments]) + '.h7n9.test')
    f = open(outfile, 'w')
    for i in range(num_sequences):
        f.write(','.join(segfiles[segments[0]][i].split(',')[:-1]))
        f.write(',*,')
        for s in range(1,len(segments)-1):
            tmp = segfiles[segments[s]][i].split(',')[1:-1]
            f.write(','.join(tmp))
            f.write(',*,')
        f.write(','.join(segfiles[segments[-1]][i].split(',')[1:]))
        f.write('\n')


