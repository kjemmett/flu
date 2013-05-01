from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import csv
import pdb

def main(segments, projpath):

    seqs = dict()
    labels = dict()
    lengths = dict()

    # load sequences
    for segment in segments:

        seqfile = os.path.join(projpath, 'data', \
                'ncbi.allsubtypes.' + segment + '.aa.sequences.fasta')
        labelfile = os.path.join(projpath, 'data', \
                'ncbi.allsubtypes.' + segment + '.aa.labels.csv')
        reader = csv.reader(open(labelfile, 'r'))

        seglabels = dict((rows[0],rows[1]) for rows in reader)
        seqs[segment] = dict()
        labels[segment] = dict()

        for seq in list(SeqIO.parse(open(seqfile, 'r'), 'fasta')):
            if segment in ['HA', 'NA', 'NS1']:
                desc = ' '.join(seq.description.split(' ')[1:-2])
            else:
                desc = ' '.join(seq.description.split(' ')[1:-1])
            seqs[segment][desc] = seq
            labels[segment][desc] = seglabels[seq.id]
            lengths[segment] = len(seq)

    # get intersection (by description)
    desc_intersection = reduce(lambda x, y: x & y, \
            [set(seqs[seg].keys()) for seg in segments])

    # host dict
    host = {'1':'mammal', '2':'avian'}

    # build new SeqRecord list
    intseqs = []
    for desc in desc_intersection:
        seqname = '_'.join([seqs[seg][desc].name for seg in segments])
        seqid = '_'.join([seqs[seg][desc].name for seg in segments])
        seqdesc = desc
        seqhost = host[str(labels[segments[0]][desc])]
        seqseq = Seq('*'.join([str(seqs[seg][desc].seq) for seg in segments]))
        intseqs.append(SeqRecord(seqseq, id=seqid, name=seqname, \
                description=seqdesc))
        intseqs[-1].host = seqhost

    seqoutfile = os.path.join(projpath, 'data', 'ncbi.allsubtypes.' + '.'.join(segments) + '.concat.sequences.fasta')
    SeqIO.write(intseqs, seqoutfile, 'fasta')

    # write concat object
    # step 2 - build jboost data files

    prefix = 'concat.' + '.'.join([seg for seg in segments])
    f = open(os.path.join(projpath, 'data', 'jboost', \
            prefix + '.full'), 'w')

    i=1
    for seq in intseqs:
        
        # numerical index
        f.write(str(i) + ',')

        # sequence as comma delimited features
        f.write(','.join([aa for aa in seq.seq]))

        # label
        f.write(',' + seq.host + ';\n')
        i = i+1
        
    f.close() 

    # write spec file
    specfile = open(os.path.join(projpath, 'data', 'jboost', \
            prefix+'.spec'), 'w')
    specfile.write('exampleTerminator=;\n')
    specfile.write('attributeTerminator=,\n')
    specfile.write('maxBadExa=20\n')
    specfile.write('INDEX\tnumber\n')
    for idx, seg in enumerate(segments):
        for i in range(1, lengths[seg] + 1):
            specfile.write('%s_x%03d\t(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,-)\n' % (seg, i))
        if idx < len(segments)-1:
            specfile.write('x\t(*)\n')
    specfile.write('labels\t(mammal, avian)')
    specfile.close()


if __name__=="__main__":

    # segments to concatenate
    projpath = '/Users/kje/work/proj/flu'
    segments = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M1', 'NS1']

    # step 1 : load and concat sequences
    main(segments, projpath)


