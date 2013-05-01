from Bio import SeqIO
import csv
import os

proj_path='/Users/kje/work/proj/flu/data'
segment='NS1'

seq_list = list(SeqIO.parse(open(os.path.join(proj_path,'ncbi.allsubtypes.'+segment+'.aa.sequences.fasta'), 'r'), 'fasta'))
reader = csv.reader(open(os.path.join(proj_path,'ncbi.allsubtypes.'+segment+'.labels.csv'), 'r'))
seq_labels = dict((rows[0],rows[1]) for rows in reader)
f = open(os.path.join(proj_path,'jboost/'+segment+'.full'), 'w')
host = {'1':'mammal', '2':'avian'}
i=1
for seq in seq_list:
    
    # numerical index
    f.write(str(i) + ',')
    # sequence as comma delimited features
    f.write(','.join([aa for aa in seq.seq]))
    # label
    
    f.write(',' + host[seq_labels[seq.id]] + ';\n')
    #f.write(',' + host['1'] + ';\n')
    i = i+1
    
f.close() 

# write spec file
specfile = open(os.path.join(proj_path,'jboost/'+segment+'.spec'), 'w')
specfile.write('exampleTerminator=;\n')
specfile.write('attributeTerminator=,\n')
specfile.write('maxBadExa=20\n')
specfile.write('INDEX\tnumber\n')
specfile.write('\n'.join(['pos%03d\t(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,-)' % i for i in range(1, len(seq_list[0])+1)]))
specfile.write('\nlabels\t(mammal, avian)\n')
specfile.close()
