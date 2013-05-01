from Bio import SeqIO
import csv
import os

data_path = '/Users/kje/work/proj/flu/data/'
proj_path = '/Users/kje/work/proj/flu/data/jboost'

#segments = ['PB2', 'PB1', 'PA', 'NP', 'NA', 'MP', 'NS']
segments = ['NS']

for segment in segments:
    seq_list = list(SeqIO.parse(open(os.path.join(data_path,'gisaid.H7N9.'+segment+'.human.test.fa'), 'r'), 'fasta'))
    f = open(os.path.join(proj_path, segment+'.h7n9.test'), 'w')

    i=1
    for seq in seq_list:
        
        # numerical index
        f.write(str(i) + ',')
        # sequence as comma delimited features
        f.write(','.join([aa for aa in seq.seq]))
        # label
        
        f.write(',' + 'mammal' + ';\n')
        i = i+1
        
    f.close() 
