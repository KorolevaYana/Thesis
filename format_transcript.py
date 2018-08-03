#!/usr/bin/python3

from collections import defaultdict
import sys
import subprocess

files_path = "../files/"
new_transcript = "transcript.fa"

transcript = sys.argv[1]
rnas = defaultdict(list)
with open(files_path + transcript) as f:
    data = f.readlines()

for i in range(0, len(data), 2):
    rna_data = data[i].strip().split('|')
    # Ensembl_id, gene_name, length, biotype, sequence
    ens_id = rna_data[1].split('.')[0]
    rnas[ens_id].append((ens_id, rna_data[5], rna_data[6], rna_data[7], data[i + 1].strip()));

new_rnas = defaultdict(tuple)
for ids, rna in rnas.items():
    freq = defaultdict(int)
    for i in rna:
        freq[i[3]] += 1
    top_biotype = sorted(freq.items(), key = lambda x: -x[1])[0][0]
    
    tmp = []
    for i in rna:
        if i[3] == top_biotype:
            tmp.append(i)
    new_rna = sorted(tmp, key = lambda x: int(x[2]), reverse = True)[0]
    new_rnas[new_rna[0]] = (new_rna[0], new_rna[1], new_rna[3], new_rna[4])

output = open(files_path + new_transcript, 'w')
for ids, rna in new_rnas.items():
    print("%s|%s|%s\n%s" % rna, file = output)
