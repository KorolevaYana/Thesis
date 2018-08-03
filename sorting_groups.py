#!/usr/bin/python3

from collections import defaultdict
import subprocess

def load_rnas(name):
    with open (name, 'r') as f:
        data = f.readlines ()
    result = []
    for i in range(0, len(data)):
        tmp = data[i].strip().split('\t', 6)
        result.append(tmp[3])
    return result

names = load_rnas('.work/GSM2396700_mESCRNAs.bed')

def load(name, is_norm):
    with open(name, 'r') as f:
        data = f.readlines ()
    result = []
    for i in range(0, len(data)):
        if (is_norm):
            tmp = data[i].strip().replace("| ", '').split(' ')
            result.append((names[int(tmp[0]) - 1], tmp[1], 
                float(tmp[2]), float(tmp[3]),
                float(tmp[4]), float(tmp[5]),
                float(tmp[6]), float(tmp[7])))
        else:
            tmp = data[i].strip().split(' ',3)
            result.append((names[int(tmp[0]) - 1], tmp[1], int(tmp[2])))
    return result

res = open('.work3/sorted.stats', 'w')
norm_res = open('.work3/norm_sorted.stats', 'w')

rna = load('.work3/graph.stats', False)
rna = sorted(rna, key = lambda x: x[2], reverse = True)

norm_rna = load('.work3/graph_norm.stats', True)
norm_rna = sorted(norm_rna, key = lambda x: x[2], reverse = True)

print("Gene_id\tchrom_pos\tgrid_norm3\tnorm3\t|\tgrid_norm2\tnorm2\t|\tgrid_norm1\tnorm1", file = norm_res)

for meow in rna:
    print("%s %s %d" % (meow[0], meow[1], meow[2]), file = res)

for meow in norm_rna:
    print("%s %s %f %f | %f %f | %f %f" 
            % (meow[0], meow[1], meow[2], meow[3], meow[4], meow[5], meow[6], meow[7]), file = norm_res)
