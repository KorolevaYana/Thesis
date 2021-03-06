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

names = load_rnas('GSM2396700_mESCRNAs.bed')

def load(name, is_norm):
    with open(name, 'r') as f:
        data = f.readlines ()
    result = []
    for i in range(0, len(data)):
        if (is_norm):
            tmp = data[i].strip().split(' ')
            result.append((names[int(tmp[0]) - 1], tmp[1], float(tmp[2]), float(tmp[3])))
        else:
            tmp = data[i].strip().split(' ',3)
            result.append((names[int(tmp[0]) - 1], tmp[1], int(tmp[2])))
    return result

res = open('.work2/sorted.stats', 'w')
norm_res = open('.work2/norm_sorted.stats', 'w')

rna = load('.work2/graph2.stats', False)
rna = sorted(rna, key = lambda x: x[2], reverse = True)

norm_rna = load('.work2/graph2_norm.stats', True)
norm_rna = sorted(norm_rna, key = lambda x: x[2], reverse = True)

for meow in rna:
    print("%s %s %d" % (meow[0], meow[1], meow[2]), file = res)

for meow in norm_rna:
    print("%s %s %f %f" % (meow[0], meow[1], meow[2], meow[3]), file = norm_res)
