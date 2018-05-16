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

def load(name):
    with open(name, 'r') as f:
        data = f.readlines ()
    result = []
    for i in range(0, len(data)):
        tmp = data[i].strip().split(' ',3)
        result.append((names[int(tmp[0]) - 1], tmp[1], int(tmp[2])))
    return result

res = open('.work2/sorted.stats', 'w')

rna = load('.work2/graph2.stats')
rna = sorted(rna, key = lambda x: x[2], reverse = True)

for meow in rna:
    print("%s %s %d" % (meow[0], meow[1], meow[2]), file = res)
