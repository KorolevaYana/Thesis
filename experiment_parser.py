#!/usr/bin/python3

from collections import defaultdict

path = '../'
experiment_file = 'GSM2396700_mESC_merged.ghits.pkbin.net.txt'

with open(path + 'files/' + experiment_file, 'r') as f:
    data = f.readlines()

dnas = defaultdict(int)
rnas = defaultdict(int)

rnas_file = open(path + 'tmp/rnas.bed', 'w')
dnas_file = open(path + 'tmp/dnas.bed', 'w')
graph_file = open(path + 'tmp/experiment.graph', 'w')

cur_rna = 0
cur_dna = 0

for s in data:
    i = s.strip().split('\t');
    # rna_chr rna_start rna_end gene_id rna_reads rna_strand dna_chr dna_start val
    rna_chr = i[0]
    rna_start = i[1]
    rna_end = i[2]
    rna_id = i[3]

    dna_chr = i[6]
    dna_start = i[7]

    val = i[8]

    new_rna = (rna_chr, rna_start, rna_end, rna_id)
    new_dna = (dna_chr, dna_start, int(dna_start) + 999)

    if rnas[new_rna] == 0:
        cur_rna += 1
        rnas[new_rna] = cur_rna
        print("%s %s %s %s" % new_rna, file = rnas_file)
    if dnas[new_dna] == 0:
        cur_dna += 1
        dnas[new_dna] = cur_dna
        print("%s %s %s" % new_dna, file = dnas_file)

    print("%d %d %s" % (rnas[new_rna], dnas[new_dna], val), file = graph_file)

rnas_file.close()
dnas_file.close()
graph_file.close()
