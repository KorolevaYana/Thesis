#!/usr/bin/python3

from collections import defaultdict
import subprocess


ncRNA = {'Mt_rRNA',
         'Mt_tRNA',
         'miRNA', 
         'misc_RNA',
         'rRNA',
         'scRNA',
         'snRNA',
         'snoRNA',
         'ribozyme',
         'sRNA',
         'scaRNA'}
lncRNA = {'lincRNA'}
pr_coding = {'protein_coding',
             'IG_X',
             'TR_X'}


def load(name):
    with open(name, 'r') as f:
        data = f.readlines()
    result = []
    rev = {}
    for i in range(0, len(data), 2):
        rev[data[i][1:].strip()] = len(result)
        result.append((data[i][1:].strip(), data[i + 1].strip()))
    return result, rev


rna, rna_rev = load('spliced_rnas.fa')
dna, dna_rev = load('.work/dnas.fa')


res_ncRNA = open('ncRNAs.stats', 'w')
res_lncRNA = open('lncRNAs.stats', 'w')
res_prot_c = open('prot_c.stats', 'w')
res_others = open('others.stats', 'w')
res_none = open('none.stats', 'w')

graph = defaultdict(list)
values = defaultdict(float)
g_num = defaultdict(int)
total = 0
with open('.work/GSM2396700_mESCgraph.bed', 'r') as f:
    for line in f.readlines():
        r, d, value, g = line.strip().split(' ')
        r = int(r) - 1
        d = int(d) - 1
        graph[r].append(d)
        values[(r,d)] = float(value)
        g_num[(r,d)] = int(g)
        total += 1

count = 0
for r in graph:
    if rna[r][0] != "":
        tmp = rna[r][0].strip().split('|');
        rna_name = tmp[1]
        rna_type = tmp[2]
        if rna_type in ncRNA:
            result = res_ncRNA
        elif rna_type in lncRNA:
            result = res_lncRNA
        elif rna_type in pr_coding:
            result = res_prot_c
        else:
            result = res_others
        with open('.work3/rna.temp', 'w') as f:
            print(">%s" % rna[r][0], file = f)
            print(rna[r][1], file = f)
        with open('.work3/dna.temp', 'w') as f:
            for d in graph[r]:
                print(">%s" % dna[d][0], file = f)
                print("%s" % dna[d][1], file = f)
                count += 1
        subprocess.run(["./triplexator-1.3.2-Linux/bin/triplexator", "-rm", "2", "-p", "4", "-l", "12", "-o", ".work3/res.temp", "-ss", ".work3/rna.temp", "-ds", ".work3/dna.temp"])
        with open(".work3/res.temp.summary") as f:
            data = f.readlines()
        triplexes = 0
        norm_res1 = 1
        norm_res2 = 1
        norm_res3 = 1
        grid_norm_res1 = 1
        grid_norm_res2 = 1
        grid_norm_res3 = 1
        for line in data[1:]:  # skip headers
            tokens = line.split()
            # chr1:158318000-158318999  chrX:159198334-159217024  22  4.02e-08  0 0 0 0 22  4.02e-08
            d = dna_rev[tokens[0]]
            r = rna_rev[tokens[1]]
            gr_num = g_num[(r,d)]
            if gr_num == 1:
                grid_norm_res1 *= 1 - values[(r, d)] * float(tokens[3])
                norm_res1 *= 1 - float(tokens[3])
            elif gr_num == 2:
                grid_norm_res2 *= 1 - values[(r, d)] * float(tokens[3])
                norm_res2 *= 1 - float(tokens[3])
            if gr_num == 3:
                grid_norm_res3 *= 1 - values[(r, d)] * float(tokens[3])
                norm_res3 *= 1 - float(tokens[3])
            triplexes += int(tokens[2])
        print("%d %s %.10e %.10e | %.10e %.10e | %.10e %.10e" % (r + 1, rna_name, 
            1 - grid_norm_res3, 1 - norm_res3, 
            1 - grid_norm_res2, 1 - norm_res2, 
            1 - grid_norm_res1, 1 - norm_res1), file = result)
        print("%d/%d" % (count, total))
    else:
        print(r, file = res_none)

res_ncRNA.close()
res_lncRNA.close()
res_prot_c.close()
res_others.close()
