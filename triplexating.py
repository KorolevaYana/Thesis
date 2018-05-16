#!/usr/bin/python3

from collections import defaultdict
import subprocess


def load(name):
    with open(name, 'r') as f:
        data = f.readlines()
    result = []
    rev = {}
    for i in range(0, len(data), 2):
        rev[data[i][1:].strip()] = len(result)
        result.append((data[i][1:].strip(), data[i + 1].strip()))
    return result, rev


rna, rna_rev = load('.work/rnas.fa')
dna, dna_rev = load('.work/dnas.fa')


result = open('.work2/graph2.result', 'w')
result_stats = open('.work2/graph2.stats', 'w')

graph = defaultdict(list)
total = 0
with open('GSM2396700_mESCgraph.bed', 'r') as f:
    for line in f.readlines():
        r, d, value = line.strip().split (' ', 2)
        r = int(r) - 1
        d = int(d) - 1
        graph[r].append(d)
        total += 1

count = 0
for r in graph:
    with open('.work2/rna.temp', 'w') as f:
        print(">%s" % rna[r][0], file = f)
        print(rna[r][1], file = f)
    with open('.work2/dna.temp', 'w') as f:
        for d in graph[r]:
            print(">%s" % dna[d][0], file = f)
            print("%s" % dna[d][1], file = f)
            count += 1
    subprocess.run(["./triplexator-1.3.2-Linux/bin/triplexator", "-rm", "2", "-p", "4", "-o", ".work2/res.temp", "-ss", ".work2/rna.temp", "-ds", ".work2/dna.temp"])
    with open(".work2/res.temp") as f:
        data = f.readlines()
    for line in data[1:]:  # skip headers
        tokens = line.split()
        # chr1:4797869-4876851    3236    3252    chr1:90191000-90191999  167     183     16      0               Y       +       P       0.5
        r = rna_rev[tokens[0]]
        d = dna_rev[tokens[3]]
        print("%d %d %s" % (r + 1, d + 1, line.strip()), file = result)
    print("%d %s %d" % (r + 1, rna[r][0], len(data) - 1), file = result_stats)
    print("%d/%d" % (count, total))

result.close()
result_stats.close()

