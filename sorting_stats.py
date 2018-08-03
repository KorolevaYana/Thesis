#!/usr/bin/python3

from collections import defaultdict
import subprocess
import operator

def load(name):
    with open(name, 'r') as f:
        data = f.readlines ()
    result = []
    for i in range(0, len(data)):
        tmp = data[i].strip().replace("| ", "").split(' ')
        result.append((tmp[1], float(tmp[2]), float(tmp[3]), 
                               float(tmp[4]), float(tmp[5]),
                               float(tmp[6]), float(tmp[7])))
    return result

lncRNA_res = 'lncRNAs.stats'
ncRNA_res = 'ncRNAs.stats'
pr_c_res = 'prot_c.stats'
others_res = 'others.stats'

medians_res = open('medians.stats', 'w')
norms = ['grid_norm3', 'norm3', 'grid_norm2', 'norm2', 'grid_norm1', 'norm1']

def print_stats(name):
    res = load(name)
    if len(res) > 0:
        print("TOP-20:", file = medians_res)
        res = sorted(res, key = operator.itemgetter(1, 3, 2, 4, 5, 6), reverse = True)
        for i in range(0, min(20, len(res))):
            print("%s %.10e %.10e | %.10e %.10e | %.10e %.10e" %
                    (res[i][0], res[i][1], res[i][2], res[i][3], res[i][4], res[i][5], res[i][6]),
                    file = medians_res)
        print('\n', file = medians_res)
        for i in range(1, 7):
            print(norms[i - 1], file = medians_res)
            res = sorted(res, key = lambda x: x[i], reverse = True)
            print("%.10e" % res[int(len(res)/2)][i], file = medians_res)
            print("%.10e" % res[int(len(res)/4)][i], file = medians_res)
        print('\n', file = medians_res)

print('lncRNA:', file = medians_res)
print_stats(lncRNA_res)

print('ncRNA:', file = medians_res)
print_stats(ncRNA_res)

print('pr_c:', file = medians_res)
print_stats(pr_c_res)

print('others:', file = medians_res)
print_stats(others_res)

medians_res.close()
