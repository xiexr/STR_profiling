from Bio import SeqIO
from collections import OrderedDict
from collections import Counter
import re
import pandas as pd
import numpy as np
from scipy.stats import zscore

p=re.compile(r'\d+')
dt=OrderedDict()
for rec in SeqIO.parse('Nipponbare.fa','fasta'):
    for i in range(0, len(rec.seq), 100000):
        if len(rec.seq) - i < 100000:
            for j in range(1,7):
                dt.setdefault(f'chr{str(rec.id).rjust(2, "0")}', {}).setdefault((i + 1, len(rec.seq)), {}).setdefault(j,0)
                dt.setdefault(f'chr{str(rec.id).rjust(2, "0")}', {}).setdefault((i + 1, len(rec.seq)), {}).setdefault('total',0)
        else:
            for j in range(1, 7):
                dt.setdefault(f'chr{str(rec.id).rjust(2, "0")}', {}).setdefault((i + 1, i + 100000), {}).setdefault(j,0)
                dt.setdefault(f'chr{str(rec.id).rjust(2, "0")}', {}).setdefault((i + 1, i + 100000), {}).setdefault('total', 0)
with open(f'Lastresult','r') as r:
    for line in r:
        if line.startswith('STR') or line.startswith('lob') or line.startswith('vg'):
            fline=line.strip().split(',')
            chromo='chr'+p.search(fline[0]).group()[:2]
            pos=int(p.search(fline[0]).group()[2:])
        elif line.startswith('period'):
            fline=line.strip().split(',')
            periodlen=len(fline[1])
            for k,v in dt.get(chromo).items():
                if k[0]<=pos<k[1]:
                    v[periodlen]+=1
                    v['total']+=1
                    break

with open('result.txt','w') as w:
    w.write('chrom\tstart\tend\tp1\tp2\tp3\tp4\tp5\tp6\tAll\n')
    for k,v in dt.items():
        for k1,v1 in v.items():
            w.write(f'{k}\t{k1[0]}\t{k1[1]}\t{v1.get(1)}\t{v1.get(2)}\t{v1.get(3)}\t{v1.get(4)}\t{v1.get(5)}\t{v1.get(6)}\t{v1.get("total")}\n')

dat = pd.read_csv('result.txt', sep='\t')
cents = pd.read_csv('cento.txt', sep='\t')

def neatPileup(dat, cents, window=11, smooth=True):
    if smooth:
        for k in dat['chrom'].unique():
            dat.loc[dat['chrom'] == k, dat.columns[3:]] = dat.loc[dat['chrom'] == k, dat.columns[3:]].apply(
                lambda vals: pd.Series(vals).rolling(window=window, min_periods=1, center=True).mean())
    for k in dat['chrom'].unique():
        b=dat.loc[dat['chrom'] == k, dat.columns[2]].max()
        dat.loc[dat['chrom'] == k, 'max_chrom'] = b
    # Calculate normalized distance from centromere
    dat['cdist.norm'] = dat.apply(lambda row: calculate_normalized_distance(row, cents), axis=1)

    return dat


def calculate_normalized_distance(row, cents):
    k = row['chrom']
    midpoints = (row['start'] + row['end']) / 2
    max_point=row['max_chrom']
    c_mid = (cents.loc[cents['chromosome'] == k, 'centostart'].min() + cents.loc[cents['chromosome'] == k, 'centoend'].max()) / 2
    c_dists = midpoints - c_mid
    plen = c_mid
    qlen = max_point - c_mid
    c_dists_n = np.where(c_dists < 0, c_dists / plen, c_dists / qlen)
    c_dists_n = np.where(c_dists_n < -1, -1, np.where(c_dists_n > 1, 1, c_dists_n))
    return c_dists_n

dat = neatPileup(dat, cents, smooth=True)


dat.to_csv('STR_densitytocento.dat.txt', sep='\t', index=False)