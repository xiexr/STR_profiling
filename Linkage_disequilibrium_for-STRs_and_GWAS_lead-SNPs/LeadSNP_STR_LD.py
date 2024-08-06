import pandas as pd
from multiprocessing import Pool
import numpy as np
from scipy import stats
from scipy.stats import pearsonr
import os
import sys
import numpy as np
import re
###read_snp_STR_pair
snp_str_pair={}
STRs = {}
def cal_snp_str_ld(df_snp = None, df_str = None, snp=None, STR=None):

    x = df_str[STR]


    snp_chrom, snp_pos = snp.split(':')
    if snp in df_snp.columns:
        y = df_snp[snp]
        ld = cal_snp_str_ld_helper(x, y)
    else:
        ld = 'nan'
    return [snp, STR, ld]

def cal_snp_str_ld_helper(x=None, y=None):
    if isinstance(y, pd.Series):
        x_numeric = pd.to_numeric(x, errors='coerce')
        y_numeric = pd.to_numeric(y, errors='coerce')
        contingency = pd.crosstab(x_numeric, y_numeric)
        with open('test.txt','a') as w:
            print(contingency,file=w)
        cor = x_numeric.corr(y_numeric)
        ld = cor * cor
    else:
        ld = 'nan'
    return str(ld)
with open('SNP_STR_250kres.txt','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        snp = fline[0] + ":" + fline[2]
        STR = fline[3] + ":" + fline[-1]
        snp_str_pair.setdefault(snp,[]).append(STR)
        STRs.setdefault(STR,0)
############### Load SNP dosage
df_snp = pd.read_csv('res_of_SNP_trans_num.txt', sep=',')
df_snp.rename(columns={"Unnamed: 0": "ID"}, inplace=True)
df_snp = df_snp.set_index('ID')
df_snp=df_snp.T
# print(df_snp)
# Load STR dosage
df_str = pd.read_csv('STR_res.txt', sep=',',low_memory=False)
df_str.rename(columns={"Unnamed: 0": "ID"}, inplace=True)
df_str = df_str.set_index('ID')
df_str = df_str.T
x = df_str.columns.to_list()
x = {i:f'{str(int(i[3:5]))}:{i}' for i in x}
df_str.rename(columns=x, inplace=True)
df_str = df_str.reindex(index=df_snp.index)
x = df_snp.columns.to_list()

x = {i:f'{str(int(i[2:4]))}:{i}' for i in x}
df_snp.rename(columns=x, inplace=True)
df_snp.replace("N", np.nan, inplace=True)
df_snp.replace("H", np.nan, inplace=True)
df_str.replace("N", np.nan, inplace=True)
df_str.replace("H", np.nan, inplace=True)

results=[]
with open('SNP_STR_LD_value.txt','w') as w:
    for snp in snp_str_pair:
        for STR in snp_str_pair[snp]:
            res = cal_snp_str_ld(df_snp, df_str, snp, STR)
            results.append(res)

    for res in results:
        snp, STR, ld = res
        snp_chrom, snp_pos = snp.split(':')
        STR_chrom, STR_pos = STR.split(':')
        w.write('\t'.join([snp_chrom, snp_pos, STR_chrom, STR_pos, ld]) + '\n')




p=re.compile(r'\d+')
dt={}
for i in range(0,10000,1000):
    dt.setdefault((i,i+1000),[])
for i in range(10000,250000,10000):
    dt.setdefault((i,i+10000),[])
print(dt)
D={}
with open("SNP_STR_LD_value.txt", 'r') as f:
   for line in f:
      fline=line.strip().split("\t")
      if fline[-1] != "nan":
        snp_pos=p.search(fline[1]).group()[2:]
        str_pos=p.search(fline[3]).group()[2:]
        r2=float(fline[-1])
        dist = abs(int(snp_pos)-int(str_pos))
        for k,v in dt.items():
         if k[0]<=dist<k[1]:
             dt.setdefault(k,[]).append(r2)
             break
with open("output_for_SNP_STR_LD_paint.txt", 'w') as w:
    w.write(f"dist\taverLD\n")
    for k,v in dt.items():
        averageofLD=sum(v)/len(v)
        w.write(f'{k[1]}\t{averageofLD}\n')

