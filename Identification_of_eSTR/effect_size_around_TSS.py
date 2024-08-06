import pandas as pd
from collections import OrderedDict
SDAT6_PATH=  'eSTRwithFDRxiaoyu005.txt'
realstrs = pd.read_csv(SDAT6_PATH, sep = '\t')
dt={}
for i in range(-100000,110000,5000):
    dt.setdefault((i,i+5000),{})
    dt.setdefault((i,i+5000),{})
tempdt=realstrs.to_dict(orient='index')
for k,v in tempdt.items():
    lengthofSTR = len(v.get('STR site').split('-')[-1])
    for k1 in dt:
        if v.get('feature')=='upstream':
            if k1[0]<=(-int(v.get('distance').split('_')[-1]))<=k1[1]:
                dt.setdefault(k1, {}).setdefault(lengthofSTR,[]).append(abs(v.get('effect size')))
                dt.setdefault(k1,{}).setdefault('totalSTR',[]).append(abs(v.get('effect size')))
                break
        elif v.get('feature') == 'downstream':
            if k1[0] <= int(v.get('distance').split('_')[-1])+10000 <= k1[1]:
                dt.setdefault(k1, {}).setdefault(lengthofSTR,[]).append(abs(v.get('effect size')))
                dt.setdefault(k1, {}).setdefault('totalSTR', []).append(abs(v.get('effect size')))
                break
        else:
            if k1[0] <= round(v.get('relative position in gene')*10000) <= k1[1]:
                dt.setdefault(k1, {}).setdefault(lengthofSTR, []).append(abs(v.get('effect size')))
                dt.setdefault(k1, {}).setdefault('totalSTR', []).append(abs(v.get('effect size')))
                break

################################
size=[]
meaneffsize={}
for i in range(1,7):
    meaneffsize.setdefault(i,[])
meaneffsize.setdefault('totalSTR',[])
for k,v in dt.items():
    size.append(k[1])
    for k1,v1 in v.items():
        meaneff=sum(v1)/len(v1)
        meaneffsize.setdefault(k1,[]).append(meaneff)
lastdt={}
totaldt={}
for i in range(len(size)):
    totaldt.setdefault('totalSTR',{}).setdefault(size[i],meaneffsize.get('totalSTR')[i])
for u in range(1, 7):
    for i in range(len(size)):
        lastdt.setdefault(f'{u}bp',{}).setdefault(size[i],meaneffsize.get(u)[i])
df=pd.DataFrame.from_dict(lastdt)
df1=pd.DataFrame.from_dict(totaldt)
df.to_csv('output.csv',sep='\t')
df1.to_csv('totalSTRoutput.csv',sep='\t')