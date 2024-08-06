from collections import OrderedDict
import sys
import pandas as pd
xx,file=sys.argv
dt=OrderedDict()
with open('Lastresult.txt','r') as r1:
    with open(file,'r') as r2:
        for line in r1:
            if line.startswith('STR'):
                fline=line.strip().split(',')
                STRname=fline[0]
            elif line.startswith('period'):
                fline=line.strip().split(',')
                pdlen=len(fline[1])
                dt.setdefault(STRname,{}).setdefault('period',pdlen)
        for line in r2:
            fline=line.strip().split('\t')
            dt.setdefault(fline[0],{}).setdefault('pos',fline[-1])
tempdt=OrderedDict()
for k,v in dt.items():
    tempdt.setdefault(v.get('pos'),{}).setdefault(v.get('period'),0)
    tempdt[v.get('pos')][v.get('period')]+=1
df=pd.DataFrame.from_dict(tempdt)
df.to_csv('STRstatisticwitheachlen.csv',header=True)