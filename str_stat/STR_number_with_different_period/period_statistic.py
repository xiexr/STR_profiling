import sys
from collections import OrderedDict
xx,file=sys.argv
dt=OrderedDict()
with open(file,'r') as r:
    for line in r:
        if line.startswith('period'):
            fline=line.strip().split(',')
            dt.setdefault(len(fline[1]),0)
            dt[len(fline[1])]+=1
newdt=sorted(dt.items(),key=lambda x:x[0],reverse=False)
with open('statistic.txt','w') as w:
    for k,v in newdt:
        w.write(f'{k}\t{v}\n')
