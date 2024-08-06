dt={}
with open('read.txt','r') as r:
    for line in r:
        fline=line.strip().split('\t')
        dt.setdefault(fline[0],[]).append(int(fline[1]))
for k,v in dt.items():
    print(k,sum(v)/len(v))