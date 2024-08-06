import pickle
from collections import OrderedDict
import sys
xx,file=sys.argv
dt=OrderedDict()
with open(file,'r') as r:
    for line in r:
        if line.startswith('\n'):
            continue
        elif line.startswith('gatkvg') or line.startswith('lobSTR') or line.startswith('hipSTR'):
            key=line.strip().split(',')[0]
        elif line[0].isdigit():
            fline = line.strip().split(',')
            templs = []
            freq = []
            templs = fline[0].split('-')[::2]
            templs = list(map(int, templs))
            total = sum(templs)
            freq = [(i / total) ** 2 for i in templs]
            het = round(1 - sum(freq), 4)
            Nfreq = round(int(fline[-3]) / int(fline[1]), 4)
            dt.setdefault(key,[]).extend([str(het),str(Nfreq)])
with open('STRhetall.txt','w') as w:
    w.write(
        'Name,All,AllNfreq,AllIndica,AllIndicaNfreq,AllJaponica,AllJaponicaNfreq,IndicaI,IndicaINfreq,IndicaII,IndicaIINfreq,IndicaIII,IndicaIIINfreq,IndicaIntermediate,IndicaIntermediateNfreq,Aus,AusNfreq,TemperateJaponica,TemperateJaponicaNfreq,TropicalJaponica,TropicalJaponicaNfreq,JaponicaIntermediate,JaponicaIntermediateNfreq,Aromatic,AromaticNfreq,Intermediate,IntermediateNfreq\n')
    for k,v in dt.items():
        w.write(f'{key},{",".join(v)}\n')


templs=["All","All_Indica","All_Japonica","Aus","Aromatic","Intermediate"]
dt=OrderedDict()
with open('STRhetall.txt', 'r') as r:
    for line in r:
        if line.startswith('Name') or 'unknown' in line:
            continue
        fline = line.strip().split(',')
        thresholds = [2, 4, 6, 16, 24, 26]
        for idx, threshold in enumerate(thresholds):
            if float(fline[threshold]) < 0.6:
                dt.setdefault(templs[idx], []).append((fline[0], fline[threshold - 1]))

for k,v in dt.items():
    templs = []
    for each in v:
        templs.append(each)
    templs.sort(key=lambda x:x[1],reverse=True)
    if len(templs)%4:
        num=(len(templs)//4)+1
    else:
        num=(len(templs)//4)
    dt[k]=templs[:num]

with open(f'STRhetall_filter06_leave6_qian_25.txt','w') as w:
    w.write('Het\tSpecies\n')
    for k,v in dt.items():
        for each in v:
            w.write(f'{each[1]}\t{k}\n')