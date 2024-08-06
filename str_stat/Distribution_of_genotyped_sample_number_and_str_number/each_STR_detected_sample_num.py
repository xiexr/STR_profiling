STRnamewithNumber={}
with open('Lastresult.txt','r') as r:
    for line in r:
        if line.startswith('\n'):
            continue
        fline=line.strip().split(',')
        if line.startswith('lobSTR') or line.startswith('STR') or line.startswith('vg'):
            title=fline[0]
        if not title.startswith('vg'):
            if len(fline)>2:
                if fline[2]=='All':
                    realSTRnumber=int(fline[1])-int(fline[-1])-int(fline[-3])
                    STRnamewithNumber.setdefault(title,realSTRnumber)
with open('totalstatisticofeachSTRfrequencywithoutvg.txt','w') as w:
    for k,v in STRnamewithNumber.items():
        w.write(f'{k}\t{v}\n')

import sys
from collections import Counter
xx,vgfiles,vgoriginal=sys.argv
vgdt={}
with open(vgfiles,'r') as r:
    for line in r:
        if line.startswith('\n'):
            continue
        elif line.startswith('vg'):
            fline=line.strip().split(',')
            vgdt.setdefault(fline[0],'')
indexcorrespondsample={}
rice3k={}
with open('samplename.txt','r') as r:
    for line in r:
        fline=line.strip().split(',')
        rice3k.setdefault(fline[0],'')
totalvg={}
with open(vgoriginal,'r') as r:
    samples=r.readline().strip().split(',')[1:]
    indexcorrespondsample=dict(enumerate(samples))
    tempdt=indexcorrespondsample.copy()
    for k,v in tempdt.items():
        if not v in rice3k:
            del indexcorrespondsample[k]
    print(len(indexcorrespondsample))
    for line in r:
        fline=line.strip().split(',')
        if fline[0] in vgdt:
            allsamples=fline[1:]
            for k in indexcorrespondsample:
                totalvg.setdefault(fline[0],[]).append(allsamples[int(k)])
with open(f'statisticofeachSTRfrequencychr{vgfiles[17:-26]}.txt','w') as w:
    for k,v in totalvg.items():
        statistic=Counter(v)
        NandH=0
        for k1,v1 in statistic.items():
            if k1 == 'N':
                NandH += v1
            elif k1=='H':
                NandH += v1
        restfind=3046-NandH
        w.write(f'{k}\t{restfind}\n')


