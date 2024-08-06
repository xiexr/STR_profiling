import sys
import re
from pathlib import Path
xx,totaltable,lobstr,hipstr,gatk=sys.argv

lobdt={}
hipdt={}
vgdt={}
with open(totaltable,'r') as r:
    for line in r:
        if line.startswith('\n'):
            continue
        fline = line.strip().split(',')
        if line.startswith('lobSTR'):
            lobdt.setdefault(str(int(fline[0][8:])),'')
        elif line.startswith('STR'):
            hipdt.setdefault(str(int(fline[0][5:])),'')
        elif line.startswith('vg'):
            vgdt.setdefault(fline[0],'')
# print(vgdt)
samplewithindex={}
totaltime={}
with open(lobstr,'r') as r:
    for line in r:
        if line.startswith('##'):
            continue
        samplefline=line.strip().split('\t')[9:]
        if line.startswith('#CHROM'):
            for index,sample in enumerate(samplefline):
                samplewithindex.setdefault(index,sample)
                totaltime.setdefault(sample, 0)
            # print(len(totaltime))
        else:
            pos=line.strip().split('\t')[1]
            if pos in lobdt:
                for index,sample in enumerate(samplefline):
                    if re.search(r'(\d)/\1',sample.split(':')[0]):
                        samplename=samplewithindex.get(index)
                        totaltime[samplename] += 1
samplewithindex={}
with open(hipstr,'r') as r:
    for line in r:
        if line.startswith('##'):
            continue
        samplefline=line.strip().split('\t')[9:]
        if line.startswith('#CHROM'):
            for index,sample in enumerate(samplefline):
                samplewithindex.setdefault(index,sample)
                totaltime.setdefault(sample, 0)
        else:
            pos=line.strip().split('\t')[1]
            if pos in hipdt:
                for index,sample in enumerate(samplefline):
                    if re.search(r'(\d)\|\1',sample.split(':')[0]):
                        samplename=samplewithindex.get(index)
                        totaltime[samplename] += 1
samplewithindex={}
with open(gatk,'r') as r:
    title=r.readline()
    for index, sample in enumerate(title.strip().split(',')[1:]):
        samplewithindex.setdefault(index, sample)
        totaltime.setdefault(sample, 0)
    for line in r:
        samplefline=line.strip().split(',')
        pos=samplefline[0]
        # print(pos)
        if pos in vgdt:
            for index, sample in enumerate(samplefline[1:]):
                # print(sample)
                if (not 'N' in sample) and (not 'H' in sample):
                    samplename = samplewithindex.get(index)
                    totaltime[samplename] += 1
with open(f'statisticofeachsamlefrequency{totaltable[13:-4]}.txt','w') as w:
    for k,v in totaltime.items():
        w.write(f'{k}\t{v}\n')
#
for file_name in Path('statisticofeachsamlefrequency*'):
    with open('totalstatisticofeachsamplefrequency.txt','w') as w:
        with open(file_name,'r') as r:
            for line in r:
                w.write(line)


rice3k={}
with open('samplename.txt','r') as r:
    for line in r:
        fline=line.strip().split(',')
        rice3k.setdefault(fline[0],'')
each3k={}
with open('totalstatisticofeachsamplefrequency.txt','r') as r:
    with open('totalresultofeachsample.txt','w') as w:
        for line in r:
            fline=line.strip().split('\t')
            if fline[0] in rice3k:
                each3k.setdefault(fline[0],0)
                each3k[fline[0]] += int(fline[1])
        for k,v in each3k.items():
            w.write(f'{k}\t{v}\n')
