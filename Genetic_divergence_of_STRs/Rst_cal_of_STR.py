import sys
from math import log2
from statistics import variance
import re
p=re.compile(r'\d+')
xx,popu1,popu2=sys.argv
dt={}
hetdt={}
with open('Each_population_Het.txt','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        if float(fline[1]) >0.1:
            hetdt.setdefault(fline[0],'')
with open('Lastresult.txt','r') as r:
    for line in r:
        if line.startswith('\n'):
            continue
        elif line.startswith('lobSTR') or line.startswith('hipSTR') or line.startswith('gatkvg'):
            title=line.strip().split(',')[0]
        elif line.startswith('period'):
            fline=line.strip().split(',')
            cpnum=fline[-1].split('-')
            cpnum_of_pop1=[]
            cpnum_of_pop2=[]
            temptotalnum1=0
            temptotalnum2=0
            i=0
        elif line[0].isdigit():
            if title in hetdt:
                fline = line.strip().split(',')
                if fline[2]==popu1:
                    for index,each in enumerate(fline[0].split('-')[::2]):
                        cpnum_of_pop1.extend([int(cpnum[index])*2]*int(each))
                    i+=1
                elif fline[2]==popu2:
                    for index,each in enumerate(fline[0].split('-')[::2]):
                        cpnum_of_pop2.extend([int(cpnum[index])*2]*int(each))
                    i += 1

                if i==2:
                    if len(cpnum_of_pop1) < 2 or len(cpnum_of_pop2) < 2:
                        dt.setdefault(title, []).extend(['0', str(len(cpnum_of_pop1)), str(len(cpnum_of_pop2))])
                        i=0
                        continue
                    S_w = (variance(cpnum_of_pop1) + variance(cpnum_of_pop2)) / 2.0
                    S_bar = variance(cpnum_of_pop1 + cpnum_of_pop2)
                    if S_bar != 0:
                        Rst = (S_bar - S_w) / S_bar
                    else:
                        Rst = 0
                    if Rst < 0:
                        Rst = 0
                    dt.setdefault(title, []).extend([ str(Rst), str(len(cpnum_of_pop1)), str(len(cpnum_of_pop2))])
                    i = 0
with open(f'cal_Rst_res_{popu1}_{popu2}.txt','w') as w:
    for k,v in dt.items():
        w.write(f'{k}\t{v[-2]}\t{v[-1]}\t{v[0]}\n')
