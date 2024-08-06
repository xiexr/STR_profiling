from collections import OrderedDict
import pandas as pd
dt=OrderedDict()
import sys
import re
fairedt={}
memedt={}
p=re.compile(r'\d+')

with open('overlappedwithfaire.txt','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        number=p.search(fline[0]).group(0)
        newname='STR'+number
        fairedt.setdefault(newname,{})
with open('overlappedwithmemeres.txt','r') as r:#STR名字与meme相关
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        number=p.search(fline[0]).group(0)
        newname='STR'+number
        memedt.setdefault(newname,{})
with open('STRdensityoutput10k.txt','r') as r:
    with open('overlappedwithfairedensityoutput_10k.txt','w') as w1:#跟faire重叠的STR的位置信息
        with open('overlappedwithmemedensityoutput_10k.txt','w') as w2:#跟meme重叠的STR的位置信息
            for line in r:
                fline=line.strip().split('\t')
                if fline[0] in fairedt:
                    w1.write(line)
                if fline[0] in memedt:
                    w2.write(line)

lst=[]
xx,file3=sys.argv
for j in range(-10000,10001,100):
    dt.setdefault((j,j+100),[])
with open(file3,'r') as r:
    for line in r:
        fline=line.strip().split('\t')
        number=fline[-2].split(':')[-1]
        if number.startswith('-') or number.isdigit():
            number = int(number)
            if -10000<number<10001:
                lst.append(number)
lst.sort()
with open(f'{file3[:-4]}_not_uniform.txt','w') as w:
    for each in lst:
        w.write(f'{file3[14:-18]}\t{each}\n')