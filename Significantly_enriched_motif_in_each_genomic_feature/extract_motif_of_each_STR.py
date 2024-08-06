import re
from collections import Counter
p=re.compile(r'\d+')
perioddt={}
i=0
with open('Lastresult.txt','r') as r:
    for line in r:
        fline=line.strip().split(',')
        if line.startswith('\n'):
            continue
        elif line.startswith('lob') or line.startswith('STR') or line.startswith('vg'):
            title=p.search(fline[0]).group()
        elif line.startswith('period'):
            perioddt.setdefault(title,fline[1])

totaldt={}
countdt={}

with open('STRoutput5kb.txt','r') as r:
    for line in r:
        fline=line.strip().split('\t')
        title=p.search(fline[0]).group()
        period=perioddt.get(title)
        totaldt.setdefault(title,{}).setdefault(period,fline[-1])
        if fline[-1]=='exon':
            fline[-1]='CDS'
        countdt.setdefault(fline[-1],[]).append(period)

with open('result.txt','w') as w:
    w.write(f'title\tmotif_geno\tgfeature\n')
    for k,v in totaldt.items():
        for k1,v1 in v.items():
            w.write(f'{k}\t{k1}\t{v1}\n')
