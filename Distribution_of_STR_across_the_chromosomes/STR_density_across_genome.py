import re
p=re.compile(r'\d+')
dt={}
with open('Lastresult.txt','r') as r:
    for line in r:
        fline=line.strip().split(',')
        if line.startswith('\n'):
            continue
        elif line.startswith('vg') or line.startswith('STR') or line.startswith('lob'):
            chrom=str(int(p.search(fline[0]).group()[:2]))
            STRname=fline[0]
            pos=fline[2]
            dt.setdefault(chrom,{}).setdefault(STRname,pos)
with open('res.txt','w') as w:
    for k,v in dt.items():
        for k1,v1 in v.items():
            w.write(f'{k},{k1},{v1}\n')