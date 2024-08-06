from pathlib import Path
import re
u=re.compile(r'\d+')
dt={}
with open('Last_GWAS_Rst_res.txt','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        dt.setdefault(fline[0],fline[4])
p=Path()
totaldt={}
for files in p.glob('cal_*.txt'):
    title=files.stem.split('_')[-2:]
    title='_'.join(title)
    with open(files,'r') as r:
        for line in r:
            fline=line.strip().split('\t')
            STRname='STR'+u.search(fline[0]).group(0)
            if dt.get(STRname):
                if float(fline[-1])>=0.25:
                    totaldt.setdefault(title,{}).setdefault(dt.get(STRname),0)
                    totaldt[title][dt.get(STRname)]+=1
with open('GWAS_pop_sangkitu.txt','w') as w:
    w.write('pop1_pop2\tGWAS_trait\tvalue\n')
    for k,v in totaldt.items():
        for k1,v1 in v.items():
            w.write(f'{k}\t{k1}\t{v1}\n')
templs=[]
for k, v in totaldt.items():
    templs.append(k)
for k, v in totaldt.items():
    for k1,v1 in v.items():
        if k1 not in templs:
            templs.append(k1)
with open('nodes_sangki.txt','w') as w:
    total='\n'.join(templs)
    w.write(f'{total}')

