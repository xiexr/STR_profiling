dt={}
with open('Lastresult.txt','r') as r:
    for line in r:
        fline = line.strip().split(',')
        if line.startswith('\n'):
            continue
        elif line.startswith('lobSTR') or line.startswith('hipSTR') or line.startswith('gatkvg'):
            title=fline[0]
        elif line.startswith('period'):
            periodlen=len(fline[1])
            allelenum=len(set(fline[-1].split('-')))
            if allelenum>7:
                print(title)
            dt.setdefault(title,[]).extend([periodlen,allelenum])
with open('res_of_each_allele.txt','w') as w:
    w.write('STRname\tperiod\tallelenum\n')
    for k,v in dt.items():
        w.write(f'{k}\t{v[0]}\t{v[1]}\n')

dt={}
with open('res_of_each_allele.txt','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        dt.setdefault(fline[1],{}).setdefault(fline[2],0)
        dt[fline[1]][fline[2]]+=1
with open('each_bp_for_each_allele.txt','w') as w:
    w.write('period\tallelenum\tcounts\n')
    for k,v in dt.items():
        for k1,v1 in v.items():
            w.write(f'{k}\t{k1}\t{v1}\n')