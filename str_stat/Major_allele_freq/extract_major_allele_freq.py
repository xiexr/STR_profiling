dt={}
with open('Lastresult.txt','r') as r:
    for line in r:
        fline=line.strip().split(',')
        if line.startswith('\n'):
            continue
        elif line.startswith('lobSTR') or line.startswith('STR') or line.startswith('vg'):
            title=fline[0]
        elif line.startswith('period'):
            periodlen=len(fline[1])
        elif fline[0][0].isdigit():
            if fline[2]=='All':
                eachfreq=[float(each.split('%')[0])*0.01 for each in fline[0].split('-')[1::2]]
                dt.setdefault(title,{}).setdefault(periodlen,max(eachfreq))
with open('majorallelefreq.txt','w') as w:
    w.write(f'period\tmajorallelefreq\tSTRname\n')
    for k,v in dt.items():
        for k1,v1 in v.items():
            w.write(f'{k1}\t{v1}\t{k}\n')