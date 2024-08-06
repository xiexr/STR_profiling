dt={}
tempdt={}
with open('Lastresult.txt','r') as r:
    for line in r:
        if line.startswith('\n'):
            continue
        elif line.startswith('period'):
            fline = line.strip().split(',')
            repeattime=fline[-1].split('-')
            period=fline[1]
        elif line[0].isdigit():
            fline = line.strip().split(',')
            if fline[1]=='All':
                for index,each in enumerate(fline[0].split('-')[::2]):
                    tempdt.setdefault(repeattime[index],0)
                    tempdt[repeattime[index]]+=int(each)
                max_key = max(tempdt, key=tempdt.get)
                dt.setdefault(period,[].append(max_key))
                tempdt={}
with open('Repeat_unit_major_allele_repeat_res.txt','w') as w:
    for k,v in dt.items():
        for each in v:
            w.write(f'{each}\t{k}\n')

