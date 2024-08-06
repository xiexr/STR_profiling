refdiffdt={}
import re
p=re.compile(r'\d+')
Flag=False
with open('Lastresult.txt','r') as r:
    for line in r:
        fline=line.strip().split(',')
        if line.startswith('\n'):
            continue
        elif line.startswith('lobSTR') or line.startswith('hipSTR') or line.startswith('gatkvg'):
            title=p.search(line.strip().split(',')[0]).group()
            Flag=False
        elif Flag:
            continue
        elif line.startswith('ref'):
            refseq=fline[1]
        elif line.startswith('seq'):
            if refseq in fline[1].split('-'):
                refseqindex=fline[1].split('-').index(refseq)
            else:
                Flag=True
        elif line.startswith('period'):
            periodlen=len(fline[1])
            cpls=fline[-1].split('-')
            cpls=list(map(int,cpls))
            refseqrptime=cpls[refseqindex]
        elif line[0].isdigit():
            fline = line.strip().split(',')
            if fline[-7]=='All':
                eachfreq=fline[0].split('-')[1::2]
                eachfreq1=[float(each.split('%')[0])*0.01 for each in eachfreq]
                majorindex=eachfreq1.index(max(eachfreq1))
                majorallelerptime=cpls[majorindex]
                differenceofrepeattime=majorallelerptime-refseqrptime
                refdiffdt.setdefault(title, {}).setdefault(periodlen,differenceofrepeattime)
with open('majoralleledifffromref.txt','w') as w:
    w.write('title\tmajordifffromref\tperiodlen\n')
    for k,v in refdiffdt.items():
        for k1,v1 in v.items():
            w.write(f'{k}\t{v1}\t{k1}\n')

