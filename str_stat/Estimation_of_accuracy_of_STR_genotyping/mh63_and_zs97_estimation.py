import re
from Bio import SeqIO
from pathlib import Path
from collections import OrderedDict
pathway=Path('/HipSTR')
p=re.compile(r'\d+')
totaldt=OrderedDict()
recdt={}
for rec in SeqIO.parse(r'/MH63RS3.fasta', 'fasta'):
    recid=str(int(p.search(str(rec.id)).group()))
    recdt.setdefault(recid,rec.seq)
threekrefdt={}
with open('hipSTRref.vcf', 'r') as r:
    r.readline()  # Skip the header line
    for line in r:
        fline = line.strip().split('\t')
        refseq = fline[3]
        altseqls = fline[4].split(',')
        # zs97=fline[-2].split(':')[0].split('|')
        mh63=fline[-1].split(':')[0].split('|')
        skip_line = False

        for each in altseqls:
            if each == '.':
                skip_line = True
                break

        if skip_line:
            continue
        if mh63[0]=='./.':#or zs97
            continue
        elif '.' in mh63[0]:
            continue
        else:
            # print(nip)
            if mh63[0]==mh63[1]:
                if mh63[0]=='0':
                    threekrefdt.setdefault(fline[0],{}).setdefault(fline[1],refseq)
                else:
                    threekrefdt.setdefault(fline[0], {}).setdefault(fline[1], altseqls[int(mh63[0])-1])

def zsormh(filename):
    dt={}
    with open(filename,'r') as r:
        for line in r:
            fline =line.strip().split('\t')
            chrom=fline[0].split('-')[0]
            varpos=str(int(p.search(fline[0].split('-')[-1]).group(0)[2:]))
            dt.setdefault(chrom,{}).setdefault(varpos,fline[9])
    return dt
mh63ref=zsormh(filename)
# print(mh63ref)
Flag=False
for files in pathway.glob('hipSTRchrom*.csv'):
    totaldt=OrderedDict()
    with open(files,'r') as r:
        refchrom=str(int(p.search(files.name).group()))
        print(refchrom)
        # print(files.name)
        for line in r:
            if line.startswith('\n'):
                continue
            elif line.startswith('hipSTR'):
                fline= line.strip().split(',')
                nipchrom=str(int(p.search(fline[0]).group()[:2]))
                nipvariantpos=str(int(p.search(fline[0]).group()[2:]))
                nipstartpos=int(fline[-2])
                nipendpos=int(fline[-1])
                if mh63ref.get(nipchrom).get(nipvariantpos):
                    anothervariantpos=int(mh63ref.get(nipchrom).get(nipvariantpos))
                    if anothervariantpos:
                        Flag=True
                        if int(nipvariantpos)>int(nipstartpos):
                            anotherstartpos=int(anothervariantpos-(int(nipvariantpos)-nipstartpos))
                        else:
                            anotherstartpos=int(anothervariantpos+(nipstartpos-int(nipvariantpos)))

                        anotherendpos=int(anothervariantpos+(nipendpos-int(nipvariantpos)))
            elif line.startswith('period'):
                if Flag:
                    Flag=False
                    if threekrefdt.get(refchrom).get(str(nipvariantpos)):
                        threekrefseq = threekrefdt.get(refchrom).get(str(nipvariantpos))
                        if anothervariantpos < anotherstartpos:
                            if anothervariantpos - 1 + len(threekrefseq) > anotherendpos + 3 - 1:
                                refseq = recdt.get(refchrom)[anothervariantpos - 3 - 1:anothervariantpos - 1 + len(threekrefseq)]
                                threekseq = recdt.get(refchrom)[anothervariantpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                            else:
                                refseq = recdt.get(refchrom)[anothervariantpos - 3 - 1:anotherendpos + 3 - 1]
                                threekseq = recdt.get(refchrom)[
                                            anothervariantpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                        else:
                            if anothervariantpos - 1 + len(threekrefseq) > anotherendpos + 3 - 1:
                                refseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anothervariantpos - 1 + len(threekrefseq)]
                                threekseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                            else:
                                refseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anotherendpos + 3 - 1]
                                threekseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                        totaldt.setdefault(f'{fline[0]}-ref', [refseq])
                        totaldt.setdefault(f'{fline[0]}-threek', [threekseq])
                        period=line.strip().split(',')[1]
                        totaldt.setdefault(f'{fline[0]}-ref', []).append(period)
                        totaldt.setdefault(f'{fline[0]}-threek', []).append(period)
    with open(f'{files.name[:-4]}_{mh63orzs97}_forTRF.txt','w') as w:
        for k,v in totaldt.items():
            w.write(f'>{k}-{v[1]}\n{v[0]}\n')

import re
from Bio import SeqIO
from pathlib import Path
from collections import OrderedDict
pathway=Path('/lobSTR')
p=re.compile(r'\d+')
totaldt=OrderedDict()
recdt={}
for rec in SeqIO.parse(r'MH63RS3.fasta', 'fasta'):
    recid=str(int(p.search(str(rec.id)).group()))
    recdt.setdefault(recid,rec.seq)
threekrefdt={}
with open('lobSTRref.vcf', 'r') as r:
    r.readline()  # Skip the header line
    for line in r:
        fline = line.strip().split('\t')
        refseq = fline[3]
        altseqls = fline[4].split(',')
        # zs97=fline[-2].split(':')[0].split('|')
        mh63=fline[-1].split(':')[0].split('/')
        skip_line = False

        for each in altseqls:
            if each == '.':
                skip_line = True
                break

        if skip_line:
            continue
        if mh63[0]=='./.':
            continue
        elif '.' in mh63[0]:
            continue
        else:
            if mh63[0]==mh63[1]:
                if mh63[0]=='0':
                    threekrefdt.setdefault(fline[0],{}).setdefault(fline[1],refseq)
                else:
                    threekrefdt.setdefault(fline[0], {}).setdefault(fline[1], altseqls[int(mh63[0])-1])

def zsormh(filename):
    dt={}
    with open(filename,'r') as r:
        for line in r:
            fline =line.strip().split('\t')
            chrom=fline[0].split('-')[0]
            varpos=str(int(p.search(fline[0].split('-')[-1]).group(0)[2:]))
            dt.setdefault(chrom,{}).setdefault(varpos,fline[9])
    return dt
mh63ref=zsormh(filename)
Flag=False
for files in pathway.glob('lobSTRchrom*.csv'):
    totaldt=OrderedDict()
    with open(files,'r') as r:
        refchrom=str(int(p.search(files.name).group()))
        print(refchrom)
        for line in r:
            if line.startswith('\n'):
                continue
            elif line.startswith('lobSTR'):
                fline= line.strip().split(',')
                nipchrom=str(int(p.search(fline[0]).group()[:2]))
                nipvariantpos=str(int(p.search(fline[0]).group()[2:]))
                nipstartpos=int(fline[-2])
                nipendpos=int(fline[-1])
                if mh63ref.get(nipchrom).get(nipvariantpos):
                    anothervariantpos=int(mh63ref.get(nipchrom).get(nipvariantpos))
                    if anothervariantpos:
                        Flag=True
                        if int(nipvariantpos)>int(nipstartpos):
                            anotherstartpos=int(anothervariantpos-(int(nipvariantpos)-nipstartpos))
                        else:
                            anotherstartpos=int(anothervariantpos+(nipstartpos-int(nipvariantpos)))

                        anotherendpos=int(anothervariantpos+(nipendpos-int(nipvariantpos)))
            elif line.startswith('period'):
                if Flag:
                    Flag=False
                    if threekrefdt.get(refchrom).get(str(nipvariantpos)):
                        threekrefseq = threekrefdt.get(refchrom).get(str(nipvariantpos))
                        if anothervariantpos < anotherstartpos:
                            if anothervariantpos - 1 + len(threekrefseq) > anotherendpos + 3 - 1:
                                refseq = recdt.get(refchrom)[anothervariantpos - 3 - 1:anothervariantpos - 1 + len(threekrefseq)]
                                threekseq = recdt.get(refchrom)[anothervariantpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                            else:
                                refseq = recdt.get(refchrom)[anothervariantpos - 3 - 1:anotherendpos + 3 - 1]
                                threekseq = recdt.get(refchrom)[
                                            anothervariantpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                        else:
                            if anothervariantpos - 1 + len(threekrefseq) > anotherendpos + 3 - 1:
                                refseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anothervariantpos - 1 + len(threekrefseq)]
                                threekseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                            else:
                                refseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anotherendpos + 3 - 1]
                                threekseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                        totaldt.setdefault(f'{fline[0]}-ref', [refseq])
                        totaldt.setdefault(f'{fline[0]}-threek', [threekseq])
                        period=line.strip().split(',')[1]
                        totaldt.setdefault(f'{fline[0]}-ref', []).append(period)
                        totaldt.setdefault(f'{fline[0]}-threek', []).append(period)
    with open(f'{files.name[:-4]}_{zs97ormh63}_forTRF.txt','w') as w:
        for k,v in totaldt.items():
            w.write(f'>{k}-{v[1]}\n{v[0]}\n')



import re
from Bio import SeqIO
from pathlib import Path
from collections import OrderedDict
pathway=Path('/gatk_res')
p=re.compile(r'\d+')
totaldt=OrderedDict()
recdt={}
for rec in SeqIO.parse(r'MH63RS3.fasta', 'fasta'):
    recid=str(int(p.search(str(rec.id)).group()))
    recdt.setdefault(recid,rec.seq)
threekrefdt={}
with open('gatkref.vcf', 'r') as r:
    r.readline()  # Skip the header line
    for line in r:
        fline = line.strip().split('\t')
        chrom=str(int(fline[0][2:4]))
        pos=str(int(fline[0][4:]))
        refseq=fline[-1]
        if refseq!='N' and refseq!='H':
            threekrefdt.setdefault(chrom,{}).setdefault(pos,refseq)
def zsormh(filename):
    dt={}
    with open(filename,'r') as r:
        for line in r:
            fline =line.strip().split('\t')
            chrom=fline[0].split('-')[0]
            # print(str(int(p.search(fline[0].split('-')[-1]).group(0)[2:])))
            varpos=str(int(p.search(fline[0].split('-')[-1]).group(0)[2:]))
            dt.setdefault(chrom,{}).setdefault(varpos,fline[9])
    return dt
mh63ref=zsormh(filename)
Flag=False
for files in pathway.glob('gatkchrom*STR.csv'):
    totaldt=OrderedDict()
    with open(files,'r') as r:
        refchrom=str(int(p.search(files.name).group()))
        print(refchrom)
        # print(files.name)
        for line in r:
            if line.startswith('\n'):
                continue
            elif line.startswith('vg'):
                fline= line.strip().split(',')
                nipchrom=str(int(p.search(fline[0]).group()[:2]))
                nipvariantpos=str(int(p.search(fline[0]).group()[2:]))
                nipstartpos=int(fline[-2])
                nipendpos=int(fline[-1])
                if mh63ref.get(nipchrom).get(nipvariantpos):
                    anothervariantpos=int(mh63ref.get(nipchrom).get(nipvariantpos))
                    if anothervariantpos:
                        Flag=True
                        if int(nipvariantpos)>int(nipstartpos):
                            anotherstartpos=int(anothervariantpos-(int(nipvariantpos)-nipstartpos))
                        else:
                            anotherstartpos=int(anothervariantpos+(nipstartpos-int(nipvariantpos)))

                        anotherendpos=int(anothervariantpos+(nipendpos-int(nipvariantpos)))
            elif line.startswith('period'):
                if Flag:
                    Flag=False
                    if threekrefdt.get(refchrom).get(str(nipvariantpos)):
                        threekrefseq = threekrefdt.get(refchrom).get(str(nipvariantpos))
                        if anothervariantpos < anotherstartpos:
                            if anothervariantpos - 1 + len(threekrefseq) > anotherendpos + 3 - 1:
                                refseq = recdt.get(refchrom)[anothervariantpos - 3 - 1:anothervariantpos - 1 + len(threekrefseq)]
                                threekseq = recdt.get(refchrom)[anothervariantpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                            else:
                                refseq = recdt.get(refchrom)[anothervariantpos - 3 - 1:anotherendpos + 3 - 1]
                                threekseq = recdt.get(refchrom)[
                                            anothervariantpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                        else:
                            if anothervariantpos - 1 + len(threekrefseq) > anotherendpos + 3 - 1:
                                refseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anothervariantpos - 1 + len(threekrefseq)]
                                threekseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                            else:
                                refseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anotherendpos + 3 - 1]
                                threekseq = recdt.get(refchrom)[anotherstartpos - 3 - 1:anothervariantpos - 1] + threekrefseq + recdt.get(
                                    refchrom)[anothervariantpos - 1 + len(threekrefseq):anotherendpos + 3 - 1]
                        totaldt.setdefault(f'{fline[0]}-ref', [refseq])
                        totaldt.setdefault(f'{fline[0]}-threek', [threekseq])
                        period=line.strip().split(',')[1]
                        totaldt.setdefault(f'{fline[0]}-ref', []).append(period)
                        totaldt.setdefault(f'{fline[0]}-threek', []).append(period)
    with open(f'{files.name[:-4]}_{mh63orzs97}_forTRF.txt','w') as w:
        for k,v in totaldt.items():
            w.write(f'>{k}-{v[1]}\n{v[0]}\n')
