import re
from Bio import SeqIO
from pathlib import Path
from collections import OrderedDict

pathway = Path('eachchro/Lastresult_gatk')
p = re.compile(r'\d+')
totaldt = OrderedDict()
recdt = {}
for rec in SeqIO.parse(r'Nipponbare.fa', 'fasta'):
    recdt.setdefault(rec.id, rec.seq)
Flag = True
for files in pathway.glob('gatkchr*STR.csv'):
    totaldt = OrderedDict()
    with open(files, 'r') as r:
        refchrom = str(int(p.search(files.name).group()))
        for line in r:
            if line.startswith('\n'):
                continue
            elif line.startswith('gatkvg'):
                fline = line.strip().split(',')
                variantpos = int(p.search(fline[0]).group()[2:])
                startpos = int(fline[-2])
                endpos = int(fline[-1])
            elif line.startswith('ref'):  # vg0900314500
                threekrefseq = line.strip().split(',')[-1]
            elif line.startswith('seq'):
                if variantpos <= startpos:
                    if variantpos - 1 + len(threekrefseq) > endpos + 3 - 1:
                        refseq = recdt.get(refchrom)[variantpos - 3 - 1:variantpos - 1 + len(threekrefseq)]
                        threekseq = recdt.get(refchrom)[variantpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                            refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                    else:
                        refseq = recdt.get(refchrom)[variantpos - 3 - 1:endpos + 3 - 1]
                        threekseq = recdt.get(refchrom)[
                                    variantpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                            refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                else:
                    if variantpos - 1 + len(threekrefseq) > endpos + 3 - 1:
                        refseq = recdt.get(refchrom)[startpos - 3 - 1:variantpos - 1 + len(threekrefseq)]
                        threekseq = recdt.get(refchrom)[startpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                            refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                    else:
                        refseq = recdt.get(refchrom)[startpos - 3 - 1:endpos + 3 - 1]
                        threekseq = recdt.get(refchrom)[startpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                            refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                totaldt.setdefault(f'{fline[0]}-ref', [refseq])
                totaldt.setdefault(f'{fline[0]}-threek', [threekseq])
            elif line.startswith('period'):
                if Flag == False:
                    Flag = True
                    continue
                period = line.strip().split(',')[1]
                totaldt.setdefault(f'{fline[0]}-ref', []).append(period)
                totaldt.setdefault(f'{fline[0]}-threek', []).append(period)
    with open(f'{files.name[:-4]}forTRF.txt', 'w') as w:
        for k, v in totaldt.items():
            w.write(f'>{k}-{v[1]}\n{v[0]}\n')




pathway = Path('/HipSTR')
p = re.compile(r'\d+')
totaldt = OrderedDict()
recdt = {}
for rec in SeqIO.parse(r'Nipponbare.fa', 'fasta'):
    recdt.setdefault(rec.id, rec.seq)
threekrefdt = {}
with open('hipstr_ref.vcf', 'r') as r:
    r.readline()  # Skip the header line
    for line in r:
        fline = line.strip().split('\t')
        refseq = fline[3]
        altseqls = fline[4].split(',')
        nip = fline[-3].split(':')[0].split('|')
        skip_line = False

        for each in altseqls:
            if each == '.':
                skip_line = True
                break

        if skip_line:
            continue
        if nip[0] == './.':
            continue
        elif '.' in nip[0]:
            continue
        else:
            if nip[0] == nip[1]:
                if nip[0] == '0':
                    threekrefdt.setdefault(fline[0], {}).setdefault(fline[1], refseq)
                else:
                    threekrefdt.setdefault(fline[0], {}).setdefault(fline[1], altseqls[int(nip[0]) - 1])


Flag = True
for files in pathway.glob('hipstrchrom*.csv'):
    totaldt = OrderedDict()
    with open(files, 'r') as r:
        refchrom = str(int(p.search(files.name).group()))
        print(refchrom)
        for line in r:
            if line.startswith('\n'):
                continue
            elif line.startswith('hipSTR'):
                fline = line.strip().split(',')
                variantpos = int(p.search(fline[0]).group()[2:])
                startpos = int(fline[-2])
                endpos = int(fline[-1])
            elif line.startswith('period'):
                if threekrefdt.get(refchrom).get(str(variantpos)):
                    threekrefseq = threekrefdt.get(refchrom).get(str(variantpos))
                    if variantpos < startpos:
                        if variantpos - 1 + len(threekrefseq) > endpos + 3 - 1:
                            refseq = recdt.get(refchrom)[variantpos - 3 - 1:variantpos - 1 + len(threekrefseq)]
                            threekseq = recdt.get(refchrom)[
                                        variantpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                                refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                        else:
                            refseq = recdt.get(refchrom)[variantpos - 3 - 1:endpos + 3 - 1]
                            threekseq = recdt.get(refchrom)[
                                        variantpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                                refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                    else:
                        if variantpos - 1 + len(threekrefseq) > endpos + 3 - 1:
                            refseq = recdt.get(refchrom)[startpos - 3 - 1:variantpos - 1 + len(threekrefseq)]
                            threekseq = recdt.get(refchrom)[startpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                                refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                        else:
                            refseq = recdt.get(refchrom)[startpos - 3 - 1:endpos + 3 - 1]
                            threekseq = recdt.get(refchrom)[startpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                                refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                    totaldt.setdefault(f'{fline[0]}-ref', [refseq])
                    totaldt.setdefault(f'{fline[0]}-threek', [threekseq])
                    period = line.strip().split(',')[1]
                    totaldt.setdefault(f'{fline[0]}-ref', []).append(period)
                    totaldt.setdefault(f'{fline[0]}-threek', []).append(period)
    with open(f'{files.name[:-4]}forTRF.txt', 'w') as w:
        for k, v in totaldt.items():
            w.write(f'>{k}-{v[1]}\n{v[0]}\n')



pathway = Path('lobstr_res')
p = re.compile(r'\d+')
totaldt = OrderedDict()
recdt = {}
for rec in SeqIO.parse(r'Nipponbare.fasta', 'fasta'):
    recdt.setdefault(rec.id, rec.seq)
threekrefdt = {}

with open('lobSTR_ref.vcf','r') as r:
    r.readline()  # Skip the header line
    for line in r:
        fline = line.strip().split('\t')
        refseq = fline[3]
        altseqls = fline[4].split(',')
        nip = fline[-3].split(':')[0].split('/')
        skip_line = False

        for each in altseqls:
            if each == '.':
                skip_line = True
                break

        if skip_line:
            continue
        if nip[0] == './.':
            continue
        elif '.' in nip[0]:
            continue
        else:
            if nip[0] == nip[1]:
                if nip[0] == '0':
                    threekrefdt.setdefault(fline[0], {}).setdefault(fline[1], refseq)
                else:
                    threekrefdt.setdefault(fline[0], {}).setdefault(fline[1], altseqls[int(nip[0]) - 1])

Flag = True
for files in pathway.glob('lobSTRchrom*.csv'):
    totaldt = OrderedDict()
    with open(files, 'r') as r:
        refchrom = str(int(p.search(files.name).group()))
        for line in r:
            if line.startswith('\n'):
                continue
            elif line.startswith('lobSTR'):
                fline = line.strip().split(',')
                variantpos = int(p.search(fline[0]).group()[2:])
                startpos = int(fline[-2])
                endpos = int(fline[-1])
            elif line.startswith('period'):
                if threekrefdt.get(refchrom).get(str(variantpos)):
                    threekrefseq = threekrefdt.get(refchrom).get(str(variantpos))
                    if variantpos < startpos:
                        if variantpos - 1 + len(threekrefseq) > endpos + 3 - 1:
                            refseq = recdt.get(refchrom)[variantpos - 3 - 1:variantpos - 1 + len(threekrefseq)]
                            threekseq = recdt.get(refchrom)[
                                        variantpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                                refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                        else:
                            refseq = recdt.get(refchrom)[variantpos - 3 - 1:endpos + 3 - 1]
                            threekseq = recdt.get(refchrom)[
                                        variantpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                                refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                    else:
                        if variantpos - 1 + len(threekrefseq) > endpos + 3 - 1:
                            refseq = recdt.get(refchrom)[startpos - 3 - 1:variantpos - 1 + len(threekrefseq)]
                            threekseq = recdt.get(refchrom)[startpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                                refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                        else:
                            refseq = recdt.get(refchrom)[startpos - 3 - 1:endpos + 3 - 1]
                            threekseq = recdt.get(refchrom)[startpos - 3 - 1:variantpos - 1] + threekrefseq + recdt.get(
                                refchrom)[variantpos - 1 + len(threekrefseq):endpos + 3 - 1]
                    totaldt.setdefault(f'{fline[0]}-ref', [refseq])
                    totaldt.setdefault(f'{fline[0]}-threek', [threekseq])
                    period = line.strip().split(',')[1]
                    totaldt.setdefault(f'{fline[0]}-ref', []).append(period)
                    totaldt.setdefault(f'{fline[0]}-threek', []).append(period)
    with open(f'{files.name[:-4]}forTRF.txt', 'w') as w:
        for k, v in totaldt.items():
            w.write(f'>{k}-{v[1]}\n{v[0]}\n')