from Bio import SeqIO
from collections import OrderedDict
from collections import Counter
dt=OrderedDict()

for rec in SeqIO.parse('Nipponbare.fa','fasta'):
    for i in range(0, len(rec.seq), 100000):
        if len(rec.seq) - i < 100000:
            seqall = rec.seq[i:len(rec.seq)]
            tempdt = Counter(seqall)
            GandC = tempdt.get('G') + tempdt.get('C')
            together = tempdt.get('G') + tempdt.get('C') + tempdt.get('A') + tempdt.get('T')
            gcratio = int(round(GandC / together, 2) * 100)
            dt.setdefault(f'Chr{str(rec.id).rjust(2, "0")}', {}).setdefault((i + 1, len(rec.seq)), {}).setdefault(
                'GCcontent', [gcratio])
            dt.setdefault(f'Chr{str(rec.id).rjust(2, "0")}', {}).setdefault((i + 1, len(rec.seq)), {}).setdefault('STRdensity',[])
            dt.setdefault(f'Chr{str(rec.id).rjust(2, "0")}', {}).setdefault((i + 1, len(rec.seq)), {}).setdefault('Genedensity',[])
        else:
            seqall = rec.seq[i:i + 100000]
            tempdt = Counter(seqall)
            GandC = tempdt.get('G') + tempdt.get('C')
            together = tempdt.get('G') + tempdt.get('C') + tempdt.get('A') + tempdt.get('T')
            gcratio = int(round(GandC / together, 2) * 100)
            dt.setdefault(f'Chr{str(rec.id).rjust(2, "0")}', {}).setdefault((i + 1, i + 100000), {}).setdefault(
                'GCcontent', [gcratio])
            dt.setdefault(f'Chr{str(rec.id).rjust(2, "0")}', {}).setdefault((i + 1, i + 100000), {}).setdefault(
                'STRdensity', [])
            dt.setdefault(f'Chr{str(rec.id).rjust(2, "0")}', {}).setdefault((i + 1, i + 100000), {}).setdefault(
                'Genedensity', [])

with open('locus.txt','r') as r:
    for line in r:
        fline = line.strip().split('\t')
        for k1,k2 in dt.get(fline[0]):
            if k1<=int(fline[1])<k2:
                dt.setdefault(fline[0],{}).setdefault((k1,k2), {}).setdefault('Genedensity',[]).append(int(fline[1]))

with open('Lastresult.txt','r') as r:
    for line in r:
        if line.startswith('\n'):
            continue
        elif line.startswith('lob') or line.startswith('vg') or line.startswith('STR'):
            chromandpos=p.search(line.strip().split(',')[0]).group(0)
            chrom=f'Chr{chromandpos[:2]}'
            pos=int(chromandpos[2:])
            for k1,k2 in dt.get(chrom):
                if k1<=pos<k2:
                    dt.setdefault(chrom,{}).setdefault(k,{}).setdefault('STRdensity',[]).append(pos)
                    break
with open('STRdensityandGenedensity.txt','w') as w:
    w.write('chromsome\tstart\tend\tgenedensity\tSTRdensity\tGCcontent\n')
    for k,v in dt.items():
        for k1,v1 in v.items():
            genedensity=len(v1.get('Genedensity'))
            STRdensity=v1.get('STRdensity')[0]
            GCcontent=v1.get('GCcontent')[0]
            w.write(f'{k}\t{k1[0]}\t{k1[1]}\t{genedensity}\t{STRdensity}\t{GCcontent}\n')

