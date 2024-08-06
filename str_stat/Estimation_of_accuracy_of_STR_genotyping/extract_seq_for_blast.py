from pathlib import Path
import re
from Bio import SeqIO
p=re.compile(r'\d+')
pathway=Path('/home/tanxy/3krice/01.bwamem/gatkres/HipSTR')
dt={}
recdt={}
for rec in SeqIO.parse(r'Nipponbare.fa', 'fasta'):
    recdt.setdefault(rec.id,rec.seq)
for files in pathway.glob('leaveonefilterhipstrmergechrom*.csv'):
    with open(files,'r') as r:
        for line in r:
            fline=line.strip().split(',')
            if line.startswith('\n'):
                continue
            elif line.startswith(('hipSTR','lobSTR','gatkvg')):
                chrom=str(int(p.search(fline[0]).group(0)[:2]))
                varpos=int(p.search(fline[0]).group(0)[2:])
                seqofvar=recdt.get(chrom)[varpos-100:varpos]
                dt.setdefault(f'{chrom}-{fline[0]}',seqofvar)

with open(f'seq_of_{files.name}_100bp.txt','w') as w:
    for k,v in dt.items():
        w.write(f'>{k}\n{v}\n')

from pathlib import Path
import re
from Bio import SeqIO


def extract_sequences(pathway,file_name , fasta_file='Nipponbare.fa', output_suffix='_100bp.txt'):
    # 编译正则表达式
    p = re.compile(r'\d+')

    # 初始化字典
    dt = {}
    recdt = {}

    # 解析FASTA文件
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        recdt[rec.id] = rec.seq

    # 遍历指定路径下的文件
    for files in Path(pathway).glob(f'{file_name}chrom*.csv'):
        with open(files, 'r') as r:
            for line in r:
                if line.startswith('\n'):
                    continue
                elif line.startswith(('hipSTR', 'lobSTR', 'gatkvg')):
                    fline = line.strip().split(',')
                    chrom = str(int(p.search(fline[0]).group(0)[:2]))
                    varpos = int(p.search(fline[0]).group(0)[2:])
                    seqofvar = recdt.get(chrom)[varpos - 100:varpos]
                    dt[f'{chrom}-{fline[0]}'] = seqofvar

    # 写入输出文件
    output_file = Path(pathway) / f'seq_of_{file_name}{output_suffix}'
    with open(output_file, 'w') as w:
        for k, v in dt.items():
            w.write(f'>{k}\n{v}\n')


# 使用示例
extract_sequences('/hipSTR','hipstr')
extract_sequences('/gatk','gatk')
extract_sequences('/lobSTR','lobstr')
