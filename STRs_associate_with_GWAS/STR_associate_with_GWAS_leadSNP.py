import re
snp_str_dt={}
p=re.compile(r'\d+')
with open('SNP_STR_LD_value.txt','r') as r:
    for line in r:
        fline=line.strip().split('\t')
        ldval=float(fline[-1])
        if float(ldval)> 0.7:
            snp_str_dt.setdefault(f'{fline[2]}:{fline[3]}',{}).setdefault(f'{fline[0]}:{fline[1]}',ldval)
tempdt={}
with open('totalGWASSNP_of_ZQF.txt','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        tempdt.setdefault(f'{fline[1]}:{fline[0]}',[fline[2],fline[3],fline[4],fline[5]])
tempdt1={}
with open('totalGWAS_Hanbing.txt','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        tempdt1.setdefault(f'{fline[1]}:{fline[0]}',[fline[2],fline[3],fline[-1]])

with open('STR_SNP_combination_res.txt','w') as w:
    w.write('STRID\tchr:STRpos\tSNPID\tchr:SNPpos\tTrait\tLRP_pvalue\tLMM_pvalue\tMLM_pvalue(-log10)\tLD\n')
    for k,v in snp_str_dt.items():
        for k1,v1 in v.items():
            if k1 in tempdt:
                STRpos=p.search(k.split(":")[-1]).group()[2:]
                SNPpos=tempdt.get(k1)[0]
                trait=tempdt.get(k1)[-1]
                LRPvalue=tempdt.get(k1)[1]
                LMMpvalue=tempdt.get(k1)[2]
                w.write(f'{k.split(":")[-1]}\t{k.split(":")[0]}:{STRpos}\t{k1.split(":")[1]}\t{k1.split(":")[0]}:{SNPpos}\t{trait}\t{LRPvalue}\t{LMMpvalue}\tNA\t{v1}\n')
    for k, v in snp_str_dt.items():
        for k1, v1 in v.items():
            if k1 in tempdt1:
                STRpos = p.search(k.split(":")[-1]).group()[2:]
                SNPpos = tempdt1.get(k1)[0]
                trait = tempdt1.get(k1)[-1]
                MLMpvalue = tempdt1.get(k1)[1]
                w.write(
                    f'{k.split(":")[-1]}\t{k.split(":")[0]}:{STRpos}\t{k1.split(":")[1]}\t{k1.split(":")[0]}:{SNPpos}\t{trait}\tNA\tNA\t{MLMpvalue}\t{v1}\n')

dt={}
with open('STR_SNP_combination_res.txt','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        dt.setdefault(fline[4],0)
        dt[fline[4]]+=1
with open('pheno_num_with_STRs.txt','w') as w:
    w.write('pheno\tnum\n')
    for k,v in dt.items():
        w.write(f'{k}\t{v}\n')