genewitheSTR={}
with open('eSTRwithFDR.csv','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        genename=fline[1]
        if float(fline[-1])<0.05:
            if 'upstream' in fline[-2]:
                if int(fline[-2].split('_')[-1]) < 5000 :
                    genewitheSTR.setdefault(genename,'')
            elif 'downstream' in fline[-2]:
                if int(fline[-2].split('_')[-1]) < 2000 :
                    genewitheSTR.setdefault(genename,'')
            else:
                genewitheSTR.setdefault(genename, '')
with open('extractgenewitheSTRforGOandKEGG.txt','w') as w:
    for k in genewitheSTR:
        w.write(f'{k}\n')