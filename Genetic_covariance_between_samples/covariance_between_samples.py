import pandas as pd
import numpy as np
from statsmodels.stats import moment_helpers
df=pd.read_csv('newfiltersnpoutput.txt',header=0,index_col=0,sep='\t')

tempstd=[]
for each in df:
    test=df[each].std()
    tempstd.append(test)

corr_matrix=df.corr()

result=moment_helpers.corr2cov(corr_matrix,std=tempstd)
print(result)
print(type(result))
with open('CovarianceofSNP.txt','w') as w:
    w.write('\t')
    for each in df:
        if each[0].isdigit():
            w.write(f'IRIS_{each}\t')
        else:
            w.write(f'{each}\t')
    w.write('\n')
    columnname=list(df)
    i=0
    for each in result:
        if columnname[i][0].isdigit():
            w.write(f'IRIS_{columnname[i]}\t')
        else:
            w.write(f'{columnname[i]}\t')
        i+=1
        for each1 in each:
            w.write(f'{each1}\t')
        w.write('\n')