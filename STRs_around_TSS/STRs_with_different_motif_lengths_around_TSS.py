from collections import OrderedDict
import pandas as pd
dt=OrderedDict()
import sys
lst=[]
xx,file2,file3=sys.argv
for i in range(1,7):
    for j in range(-10000,10000,100):
        dt.setdefault(i,{}).setdefault((j,j+100),[])
tempdt={}
with open(file3,'r') as r:
    for line in r:
        if line.startswith('STR'):
            name=line.strip().split(',')[0]
        elif line.startswith('period'):
            pdlen=len(line.strip().split(',')[1])
            tempdt.setdefault(name,pdlen)
with open(file2,'r') as r:
    for line in r:
        fline=line.strip().split('\t')
        print(line)
        number=fline[-2].split(':')[-1]
        length=tempdt.get(fline[0])
        if number.startswith('-') or number.isdigit():
            number = int(number)
            for k in dt.get(length):
                if k[0] <= number < k[1]:
                    dt.get(length).get(k).append(number)
                    break


resdt=OrderedDict()
for k,v in dt.items():
    for k1,v1 in v.items():
        resdt.setdefault(k,{}).setdefault(k1[0],len(v1))

df=pd.DataFrame.from_dict(resdt)
df=df.fillna(0)
df.to_csv('eachbpposition.csv',header=True)

df=pd.read_csv('eachbppositionmeijunyi.csv',index_col=0)
window_size = 5
step_size = 1


rolling_mean = df.rolling(window=window_size, min_periods=1).mean()

rolling_mean.to_csv('Reswindows500kbstep100.txt')