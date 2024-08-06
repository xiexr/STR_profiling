from pathlib import Path
p=Path()
dt={}
for files in p.glob('cal_*.txt'):
    with open(files,'r') as r:
        filename = files.stem
        pop1_pop2=filename.split('_')[3:]
        for line in r:
            fline=line.strip().split('\t')
            if float(fline[-1])>=0.25:
                dt.setdefault(fline[0],{}).setdefault('_'.join(pop1_pop2),'1')
import pandas as pd
df=pd.DataFrame.from_dict(dt,orient='index')
df.fillna(0, inplace=True)
df.to_csv('output_Rst_0.25_upset_plot.csv')