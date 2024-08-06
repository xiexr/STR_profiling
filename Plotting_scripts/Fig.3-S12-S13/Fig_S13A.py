import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.mlab as mlab
import  numpy as np
from scipy.stats import norm
templs=[]
dt={}
with open('RealSTRresultaddmore.csv','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        dt.setdefault(fline[1],0)
        dt[fline[1]]+=1
for k,v in dt.items():
    templs.append(v)
mu=sum(templs)/len(templs)
f, ax = plt.subplots(figsize = (8, 4))
n,bins,patches=ax.hist(templs, bins = 100,density=1,facecolor='#FFCCFF',edgecolor = 'black', color = 'white')
arr_std=np.std(templs)
y=norm.pdf(bins,mu,arr_std)
ax.plot(bins,y,'r--')
ax.axis([0,140,0.000,0.035])
ax.set_yticklabels(['0','0.5','1','1.5','2','2.5','3','3.5'],family = 'Arial',fontsize = 14)
labels=ax.get_xticklabels()
for label in labels:
    label.set_fontsize(14)
    label.set_fontname('Arial')
# sns.despine(offset = 10)
plt.text(0, 0.035, r'%',fontsize=14,fontname='Arial')
plt.ylabel('The proportion of genes tested',fontname='Arial',fontsize=14)
plt.xlabel('Sites',fontname='Arial',fontsize=14)
plt.tight_layout()
plt.savefig('Meannumberoftestedsites.png',dpi=300)
#Mean number of tested sites: 55