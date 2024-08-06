import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
samplesize=[]
with open('RealSTRresultaddmore.csv','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split('\t')
        samplesize.append(int(float(fline[5])))

##########################作图部分#########################
f, ax = plt.subplots(figsize = (8, 4))
plt.rcParams['font.family'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
ax.hist(samplesize,facecolor='#98FB98' , edgecolor = 'black', color = 'white', bins = 100)
ax.axis([0,130,0,150000])
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=10))
ax.set_yticklabels(['0','2','4','6','8','10','12'],family = 'Arial',fontsize = 14)
labels=ax.get_xticklabels()
for label in labels:
    label.set_fontsize(14)
    label.set_fontname('Arial')

plt.text(0, 150000, r'x 10$^{4}$',fontsize=14,fontname='Arial')
ax.set_ylabel('Tests',fontsize=14,fontname='Arial')
ax.set_xlabel('Sample size',fontsize=14,fontname='Arial')
plt.tight_layout()
plt.savefig('Samplesizes.png',dpi=300)

########################################临时统计一下100以上的sample数量占比##################
