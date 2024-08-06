
x1,x2,x3,x4,x5,x6=[],[],[],[],[],[]
y1,y2,y3,y4,y5,y6=[],[],[],[],[],[]
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
with open('Result.txt','r') as r:
    r.readline()
    for line in r:
        fline=line.strip().split(',')
        x1.append(float(fline[0]))
        x2.append(float(fline[0]))
        x3.append(float(fline[0]))
        x4.append(float(fline[0]))
        x5.append(float(fline[0]))
        x6.append(float(fline[0]))
        y1.append(float(fline[1]))
        y2.append(float(fline[2]))
        y3.append(float(fline[3]))
        y4.append(float(fline[4]))
        y5.append(float(fline[5]))
        y6.append(float(fline[6]))

plt.figure(figsize=(12,6),dpi=300)
plt.rcParams['font.family'] = ['Arial']

# plt.plot(x,y,color='red',linewidth=3.0,label='STR')
plt.plot(x1,y1,color='#5B9BD5',linewidth=3.0,label='1 bp')
plt.plot(x2,y2,color='#ED7D31',linewidth=3.0,label='2 bp')
plt.plot(x3,y3,color='#A5A5A5',linewidth=3.0,label='3 bp')
plt.plot(x4,y4,color='#FFC000',linewidth=3.0,label='4 bp')
plt.plot(x5,y5,color='#4472C4',linewidth=3.0,label='5 bp')
plt.plot(x6,y6,color='#70AD47',linewidth=3.0,label='6 bp')

plt.title('str',fontdict={'weight':'normal','size': 30} )
plt.xlabel('Position (bp)',fontdict={'weight':'normal','size': 25})
plt.ylabel('Number',fontdict={'weight':'normal','size': 25})
plt.tick_params(labelsize=15)
plt.legend(loc=1,fontsize=15)
plt.savefig('output.png',dpi=600)
plt.savefig('output.pdf')
plt.show()
