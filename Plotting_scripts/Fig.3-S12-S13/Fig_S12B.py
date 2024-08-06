from collections import OrderedDict
import pandas as pd, seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
SDAT3_PATH= root + 'GeneticCovariancetrans.csv'
K = pd.read_csv(SDAT3_PATH, sep = '\t', index_col = 'Unnamed: 0')
kdict=K.to_dict()
lst=[6,7,9,14,27,28,29,30,32,46,17,18,19,26,34,39,57,58,59,60,99,16,33,43,47,49,53,54,55,65,73,75,76,80,81,85,92,93,95,96,97,98,101,102,104,106,107,109,110,111,112,114,115,116,117,118,119,120,122,123,125,126,127,23,31,64,69,70,74,77,88,94,124,24,40,63,67,71,72,90,121,41,44,50,51,62,66,68,113,37,78,79,82,84,89,91,103,108,1,2,3,4,5,8,10,11,13,20,22,12,15,21,25,35,36,38,42,45,48,52,56,61,83,86,87,100,105]
lst=list(map(str,lst))
tempdt=OrderedDict()
for k in lst:
    for i in lst:
        value=kdict.get(k).get(int(i))
        tempdt.setdefault(k,{}).setdefault(i,value)
newK=pd.DataFrame.from_dict(tempdt)
newK.to_csv('GeneticCovariancetransfromindicatojap.tsv',sep='\t')
plt.figure(dpi=300)
sns.clustermap(newK, cmap=sns.dark_palette("#2ecc71", as_cmap=True,reverse=True),cbar_kws={'aspect':0.5},xticklabels=False,yticklabels=False, row_cluster=False,col_cluster=False)# tree_kws={'linewidths':0}
plt.savefig('myplottest.png')