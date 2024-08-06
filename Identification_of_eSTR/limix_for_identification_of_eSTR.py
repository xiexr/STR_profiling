import glob
import limix
import pickle
import pandas as pd
import pandas_plink
import seaborn as sns
import pybedtools as pyb
import numpy as np
import scipy
import h5py
from pylab import *
from collections import Counter
from itertools import product
from sklearn.preprocessing import scale
from limix.qtl import scan, iscan
from statsmodels.stats import moment_helpers
from statsmodels.regression.linear_model import OLS
from scipy.stats import jarque_bera, normaltest, anderson
SDAT2_PATH='DiploidUnitNumberCallsresorttrans.csv'
SDAT3_PATH='GeneticCovariancetrans.csv'
XDAT1_PATH='GenewithSTR100kbinsidepy37.pkl'
XDAT2_PATH='expressiondatatrans.csv'
XDAT3_PATH='Oryzasativa.genes.bed'
# gene expression data from Kawakatsu et al., zeros removed and log(e) transformed
X= pd.read_csv(XDAT2_PATH, sep = '\t', index_col = 'Unnamed: 0')
fXd= X.to_dict()
# STR variant calling results, units counted and combined for both strands

G= pd.read_csv(SDAT2_PATH, sep = '\t', index_col = 'Unnamed: 0')
# genetic covariance
K= pd.read_csv(SDAT3_PATH, sep = '\t', index_col = 'Unnamed: 0')

K.index = K.index.astype(int)
K.columns = K.columns.astype(int)
# dictionary, genes as keys and STRs within 100 kb as values
closed= pickle.load(open(XDAT1_PATH, 'rb'))
def CommonAlleles(G, loci):
    count = Counter(G[loci])
    countvalues = np.array(list(count.values()))
    countkeys = np.array(list(count.keys()))
    commonalleles = countkeys[(countvalues / len(G)) > 0.05]
    return G[loci].isin(commonalleles)#[7,8,9]
common_loci = {}
for loci in G.columns:
    common_loci[loci] = G[loci][CommonAlleles(G, loci)]#只保留重复次数中频率大于0.05的样品,过滤nan和频率小的


def runLimix(path, genes,X, G, closed, common, K, common_loci, genes_tested):
    n = 0
    for gene in genes:
        n += 1
        if gene in closed and gene not in genes_tested:
            res = pd.DataFrame()  #
            y = pd.Series(X[gene]).dropna()  #
            loci = closed[gene]  #
            res['loci'] = loci
            res['gene'] = [gene] * len(res)
            effsizes = []  #
            effsizes_ses = []  #
            pvs = []  #
            shapiro_ps = []  #
            sample_sizes = []  #
            for locus in loci:
                if common == True:
                    g = common_loci[locus]  #
                else:
                    g = G[locus]
                # get common indices
                common_indices = list(set(y.index).intersection(g.index))  #既有表达量数据同时某种重复次数频率大于0.05
                if len(common_indices)<7:
                    res = res[~res['loci'].isin([locus])]
                    continue
                g_reindexed = g.reindex(common_indices)  #去掉重复次数中没有达到要求的样品的series
                y_reindexed = y.reindex(common_indices)  #去掉表达量中没有达到要求的样品
                k_reindexed = K.reindex(common_indices)  #去掉相关系数中没有达到要求的样品
                k_reindexed = k_reindexed[k_reindexed.index]  #去掉另一列的样品
                # if genetic variation
                if g_reindexed.nunique() > 1:
                    # run model
                    r = scan(G=g_reindexed,  #
                             Y=y_reindexed,
                             lik='normal',
                             M=None,
                             K=k_reindexed.values,
                             verbose=False)
                    eff = r.effsizes['h2']  #
                    eff = eff[eff['effect_type'] == 'candidate']  #
                    effsize = eff.effsize.values[0]  #
                    effsize_se = eff.effsize_se.values[0]  #
                    pv = r.stats.pv20.values[0]  #
                    pvs.append(pv)  #
                    sample_sizes.append(len(g_reindexed))  #
                    effsizes.append(effsize)  #
                    effsizes_ses.append(effsize_se)  #
                else:
                    effsizes.append(nan)
                    effsizes_ses.append(nan)
                    pvs.append(nan)
                    sample_sizes.append(nan)
            res = res.reset_index(drop=True)
            res['effsize'] = effsizes
            res['effsize_se'] = effsizes_ses
            res['p'] = pvs
            res['ss'] = sample_sizes
            res.dropna().to_csv(path + 'results.' + gene + '.tsv', sep='\t', index=None)
            res = None
            effsizes = None
            effsizes_ses = None
            pvs = None
            shapiro_ps = None
            sample_sizes = None
            r = None
            loci = None
            g = None
            common_indices = None
            g_reindexed = None
            y_reindexed = None
            k_reindexed = None
            eff = None
            effsize = None
            effsize_se = None
            pv = None
path = 'sampledatarun/'
genes_tested = []
for f in glob.glob(path + '*'):
    genes_tested.append(f.split('.')[-2])
X_sample = X
X_sample.head()
runLimix(path= path,genes=X_sample.columns,X = fXd,G= G,closed = closed,common = True,K = K,common_loci  = common_loci,genes_tested = genes_tested)
frames = []
for f in glob.glob(path + '*'):
    frames.append(pd.read_csv(f, sep = '\t', index_col = None))
realstrs = pd.concat(frames)
realstrs.head()
realstrs.to_csv('RealSTRresult.csv',header=True,sep='\t')