
# coding: utf-8

# In[17]:


import pickle
import pandas as pd
from scipy import stats
from multiprocessing import Pool
from statsmodels.stats.multitest import multipletests
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import h5py
import os
from scipy.stats import t as t_dist
from functools import partial


# In[23]:


def vcorrcoef(target, X, data):
    chunk_index = data[0]
    Y = data[1]
    Xm = np.reshape(np.mean(X,axis=1),(X.shape[0],1))
    Ym = np.reshape(np.mean(Y,axis=1),(Y.shape[0],1))
    r_num = (X-Xm).dot((Y-Ym).T)
    r_den = np.sqrt(np.sum((X-Xm)**2,axis=1)[:, np.newaxis]).dot(np.sqrt(np.sum((Y-Ym)**2,axis=1))[:, np.newaxis].T)
    r = r_num/r_den

    with h5py.File("../all_marks/" + target + "/lncRNA_Peaks_corrs/lncRNA_Peaks_Correlations_" + str(chunk_index) + '.hdf5', 'w') as f:
        dset = f.create_dataset("default", data=r.T)
        
    print("done chunk " + str(chunk_index))
    return chunk_index


# In[19]:


#По хорошему надо делать свою реализацию ф-ций sf и fdr_bh векторизованные
def testSignificance(data, dof):
    
    def FDR(data):
        return np.array(multipletests(data, alpha=0.05, method='fdr_bh')[0])
    
    t = data * np.sqrt((dof/((data+1.0)*(1.0-data))).clip(0))
    print("count t")
    prob = 2 * t_dist.sf(np.abs(t), dof)
    t = None
    print("count p-values")
    return np.apply_along_axis(FDR, 1, prob)#Boolean mask


# In[24]:


def makeCorrelations(target, chunk_size, threads_num): 
    
    lncRNA = pd.read_csv("../all_marks/"+ target + "/lncRNA_matrix_filtered_norm.csv", sep="\t", index_col=0)
    chip = pd.read_csv("../all_marks/" + target + "/peaks_signal_matrix_norm.csv", sep="\t", index_col=0)
    
    print("Make matrices...")
    X = stats.mstats.rankdata(chip.values, axis=1)
    X = X.astype('float32')
    chip = None
    
    Y = stats.mstats.rankdata(lncRNA.values, axis=1)
    Y = Y.astype('float32')
    
    chunks = np.array_split(Y, int(Y.shape[0]/chunk_size) + 1)
    chunks_pos = np.cumsum([c.shape[0] for c in chunks])
    print(str(len(chunks_pos)) + " chunks, make corrs in " + str(threads_num) + " threads...")
    
    if not os.path.exists("../all_marks/" + target + "/lncRNA_Peaks_corrs/"):
        os.makedirs("../all_marks/" + target + "/lncRNA_Peaks_corrs/")
        
    func = partial(vcorrcoef, target, X)
    
    pool = Pool(processes=threads_num)
    corr = pool.map(func, enumerate(chunks))

    pool.close()
    pool.join()


# In[26]:


def makeCorrelationsMatrix(target, chunk_size):
    i = 0
    dfs = []
    lncRNA = pd.read_csv("../all_marks/"+ target + "/lncRNA_matrix_filtered_norm.csv", sep="\t", index_col=0)
    chip = pd.read_csv("../all_marks/" + target + "/peaks_signal_matrix_norm.csv", sep="\t", index_col=0)
    dof = lncRNA.shape[1] - 2
    chunk_size = 1000
    lncRNAnames = lncRNA.index
    chunks = np.array_split(lncRNA.values, int(lncRNA.shape[0]/chunk_size) + 1)
    chunks_pos = np.cumsum([c.shape[0] for c in chunks])
    
    for i in range(0, len(chunks_pos)):
        with h5py.File("../all_marks/" + target + "/lncRNA_Peaks_corrs/lncRNA_Peaks_Correlations_" + str(i) + '.hdf5', 'r') as f:
            data = f['default'][:]
            names = lncRNAnames[(chunks_pos[i - 1] if i > 0 else 0):chunks_pos[i]]

            signif_mask = testSignificance(data, dof)

            data[~signif_mask] = 0
            signif_mask = None
            zero_ind = np.where(~data.any(axis=1))[0]

            dfs.append(([n for i, n in enumerate(names) if i not in zero_ind], (data[~np.all(data == 0, axis=1)])))
    
    all_corrs = np.concatenate([i[1] for i in dfs])
    all_names = [item for sublist in [i[0] for i in dfs] for item in sublist]
        
    with h5py.File("../all_marks/" + target + "/lncRNA_Peaks_corrs/lncRNA_Peaks_Correlations_corrected_non_zero.hdf5", 'w') as f:
        f.create_dataset("corrs_matrix", data=all_corrs)
        string_dt = h5py.special_dtype(vlen=str)
        f.create_dataset("lncRNAs_names", data=np.array(all_names, dtype='object'), dtype=string_dt)


# In[1]:


targets = [("H3K4me1", "_narrow"), ("H3K4me2", "_narrow"), ("H3K4me3", "_narrow"), ("H3K79me2", ""), 
           ("H3K9ac", "_narrow"), ("H3K9me3", ""), ("H4K20me1", "")]


# In[ ]:


for target in targets:
    chunk_size = 1000
    print("Make correlations for " + target[0] + "...")
    makeCorrelations(target[0], chunk_size, 5)
    makeCorrelationsMatrix(target[0], chunk_size)


# In[ ]:


#========================================== графики ====================================================================


# In[9]:


def getCorrNoZero():
    with open("../H3K27me3/lncRNA_All_Peaks_Correlations_corrected.pickle", 'rb') as f:
        corr = pickle.load(f)
        
    corr_no_zero = [i for i in corr if len(i) != 0]
    corr_no_zero = [item for sublist in corr_no_zero for item in sublist]


# In[33]:


def configPlots():
    sns.set(color_codes=True)
    rcParams['figure.figsize'] = 11.7,8.27
    rcParams["patch.force_edgecolor"] = True


# In[177]:


def corrrelationsHist():
    configPlots()
    corr_no_zero = getCorrNoZero()
    ax = sns.distplot([c[2][0] for c in corr_no_zero], bins=16, kde=False)
    ax.set_yticks(range(0, 2500001, 250000))
    ax.set_xticks(np.arange(-0.7, 0.75, 0.1))
    #ax.set_xticklabels(["{:.1e}".format(x) for x in ax.get_xticks()], rotation=30)
    ax.set(xlabel='correlation')
    plt.show()
    
    #save
    fig = ax.get_figure()
    fig.patch.set_alpha(0)
    fig.savefig("../H3K27me3/plots/MakeLncRNAWithPeaksCorrelations_correlations_hist.png", bbox_inches='tight', pad_inches = 0)


# In[181]:


#barplot +/- correlations
def corrBarplot1():
    ax = sns.barplot(x=['negative correlations', 'positive correlations'], y=[4310512, 3391340])
    ax.set_yticks(range(0, 4500001, 500000))
    plt.show()
    
    fig = ax.get_figure()
    fig.patch.set_alpha(0)
    fig.savefig("../H3K27me3/plots/MakeLncRNAWithPeaksCorrelations_correlations_sign_barplot.png", bbox_inches='tight', pad_inches = 0)


# In[17]:


def peaksCountLog10Hist():
    ax = sns.distplot(np.log10([len(c) for c in corr if len(c) > 0]), bins=5, kde=False)
    #ax.set_xticks(range(0, 80000, 10000))
    ax.set(xlabel='correlation peaks count(count > 0)')
    plt.show()
    
    fig = ax.get_figure()
    fig.patch.set_alpha(0)
    fig.savefig("../H3K27me3/plots/MakeLncRNAWithPeaksCorrelations_peaks_count_hist.png", bbox_inches='tight', pad_inches = 0)


# In[ ]:


#TODO: построить такой график только для пиков в генах и сравнить с фантомом


# In[ ]:


def getAnno():
    return pd.read_csv("../H3K27me3/peaks/anno_merged_peaks.csv", sep="\t", header=None)


# In[35]:


def plotCorrGenesCountLog10Hist():
    df = pd.read_csv("../H3K27me3/peaks/lncRNA_peaks_gene_association.tsv", sep="\t")
    d = {}
    i = 1
    lncRNAs = df['lncRNA'].unique()
    for r in lncRNAs:
        print(str(i))
        d[r] = len(df[df['lncRNA'] == r]['gene'].unique())
        i = i + 1
    
    ax = sns.distplot(np.log10(list(d.values())), bins=8, kde=False)
    ax.set_xticks(np.arange(0, 5, 0.5))
    ax.set(xlabel='log10(correlation genes count)')
    plt.show()
    
    fig = ax.get_figure()
    fig.patch.set_alpha(0)
    fig.savefig("../H3K27me3/plots/MakeLncRNAWithPeaksCorrelations_genes_count_log_10_hist.png", bbox_inches='tight', pad_inches = 0)


# In[ ]:


#Есть ли уникальные пары?(график)


# In[ ]:


#График где выше корреляция, у пиков рядом с генами или нет? в промоторах? или никакой разницы?


# In[ ]:


#Deprecated
def corrForLncRNA(lncRNA):
    print(lncRNA[0] + ", " + str(lncRNA[2]))
    corr = []
    for i, r in chip.iterrows():
        corr.append((lncRNA[0], i, stats.spearmanr(lncRNA[1], r)))
        
    pvalues = [sr[1] for g, p, sr in corr]
    adjustedPvalues = multipletests(pvalues, alpha=0.05, method='fdr_bh')
    corr_corrected_significance = [i for i, j in zip(corr, adjustedPvalues[1]) if j < 0.05]
    
    return corr, corr_corrected_significance

