
# coding: utf-8

# In[23]:


import pandas as pd
from scipy import stats
import pickle
#%matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
from matplotlib import rcParams
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool
import datetime


# In[9]:


anno = pd.read_csv("../all_marks/H3K27ac/peaks_anno.csv", sep="\t")


# In[3]:


chip = pd.read_csv("../all_marks/H3K27ac/peaks_signal_matrix_norm.csv", sep="\t", index_col=0)


# In[5]:


rna = pd.read_csv("../all_marks/H3K27ac/rna_matrix_norm.csv", sep="\t")


# In[22]:


print("create data... " + str(datetime.datetime.now()))
data = []
for index, row in anno.iterrows():
    gene = rna.filter(like=row['feature'], axis=0)
    if(gene.shape[0] != 0):
        data.append((chip.loc[row['peak']].tolist(), gene.iloc[0].tolist()))


# In[ ]:


def makeCorrs(data):
    corrs = []
    for d in data:
        corrs.append(stats.spearmanr(d[0], d[1]))


# In[ ]:


print("Make corrs... " + str(datetime.datetime.now()))
pool = Pool(processes=20)
corr = pool.map(makeCorrs, [data[i:i + 35000] for i in range(0, len(data), 35000)])

pool.close()
pool.join()


# In[46]:


with open("../all_marks/H3K27ac/check_peaks_with_corr.pickle", 'wb') as f:
        pickle.dump(corr, f)


# In[101]:


#with open("../H3K27me3/corrs.pickle", 'rb') as f:
#    corr = pickle.load(f)


# In[102]:


#adjustedPvalues = multipletests([j for i, j in corr if ~np.isnan(i)], alpha=0.01, method='fdr_bh')


# In[104]:


#corr1 = [i for i, j in zip([k for k, v in corr if ~np.isnan(k)], adjustedPvalues[1]) if j < 0.01]


# In[107]:


#rcParams['figure.figsize'] = 11.7,8.27
#rcParams["patch.force_edgecolor"] = True
#sb.set(color_codes=True)


# In[109]:


#ax = sb.distplot(corr1, bins=5, kde=False)
#ax.set_xticks(np.arange(-1, 1, 0.16))
#ax.set_yticks(range(0, 5000, 500))
#ax.set_xticklabels(np.around(ax.get_xticks(), decimals=2), rotation=35)
#ax.set(xlabel='correlation')
#plt.show()


# In[110]:


#fig = ax.get_figure()
#fig.patch.set_alpha(0)
#fig.savefig("../H3K27me3/plots/MakeCorrespondenceBetweenGenesAndPeaks_correlation_hist.png", bbox_inches='tight', pad_inches = 0)

