
# coding: utf-8

# In[161]:


import pandas as pd
import pickle
from pybedtools import BedTool
from multiprocessing import Pool
import pickle


# In[157]:


def countMatrix(f):
    a = BedTool("../H3K27me3/bw/subsampleHimericBam_peaks.bed")
    b = BedTool("../H3K27me3/bw/sorted_t_pure_" + f + ".bdg")
    intersect = a.intersect(b, wao=True, sorted=True)
    return {f: countPeakSums(intersect)}


# In[158]:


def countPeakSums(intersect):
    signalSum = {}
    peak = ""
    for row in intersect:
        if(peak != row[3]):
            if(peak != ""): signalSum[peak] = int(partPeakSum)
            peak = row[3]
            partPeakSum = 0

        partPeakSum = partPeakSum + float(row[18])*int(row[19]) if int(row[19]) != 0 else 0
    
    return signalSum


# In[159]:


d = pd.read_csv("../H3K27me3/H3K27me3_filtered_by_biosample.csv", sep=';', index_col=False)
files = d['File accession'].tolist()


# In[ ]:


pool = Pool(processes=12)
matrix = pool.map(countMatrix, files)

pool.close()
pool.join()


# In[ ]:


with open("../H3K27me3/pre_signal_matrix.pickle", 'wb') as f:
        pickle.dump(matrix, f)

