
# coding: utf-8

# In[1]:


#1)Вызвать для всех bam'ов контроля bedtools genomecov с превращением в bedGraph со скейлом(norm_factor)
#2)Вызвать для всех bam'ов сигнала bedtools genomecov с превращением в bedGraph
#3)Используя unionbedg объединить 2 верхних bedgraph в один
#4)Вычесть из колонки сигнала колонку контроля, схлопнуть позиции(если надо будет) - получить новый bedGraph - "очищенный" сигнал,
#годный для дальнейшей нормализации(подсчет  плотности в "генах" и последующая нормализация DESeq2)


# In[1]:


import pandas as pd
import os
from pybedtools import BedTool
from multiprocessing import Pool


# In[18]:


d = pd.read_csv("../H3K27me3/H3K27me3_filtered_by_biosample.csv", sep=';', index_col=False)


# In[19]:


data = []
for index, row in d.iterrows():
    data.append((row['File accession'], row['Control accession'], row['norm_factor']))


# In[21]:


def makePureSignal(row):
    signal = row[0]
    control = row[1]
    norm_factor = row[2]
    
    bam = BedTool('../H3K27me3/bw/' + signal + '.bam')
    bam.genome_coverage(bg=True, output="../H3K27me3/bw/" + signal + ".bdg")
    
    bam_control = BedTool('../H3K27me3/control/' + control + '.bam')
    bam_control.genome_coverage(bg=True, scale=norm_factor, output="../H3K27me3/control/" + control + ".bdg")
    
    s = BedTool('../H3K27me3/bw/' + signal + '.bdg')
    c = BedTool('../H3K27me3/control/' + control + '.bdg')
    BedTool().union_bedgraphs(i=[s.fn, c.fn], output="../H3K27me3/bw/" + signal + "_" + control + ".bdg")
    
    union = pd.read_csv("../H3K27me3/bw/" + signal + "_" + control + ".bdg", sep='\t', index_col=False, header=None)
    diff = union[3] - union[4]
    diff = diff.apply(lambda x: 0 if x < 0 else x)
    union["diff"] = diff
    union.drop([3, 4], axis=1, inplace=True)
    without_zero = union[union['diff'] != 0]
    os.remove("../H3K27me3/bw/" + signal + "_" + control + ".bdg")
    without_zero.to_csv("../H3K27me3/bw/pure_" + signal + ".bdg", sep='\t') 


# In[ ]:


pool = Pool(processes=16)
r = pool.map(makePureSignal, data)

pool.close()
pool.join()

