
# coding: utf-8

# In[38]:


import nbimporter
import load_data as ld


# In[39]:


from importlib import reload
reload(ld)


# In[27]:


binLength = 450
binStep = 150
dataPath = "../H3K27me3/H3K27me3_filtered_by_biosample.csv"
blackListPath = "../wgEncodeDacMapabilityConsensusExcludable.bed"
binValueBoundary = 5#если значение бина меньше этого - обнуляем(считаем, что это шум)


# In[ ]:


#ld.loadWholeGenome(dataPath, binLength, binStep, "../H3K27me3/WGBins.pickle", False)
ld.getBinsCoords(binLength, binStep, "../H3K27me3/WGCoordsBins.pickle")
ld.loadWholeGenome(dataPath, binLength, binStep, "../H3K27me3/WGControlBins.pickle", True)

