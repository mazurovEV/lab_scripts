
# coding: utf-8

# In[2]:


import pandas as pd
import rpy2.robjects as robjects
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams


# In[11]:


sns.set(color_codes=True)
rcParams['figure.figsize'] = 11.7,8.27
rcParams["patch.force_edgecolor"] = True


# In[29]:


r = robjects.r
r.source("PeakAnnotation_R.r")


# In[ ]:


r.peaksAnno("../all_marks/H3K27ac/merged_peaks_first_in_biosample.bed", "../all_marks/H3K27ac/peaks_anno.csv")


# In[16]:


anno = pd.read_csv("../all_marks/H3K27ac/peaks_anno.csv", sep="\t")


# In[17]:


anno.head()


# In[10]:


anno = anno[anno['insideFeature'].notnull()]
anno = anno[anno.feature.str.startswith('ENSG')]


# In[11]:


anno.to_csv("../all_marks/H3K27me3/peaks_anno.csv", sep="\t", index=None)


# In[13]:


pos = anno['insideFeature'].value_counts()
ax = sns.barplot(pos.index, pos.values)
ax.set(xlabel='positions', ylabel='count')
#ax.set_yticks(range(0, 700001, 30000))
plt.show()

