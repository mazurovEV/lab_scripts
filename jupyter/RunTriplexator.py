
# coding: utf-8

# In[1]:


import pandas as pd
import subprocess
import pickle


# In[2]:


with open("../H3K27me3/lncRNA_Peaks_Correlations_corrected.pickle", 'rb') as f:
    corr = pickle.load(f)


# In[4]:


corr_sorted = sorted(corr, key=lambda x: x[2][0])


# In[18]:


minus_corr_top_20_lncRNA = set([g for g, p, c in corr_sorted[0:20]])


# In[19]:


plus_corr_top_20_lncRNA = set([g for g, p, c in corr_sorted[-20:]])


# In[36]:


def run(lncRNA, output_folder, filename):
    dataframes = []
    noTriplexList = []
    i = 0
    for l in lncRNA:
        print('start ' + str(i) + " lncRNA = " + l)
        result = runTriplexator(l, output_folder)
        if(len(result.index) != 0): 
            dataframes.append(result)
        else:
            noTriplexList.append(l)
        i = i + 1

    data = pd.concat(dataframes)
    data = data.reset_index(drop=True)
    
    data.to_csv('../H3K27me3/' + filename)


# In[34]:


#Почему -m R ?
def runTriplexator(l, output):
    cmd = 'triplexator -l 15 -e 20 -c 2 -fr off -g 20 -fm 0 -of 1 -od ' + output + '  -o ' + l + '.tsv -po -rm 2 -p 3 -ss ../H3K27me3/lnc_gff/' + l + '_transcript_lncRNA.fa -ds ../H3K27me3/peaks_coords/seqs_for_lnc_' + l + '.fa'     
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    p.wait()
    print(l + ' ' + str(p.returncode))
    return parseTriplexatorResult(l, 'ENST')


# In[63]:


def parseTriplexatorResult(l, filterName):
    data = pd.read_csv("../H3K27me3/triplexator_output/" + l + '.tsv', sep='\t')
    data = data[data['# Sequence-ID'].str.startswith(filterName)]
    if(len(data.index) == 0): return pd.DataFrame()
    dataGroup = data.groupby(['# Sequence-ID', 'Duplex-ID'])[["Score"]].sum()
    
    x = pd.DataFrame({'lncRNA':[l]*len(dataGroup.index.tolist()), 
                      'Transcript': [i for i, j in dataGroup.index.tolist()],
                  'Bin' : [j for i, j in dataGroup.index.tolist()], 
                 'Score' : dataGroup['Score'].tolist()})
    x = x.reindex_axis(['lncRNA','Transcript', 'Bin', 'Score'], axis=1)
    x.to_csv(l + ".csv")
    
    return x


# In[ ]:


run(minus_corr_top_20_lncRNA, "../H3K27me3/triplexator_output_minus_top/", "minus_top_lncRNA_triplexes.csv")


# In[ ]:


run(plus_corr_top_20_lncRNA, "../H3K27me3/triplexator_output_plus_top/", "plus_top_lncRNA_triplexes.csv")

