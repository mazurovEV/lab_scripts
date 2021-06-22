
# coding: utf-8

# In[47]:


from bs4 import BeautifulSoup
from collections import defaultdict
import pandas as pd
import sys
import os
import argparse


# In[48]:


parser = argparse.ArgumentParser(description='Process TDF results')


# In[49]:


parser.add_argument('--path', type=str, nargs='+', help='path to the tdf folder(without last slash); script will use index.html from this folder and create data.tsv file in this folder')


# In[50]:


parser.add_argument('--folder', type=str, help='path to the folder with tdf folders(without last slash)')


# In[51]:


args = parser.parse_args()
print(args)

# In[53]:


rnas = ['ENSG00000204054', 'ENSG00000204054_p', 'ENSG00000204054_pp',
        'ENSG00000223478', 'ENSG00000223478_p', 'ENSG00000223478_pp', 
        'ENSG00000231312', 'ENSG00000231312_p', 'ENSG00000231312_pp',
        'ENSG00000233527', 'ENSG00000233527_p', 'ENSG00000233527_pp', 
        'ENSG00000246067', 'ENSG00000246067_p', 'ENSG00000246067_pp']
def custom():
    for rna in rnas:
        with open("../H3K27me3/fantom_tdf/region_test/" + rna + "/index.html") as fp:
            sp = BeautifulSoup(fp, "lxml")

        data = html2frame(sp)
        print("Significant rows: " + str(data[data['p-value'] < 0.05].shape[0]))

        #data.to_csv("../H3K27me3/fantom_tdf/gencode_ann/" + rna + ".tsv", sep="\t", index=None)


# In[61]:


def html2frame(sp):
    table = sp.find_all('tbody')[0]
    d = defaultdict(list)
    for row in table.find_all('tr'):
        columns = row.find_all('td')
        d['DBD'].append(str(columns[1].string))
        d['tr_with DBS'].append(str(columns[2].string))
        d['tr_without DBS'].append(str(columns[3].string))
        d['ntr_with DBS(average)'].append(str(columns[4].string))
        d['ntr_with DBS(std)'].append(str(columns[5].string))
        d['p-value'].append(float(columns[6].string))
        d['z-score'].append(float(columns[7].string))
    
    return pd.DataFrame(d)


# In[64]:


def parseResults(path):
    try:
        with open(path + "/index.html") as fp:
            sp = BeautifulSoup(fp, "lxml")

        data = html2frame(sp)
        print("Significant rows count for " + path + ": " + str(data[data['p-value'] < 0.05].shape[0]))

        data.to_csv(path + "/data.tsv", sep="\t", index=None)
    except EnvironmentError: # parent of IOError, OSError *and* WindowsError where available
        print('Oops, ' + path + '/index.html not found')


# In[ ]:


if(args.path is not None):
    
    for path in args.path:
        parseResults(path)
        
elif(args.folder is not None):
    folders = os.listdir(args.folder)
    
    for path in folders:
        parseResults(args.folder + "/" + path)

