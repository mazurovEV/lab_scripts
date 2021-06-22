
# coding: utf-8

# In[22]:


#TODO: решить вопрос, должны ли экспериментаторские файлы быть "нормализованы"(нужно ли вычитать контроль?) - иначе
#мы находим пики там, где они не были найдены(Типа MACS не нашел пиков)!
#Вместо предыдущей свистопляски просто считаем кол-во ридов в какждом пике
#(никаких вычитаний контроля и прочего, мы использовали контроль, когда считали пики в каждом образце)
#bedtools coverage -a A.bed -b B.bed, где A это пики, а B это эксперимент
import pandas as pd
import os
from pybedtools import BedTool
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import pickle
from multiprocessing import Pool


# In[23]:


targets = [("H3K4me1", "_narrow"), ("H3K4me2", "_narrow"), ("H3K4me3", "_narrow"), ("H3K79me2", ""), 
           ("H3K9ac", "_narrow"), ("H3K9me3", ""), ("H4K20me1", "")]


# In[27]:


#Список путей к первой реплкие каждого биосемпла
def getFiles(target):
    files = []
    for dirpath, dirnames, filenames in os.walk("/data/mazurovev/" + target[0] + "/"):
        for filename in sorted([f for f in filenames if not f.endswith("peaks.bed")])[:1]:
            #Делаем путем к файлу с пиками
            files.append(dirpath + "/" + filename)
    return sorted(files)


# In[28]:


def makeCoverage(bed_path):
    exp = BedTool(bed_path)
    print("/".join(bed_path.split("/")[-2:]))
    return ("/".join(bed_path.split("/")[-2:]), peaks.coverage(a=peaks.fn, b=exp.fn).to_dataframe()[['score', 'name']])


# In[29]:


def makeSignalMatrix(target):
    
    pool = Pool(processes=3)
    data = pool.map(makeCoverage, getFiles(target))

    pool.close()
    pool.join()
    
    for_data_frame = {}
    for n, d in data:
        for_data_frame[n] = d['score']
        
    matrix = pd.DataFrame(for_data_frame)
    
    matrix.index = ["peak_" + str(i) for i in matrix.index]
    
    matrix.to_csv("../all_marks/" + target[0] + "/peaks_signal_matrix.csv", sep="\t")


# In[10]:


for target in targets:
    print("Make Coverage in Peaks for " + target[0] + "...")
    peaks = BedTool("../all_marks/" + target[0] + "/merged_peaks_first_in_biosample.bed")
    df = makeSignalMatrix(target)

