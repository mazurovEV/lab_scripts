
# coding: utf-8

# In[1]:


import pandas as pd
import requests
import os
from subprocess import Popen, PIPE
from multiprocessing import Pool
import sys


# In[5]:


chip = pd.read_csv("../all_marks/ChiP-Seq_for_all_marks_all_replicas.tsv", sep = "\t")


# In[59]:


def countDensity(data):
    bam_file_path = data[0]
    bam_file_name = data[1]
    callWithNiceOutput('bedtools bamtobed -i ' + bam_file_path + "/" + bam_file_name + '.bam > ' +  bam_file_path + "/" + bam_file_name + '.bed')
    #print('bedtools bamtobed -i ' + bam_file_path + "/" + bam_file_name + '.bam > ' +  bam_file_path + "/" + bam_file_name + '.bed')
    callWithNiceOutput('sort-bed --max-mem 3G --tmpdir ' + bam_file_path + " " + bam_file_path + "/" + bam_file_name + '.bed > ' + bam_file_path + "/" + bam_file_name + '.sort.bed')
    #print('sort-bed --max-mem 3G --tmpdir ' + bam_file_path + " " + bam_file_path + "/" + bam_file_name + '.bed > ' + bam_file_path + "/" + bam_file_name + '.sort.bed')
    os.remove(bam_file_path + "/" + bam_file_name + ".bam")
    os.remove(bam_file_path + "/" + bam_file_name + ".bed")
    os.rename(bam_file_path + "/" + bam_file_name + ".sort.bed", bam_file_path + "/" + bam_file_name + ".bed")
    print("create " + bam_file_path + "/" + bam_file_name + '.bed sorted file and remove .bam')


# In[17]:


def bed2bam(data):
    bed_file_path = data[0]
    bed_file_name = data[1]
    
    callWithNiceOutput('bedtools bedtobam -i ' + bed_file_path + "/" + bed_file_name + '.bed -g ../hg38/chrom.sizes > ' +  bed_file_path + "/" + bed_file_name + '.bam')


# In[3]:


def callWithNiceOutput(cmd):
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, bufsize=-1)
    
    for line in proc.stdout:
        print("output: ")
        print(line)
        sys.stdout.flush()
    
    for line in proc.stderr:
        print("error: ")
        print(line)
        sys.stdout.flush()


# In[112]:


#Если пришлось останавливать процесс и у нас есть уже bed(не полный или не отсортированный), но не удален bam — удаляем такой bed
def cleanCorruptedBeds():
    chip = pd.read_csv("../ChiP-Seq_for_all_marks_all_replicas.tsv", sep = "\t")
    
    for i, row in chip.iterrows():
        target = row['Experiment target'].split('-')[0]
        file_wo_type = '/data/mazurovev/' + target + "/" + row['Biosample term id'] + "/" + row['File accession']
        if((os.path.exists(file_wo_type + '.bed') and os.path.exists(file_wo_type + '.bam')) or (os.path.exists(file_wo_type + '.bed') and os.path.getsize(file_wo_type + '.bed') < 100 * 1024)):
            os.remove(file_wo_type + '.bed')
            print("remove " + file_wo_type + '.bed')
            
        if(os.path.exists(file_wo_type + 'sort.bed')):
            os.remove(file_wo_type + 'sort.bed')
            print("remove " + file_wo_type + 'sort.bed')
            
        if(os.path.exists(file_wo_type + '.sort.bed')):
            os.remove(file_wo_type + '.sort.bed')
            print("remove " + file_wo_type + '.sort.bed')
            
        if(os.path.exists('/data/mazurovev/controls/' + row['Control'] + 'sort.bed')):
            os.remove('/data/mazurovev/controls/' + row['Control'] + 'sort.bed')
            print("remove " + '/data/mazurovev/controls/' + row['Control'] + 'sort.bed')


# In[ ]:


#cleanCorruptedBeds()


# In[8]:


data = []


# In[130]:


for i, row in chip.iterrows():
    target = row['Experiment target'].split('-')[0]
    dir_path = '/data/mazurovev/' + target + "/" + row['Biosample term id']
    
    if not os.path.exists(dir_path + "/" + row['File accession'] + '.bed'):
        data.append((dir_path, row['File accession']))


# In[9]:


for control_name in chip['Control'].unique():
    if not os.path.exists('/data/mazurovev/controls/' + control_name + '.bam'):
        data.append(('/data/mazurovev/controls', control_name))


# In[18]:


pool = Pool(processes=28)
pool.map(bed2bam, data)
pool.close()
pool.join()

