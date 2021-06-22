
# coding: utf-8

# In[4]:


import pandas as pd
import requests
import os
from subprocess import Popen, PIPE
from multiprocessing import Pool
import sys


# In[7]:


os.chdir("/home/mazurovev/H3K27me3/chIP-Seq_control_bams")


# In[8]:


chipData = pd.read_csv("../H3K27me3_filtered_by_biosample.csv", sep=';')


# In[9]:


def countDensity(bam_file_name):
    callWithNiceOutput('../../IGVTools/igvtools sort ' + bam_file_name + '.bam ' + bam_file_name + '_sort.bam')
    os.remove(bam_file_name + ".bam")
    callWithNiceOutput('../../IGVTools/igvtools index ' + bam_file_name + '_sort.bam')
    callWithNiceOutput('../../IGVTools/igvtools count -e 200 ' + bam_file_name + '_sort.bam ' + bam_file_name + '.wig ../../IGVTools/genomes/hg19.chrom.sizes')
    os.remove(bam_file_name + '_sort.bam')
    os.remove(bam_file_name + '_sort.bam.bai')
    callWithNiceOutput('wigToBigWig ' + bam_file_name + '.wig ../../IGVTools/genomes/hg19.chrom.sizes ' + bam_file_name + '.bw')
    os.remove(bam_file_name + ".wig")


# In[11]:


def callWithNiceOutput(cmd):
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, bufsize=-1)
    
    for line in proc.stdout:
        print(line)
        sys.stdout.flush()


# In[ ]:


pool = Pool(processes=8)
pool.map(countDensity, chipData['Control accession'].tolist())

