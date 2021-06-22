
# coding: utf-8

# In[18]:


import pickle
import pandas as pd
from pybedtools import BedTool
import matplotlib.pyplot as plt
import seaborn as sb
from collections import defaultdict
import gffutils
from BCBio import GFF
import subprocess
import os
from multiprocessing import Pool


# In[19]:


with open("../H3K27me3/lncRNA_Peaks_Correlations_corrected.pickle", 'rb') as f:
    corr = pickle.load(f)


# In[20]:


peaks = pd.read_csv("../H3K27me3/bw/subsampleHimericBam_peaks.bed", sep="\t", skiprows=1, header=None)


# In[21]:


rnas = set([g for g, p, c in corr])


# In[22]:


gene_key_corrs = defaultdict(list)
for g, p, c in corr:
    gene_key_corrs[g].append((p, c))


# In[23]:


corrs_with_sign = []
for rna in rnas:
    corrs = gene_key_corrs[rna]
    plus_corrs = len([(p, c) for p, c in corrs if c[0] > 0])
    minus_corrs = len([(p, c) for p, c in corrs if c[0] < 0])
    corrs_with_sign.append((rna, plus_corrs, minus_corrs))


# In[24]:


corrs_with_sign_sorted_plus = sorted(corrs_with_sign, key=lambda x: x[1], reverse=True)
corrs_with_sign_sorted_minus = sorted(corrs_with_sign, key=lambda x: x[2], reverse=True)


# In[25]:


corrs = list(set(corrs_with_sign_sorted_plus[0:20]).union(set(corrs_with_sign_sorted_minus[0:20])))


# In[ ]:


def runTDF(corr):
    g = corr[0]
    cmd = "rgt-TDF promotertest -rm 2 -r ../H3K27me3/tdf_lncRNA_gene_sequences/" + g + "_gene_lncRNA.fa -rn " + g + " -bed ../H3K27me3/peaks_for_tdf/plus_" + g + "_peaks.bed -bg ../H3K27me3/peaks_for_tdf/ne_" + g + "_peaks.bed -o ../H3K27me3/promoter_test/ -organism hg19 -l 12 -ccf 20"
    print(cmd)
    return subprocess.check_output(cmd, shell=True) 


# In[16]:


not_done = ["ENSG00000197291.4"]


# In[ ]:


def getConvertIds():
    lncRNAMatrix = pd.read_csv("../H3K27me3/lncRNA_matrix_filtered.csv", sep="\t", index_col=0)

    #hg19 айдишники
    lncRNAGenes = list(lncRNAMatrix.index)

    in_file = "../H3K27me3/gencode.v29lift37.long_noncoding_RNAs.gff3"
    in_handle = open(in_file)

    limit_info = dict(
        gff_id = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 
     'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9'],
        gff_type = ["gene"])

    lncRNAgenes = []
    for rec in GFF.parse(in_handle, limit_info=limit_info):
        lncRNAgenes.append((rec.id, rec.features)) 

    in_handle.close()

    tmp = [k.split('.')[0] for k in lncRNAGenes]

    #38 айдишники
    lncRNAgenesIds = [g.id for i, j in lncRNAgenes for g in j if g.id.split('.')[0] in tmp]

    lncRNAgenesIds.sort(key = lambda x: x.split('.')[0])
    lncRNAGenes.sort(key = lambda x: x.split('.')[0])

    #в gff лежат названия айдишников для версии ghr38, а в данных hg19
    convertIds = {i:j for i, j in zip(lncRNAGenes, lncRNAgenesIds)} 
    
    return convertIds


# In[79]:


def generateFilesForTDF(corrs):
    db = gffutils.FeatureDB('../H3K27me3/long_noncoding_RNAs.db', keep_order=True)
    convertIds = getConvertIds()
    
    for g, p, m in corrs:
        print("lncRNA: " + g)
        corrs = gene_key_corrs[g]
        plus_corrs = [p for p, c in corrs if c[0] > 0]
        minus_corrs = [p for p, c in corrs if c[0] < 0]

        #Сделать TDF положительные против отрицательных
        gene_plus_peaks = peaks[peaks[3].isin(plus_corrs)]
        gene_plus_peaks.to_csv("../H3K27me3/peaks_for_tdf/plus_" + g + "_peaks.bed", sep="\t", index=False, header=None)

        gene_minus_peaks = peaks[peaks[3].isin(plus_corrs)]
        gene_minus_peaks.to_csv("../H3K27me3/peaks_for_tdf/minus_" + g + "_peaks.bed", sep="\t", index=False, header=None)

        ne_gene_peaks = peaks[~peaks[3].isin([p for p, c in corrs])]
        ne_gene_peaks.to_csv("../H3K27me3/peaks_for_tdf/ne_" + g + "_peaks.bed", sep="\t", index=False, header=None)

        #Делаем геномную(!) последовательность нужной lncRNA
        gene = db[convertIds[g]]
        with open("../H3K27me3/tdf_lncRNA_gene_sequences/" + g + ".gff3", "a") as myfile:
            myfile.write(str(gene) + "\n")
            for i in db.children(gene, featuretype='exon', order_by='start'):
                myfile.write(str(i) + "\n")

        cmd = 'gffread ../H3K27me3/tdf_lncRNA_gene_sequences/' + g + '.gff3 -g ../hg19/allChr.fa -w ../H3K27me3/tdf_lncRNA_gene_sequences/' + g + '_gene_lncRNA.fa' 
        print(cmd + "\n")
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        p.wait()
        print(p.returncode)


# In[ ]:


def generate():
    generateFilesForTDF(corrs_with_sign_sorted_plus[0:20])
    generateFilesForTDF(list(set(corrs_with_sign_sorted_minus[0:20]) - set(corrs_with_sign_sorted_plus[0:20])))


# In[ ]:


#pool = Pool(processes=1)
#r = pool.map(runTDF, not_done)

#pool.close()
#pool.join()
runTDF(("ENSG00000197291.4", 1, 1))

