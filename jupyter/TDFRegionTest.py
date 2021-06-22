
# coding: utf-8

# In[2]:


import subprocess
from pybedtools import BedTool


# In[ ]:


rna_for_region_test = ["ENSG00000185044.10", "ENSG00000204055.4", "ENSG00000206195.6", "ENSG00000257176.1"]


# In[3]:


def runTDFregionTest(rna):
    cmd = "rgt-TDF regiontest -rm 2 -r ../H3K27me3/tdf_lncRNA_gene_sequences/" + rna + "_gene_lncRNA.fa -rn " + rna + " -bed ../H3K27me3/peaks_for_tdf/plus_" + rna + "_peaks.bed -f ../H3K27me3/region_test/bg_" + rna + ".bed -o ../H3K27me3/region_test/ -organism hg19 -l 12 -ccf 50 -n " + str(1000)
    print(cmd)
    return subprocess.check_output(cmd, shell=True)


# In[ ]:


#bedtools subtract -a A.bed -b B.bed
#A - B = substract
def generateBackgroundForRegionTest(rna):
    target = BedTool('../H3K27me3/peaks_for_tdf/plus_' + rna + '_peaks.bed')
    hg19 = BedTool('../hg19/allChr.bed')
    hg19.subtract(target, output='../H3K27me3/region_test/bg_' + rna + ".bed")


# In[ ]:


for rna in rna_for_region_test:
    runTDFregionTest(rna)

