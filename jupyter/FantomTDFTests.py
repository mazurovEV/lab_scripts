
# coding: utf-8

# In[139]:


import pandas as pd
import subprocess
from Bio import SeqIO
import pickle
from pybedtools import BedTool
from multiprocessing import Pool
from itertools import chain


# In[24]:


rnas_pp = ['ENSG00000231312']
rnas_mp = ['ENSG00000214548', 'ENSG00000246273']
rnas_mm = ['ENSG00000231312', 'ENSG00000246273', 'ENSG00000260032']


# In[140]:


rnas_with_prefixes = [('ENSG00000231312', '', ''), ('ENSG00000231312', '', '_p'), ('ENSG00000231312', '_p', 'p'),
                     ('ENSG00000231312', '', '_m'), ('ENSG00000231312', '_m', 'm'),
                     ('ENSG00000214548', '', ''), ('ENSG00000214548', '', '_p'), ('ENSG00000214548', '_m', 'p'),
                     ('ENSG00000246273', '', ''), ('ENSG00000246273', '', '_p'), ('ENSG00000246273', '_m', 'p'),
                     ('ENSG00000246273', '', '_m'), ('ENSG00000246273', '_m', 'm'),
                     ('ENSG00000260032', '', ''), ('ENSG00000260032', '', '_m'), ('ENSG00000260032', '_m', 'm')]


# In[ ]:


rnas_with_prefixes_recalc = [('ENSG00000214548', '', '_p'), ('ENSG00000214548', '_m', 'p'), 
                             ('ENSG00000246273', '_m', 'p'), ('ENSG00000246273', '', '_p')]


# In[53]:


rnas_with_prefixes_recalc_another_ccf = [('ENSG00000214548', '_m', 'p')]


# In[45]:


def generateTranscriptFiles(rnas):
    fasta_sequences = SeqIO.parse(open("/home/ymedvedeva/FANTOM6/data/private/F6_CAT.DMFB.target.all_transcript.fa"), 'fasta')
    for fasta in fasta_sequences:
        if(fasta.id.split('|')[1] in ['ENSG00000214548', 'ENSG00000246273', 'ENSG00000260032']):
            print(fasta.id.split('|')[1])
            with open("../fantom6/transcripts/" + fasta.id.split('|')[1] + ".fa", "w") as output_handle:
                SeqIO.write(fasta, output_handle, "fasta")


# In[ ]:


#generateFilesForTests(['ENSG00000231312', 'ENSG00000214548', 'ENSG00000246273', 'ENSG00000260032'])


# In[50]:


#generate proms and background files for promoter and region tests
def generateFilesForTests(rnas):
    fantom = pd.read_csv("../fantom6/oligo_DE_Summary_gene_filtered.tsv", sep="\t")
    
    #В качестве бэкграунда берем все деги(а не просто все фантомовские гены)
    fantom_all_genes = set(fantom['geneID'].unique())
    
    #Уже проверено, что гены из общего списка
    with open("../H3K27me3/our_fantom_genes_association.pickle", 'rb') as f:
        our_fantom_genes = pickle.load(f)
        
    fantom_ann = pd.read_csv("/home/ymedvedeva/FANTOM6/data/public/FANTOM5-CAT/F6_CAT.gene.info.tsv", sep="\t")
    fantom_ann = fantom_ann[['cntg', 'geneStart', 'geneEnd', 'geneID', 'geneName', 'strnd']]
    
    for rna in rnas:
        #Все деги
        genes = set(our_fantom_genes[rna]['fantom'])
        #Все деги c положительным fc
        genes_p = set(our_fantom_genes[rna]['fantom_plus'])
        #Все деги c положительным fc
        genes_m = set(our_fantom_genes[rna]['fantom_minus'])
        #Все деги c положительным fc и положительной корреляцией
        genes_pp = set(our_fantom_genes[rna]['fantom_plus']).intersection(set(our_fantom_genes[rna]['our_plus']))
        #Все деги c положительным fc и отрицательной корреляцией
        genes_pm = set(our_fantom_genes[rna]['fantom_minus']).intersection(set(our_fantom_genes[rna]['our_plus']))
        #Все деги c отрицательным fc и положительной корреляцией
        genes_mp = set(our_fantom_genes[rna]['fantom_plus']).intersection(set(our_fantom_genes[rna]['our_minus']))
        #Все деги c отрицательным fc и отрицательной корреляцией
        genes_mm = set(our_fantom_genes[rna]['fantom_minus']).intersection(set(our_fantom_genes[rna]['our_minus']))

        background_genes = fantom_all_genes - genes
        background_genes_p = fantom_all_genes - genes_p
        background_genes_m = fantom_all_genes - genes_m
        background_genes_pp = fantom_all_genes - genes_pp
        background_genes_pm = fantom_all_genes - genes_pm
        background_genes_mp = fantom_all_genes - genes_mp
        background_genes_mm = fantom_all_genes - genes_mm

        genes_bed = fantom_ann[fantom_ann['geneID'].isin(genes)]
        background_genes_bed = fantom_ann[fantom_ann['geneID'].isin(background_genes)]

        genes_bed_p = fantom_ann[fantom_ann['geneID'].isin(genes_p)]
        background_genes_bed_p = fantom_ann[fantom_ann['geneID'].isin(background_genes_p)]
        
        genes_bed_m = fantom_ann[fantom_ann['geneID'].isin(genes_m)]
        background_genes_bed_m = fantom_ann[fantom_ann['geneID'].isin(background_genes_m)]

        genes_bed_pp = fantom_ann[fantom_ann['geneID'].isin(genes_pp)]
        background_genes_bed_pp = fantom_ann[fantom_ann['geneID'].isin(background_genes_pp)]
        
        genes_bed_pm = fantom_ann[fantom_ann['geneID'].isin(genes_pm)]
        background_genes_bed_pm = fantom_ann[fantom_ann['geneID'].isin(background_genes_pm)]
        
        genes_bed_mp = fantom_ann[fantom_ann['geneID'].isin(genes_mp)]
        background_genes_bed_mp = fantom_ann[fantom_ann['geneID'].isin(background_genes_mp)]
        
        genes_bed_mm = fantom_ann[fantom_ann['geneID'].isin(genes_mm)]
        background_genes_bed_mm = fantom_ann[fantom_ann['geneID'].isin(background_genes_mm)]

        proms = getProms(genes_bed)
        background_proms = getProms(background_genes_bed)

        proms.to_csv("../H3K27me3/fantom_tdf/" + rna + "_genes_proms.bed", sep="\t", header=None, index=None)
        background_proms.to_csv("../H3K27me3/fantom_tdf/" + rna + "_background_genes_proms.bed", sep="\t", header=None, index=None)

        proms_p = getProms(genes_bed_p)
        background_proms_p = getProms(background_genes_bed_p)

        proms_p.to_csv("../H3K27me3/fantom_tdf/" + rna + "_genes_p_proms.bed", sep="\t", header=None, index=None)
        background_proms_p.to_csv("../H3K27me3/fantom_tdf/" + rna + "_background_genes_p_proms.bed", sep="\t", header=None, index=None)
        
        proms_m = getProms(genes_bed_m)
        background_proms_m = getProms(background_genes_bed_m)

        proms_m.to_csv("../H3K27me3/fantom_tdf/" + rna + "_genes_m_proms.bed", sep="\t", header=None, index=None)
        background_proms_m.to_csv("../H3K27me3/fantom_tdf/" + rna + "_background_genes_m_proms.bed", sep="\t", header=None, index=None)

        proms_pp = getProms(genes_bed_pp)
        background_proms_pp = getProms(background_genes_bed_pp)

        proms_pp.to_csv("../H3K27me3/fantom_tdf/" + rna + "_genes_pp_proms.bed", sep="\t", header=None, index=None)
        background_proms_pp.to_csv("../H3K27me3/fantom_tdf/" + rna + "_background_genes_pp_proms.bed", sep="\t", header=None, index=None)
        
        proms_pm = getProms(genes_bed_pm)
        background_proms_pm = getProms(background_genes_bed_pm)

        proms_pm.to_csv("../H3K27me3/fantom_tdf/" + rna + "_genes_pm_proms.bed", sep="\t", header=None, index=None)
        background_proms_pm.to_csv("../H3K27me3/fantom_tdf/" + rna + "_background_genes_pm_proms.bed", sep="\t", header=None, index=None)
        
        proms_mp = getProms(genes_bed_mp)
        background_proms_mp = getProms(background_genes_bed_mp)

        proms_mp.to_csv("../H3K27me3/fantom_tdf/" + rna + "_genes_mp_proms.bed", sep="\t", header=None, index=None)
        background_proms_mp.to_csv("../H3K27me3/fantom_tdf/" + rna + "_background_genes_mp_proms.bed", sep="\t", header=None, index=None)
        
        proms_mm = getProms(genes_bed_mm)
        background_proms_mm = getProms(background_genes_bed_mm)

        proms_mm.to_csv("../H3K27me3/fantom_tdf/" + rna + "_genes_mm_proms.bed", sep="\t", header=None, index=None)
        background_proms_mm.to_csv("../H3K27me3/fantom_tdf/" + rna + "_background_genes_mm_proms.bed", sep="\t", header=None, index=None)

        #for region test
        generateBackgroundForRegionTest(rna, "../H3K27me3/fantom_tdf/" + rna + "_background_genes_proms.bed")
        generateBackgroundForRegionTest(rna + "_p", "../H3K27me3/fantom_tdf/" + rna + "_background_genes_p_proms.bed")
        generateBackgroundForRegionTest(rna + "_m", "../H3K27me3/fantom_tdf/" + rna + "_background_genes_m_proms.bed")
        generateBackgroundForRegionTest(rna + "_pp", "../H3K27me3/fantom_tdf/" + rna + "_background_genes_pp_proms.bed")
        generateBackgroundForRegionTest(rna + "_pm", "../H3K27me3/fantom_tdf/" + rna + "_background_genes_pm_proms.bed")
        generateBackgroundForRegionTest(rna + "_mp", "../H3K27me3/fantom_tdf/" + rna + "_background_genes_mp_proms.bed")
        generateBackgroundForRegionTest(rna + "_mm", "../H3K27me3/fantom_tdf/" + rna + "_background_genes_mm_proms.bed")


# In[4]:


def getProms(bed):
    bed['startProm'] = bed['geneStart'] - 500
    bed['startProm'] = bed['startProm'].apply(lambda x: 0 if x < 0 else x)
    bed['endProm'] = bed['geneStart'] + 500
    
    return bed[['cntg', 'startProm', 'endProm', 'geneID', 'geneName', 'strnd']]


# In[24]:


def runTDFpromoterTest(rna, prefix):
    cmd = "rgt-TDF promotertest -rm 2 -r ../fantom6/transcripts/" + rna + ".fa -rn " + rna + prefix + " -bed ../H3K27me3/fantom_tdf/" + rna + "_genes" + prefix + "_proms.bed -bg ../H3K27me3/fantom_tdf/" + rna + "_background_genes" + prefix + "_proms.bed -o ../H3K27me3/fantom_tdf/promoter_test/ -organism hg19 -l 12 -ccf 20"
    print(cmd)
    return subprocess.check_output(cmd, shell=True)


# In[141]:


def runTDFregionTest(rna, prefix_our, prefix_fantom):
    cmd = "rgt-TDF regiontest -r ../fantom6/transcripts/" + rna + ".fa -rn " + rna + prefix_our + prefix_fantom + " -bed ../H3K27me3/fantom_tdf/" + rna + "_genes" + prefix_our + prefix_fantom + "_proms.bed -f ../H3K27me3/fantom_tdf/bg_" + rna + prefix_our + prefix_fantom + ".bed -o ../H3K27me3/fantom_tdf/region_test_2/ -organism hg38 -l 12 -e 10 -obed -rt -mp 6 -ccf 50 -n" + str(1000)
    print(cmd)
    return subprocess.check_output(cmd, shell=True)


# In[5]:


#bedtools subtract -a A.bed -b B.bed
#A - B = substract
#геном - промоторы нетаргетных дегов, тогда в качестве бэкграунда будут браться только промоторы нетаргетных дегов(т.е. то, что вычли)
def generateBackgroundForRegionTest(rna, nontarget_promoters):
    target = BedTool(nontarget_promoters)
    hg19 = BedTool('../hg38/allChr.bed')
    hg19.subtract(target, output='../H3K27me3/fantom_tdf/bg_' + rna + ".bed")


# In[42]:


def runTDF(rna_prefix):
    print(runTDFregionTest(rna_prefix[0], rna_prefix[1], rna_prefix[2]))


# In[ ]:


pool = Pool(processes=5)
res = pool.map(runTDF, rnas_with_prefixes)

pool.close()
pool.join()


# In[107]:


def tdfDataToHtml():
    import os
    folders = os.listdir("../H3K27me3/fantom_tdf/region_test")

    for path in folders:
        print(path)
        try:
            #os.chmod("../H3K27me3/fantom_tdf/region_test/" + path + "/data.tsv", 0o777)
            getHtmlTableWithTDFResults("../H3K27me3/fantom_tdf/region_test/" + path)
        except FileNotFoundError: 
            print('Oops, ' + path + '/data.tsv not found')


# In[137]:


def getHtmlTableWithTDFResults(path):
    def color_negative_red(row):
        val = row[5]
        try:
            color = 'red' if val < 0.05 else 'black'
        except Exception:
            color = 'black'
        return ['color: black']*5 + ['color: %s' % color] + ['color: black']
    
    def processTable(df):
        new = df['DBD'].str.split("-", n = 1, expand = True) 
        new[0] = pd.to_numeric(new[0])
        new[1] = pd.to_numeric(new[1])
        
        df['startDBD'] = new[0]
        df['endDBD'] = new[1]
        
        df.sort_values(by=['startDBD'], ascending=True, inplace=True)
        df.drop(columns=['startDBD', 'endDBD'], inplace=True)
        df = df[df['p-value'] < 0.05]
        df = df.reset_index(drop=True)
        
        return df
    
    df = processTable(pd.read_csv(path + "/data.tsv", sep="\t"))
    #s = df.style.apply(color_negative_red, axis=1)
    s = df.style
    
    with open(path + "/data_table.html", 'w') as f:
        for index, item in enumerate(s.render().split("\n")):
            if(index == 0):
                f.write('<style  type="text/css" >\n')
                f.write('\n')
                f.write('table, th, td {\n')
                f.write('border: 1px solid black;\n')
                f.write('border-collapse: collapse;\n')
                f.write('font-size: 11px;\n')
                f.write('}\n')
                f.write('\n')
                f.write('th,\n')
                f.write('td {\n')
                f.write('border: 1px solid black;\n')
                f.write('width: 100px;\n')
                f.write('height: 25px;\n')
                f.write('text-align:center;\n')
                f.write('font-family: Montserrat;\n')
                f.write('overflow: hidden;\n')
                f.write('}\n')
                f.write('\n')
                f.write('tr:nth-child(even) {\n')
                f.write('background-color: #ffe6e6\n')
                f.write('}\n')
                f.write('\n')
            else:
                f.write("%s\n" % item)


# In[135]:


def processTable(df):
        new = df['DBD'].str.split("-", n = 1, expand = True) 
        new[0] = pd.to_numeric(new[0])
        new[1] = pd.to_numeric(new[1])
        
        df['startDBD'] = new[0]
        df['endDBD'] = new[1]
        
        df.sort_values(by=['startDBD'], ascending=True, inplace=True)
        df.drop(columns=['startDBD', 'endDBD'], inplace=True)
        df = df[df['p-value'] < 0.05]
        df = df.reset_index(drop=True)
        
        return df

