
# coding: utf-8

# In[55]:


import pyBigWig
import pandas as p
import glob
import numpy as np
from pybedtools import BedTool
import os
from functools import partial
from multiprocessing import Pool
import pickle
import shutil


# In[11]:


def getChroms():
    return {'chrY': 59373566, 'chr1': 249250621, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878,
 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 
 'chr19': 59128983, 'chr2': 243199373, 'chr20': 63025520, 'chr21': 48129895, 'chr22': 51304566, 
 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663,
 'chr8': 146364022, 'chr9': 141213431,'chrX': 155270560, 'chrM': 16571}


# In[7]:


def getRNAseq(biosample):
    return glob.glob("RNA-seq/" + biosample + "_*.tsv")[0]


# In[30]:


def getSamplesSizes(path):
    metaFile = p.read_csv(path, sep=';', index_col=False)
    return metaFile['Reads count'].tolist()


# In[31]:


def getControlSamplesSizes(path):
    metaFile = p.read_csv(path, sep=';', index_col=False)
    return metaFile['Control reads count'].tolist()


# In[21]:


def getListOfBlackZones(chrom):
    blackList = BedTool('../wgEncodeDacMapabilityConsensusExcludable.bed')
    blackListChrom = blackList.filter(lambda b: b.chrom == chrom)
    return [(i.start, i.end) for i in blackListChrom]


# In[25]:


def getBigWigFile(path):
    return pyBigWig.open(path)


# In[25]:


def getBinsCoords(binLength, binStep, pathForSave):
    WGBinsCoords = {}
    for name, length in getChroms().items():
        print("Count coors for  " + name)
        blackZones = getListOfBlackZones(name)
        chromBinsCoords = []
        i = 0
        while i < length:
            j = i + binLength
            if(j > length):
                j = length

            #Если бин включен или пересекает blacklist, то не добавляем его
            #Можем потерять модификации возле стоп зон? Если они и есть, то нужны ли они для нашей задачи?
            inBlacklist = False
            for start, end in blackZones:
                if((j >= start) & (i <= end)):
                    inBlacklist = True
                    break

            if not inBlacklist:
                chromBinsCoords.append((i, j))
                inBlacklist = False

            i = i + binStep
        WGBinsCoords[name] = chromBinsCoords
        
    with open(pathForSave, 'wb') as f:
        pickle.dump(WGBinsCoords, f)


# In[52]:


def loadData(data, chrom, chrLen, binLength, binStep):
    
    blackZones = getListOfBlackZones(chrom)
    f = partial(loadSample, chrom=chrom, chrLen=chrLen, binLength=binLength, binStep=binStep, blackZones=blackZones)
    pool = Pool(processes=32)
    chromData = pool.map(f, data)
    
    pool.close()
    pool.join()
    
    tmp = []
    samplesNames = []
    for sampleName, bins in chromData:
        samplesNames.append(sampleName)
        tmp.append(bins)
    
    filename = "./tmpChromBins/" + chrom + ".pkl"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, 'wb') as f:
        pickle.dump({"data": np.array(tmp), "samplesNames": samplesNames}, f)


# In[45]:


def loadSample(sample, chrom, chrLen, binLength, binStep, blackZones):
    bins = []
    
    bw = getBigWigFile(sample)
        
    if(chrom not in bw.chroms()):#Речь про Y или M хромосому скорее всего
        print("Sample " + sample + " not contain " + chrom + " chromosome")
    else:
        i = 0
        while i < chrLen:
            j = i + binLength
            if(j > chrLen):
                j = chrLen

            #Если бин включен или пересекает blacklist, то не добавляем его
            #Можем потерять модификации возле стоп зон? Если они и есть, то нужны ли они для нашей задачи?
            inBlacklist = False
            for start, end in blackZones:
                if((j >= start) & (i <= end)):
                    inBlacklist = True
                    break

            if not inBlacklist:
                stat = bw.stats(chrom, i, j, exact=True)[0]
                bins.append(stat)
                inBlacklist = False

            i = i + binStep
                
    return (sample.split(os.sep)[-1], bins)


# In[44]:


def loadWholeGenome(path, binLength, binStep, pathForSave, isControlData):
    WGBins = {}
    metaFile = p.read_csv(path, sep=';', index_col=False)
    paths = ['../H3K27me3/control/' + i + '.bw' for i in sorted(metaFile['Control accession'])] if isControlData else ['../H3K27me3/bw/' + i + '.bw' for i in sorted(metaFile['File accession'])]
                        
    for name, length in getChroms().items():
        print("Chrom " + name + ", length = " + str(length))
        loadData(paths, name, length, binLength, binStep)
    
    print("Build one file from chroms bins...")
    
    for name, length in getChroms().items():
        with open("./tmpChromBins/" + name + ".pkl", 'rb') as f:
            data = pickle.load(f)
        f.close()
        print("Add data from " + name)
        WGBins[name] = data
        data = None
   
    with open(pathForSave, 'wb') as f:
        pickle.dump(WGBins, f)
    
    print("Finish!")
    shutil.rmtree("./tmpChromBins/")

