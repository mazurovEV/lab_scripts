
# coding: utf-8

# In[1]:


import rpy2.robjects as robjects
import pandas as pd
from multiprocessing import Pool
import pickle
import os
from subprocess import Popen, PIPE
import shutil
from functools import partial
import re


# In[24]:


r = robjects.r
r.source("bin_data_for_diagnostics_plot.R")


# In[84]:


def createChroms(path, signal, control, signal_path, control_path):
    #0. Создаем папку
    if not os.path.exists(path):
        os.makedirs(path)
        os.chmod(path, 0o777)
    #Все bed файлы должны быть отсортированны
    print(signal)
    try:
        callWithNiceOutput("for chr in `bedextract --list-chr " + signal_path + signal + "`; do bedextract $chr " + signal_path + signal + "> " + path + signal + ".$chr; done")
        callWithNiceOutput("for chr in `bedextract --list-chr " + control_path + control + "`; do bedextract $chr " + control_path + control + "> " + path + control + ".$chr; done")
        #subprocess.check_output("for chr in `bedextract --list-chr " + signal_path + signal + "`; do bedextract $chr " + signal_path + signal + "> " + path + signal + ".$chr; done", shell=True)
        #subprocess.check_output("for chr in `bedextract --list-chr " + control_path + control + "`; do bedextract $chr " + control_path + control + "> " + path + control + ".$chr; done", shell=True)
        
        print("remove ugly chroms")
        files = [f for f in os.listdir(path) if not re.match(r'(' + signal + '|' + control +')\.chr[M|X|Y|[0-9][0-9]?$', f)]
        for f in files:
            os.remove(path + f) 
    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))


# In[ ]:


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


# In[22]:


#est : the estimated normalizing factor.

#binsize.est : the estimated binsize.

#r.seq.depth: sequencing depth ratio.

#pi0: the estimated proportion of background reads among ChIP sample.

def makeNormalizationFactorAndPlot(signal_path, control_path, signal_name, control_name):
    signal = signal_name + ".bed"
    control = control_name + ".bed"
    path = "/data/mazurovev/norm_factors/f_" + signal_name + "/"
    
    #createChroms(path, signal, control, signal_path, control_path)
    print("start count factor for : " + signal)
    res = r.countFactorAndDiagnosticPlot(path, path, signal + ".chr", control + ".chr",
                                         "/data/mazurovev/norm_factors/plots/" + signal_name)
    #python_res = dict(zip(res.names, map(list,list(res))))
    #python_res['chip'] = data[2]
    
    #with open('/data/mazurovev/norm_factors/plots/' + signal_name + ".txt", 'w') as f:
    #    print(python_res, file=f)
        
    print("done with " + signal_name)
    #shutil.rmtree(path)
    #print("delete " + path)
    #return python_res


# In[86]:


def run(data_path, res_path):
    d = pd.read_csv(data_path, sep='\t')

    data = []
    exist = [f.split('.')[0] for f in os.listdir("/data/mazurovev/norm_factors/plots") if f.split('.')[-1] == 'txt']
    for index, row in d.iterrows():
        target = row['Experiment target'].split('-')[0]
        file_path = '/data/mazurovev/' + target + "/" + row['Biosample term id'] + "/"
        control_path = '/data/mazurovev/controls/'
        if(row['File accession'] not in exist):
            data.append((file_path, control_path, row['File accession'], row['Control']))
 
    pool = Pool(processes=25)
    factors = pool.starmap(makeNormalizationFactorAndPlot, data)
    #print(factors)
    pool.close()
    pool.join()
    
    #with open(res_path, 'wb') as f:
    #    pickle.dump(factors, f)


# In[ ]:


def addNormFactorsToFile(pickle_path, csv_path):
    d = pd.read_csv(csv_path, sep=';', index_col=False)
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f)
        
    factors_col = []
    for index, row in d.iterrows():
        factors_col.append([i['est'] for i in data if i['chip'] == row['File accession']][0][0])
        
    d['norm_factor'] = factors_col
    d.to_csv(csv_path, sep=';', index=False)


# In[ ]:


#Считаем факторы для H3K27me3
#run("../H3K27me3/H3K27me3_filtered_by_biosample.csv", "../H3K27me3/normalizationFactors.pickle", "../H3K27me3/bw/", "../H3K27me3/control/")
#addNormFactorsToFile("../H3K27me3/normalizationFactors.pickle", "../H3K27me3/H3K27me3_filtered_by_biosample.csv")

#Считаем факторы для H3K4me1
run("../ChiP-Seq_for_all_marks_all_replicas.tsv", "/data/mazurovEV/norm_factors/normalizationFactors.pickle")


# In[2]:


exist = [f.split('.')[0] for f in os.listdir("/data/mazurovev/norm_factors/plots") if f.split('.')[-1] == 'txt']


# In[4]:


len(exist)

