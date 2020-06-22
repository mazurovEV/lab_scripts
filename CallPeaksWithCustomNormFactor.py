
# coding: utf-8

# In[35]:


import pandas as pd
import subprocess
import re
import os
import matplotlib.pyplot as plt
import shutil
from multiprocessing import Pool


# In[36]:


#!!TODO: возможно надо делать меньше порог на p-value


# In[37]:


d = pd.read_csv("../all_marks/ChiP-Seq_for_all_marks_all_replicas.tsv", sep='\t', index_col=False)


# In[187]:


#macs2 filterdup -i CTCF_ChIP_200K.bed.gz --keep-dup=1 -o CTCF_ChIP_200K_filterdup.bed
#macs2 filterdup -i CTCF_ChIP_200K_control.bed.gz --keep-dup=1 -o CTCF_ChIP_200K_control_filterdup.bed
#macs2 predictd -i CTCF_ChIP_200K_filterdup.bed -g hs -m 5 50
#macs2 pileup -i CTCF_ChIP_200K_filterdup.bed -o CTCF_ChIP_200K_filterdup.pileup.bdg --extsize 254
#macs2 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 127 -o d_bg.bdg
#macs2 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 500 -o 1k_bg.bdg
#macs2 bdgopt -i 1k_bg.bdg -m multiply -p 0.254 -o 1k_bg_norm.bdg
#macs2 pileup -i CTCF_Control_200K_filterdup.bed -B --extsize 5000 -o 10k_bg.bdg
#macs2 bdgopt -i 1k_bg.bdg -m multiply -p 0.0254 -o 1k_bg_norm.bdg
#macs2 bdgcmp -m max -t 1k_bg_norm.bdg -c 10k_bg_norm.bdg -o 1k_10k_bg_norm.bdg
#macs2 bdgcmp -m max -t 1k_10k_bg_norm.bdg -c d_bg.bdg -o d_1k_10k_bg_norm.bdg
#macs2 bdgopt -i d_1k_10k_bg_norm.bdg -m max -p .0188023 -o local_bias_raw.bdg
#macs2 bdgopt -i local_bias_raw.bdg -m multiply -p .99858 -o local_lambda.bdg
#macs2 bdgcmp -t CTCF_ChIP_200K_filterdup.pileup.bdg -c local_lambda.bdg -m qpois -o CTCF_ChIP_200K_qvalue.bdg
#macs2 bdgbroadcall -i CTCF_ChIP_200K_qvalue.bdg -c 1.301 -l 245 -g 100 -o CTCF_ChIP_200K_peaks.bed
#or& for narrow peaks:
#macs2 bdgpeakcall -i CTCF_ChIP_200K_qvalue.bdg -c 1.301 -l 245 -g 100 -o CTCF_ChIP_200K_peaks.bed
def callPeaks(data):
    signal = data[0]
    signal_path = data[1]
    control = data[2]
    control_path = data[3]
    norm_factor = data[4]
    tmp_dir = "/data/mazurovev/tmp_peaks/tmp_call_peaks_for_" + signal
    os.mkdir(tmp_dir)
    print(signal)
    try:
        output = subprocess.check_output("macs2 filterdup -i " + signal_path + signal + ".bed" + " --keep-dup=1 -o " + tmp_dir + "/" + signal + "_filterdup.bed", shell=True, stderr=subprocess.STDOUT).decode()
        m = re.match(r".+tags after filtering in alignment file: (\d+).+", output, re.MULTILINE|re.DOTALL)
        signal_reads_count = m.group(1)
        
        output = subprocess.check_output("macs2 filterdup -i " + control_path + control + ".bed" + " --keep-dup=1 -o " + tmp_dir + "/" + control + "_filterdup.bed", shell=True, stderr=subprocess.STDOUT).decode()
        m = re.match(r".+tags after filtering in alignment file: (\d+).+", output, re.MULTILINE|re.DOTALL)
        control_reads_count = m.group(1)
        
        output = subprocess.check_output("macs2 predictd -i " + tmp_dir + "/" + signal + "_filterdup.bed" + " -g hs -m 5 50", shell=True, stderr=subprocess.STDOUT).decode()
        m = re.match(r".+predicted fragment length is (\d+).+", output, re.MULTILINE|re.DOTALL)
        extsize = m.group(1)
        
        output = subprocess.check_output("macs2 pileup -i " + tmp_dir + "/" + signal + "_filterdup.bed" + " -f BED -o " + tmp_dir + "/" + signal + "_filterdup.pileup.bdg --extsize " + extsize, shell=True, stderr=subprocess.STDOUT).decode()
        
        output = subprocess.check_output("macs2 pileup -i " + tmp_dir + "/" + control + "_filterdup.bed" + " -B -f BED --extsize " + str(int(extsize)//2) + " -o " + tmp_dir + "/" + control + "_d_bg.bdg", shell=True, stderr=subprocess.STDOUT).decode()
        
        output = subprocess.check_output("macs2 pileup -i " + tmp_dir + "/" + control + "_filterdup.bed" + " -B -f BED --extsize 500 -o " + tmp_dir + "/" + control + "_1k_bg.bdg", shell=True, stderr=subprocess.STDOUT).decode()
        
        output = subprocess.check_output("macs2 bdgopt -i " + tmp_dir + "/" + control + "_1k_bg.bdg" + " -m multiply -p " + str(int(extsize)/1000) + " -o " + tmp_dir + "/" + control + "_1k_bg_norm.bdg", shell=True, stderr=subprocess.STDOUT).decode()
        
        output = subprocess.check_output("macs2 pileup -i " + tmp_dir + "/" + control + "_filterdup.bed" + " -B -f BED --extsize 5000 -o " + tmp_dir + "/" + control + "_10k_bg.bdg", shell=True, stderr=subprocess.STDOUT).decode()
        
        output = subprocess.check_output("macs2 bdgopt -i " + tmp_dir + "/" + control + "_10k_bg.bdg" + " -m multiply -p " + str(int(extsize)/10000) + " -o " + tmp_dir + "/" + control + "_10k_bg_norm.bdg", shell=True, stderr=subprocess.STDOUT).decode()
        
        genome_bg = int(control_reads_count)*int(extsize)/2700000000
        
        output = subprocess.check_output("macs2 bdgcmp -m max -t " + tmp_dir + "/" + control + "_1k_bg_norm.bdg" + " -c " + tmp_dir + "/" + control + "_10k_bg_norm.bdg" + " -o " + tmp_dir + "/" + control + "_1k_10k_bg_norm.bdg", shell=True, stderr=subprocess.STDOUT).decode()
        
        output = subprocess.check_output("macs2 bdgcmp -m max -t " + tmp_dir + "/" + control + "_1k_10k_bg_norm.bdg" + " -c " + tmp_dir + "/" + control + "_d_bg.bdg" + " -o " + tmp_dir + "/" + control + "_d_1k_10k_bg_norm.bdg", shell=True, stderr=subprocess.STDOUT).decode()
        
        output = subprocess.check_output("macs2 bdgopt -i " + tmp_dir + "/" + control + "_d_1k_10k_bg_norm.bdg" + " -m max -p " + str(genome_bg) + " -o " + tmp_dir + "/" + control + "_local_bias_raw.bdg", shell=True, stderr=subprocess.STDOUT).decode()
        
        output = subprocess.check_output("macs2 bdgopt -i " + tmp_dir + "/" + control + "_local_bias_raw.bdg" + " -m multiply -p " + str(norm_factor) + " -o " + tmp_dir + "/" + control + "_local_lambda.bdg", shell=True, stderr=subprocess.STDOUT).decode()
        
        output = subprocess.check_output("macs2 bdgcmp -t " + tmp_dir + "/" + signal + "_filterdup.pileup.bdg" + " -c " + tmp_dir + "/" + control + "_local_lambda.bdg -m qpois" + " -o " + tmp_dir + "/" + signal + "_qvalue.bdg", shell=True, stderr=subprocess.STDOUT).decode()
        
        #broad or narrow peaks
        output = subprocess.check_output("macs2 bdgbroadcall -i " + tmp_dir + "/" + signal + "_qvalue.bdg" + " -c 1.301 " + " -o " + signal_path + signal + "_peaks.bed", shell=True, stderr=subprocess.STDOUT).decode()
        #output = subprocess.check_output("macs2 bdgpeakcall -i " + tmp_dir + "/" + signal + "_qvalue.bdg" + " -c 1.301 -l " + extsize + " -o " + signal_path + signal + "_narrow_peaks.bed", shell=True, stderr=subprocess.STDOUT).decode()
        print("done " + signal)
    except subprocess.CalledProcessError as e:
        output = e.output.decode()
        
    shutil.rmtree(tmp_dir)
    return output


# In[38]:


data = []
for index, row in d.iterrows():
    if(row['Experiment target'] == 'H4K20me1-human'):
        signal_path = "/data/mazurovev/" + row['Experiment target'].split('-')[0] + "/" + row['Biosample term id'] + "/"
        control_path = "/data/mazurovev/controls/"
        data.append((row['File accession'], signal_path, row['Control'], control_path, row['norm_factor']))


# In[39]:


len(data)


# In[ ]:


pool = Pool(processes=10)
outputs = pool.map(callPeaks, data)

pool.close()
pool.join()


# In[21]:


#t = "/data/mazurovev/H3K27ac/UBERON:0001515/ENCFF946NGS.bed"
#c = "/data/mazurovev/controls/ENCFF901HZC.bed"
#subprocess.check_output("macs2 callpeak -t " + t + " -c " + c + " -g hg38 --outdir /data/mazurovev/H3K27ac/UBERON:0001515 -n test 2> /data/mazurovev/H3K27ac/UBERON:0001515/test.log", shell=True, stderr=subprocess.STDOUT).decode()

