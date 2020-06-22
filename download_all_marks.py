import pandas as pd
import os
import wget
from multiprocessing import Pool
import gzip
import shutil

#chip = pd.read_csv("../ChiP-Seq_for_all_marks_all_replicas.tsv", sep = "\t")
#rna = pd.read_csv("../RNA-Seq_for_all_marks_all_replicas.tsv", sep="\t")
methyl = pd.read_csv("../Methylation_all_replicas.tsv", sep="\t")

data = []
#for index, row in chip.iterrows():
#        target = row['Experiment target'].split('-')[0]
#        dir_path = '/data/mazurovev/' + target + "/" + row['Biosample term id']
#        
#        if not os.path.exists(dir_path + "/" + row['File accession'] + '.bam'):
#            data.append((dir_path, row['File accession'], row['File download URL']))
            
for index, row in methyl.iterrows():
        dir_path = "/data/mazurovev/methylation/" + row['Biosample term id']
        
        if not os.path.exists(dir_path + "/" + row['File accession'] + '.bed.gz'):
            data.append((dir_path, row['File accession'], row['File download URL']))
            
        if not os.path.exists(dir_path + "/" + row['File accession'] + '.bed'):
            with gzip.open(dir_path + '/' + row['File accession'] + ".bed.gz", 'rb') as f_in:
                with open(dir_path + '/' + row['File accession'] + ".bed", 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        
        #os.rename(dir_path + '/' + row['File accession'] + ".bed", dir_path + '/' + row['File accession'] + ".bed.gz")
                    
#for index, row in rna.iterrows():
#        dir_path = '/data/mazurovev/RNA-Seq/' + row['Biosample term id']
#        data.append((dir_path, row['File accession'], row['File download URL']))


#for control_name in chip['Control'].unique():
#   url = 'https://www.encodeproject.org/files/' + control_name + '/@@download/' + control_name + '.bam'
#    
#    if not os.path.exists('/data/mazurovev/controls' + control_name + ".bam"):       
#        data.append(('/data/mazurovev/controls', control_name, url))

def download(data):
    dir_path = data[0]
    file_path = data[1]
    file_url = data[2]
                    
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
                                            
    wget.download(file_url, dir_path + '/' + file_path + ".bed")
    print(file_path + ".bed downloaded")

print(str(len(data)))
p = Pool(processes=30)
      
tmp = p.map(download, data)
p.close()
p.join()
      
