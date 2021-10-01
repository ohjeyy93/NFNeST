import os
import gzip

files = [f for f in os.listdir('.') if os.path.isfile(f)]
for f in files:
    if f.startswith("SRR"):
        #print(f[-3:-1])
        if "_R1" in f:
            #print(f.replace("_1", ""))
            #print(f)
            os.rename(f,f.replace("_R1", "_1"))
            #print(f[0:-6]+"_R1.fastq")
        #if f[-8:-6]==("_2"):
        if "_R2" in f:
            os.rename(f,f.replace("_R2", "_2"))
            #os.rename(f,f.strip("_2"))
            #print(f)
            #os.rename(f,f.strip("_2"))
            #print(f[0:-6]+"_R2.fastq")
        #print(f[0:-6]+"R2")
       
        