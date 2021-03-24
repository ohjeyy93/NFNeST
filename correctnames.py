import argparse
import pandas as pd
import subprocess
import csv
from collections import OrderedDict
import json
import os
import matplotlib.pyplot
import pylab
import math

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-d1', dest='directory', type=str, help="name of directory", default=None)
parser.add_argument('-c1', dest='corrected', type=str, help="name of directory", default=None)
args = parser.parse_args() 

incorrected=[]
corrected=[]
df=pd.read_excel(args.corrected)
for index,row in df.iterrows():
    if type(row[0])==str or type(row[2])==str:
        incorrected+=[row[0]]
        corrected+=[row[2]]


for filename in os.listdir(args.directory):
    for x in range(len(incorrected)):
        if filename.replace("_L001_R1_001.fastq.gz","") == incorrected[x]:
            #print(filename)
            old_file = os.path.join(args.directory, filename)
            new_file = os.path.join(args.directory, corrected[x]+"_L001_R1_001.fastq.gz")
            os.rename(old_file,new_file)
        if filename.replace("_L001_R2_001.fastq.gz","") == incorrected[x]:
            #print(filename)
            old_file = os.path.join(args.directory, filename)
            new_file = os.path.join(args.directory, corrected[x]+"_L001_R2_001.fastq.gz")
            os.rename(old_file,new_file)


#matplotlib.pyplot.scatter(wholefilter,wholemajor)
#matplotlib.pyplot.savefig('test.png')