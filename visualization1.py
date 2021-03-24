import argparse
import pandas as pd
import subprocess
import csv
from collections import OrderedDict
import json
import os
import matplotlib.pyplot
import pylab

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-d1', dest='directory', type=str, help="name of directory", default=None)
args = parser.parse_args() 

wholemajor=[]
wholefilter=[]

for filename in os.listdir(args.directory):
    wholelist1=[]
    with open(os.path.join(args.directory, filename), "r") as f1:
        for line in f1:
            templist1=[]
            for word in line.split(","):
                #print(word)
                templist1+=[word.strip("\n")]
            #print(templist1)
            wholelist1+=[templist1]
    numberofmajor=0
    numberoffilter=0
    sumoffilter=0
    averagefilter=0
    for k in range(1,len(wholelist1)):
        #print(wholelist1[x])
        #print(wholelist1[x][18])
        if wholelist1[k][10]=="Major":
            numberofmajor+=1
        if wholelist1[k][15]!="NA":
            numberoffilter+=1
            sumoffilter+=float(wholelist1[k][15])
    if numberoffilter!=0:
        averagefilter=sumoffilter/numberoffilter
        wholemajor+=[numberofmajor]
        wholefilter+=[averagefilter]
        #print(numberofmajor)
        #print(averagefilter)


matplotlib.pyplot.scatter(wholefilter,wholemajor)
matplotlib.pyplot.savefig('test.png')