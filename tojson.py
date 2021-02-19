import argparse
import pandas as pd
import subprocess
import csv
from collections import OrderedDict
import json

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-s1', dest='study', type=str, help="name of study", default="xxxx")
parser.add_argument('-y1', dest='year', type=str, help="year", default="xxxx")
parser.add_argument('-n1', dest='numsam', type=str, help="name of sample", default="xxxx")
parser.add_argument('-c1', dest='country', type=str, help="country", default="xxxx")
parser.add_argument('-s2', dest='site', type=str, help="name of site", default="xxxx")
parser.add_argument('-s3', dest='sampleID', type=str, help="ID of sample", default="xxxx")
parser.add_argument('-g1', dest='genus', type=str, help="name of genus", default="xxxx")
parser.add_argument('-s4', dest='sampletype', type=str, help="type of sample", default="xxxx")
parser.add_argument('-m1', dest='markers', type=str, help="name of markers", default="PfK13,PfCRT,PfMDR,MT,CytB,PfDHPS,PfDHFR")
parser.add_argument('-r1', dest='replicate', type=str, help="name of replicate", default="xxxx")
parser.add_argument('-v1', dest='variantCalls', type=str, help="variantCalls", default="xxxx")
parser.add_argument('-t1', dest='treatment', type=str, help="day of treatment", default="xxxx")
parser.add_argument('-a1', dest='analysistime', type=str, help="day of treatment", default="now")
parser.add_argument('-f1', dest='file', type=str, help="name of filterd otuput", default=None)
args = parser.parse_args() 

jsonFile = open('teststudy.json', 'a')

json_dict = OrderedDict({'Study': args.study,
            'Number of samples': args.numsam, 'Date': args.analysistime})

wholelist1=[]
#from the final output file
if args.file != None:
    with open(args.file, "r") as f1:
        for line in f1:
            templist1=[]
            for word in line.split(","):
                #print(word)
                templist1+=[word.strip("\n")]
            #print(templist1)
            wholelist1+=[templist1]

#print(wholelist1)
#name of samples
index=[str(args.file).strip(".csv")]
json_dict['Sample'] = {}

for x in index:
    json_dict['Sample'][x] = {}
    json_dict['Sample'][x]['Year'] = args.study
    json_dict['Sample'][x]['Country'] = args.country
    json_dict['Sample'][x]['Site'] = args.site
    json_dict['Sample'][x]['TreatmentDay'] = args.treatment
    json_dict['Sample'][x]['SampleID'] = args.sampleID
    json_dict['Sample'][x]['Genus&Species'] = args.genus
    json_dict['Sample'][x]['SampleType'] = args.sampletype
    json_dict['Sample'][x]['Markers'] = args.markers
    json_dict['Sample'][x]['Replicate'] = args.replicate
    json_dict['Sample'][x]['VariantCalls'] = {}
    for k in range(1,len(wholelist1)):
        #print(wholelist1[x])
        #print(wholelist1[x][18])
        json_dict['Sample'][x]['VariantCalls'][wholelist1[k][0]+":"+str(wholelist1[k][5])+str(wholelist1[k][7])+str(wholelist1[k][6])] = {
                        'Ref' : wholelist1[k][5], 'Pos': wholelist1[k][7], 
                        'Alt': wholelist1[k][6] , 'Call' :  wholelist1[k][10],
                        'AF' : wholelist1[k][9], 'Confidence': wholelist1[k][15],
                        'Status' : wholelist1[k][21]}

json.dump(json_dict, jsonFile, indent=4)
jsonFile.close()