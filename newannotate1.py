import os
import copy
import logging 
import itertools
from collections import namedtuple
from collections import OrderedDict
from nest.parsers.vcfReader import Reader 
from nest.parsers.vcfwriter import Writer

class Merge:

    def __init__(self, tmp_dir, vcf_dict, ref_path):
         ## Edit: (10/10/19): Added fasta path as requirement for initialization of merge
         ## Reason: Reader was modified to require fasta path, see vcfReader for more info
         self.out_path = tmp_dir
         self.vcf_dict = vcf_dict
         self.ref_path = ref_path

    def merge(self, vcflist1):
        dict1={}
        dict2={}
        dict3={}
        #print(vcflist1[0][0:-14])
        with open("final_"+vcflist1[0][0:-14]+".vcf", "w") as w1:
            for item in vcflist1:
                with open(item,"r") as w2:
                    for line in w2:
                        #print(line)
                        count=0
                        key=""
                        key2=""
                        for word in line.split():
                            #print("what")
                            #print(line)
                            #print(word)
                            count+=1
                            if count==1:
                                key=word
                                key2=word
                            if count==2:
                                key+=","+word
                                key2+=","+word
                            if count==4:
                                key+=","+word
                        #print(key)
                        if key in dict1 and "#" not in key:
                            countw=0
                            dict1[key]=dict1[key][0:dict1[key].find("VARTYPE=SNP")+11]
                            #print(dict1[key])
                            for word in line.split(";"):
                                countw+=1
                                #print(word)
                                #print(word[0:word.find("=")])
                                if word[0:word.find("=")] not in dict1[key]:
                                    #print(word)
                                    if key[0:key.find(",")] not in word and key[key.find(",")::] not in word:
                                        #if key=="PfCRT,211,G" and countw==len(line.split(";"))-4:
                                        #    print(word)
                                        #print(word)
                                        #print(len(line.split(";")))
                                        #print(item)
                                        #print(word)
                                        #print(dict1[key])
                                        #print(word)
                                        #print(countw)
                                        #print(len(line.split(";")))
                                        if countw<len(line.split(";")):
                                            #print(word)
                                            dict1[key]=dict1[key]+";"+word.strip("\n")
                                if countw==len(line.split(";")):
                                    #print(word)
                                    #print("True")
                                    dict1[key]+word[12::]
                                    if "_1." in item:
                                        dict2[key]="samtools"
                                        dict3[key]=1
                                    if "_2." in item:
                                        if key in dict2:
                                        #if dict2[key]=="samtools":
                                            dict2[key]="samtools,GATK"
                                            dict3[key]=2
                                        else:
                                            dict2[key]="GATK"
                                            dict3[key]=1
                                    if "_3." in item:
                                        if key in dict2:
                                            if dict2[key]=="samtools,GATK":
                                                dict2[key]="samtools,GATK,Freebayes"
                                                dict3[key]=3
                                            if dict2[key]=="samtools":
                                                dict2[key]="samtools,Freebayes"
                                                dict3[key]=2
                                            if dict2[key]=="GATK":
                                                dict2[key]="GATK,Freebayes"
                                                dict3[key]=2
                                        if key not in dict2:
                                            dict2[key]="Freebayes"    
                                            dict3[key]=1
                                    dict1[key]=dict1[key]+";"+word+"\n"
                        if key not in dict1:
                            if "#" not in line:
                                #print(line)
                                #print(line[0:line.find("VARTYPE=SNP")+12])
                                #dict1[key]=line[0:line.find("VARTYPE=SNP")+11]
                                #print(line)
                                dict1[key]=line
                            if "#" in line:
                                dict1[key]=line
                            #print(key)
                            if "_1." in item:
                                dict2[key]="samtools"
                                dict3[key]=1
                            if "_2." in item:
                                dict2[key]="GATK"
                                dict3[key]=1
                            if "_3." in item:
                                dict2[key]="Freebayes"
                                dict3[key]=1  
            for key in dict1: 
                #print(key)
                #    print(dict1[key]) 
                #print(key)
                #if key in dict2:
                #    print(dict2[key])
                #if key not in dict2:
                #    print(key)
                if "#" in dict1[key]:
                    w1.write(dict1[key])
                #else:
                    #print(dict1[key][:-1])
                #    w1.write(dict1[key][:-1]+";Confidence="+str(dict3[key])+";Sources="+dict2[key]+"\n")
            w1.write("""##INFO=<ID=Confidence,Number=1,Type=Integer,Description="Number of variant callers that identified this variant">\n""")
            w1.write("""##INFO=<ID=Sources,Number=1,Type=String,Description="Variant callers that called the variant">\n""")
            for key in dict1:
                if "#" not in dict1[key]:
                    #print(dict1[key][:-1])
                    w1.write(dict1[key][:-1]+";Confidence="+str(dict3[key])+";Sources="+dict2[key]+"\n")
        return w1