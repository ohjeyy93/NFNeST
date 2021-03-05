import argparse
import pandas as pd
import subprocess
import csv
import itertools

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-v1', dest='unfiltered', type=str, help="name of unfilterd merged vcf file")
parser.add_argument('-v2', dest='filtered', type=str, help="name of filtered merged vcf file")
args = parser.parse_args()

AAdic= {"Ala":"A", "Arg":"R", "Asn":"N", "Asp":"D", "Asx":"B", "Cys":"C", "Glu":"E", "Gln":"Q", "Glx":"Z", "Gly":"G", "His":"H", "Ile":"I", "Leu":"L", "Lys":"K", "Met":"M", "Phe":"F", "Pro":"P", "Ser":"S", "Thr":"T", "Trp":"W", "Tyr":"Y", "Val":"V", "STOP":"STOP"}

currentcodon=""
keywords = [''.join(i) for i in itertools.product(["T","G","A","C"], repeat = 3)]
codondic={}
for x in keywords:
    currentcodon=""
    for y in x:
        currentcodon+=y
    #print(currentcodon)
    if currentcodon == "TTT" or currentcodon == "TTC":
        #print("True")
        codondic[currentcodon]="Phe"
    if currentcodon == "TTA" or currentcodon == "TTG": 
        codondic[currentcodon]="Leu"
    if currentcodon == "CTT" or currentcodon =="CTC" or currentcodon =="CTA" or currentcodon =="CTG":
        codondic[currentcodon]="Leu"
    if currentcodon == "ATT" or currentcodon =="ATC" or currentcodon =="ATA":
        codondic[currentcodon]="Ile"
    if currentcodon == "ATG":
        codondic[currentcodon]="Met"
    if currentcodon == "GTT" or currentcodon =="GTC" or currentcodon =="GTA" or currentcodon =="GTG":
        codondic[currentcodon]="Val"
    if currentcodon == "TCT" or currentcodon =="TCC" or currentcodon =="TCA" or currentcodon =="TCG":
        codondic[currentcodon]="Ser"
    if currentcodon == "CCT" or currentcodon =="CCC" or currentcodon =="CCA" or currentcodon =="CCG": 
        codondic[currentcodon]="Pro"
    if currentcodon == "ACT" or currentcodon =="ACC" or currentcodon =="ACA" or currentcodon =="ACG":
        codondic[currentcodon]="Thr"
    if currentcodon == "GCT" or currentcodon =="GCC" or currentcodon =="GCA" or currentcodon =="GCG":
        codondic[currentcodon]="Ala"
    if currentcodon == "TAT" or currentcodon =="TAC":
        codondic[currentcodon]="Tyr"
    if currentcodon == "CAT" or currentcodon =="CAC":
        codondic[currentcodon]="His"
    if currentcodon == "CAA" or currentcodon =="CAG":
        codondic[currentcodon]="Gln"
    if currentcodon == "AAT" or currentcodon =="AAC":
        codondic[currentcodon]="Asn"
    if currentcodon == "AAA" or currentcodon =="AAG":
        codondic[currentcodon]="Lys"
    if currentcodon == "GAT" or currentcodon =="GAC":
        codondic[currentcodon]="Asp"
    if currentcodon == "GAA" or currentcodon =="GAG":
        codondic[currentcodon]="Glu"
    if currentcodon == "TGT" or currentcodon== "TGC":
        #print("I am correct")
        codondic[currentcodon]="Cys"
    if currentcodon == "TGG":
        codondic[currentcodon]="Trp"
    if currentcodon == "AGA" or currentcodon =="AGG" or currentcodon =="CGT" or currentcodon =="CGC" or currentcodon =="CGA" or currentcodon =="CGG":
        codondic[currentcodon]="Arg"
    if currentcodon == "AGT" or currentcodon == "AGC":
        codondic[currentcodon]="Ser"
    if currentcodon == "GGT" or currentcodon =="GGC" or currentcodon =="GGA" or currentcodon =="GGG":
        codondic[currentcodon]="Gly"
    if currentcodon == "TGA" or currentcodon == "TAA" or currentcodon == "TAG":
        codondic[currentcodon]="STOP"


with open(args.filtered, "r") as f1:
    #with open(args.fasta, "r") as f2:
    tempGene1=""
    tempPOS1=""
    Genelist1=[]
    POSlist1=[]
    GenePOSSet1=[]
    for line in f1:
        if line.find("ANN[*].EFFECT")==-1:
            count=0
            for word in line.split():
                #print(line)
                if count==0:
                    tempGene1=word
                    #print(word)
                if count==1:
                    tempPOS1=word
                if count==2:
                    if [tempGene1,tempPOS1] not in GenePOSSet1:
                        Genelist1+=[tempGene1]
                        POSlist1+=[tempPOS1]
                        GenePOSSet1+=[[tempGene1,tempPOS1]]
                count+=1



    filterlist1=[]
    f2= open(args.unfiltered, "r")
    f3= open(args.filtered, "r")
    wholefilnumsx=[]
    wholefilnums=[]
    wholenumset=[]
    for line in f3: 
        count=0
        tempnum=""
        tempGene2=""
        #print(line)
        if line.find("ANN[*].EFFECT")==-1:
            for word in line.split():
                count+=1
                if count==1:
                    tempGene2=word
                if count==2:
                    tempnum=word
                #print(tempnum)
                    #print(word)
                #print(tempnum)
                if word.startswith("c.") or word.startswith("n."):
                    #print(word)
                    if word.find("*")!=-1 and [tempGene2,tempnum] not in wholenumset:
                        wholefilnumsx+=[tempnum]
                        #print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                        wholefilnums+=[word[word.find("c.*"):word.find(tempnum)+len(tempnum)]]
                        wholenumset+=[[tempGene2,tempnum]]
                    if word.find("*")==-1 and [tempGene2,tempnum] not in wholenumset:
                        #print(word[word.find("c.*"):word.find(tempnum)+len(tempnum)])
                        wholefilnumsx+=[tempnum]
                        wholefilnums+=[word[word.find("c."):word.find(tempnum)+len(tempnum)]]
                        wholenumset+=[[tempGene2,tempnum]]
                    if word.find("*")!=-1 and [tempGene2,tempnum] not in wholenumset:
                        wholefilnumsx+=[tempnum]
                        wholefilnums+=[word[word.find("n.*"):word.find(tempnum)+len(tempnum)]]
                        wholenumset+=[[tempGene2,tempnum]]
                    if word.find("*")==-1 and [tempGene2,tempnum] not in wholenumset:
                        wholefilnumsx+=[tempnum]
                        wholefilnums+=[word[word.find("n."):word.find(tempnum)+len(tempnum)]]
                        wholenumset+=[[tempGene2,tempnum]]

    #print(wholefilnumsx)
    #print(wholefilnums)


    Genelist2=[]
    POSlist2=[]
    GenePosSet1=[]
    Geneset1=[]
    customAA1=[]
    for k in range(len(Genelist1)):
        Geneset1+=[[Genelist1[k],POSlist1[k]]]

    f2= open(args.unfiltered, "r")
    for line in f2:
        tempcount=0
        count=0
        for word in line.split():
            #print(line)
            #print(Genelist1[x])
            #print(word)
            #if count==0 and word == Genelist1[x]:
                #print(word+"test1111")
            if count==0:
                word1=""
                if word.startswith("#")==False:
                    word1=word
                    #print(word1)
            #print(word1)
            #print(count)
            #print(word1)
            #temppos3=""
            #tempgene1=""
            if count==1:
                #print(word1)
                #print(word1)
                if word1!="":
                    #filtercount=-1
                    #print(word1)
                    #print(word)
                    #print(POSlist1)
                    #print(Genelist1)
                    for x in range(len(Genelist1)):
                        #print(str(int(POSlist1[x])-1))
                        #if word1== Genelist1[x]:print(word1,word)
                        #print(POSlist1)
                        #print(word)
                        #if word.startswith("c."):
                            #temppos2=""
                            #for x in word:
                                #print(word)
                                #if x.isdigit():
                                    #print(x)
                                    #temppos2+=x
                            #print(temppos2)
                        #print(word1)
                        #print(Genelist1[x])
                        #print(POSlist1[x])
                        #print(word)
                        #print(Genelist1[x],POSlist1[x])
                        #print(word1,word)
                        #print(POSlist1)
                        #print(type(POSlist1[x]))
                        #if word1 == Genelist1[x]:print("no")
                        #print(word)
                        #print(POSlist1[x])
                        #if word == POSlist1[x]: print("yes")
                        if word1== Genelist1[x] and word== POSlist1[x] and line.startswith("#")==False:
                            #print(word1)
                            #print(word)
                            #print(word)
                            #print(word)
                            #filtercount+=1
                            #print(POSlist1[x])
                            #print("test"+POSlist1[x])
                            #print(Genelist1[x],POSlist1[x])
                            Genelist2+=[Genelist1[x]]
                            POSlist2+=[POSlist1[x]]
                            GenePosSet1+=[[Genelist1[x],POSlist1[x]]]
                            #print(line.find("CUSTOM&annotation")+18)
                            #print(line)
                            if line.find("CUSTOM&annotation")!=-1:
                                #print(line.find("CUSTOM&annotation"))
                                #print(line[2::])
                                #print(line[line.find("CUSTOM&annotation")+18::])
                                #print(line[line.find("CUSTOM&annotation")+18::].find("|"))
                                customAA1+=[line[line.find("CUSTOM&annotation")+18:line.find("CUSTOM&annotation")+18+line[line.find("CUSTOM&annotation")+18::].find("|")]]
                            else:
                                customAA1+=["None"]
                            #temppos3=POSlist1[x]
            #print(len(QDlist1))
            #print(GenePosSet1)
            #print(AFlist1)
            if count==7:
                if word1 != "":
                    #print(word)
                    if line.startswith("#")==False:
                        if word.startswith("AC=") or word.startswith("AB="):
                            #print(word[word.find("c."):word.find("c.")+10])
                            #print(line)
                            for x in range(len(wholefilnums)):
                                #print(x)
                                #print(wholefilnumsx[count])
                                #print(POSlist2)
                                #print(wholefilnumsx[count2])
                                #print([word1,str(int(wholefilnumsx[x]))])
                                #print(word)
                                #find part of c.138 in c.1381 in different gene.....
                                #print(line.find(wholefilnums[x]))
                                #if word.find(wholefilnums[x]) != -1:
                                    #print(wholefilnums[x])
                                    #print(word.find(wholefilnums[x]))
                                    #print(len(wholefilnums[x]))
                                    #print(wholefilnums[x])
                                    #print(word[word.find(wholefilnums[x])+len(wholefilnums[x])])
                                    #print(word[word.find(wholefilnums[x])+len(wholefilnums[x])-len(wholefilnumsx[x])-1])
                                if word.find(wholefilnums[x])!=-1 and str(word.find(wholefilnums[x])+len(wholefilnums[x])+1).isdigit()==False and [word1,wholefilnumsx[x]] in Geneset1 and [word1,str(int(wholefilnumsx[x]))] not in GenePosSet1:
                                    #print(word1)
                                    #print(wholefilnumsx[x])
                                    #print(wholefilnums[x])
                                    #print(wholefilnumsx[x])
                                    #print(word1)
                                    #print(word1)
                                    #print(wholefilnumsx[x])
                                    Genelist2+=[word1]
                                    #print(wholefilnumsx[x])
                                    POSlist2+=[wholefilnumsx[x]]
                                    #"QD": QDlist1, "SOR":SORlist1, "MQ":MQlist1, "MQRankSum":MQRankSumlist1,
                                    GenePosSet1+=[[word1,str(int(wholefilnumsx[x]))]]
                                    if line.find("CUSTOM&annotation")!=-1:
                                #print(line.find("CUSTOM&annotation"))
                                #print(line[2::])
                                #print(line[line.find("CUSTOM&annotation")+18::])
                                #print(line[line.find("CUSTOM&annotation")+18::].find("|"))
                                        customAA1+=[line[line.find("CUSTOM&annotation")+18:line.find("CUSTOM&annotation")+18+line[line.find("CUSTOM&annotation")+18::].find("|")]]
                                    else:
                                        customAA1+=["None"]
                                #print(Genelist1[count])
                                #print(wholefilnumsx[count])
            count+=1

#print(Genelist1)
#print(POSlist1)
#print(Genelist2)
#print(POSlist2)
#print(customAA1)
d3={"Gene":Genelist2, "BasePOS":POSlist2, "customAA":customAA1}
df3=pd.DataFrame(data=d3)
df3.to_csv('custom1.csv', index=False,) 