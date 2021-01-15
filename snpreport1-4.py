import argparse
import pandas as pd
import subprocess
import csv

parser = argparse.ArgumentParser(description='name')
parser.add_argument('-v1', dest='unfiltered', type=str, help="name of unfilterd merged vcf file")
parser.add_argument('-v2', dest='filtered', type=str, help="name of filtered merged vcf file")
parser.add_argument('-b1', dest='bam', type=str, help="name of bam file")
parser.add_argument('-o1', dest='name', type=str, help="name of output")
args = parser.parse_args()

AAdic= {"Ala":"A", "Arg":"R", "Asn":"N", "Asp":"D", "Asx":"B", "Cys":"C", "Glu":"E", "Gln":"Q", "Glx":"Z", "Gly":"G", "His":"H", "Ile":"I", "Leu":"L", "Lys":"K", "Met":"M", "Phe":"F", "Pro":"P", "Ser":"S", "Thr":"T", "Trp":"W", "Tyr":"Y", "Val":"V"}


with open(args.filtered, "r") as f1:
    Genelist1=[]
    POSlist1=[]
    Reflist1=[]
    Altlist1=[]
    AAlist1=[]
    AAlist2=[]
    for line in f1:
        count=0
        if "p." in line:
            for word in line.split():
                #print(word)
                if count==0:
                    Genelist1+=[word]
                if count==1:
                    POSlist1+=[word]
                if count==2:
                    Reflist1+=[word]
                if count==3:
                    Altlist1+=[word]
                if word.startswith("p."):
                    if word[2:5] in AAdic:
                        AAlist1+=[AAdic[word[2:5]]]
                    else:
                        Genelist1=Genelist1[:-1]
                        POSlist1=POSlist1[:-1]
                        Reflist1=Reflist1[:-1]
                        Altlist1=Altlist1[:-1]
                        break
                    if word[-3:] in AAdic:
                        AAlist2+=[AAdic[word[-3:]]]
                    else:
                        Genelist1=Genelist1[:-1]
                        POSlist1=POSlist1[:-1]
                        Reflist1=Reflist1[:-1]
                        Altlist1=Altlist1[:-1]
                        AAlist1=AAlist1[:-1]
                        break
                count+=1
            #print(Genelist1)
            #print(POSlist1)
            #print(Reflist1)
            #print(Altlist1)
            #print(AAlist1)
            #print(AAlist2)
    process=subprocess.Popen(['samtools', 'depth', '-a', args.name+"_SR.bam"], stdout=subprocess.PIPE)
    stdout, stderr = process.communicate()
    reader = csv.reader(stdout.decode('ascii').splitlines(),
                            delimiter='\t', skipinitialspace=True)
    depthlist1=[]
    for row in reader:
        depthlist1+=[row]
    depthlist1up=[]
    #for x in range(len(Genelist1)):
    #    if Genelist1[x] in depthlist1:
    #        print(Genelist1[x])
    #    if POSlist1[x] in depthlist1:
    #        print(POSlist1[x])
    for x in range(len(Genelist1)):
        for y in range(len(depthlist1)):
            if Genelist1[x] in depthlist1[y] and POSlist1[x] in depthlist1[y]:
                #print(depthlist1[y])
                depthlist1up+=[depthlist1[y][2]]
                break

    
    #print(len(Genelist1))
    #print(len(POSlist1))
    #print(len(depthlist1up))
    #print(depthlist1up)
    #with open(args.unfiltered, "r") as f2:
    #    for line in f2:
    #        print(line)
   # print(len(Genelist1))
    Genelist2=[]
    POSlist2=[]
    AFlist1=[]
    QDlist1=[]
    SORlist1=[]
    MQlist1=[]
    MQRankSumlist1=[]
    filterlist1=[]
    f2= open(args.unfiltered, "r")
    for line in f2:
        count=0
        for word in line.split():
            #print(line)
            #print(Genelist1[x])
            #print(word)
            #if count==0 and word == Genelist1[x]:
                #print(word+"test1111")
            if count==0:
                word1=word
            if count==1:
                filtercount=-1
                for x in range(len(Genelist1)):
                    if word1== Genelist1[x] and word== str(int(POSlist1[x])-1) and line.startswith("#")==False:
                        filtercount+=1
                        if line.startswith("#")==False and line.find("AF=")!=-1 and line.find("AN=") !=-1 and line.find("QD=")!=-1 and line.find("ReadPosRankSum=") !=-1:
                            AFlist1+=[float(line[line.find("AF=")+3:line.find("AN=")-1])]
                            QDlist1+=[float(line[line.find("QD=")+3:line.find("ReadPosRankSum=")-1])]
                            Genelist2+=[Genelist1[x]]
                            POSlist2+=[POSlist1[x]]
                        elif line.startswith("#")==False and line.find("AF=")!=-1 and line.find("QD=")!=-1 and line.find("ReadPosRankSum=") ==-1:
                            AFlist1+=[float(line[line.find("AF=")+3:line.find("AN=")-1])]
                            QDlist1+=[float(line[line.find("QD=")+3:line.find("SOR=")-1])]
                            Genelist2+=[Genelist1[x]]
                            POSlist2+=[POSlist1[x]]
                        if line.startswith("#")==False and line.find("MQ=")!=-1 and line.find("MQRankSum") != -1:
                            MQlist1+=[float(line[line.find("MQ=")+3:line.find("MQRankSum=")-1])]
                        elif line.startswith("#")==False and line.find("MQ=")!=-1 and line.find("MQRankSum") == -1:
                            MQlist1+=[float(line[line.find("MQ=")+3:line.find("QD=")-1])]
                            MQRankSumlist1+=["None"]
                        if line.startswith("#")==False and line.find("SOR=")!=-1:
                            SORlist1+=[float(line[line.find("SOR=")+4:line.find("ANN=")-1])]
                        if line.startswith("#")==False and line.find("MQRankSum") != -1:
                            MQRankSumlist1+=[float(line[line.find("MQRankSum=")+10:line.find("QD=")-1])]
                        if line.startswith("#")==False and line.find("AF=")!=-1 and line.find("AN=") !=-1 and line.find("QD=")!=-1:
                            if QDlist1[filtercount] != "None" and QDlist1[filtercount]<15 or  SORlist1[filtercount]!="None"and  SORlist1[filtercount] >1 or MQlist1[filtercount]!="None" and MQlist1[filtercount]<35 or MQRankSumlist1[filtercount]!="None" and MQRankSumlist1[filtercount] <-2:
                                filterlist1+=["bad"]
                            else: filterlist1+=["good"]
            count+=1
                        #print(word+"test2222")
    #for x in range(len(Genelist1)):
        #print(Genelist1[x])
    #    for line in f2:
            #count=0
            #print(line)
    #        for word in line.split():
                #if count==0:
                    #print(word + "test")
    #            if Genelist1[x] == word:
    #                print(Genelist1[x] +"test222")
                    #    print(line)
                    #    break
            #    count+=1
                #if count==1:
                #    if POSlist1[x] == word:
                #        print(line)
    Reflist2=[]
    Altlist2=[]
    AAlist21=[]
    AAlist22=[]
    depthlist2up=[]
    for x in range(len(Genelist2)):
        for y in range(len(Genelist1)):
            if Genelist2[x]==Genelist1[y] and POSlist2[x]==POSlist1[y]:
                Reflist2+=[Reflist1[y]]
                Altlist2+=[Altlist1[y]]
                AAlist21+=[AAlist1[y]]
                AAlist22+=[AAlist2[y]]
                depthlist2up+=[depthlist1up[y]]
    mutationlist1=[]
    for x in range(len(AFlist1)):
        if AFlist1[x]==1:mutationlist1+=["Major"]
        if AFlist1[x]==0.5:mutationlist1+=["Minor"]
        if AFlist1[x]==0:mutationlist1+=["Wildtype"]

    #print(len(AFlist1))
    print(len(QDlist1))
    print(len(MQlist1))
    print(len(SORlist1))
    print(len(MQRankSumlist1))
    print(len(filterlist1))
    #d1={"Gene":Genelist1, "POS":POSlist1, "Ref":Reflist1, "Alt":Altlist1, "AAref":AAlist1, "AAalt":AAlist2, "Depth":depthlist1up}
    #df=pd.DataFrame(data=d1)
    #d2={"Gene":Genelist2, "BasePOS":POSlist2, "Ref":Reflist2, "Alt":Altlist2, "AAref":AAlist21, "AAalt":AAlist22, "Depth":depthlist2up, "AF": AFlist1, "Mutation": mutationlist1, "QD": QDlist1}
    #df2=pd.DataFrame(data=d2)
    d3={"Gene":Genelist2, "BasePOS":POSlist2, "Ref":Reflist2, "Alt":Altlist2, "AAref":AAlist21, "AAalt":AAlist22, "Depth":depthlist2up, "AF": AFlist1, "Mutation": mutationlist1, "QD": QDlist1, "SOR":SORlist1, "MQ":MQlist1, "MQRankSum":MQRankSumlist1, "Filter":filterlist1}
    df3=pd.DataFrame(data=d3)
    #df.to_csv('snpreport4.csv', index=False,) 
    df3.to_csv(args.name+'.csv', index=False,) 
