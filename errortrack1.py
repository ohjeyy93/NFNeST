from os import listdir
import pandas as pd
import csv

wholelist1=[["samplename","error", "CTvalue"]]
with open(".nextflow.log", "r") as logfile:
    for line in logfile:
        if line.find("error") != -1:
            if line.find("INFO") != -1:
                #print(line[line.find("["):line.find("]")+1])
                #print((line[line.find("]")+1::]))
                #print((line[line.find("]")+1::])[(line[line.find("]")+1::]).find("["):(line[line.find("]")+1::]).find("]")+1])
                workdir=(line[line.find("]")+1::])[(line[line.find("]")+1::]).find("["):(line[line.find("]")+1::]).find("]")+1]
                #print(workdir)
                #print(workdir[1:workdir.find("/")])
                #print(workdir[workdir.find("/")+1:-1])
                workdir1=workdir[1:workdir.find("/")]
                workdir2=workdir[workdir.find("/")+1:-1]
                for x in listdir("work/"+workdir1):
                    if x.startswith(workdir2):
                        workdir2=x
                with open("work/"+workdir1+"/"+workdir2+"/.command.sh", "r") as commandsh:
                    for line in commandsh:
                        if line.find("samtools index") !=-1:
                            samplename=line[14::]
                            #print(samplename)
                        if line.find("bbduk.sh") !=-1:
                            samplename=line[line.find("in=")+3:line.find("_R1")]
                with open("work/"+workdir1+"/"+workdir2+"/.command.err", "r") as commanderror:
                    for line in commanderror:
                        if line.find("KeyError: 'fields")!=-1:
                            error="noVCF?fields"
                            #print(error)
                            #print(samplename,error)
                        if line.find("KeyError: 'info")!=-1:
                            error="noVCF?info"
                            #print(samplename,error)
                        if line.find("There appear to be different numbers of reads in the paired input files.")!=-1:
                            error="different pair read"
                            #print(samplename,error)
                samplename=samplename[:samplename.find("Fxxx0")+1]+"1230"
                xls=pd.ExcelFile("ANG_2019_TES_master.xlsx")
                df1 = pd.read_excel(xls, 'TES PCR gel results')
                #print(df1.columns)
                #df1.to_csv("test.csv", index=False)
                #print(df1["Date completed"].tolist())
                #print(df1["Unnamed: 2"].tolist())
                #print(df1['Date completed.1'].tolist())
                #print(df1["Unnamed: 15"].tolist())
                #print(samplename)
                for x in range(len(df1["Date completed"].tolist())):
                    #print((df1["Date completed"].tolist())[x])
                    #print(samplename)
                    if (df1["Date completed"].tolist())[x]==samplename:
                        #print("Test")
                        CTvalue=(df1["Unnamed: 2"].tolist())[x]
                        wholelist1+=[samplename,error,CTvalue]
                for x in range(len(df1["Date completed.1"].tolist())):
                    if (df1["Date completed.1"].tolist())[x]==samplename:
                        #print("Test")
                        #print((df1["Date completed"].tolist())[x])
                        #print(samplename)
                        CTvalue=(df1["Unnamed: 15"].tolist())[x]
                        wholelist1+=[[samplename,error,CTvalue]]
#print(wholelist1)
with open("out.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(wholelist1)

                