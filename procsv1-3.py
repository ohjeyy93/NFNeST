import csv
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

Druglist1 = ["Chloroquine", "Piperaquine", "Piperaquine", "Chloroquine", "Mefloquine", "Pyrimethamine", "Proguanil", "Sulfadoxine", "Artemisinin and its derivatives"]
Genelist1 = ["Pfcrt", "Pfcrt", "Pfcrt", "Pfmdr1 (in combination with Pfcrt mutations only)", "Pfmdr1", "Pfdhfr", "Pfdhfr", "Pfdhps","PfK13"]
Mutationlist1 = ["K76T + different sets of mutations at other codons (including C72S, M74I, N75E, A220S, Q271E, N326S, I356T and R371I)", "Detected in vivo: T93S, H97Y, F145I, I218F and C350R", "Detected in vitro: T93S, H97Y, F145I, I218F, M343L and G353V", "N51I, C59R, S108N and I164L", "N86Y, Y184F, S1034C, N1042D and D1246Y", "Pfmdr1 increased copy number", "A16V, N51I, C59R, S108N and I164L", "S436A/F, A437G, K540E, A581G and A613T/S",  "P441L, G449A, C469F/Y, A481V, R515K, P527H, N537I/D, G538V, V568G, R622I, A675V"] 
#print(len(Druglist1))
#print(len(Genelist1))
#print(len(Mutationlist1))
d= {'Gene':Genelist1,'Drug':Druglist1, "Mutation":Mutationlist1}
df = pd.DataFrame(data=d)
df.to_csv('report1-3-2.csv', index=False)

ChrList1 = ["Pfdhps", "Pfdhfr", "Pfmdr1", "PfK13", "Pfcrt"]
Genelist2 = ["Pfdhps", "Pfdhfr", "Pfmdr1", "PfK13", "Pfcrt"]
Distancelist1 = [1487,1827, 2157,  2181, 2586]
Pfcrtgene1=['K76T', "C72S", "M74I", "N75E", "T93S", "H97Y", "F145I", "I218F", "A220S", "Q271E", "N326S", "M343L","C350R", "G353V", "I356T", "R371I"]
Pfcrtlist1=[76,72,74,75,93,97, 145,218, 220,271,326, 343, 350,353, 356,371,]
Pfmdr1gene1=['N51I', "C59R", "N86Y", "S108N", "I164L", "Y184F", "S1034C", "N1042D", "D1246Y"]
Pfmdr1=[51,59,86, 108, 164, 184, 1034, 1042, 1246]
Pfk13gene1=["P441L", "G449A", "C469F/Y", "A481V", "R515K", "P527H", "N537I/D", "G538V", "V568G", "R622I", "A675V"]
Pfk13=[441,449,469,481,515,527,537,538,568,622,675]
Pfdhfrgene1=["A16V", "N51I", "C59R", "S108N","I164L"]
Pfdhfr=[16, 51, 59, 108, 164]
Pfdhpsgene1=["S436A/F", "A437G", "K540E", "A581G", "A613T/S"]
Pfdhps=[436,437,540,581,613]

Sumlist1=[]
Sumlen=len(Pfcrtgene1)+len(Pfmdr1gene1)+len(Pfk13gene1)+len(Pfdhfrgene1)+len(Pfdhpsgene1)
print(Sumlen)
ChrList1 = []
Genelist3 = []
RefAA1 = []
AAPos1= []
AltAA1= []
for x in range(len(Pfcrtgene1)):
    ChrList1+=["Pfcrt"]
    Genelist3+=["Pfcrt"]
    RefAA1+=[Pfcrtgene1[x][0]]
    if "/" in Pfcrtgene1[x]:
        ChrList1+=["Pfcrt"]
        Genelist3+=["Pfcrt"]
        RefAA1+=[Pfcrtgene1[x][0]]
        AAPos1+=[Pfcrtgene1[x][1:len(Pfcrtgene1[x])-3]]
        AAPos1+=[Pfcrtgene1[x][1:len(Pfcrtgene1[x])-3]]
        AltAA1+=[Pfcrtgene1[x][-3]]
        AltAA1+=[Pfcrtgene1[x][-1]]
    if "/" not in Pfcrtgene1[x]:
        AAPos1+=[Pfcrtgene1[x][1:len(Pfcrtgene1[x])-1]]
    if "/" not in Pfcrtgene1[x]:
        AltAA1+=[Pfcrtgene1[x][-1:]]
for x in range(len(Pfmdr1gene1)):
    ChrList1+=["Pfmdr1"]
    Genelist3+=["Pfmdr1"]
    RefAA1+=[Pfmdr1gene1[x][0]]
    if "/" in Pfmdr1gene1[x]:
        ChrList1+=["Pfmdr1"]
        Genelist3+=["Pfmdr1"]
        RefAA1+=[Pfmdr1gene1[x][0]]
        AAPos1+=[Pfmdr1gene1[x][1:len(Pfmdr1gene1[x])-3]]
        AltAA1+=[Pfmdr1gene1[x][-3]]
        AltAA1+=[Pfmdr1gene1[x][-1]]
        AAPos1+=[Pfmdr1gene1[x][1:len(Pfmdr1gene1[x])-3]]
    if "/" not in Pfcrtgene1[x]:
        AAPos1+=[Pfmdr1gene1[x][1:len(Pfmdr1gene1[x])-1]]
    if "/" not in Pfmdr1gene1[x]:
        AltAA1+=[Pfmdr1gene1[x][-1:]]
for x in range(len(Pfk13gene1)):
    ChrList1+=["Pfk13"]
    Genelist3+=["Pfk13"]
    RefAA1+=[Pfk13gene1[x][0]]
    if "/" in Pfk13gene1[x]:
        ChrList1+=["Pfk13"]
        Genelist3+=["Pfk13"]
        RefAA1+=[Pfk13gene1[x][0]]
        AAPos1+=[Pfk13gene1[x][1:len(Pfk13gene1[x])-3]]
        AltAA1+=[Pfk13gene1[x][-3]]
        AltAA1+=[Pfk13gene1[x][-1]]
        AAPos1+=[Pfk13gene1[x][1:len(Pfk13gene1[x])-3]]
    if "/" not in Pfk13gene1[x]:
        AAPos1+=[Pfk13gene1[x][1:len(Pfk13gene1[x])-1]]
    if "/" not in Pfk13gene1[x]:
        AltAA1+=[Pfk13gene1[x][-1:]]
for x in range(len(Pfdhfrgene1)):
    ChrList1+=["Pfdhfr"]
    Genelist3+=["Pfdhfr"]
    RefAA1+=[Pfdhfrgene1[x][0]]
    if "/" in Pfdhfrgene1[x]:
        ChrList1+=["Pfdhfr"]
        Genelist3+=["Pfdhfr"]
        RefAA1+=[Pfdhfrgene1[x][0]]
        AAPos1+=[Pfdhfrgene1[x][1:len(Pfdhfrgene1[x])-3]]
        AltAA1+=[Pfdhfrgene1[x][-3]]
        AltAA1+=[Pfdhfrgene1[x][-1]]
        AAPos1+=[Pfdhfrgene1[x][1:len(Pfdhfrgene1[x])-3]]
    if "/" not in Pfdhfrgene1[x]:
        AAPos1+=[Pfdhfrgene1[x][1:len(Pfdhfrgene1[x])-1]]
    if "/" not in Pfdhfrgene1[x]:
        AltAA1+=[Pfdhfrgene1[x][-1:]]
for x in range(len(Pfdhpsgene1)):
    ChrList1+=["Pfdhps"]
    Genelist3+=["Pfdhps"]
    RefAA1+=[Pfdhpsgene1[x][0]]
    if "/" in Pfdhpsgene1[x]:
        ChrList1+=["Pfdhps"]
        Genelist3+=["Pfdhps"]
        RefAA1+=[Pfdhpsgene1[x][0]]
        AAPos1+=[Pfdhpsgene1[x][1:len(Pfdhpsgene1[x])-3]]
        AltAA1+=[Pfdhpsgene1[x][-3]]
        AltAA1+=[Pfdhpsgene1[x][-1]]
        AAPos1+=[Pfdhpsgene1[x][1:len(Pfdhpsgene1[x])-3]]
    if "/" not in Pfdhpsgene1[x]:
        AAPos1+=[Pfdhpsgene1[x][1:len(Pfdhpsgene1[x])-1]]
    if "/" not in Pfdhpsgene1[x]:
        AltAA1+=[Pfdhpsgene1[x][-1:]]

print(ChrList1)
print(Genelist3)
print(RefAA1)
print(AAPos1)
print(AltAA1)
print(len(ChrList1))
print(len(Genelist3))
print(len(RefAA1))
print(len(AAPos1))
print(len(AltAA1))



d2= {'Gene': Genelist2, 'Distance': Distancelist1}
d3= {'Pfcrt': Pfcrtgene1, 'Distance': Pfcrtlist1}
d4= {'Pfmdr1': Pfmdr1gene1, 'Distance': Pfmdr1}
d5= {'Pfk13': Pfk13gene1, 'Distance': Pfk13}
d6= {'Pfdhfr': Pfdhfrgene1, 'Distance': Pfdhfr}
d7= {'Pfdhps': Pfdhpsgene1, 'Distance': Pfdhps}

d8= {'Chr': ChrList1, 'Gene': Genelist3, 'RefAA': RefAA1, 'AAPos': AAPos1, 'AltAA': AltAA1}

#df2 = pd.DataFrame(data=d2)
df3 = pd.DataFrame(data=d3)
df4 = pd.DataFrame(data=d4)
df5 = pd.DataFrame(data=d5)
df6 = pd.DataFrame(data=d6)
df7 = pd.DataFrame(data=d7)
df8 = pd.DataFrame(data=d8)

df8.to_csv('voi2.csv', index=False,) 
