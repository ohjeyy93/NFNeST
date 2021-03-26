from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import statistics
import csv
import math

#######################2017
####K13#####

possum1K132017=[]
negsum1K132017=[]
posmean1K132017=[]
negmean1K132017=[]

possum1Pfcrt2017=[] 
negsum1Pfcrt2017=[]
posmean1Pfcrt2017=[]
negmean1Pfcrt2017=[]

possum1Pfcrt_newprimer2017=[]
negsum1Pfcrt_newprimer2017=[]
posmean1Pfcrt_newprimer2017=[]
negmean1Pfcrt_newprimer2017=[]

possum1MDR2017=[]
negsum1MDR2017=[]
posmean1MDR2017=[]
negmean1MDR2017=[]

possum1MDR_newprimer2017=[]
negsum1MDR_newprimer2017=[]
posmean1MDR_newprimer2017=[]
negmean1MDR_newprimer2017=[]

possum1CytoB2017=[]
negsum1CytoB2017=[]
posmean1CytoB2017=[]
negmean1CytoB2017=[]

possum1dhps2017=[]
negsum1dhps2017=[]
posmean1dhps2017=[]
negmean1dhps2017=[]

possum1dhfr2017=[]
negsum1dhfr2017=[]
posmean1dhfr2017=[]
negmean1dhfr2017=[]

possum1cpmp2017=[]
negsum1cpmp2017=[]
posmean1cpmp2017=[]
negmean1cpmp2017=[]

possum1Pfs472017=[]
negsum1Pfs472017=[]
posmean1Pfs472017=[]
negmean1Pfs472017=[]

#######################2018

possum1K132018=[]
negsum1K132018=[]
posmean1K132018=[]
negmean1K132018=[]

possum1Pfcrt2018=[]
negsum1Pfcrt2018=[]
posmean1Pfcrt2018=[]
negmean1Pfcrt2018=[]

possum1Pfcrt_newprimer2018=[]
negsum1Pfcrt_newprimer2018=[]
posmean1Pfcrt_newprimer2018=[]
negmean1Pfcrt_newprimer2018=[]

possum1MDR2018=[]
negsum1MDR2018=[]
posmean1MDR2018=[]
negmean1MDR2018=[]

possum1MDR_newprimer2018=[]
negsum1MDR_newprimer2018=[]
posmean1MDR_newprimer2018=[]
negmean1MDR_newprimer2018=[]

possum1CytoB2018=[]
negsum1CytoB2018=[]
posmean1CytoB2018=[]
negmean1CytoB2018=[]

possum1dhps2018=[]
negsum1dhps2018=[]
posmean1dhps2018=[]
negmean1dhps2018=[]

possum1dhfr2018=[]
negsum1dhfr2018=[]
posmean1dhfr2018=[]
negmean1dhfr2018=[]

possum1cpmp2018=[]
negsum1cpmp2018=[]
posmean1cpmp2018=[]
negmean1cpmp2018=[]

possum1Pfs472018=[]
negsum1Pfs472018=[]
posmean1Pfs472018=[]
negmean1Pfs472018=[]


#######################2019

possum1K132019=[]
negsum1K132019=[]
posmean1K132019=[]
negmean1K132019=[]

possum1Pfcrt2019=[]
negsum1Pfcrt2019=[]
posmean1Pfcrt2019=[]
negmean1Pfcrt2019=[]

possum1Pfcrt_newprimer2019=[]
negsum1Pfcrt_newprimer2019=[]
posmean1Pfcrt_newprimer2019=[]
negmean1Pfcrt_newprimer2019=[]

possum1MDR2019=[]
negsum1MDR2019=[]
posmean1MDR2019=[]
negmean1MDR2019=[]

possum1MDR_newprimer2019=[]
negsum1MDR_newprimer2019=[]
posmean1MDR_newprimer2019=[]
negmean1MDR_newprimer2019=[]

possum1CytoB2019=[]
negsum1CytoB2019=[]
posmean1CytoB2019=[]
negmean1CytoB2019=[]

possum1dhps2019=[]
negsum1dhps2019=[]
posmean1dhps2019=[]
negmean1dhps2019=[]

possum1dhfr2019=[]
negsum1dhfr2019=[]
posmean1dhfr2019=[]
negmean1dhfr2019=[]

possum1cpmp2019=[]
negsum1cpmp2019=[]
posmean1cpmp2019=[]
negmean1cpmp2019=[]

possum1Pfs472019=[]
negsum1Pfs472019=[]
posmean1Pfs472019=[]
negmean1Pfs472019=[]


with open('test.csv', "r") as csv1:
    reader = csv.reader(csv1)
    for row in reader:
        #print(row[1])
        if row[1]=="2017" and row[4]=="Negative Amplification":
            negsum1K132017+=[float(row[2])]
        if row[1]=="2017" and row[4]=="Positive Amplification":
            possum1K132017+=[float(row[2])]
        if row[1]=="2017" and row[5]=="Negative Amplification":
            negsum1Pfcrt2017+=[float(row[2])]
        if row[1]=="2017" and row[5]=="Positive Amplification":
            possum1Pfcrt2017+=[float(row[2])]
        if row[1]=="2017" and row[6]=="Negative Amplification":
            negsum1Pfcrt_newprimer2017+=[float(row[2])]
        if row[1]=="2017" and row[6]=="Positive Amplification":
            possum1Pfcrt_newprimer2017+=[float(row[2])]
        if row[1]=="2017" and row[7]=="Negative Amplification":
            negsum1MDR2017+=[float(row[2])]
        if row[1]=="2017" and row[7]=="Positive Amplification":
            possum1MDR2017+=[float(row[2])]
        if row[1]=="2017" and row[8]=="Negative Amplification":
            negsum1MDR_newprimer2017+=[float(row[2])]
        if row[1]=="2017" and row[8]=="Positive Amplification":
            possum1MDR_newprimer2017+=[float(row[2])]
        if row[1]=="2017" and row[9]=="Negative Amplification":
            negsum1CytoB2017+=[float(row[2])]
        if row[1]=="2017" and row[9]=="Positive Amplification":
            possum1CytoB2017+=[float(row[2])]
        if row[1]=="2017" and row[10]=="Negative Amplification":
            negsum1dhps2017+=[float(row[2])]
        if row[1]=="2017" and row[10]=="Positive Amplification":
            possum1dhps2017+=[float(row[2])]
        if row[1]=="2017" and row[11]=="Negative Amplification":
            negsum1dhfr2017+=[float(row[2])]
        if row[1]=="2017" and row[11]=="Positive Amplification":
            possum1dhfr2017+=[float(row[2])]
        if row[1]=="2017" and row[12]=="Negative Amplification":
            negsum1cpmp2017+=[float(row[2])]
        if row[1]=="2017" and row[12]=="Positive Amplification":
            possum1cpmp2017+=[float(row[2])]
        if row[1]=="2017" and row[13]=="Negative Amplification":
            negsum1Pfs472017+=[float(row[2])]
        if row[1]=="2017" and row[13]=="Positive Amplification":
            possum1Pfs472017+=[float(row[2])]###

        if row[1]=="2018" and row[4]=="Negative Amplification":
            negsum1K132018+=[float(row[2])]
        if row[1]=="2018" and row[4]=="Positive Amplification":
            possum1K132018+=[float(row[2])]
        if row[1]=="2018" and row[5]=="Negative Amplification":
            negsum1Pfcrt2018+=[float(row[2])]
        if row[1]=="2018" and row[5]=="Positive Amplification":
            possum1Pfcrt2018+=[float(row[2])]
        if row[1]=="2018" and row[6]=="Negative Amplification":
            negsum1Pfcrt_newprimer2018+=[float(row[2])]
        if row[1]=="2018" and row[6]=="Positive Amplification":
            possum1Pfcrt_newprimer2018+=[float(row[2])]
        if row[1]=="2018" and row[7]=="Negative Amplification":
            negsum1MDR2018+=[float(row[2])]
        if row[1]=="2018" and row[7]=="Positive Amplification":
            possum1MDR2018+=[float(row[2])]
        if row[1]=="2018" and row[8]=="Negative Amplification":
            negsum1MDR_newprimer2018+=[float(row[2])]
        if row[1]=="2018" and row[8]=="Positive Amplification":
            possum1MDR_newprimer2018+=[float(row[2])]
        if row[1]=="2018" and row[9]=="Negative Amplification":
            negsum1CytoB2018+=[float(row[2])]
        if row[1]=="2018" and row[9]=="Positive Amplification":
            possum1CytoB2018+=[float(row[2])]
        if row[1]=="2018" and row[10]=="Negative Amplification":
            negsum1dhps2018+=[float(row[2])]
        if row[1]=="2018" and row[10]=="Positive Amplification":
            possum1dhps2018+=[float(row[2])]
        if row[1]=="2018" and row[11]=="Negative Amplification":
            negsum1dhfr2018+=[float(row[2])]
        if row[1]=="2018" and row[11]=="Positive Amplification":
            possum1dhfr2018+=[float(row[2])]
        if row[1]=="2018" and row[12]=="Negative Amplification":
            negsum1cpmp2018+=[float(row[2])]
        if row[1]=="2018" and row[12]=="Positive Amplification":
            possum1cpmp2018+=[float(row[2])]
        if row[1]=="2018" and row[13]=="Negative Amplification":
            negsum1Pfs472018+=[float(row[2])]
        if row[1]=="2018" and row[13]=="Positive Amplification":
            possum1Pfs472018+=[float(row[2])]###


        if row[1]=="2019" and row[4]=="Negative Amplification":
            negsum1K132019+=[float(row[2])]
        if row[1]=="2019" and row[4]=="Positive Amplification":
            possum1K132019+=[float(row[2])]
        if row[1]=="2019" and row[5]=="Negative Amplification":
            negsum1Pfcrt2019+=[float(row[2])]
        if row[1]=="2019" and row[5]=="Positive Amplification":
            possum1Pfcrt2019+=[float(row[2])]
        if row[1]=="2019" and row[6]=="Negative Amplification":
            negsum1Pfcrt_newprimer2019+=[float(row[2])]
        if row[1]=="2019" and row[6]=="Positive Amplification":
            possum1Pfcrt_newprimer2019+=[float(row[2])]
        if row[1]=="2019" and row[7]=="Negative Amplification":
            negsum1MDR2019+=[float(row[2])]
        if row[1]=="2019" and row[7]=="Positive Amplification":
            possum1MDR2019+=[float(row[2])]
        if row[1]=="2019" and row[8]=="Negative Amplification":
            negsum1MDR_newprimer2019+=[float(row[2])]
        if row[1]=="2019" and row[8]=="Positive Amplification":
            possum1MDR_newprimer2019+=[float(row[2])]
        if row[1]=="2019" and row[9]=="Negative Amplification":
            negsum1CytoB2019+=[float(row[2])]
        if row[1]=="2019" and row[9]=="Positive Amplification":
            possum1CytoB2019+=[float(row[2])]
        if row[1]=="2019" and row[10]=="Negative Amplification":
            negsum1dhps2019+=[float(row[2])]
        if row[1]=="2019" and row[10]=="Positive Amplification":
            possum1dhps2019+=[float(row[2])]
        if row[1]=="2019" and row[11]=="Negative Amplification":
            negsum1dhfr2019+=[float(row[2])]
        if row[1]=="2019" and row[11]=="Positive Amplification":
            possum1dhfr2019+=[float(row[2])]
        if row[1]=="2019" and row[12]=="Negative Amplification":
            negsum1cpmp2019+=[float(row[2])]
        if row[1]=="2019" and row[12]=="Positive Amplification":
            possum1cpmp2019+=[float(row[2])]
        if row[1]=="2019" and row[13]=="Negative Amplification":
            negsum1Pfs472019+=[float(row[2])]
        if row[1]=="2019" and row[13]=="Positive Amplification":
            possum1Pfs472019+=[float(row[2])]###


posmean1K132017=statistics.mean(possum1K132017)
negmean1K132017=statistics.mean(negsum1K132017)

posmean1Pfcrt2017=statistics.mean(possum1Pfcrt2017)
negmean1Pfcrt2017=statistics.mean(negsum1Pfcrt2017)

posmean1Pfcrt_newprimer2017=statistics.mean(possum1Pfcrt_newprimer2017)
negmean1Pfcrt_newprimer2017=statistics.mean(negsum1Pfcrt_newprimer2017)

posmean1MDR2017=statistics.mean(possum1MDR2017)
negmean1MDR2017=statistics.mean(negsum1MDR2017)

posmean1MDR_newprimer2017=statistics.mean(possum1MDR_newprimer2017)
negmean1MDR_newprimer2017=statistics.mean(negsum1MDR_newprimer2017)

posmean1CytoB2017=statistics.mean(possum1CytoB2017)
negmean1CytoB2017=statistics.mean(negsum1CytoB2017)

posmean1dhps2017=statistics.mean(possum1dhps2017)
negmean1dhps2017=statistics.mean(negsum1dhps2017)

posmean1dhfr2017=statistics.mean(possum1dhfr2017)
negmean1dhfr2017=statistics.mean(negsum1dhfr2017)

posmean1cpmp2017=statistics.mean(possum1cpmp2017)
negmean1cpmp2017=statistics.mean(negsum1cpmp2017)

posmean1Pfs472017=statistics.mean(possum1Pfs472017)
negmean1Pfs472017=statistics.mean(negsum1Pfs472017)
#############

posmean1K132018=statistics.mean(possum1K132018)
negmean1K132018=statistics.mean(negsum1K132018)

posmean1Pfcrt2018=statistics.mean(possum1Pfcrt2018)
negmean1Pfcrt2018=statistics.mean(negsum1Pfcrt2018)

posmean1Pfcrt_newprimer2018=statistics.mean(possum1Pfcrt_newprimer2018)
negmean1Pfcrt_newprimer2018=statistics.mean(negsum1Pfcrt_newprimer2018)

posmean1MDR2018=statistics.mean(possum1MDR2018)
negmean1MDR2018=statistics.mean(negsum1MDR2018)

posmean1MDR_newprimer2018=statistics.mean(possum1MDR_newprimer2018)
negmean1MDR_newprimer2018=statistics.mean(negsum1MDR_newprimer2018)

posmean1CytoB2018=statistics.mean(possum1CytoB2018)
negmean1CytoB2018=statistics.mean(negsum1CytoB2018)

posmean1dhps2018=statistics.mean(possum1dhps2018)
negmean1dhps2018=statistics.mean(negsum1dhps2018)

posmean1dhfr2018=statistics.mean(possum1dhfr2018)
negmean1dhfr2018=statistics.mean(negsum1dhfr2018)

posmean1cpmp2018=statistics.mean(possum1cpmp2018)
negmean1cpmp2018=statistics.mean(negsum1cpmp2018)

posmean1Pfs472018=statistics.mean(possum1Pfs472018)
negmean1Pfs472018=statistics.mean(negsum1Pfs472018)

posmean1K132019=statistics.mean(possum1K132019)
negmean1K132019=statistics.mean(negsum1K132019)

#####
posmean1Pfcrt2019=None
negmean1Pfcrt2019=statistics.mean(negsum1Pfcrt2019)

posmean1Pfcrt_newprimer2019=statistics.mean(possum1Pfcrt_newprimer2019)
negmean1Pfcrt_newprimer2019=statistics.mean(negsum1Pfcrt_newprimer2019)

posmean1MDR2019=statistics.mean(possum1MDR2019)
negmean1MDR2019=statistics.mean(negsum1MDR2019)

posmean1MDR_newprimer2019=statistics.mean(possum1MDR_newprimer2019)
negmean1MDR_newprimer2019=statistics.mean(negsum1MDR_newprimer2019)

posmean1CytoB2019=statistics.mean(possum1CytoB2019)
negmean1CytoB2019=statistics.mean(negsum1CytoB2019)

posmean1dhps2019=statistics.mean(possum1dhps2019)
negmean1dhps2019=statistics.mean(negsum1dhps2019)

posmean1dhfr2019=statistics.mean(possum1dhfr2019)
negmean1dhfr2019=statistics.mean(negsum1dhfr2019)

posmean1cpmp2019=statistics.mean(possum1cpmp2019)
negmean1cpmp2019=statistics.mean(negsum1cpmp2019)

posmean1Pfs472019=statistics.mean(possum1Pfs472019)
negmean1Pfs472019=statistics.mean(negsum1Pfs472019)




posmean1K132017=statistics.mean(possum1K132017)
negmean1K132017=statistics.mean(negsum1K132017)

posmean1Pfcrt2017=statistics.mean(possum1Pfcrt2017)
negmean1Pfcrt2017=statistics.mean(negsum1Pfcrt2017)

posmean1Pfcrt_newprimer2017=statistics.mean(possum1Pfcrt_newprimer2017)
negmean1Pfcrt_newprimer2017=statistics.mean(negsum1Pfcrt_newprimer2017)

posmean1MDR2017=statistics.mean(possum1MDR2017)
negmean1MDR2017=statistics.mean(negsum1MDR2017)

posmean1MDR_newprimer2017=statistics.mean(possum1MDR_newprimer2017)
negmean1MDR_newprimer2017=statistics.mean(negsum1MDR_newprimer2017)

posmean1CytoB2017=statistics.mean(possum1CytoB2017)
negmean1CytoB2017=statistics.mean(negsum1CytoB2017)

posmean1dhps2017=statistics.mean(possum1dhps2017)
negmean1dhps2017=statistics.mean(negsum1dhps2017)

posmean1dhfr2017=statistics.mean(possum1dhfr2017)
negmean1dhfr2017=statistics.mean(negsum1dhfr2017)

posmean1cpmp2017=statistics.mean(possum1cpmp2017)
negmean1cpmp2017=statistics.mean(negsum1cpmp2017)

posmean1Pfs472017=statistics.mean(possum1Pfs472017)
negmean1Pfs472017=statistics.mean(negsum1Pfs472017)

positive2017=[possum1K132017,possum1Pfcrt2017, \
possum1Pfcrt_newprimer2017,possum1MDR2017, possum1MDR_newprimer2017, \
possum1CytoB2017,possum1dhps2017,possum1dhfr2017, \
possum1cpmp2017,possum1Pfs472017]

negative2017=[negsum1K132017,negsum1Pfcrt2017,negsum1Pfcrt_newprimer2017,negsum1MDR2017,negsum1MDR_newprimer2017,negsum1CytoB2017,negsum1dhps2017,negsum1dhfr2017,negsum1cpmp2017,negsum1Pfs472017]

positive12017=[statistics.mean(possum1K132017),statistics.mean(possum1Pfcrt2017), \
statistics.mean(possum1Pfcrt_newprimer2017),statistics.mean(possum1MDR2017),statistics.mean(possum1MDR_newprimer2017), \
statistics.mean(possum1CytoB2017),statistics.mean(possum1dhps2017),statistics.mean(possum1dhfr2017), \
statistics.mean(possum1cpmp2017),statistics.mean(possum1Pfs472017)]

negative12017=[negmean1K132017,negmean1Pfcrt2017,negmean1Pfcrt_newprimer2017,negmean1MDR2017,negmean1MDR_newprimer2017,negmean1CytoB2017,negmean1dhps2017,negmean1dhfr2017,negmean1cpmp2017,negmean1Pfs472017]


CT = np.random.rand((10*18))*40
Parasitemia =   np.random.rand((10*18))*3000

#t = np.arange(0.01, 10.0, 0.01)
#data1 = np.exp(t)
#data2 = np.sin(2 * np.pi * t)

fig, ax1 = plt.subplots()

negative2017log=[]
for x in negative2017:
    negative2017log+=[-1.49*np.log(x)+37.641]

combined2017=positive2017+negative2017

#print(negative2017log)

genes = [
    'K13', 'Pfcrt', 'Pfcrt_newprimer', 'MDR', 'MDR_newprimer', 'CytoB',
    'dhps', 'dhfr', 'cpmp', 'Pfs47'
]

#print(len(possum1cpmp2017))
#print(len(negsum1cpmp2017))

ax1.set_xlabel('TES 2017')
ax1.set_ylabel('Cycle Threshold - Falciparum', color='tab:gray')
res1 = ax1.boxplot(
    positive2017, positions = np.arange(10)-0.25, widths=0.4,
    patch_artist=True,
)
for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(res1[element], color='k')

for patch in res1['boxes']:
    patch.set_facecolor('tab:blue')



ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel('Parasitemia (p/\u03BCL)', color='tab:gray')
res2 = ax1.boxplot(
    negative2017, positions = np.arange(10)+0.25, widths=0.4,
    patch_artist=True,
)
##from https://stackoverflow.com/a/41997865/2454357
for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(res2[element], color='k')

for patch in res2['boxes']:
    patch.set_facecolor('tab:orange')

ax1.set_xlim([-0.55, 11.55])
ax1.set_xticks(np.arange(10))
ax1.set_xticklabels(genes)
plt.legend([res1["boxes"][0], res2["boxes"][0]], ["Positive Amplification", "Negative Amplification"], loc='upper right')

print(min(min(combined2017)))
print(max(max(combined2017)))

#print(-1.49*np.log(min(min(combined2017)))+37.641)
#print(-1.49*np.log(max(max(combined2017)))+37.641)
min2017=min(min(combined2017))
max2017=max(max(combined2017))
#print(math.exp((-(min2017)+39.258)/1.722))
#print(math.exp((-(max2017)+39.258)/1.722))
#print([math.exp((-(min2017)+39.258)/1.722),math.exp((-(max2017)+39.258)/1.722)])
#ax2.set_ylim([math.exp((-(min2017)+39.258)/1.722),math.exp((-(max2017)+39.258)/1.722)])
#ax2.set_ylim([-1.49*np.log(max(max(combined2017)))+37.641,-1.49*np.log(min(min(combined2017)))+37.641])
#ax2.set_yticks(np.arange(100))

yticks2017=[]
yticks20172=[]
for x in range(int(min(min(combined2017)))-2,int(max(max(combined2017)))+2,2):
    yticks2017+=[x]
    yticks20172+=[str((math.exp((-x+39.258)/1.722)))]

ax1.set_ylim(int(min(min(combined2017)))-2,int(max(max(combined2017)))+2)
ax1.set_yticks(yticks2017)

ax2.set_ylim(int(min(min(combined2017)))-2,int(max(max(combined2017)))+2)
ax2.set_yticks(yticks2017)
ax2.set_yticklabels(yticks20172)


fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig("testplot1.png")
plt.show()
plt.close()

positive2018=[possum1K132018,possum1Pfcrt2018, \
possum1Pfcrt_newprimer2018,possum1MDR2018, possum1MDR_newprimer2018, \
possum1CytoB2018,possum1dhps2017,possum1dhfr2018, \
possum1cpmp2018,possum1Pfs472018]

negative2018=[negsum1K132018,negsum1Pfcrt2018,negsum1Pfcrt_newprimer2018,negsum1MDR2018,\
negsum1MDR_newprimer2018,negsum1CytoB2018,negsum1dhps2018,negsum1dhfr2018,negsum1cpmp2018,negsum1Pfs472018]

combined2018=positive2018+negative2018

#print(positive2018)

fig, ax3 = plt.subplots()

ax3.set_xlabel('TES 2018')
ax3.set_ylabel('Cycle Threshold - Falciparum', color='tab:gray')
res1 = ax3.boxplot(
    positive2018, positions = np.arange(10)-0.25, widths=0.4,
    patch_artist=True,
)
for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(res1[element], color='k')

for patch in res1['boxes']:
    patch.set_facecolor('tab:blue')



ax4 = ax3.twinx()  # instantiate a second axes that shares the same x-axis
ax4.set_ylabel('Parasitemia (p/\u03BCL)', color='tab:gray')
res2 = ax3.boxplot(
    negative2018, positions = np.arange(10)+0.25, widths=0.4,
    patch_artist=True,
)
##from https://stackoverflow.com/a/41997865/2454357
for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(res2[element], color='k')

for patch in res2['boxes']:
    patch.set_facecolor('tab:orange')

ax3.set_xlim([-0.55, 11.55])
ax3.set_xticks(np.arange(10))
ax3.set_xticklabels(genes)
plt.legend([res1["boxes"][0], res2["boxes"][0]], ["Positive Amplification", "Negative Amplification"], loc='upper right')

#print(min(min(combined2017)))
#print(max(max(combined2017)))

#print(-1.49*np.log(min(min(combined2017)))+37.641)
#print(-1.49*np.log(max(max(combined2017)))+37.641)
#ax4.set_ylim([math.exp(-((min(min(combined2018)))-39.258)/1.722),math.exp(-((max(max(combined2018)))-39.258)/1.722)])
yticks2018=[]
yticks20182=[]
for x in range(int(min(min(combined2018)))-2,int(max(max(combined2018)))+2,2):
    yticks2018+=[x]
    yticks20182+=[str((math.exp((-x+39.258)/1.722)))]

ax3.set_ylim(int(min(min(combined2018)))-2,int(max(max(combined2018)))+2)
ax3.set_yticks(yticks2018)

ax4.set_ylim(int(min(min(combined2018)))-2,int(max(max(combined2018)))+2)
ax4.set_yticks(yticks2018)
ax4.set_yticklabels(yticks20182)
#ax2.set_ylim([-1.49*np.log(max(max(combined2017)))+37.641,-1.49*np.log(min(min(combined2017)))+37.641])
#ax2.set_yticks(np.arange(100))


fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig("testplot2.png")
plt.show()
plt.close()

positive2019=[possum1K132019,possum1Pfcrt2019, \
possum1Pfcrt_newprimer2019,possum1MDR2019, possum1MDR_newprimer2019, \
possum1CytoB2019,possum1dhps2019,possum1dhfr2019, \
possum1cpmp2019,possum1Pfs472019]

negative2019=[negsum1K132019,negsum1Pfcrt2019,negsum1Pfcrt_newprimer2019,negsum1MDR2019,\
negsum1MDR_newprimer2019,negsum1CytoB2019,negsum1dhps2019,negsum1dhfr2019,negsum1cpmp2019,negsum1Pfs472019]

combined2019=positive2019+negative2019

#print(positive2018)

fig, ax5 = plt.subplots()

ax5.set_xlabel('TES 2019')
ax5.set_ylabel('Cycle Threshold - Falciparum', color='tab:gray')
res1 = ax5.boxplot(
    positive2019, positions = np.arange(10)-0.25, widths=0.4,
    patch_artist=True,
)
for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(res1[element], color='k')

for patch in res1['boxes']:
    patch.set_facecolor('tab:blue')



ax6 = ax5.twinx()  # instantiate a second axes that shares the same x-axis
ax6.set_ylabel('Parasitemia (p/\u03BCL)', color='tab:gray')
res2 = ax5.boxplot(
    negative2019, positions = np.arange(10)+0.25, widths=0.4,
    patch_artist=True,
)
##from https://stackoverflow.com/a/41997865/2454357
for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(res2[element], color='k')

for patch in res2['boxes']:
    patch.set_facecolor('tab:orange')

ax5.set_xlim([-0.55, 11.55])
ax5.set_xticks(np.arange(10))
ax5.set_xticklabels(genes)
plt.legend([res1["boxes"][0], res2["boxes"][0]], ["Positive Amplification", "Negative Amplification"], loc='upper right')

#print(min(min(combined2017)))
#print(max(max(combined2017)))

#print(-1.49*np.log(min(min(combined2017)))+37.641)
#print(-1.49*np.log(max(max(combined2017)))+37.641)
#ax6.set_ylim([math.exp(-((min(min(combined2019)))-39.258)/1.722),math.exp(-((max(max(combined2019)))-39.258)/1.722)])
yticks2019=[]
yticks20192=[]
for x in range(22,int(max(max(combined2019)))+2,2):
    yticks2019+=[x]
    yticks20192+=[str((math.exp((-x+39.258)/1.722)))]

ax5.set_ylim(22,int(max(max(combined2019)))+2)
ax5.set_yticks(yticks2019)

ax6.set_ylim(22,int(max(max(combined2019)))+2)
ax6.set_yticks(yticks2019)
ax6.set_yticklabels(yticks20192)
#ax2.set_ylim([-1.49*np.log(max(max(combined2017)))+37.641,-1.49*np.log(min(min(combined2017)))+37.641])
#ax2.set_yticks(np.arange(100))


fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig("testplot3.png")
plt.show()
plt.close()
