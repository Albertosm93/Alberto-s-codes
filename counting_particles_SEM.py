# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 09:54:42 2016

@author: py15asm
"""

import os, os.path
import glob
import numpy as np
from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import pylab
import math

#path of the two csv files, that should be called washed and nwashed

path='C:\\Users\\py15asm\\PHD INP\\Eng SEM\\Images of the filters to count\\161024 analysed filter\\areas.csv'
title='APS vs SEM counting 161005 afternoon'
nbin=20

#sampled air in /cm3
air=2610000

#surfice of the filter
surf=3.1415*(47000/2)**2


rowdata=np.genfromtxt(path,delimiter=',',skip_header=1,dtype=np.float64)


lowm=np.array(rowdata[:,1])
medm=np.array(rowdata[:,6])

dlowm=np.sqrt((4/3.1415)*lowm)
dmedm=np.sqrt((4/3.1415)*medm)
alowm=150.6*131.8
amedm=37.5*32.8



lowmhist,lowmedges=np.histogram(dlowm,bins=nbin,range=(0.1,12))
medmhist,medmedges=np.histogram(dmedm,bins=nbin,range=(0.1,12))

#number of images analysed
    #lowmag
imageslowm=0
counting=np.array(rowdata[:,0])
for i in range(1,len(counting)):
    
    if counting[i]>counting[i-1] or math.isnan(counting[i]):
        continue
    imageslowm=imageslowm+1
    print imageslowm
    
imagesmedm=0
counting=np.array(rowdata[:,5])
for i in range(1,len(counting)):
    
    if counting[i]>counting[i-1] or math.isnan(counting[i]):
        continue
    imagesmedm=imagesmedm+1
    print imagesmedm
histo=lowmhist+medmhist
lowmreal=lowmhist/(imageslowm*alowm)   
medmreal=medmhist/(imagesmedm*amedm)  

#calculating dLogDp


dlogDp=np.zeros(nbin)
bincenter=np.zeros(nbin)
for i in range(0,len(lowmedges)-1):
    bincenter[i]=(lowmedges[i]+lowmedges[i+1])/2
    dlogDp[i]=math.log10(lowmedges[i+1])-math.log10(lowmedges[i])

hist=(lowmreal+medmreal)*(surf/air)


#error
total=np.sum(histo)
fraction=histo/float(total)
relerror=1.96*np.sqrt(fraction*(1-fraction)/total)
error=relerror*hist/dlogDp


hist=hist/dlogDp

#SEM
#SMPS
SMPSlocation='W:\\Formatted Correctly\\161027\\SMPS\\161027 1144.csv'
dataSMPS=np.genfromtxt(SMPSlocation,delimiter=',',skip_header=26,dtype=np.float64)
    
SMPS=dataSMPS[0:-1,9:119]
SMPSdavg=SMPS.mean(axis=0)
SMPS,=plt.plot(SMPSbincenter,SMPSdavg,label='SMPS')

#APS
APSfile='W:\\Formatted Correctly\\161005\APS\\161005 1102.csv'
dataAPS=np.genfromtxt(APSfile,delimiter=',',skip_header=7,dtype=np.float64)

#number of intervals you want to divide the APS (always 1)
n=1

APSdlogDp=[0.015497598,0.031669268,0.030905778,0.031440731,0.031123233,0.031282722,0.031265457,0.031589861,0.031250977,0.031251443,0.031106434,0.03120896,0.031146818,0.031267074,0.031510766,0.031014078,0.031231237,0.031261211,0.03135914,0.031291186,0.031280359,0.031125858,0.031350935,0.031232579,0.031276815,0.031166742,0.031303694,0.031267966,0.031200335,0.031199862,0.031332206,0.031206334,0.031280396,0.031194313,0.031257694,0.031292358,0.031234397,0.031266032,0.031249935,0.031236754,0.031261053,0.031386478,0.031106434,0.03120896,0.031146818,0.031267074,0.031510766,0.031014078,0.031231237,0.031261211,0.03135914]
APSbincenter=np.array([0.5325,0.5625,0.6045,0.6495,0.698,0.75,0.806,0.8665,0.9315,1.001,1.0755,1.1555,1.2415,1.334,1.434,1.541,1.6555,1.779,1.912,2.055,2.2085,2.373,2.55,2.7405,2.945,3.1645,3.4005,3.6545,3.927,4.2195,4.5345,4.873,5.2365,5.627,6.0465,6.498,6.983,7.504,8.064,8.6655,9.312,10.0085,10.755,11.555,12.415,13.34,14.34,15.41,16.555,17.79,19.12])

nint=dataAPS[-1,0]
import math
u=math.modf(nint/n)
length=u[1]

plt.xscale('log')
#plt.xscale('log')
plt.errorbar(bincenter,hist,yerr=error,linestyle="",marker='o',markersize=4,label='SEM, two resolutions')

for i in range(0,len(histo)):
    plt.annotate(histo[i],(bincenter[i],0.5*hist[i]))
    
for i in range(0,n):     #this doesn't includes the end, that is why I wrote n rather than n-1
    APS=dataAPS[i*length:(i+1)*length,5:56]
    APSavg=APS.mean(axis=0)
    APSdavg=np.divide(APSavg,APSdlogDp)
    APSdavg=APSdavg[:43]
    APSbincenter=APSbincenter[:43]
    plt.errorbar(APSbincenter,APSdavg,label='APS')
    plt.xlabel('size ($\mu m$)')
    plt.ylabel('dN/dlogDp (cm$^{-3}$)')
    pylab.xlim([0.3,12])
    pylab.ylim([0.001,max(hist)+10])
    

#Adding a second line with the one resolution infomration
path2='C:\\Users\\py15asm\\PHD INP\\Eng SEM\\Images of the filters to count\\161024 analysed filter\\areas1resolution.csv'
rowdata2=np.genfromtxt(path2,delimiter=',',skip_header=1,dtype=np.float64)
lowm2=np.array(rowdata2[:,1])


dlowm2=np.sqrt((4/3.1415)*lowm2)

lowmhist2,lowmedges2=np.histogram(dlowm2,bins=nbin,range=(0.1,12))


#number of images analysed
    #lowmag
imageslowm2=0
counting2=np.array(rowdata2[:,0])
for i in range(1,len(counting2)):
    
    if counting2[i]>counting2[i-1] or math.isnan(counting2[i]):
        continue
    imageslowm2=imageslowm2+1
    print imageslowm2
    

histo2=lowmhist2
lowmreal2=lowmhist2/(imageslowm2*alowm)   
#calculating dLogDp
dlogDp2=np.zeros(nbin)
bincenter2=np.zeros(nbin)
for i in range(0,len(lowmedges2)-1):
    bincenter2[i]=(lowmedges2[i]+lowmedges2[i+1])/2
    dlogDp2[i]=math.log10(lowmedges2[i+1])-math.log10(lowmedges2[i])

hist2=(lowmreal2)*(surf/air)
hist2=hist2/dlogDp2

#error
total2=np.sum(histo2)
fraction2=histo2/float(total2)
relerror2=1.96*np.sqrt(fraction2*(1-fraction2)/total2)
error2=relerror2*hist2/dlogDp2

plt.errorbar(bincenter2,hist2,yerr=error2,linestyle="",marker='o',markersize=4,label='SEM, one resolution')

#plot        
 #marker='o',width=0.5,,log=True
plt.title(title)
plt.legend()
pylab.savefig(title)

plt.show()
plt.clf
