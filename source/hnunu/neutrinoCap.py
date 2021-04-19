from scipy.interpolate import interp1d
import numpy as np
import math
import sys
from scipy.optimize import curve_fit
from hnunu.unit_const import *
from hnunu.albinoBasics import *

def func1(r,a,b):
    return a+b/r**2

def func2(r,a,b):
    return a+b/r**2

# extend the neutrino capture rate
def exNuCap(Nucapfile,MaxR,nNew):
    Enew,dEnew,EnewCgs,ECgs,dECgs,dEnewCgs=getNewE(nNew)
    print(Nucapfile)
    dR=0.0
    inFileNC=open(Nucapfile+'nucap.txt','r')
    #    inFilePotential = open(dir + '/' + str(nNew) + 'NE/'  +whichFileP+ str(nf) + note + '.txt','r')  # plot linear and cubic spline
    nline=sum(1 for line in open(Nucapfile + 'nucap.txt'))
    lineNC = inFileNC.readlines()
    rOld=[0.]*(nline-1)
    nrate=[0.]*(nline-1)
    nrateOsc=[0.]*(nline-1)
    prate=[0.]*(nline-1)
    prateOsc=[0.]*(nline-1)
    outfilenprate=open(Nucapfile + 'nprate.txt','w')
    outfilenprateOsc=open(Nucapfile + 'nprateOsc.txt','w')
    for i in range(1,nline):
        coordinatesNC=lineNC[i]
        pointNC=coordinatesNC.split()
        rOld[i-1]=float(pointNC[0])
        nrate[i-1]=float(pointNC[1])
        nrateOsc[i-1]=float(pointNC[2])
        prate[i-1]=float(pointNC[3])
        prateOsc[i-1]=float(pointNC[4])
        outfilenprate.write(str(rOld[i-1])+"\t"+str(nrate[i-1])+"\t"+str(prate[i-1])+"\n")
        outfilenprateOsc.write(str(rOld[i-1])+"\t"+str(nrateOsc[i-1])+"\t"+str(prateOsc[i-1])+"\n")
    print("rold",min(rOld)/1e5,max(rOld)/1e5)
    inFileNC.close()
    nrateSample=nrate[-2:]
    prateSample=prate[-2:]
    rOldSample=rOld[-2:]
    print(rOldSample,prateSample,nrateSample)
    p0=[1.8,1.3e15]
    n0=[2.4,9.2e14]               #initial guess is necessary
    poptn, pcovn = curve_fit(func1, rOldSample, nrateSample,p0=n0)
    poptp, pcovp = curve_fit(func1, rOldSample, prateSample,p0=p0)
    poptnOsc, pcovnOsc = curve_fit(func2, rOld, nrateOsc,p0=p0)
    poptpOsc, pcovpOsc = curve_fit(func2, rOld, prateOsc,p0=p0)
    dR=rOld[-1]-rOld[-2]
    print(rOld[-1],rOld[-2],dR)
    #continuting writing file of neutrino capture rate
    print(poptn, pcovn,poptp, pcovp)
    print(int((MaxR-rOld[-1])/dR),poptnOsc, pcovnOsc)
    for j in range(0,int((MaxR-rOld[-1])/dR)):
        R=rOld[-1]+dR*j
        outfilenprate.write(str(R)+"\t"+str(func2(R,*poptn))+"\t"+str(func2(R,*poptp))+"\n")
        outfilenprateOsc.write(str(R)+"\t"+str(func2(R,*poptnOsc))+"\t"+str(func2(R,*poptpOsc))+"\n")
        #print(R)
    print(MaxR,rOld[-1],dR,((MaxR-rOld[-1])/dR))
    outfilenprate.close()
    outfilenprateOsc.close()