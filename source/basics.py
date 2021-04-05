__author__ = 'yzhu14'
import numpy as np
from unit_const import *
from scipy.interpolate import interp1d
import numpy as np
import math
import sys
import csv
import itertools as IT
from numpy import linalg as LA



#read in interpolated density with nNew bins from albino's density, interpolate into more x of newPnum
#for brett's input, Jan 2018
nNew= 30
nCommentLines=4
neutrinotype= ['a','e','x']
potype=['nnu','vnu','ye','rho']
dir= "/Users/yzhu14/Research/forBrettPaper/potentials-RunBA-025deg-scat-04/"
for pt in potype:
    if pt=="ye" or pt=="rho":
        inFileyr=dir+pt+"-nsnssc04-000Mo-025deg.dat"
        print(inFileyr)
        oldPnumyr=sum(1 for line in open(inFileyr))
        inFileYR=open(inFileyr,"r")
        linesyr=inFileYR.readlines()
        R1=[0.]*(oldPnumyr-nCommentLines)
        YR=[0.]*(oldPnumyr-nCommentLines)
        print(oldPnumyr)
        for nt in range(nCommentLines,oldPnumyr):
            coordinatesLinesyr=linesyr[nt]
            pointsCoordinatesyr=coordinatesLinesyr.split()
            R1[nt-nCommentLines]=float(pointsCoordinatesyr[0])
            YR[nt-nCommentLines]=float(pointsCoordinatesyr[1])
            print(nt,float(pointsCoordinatesyr[0]),float(pointsCoordinatesyr[1]))
            #if nt==0:
            #print("0",float(pointsCoordinatesD[1+k+it*nNew]))
        f1=interp1d(R1,YR)
        outFile= open(dir+pt+"-nsnssc04.dat",'w')
        dR=float(R1[-1]-R1[0])/(oldPnum-nCommentLines)
        dR/=10
        r=float(R1[0])
        while r<R1[-2]:
            r+=dR
            outFile.write(str(r)+"\t"+str(f1(r))+"\n")
        outFile.close()
    else:
        for it in neutrinotype:
            for k in range(0,nNew):
                inFilename=dir+pt+str(it)+"-"+str(format(k+1,'02d'))+"-nsnssc04-000Mo-025deg.dat"
                oldPnum=sum(1 for line in open(inFilename))
                inFile=open(inFilename,"r")
                linesD=inFile.readlines()
                RR=[0.]*(oldPnum-nCommentLines)
                Density=[0.]*(oldPnum-nCommentLines)
                for nt in range(nCommentLines,oldPnum):
                    coordinatesLinesD=linesD[nt]
                    pointsCoordinatesD=coordinatesLinesD.split()
                    RR[nt-nCommentLines]=float(pointsCoordinatesD[0])
                    Density[nt-nCommentLines]=float(pointsCoordinatesD[1])
                    #if nt==0:
                    #print("0",float(pointsCoordinatesD[1+k+it*nNew]))
                f=interp1d(RR,Density)
                outFile= open(dir+pt+str(it)+"-"+str(format(k+1,'02d'))+"-nsnssc04.dat",'w')
                dR=float(RR[-1]-RR[0])/(oldPnum-nCommentLines)
                dR/=10
                r=float(RR[0])
                while r<RR[-2]:
                    r+=dR
                    outFile.write(str(r)+"\t"+str(f(r))+"\n")
                outFile.close()
