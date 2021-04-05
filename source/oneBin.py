from unit_const import *
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
import os
######################################
neutrinotype = ['e','a','x']
# Parameters for Hydro-simulation from Albino Perego
nEnergyBins = 30
numberofneutrinotype = 3
maxE = 100
minE = 0.01
E = np.linspace(minE,maxE,nEnergyBins)  # original bins
dE = np.array([np.exp(np.log(minE*M)+(i+0.5)*((np.log(maxE*M)-np.log(minE*M))/(nEnergyBins-1)))-np.exp(np.log(minE*M)+(i-0.5)*((np.log(maxE*M)-np.log(minE*M))/(nEnergyBins-1))) for i in range(nEnergyBins)])
dE/=M #MeV
##################################
def getNewE(Emin,Emax,nNew):
    if nNew ==1:
        Enew =np.array([np.exp(np.log(Emin * M) + (np.log(Emax * M) - np.log(Emin* M))/2)])
        dEnew = np.array([0.])
    else:
        Enew = np.array([np.exp(np.log(Emin * M) + i * ((np.log(Emax * M) - np.log(Emin * M)) / (nNew - 1))) for i in range(nNew)])
        dEnew = np.array([np.exp(np.log(Emin*M)+(i+0.5)*((np.log(Emax*M)-np.log(Emin*M))/(nNew-1)))-np.exp(np.log(Emin*M)+(i-0.5)*((np.log(Emax*M)-np.log(Emin*M))/(nNew-1))) for i in range(nNew)])
    dEnew/=M
    Enew /=M
    EnewCgs=Enew*MeV
    ECgs =E*MeV #erg
    dECgs = dE*MeV #erg
    dEnewCgs = dEnew*MeV #erg
    return Enew,dEnew,EnewCgs,ECgs,dECgs,dEnewCgs

def onebin(dir,nNew,ntestpoints):
    Enew,dEnew,EnewCgs,ECgs,dECgs,dEnewCgs=getNewE(minE,maxE,nNew)
    outFilePotentiale = open(dir + '/' + str(nNew) + 'NE/'  +'vnue-nsnssc03-000Mo-022deg.dat',
                            'w')  # plot linear and cubic spline
    outFileDensitye = open(dir +'/' + str(nNew) + 'NE/' +'nnue-nsnssc03-000Mo-022deg.dat',
                          'w')  # plot linear and cubic spline
    outFilePotentiala = open(dir + '/' + str(nNew) + 'NE/'+'vnua-nsnssc03-000Mo-022deg.dat',
                         'w')  # plot linear and cubic spline
    outFileDensitya = open(dir +'/' + str(nNew) + 'NE/' + 'nnua-nsnssc03-000Mo-022deg.dat',
                      'w')  # plot linear and cubic spline
    outFilePotentialx = open(dir + '/' + str(nNew) + 'NE/'  +'vnux-nsnssc03-000Mo-022deg.dat',
                        'w')  # plot linear and cubic spline
    outFileDensityx = open(dir +'/' + str(nNew) + 'NE/' + 'nnux-nsnssc03-000Mo-022deg.dat',
                      'w')  # plot linear and cubic spline
    pe = [0] * nEnergyBins
    pnewe = [0]*numberofneutrinotype*nNew
    de = [0] * nEnergyBins
    dnewe = [0]*numberofneutrinotype*nNew
    pa = [0] * nEnergyBins
    pnewa = [0]*numberofneutrinotype*nNew
    da = [0] * nEnergyBins
    dnewa = [0]*numberofneutrinotype*nNew
    px = [0] * nEnergyBins
    pnewx = [0]*numberofneutrinotype*nNew
    dx = [0] * nEnergyBins
    dnewx = [0]*numberofneutrinotype*nNew

    for n in range(4,ntestpoints):  # all points in the test trajectory
        Rcore=0. #from the core, use this for neutrino caputure rate extension
        Rs=0. #from the stating point, usually used for neutrino oscillation
        totalPotential=[0.]*numberofneutrinotype
        totalDensity=[0.]*numberofneutrinotype
        newtotalPotential=[0.]*numberofneutrinotype
        newtotalDensity=[0.]*numberofneutrinotype
        print(n)
        for j in range(0, nEnergyBins):
                inFilePe = open(dir +'/vnue-'+'{:02}'.format(1+j)+'-nsnssc03-000Mo-022deg.dat', 'r')
                inFileNe = open(dir +'/nnue-'+'{:02}'.format(1+j)+'-nsnssc03-000Mo-022deg.dat', 'r')
                inFilePa = open(dir +'/vnua-'+'{:02}'.format(1+j)+'-nsnssc03-000Mo-022deg.dat', 'r')
                inFileNa = open(dir +'/nnua-'+'{:02}'.format(1+j)+'-nsnssc03-000Mo-022deg.dat', 'r')
                inFilePx = open(dir +'/vnux-'+'{:02}'.format(1+j)+'-nsnssc03-000Mo-022deg.dat', 'r')
                inFileNx = open(dir +'/nnux-'+'{:02}'.format(1+j)+'-nsnssc03-000Mo-022deg.dat', 'r')
                linesPe = inFilePe.readlines()  # yzhu: START FROM 0
                linesNe = inFileNe.readlines()  # yzhu: START FROM 0
                linesPa = inFilePa.readlines()  # yzhu: START FROM 0
                linesNa = inFileNa.readlines()  # yzhu: START FROM 0
                linesPx = inFilePx.readlines()  # yzhu: START FROM 0
                linesNx = inFileNx.readlines()  # yzhu: START FROM 0
                inFileNe.close()
                inFilePe.close()
                inFileNa.close()
                inFilePa.close()
                inFileNx.close()
                inFilePx.close()
                coordinatesLinesPe = linesPe[n]
                coordinatesLinesNe = linesNe[n]
                coordinatesLinesPa = linesPa[n]
                coordinatesLinesNa = linesNa[n]
                coordinatesLinesPx = linesPx[n]
                coordinatesLinesNx = linesNx[n]
                pointsCoordinatesPe = coordinatesLinesPe.split()
                pointsCoordinatesNe = coordinatesLinesNe.split()
                pointsCoordinatesPa = coordinatesLinesPa.split()
                pointsCoordinatesNa = coordinatesLinesNa.split()
                pointsCoordinatesPx = coordinatesLinesPx.split()
                pointsCoordinatesNx = coordinatesLinesNx.split()
                pe[j] = float(pointsCoordinatesPe[1])/dECgs[j]
                de[j] = float(pointsCoordinatesNe[1])/dECgs[j]
                pa[j] = float(pointsCoordinatesPa[1])/dECgs[j]
                da[j] = float(pointsCoordinatesNa[1])/dECgs[j]
                px[j] = float(pointsCoordinatesPx[1])/dECgs[j]
                dx[j] = float(pointsCoordinatesNx[1])/dECgs[j]
                totalPotential[0]=totalPotential[0]+pe[j]*dECgs[j] #/dE
                totalDensity[0]=totalDensity[0]+de[j]*dECgs[j]
                totalPotential[1]=totalPotential[1]+pa[j]*dECgs[j] #/dE
                totalDensity[1]=totalDensity[1]+da[j]*dECgs[j]
                totalPotential[2]=totalPotential[2]+px[j]*dECgs[j] #/dE
                totalDensity[2]=totalDensity[2]+dx[j]*dECgs[j]
        outFileDensitye.write(str(float(pointsCoordinatesNx[0])) + "\t" +str(float(totalDensity[0]))+"\n")
        outFilePotentiale.write(str(float(pointsCoordinatesNx[0])) + "\t" +str(float(totalPotential[0]))+"\n")
        outFileDensitya.write(str(float(pointsCoordinatesNx[0])) + "\t" +str(float(totalDensity[1]))+"\n")
        outFilePotentiala.write(str(float(pointsCoordinatesNx[0])) + "\t" +str(float(totalPotential[1]))+"\n")
        outFileDensityx.write(str(float(pointsCoordinatesNx[0])) + "\t" +str(float(totalDensity[2]))+"\n")
        outFilePotentialx.write(str(float(pointsCoordinatesNx[0])) + "\t" +str(float(totalPotential[2]))+"\n")
    outFilePotentiale.close()
    outFileDensitye.close()
    outFilePotentiala.close()
    outFileDensitya.close()
    outFilePotentialx.close()
    outFileDensityx.close()
#######################################################################

dir = "/Users/yzhu14/Research/forBrettPaper/potentials-22deg-RunAB-sc/"
#potentials-RunBA-025deg-scat-04/"
newElist=[1]

########################################################################
for nNew in newElist:
        ntestpointsN = sum(1 for line in open(dir +'/nnux-'+'01'+'-nsnssc03-000Mo-022deg.dat'))
        ntestpointsP = sum(1 for line in open(dir +'/vnux-'+'01'+'-nsnssc03-000Mo-022deg.dat'))
        ntestpoints =min(ntestpointsN,ntestpointsP)
        print(ntestpoints)
        onebin(dir,nNew,ntestpoints)
