import matplotlib.pyplot as plt
import numpy as np
import os

nFlavor = 3
nE = 30
fitcoe =np.array([0.,0.,0.,0.,0.])*nFlavor

#filename = ['vnue-01-nsnssc03-000Mo-022deg','vnua-01-nsnssc03-000Mo-022deg','vnux-01-nsnssc03-000Mo-022deg']
dir = "/Users/yzhu14/Research/forBrettPaper/potentials-22deg-RunAB-sc/"
#dir = '/Users/yzhu14/Research/forBrettPaper/potentials-RunBA-025deg-scat-04/'
outfile = open(dir +'p.dat','w')
i=0
filename = ['pe','pa','px']
for file in filename:
    file1 = dir + file+".dat"
    data = np.loadtxt(file1,skiprows=0)
    x4 = data[:,1]
    x3 = data[:,2]
    x2 = data[:,3]
    x1 = data[:,4]
    x0 = data[:,5]
    #if i==0:
    #print(x4,x3,x2,x1,x0)
    print(np.sum(x4),np.sum(x3),np.sum(x2),np.sum(x1),np.sum(x0))
    #print(i,np.array([np.sum(x4),np.sum(x3),np.sum(x2),np.sum(x1),np.sum(x0)]))
    i+=1
for x1 in np.linspace(300, 1300, 1e3):#
        plt.semilogy(x1, 1/(), 'b+')
