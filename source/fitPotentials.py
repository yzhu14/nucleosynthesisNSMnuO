import matplotlib.pyplot as plt
import numpy as np
import os
import math
n=4
nFlavor = 3
z= np.array([0.])*30
#filename= ['vnue-nsnssc03-000Mo-022deg.dat','vnua-nsnssc04-000Mo-025deg.dat','vnux-nsnssc04-000Mo-025deg.dat']
#filename= 'rho-nsnssc04-000Mo-025deg.dat'
#filename = ['nnue-01-nsnssc03-000Mo-022deg','nnua-01-nsnssc03-000Mo-022deg','nnux-01-nsnssc03-000Mo-022deg']
#filename = ['vnue-nsnssc03-000Mo-022deg','vnua-nsnssc03-000Mo-022deg','vnux-nsnssc03-000Mo-022deg']
#dir = '/Users/yzhu14/Research/forBrettPaper/potentials-RunBA-025deg-scat-04/'
#fitcoe =np.array([0.,0.,0.,0.,0.])*nFlavor
dir = '/Users/yzhu14/Research/forBrettPaper/potentials-RunBA-025deg-scat-04/1NE/'
dir = "/Users/yzhu14/Research/forBrettPaper/potentials-22deg-RunAB-sc/"
dir = '/Users/yzhu14/Research/forBrettPaper/'
# "160-600kmACC12Step3NE30NH011818sc.dat2p" using ($1/1e5):(abs(($7-$8)/$2)) with lines ls 2 ti "Hnn/Ve",\
nfit=1
outfile = open(dir +'pe.dat','w')
#for file in range(0,30):
#    file1 = dir + 'vnue-'+'{:02}'.format(1+file)+'-nsnssc03-000Mo-022deg.dat'
filename=['160-600kmACC12Step3NE30NH011818sc.dat2p']
for file in filename:
    file1 = dir + file
    data = np.loadtxt(file1,skiprows=n)
    pe = data[:,6]
    pa = data[:,7]
    ve = data[:,1]
    x = data[:,0]
    ya = pa/ve
    ye = pe/ve
    #xx = np.zeros(math.floor(len(x)/nfit))
    #yy = np.zeros(math.floor(len(x)/nfit))
    #for j in range(math.floor(len(x)/nfit)):
    #    xx[j]=x[nfit*j]
    #    yy[j]=y[nfit*j]
    xx = x/1e5
    yy = ye*1e19
    yIn = 1/ye
    for x1, y1 in zip(xx, yIn):
        plt.plot(x1, 1/y1, 'ro')
    zz = np.polyfit(xx, yIn, 2)
    fe = np.poly1d(zz)
    outfile.write(str(file)+"\t"+str(float(zz[0]))+"\t"+str(float(zz[1]))+"\n")#+str(float(zz[2]))+"\n")
    print(zz)
    for x1 in np.linspace(100, 1000, 1e3):#
        plt.plot(x1, 1/fe(x1), 'b+')
    yy = ya
    yIn = 1/ya
    for x1, y1 in zip(xx, yIn):
        plt.plot(x1, 1/y1, 'yo')
    zz = np.polyfit(xx, yIn, 2)
    fa = np.poly1d(zz)
    outfile.write(str(file)+"\t"+str(float(zz[0]))+"\t"+str(float(zz[1]))+"\n")#+str(float(zz[2]))+"\n")
    print(zz)
    for x1 in np.linspace(100, 1000, 1e3):#
        plt.plot(x1, 1/fa(x1), 'g+')
plt.axis([xx[0], xx[-1],0.01/yIn[-1],1/yIn[0]])
plt.show()
outfile.close()
