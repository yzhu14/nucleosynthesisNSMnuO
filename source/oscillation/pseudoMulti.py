import matplotlib as mpl
import matplotlib.pyplot as plt
from albinoBasics import *
from neutrinoCap import *
from mpl_toolkits.mplot3d import Axes3D
from math import *
import math


def picksample(dir,fileID,nsample):
    for id in fileID:
        outFileSample = open(dir + '/tracer/'+ str(nsample) +'sample'+str(id)+'.txt',
                             'w')
        fileName="/tracer/tracer_"+str(id)+".txt" #","06174.txt
        nfile = sum(1 for line in open(dir+fileName))
        infile = open(dir+fileName,"r")
        lineS = infile.readlines()
        x = np.zeros((nfile,1))
        y = np.zeros((nfile,1))
        z = np.zeros((nfile,1))
        xn = np.zeros((nsample,1))
        yn = np.zeros((nsample,1))
        zn = np.zeros((nsample,1))
        stepAdjust=int(math.ceil(nfile/nsample))
        #print("stepAdjust",stepAdjust,nsample)
        for i in range(nCommentLines,nfile):
            corrS=lineS[i]
            pointPS=corrS.split()
            x[i]=pointPS[1]
            y[i]=pointPS[2]
            z[i]=pointPS[3]
            #print(x[i],y[i],z[i])
            if (i-nCommentLines-1) % stepAdjust == 0:
                #print(i,i/stepAdjust,int(i/stepAdjust))
                n=int(i/stepAdjust)
                xn[n]=float(pointPS[1])
                yn[n]=float(pointPS[2])
                zn[n]=float(pointPS[3])
                #print(n,i,float(pointPS[1]),float(pointPS[2]),float(pointPS[3]))
                outFileSample.write(("\t%3.10f\n" % (float(xn[n]))))
                outFileSample.write(("\t%3.10f\n" % (float(yn[n]))))
                outFileSample.write(("\t%3.10f\n" % (float(zn[n]))))
        outFileSample.close()
        plotxyz(id,xn,yn,zn,x,y,z)

def plotxyz(id,xn,yn,zn,x,y,z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(xn[1:], yn[1:], zn[1:], label='Tracer',s=100, c='g')
    ax.scatter(x[1:], y[1:], z[1:], label='Sample')
    ax.legend()
    ax.set_xlim(-500,500)
    ax.set_ylim(-500,500)
    ax.set_zlim(-0,1200)
    plt.savefig(str(id)+'sampleXYZ.pdf')


def sampletra(samplex,dir,fileID):
    step=150
    stepsize=0.5
    for id in fileID:
        fileName=dir + '/tracer/'+ str(id)+'ntest0.txt' #","06174.txt
        nfile = sum(1 for line in open(fileName))
        infile = open(fileName,"r")
        lineS = infile.readlines()
        x = np.zeros((nfile,1))
        y = np.zeros((nfile,1))
        z = np.zeros((nfile,1))
        vector = np.zeros((nfile,3))
        print(vector[1])
        for i in range(0,nfile):
            corrS=lineS[i]
            pointPS=corrS.split()
            x[i]=float(pointPS[6])
            y[i]=float(pointPS[7])
            z[i]=float(pointPS[8])
        x=np.unique(x)
        y=np.unique(y)
        z=np.unique(z)
        if len(x) == len(y) and len(y) == len(z):
            for i in range(len(x)):
                outsampletra=open(dir + '/tracer/'+ str(id)+"sample"+str(i)+'.txt','w')
                #print(vector)
                for j in range(2*step):
                    r=stepsize*(j-step)
                    xstep=float(samplex[0])+r*float(x[i])
                    ystep=float(samplex[1])+r*float(y[i])
                    zstep=float(samplex[2])+r*float(z[i])
                    outsampletra.write("   "+str(xstep)+"\n"+"   "+str(ystep)+"\n"+"   "+str(zstep)+"\n")
                    if j ==(1):
                        print(i,math.sqrt(xstep**2+ystep**2+zstep**2))
                    if j ==(2*step-1):
                        print(i,math.sqrt(xstep**2+ystep**2+zstep**2))
                outsampletra.close()



def plotangle(dir,fileID):
    for id in fileID:
        fileName=dir + '/tracer/'+ str(id)+'ntest0.txt' #","06174.txt
        nfile = sum(1 for line in open(fileName))
        infile = open(fileName,"r")
        lineS = infile.readlines()
        x = np.zeros((nfile,1))
        y = np.zeros((nfile,1))
        z = np.zeros((nfile,1))
        for i in range(0,nfile):
            corrS=lineS[i]
            pointPS=corrS.split()
            x[i]=pointPS[6]
            y[i]=pointPS[7]
            z[i]=pointPS[8]
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.scatter(x[1:], y[1:], z[1:], label='UV')
        ax.legend()
        ax.set_xlim(-2,2)
        ax.set_ylim(-1.5,1.5)
        ax.set_zlim(-0,2)
        plt.savefig(str(id)+'angle.pdf')

######################################### getTra() #################################################################
#   Yonglin Zhu 2017
#   input: nf - file#, NE - # of energy bins, ntestpoints - # of line from input files, nNew - new # of energy bins
#   output: potential/density files with NE energy bins, Ye/rho files
#   check: R_original(core or not), potential/density units, potential/density normalization, interplote algorithm
#   to do: checkmorebins()
#
####################################################################################################################
def getTra():
    inFileSample = open(dir + '/tracer/'+ str(nsample) +'sample'+str(id)+'.txt','r')
    nfile = sum(1 for line in open(dir + '/tracer/'+ str(nsample) +'sample'+str(id)+'.txt'))
    lineS = inFileSample.readlines()

def getTra(a,b,c):
    #Target Point a b c,
    nTarget=len(a)
    #read in neutrino surface
    dirNeu="/Users/yzhu14/Dropbox/paperOne/results/NeuSphere/nu_surf_elan_e08_taueff.txt"
    nNeufile = sum(1 for line in open(dirNeu))
    infileNeu= open(dirNeu,"r")
    lineNeu = infileNeu.readlines()
    rp = np.zeros((nNeufile,1))
    xp = np.zeros((nNeufile,1))
    yp = np.zeros((nNeufile,1))
    zp = np.zeros((nNeufile,1))
    for i in range(1,nNeufile):
        corrNeu=lineNeu[i]
        pointsNeu=corrNeu.split()
        rp[i]=float(pointsNeu[0])
        zp[i]=float(pointsNeu[1])
    #        print("rp=zp",rp[i],zp[i])
    ###starting point zp/rp  nZ
    #######clean this z colunm before count line, becasue z is 0 after certain distance.
    #trim the leading or trailing zeros form the 1-d array
    zp=np.trim_zeros(zp)
    nZ=len(zp)
    rp=rp[1:nZ+1]
    #print("length",len(rp),len(zp))
    if len(rp)==len(zp):
        print("True")
        #print("print",rp,zp)
    ##########
    #Starting point nA,nR,
    #target point abc,len(a)
    print(zp)
    nR=4
    nA=4
    print("rp",len(rp))
    rRange=rp[-1]-rp[0]
    rStep=nZ/nR
    stepSize=5 #km
    step=10000
    tradr=np.zeros((nTarget,nA,nR,3))
    #print("?",tradr[1][2][1])
    outFileFT=open(dir+fileID+"/AP.txt",'w')
    print(nTarget,nA,nR)
    for j in range(1,nTarget):
        for A in range(0,nA):#angel
            for R in range(0,nR):#radius
                print(zp[int(1+R*rStep/2)],int(R*rStep/2),"??")
                print("int J",int(j))
                traX=float(a[j]-rp[int(R*rStep/2)]*math.sin(A*2*pi/nA))
                traY=float(b[j]-rp[int(R*rStep/2)]*math.cos(A*2*pi/nA))
                traZ=float(c[j]-zp[int(R*rStep/2)])
                traR=math.sqrt(traX*traX+traY*traY+traZ*traZ)
                #print("z=",c[j],traZ,zp[int(R*rStep/2)])

                if traR!=0:
                    tradr[j][A][R]=np.array([traX/traR,traY/traR,traZ/traR])
                    print(tradr[j][A][R],zp[int(R*rStep/2)])
                    outFile=open(dir+fileID+"/"+str(int(j))+"_"+str(A)+"_"+str(R)+".txt",'w')
                    outFileAP=open(dir+fileID+"/AP"+str(int(j))+"_"+str(A)+"_"+str(R)+".txt",'w')
                    outFileFT.write(str(float(a[j]))+"\n"+
                                    str(float(rp[int(R*rStep/2)]*math.sin(A*2*pi/nA)))+"\n"+
                                    str(float(b[j]))+"\n"+
                                    str(float(rp[int(R*rStep/2)]*math.cos(A*2*pi/nA)))+"\n"+
                                    str(float(c[j]))+"\n"+
                                    str(float(zp[int(R*rStep/2)]))+"\n")
                    for k in range(0,step):
                        #print(k)
                        if k*stepSize*tradr[j][A][R][2]<0:
                            print("error",k*stepSize*tradr[j][A][R][2])
                        if float(k*stepSize*tradr[j][A][R][2]+zp[int(R*rStep/2)])!=0:
                            outFile.write(str(float(k*stepSize*tradr[j][A][R][0]+rp[int(R*rStep/2)]*math.sin(A*2*pi/nA)))+
                                          "\t"+str(float(k*stepSize*tradr[j][A][R][1]+rp[int(R*rStep/2)]*math.cos(A*2*pi/nA)))+
                                          "\t"+str(float(k*stepSize*tradr[j][A][R][2]+zp[int(R*rStep/2)]))+"\n")
                            outFileAP.write("   "+str(float(k*stepSize*tradr[j][A][R][0]+rp[int(R*rStep/2)]*math.sin(A*2*pi/nA)))+
                                            "\n   "+str(float(k*stepSize*tradr[j][A][R][1]+rp[int(R*rStep/2)]*math.cos(A*2*pi/nA)))+
                                            "\n   "+str(float(k*stepSize*tradr[j][A][R][2]+zp[int(R*rStep/2)]))+"\n")

                            #print("??")

    outFile.close()
    outFileAP.close()
    outFileFT.close()
    #generate trajectory



    ########################################################################
def get_tracer_id(tracers):
    if tracers=="xy70":
        fileID =["14746", "14782","15364","15377","15378","15379","15388","15400",
                 "15407","15412","15414","17090","17096","17704","17708","17721",
                 "17779","17781","17783","17786","17789","17792","22290"]
        dir="/Users/yzhu14/Research/hnunu/tests/data/neuCap/xy70/forNueosc"
    elif tracers=="dirk":
        fileID=["61778"]
        #fileID=['57221',"61778",'66394','78323','79049',"80224"]
        dir="/Users/yzhu14/Research/hnunu/tests/data/neuCap/dirksample/forNueosc"
    else:
        print("Wrong Dir!")
    return fileID,dir
#######################################################################

tracers="dirk"
fileID,dir=get_tracer_id(tracers)
newElist=[16]
start = "0"
dir="/Users/yzhu14/Research/hnunu/tests/data/tracers/dirksample"
######### Main() ########
outfilelist=open(dir+"tracerList.txt","w")
deleteTracer= np.zeros((10000,1))
zmax=500
nsample=10
#for id in fileID:
#fileID=["61778"]
fileID=["17090"]
#picksample(dir,fileID,nsample)
#plotangle(dir,fileID)
#samplex=[-42.072,-90.95,124.90]
samplex=[-67.2711,-13.720,78.231]
sampletra(samplex,dir,fileID)


