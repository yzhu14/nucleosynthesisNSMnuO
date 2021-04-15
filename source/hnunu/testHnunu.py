# read in file and get the xyz vector
from shared.albinoBasics import get_tracer_id
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from shared.albinoBasics import *

######################################### getTra() #################################################################
#   Yonglin Zhu 2017
#   input: nf - file#, NE - # of energy bins, ntestpoints - # of line from input files, nNew - new # of energy bins
#   output: potential/density files with NE energy bins, Ye/rho files
#   check: R_original(core or not), potential/density units, potential/density normalization, interplote algorithm
#   to do: checkmorebins()
#
####################################################################################################################
def getTra(a, b, c):
    """[summary]

    Args:
        a ([type]): [description]
        b ([type]): [description]
        c ([type]): [description]
    """
    # Target Point a b c,
    nTarget = len(a)
    # read in neutrino surface
    dirNeu = (
        "../../data/neuSphere/nu_surf_elan_e08_taueff.txt"
    )
    nNeufile = sum(1 for line in open(dirNeu))
    infileNeu = open(dirNeu, "r")
    lineNeu = infileNeu.readlines()
    rp = np.zeros((nNeufile, 1))
    xp = np.zeros((nNeufile, 1))
    yp = np.zeros((nNeufile, 1))
    zp = np.zeros((nNeufile, 1))
    for i in range(1, nNeufile):
        corrNeu = lineNeu[i]
        pointsNeu = corrNeu.split()
        rp[i] = float(pointsNeu[0])
        zp[i] = float(pointsNeu[1])
    #        print("rp=zp",rp[i],zp[i])
    ###starting point zp/rp  nZ
    #######clean this z colunm before count line, becasue z is 0 after certain distance.
    # trim the leading or trailing zeros form the 1-d array
    zp = np.trim_zeros(zp)
    nZ = len(zp)
    rp = rp[1 : nZ + 1]
    # print("length",len(rp),len(zp))
    if len(rp) == len(zp):
        print("True")
    # print("print",rp,zp)
    ##########
    # Starting point nA,nR,
    # target point abc,len(a)
    print(zp)
    nR = 4
    nA = 4
    print("rp", len(rp))
    rRange = rp[-1] - rp[0]
    rStep = nZ / nR
    stepSize = 5  # km
    step = 10000
    tradr = np.zeros((nTarget, nA, nR, 3))
    # print("?",tradr[1][2][1])
    outFileFT = open(dirtracer + fileID + "/AP.txt", "w")
    print(nTarget, nA, nR)
    for j in range(1, nTarget):
        for A in range(0, nA):  # angel
            for R in range(0, nR):  # radius
                print(zp[int(1 + R * rStep / 2)], int(R * rStep / 2), "??")
                print("int J", int(j))
                traX = float(a[j] - rp[int(R * rStep / 2)] * math.sin(A * 2 * pi / nA))
                traY = float(b[j] - rp[int(R * rStep / 2)] * math.cos(A * 2 * pi / nA))
                traZ = float(c[j] - zp[int(R * rStep / 2)])
                traR = math.sqrt(traX * traX + traY * traY + traZ * traZ)
                # print("z=",c[j],traZ,zp[int(R*rStep/2)])

                if traR != 0:
                    tradr[j][A][R] = np.array([traX / traR, traY / traR, traZ / traR])
                    print(tradr[j][A][R], zp[int(R * rStep / 2)])
                    outFile = open(
                        dirtracer
                        + fileID
                        + "/"
                        + str(int(j))
                        + "_"
                        + str(A)
                        + "_"
                        + str(R)
                        + ".txt",
                        "w",
                    )
                    outFileAP = open(
                        dirtracer
                        + fileID
                        + "/AP"
                        + str(int(j))
                        + "_"
                        + str(A)
                        + "_"
                        + str(R)
                        + ".txt",
                        "w",
                    )
                    outFileFT.write(
                        str(float(a[j]))
                        + "\n"
                        + str(float(rp[int(R * rStep / 2)] * math.sin(A * 2 * pi / nA)))
                        + "\n"
                        + str(float(b[j]))
                        + "\n"
                        + str(float(rp[int(R * rStep / 2)] * math.cos(A * 2 * pi / nA)))
                        + "\n"
                        + str(float(c[j]))
                        + "\n"
                        + str(float(zp[int(R * rStep / 2)]))
                        + "\n"
                    )
                    for k in range(0, step):
                        # print(k)
                        if k * stepSize * tradr[j][A][R][2] < 0:
                            print("error", k * stepSize * tradr[j][A][R][2])
                        if (
                            float(
                                k * stepSize * tradr[j][A][R][2]
                                + zp[int(R * rStep / 2)]
                            )
                            != 0
                        ):
                            outFile.write(
                                str(
                                    float(
                                        k * stepSize * tradr[j][A][R][0]
                                        + rp[int(R * rStep / 2)]
                                        * math.sin(A * 2 * pi / nA)
                                    )
                                )
                                + "\t"
                                + str(
                                    float(
                                        k * stepSize * tradr[j][A][R][1]
                                        + rp[int(R * rStep / 2)]
                                        * math.cos(A * 2 * pi / nA)
                                    )
                                )
                                + "\t"
                                + str(
                                    float(
                                        k * stepSize * tradr[j][A][R][2]
                                        + zp[int(R * rStep / 2)]
                                    )
                                )
                                + "\n"
                            )
                            outFileAP.write(
                                "   "
                                + str(
                                    float(
                                        k * stepSize * tradr[j][A][R][0]
                                        + rp[int(R * rStep / 2)]
                                        * math.sin(A * 2 * pi / nA)
                                    )
                                )
                                + "\n   "
                                + str(
                                    float(
                                        k * stepSize * tradr[j][A][R][1]
                                        + rp[int(R * rStep / 2)]
                                        * math.cos(A * 2 * pi / nA)
                                    )
                                )
                                + "\n   "
                                + str(
                                    float(
                                        k * stepSize * tradr[j][A][R][2]
                                        + zp[int(R * rStep / 2)]
                                    )
                                )
                                + "\n"
                            )

                # print("??")

    outFile.close()
    outFileAP.close()
    outFileFT.close()
    # generate trajectory


# def getxyz():#as input for albino's code
def pickSample(stepAjust):
    for id in fileID:
        fileName = "/tracer/", "" + str(id) + ".txt"  # ","06174.txt
        nfile = sum(1 for line in open(dirtracer + fileName))
        infile = open(dirtracer + fileName, "r")
        lineS = infile.readlines()
        x = np.zeros((nfile, 1))
        y = np.zeros((nfile, 1))
        z = np.zeros((nfile, 1))
        xn = np.zeros((1 + int(nfile / stepAdjust), 1))
        yn = np.zeros((1 + int(nfile / stepAdjust), 1))
        zn = np.zeros((1 + int(nfile / stepAdjust), 1))
        print(nfile / stepAdjust, int(nfile / stepAdjust))
        for i in range(commentline, nfile):
            corrS = lineS[i]
            pointPS = corrS.split()
            x[i] = pointPS[1]
            y[i] = pointPS[2]
            z[i] = pointPS[3]
            # print(x[i],y[i],z[i])
            if i % stepAdjust == 0:
                # print(i,i/stepAdjust,int(i/stepAdjust))
                n = int(nfile / stepAdjust)
                xn[int(i / stepAdjust)] = float(pointPS[1])
                yn[int(i / stepAdjust)] = float(pointPS[2])
                zn[int(i / stepAdjust)] = float(pointPS[3])
                print(float(pointPS[1]), float(pointPS[2]), float(pointPS[3]))
        fig = plt.figure()
        ax = fig.gca(projection="3d")
        ax.scatter(xn[1:], yn[1:], zn[1:], label="parametric scatter")
        ax.legend()
        ax.set_xlim(-200, 100)
        ax.set_ylim(-200, 200)
        ax.set_zlim(-200, 400)
        plt.show()


def generate_tracer_nupotential(id, zmax, dirtracer,diroutput):
    """
        Generate tracer that nu_potential code take as input files.

    Args:
        id ([type]): id of tracers from Albino's simulation.
        zmax ([type]): Max z-axis for generated trajectories.

    Returns:
        [type]: [description]
    """
    fileName = "/tracer_" + str(id)  # ","06174.txt
    nline_tracer = sum(1 for line in open(dirtracer + fileName + ".txt"))
    print("\nChecking Number of lines in file:\t",nline_tracer)
    infile = open(dirtracer + fileName + ".txt", "r")
    outfilexyz = open(diroutput + fileName + "xyz.txt", "w")  # as input for albino's code
    linexyz = infile.readlines()
    xg = np.zeros((nline_tracer, 1))
    yg = np.zeros((nline_tracer, 1))
    zg = np.zeros((nline_tracer, 1))
    # read all coordinate first and trim them later
    for i in range(2, nline_tracer):
        corr = linexyz[i]
        points = corr.split()
        xg[i] = float(points[1])
        yg[i] = float(points[2])
        zg[i] = float(points[3])
        if float(zg[i]) > zmax:
            print("trim Z above", zmax, id, i)
            break
        else:
            outfilexyz.write(
                "   "
                + str(float(xg[i]))
                + "\n   "
                + str(float(yg[i]))
                + "\n   "
                + str(float(abs(zg[i]))) # mirror tracers that below the merger core
                + "\n"
            )
    print(str(id), str(zg[-1]))
    outfilexyz.close()
    # get the new number of lines of tracer file upto zmax
    nline_tracer = int(sum(1 for line in open(diroutput + fileName + "xyz.txt")) / 3)

    return nline_tracer

if __name__ == "__main__":

    diroutput = "../inoutput/neuField/"
    commentline = 2
    stepAdjust = 50
    dirkID, dirtracer= get_tracer_id('dirk')
    ###########
    # pickSample(stepAdjust)
    print(dirtracer)
    outfilelist = open(dirtracer + "tracerList.txt", "w")
    deleteTracer = np.zeros((10000, 1))
    zmax = 500
    for id in dirkID:
        nline = generate_tracer_nupotential(id, zmax,dirtracer,diroutput)
        outfilelist.write(str(id) + "\t" + str(zmax) + "\t" + str(nline) + "\n")
    outfilelist.close()
