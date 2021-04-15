# read in file and get the xyz vector
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from hnunu.albinoBasics import *


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
            # mirror tracers that below the merger core
            outfilexyz.write(("\t%3.10f\n" % (float(xg[i]))))
            outfilexyz.write(("\t%3.10f\n" % (float(yg[i]))))
            outfilexyz.write(("\t%3.10f\n" % (float(abs(zg[i])))))

    print(str(id), str(zg[-1]))
    outfilexyz.close()
    # get the new number of lines of tracer file upto zmax
    nline_tracer = int(sum(1 for line in open(diroutput + fileName + "xyz.txt")) / 3)

    return nline_tracer
