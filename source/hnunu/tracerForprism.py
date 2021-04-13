__author__ = "yzhu14"
from scipy.interpolate import interp1d
import numpy as np
import math
import itertools as IT
from decimal import Decimal
from hnunu.unit_const import *
from hnunu.albinoBasics import get_tracer_id

# tracerlist=['57221',"61778",'66394','78323','79049',"80224"]

# "/Users/yzhu14/Research/nucleosynthesisNSMnuO/testTracerofDirk/tracers/"#s/physics/yzhu14/nucleosyn/thermo/tracer/"
tmax = 10  # GK
initialS = 29
dt = 1e-3
ds = 1e-3
tempTol = 1e-10  # erg
commentline = 2  # of commentline in tracers
#############################################################
######################################### getLastLine #################################################################
#   Yonglin Zhu 2017
#   input:ntestpoints - # of line from input files, id, dirOut
#   output: .run file for Thermodynamic.cpp
#
####################################################################################################################
def prep_thermodynamic_input(line,ntestpoints, id, dirOut):
    """[summary]
        Take the last line of the tracer and initial parameter
        to gengerate Thermodynamic.cpp input file (id+new.run)

    Args:
        ntestpoints ([type]): [description]
        id ([type]): [description]
        dirOut ([type]): [description]
    """
    #   time        x coord        y coord        z coord         radius           dens           temp             ye            |v|             vx             vy             vz
    #    sec             km             km             km             km         g/cm^3             GK              -           cm/s           cm/s           cm/s           cm/s
    thermo = line[ntestpoints - 1]
    parameter = thermo.split()
    tfin = float(parameter[0])
    t9fin = float(parameter[6])
    rhofin = float(parameter[5])
    rfin = float(
        math.sqrt(
            float(parameter[1]) ** 2
            + float(parameter[2]) ** 2
            + float(parameter[3]) ** 2
        )
    )
    yefin = float(parameter[7])
    # print(id,ntestpoints,tfin,t9fin,rhofin,rfin,
    print(id, "ye", yefin)
    outforThermorun = open(dirOut + str(id) + "new.run", "w")
    outforThermorun.write(
        str(str(id) + "ex.txt")
        + "\n"
        + str(tfin)
        + "\n"
        + str(t9fin)
        + "\n"
        + str(rhofin)
        + "\n"
        + str(rfin)
        + "\n"
        + str(yefin)
        + "\n"
        + str(initialS)
        + "\n"
        + str(tmax)
        + "\n"
        + str(dt)
        + "\n"
        + str(ds)
        + "\n"
        + str(tempTol)
        + "\n"
    )
    outforThermorun.close()


######################################### etend_tracer #################################################################
#   Yonglin Zhu 2017
#   combine
#   input:ntestpoints - # of line from input files, ntestpoints - # of line from input files,
#
####################################################################################################################
def etend_tracer(R,time,line,ntestpoints,dirIN,dirOut):
    """[summary]
    Part two use the output of Thermodynamic.cpp and combine it with original tracer, fit a radius to construct input for prism

    """
    outforFulltracer = open(dirOut + str(id) + "newprism.txt", "w")
    for n in range(commentline, ntestpoints):
        thermo = line[n]
        parameter = thermo.split()
        time[n - commentline] = float(parameter[0])
        t9 = float(parameter[6])
        rho = float(parameter[5])
        ye = float(parameter[7])
        R[n - commentline] = float(
            math.sqrt(
                float(parameter[1]) ** 2
                + float(parameter[2]) ** 2
                + float(parameter[3]) ** 2
            )
        )
        if n == commentline:
            print(
                id,
                str(time[n - commentline])
                + " "
                + str(t9)
                + " "
                + str(rho)
                + " "
                + str(float(R[n - commentline]) * km),
            )
        outforFulltracer.write(
            (
                "%3.6e %3.6e %3.6e %3.6e\n"
                % (
                    Decimal(Decimal(str(time[n - commentline]))),
                    float(str(t9)),
                    Decimal(str(rho)),
                    Decimal(str(float(R[n - commentline]) * km)),
                )
            )
        )
    #############fit the radius#######
    r = np.polyfit(time[-5:], R[-5:], 1)  # use last few points to fit
    f = np.poly1d(r)
    ####################################
    intracerEx = open(dirIN + str(id) + "ex.txt", "r")
    exline = intracerEx.readlines()
    exstep = sum(1 for line in open(dirIN + str(id) + "ex.txt"))
    for i in range(0, exstep):
        exthermo = exline[i]
        expara = exthermo.split()
        exR = f(float(expara[0])) * km
        outforFulltracer.write(
            (
                "%3.6e %3.6e %3.6e %3.6e\n"
                % (
                    Decimal(str(expara[0])),
                    float(str(expara[1])),
                    Decimal(str(expara[2])),
                    Decimal(str(exR)),
                )
            )
        )
    outforFulltracer.close()


