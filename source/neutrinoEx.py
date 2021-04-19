# energy Bin interpolation from 8 bins from Albino's simulation(nu_potential_v2.0)
# For neutrino capture rate calculation
#March 2017
from scipy.interpolate import interp1d
import numpy as np
import math
import sys
from scipy.optimize import curve_fit
from hnunu.albinoBasics import *
from hnunu.neutrinoCap import *

#def get_tracer_id(tracers):
#    if tracers=="xy70":
#        fileID =["14746", "14782","15364","15377","15378","15379","15388","15400",
#                 "15407","15412","15414","17090","17096","17704","17708","17721",
#                 "17779","17781","17783","17786","17789","17792","22290"]
#        dir="/Users/yzhu14/Research/hnunu/tests/data/neuCap/xy70/forNueosc"
#    elif tracers=="dirk":
#        fileID=['57221',"61778",'66394','78323','79049',"80224"]
#        dir="/Users/yzhu14/Research/hnunu/tests/data/neuCap/dirksample/forNueosc"
#    else:
#        print("Wrong Dir!")
#    return fileID,dir
#

note="log1"
# dir="/Users/yzhu14/Research/hnunu/tests/data/neuCap/dirksample/forNueosc"
# #dir="/Users/yzhu14/Research/hnunu/tests/data/neuCap/xy70/forNueosc"
# #dir="/Users/yzhu14/Research/nucleosynthesisNSMnuO/complete_tracers_dirk_sample/density/"
# fileID =["14746","14782","15364","15377","15378","15379","15388","15400",
#          "15407","15412","15414","17090","17096","17704","17708","17721",
#          "17779","17781","17783","17786","17789","17792","22290"]
# fileID=['57221',"61778",'66394','78323','79049',"80224"]

whichFileP = 'v_potential'
whichFileN = 'v_density'
tracers="dirk"
###################################
fileID,dirt=get_tracer_id(tracers)
newElist=[16]
for nNew in newElist:
    for k in fileID:
        Nucapfile=dirt + '/' + str(nNew) + 'NE/'  + str(k) + note
        exNuCap(Nucapfile,1000e5,nNew)
    #outfileNC.write(str(pointsDensity[0])+"\t"+str(noOscillationRateMatter*Cnunu)+'\t'+str(currentRateMatter*Cnunu)+"\t"+str(noOscillationRateAntiMatter*Cnunu)
