import os
import sys
import traceback
import numpy as np
import pandas as pd
from prism_output import *
from network import *
dirt = '/Users/yzhu14/Research/nucleosynthesis/forYonglin/'
stable_path = dirt + 'input/data/'
channellist = ['ad', 'ncap', 'gn', 'nif', 'bdf', 'spf', 'bmd']
daytosec = 86400
mintime = 1e-1*daytosec
maxtime = 1e3*daytosec
betafiles = ['gtffsdn']
datelist = ['101819']
fission = np.array(["asyKarpov", "kodamaKarpov", "asyXuren", "kodamaXuren"])
#fissionlabel = (["Kodama, Xuren", "Kodama, Karpov", "Symmetric, Xuren", "Symmetric, Karpov"])
fissionlabel = (["Symmetric, Karpov","Kodama, Karpov","Symmetric, Xuren","Kodama, Xuren"])

def get_notelist():
    minye = 0
    maxye = 30
    yelist = [2, 12, 14, 16, 18, 21, 24, 28]
    notelist = []
    for i in yelist:
        notelist.append("ye"+'{:02}'.format(i))
    return notelist
notelist = get_notelist()
masslist = ["frdm2012", "hfb27", "hfb22", "ws3v6",
                "dz33", "dft_sly4", "dft_unedf1", "etfsi"]

def get_barrier(mass):
    if mass in ['hfb22','hfb27']:
       barrier = 'HFB'
    elif mass in ['etfsi']:
       barrier = 'ETFSI'
    else:
       barrier = 'FRLDM'
    return barrier

def getchannel(rn):
  if rn == 'ad':
      ichannel = 0
  elif rn == 'ncap':
      ichannel = 1
  elif rn == 'gn':
      ichannel = 2
  elif rn == 'nif':
      ichannel = 3
  elif rn == 'bdf':
      ichannel = 4
  elif rn == 'spf':
      ichannel = 5 #
  else:
      ichannel = 6 #beta decay
  return ichannel

cnlist = []
for i, rn in enumerate(channellist):
    cnlist.append(getchannel(rn))


def masslabel(masslist, massmodel):
    ml = masslist.index(massmodel)
    masslabels = ["FRDM2012", "HFB27", "HFB22",
                "WS3", "DZ33", "SLY4", "UNEDF1", "ETFSI"]
    return masslabels[ml]
##########################################################################################################################################
def yelabel(notelist, note):
    nl = notelist.index(note)
    notelabels = [r"$Y_e = 0.02$", r"$Y_e = 0.12$", r"$Y_e = 0.14$", r"$Y_e = 0.16$",
                  r"$Y_e = 0.18$", r"$Y_e = 0.21$", r"$Y_e = 0.24$", r"$Y_e = 0.28$"]
    return notelabels[nl]
####################################################################################################################


#################################################################################################
#################################################################################################
