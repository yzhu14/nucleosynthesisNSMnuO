#!/usr/bin/env python

'''
   Created: 03/16/2018: Analysis script for Yonglin.
   Updated: ??/??/2018: ?
   Author : Matthew Mumpower {matthew@mumpower.net}
   Description: Object containers for PRISM output.
   Intended: Python 2.7.11+


to do : time: t/tsec/tday
'''
#['ad', 'ncap', 'gn', 'nif', 'bdf', 'spf', 'bmd']
# ========================
import numpy as np


class AttrDict(dict):
  '''
     Treats dict key like an attribute
  '''

  def __getattr__(self, k):
    return self[k]

  def __setattr__(self, k, v):
    self[k] = v

  def _asdict(self):
    return self


class AbTime(object):
  '''
     Holds abundance info at each timestep.
  '''
  def __init__(self, filename):
    '''
       Initialize the object from the path and file name.
    '''
    self.data = dict()
    self.read(filename)

  def read(self, filename):
    '''
       Read the data from file
    '''
    f = file(filename)
    f.readline() # read header line
    # Loop over time data
    for line in f:
      li = line.split()
      if "timestep" in line:
        m  = int(li[1])
        self.data[m] = AttrDict()
        self.data[m].Y = AttrDict()
        self.data[m].Ya = {}
        self.data[m].t = float(li[3])
        self.data[m].T9 = float(li[5])
        self.data[m].rho = float(li[7])
      elif len(li) == 3:
        Z = int(li[0])
        A = int(li[1])
        N = A - Z
        Y = float(li[2])
        self.data[m].Y[Z,N] = Y
        if self.data[m].Ya.get(A,0) != 0:
          self.data[m].Ya[A] = self.data[m].Ya[A] + Y
        else:
          self.data[m].Ya[A] = Y

  def find_m_closest_to_t(self, t):
    '''
       Find a timestep for which the time in the simulation is closest to the given time (t) in seconds.
    '''
    tlst = []
    for m in sorted(self.data.keys()):
      tlst.append(abs(t-self.data[m].t))
    return tlst.index(min(tlst))

  def find_m_closest_to_T9(self, T9):
    '''
       Find a timestep for which the temperature [GK] in the simulation is closest to the given T9.
    '''
    t9lst = []
    for m in sorted(self.data.keys()):
      t9lst.append(abs(T9-self.data[m].T9))
    return t9lst.index(min(t9lst))


class Flows(object):
  '''
     Holds flow info at each timestep.
  '''
  def __init__(self, filename):
    '''
       Initialize the object from the path and file name.
    '''
    self.data = dict()
    self.read(filename)

  def read(self, filename):
    '''
       Read the data from file
    '''
    f = file(filename)
    self.names = f.readline().split()
    # Remove Z, A headers
    self.names.pop(0)
    self.names.pop(0)

    # Loop over flow data
    for line in f:
      li = line.split()
      if "timestep" in line:
        m  = int(li[1])
        self.data[m] = AttrDict()
        self.data[m].f = AttrDict()
        self.data[m].t = float(li[3])
        self.data[m].T9 = float(li[5])
        self.data[m].rho = float(li[7])
        self.data[m].top = np.array([0,0,0.]*5).reshape(5, 3)
        topbmd = np.array([0,0,0.]*5).reshape(5, 3)
        topspf = np.array([0,0,0.]*5).reshape(5, 3)
      elif len(li)-2 == len(self.names):
        Z = int(li[0])
        A = int(li[1])
        N = A - Z
        self.data[m].f[Z,N] = []
        for i in range(0, len(self.names)):
          self.data[m].f[Z,N].append(float(li[i+2]))
        #################
        for j in range(5):
            if j == 5:
                if self.data[m].f[Z,N][j]>topbmd[j][2]:
                    topbmd[j][0]=int(Z)
                    topbmd[j][1]=int(N)
                    topbmd[j][2]=self.data[m].f[Z,N][j]
      self.data[m].top = topbmd

  def find_m_closest_to_t(self, t):
    '''
       Find a timestep for which the time in the simulation is closest to the given time (t) in seconds.
    '''
    tlst = []
    for m in sorted(self.data.keys()):
      tlst.append(abs(t-self.data[m].t))
    return tlst.index(min(tlst))

  def find_m_closest_to_T9(self, T9):
    '''
       Find a timestep for which the temperature [GK] in the simulation is closest to the given T9.
    '''
    t9lst = []
    for m in sorted(self.data.keys()):
      t9lst.append(abs(T9-self.data[m].T9))
    return t9lst.index(min(t9lst))


class Flowsin(object):
  '''
     Holds flow info at each timestep.
  '''
  def __init__(self, filename):
    '''
       Initialize the object from the path and file name.
    '''
    self.data = dict()
    self.read(filename)

  def read(self, filename):
    '''
       Read the data from file
    '''
    f = file(filename)
    self.names = f.readline().split()
    # Remove Z, A headers
    self.names.pop(0)
    self.names.pop(0)
    # Loop over flow data
    flowchannels=["ad","ncap","gn","nif","bdf","spf","bmd"]
    gettimestep = False
    for line in f:
      li = line.split()
      if "0.0" in line:
        m  = int(li[0])
        self.data[m] = AttrDict()
        self.data[m].tsec = float(li[1])
        self.data[m].tday = float(li[1])/86400
        gettimestep = True
#        oZ = int(li[2])
#        oN = int(li[3])
#        oA = oN + oZ
        self.data[m].fout= AttrDict()
        for i in range(len(flowchannels)):
            self.data[m].fout[str(flowchannels[i])] = float(li[i+4])
        self.data[m].totalout = float(li[-2])
        self.data[m].totalin =  float(li[-1])
        self.data[m].p = AttrDict()
      elif 'parent' in line: 
        pZ = int(li[1])
        pN = int(li[2])
        pA = pN + pZ
        rn = str(li[3])
        self.data[m].p[pZ,pN] = AttrDict() 
        self.data[m].p[pZ,pN].fin = float(li[4])
        self.data[m].p[pZ,pN].ch = rn 

  def find_m_closest_to_t(self, t):
    '''
       Find a timestep for which the time in the simulation is closest to the given time (t) in seconds.
    '''
    tlst = []
    for m in sorted(self.data.keys()):
      tlst.append(abs(t-self.data[m].tsec))
    return tlst.index(min(tlst))


class Topnq(object):
  '''

  '''
  def __init__(self, filename):
    '''
       Initialize the object from the path and file name.
    '''
    self.data = dict()
    self.read(filename)

  def read(self, filename):
    '''
       Read the top heating data from file
    '''
    f = file(filename)
    self.names = f.readline().split()
    m=0
    for line in f:
      li = line.split()#
      if len(li)>5:
        m+=1
        top = 1
        self.data[m] = AttrDict()
        self.data[m].t = AttrDict()
        self.data[m].h = AttrDict()
        self.data[m].q = AttrDict()
        self.data[m].tsec = float(li[1])
        self.data[m].tday = float(li[3])
      elif len(li)== 4 or len(li)== 5:
        Z = int(li[0])
        A = int(li[1])
        N = A - Z
        self.data[m].t[top]= AttrDict()
        self.data[m].t[top].h=float(li[2]) 
        self.data[m].t[top].q=float(li[3]) 
        self.data[m].t[top].z= Z
        self.data[m].t[top].a= A
        self.data[m].t[top].n= N 
        top+=1
        H = float(li[2])
        self.data[m].h[Z,N] = H
        self.data[m].q[Z,N]=float(li[3]) 


  def find_m_closest_to_t(self, t):
    '''
       Find a timestep for which the time in the simulation is closest to the given time (t) in seconds.
    '''
    tlst = []
    for m in sorted(self.data.keys()):
      tlst.append(abs(t-self.data[m].tsec))
    return tlst.index(min(tlst))


class Topn(object):
  '''

  '''
  def __init__(self, filename):
    '''
       Initialize the object from the path and file name.
    '''
    self.data = dict()
    self.read(filename)

  def read(self, filename):
    '''
       Read the top heating data from file
    '''
    f = file(filename)
    self.names = f.readline().split()
    m=0
    for line in f:
      li = line.split()#
      if len(li)>3:
        m+=1
        top = 1
        self.data[m] = AttrDict()
        self.data[m].t = AttrDict()
        self.data[m].h = AttrDict()

        self.data[m].tsec = float(li[0])
        self.data[m].tday = float(li[2])
      elif len(li)==2:
        top = 1
        self.data[m] = AttrDict()
        self.data[m].t = AttrDict()
        self.data[m].h = AttrDict()

        self.data[m].tsec = float(li[0])
        self.data[m].tday = float(li[1])
      elif len(li)== 3:
        Z = int(li[0])
        A = int(li[1])
        N = A - Z
        self.data[m].t[top]= AttrDict()
        self.data[m].t[top].h=float(li[2]) 
        self.data[m].t[top].z= Z
        self.data[m].t[top].a= A
        self.data[m].t[top].n= N 
        top+=1
        H = float(li[2])
        self.data[m].h[Z,N] = H

  def find_m_closest_to_t(self, t):
    '''
       Find a timestep for which the time in the simulation is closest to the given time (t) in seconds.
    '''
    tlst = []
    for m in sorted(self.data.keys()):
      tlst.append(abs(t-self.data[m].tsec))
    return tlst.index(min(tlst))

class heat(object):
  '''

  '''
  def __init__(self, filename):
    '''
       Initialize the object from the path and file name.
    '''
    self.data = dict()
    self.read(filename)

  def read(self, filename):
    '''
       Read the data from file
    '''
    f = file(filename)
    self.names = f.readline().split()
    # Remove Z, A headers
    self.names.pop(0)
    self.names.pop(0)
    # Loop over flow data
    for line in f:
      li = line.split()
      if "Z" in line:
        continue
      elif "timestep" in line:
        m  = int(li[1])
        self.data[m] = AttrDict()
        self.data[m].ch = AttrDict()
        self.data[m].t = float(li[3])
        self.data[m].T = float(li[5])
        self.data[m].rho = float(li[7])
        self.data[m].h = AttrDict()
      elif len(li)-2 == len(self.names):
        Z = int(li[0])
        A = int(li[1])
        N = A-Z
        self.data[m].h[Z,N] =AttrDict() 
        self.data[m].h[Z,N].ch =AttrDict()
        for i in range(0, len(self.names)):
          self.data[m].h[Z,N].ch[i] = float(li[i+2])

  def find_m_closest_to_t(self, t):
    '''
       Find a timestep for which the time in the simulation is closest to the given time (t) in seconds.
    '''
    tlst = []
    for m in sorted(self.data.keys()):
        tlst.append(abs(t-self.data[m].t))

    return tlst.index(min(tlst))

  def find_max_heat(self,z,n,nchannel,mintsec,maxtsec):
    maxheat = 0
    for m in sorted(self.data.keys()):
        if self.data[m].t>mintsec and self.data[m].t<maxtsec:
            if self.data[m].h[z,n].ch[nchannel] > maxheat:
                maxheat = self.data[m].h[z,n].ch[nchannel]
    return maxheat


class totalHeat(object):
  '''

  '''
  def __init__(self, filename):
    '''
       Initialize the object from the path and file name.
    '''
    self.data = dict()
    self.read(filename)

  def read(self, filename):
    '''
       Read the data from file
    '''
    f = file(filename)
    self.names = f.readline().split()
    # Remove Z, A headers
    self.names.pop(0)
    self.names.pop(0)
    # Loop over flow data
    m=0
    for line in f:
      li = line.split()#
      m+=1
      if len(li) == 8: # time and 7 channels
        #m = int(li[0])-1
        self.data[m] = AttrDict()
        self.data[m].H = AttrDict()
        self.data[m].tsec = float(li[0])
        self.data[m].tday = float(li[0])/86400
        self.data[m].H = []
        for i in range(7):
          self.data[m].H.append(float(li[i+1]))
      else:
        continue
  def find_m_closest_to_t(self, t):
    '''
       Find a timestep for which the time in the simulation is closest to the given time (t) in seconds.
    '''
    tlst = []
    for m in sorted(self.data.keys()): tlst.append(abs(t-self.data[m].tsec))
    return tlst.index(min(tlst))




class tau(object):
  '''

  '''
  def __init__(self, filename):
    '''
       Initialize the object from the path and file name.
    '''
    self.data = dict()
    self.read(filename)

  def read(self, filename):
    '''
       Read the data from file
    '''
    f = file(filename)
    self.names = f.readline().split()
    # Remove Z, A headers
    self.names.pop(0)
    self.names.pop(0)
    # Loop over flow data
    m=0
    for line in f:
      li = line.split()#
      m+=1
      if len(li) == 8: # time and 7 channels
        #m = int(li[0])-1
        self.data[m] = AttrDict()
        self.data[m].tau = AttrDict()
        self.data[m].tsec = float(li[0])
        self.data[m].tday = float(li[0])/86400
        self.data[m].tau = []
        for i in range(7):
          self.data[m].tau.append(float(li[i+1]))
      else:
        continue
  def find_m_closest_to_t(self, t):
    '''
       Find a timestep for which the time in the simulation is closest to the given time (t) in seconds.
    '''
    tlst = []
    for m in sorted(self.data.keys()):
        tlst.append(abs(t-self.data[m].tsec))
    return tlst.index(min(tlst))

