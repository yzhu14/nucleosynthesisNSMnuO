#!/usr/bin/env python

'''
   Created: 04/03/2018: Analysis script for Yonglin.
   Author : Yonglin Zhu
   Description: Analysis of prism network
   Intended: Python 2.7.11+
'''

# ========================
#!/usr/bin/env python

'''
   Created: 04/03/2018: Analysis script for Yonglin.
   Author : Yonglin Zhu
   Description: Analysis of prism network
   Intended: Python 2.7.11+
'''

# ========================


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
'''

'''
class network(object):
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
    # Loop over  data
    for line in f:
      li = line.split()
      z = int(li[0])
      self.data[z] = AttrDict()
      self.data[z].z = z
      self.data[z].nmin  = int(li[1])
      self.data[z].nmax  = int(li[2])

class ab(object):
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

    # Loop over flow data
    for line in f:
      li = line.split()
      z = int(li[0])
      a = int(li[1])
      self.data[a] = AttrDict()
      self.data[a].y = AttrDict()
      self.data[a].y[z,a] = float(li[2])
class mfZ(object):
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

    # Loop over flow data
    sum =0
    for line in f:
      li = line.split()
      z = int(li[0])
      self.data[z] = AttrDict()
      self.data[z].z = z
      self.data[z].mfz  = float(li[1])
      sum +=float(li[1])
    if sum -1.0>0.02:
      print "error in mass fraction sum"
      exit()

class abZ(object):
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

    # Loop over flow data
    for line in f:
      li = line.split()
      z = int(li[0])
      self.data[z] = AttrDict()
      self.data[z].z = z
      self.data[z].yz  = float(li[1])

class abA(object):
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

    # Loop over flow data
    for line in f:
      li = line.split()
      a = int(li[0])
      self.data[a] = AttrDict()
      self.data[a].f = AttrDict()
      self.data[a].a = a
      self.data[a].ya  = float(li[1])
