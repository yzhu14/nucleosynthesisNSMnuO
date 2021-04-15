import os
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.widgets import Slider
import matplotlib.animation as animation
from prism_output import Flows,totalHeat,AbTime
from prism_output import Topn,Topnq
import sys, traceback
import json
from scipy.io import FortranFile
from network import network
from network import abA

#plt.rc('text', usetex=True, fontsize=14)
plt.rc('font', family='serif')
plt.rc("ytick",direction="in")
plt.rc("xtick",direction="in")

## TODO:
daytosec = 86400
res = 10
daytosec = 86400
mintime = 0.1e2
maxtime = 2e2
numFrames = 0
nf = 0
# flows
dirt = '/Users/yzhu14/Research/nucleosynthesis/forYonglin/'
betafiles=['gtffsdn']
date = '101819'
td = 10 # day

fission=(["asyKarpov","kodamaKarpov","kodamaXuren","asyXuren"])#"gt_sdo_asy"])#,"gtff_mkt_asy","gtff_mkt_kodama","gtff_sdn_asy","gtff_sdn_kodama"])
channels = 'bmd'#'fission'
yelist = [2,12,14,16,18,21,24,28]#,29,30]
yelist =[16]
masslist = ["etfsi","frdm2012","dz33","dft_sly4","dft_unedf1","hfb22","hfb27","ws3v6"]
masslist = ["dz33","dft_sly4","dft_unedf1"]

notelist = []
for i in yelist:#range(20,30):
  notelist.append("ye"+'{:02}'.format(i))
#####################################################################################################################3
# event=[]
# for note in notelist:
#   for betaf in betafiles:
#     for massmodel in masslist:
#       for fissionfile in fission:
#          eventf= massmodel+betaf+date+fissionfile + note
#          event.append(eventf)
#          dirout=dirt+"output/"+date+"/"+betaf+"/"
#          #print dirout+eventf


def get_barrier(mass):
    if mass in ['hfb22','hfb27']:
       barrier = 'HFB'
    elif mass in ['etfsi']:
       barrier = 'ETFSI'
    else:
       barrier = 'FRLDM'
    return barrier


def getcloset(totalheating,tsec):
    ts = totalheating.find_m_closest_to_t(tsec)
    for i in range(10):
        if totalheating.data[ts+i-5].tsec - tsec <0.01*tsec:
            gettsec = ts + i-5
    return gettsec


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

#####################
nf =0 #number of frame
# =========================================================
#				Output Data Parsing
# =========================================================


# Make the temporary directory
#topheatinglatef = topn(topheatinglate_path)

maxn = 200
maxz = 110
minn = 10
minz = 20

fig = plt.figure(figsize=(7,4))
ax = fig.add_axes([0.1, 0.15, 0.57, 0.7])

ax.grid(True, ":")
ax.set_xlim(60,150)
ax.set_ylim(44,95)
ax.set_ylabel('Z')
ax.set_xlabel('N')
ax.tick_params(axis='both', which='both', top='on', right='on')

ax.plot((10, 70), (20, 20), 'k-', lw=2, alpha=0.25)
ax.plot((20, 80), (28, 28), 'k-', lw=2, alpha=0.25)
ax.plot((60, 120), (50, 50), 'k-', lw=2, alpha=0.25)
ax.plot((120, 190), (82, 82), 'k-', lw=2, alpha=0.25)

ax.plot((28,28), (10,35), 'k-', lw=2, alpha=0.25)
ax.plot((50, 50), (10, 50), 'k-', lw=2, alpha=0.25)
ax.plot((82, 82), (30, 65), 'k-', lw=2, alpha=0.25)
ax.plot((126, 126), (55, 90), 'k-', lw=2, alpha=0.25)
ax.plot((184, 184), (75, 115), 'k-', lw=2, alpha=0.25)


# =========================================================
#					Initialize
# =========================================================


# Uncomment these lines to use a slider
#tax = plt.axes([0.15, 0.93, 0.70, 0.04], axisbg='white')
#tslide = Slider(tax, 'topheating step', 0, nf, valinit=0, valfmt='%1.0f')


# Over-determined range of nuclei
n = np.arange(-0.5, maxn+0.5)
z = np.arange(-0.5, maxz+0.5)
nn, zz = np.meshgrid(n, z)

y = np.array([[np.nan]*maxn]*maxz)
yb = np.array([[np.nan]*maxn]*maxz)

stable_path = dirt + 'input/data/'

# list of stable isotopes
stable = np.genfromtxt(stable_path+'/stable.txt')
for j,row in enumerate(stable):
    y[int(row[0]), int(row[1])] = 1
ym = np.ma.array(y,mask=np.isnan(y))
ax.pcolormesh(n, z, ym, cmap='binary', vmin=0, vmax=1)
cbmdfile =np.genfromtxt(stable_path+'/beta1day.txt')
for j,row in enumerate(cbmdfile):
    yb[int(row[0]), int(row[1])] = 1
yb = np.ma.array(yb,mask=np.isnan(yb))
#ax.pcolormesh(n, z, yb, marker='.',cmap='Blues', vmin=0, vmax=1)
heatmap = ax.pcolormesh(nn, zz, y, cmap='OrRd',vmin=1e3, vmax=1e6,alpha=.53)#,orientation='horizontal')

def ploteach(filef,td):
    getplot = False
    topheatinglatef = Topn(topheatinglate_path+filef)
    y = np.array([[np.nan]*maxn]*maxz, dtype=float)
    for ii in topheatinglatef.data.keys():
      if not getplot and abs(topheatinglatef.data[ii].tday-td)<0.05:
        for (Z,N) in topheatinglatef.data[ii].h.keys():
             y[Z,N] = float(topheatinglatef.data[ii].h[Z,N])
        ym = np.ma.array(y,mask=np.isnan(y))
        heatmap.set_array(ym.ravel())
        cax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
        fig.colorbar(heatmap, cax=cax, label="%s%s" %(filef,topheatinglatef.data[ii].tday))
        ax.xaxis.grid(True, ls=':')
        ax.yaxis.grid(True, ls=':')
        getplot=True


def countTop(dirout,event,td,channels):
  y = np.array([[np.nan]*maxn]*maxz, dtype=float)
  cn = getchannel(channels)
  for filef in event:
    try:
      getheating= False
      totalheating = totalHeat(dirout+filef+"cum_heating_rates.txt")
      topheatinglatef = Topn(topheatinglate_path+filef)
      ts = getcloset(totalheating,td*daytosec)
      for ii in topheatinglatef.data.keys():
        if not getheating and abs(topheatinglatef.data[ii].tday-td)<0.05*td:
          ts = getcloset(totalheating,topheatinglatef.data[ii].tday)
          for (Z,N) in topheatinglatef.data[ii].h.keys():
              if totalheating.data[ts].H[cn]>0:                 
                  if np.isnan(y[Z,N]):
                      y[Z,N] = topheatinglatef.data[ii].h[Z,N]/totalheating.data[ts].H[cn]
                  else:
                      y[Z,N] += topheatinglatef.data[ii].h[Z,N]/totalheating.data[ts].H[cn]
          getheating=True
    except:
      traceback.print_exc()
      print "Skip: %s"%(filef)
      continue 
  ym = np.ma.array(y,mask=np.isnan(y))
  heatmap.set_array(ym.ravel())
  cax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
  fig.colorbar(heatmap, cax=cax, label="Top Counts of %s"%(len(event)))
  ax.xaxis.grid(True, ls=':')
  ax.yaxis.grid(True, ls=':')
  plt.savefig("top%s_count%s_100day" % (channels,filef)+".pdf",bbox_inches="tight")

def flowin(filef,td):
    getplot = False
    topheatinglatef = Topn(topheatinglate_path+filef)
    y = np.array([[np.nan]*maxn]*maxz, dtype=float)
    for ii in topheatinglatef.data.keys():
      if not getplot and abs(topheatinglatef.data[ii].tday-td)<0.05:
        for (Z,N) in topheatinglatef.data[ii].h.keys():
             y[Z,N] = float(topheatinglatef.data[ii].h[Z,N])
        ym = np.ma.array(y,mask=np.isnan(y))
        heatmap.set_array(ym.ravel())
        cax = fig.add_axes([0.82, 0.15, 0.02, 0.7])
        fig.colorbar(heatmap, cax=cax, label="%s%s" %(filef,topheatinglatef.data[ii].tday))
        ax.xaxis.grid(True, ls=':')
        ax.yaxis.grid(True, ls=':')
        plt.savefig("top%s_%s_100day%s" % (channels,filef,ii)+".pdf",bbox_inches="tight")
        getplot=True
#################################################################
def findcommon(event,td,channels,topn):
  topheatinglate_path = dirout+"/"+"top"+channels+"_"
  clist = []
  for ei,ev in enumerate(event):
     newl = []
     topheatinglatef = Topn(topheatinglate_path+ev)
     tn=0
     for ts in topheatinglatef.data.keys():
       if abs(topheatinglatef.data[ts].tday-td)<td*0.05:
          for (Z,N) in topheatinglatef.data[ts].h.keys():
              while(tn < topn):
                  newl.append((Z,N))
                  tn +=1# = float(topheatinglatef.data[ts].h[Z,N])
          if not clist:
             clist = newl
          else:
             clist = list(set(newl).intersection(clist))
  return clist
#################################################################
##################################################################
masslist = ["dz33"]
masslist = ["dz33","dft_sly4","dft_unedf1"]

for betaf in betafiles:
    for fissionfile in fission:
        for note in notelist:
            event = []
            cnlist = []
            czlist = []
            for massmodel in masslist:
                event.append(massmodel+betaf+date+fissionfile + get_barrier(massmodel)+ note)
                filef= massmodel+betaf+date+fissionfile + get_barrier(massmodel)+ note
                dirout=dirt+"output/"+date+"/"+betaf+"/"
                topheatinglate_path = dirout+"/"+"top"+channels+"_"
            
                topheatinglatef = Topnq(topheatinglate_path+filef+'q')
                y = np.array([[np.nan]*maxn]*maxz, dtype=float)
                ii = topheatinglatef.find_m_closest_to_t(td*daytosec)
                fileout = open("toph_"+str(channels)+str(date)+str(fissionfile)+"_"+str(td)+str(note)+str(massmodel)+'day.txt','w')
                fileout.write(str(filef)+"\t"+str(td)+"\n")
                for top in topheatinglatef.data[ii].t.keys():
                    Z = topheatinglatef.data[ii].t[top].z
                    N = topheatinglatef.data[ii].t[top].n
                    y[Z,N] = float(topheatinglatef.data[ii].t[top].h)
                    fileout.write(str(Z)+"\t"+str(N)+"\t"+str(y[Z,N])+"\n")

                fileout.write('\n')
                fileout.close()
                ym = np.ma.array(y,mask=np.isnan(y))
                heatmap.set_array(ym.ravel())
                cax = fig.add_axes([0.10, 0.92, 0.56, 0.02])
                fig.colorbar(heatmap, cax=cax,orientation='horizontal')#label="%s%s" %(filef,topheatinglatef.data[ii].tday))
                ax.xaxis.grid(True, ls=':')
                ax.yaxis.grid(True, ls=':')
                plt.savefig("toph_"+str(channels)+str(date)+str(fissionfile)+"_"+str(td)+str(note)+str(massmodel)+'day.pdf')
            #plt.title(note)
    #    comlist = findcommon(event,1,'bmd',10)
    #    print (date+fissionfile + note),comlist
    #    for z,n in comlist:
    #       cnlist.append(n)
    #       czlist.append(z)
    #    ax.scatter(cnlist,czlist,marker='s',s=10,facecolors='none',edgecolors='green')#marker='o',fillstyle='none')
