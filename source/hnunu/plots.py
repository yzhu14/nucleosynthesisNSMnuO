from basic import *
from light import *
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from matplotlib.colors import LogNorm
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import os
from pylab import *

plt.rc('font', family='serif')
plt.rc("ytick", direction="in")
plt.rc("xtick", direction="in")
rc('text', usetex=True)
rc('axes', linewidth=1.1)
rc('font', weight='bold')
eventmark = ["frdm2012gtffsdn101819asyKarpovFRLDMye28",
        "frdm2012gtffsdn101819asyKarpovFRLDMye16",
        "hfb22gtffsdn101819asyKarpovHFBye16",
        "hfb27gtffsdn101819kodamaKarpovHFBye16",    
        "dz33gtffsdn101819asyKarpovFRLDMye16",
        "dft_unedf1gtffsdn101819kodamaKarpovFRLDMye16",
        "dft_unedf1gtffsdn101819kodamaXurenFRLDMye16",
        "dft_unedf1gtffsdn101819asyKarpovFRLDMye24",
        "dft_sly4gtffsdn101819asyKarpovFRLDMye18",
        "dft_sly4gtffsdn101819asyKarpovFRLDMye21"]
notelist = get_notelist()
mixmark = ['dz33kodamaKarpov101819mix','dft_unedf1asyKarpov101819mix']
channellabel1 = ['Alpha Decay', 'Neutron Capture', 'Neutron Emission',
                 'Neutron-induced Fission', 'Beta-delayed Fission', 'SP Fission', 'Beta Decay']
channellabel = ['Beta Decay', 'Alpha Decay', 'SP Fission', '_nolegend_']
markerlist = ['.', 'x', 'X', '*', '+']
# '+'       plus marker
# '*'       star marker
# 'x'       x marker
# '.'       point marker

linestyles = ['-', '--', '-.',':']
# '-'       solid line style
# '--'      dashed line style
# '-.'      dash-dot line style
# ':'       dotted line style
mixcolor = ['#e41a1c','#377eb8']
def getcolor(eventf):
    colorlist11 = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
               '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99']
    colorlist4 = ['#b80d48', '#f29724', '#2b6a6c', '#ff3300','#ffd700', '#414ac5', '#00ffff', '#02a229', '#0f0f0f','#404040']

    try:
        return colorlist4[masslist.index(eventf)]
    except:
        try:
       	     return colorlist11[notelist.index(eventf[-4:])]
        except:
             return colorlist11[fission.tolist().index(eventf)]
    
def geteventcolor(eventf):
    eventcolor = ['#fa8a76','#e15383','#ffd700','#f29724','#41ab5d','#c6dbef','#6baed6','#2171b5','#d1afe8','#63589f',"#006d2c",'#08306b']#'#fee08b','#bf812d']

    #eventcolor = ['#fa8a76','#e15383','#ffd700','#f29724','#6cc08b','#034e7b','#3690c0','#a6bddb','#d1afe8','#63589f',"#6a51a3","#3f007d"]#'#fee08b','#bf812d']
    return eventcolor[eventmark.index(eventf)]
def get_xaxis_tday(dirout, eventf):
    heatcumtemp = totalHeat(dirout+eventf+"cum_heating_rates.txt")
    x = []
    for ts in sorted(heatcumtemp.data.keys()):
        t = heatcumtemp.data[ts].tsec
        if t > mintime and t < maxtime:
            x.append(heatcumtemp.data[ts].tsec/daytosec)
    return x

def plot_mix(eventmark,x,xmix,yeventeffmix,labellist,ax1,ymineff,ymaxeff,balpha,linecolor,bandcolor):
    for i,ev in enumerate(eventmark):
        ax1.loglog(xmix, yeventeffmix[i], '-', color=linecolor[i],label=labellist[i])
    ax1.set_ylim(min(ymineff), max(ymaxeff))
    ax1.set_xlim(mintime/daytosec, maxtime/daytosec)
    ax1.fill_between(x, ymineff, ymaxeff, facecolor=bandcolor, alpha=balpha)
    return ax1

def plot_events(eventmark,revent, x, yeventeff, labellist, ax1, ymineff, ymaxeff, balpha, linecolor,bandcolor):
    print(len(eventmark),revent)

    for i,ev in enumerate(eventmark):
        j = revent.index(ev)
        ax1.loglog(x, yeventeff[j], '-', color=linecolor[i],label=labellist[i])

    ax1.set_ylim(min(ymineff), max(ymaxeff))
    ax1.set_xlim(mintime/daytosec, maxtime/daytosec)
    ax1.fill_between(x, ymineff, ymaxeff, facecolor=bandcolor, alpha=balpha)
    #  ax1.fill_between(x,ymin,ymax,facecolor='red',alpha=0.48)
    return ax1
####################################################################################################################
def get_events(datelist, notelist, betafiles, masslist, fission, dirt):
    event = []
    for date in datelist:
        for note in notelist:
            for betaf in betafiles:
                for massmodel in range(0, len(masslist)):
                    for fissionfile in fission:
                        eventf = masslist[massmodel] + \
                            betaf+date+fissionfile + \
                            get_barrier(masslist[massmodel])+note
                        event.append(eventf)
                        dirout = dirt+"output/"+date+"/"+betaf+"/"
    return event, dirout
####################################################################################################################
def get_style(rn):
    if rn == 0:
        return linestyles[0]
    elif rn == 5:
        return linestyles[1]
    elif rn == 6:
        return linestyles[2]
    else:
        return linestyles[3]

def get_channel_label(rn):
    if rn == 0:
        return channellabel1[0]
    elif rn == 5:
        return channellabel1[5]
    elif rn == 6:
        return channellabel1[6]
    else:
        return '_nolegend_'


def get_frac(dirt, eventf):
    try:
        filefrac = pd.read_csv("fracw"+eventf+".csv")
        x = filefrac.iloc[:, 0]
        frac = filefrac.iloc[:, 1:]
        return x, frac
    except:
        print('calculating fractions')
        x = []
        y = []
        heatcum = totalHeat(dirt+eventf+"cum_heating_rates.txt")
        for ts in sorted(heatcum.data.keys()):
            t = heatcum.data[ts].tsec
            if t > mintime and t < maxtime:
                x.append(heatcum.data[ts].tsec/daytosec)

        yevent = [0]
        yeventeff = [0]
        sumh = [0.0]
        yfrac = dict()
        yfracmax = [0.0]*7
        yfracmin = [0.0]*7
        filefrac = open("frac"+eventf+".csv", 'w')
        frac = np.zeros((len(x), 7))
        for j, t in enumerate(x):
            filefrac.write(str(t)+",")
            sumh = 0.0
            ts = heatcum.find_m_closest_to_t(t*daytosec)
            for rn in range(7):
                factor = thermalf(t,t_tha,t_thg,t_thf,t_the,rn)
                sumh += factor*heatcum.data[ts].H[rn]
            temp = []
            fractemp = []
            factortemp = []
            for rn in range(7):
                factor = thermalf(t,t_tha,t_thg,t_thf,t_the,rn)
                factortemp.append(factor)
                channelfrac = factor*heatcum.data[ts].H[rn]/sumh
                temp.append(factor*heatcum.data[ts].H[rn])
                fractemp.append(channelfrac)
                filefrac.write(str(channelfrac)+",")
            for ft in factortemp:
                filefrac.write(str(ft)+",")
            frac[j] = np.array(fractemp)
            yevent.append(sumh)
            filefrac.write(str(sumh)+",")

            filefrac.write('\n')
        filefrac.close()
        return x, pd.DataFrame(frac)

def plot_frac(dirt,eventf,plt, ax1,ax2,ax3,colors,alphas,linew):
    x, frac = get_frac(dirt, eventf)
    i = 0  # ad
    ax1.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color=colors, alpha=alphas,linewidth=linew)
    i = 5  # spf
    ax2.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color=colors, alpha=alphas,linewidth=linew)
    i = 6
    ax3.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color=colors, alpha=alphas,linewidth=linew)
    return plt, ax1,ax2,ax3

def plotfrac(filename, plt, ax):
    data = pd.read_csv(filename, header=None)
    ax.set_xscale("log")
    ax.stackplot(data[0], data[7], data[1], data[6], data[2]+data[3] +
                 data[4]+data[5], colors=colorlist4, labels=channellabel, alpha=0.35)
    ax.margins(0, 0)
    #plt.legend(loc="upper left", prop={'size': 8}, shadow=False, fancybox=False)
    plt.xlabel('Time (Day)')
    ax.set_xlim(1e-1, 1e3)
    plt.legend(loc='lower left', fontsize=7,
               ncol=1, fancybox=False, shadow=False)
    plt.ylabel('Average Fraction of Effective Heating Rate')
    return plt, ax
