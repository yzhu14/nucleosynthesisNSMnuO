from plots import *
from light import *
from ab import *
from brokenaxes import brokenaxes
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

font = {'size': 6,}


eventlabel = [1,2,3,4,5,6,7,8,9,10]
eventcolor = [geteventcolor(i) for i in eventmark]
filename = 'frac3event'
yticks = [0.1, 0.5,0.9]
xticks =[0.1, 1 , 10, 100, 1000]
fontsize = 14
fontweight = 'bold'
fontproperties = {'weight' : fontweight, 'size' : fontsize}
#####
plt.rcParams["figure.subplot.right"] = 0.99
dirout = dirt + "output/"+datelist[0] +"/"+betafiles[0]+"/"
fig = plt.figure()#constrained_layout=True)
plt.yticks([])
plt.xticks([])

gs = gridspec.GridSpec(3, 1, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
plt.ylabel("Fraction of Effective Heating Rate",fontsize=15)

ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
ax3.set_xlim(0.1,1e3)
ax1.set_ylim(9e-3, 1.005)
ax2.set_ylim(9e-3, 1.005)
ax3.set_ylim(9e-3, 1.005)

ax1.set_yticks(yticks)
ax2.set_yticks(yticks)
ax3.set_yticks(yticks)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax2.tick_params(axis='both', which='major', labelsize=14)
ax3.tick_params(axis='both', which='major', labelsize=14)
ax3.set_xticks(xticks)

for e,eventf in enumerate(eventmark):
    print(eventf)
    x, frac = get_frac(dirout, eventf)
    i = 0  # ad
    ax1.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color=eventcolor[e], alpha=0.94)
    i = 5  # spf
    ax2.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color=eventcolor[e], alpha=0.94)
    i = 6
    ax3.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color=eventcolor[e], alpha=0.94)
ax1.set_xticklabels(ax1.get_xticks(), fontproperties)
ax1.set_yticklabels(ax1.get_yticks(), fontproperties)
ax1.text(0.11, 0.88,r"\textbf{Alpha Decay}", fontsize=10)
ax2.text(0.11,0.88,r"\textbf{SP Fission}", fontsize=10)
ax3.text(0.11,0.1,r"\textbf{Beta Decay}", fontsize=10)
ax3.set_xlabel("Time (Day)", fontsize=15)

plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax1.get_xticklabels(), visible=False)

ax1.yaxis.set_label_position("left")
ax1.yaxis.tick_left()
ax2.yaxis.set_label_position("left")
ax2.yaxis.tick_left()
ax3.yaxis.set_label_position("left")
ax3.yaxis.tick_left()
plt.tight_layout()

h = [plt.plot([], [], eventcolor[i], ls='-', alpha=0.95)[0] for i in range(10)]
plt.subplots_adjust(left = 0.1,right = 0.95,hspace=0.03,top = 0.89,bottom=0.11)     # the top of the subplots of the figure

plt.legend(handles=h, labels=eventlabel, loc='upper center',bbox_to_anchor=(0.5, 3.535), 
        fontsize=11, ncol=5, fancybox=False, shadow=False)

plt.savefig(filename + '_tfmkt.eps',format='eps',dpi=1000,rasterized=True)
plt.close()
