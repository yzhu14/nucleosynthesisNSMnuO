from plots import *
from light import *
from ab import *
from brokenaxes import brokenaxes
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

font = {'size': 6,}

plt.rcParams["figure.subplot.right"] = 0.99
dirout = dirt + "output/"+date +"/"+betafiles[0]+"/"
dirmix = '/Users/yzhu14/academicmap/code/fixednuc/mixedmodels/'
masslist = ['dz33','dft_unedf1']
fission = ['kodamaKarpov','asyKarpov']
colorlist4 = [getcolor(note) for note in notelist]+['k']
mixlabel = [11,12]
yelabels = [yelabel(notelist,note) for note in notelist]
for e,ev in enumerate(mixmark):
        fig = plt.figure()
        plt.yticks([])
        plt.xticks([])
        plt.grid()
        gs = gridspec.GridSpec(3, 1, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
        plt.ylabel("Fraction of Effective Heating Rate",fontsize=15)
        ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
        ax1.set_ylim(9e-3, 1.005)
        ax2.set_ylim(9e-3, 1.005)
        ax3.set_ylim(9e-3, 1.005)
        ax1.yaxis.set_label_position("left")
        ax1.yaxis.tick_left()
        ax2.yaxis.set_label_position("left")
        ax2.yaxis.tick_left()
        ax3.yaxis.set_label_position("left")
        ax3.yaxis.tick_left()
        yticks = [0.1, 0.5,0.9]
        ax1.set_yticks(yticks)
        ax2.set_yticks(yticks)
        ax3.set_yticks(yticks)
        ax1.tick_params(axis='both', which='major', labelsize=14)
        ax2.tick_params(axis='both', which='major', labelsize=14)
        ax3.tick_params(axis='both', which='major', labelsize=14)
        xticks =[0.1,1,10,100,1000]
        ax3.set_xticks(xticks)

        plt, ax1,ax2,ax3 = plot_frac(dirmix,ev,plt, ax1,ax2,ax3,'k',1,2)

        for n,note in enumerate(notelist):
            eventf = masslist[e]+'gtffsdn'+date+fission[e]+get_barrier(masslist[e])+note
            plt, ax1,ax2,ax3 = plot_frac(dirout,eventf,plt, ax1,ax2,ax3,getcolor(note),0.53,1)

        fontsize = 14
        fontweight = 'bold'
        fontproperties = {'weight' : fontweight, 'size' : fontsize}
        ax1.set_xticklabels(ax1.get_xticks(), fontproperties)
        ax1.set_yticklabels(ax1.get_yticks(), fontproperties)

        ax1.text(11e-2, 0.88,r"\textbf{Alpha Decay}", fontsize=10)
        ax2.text(11e-2,0.88,r"\textbf{SP Fission}", fontsize=10)
        ax3.text(11e-2,0.1,r"\textbf{Beta Decay}", fontsize=10)
        ax3.set_xlabel("Time (Day)", fontsize=15)
 

        ax3.set_xlim(0.099,1000.01)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.tight_layout()
        print mixlabel[e]-1
        print len(colorlist4)
        h = [plt.plot([], [], colorlist4[i], ls='-', alpha=0.95)[0] for i in range(9)]
        plt.subplots_adjust(left = 0.1,right = 0.95,hspace=0.03,top = 0.89, bottom=0.11)     # the top of the subplots of the figure

        plt.legend(handles=h, labels=yelabels+[mixlabel[e]], loc='upper center',bbox_to_anchor=(0.47, 3.54), 
                fontsize=11, ncol=5, fancybox=False, shadow=False)

        plt.savefig('./figs/frac'+str(mixlabel[e])+'mix.eps',format='eps',dpi=1000,rasterized=True)
        plt.close()