from plots import *
from light import *
from ab import *
from brokenaxes import brokenaxes
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator


colorlist4 = ['#b80d48', '#f29724', '#2b6a6c', '#ff3300','#ffd700', '#414ac5', '#00ffff', '#02a229']

if __name__ == "__main__":
    eventlabel = [1,2, 4, 5, 6, 7, 8,10,11,12]
    filename = 'frac3event'
    eventmark = []
    plt.rcParams["figure.subplot.right"] = 0.99
    dirout = dirt + "output/"+datelist[0] +"/"+betafiles[0]+"/"
    fission = np.array(["asyKarpov", "kodamaKarpov", "asyXuren", "kodamaXuren"])

    fissionlabel = (["Symmetric, KZ","Kodama, KZ",  "Symmetric, XR", "Kodama, XR",])
    for m,mass in enumerate(masslist):    
        for n,note in enumerate(notelist):
            if note=='ye16':
                fig = plt.figure()#constrained_layout=True)
                plt.yticks([])
                plt.xticks([])
                gs = gridspec.GridSpec(3, 1, figure=fig)
                ax1 = fig.add_subplot(gs[0, 0])
                ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
                plt.ylabel("Fraction of Effective Heating Rate",fontsize=15)
                ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
                ax3.set_xlim(0.1,1001.01)
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

                colors = []
                for f,fissionfile in enumerate(fission):
                    print fissionfile,fissionlabel[f]
                    eventf = mass+'gtffsdn'+date+fissionfile+get_barrier(mass)+note
                    x, frac = get_frac(dirout, eventf)
                    i = 0  # ad
                    ax1.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color=colorlist4[2*f], alpha=0.64)
                    i = 5  # spf
                    ax2.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color=colorlist4[2*f], alpha=0.64)
                    i = 6
                    ax3.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color=colorlist4[2*f], alpha=0.64)
                    colors.append(2*f)
                # ################# no cf yields
                # dirout2 = dirt + "output/122819/"+betafiles[0]+"/"

                # for f,fissionfile in enumerate(fission):
                #     eventf = mass+'gtffsdn'+'122819'+fissionfile+get_barrier(mass)+note
                #     x, frac2 = get_frac(dirout2, eventf)
                #     i = 0  # ad
                #     ax1.semilogx(x, frac2.iloc[:, i], linestyle=linestyles[2], color=colorlist4[-(1+f)], alpha=0.64)
                #     i = 5  # spf
                #     ax2.semilogx(x, frac2.iloc[:, i], linestyle=linestyles[2], color=colorlist4[-(1+f)], alpha=0.64)
                #     i = 6
                #     ax3.semilogx(x, frac2.iloc[:, i], linestyle=linestyles[2], color=colorlist4[-(1+f)], alpha=0.64)
                #     colors.append(-(1+f))

                fontsize = 14
                fontweight = 'bold'
                fontproperties = {'weight' : fontweight, 'size' : fontsize}
                ax1.set_xticklabels(ax1.get_xticks(), fontproperties)
                ax1.set_yticklabels(ax1.get_yticks(), fontproperties)

                ax1.text(11e-2, 0.88,r"\textbf{Alpha Decay}", fontsize=10)
                ax2.text(11e-2,0.88,r"\textbf{SP Fission}", fontsize=10)
                ax3.text(11e-2,0.1,r"\textbf{Beta Decay}", fontsize=10)
                ax3.set_xlabel("Time (Day)", fontsize=15)
        
                plt.setp(ax2.get_xticklabels(), visible=False)
                plt.setp(ax1.get_xticklabels(), visible=False)
                plt.tight_layout()
                h = [plt.plot([], [], colorlist4[i], ls='-', alpha=0.95)[0] for i in colors]
                newfissionlabel = [i+str('no254') for i in fissionlabel]
                plt.subplots_adjust(left = 0.1,right = 0.95,hspace=0.03,top = 0.89,bottom=0.11)     # the top of the subplots of the figure
                plt.legend(handles=h, labels=fissionlabel+newfissionlabel, loc='upper center',bbox_to_anchor=(0.5, 3.54), 
                        fontsize=10, ncol=2, fancybox=False, shadow=False)
                plt.savefig('./figs/fissionfrac'+mass+note + '.eps',format='eps',dpi=1000,rasterized=True)
                plt.close()
