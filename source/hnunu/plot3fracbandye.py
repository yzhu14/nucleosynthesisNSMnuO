from light import *
from plot import *
from brokenaxes import brokenaxes
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
matplotlib.use('TKAgg')

font = {'size': 6,}
if __name__ == "__main__":
    plt.rcParams["figure.subplot.right"] = 0.99
    yelabels = [yelabel(i) for i in notelist]
    colors = = [color(i) for i in notelist]
    for date in datelist:
        for massmodel in range(0, len(masslist)):
            for fissionfile in fission:
                fig = plt.figure()
                plt.ylabel("Fraction of Effective Heating Rate",fontsize=15)
                plt.yticks([])
                plt.xticks([])
                gs = gridspec.GridSpec(3, 1, figure=fig)
                ax1 = fig.add_subplot(gs[0, 0])
                ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
                ax3 = fig.add_subplot(gs[2, 0], sharex=ax2)
                ax1.set_ylim(9e-3, 1.005)
                ax2.set_ylim(9e-3, 1.005)
                ax3.set_ylim(9e-3, 1.005)
                ax1.yaxis.set_label_position("right")
                ax1.yaxis.tick_right()
                ax2.yaxis.set_label_position("right")
                ax2.yaxis.tick_right()
                ax3.yaxis.set_label_position("right")
                ax3.yaxis.tick_right()
                yticks = [0.1, 0.5,0.9]
                ax1.set_yticks(yticks)
                ax2.set_yticks(yticks)
                ax3.set_yticks(yticks)
                ax1.tick_params(axis='both', which='major', labelsize=14)
                ax2.tick_params(axis='both', which='major', labelsize=14)
                ax3.tick_params(axis='both', which='major', labelsize=14)
                xticks =[0.1,1,10,100,1000]
                ax3.set_xticks(xticks)
                labelshow = True
                for betaf in betafiles:
                    for note in notelist:
                        eventf = masslist[massmodel]+betaf+date + \
                            fissionfile + get_barrier(masslist[massmodel]) + note
                        print(eventf)
                        dirout = dirt+"output/"+date+"/"+betaf+"/"
                        x, frac = get_frac(dirout, eventf)
                        i = 0  # ad
                        ax1.semilogx(x, frac.iloc[:, i], label=get_channel_label(
                            i), linestyle=get_style(i), color=getcolor(note), alpha=0.84)
                        i = 5  # spf
                        ax2.semilogx(x, frac.iloc[:, i], label=get_channel_label(
                            i), linestyle=get_style(i), color=getcolor(note), alpha=0.84)
                        i = 6
                        ax3.semilogx(x, frac.iloc[:, i], label=get_channel_label(
                            i), linestyle=get_style(i), color=getcolor(note), alpha=0.84)
                    filename = 'frac'+fissionfile+masslist[massmodel]+'.csv'
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
                    h = [plt.plot([], [], colorlist4[i], ls='-', alpha=0.5)[0] for i in range(8)]
                    plt.legend(handles=h, labels=yelabel, loc='upper center',bbox_to_anchor=(0.5, 3.500), 
                            fontsize=5.5, ncol=8, fancybox=False, shadow=False)
                    plt.savefig(filename + '.eps',format='eps',dpi=1000,rasterized=True)
                    plt.close()
