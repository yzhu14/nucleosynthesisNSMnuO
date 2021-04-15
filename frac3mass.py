from plots import *
from light import *
from brokenaxes import brokenaxes
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
matplotlib.use('TKAgg')

font = {'size': 6,}
if __name__ == "__main__":
    plt.rcParams["figure.subplot.right"] = 0.99
    masslabels = [masslabel(masslist,i) for i in masslist]
    colors = [getcolor(i) for i in masslist]
    notelist = ['ye16']
    for date in datelist:
        for note in notelist:
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
                    for massmodel in range(0, len(masslist)):
                        eventf = masslist[massmodel]+betaf+date + \
                            fissionfile + get_barrier(masslist[massmodel]) + note
                        print(eventf)
                        dirout = dirt+"output/"+date+"/"+betaf+"/"
                        x, frac = get_frac(dirout, eventf)
                        i = 0  # ad
                        ax1.semilogx(x, frac.iloc[:, i], label=get_channel_label(
                            i), linestyle=get_style(i), color=getcolor(masslist[massmodel]), alpha=0.94)
                        i = 5  # spf
                        ax2.semilogx(x, frac.iloc[:, i], label=get_channel_label(
                            i), linestyle=get_style(i), color=getcolor(masslist[massmodel]), alpha=0.94)
                        i = 6
                        ax3.semilogx(x, frac.iloc[:, i], label=get_channel_label(
                            i), linestyle=get_style(i), color=getcolor(masslist[massmodel]), alpha=0.94)
                    ##########################################
                    dirtf = '/Users/yzhu14/Research/nucleosynthesis/forYonglin/output/040820/'
                    eventtf = 'tfgtffsdn040820asyKarpovFRLDMTFmktetfsiye16'
                    x, frac = get_frac(dirtf, eventtf)
                    i = 0  # ad
                    ax1.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color='black', alpha=0.94)
                    i = 5  # spf
                    ax2.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color='black', alpha=0.94)
                    i = 6
                    ax3.semilogx(x, frac.iloc[:, i], linestyle=get_style(i), color='black', alpha=0.94)

                    ############################################
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
                    colors.append('black')
                    masslabels.append('TFmkt')
                    h = [plt.plot([], [], colors[i], ls='-', alpha=0.5)[0] for i in range(len(masslabels))]
                    plt.subplots_adjust(left = 0.1,right = 0.95,hspace=0.03,top = 0.89,bottom=0.11)     # the top of the subplots of the figure
                    plt.legend(handles=h, labels=masslabels, loc='upper center',bbox_to_anchor=(0.5, 3.53), 
                            fontsize=11, ncol=5, fancybox=False, shadow=False)
                    plt.savefig('./figs/frac'+fissionfile+note+'.eps',format='eps',dpi=1000,rasterized=True)
                    plt.close()
