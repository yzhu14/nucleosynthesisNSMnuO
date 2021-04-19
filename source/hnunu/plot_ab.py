from light import *
from prism_plots import *
from abundance import *
from brokenaxes import brokenaxes
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
#####################################################
if __name__ == "__main__":
    print(3)
    #################################################################################
    fig = plt.figure()
    ddt = 0.5
    gs = gridspec.GridSpec(1, 1, figure=fig)
    ax1 = fig.add_subplot(gs[0, :])

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    dir_prismoutput = '../../inoutput/prismOut/'
    eventmark = ['_base57221']
    eventcolor = ['#fa8a76','#e15383','#ffd700','#f29724','#6cc08b','#034e7b','#3690c0','#a6bddb','#d1afe8','#63589f']
    for i,eventf in enumerate(eventmark):
        prismya = []
        solarya = []
        aa = []
        prismf = "abA"+eventf
        prismabA = abA(dir_prismoutput+ prismf)
        prismyz = []
        solaryz = []
        zz = []
        prismfz = "ab"+eventf
        prismabz = ab(dir_prismoutput+prismfz)
        ax1.set_ylabel("Abundance")
        ax1.legend(loc='upper center',fontsize = 12, bbox_to_anchor=(0.5,1.38),ncol=5, fancybox=False, shadow=False)
        print prismabA
        ax1.semilogy(aa, solarya, 'k+', alpha=0.16)
        ax1.set_xlim(minA, maxA)

        d = .015  # how big to make the diagonal lines in axes coordinates
        ax1.tick_params(axis='both', which='major', labelsize=14)

        kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)

        plt.tight_layout()
        ax1.set_xlim(120, 249)

        plt.savefig(eventf+'_ab.eps',format='eps',dpi=1000,rasterized=True)
        plt.close()


