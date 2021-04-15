from light import *
from plots import *
from ab import *
from brokenaxes import brokenaxes
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
#####################################################
#################################################################################
if __name__ == "__main__":
    fileab = open("abundance.csv", 'w')

    ddt = 0.5
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, figure=fig)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1], sharey=ax2)
    ax1.set_xlabel("A", fontsize=10)
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.set_yscale('log')
    # ax3.set_yscale('log')
    ax3.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_ylim(2e-7, 3e-3)
    ax3.set_ylim(1.1e-4, 3e-1)
    ax2.set_ylim(1.1e-4, 3e-1)
    plotsolar = False
    showlabel = True
    etf = 0
    maxmfz = dict()
    minmfz = dict()
    maxya = dict()
    minya = dict()
    marklabel = [1,2,3,4,5,6,7,8,9,10]
    eventmark = mixmark 
    mixlabel = [11,12]
    eventcolor = ['#fee08b','#bf812d']

    dirmix = '/Users/yzhu14/academicmap/code/fixednuc/mixedmodels/'
    mfday =1
    fitmix = [[4.03150890e-16, 2.82572247e-14, 1.82624236e+03, 2.42513800e-14,4.85836134e+02, 2.74045833e+02, 1.45539315e-18, 1.16426023e+02],
    [1.04512508e+02, 1.37117780e-18, 1.18310697e-18, 1.94291911e-16, 7.09143652e+02, 2.08542758e-20, 8.92243015e+02, 5.69904325e+02]]
    for i,eventf in enumerate(eventmark):
        filemf = open(eventf+"massfraction"+str(mfday)+".csv", 'w')
        try:
            m = i
            labels = mixlabel[i]
            prismya = []
            solarya = []
            aa = []
            prismf = "abA_"+eventf
            prismabA = abA(dirmix+prismf)
            heatcum = totalHeat(
                dirmix+eventf+'cum_heating_rates.txt')
            ts = heatcum.find_m_closest_to_t(1*daytosec)
            aa, prismya, solarya = sumscalenew(
                ['ye21'], ['asyKarpov'], betafiles, masslist, prismabA, solarA, 198, 188)
            solarya =np.array(solarya)/sum(fitmix[0])
            print sum(fitmix[0]),sum(fitmix[1])
            if showlabel:
                # note[-4:])
                ax1.semilogy(aa, prismya, linewidth=1.0,
                            linestyle=None, label=labels, alpha=0.96,color=eventcolor[m])
            else:
                ax1.semilogy(aa, prismya, linewidth=1.0,
                            linestyle=None, alpha=0.96,color=eventcolor[m])
            for i, a in enumerate(aa):
                try:
                    # prismya list; maxya dict()
                    maxya[a] = max(prismya[i], maxya[a])
                except:
                    maxya[a] = dict()
                    maxya[a] = prismya[i]
                try:
                    # prismya list; maxya dict()
                    minya[a] = min(prismya[i], minya[a])
                except:
                    minya[a] = dict()
                    minya[a] = prismya[i]

            prismyz = []
            solaryz = []
            zz = []
            prismfz = "ab_"+eventf
            prismabz = ab(dirmix+prismfz)
            zz, prismyz, prismmf = mfsumscale(prismabz, solarZ, 250, 0)
            for j, z in enumerate(zz):
                try:
                    maxmfz[z] = max(prismmf[z], maxmfz[z])
                except:
                    maxmfz[z] = dict()
                    maxmfz[z] = prismmf[z]
                try:
                    minmfz[z] = min(prismmf[z], minmfz[z])
                except:
                    minmfz[z] = dict()
                    minmfz[z] = prismmf[z]
            ax2.scatter(zz, prismyz, s=24, marker='*',alpha=0.96,color=eventcolor[m])
            ax3.scatter(zz, prismyz, s=24, marker='*',alpha=0.96,color=eventcolor[m])
            for k,z in enumerate(zz):
                filemf.write(str(z)+"\t"+str(prismyz[k])+"\n")
        except:
            traceback.print_exc(file=sys.stdout)
    etf += 1
    showlabel = False
    plota = []
    plotminya = []
    plotmaxya = []
    for a in minya.keys():
        plota.append(a)
        plotminya.append(minya[a])
        plotmaxya.append(maxya[a])
    laplotz = []
    laplotminmfz = []
    laplotmaxmfz = []
    acplotz = []
    acplotminmfz = []
    acplotmaxmfz = []
    for z in minmfz.keys():
        if (z > 56 and z < 72):
            laplotz.append(z)
            laplotminmfz.append(minmfz[z])
            laplotmaxmfz.append(maxmfz[z])
        elif (z > 88 and z < 104):
            acplotz.append(z)
            acplotminmfz.append(minmfz[z])
            acplotmaxmfz.append(maxmfz[z])
    ax1.fill_between(plota, plotminya, plotmaxya,
                    color='#d9d9d9', alpha=0.21)
    ax1.semilogy(aa, solarya, 'k+', alpha=0.6)
    ax1.set_xlim(minA, maxA)

    for i, a in enumerate(laplotz):
        ax2.plot(
            (a, a), (laplotminmfz[i], laplotmaxmfz[i]), color='#d9d9d9', alpha=0.15,zorder=0)
    for i, a in enumerate(acplotz):
        ax3.plot(
            (a, a), (acplotminmfz[i], acplotmaxmfz[i]), color='#d9d9d9', alpha=0.15,zorder=0)
    ax1.set_ylabel("Abundance")
    ax2.set_ylabel("Mass Fraction(Z)")
    ax3.set_xlabel("Z (Actinides)", fontsize=10)
    ax2.set_xlabel("Z (Lanthanides)", fontsize=10)
    ax1.legend(loc='upper center',fontsize = 12, bbox_to_anchor=(0.5,1.3),ncol=5, fancybox=False, shadow=False)
    ax2.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax2.yaxis.tick_left()
    ax3.set_xlim(89-ddt, 103+ddt)
    ax2.tick_params(labelright='off')
    ax3.yaxis.tick_right()
    ax2.set_xlim(57-ddt, 71+ddt)
    d = .015  # how big to make the diagonal lines in axes coordinates
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    ax3.tick_params(axis='both', which='major', labelsize=14)
    kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
    ax2.plot((1-d, 1+d), (-d, +d), **kwargs)
    ax2.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
    kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
    ax3.plot((-d, +d), (1-d, 1+d), **kwargs)
    ax3.plot((-d, +d), (-d, +d), **kwargs)

    plt.tight_layout()
    ax1.set_xlim(120, 249)

    plt.subplots_adjust(hspace=0.3,top = 0.88,bottom=0.1) 
    plt.savefig('./figs/101819abAmix.eps',format='eps',dpi=1000,rasterized=True)
    plt.close()
    filemf.close()
    fileab.close()
