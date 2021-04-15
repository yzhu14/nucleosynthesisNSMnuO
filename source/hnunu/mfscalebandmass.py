from light import *
from plots import *
from ab import *
from brokenaxes import brokenaxes
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
#####################################################
if __name__ == "__main__":
    fileab = open("abundance.csv", 'w')

    dirtf = '/Users/yzhu14/Research/nucleosynthesis/forYonglin/output/040820/'

    eventtf = 'tfgtffsdn040820asyKarpovFRLDMTFmktetfsiye16'
    tflabel = 'TF MKT'
    eventtf = 'tfgtffsdn040820asyKarpovFRLDMTFetfsiye16'
    tflabel = 'TF'
    addone = True

    ddt = 0.5
    for note in ['ye16']:
        for fissionfile in fission:
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
            ax1.set_ylim(1e-7, 1e-2)
            ax3.set_ylim(1.1e-4, 3e-1)
            ax2.set_ylim(1.1e-4, 3e-1)
            plotsolar = False
            showlabel = True
            etf = 0
            maxmfz = dict()
            minmfz = dict()
            maxya = dict()
            minya = dict()
            for betaf in betafiles:
                for massmodel in masslist:
                    eventf = massmodel+betaf+date+fissionfile +get_barrier(massmodel)+ note
                    filemf = open(eventf+"massfraction.csv", 'w')
                    try:
                        m = fission.tolist().index(fissionfile)
                        labels = masslabel(masslist, massmodel)
                        #masslabel[masslist.index(massmodel)]
                        prismya = []
                        solarya = []
                        aa = []
                        ######### add one
                        if addone:
                            prismabA = abA(dirtf+"abA"+eventtf)
                            heatcum = totalHeat(dirtf+eventtf+'cum_heating_rates.txt')
                            ts = heatcum.find_m_closest_to_t(1*daytosec)
                            aa, prismya, solarya = sumscalenew(
                                ['ye21'], ['asyKarpov'], betafiles, masslist, prismabA, solarA, 198, 188)
                            ax1.semilogy(aa, prismya, linewidth=1.0,
                                        linestyle=None, label=tflabel, alpha=0.53,color='blue')
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
                            prismabz = ab(dirtf+"ab"+eventtf)
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
                            ax2.scatter(zz, prismyz, s=24, marker='*',color='blue')
                            ax3.scatter(zz, prismyz, s=24, marker='*',color='blue')
                            for k,z in enumerate(zz):
                                filemf.write(str(z)+"\t"+str(prismyz[k])+"\n")
                            addone = False
                        prismf = "abA"+eventf
                        prismabA = abA(dirt+"output/"+date+"/"+betaf+"/"+prismf)
                        heatcum = totalHeat(
                            dirt+"output/"+date+"/"+betaf+"/"+eventf+'cum_heating_rates.txt')
                        ts = heatcum.find_m_closest_to_t(1*daytosec)
                        aa, prismya, solarya = sumscalenew(
                            ['ye21'], ['asyKarpov'], betafiles, masslist, prismabA, solarA, 198, 188)
                        solarya =np.array(solarya)/sum(fitmix[0])

                        if showlabel:
                            # note[-4:])
                            ax1.semilogy(aa, prismya, linewidth=1.0,
                                        linestyle=None, label=labels, alpha=0.53,color=getcolor(massmodel))
                        else:
                            ax1.semilogy(aa, prismya, linewidth=1.0,
                                        linestyle=None, alpha=0.52,color=getcolor(massmodel))
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
                        prismfz = "ab"+eventf
                        prismabz = ab(dirt+"output/"+date+"/"+betaf+"/"+prismfz)
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
                        ax2.scatter(zz, prismyz, s=24, marker='*',color=getcolor(massmodel))
                        ax3.scatter(zz, prismyz, s=24, marker='*',color=getcolor(massmodel))
                        for k,z in enumerate(zz):
                            filemf.write(str(z)+"\t"+str(prismyz[k])+"\n")

                    except:
                        traceback.print_exc(file=sys.stdout)
                    etf += 1
                    filemf.close()

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
                            color='#d9d9d9', alpha=0.1)
            ax1.semilogy(aa, solarya, 'k+', alpha=0.16)
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
            ax1.legend(loc='upper center',fontsize = 10, bbox_to_anchor=(0.5,1.35),ncol=4, fancybox=False, shadow=False)
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
            plt.subplots_adjust(hspace=0.31,top = 0.89,bottom=0.11) 
            plt.savefig('./figs/101819abAraw'+fissionfile+note+'tf16.eps',format='eps',dpi=1000,rasterized=True)
            plt.close()
    fileab.close()
