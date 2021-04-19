from light import *
from plots import *
#r"$\dot { Q } \propto t^{-1}$"

if __name__ == "__main__":
    for mass in masslist:
        for fissionfile in fission:
            fig, ax1 = plt.subplots()
            plt = heatandband([fissionfile], datelist, notelist,
                        betafiles, [mass], '')
            ax1.tick_params(axis='both', which='major', labelsize=14)
            plt.legend(loc="best", prop={'size': 12}, shadow=False, fancybox=False)
            plt.grid(True)
            plt.xlabel('Time (Day)',fontsize=14)
            plt.ylabel(r'Effective Heating Rate (erg$\,s^{-1}$)',fontsize=14)
            plt.subplots_adjust(top = 0.95,bottom=0.11)

            plt.savefig('./figs/heatband'+fissionfile+mass+'.eps',format='eps',dpi=1000, rasterized=True)
            plt.close()

