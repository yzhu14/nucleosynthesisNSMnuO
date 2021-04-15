from light import *
from plots import *

if __name__ == "__main__":
    for note in notelist:
            fig, ax1 = plt.subplots()
            plt = heatandband(fission, datelist, [note],
                            betafiles, masslist, '')
            plt.legend(loc="best", prop={'size': 10}, shadow=False, fancybox=False)
            plt.grid(True)
            plt.xlabel('Time (Day)',fontsize=14)
            plt.ylabel('Effective Heating Rate (erg/sec)',fontsize=14)
        ########################################
            plt.grid(True)
            plt.savefig('./figs/heatingband'+note+'.eps',format='eps',dpi=1000,rasterized=True)
            plt.close()
