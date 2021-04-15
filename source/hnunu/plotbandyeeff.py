from light import *
from plots import *

if __name__ == "__main__":
    addone = []#['/Users/yzhu14/Research/nucleosynthesis/forYonglin/output/040820/','tfgtffsdn040820asyKarpovFRLDMTFmktetfsiye16']
    fission = ['asyKarpov']
    for note in ['ye24']:
        for fissionfile in fission:
            fig, ax1 = plt.subplots()
            plt = heatbandeff([fissionfile], datelist, [note],
                            betafiles, masslist, '')
            plt.legend(loc="best", prop={'size': 10}, shadow=False, fancybox=False)
            plt.grid(True)
            plt.xlabel('Time (Day)',fontsize=14)
            plt.ylabel(r'Effective Heating Rate (erg$\,s^{-1}$)',fontsize=14)
        ########################################
            plt.grid(True)
            plt.savefig('./figs/heatingband'+fissionfile+note+'eff.eps',format='eps',dpi=1000,rasterized=True)
            plt.close()
