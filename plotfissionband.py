from light import *
from plots import *


if __name__ == "__main__":
    fig, ax1 = plt.subplots()
    print fission
    print masslist
    i = 0
    eventcolor = ['#fa8a76','#e15383','#ffd700','#f29724','#6cc08b','#034e7b','#3690c0','#a6bddb','#d1afe8','#63589f','#fee08b','#bf812d']
    fission = ['asyKarpov','asyXuren']
    for fissionfile in fission:
        plt = heatonlyband(
            plt,
            fig,
            ax1,
            [fissionfile],
            datelist,
            notelist,
            betafiles,
            masslist,
            fissionfile,
            eventcolor[2*i])
        i+=1
        plt.legend(loc="best", prop={"size": 10}, shadow=False, fancybox=False)
        plt.grid(True)
        fig, ax1 = heatandbandevent11(eventmark, fig, ax1, fission, datelist,notelist, betafiles, masslist, marklabel,eventcolor)
        # plt.xlabel('Time (Day)',fontsize=14)
        # plt.ylabel('Effective Heating Rate (erg/sec)',fontsize=14)
    #plt.show()
    plt.savefig("./figs/heatsymfission.png", dpi=1000,rasterized=True)
    plt.close()

