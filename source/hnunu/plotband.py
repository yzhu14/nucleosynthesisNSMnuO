from light import *
from ab import *
from plots import *
#plt.rc('text', usetex=True, fontsize=14)
plt.rc('font', family='serif')
plt.rc("ytick", direction="in")
plt.rc("xtick", direction="in")

####################################################################################################################
if __name__ == "__main__":
    #marklabel = [r"HFB22, $Y_e = 0.16$",r"FRDM2012, $Y_e = 0.28$",r'DFT_SLY4, $Y_e = 0.21$',r'DFT_UNEDF1, $Y_e = 0.02$']
    marklabel = [1,2,3,4,5,6,7,8,9,10]
    eventcolor = [geteventcolor(i) for i in eventmark]
    #eventcolor = ['#fa8a76','#e15383','#ffd700','#f29724','#6cc08b','#034e7b','#3690c0','#a6bddb','#d1afe8','#63589f']
    fig, ax1 = plt.subplots()
    fig, ax1 = heatandbandevent11(eventmark, fig, ax1, fission, datelist,notelist, betafiles, masslist, marklabel,eventcolor)
    plt.grid(True)
    ax1.set_xlim(0.1,1e3)
    ax1.yaxis.set_label_position("left")
    plt.subplots_adjust(top = 0.95,bottom=0.10,right= 0.95) 
    plt.legend(loc="upper right", prop={'size': 10},ncol=5, shadow=False, fancybox=False)
    ax1.yaxis.tick_left()
    plt.xlabel('Time (Day)',fontsize=14)
    plt.ylabel('Effective Heating Rate (erg/sec)',fontsize=14)
    #######
    latex = np.array([i for i in np.linspace(25,150)])
    latex1 = latex[:150]
    plt.plot(latex1,1e13*np.array(latex1)**(-4/3),'r--')
    x1 = latex1[0]
    plt.text(x1*3,4e12*(x1)**(-4/3),r"L$\propto t^{-4/3}$",fontsize = 12,color='r')
    latex2 = latex[-150:]
    x2 = latex2[0]
    plt.plot(latex2,3e11*np.array(latex2)**(-7/3),'r-')
    plt.text(x2*0.6,1e10*(x2)**(-7/3),r"L$\propto t^{-7/3}$",fontsize = 12,color ='r')
    a = plt.axes([.17, .12, .35, .35])
    mfday = 1
    for e,eventf in enumerate(eventmark):
        mfz = mfZ('../massfraction/'+eventf+'massfraction'+str(mfday)+'.csv')
        mfla = 0.
        mfac = 0.
        for z in mfz.data.keys():
            if (z > 56 and z < 72):
                mfla += mfz.data[z].mfz
            elif (z > 88 and z < 104):    
                mfac += mfz.data[z].mfz
        print eventf,mfla,mfac
        a.bar(e,mfla+mfac,color=eventcolor[e])
        a.text(e-0.25,mfla+mfac+0.02,str(e+1))
        #plt.scatter(e,mfac,color=eventcolor[e],marker='*')
    a.yaxis.set_label_position("left")
    a.yaxis.tick_right()
    plt.ylabel(r'$w_{(Ln+An)}(t_0)$')
    plt.xticks([])
    plt.grid(True)
    plt.ylim(0.001, 0.5)
    plt.yticks([0.1,0.2,0.3,0.4])
    plt.savefig('./figs/heatbandeventmfbar.eps',format='eps',dpi=1000,rasterized=True)
    plt.close()
