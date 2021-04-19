from light import *
from plots import *

daytosec = 86400
if __name__ == "__main__":
    #################################################################################################
    ##################################################################################################
    fig, ax1 = plt.subplots()
    dirmix = '/Users/yzhu14/academicmap/code/fixednuc/mixedmodels/'
    eventmark = ['dz33kodamaKarpov101819mix','dft_unedf1asyKarpov101819mix']
    mixlabel = [11,12]
    eventcolor = ['#fa8a76','#e15383','#ffd700','#f29724','#41ab5d','#c6dbef','#6baed6','#2171b5','#d1afe8','#63589f',"#006d2c",'#08306b']#'#fee08b','#bf812d']

    ##################################################################################################
    revent, x, yeventeff, ymineff, ymaxeff = heatandbandevent(eventmark,fission, datelist, notelist, betafiles, masslist)
    xmix, yeventeffmix=heatandbandmix(eventmark,dirmix)
    ax1= plot_mix(eventmark,x,xmix,yeventeffmix,mixlabel,ax1,ymineff,ymaxeff,0.01,eventcolor[-2:],'#d9d9d9')
    latex = x[110:200]
    plt.plot(latex,9e10*np.array(latex)**(-1),'r-.')
    x1 = latex[0]
    plt.text(x1*2,7e10*(x1)**(-1),r"L$\propto t^{-1}$",fontsize = 12,color='r')
    ##################################################################################################
    latex = x[-130:-50]
    plt.plot(latex,9e10*np.array(latex)**(-1.1),'r--')
    x1 = latex[0]
    plt.text(x1*2,8e10*(x1)**(-1.1),r"L$\propto t^{-1.1}$",fontsize = 12,color='r')
    latex = x[-80:-20]
    x1 = latex[10]
    plt.plot(latex,3e12*np.array(latex)**(-7/3),'r-')
    plt.text(x1*0.6,1e11*(x1)**(-7/3),r"L$\propto t^{-7/3}$",fontsize = 12,color ='r')
    ################################################################################################
    ax1.yaxis.set_label_position("left")
    plt.subplots_adjust(top = 0.95,bottom=0.10,right= 0.95) 
    plt.legend(loc="upper right", prop={'size': 10},ncol=5, shadow=False, fancybox=False)
    ax1.yaxis.tick_left()
    plt.xlabel('Time (Day)',fontsize=14)
    plt.ylabel('Effective Heating Rate (erg/sec)',fontsize=14)
    plt.grid(True)
    ax1.set_xlim(0.1,1000.01)
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
        a.bar(e,mfla+mfac,color=eventcolor[e+10])
        a.text(e-0.05,mfla+mfac+0.02,str(mixlabel[e]))
    a.yaxis.set_label_position("left")
    a.yaxis.tick_right()
    plt.ylabel(r'$w_{(Ln+An)}(t_0)$')
    plt.xticks([])
    plt.grid(True)
    plt.ylim(0.001, 0.5)
    plt.yticks([0.1,0.2,0.3,0.4])
    plt.savefig('./heatbandmixbar.eps',format='eps',dpi=1000,rasterized=True)
    plt.close()
