from light import *
from ab import *
from plots import *
#plt.rc('text', usetex=True, fontsize=14)
plt.rc('font', family='serif')
plt.rc("ytick", direction="in")
plt.rc("xtick", direction="in")
meject = 0.05*1.98847e33 #g
def heatandbandevent12(eventmark, fig, ax1, fission, datelist, notelist, betafiles, masslist, labellist, linecolor, bandcolor='#d9d9d9'):
    event, dirout = get_events(datelist, notelist, betafiles, masslist, fission, dirt)
    #event = eventmark
    print(len(event))
    y = []
    ymin = []
    ymax = []
    ymineff = []
    ymaxeff = []
    yevent = [0]*len(event)
    yeventeff = [0]*len(event)
    x = get_xaxis_tday(dirout, event[0])
    heatcum = [0]*len(event)
    revent = []
    j = 0
    for i in range(len(event)):
        try:
            statinfo = os.stat(dirout+event[i]+'cum_heating_rates.txt')
            if statinfo.st_size > 1:
                heatcum[j] = totalHeat(dirout+event[i]+'cum_heating_rates.txt')
                yeventeff[j] = []
                yevent[j] = []
                revent.append(event[i])
                j += 1
        except:
            continue

    for t in x:
        sumheff = [0.0]*len(revent)
        sumh = [0.0]*len(revent)
        for i in range(len(revent)):
            ts = heatcum[i].find_m_closest_to_t(t*daytosec)
            for rn in cnlist:  # heatcum[i].data[ts].H.keys():
                factor = thermalf(t,t_tha,t_thg,t_thf,t_the,rn)
                sumh[i] += heatcum[i].data[ts].H[rn]*meject
                sumheff[i] += heatcum[i].data[ts].H[rn]*factor*meject
            yevent[i].append(sumh[i])
            yeventeff[i].append(sumheff[i])
        ymin.append(min(sumh))
        ymax.append(max(sumh))
        ymineff.append(min(sumheff))
        ymaxeff.append(max(sumheff))
        y.append(np.mean(sumh))
    ax1 = plot_events(eventmark,revent, x, yeventeff, labellist, ax1,
                      ymineff, ymaxeff, 0.05, linecolor,'#d9d9d9')
    return fig, ax1
####################################################################################################################
if __name__ == "__main__":
    #marklabel = [r"HFB22, $Y_e = 0.16$",r"FRDM2012, $Y_e = 0.28$",r'DFT_SLY4, $Y_e = 0.21$',r'DFT_UNEDF1, $Y_e = 0.02$']
    marklabel = [1,2,3,4,5,6,7,8,9,10,11,12]
    eventcolor = ['#fa8a76','#e15383','#ffd700','#f29724','#41ab5d','#c6dbef','#6baed6','#2171b5','#d1afe8','#63589f',"#006d2c",'#08306b']#'#fee08b','#bf812d']
    eventcolor = ['#fa8a76','#e15383','#ffd700','#f29724','#41ab5d','#c6dbef','#6baed6','#2171b5','#d1afe8','#63589f','#b2df8a','#08306b']
    #eventcolor = ['#fa8a76','#e15383','#ffd700','#f29724','#6cc08b','#034e7b','#3690c0','#a6bddb','#d1afe8','#63589f']
    fig, ax1 = plt.subplots()
    fig, ax1 = heatandbandevent12(eventmark, fig, ax1, fission, datelist,notelist, betafiles, masslist, marklabel,eventcolor)
    plt.grid(True)
    ax1.tick_params(axis='both', which='major', labelsize=14)

    ax1.set_xlim(0.1,1e3)
    ax1.yaxis.set_label_position("left")
    plt.subplots_adjust(top = 0.95,bottom=0.10,right= 0.95) 
    ax1.yaxis.tick_left()
    plt.xlabel('Time (Day)',fontsize=14)
    plt.ylabel(r'Effective Heating Rate (erg$\,s^{-1}$)',fontsize=14)
    # #######
    # latex = np.array([i for i in np.linspace(25,150)])
    # latex1 = latex[:150]
    # plt.plot(latex1,1e13*np.array(latex1)**(-4/3),'r--')
    # x1 = latex1[0]
    # plt.text(x1*3,4e12*(x1)**(-4/3),r"L$\propto t^{-4/3}$",fontsize = 12,color='r')
    # latex2 = latex[-150:]
    # x2 = latex2[0]
    # plt.plot(latex2,3e11*np.array(latex2)**(-7/3),'r-')
    # plt.text(x2*0.6,1e10*(x2)**(-7/3),r"L$\propto t^{-7/3}$",fontsize = 12,color ='r')
    ######
   
    #ax.set_rasterized(True)
########
    dirmix = '/Users/yzhu14/academicmap/code/fixednuc/mixedmodels/'
    eventmark = ['dz33kodamaKarpov101819mix','dft_unedf1asyKarpov101819mix']
    mixlabel = [11,12]

    revent, x, yeventeff, ymineff, ymaxeff = heatandbandevent(eventmark,fission, datelist, notelist, betafiles, masslist)
    xmix, yeventeffmix=heatandbandmix(eventmark,dirmix)
    ax1= plot_mix(eventmark,x,xmix,yeventeffmix,mixlabel,ax1,ymineff,ymaxeff,0.01,eventcolor[-2:],'#d9d9d9')
    ax1.set_ylim(meject*3e3,meject*1.2e12)
# ##########

#     dirtf = '/Users/yzhu14/Research/nucleosynthesis/forYonglin/output/040820/'
#     eventtf = ['tfgtffsdn040820asyKarpovFRLDMTFetfsiye16']
#     tflabel = ['TF']
#     revent, x, yeventeff, ymineff, ymaxeff = heatandbandevent(eventtf,fission, datelist, notelist, betafiles, masslist)
#     xtf, yeventefftf = heatandbandmix(eventtf,dirtf)
#     ax1 = plot_mix(eventtf,x,xtf,yeventefftf,tflabel,ax1,ymineff,ymaxeff,0.01,'blue','black')
#     ax1.set_ylim(meject*3e3,meject*1.2e12)
##########
    latex = x[110:200]
    plt.plot(latex,9e8*np.array(latex)**(-1)*meject,'r-.')
    x1 = latex[0]
    plt.text(x1,1e8*(x1)**(-1)*meject,r"$\dot { Q } \propto t^{-1}$",fontsize = 12,color='r')
    ##################################################################################################
    latex = x[-130:-50]
    plt.plot(latex,1.2e11*np.array(latex)**(-1.1)*meject,'r--')
    x1 = latex[0]
    plt.text(x1*2,9.5e10*(x1)**(-1.1)*meject,r"$\dot { Q } \propto t^{-1.1}$",fontsize = 12,color='r')
    latex = x[-100:-20]
    x1 = latex[10]
    plt.plot(latex,9e9*np.array(latex)**(-2.3)*meject,'r-')
    plt.text(x1*0.6,1e9*(x1)**(-2.3)*meject,r"$\dot { Q } \propto t^{-2.3}$",fontsize = 12,color ='r')
    # for e,eventf in enumerate(eventmark):
    #     mfz = mfZ('../massfraction/'+eventf+'massfraction'+str(mfday)+'.csv')
    #     mfla = 0.
    #     mfac = 0.
    #     for z in mfz.data.keys():
    #         if (z > 56 and z < 72):
    #             mfla += mfz.data[z].mfz
    #         elif (z > 88 and z < 104):    
    #             mfac += mfz.data[z].mfz
    #     print eventf,mfla,mfac
    #     a.bar(e,mfla+mfac,color=eventcolor[e+10])
    #     a.text(e-0.05,mfla+mfac+0.02,str(mixlabel[e]))
    plt.legend(bbox_to_anchor=(0.33,0.79), prop={'size': 12},ncol=4, shadow=False, fancybox=False)
    a = plt.axes([.17, .12, .4, .35])
    mfday = 1
    eventmark = ["frdm2012gtffsdn101819asyKarpovFRLDMye28",
        "frdm2012gtffsdn101819asyKarpovFRLDMye16",
        "hfb22gtffsdn101819asyKarpovHFBye16",
        "hfb27gtffsdn101819kodamaKarpovHFBye16",    
        "dz33gtffsdn101819asyKarpovFRLDMye16",
        "dft_unedf1gtffsdn101819kodamaKarpovFRLDMye16",
        "dft_unedf1gtffsdn101819kodamaXurenFRLDMye16",
        "dft_unedf1gtffsdn101819asyKarpovFRLDMye24",
        "dft_sly4gtffsdn101819asyKarpovFRLDMye18",
        "dft_sly4gtffsdn101819asyKarpovFRLDMye21"]
    for e,eventf in enumerate(eventmark+['dz33kodamaKarpov101819mix','dft_unedf1asyKarpov101819mix']):
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
    plt.ylabel(r'$X_{(Ln+An)}(t_{\kappa})$')
    plt.xticks([])
    plt.grid(True)
    plt.ylim(0.001, 0.5)
    plt.yticks([0.1,0.2,0.3,0.4])
    plt.savefig('./figs/heatbandeventmfbar.eps',format='eps',dpi=1000,rasterized=True)
    plt.close()
