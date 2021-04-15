from basic import *
from plots import *
day = 1
mdot = 1
c = 1
mej = 5e-2*mdot# Typical ejecta mass
vej = 0.15*c# typical velocity of ejecta
m5 = mej/(0.01*mdot)
v1 = vej/0.2*c
t_the = 12.9*np.power(m5,2./3.)*np.power(v1,-2)*day
t_tha = 2*t_the
t_thf = 4*t_the
t_thg = 0.3*np.sqrt(m5)*np.power(v1,-1)*day
meject = 0.05*1.98847e33 #g
####################################################################################################################
def f_beta(t,t_the,t_thg):
    edot_g = 0.25
    edot_b = 0.4
    f_b = 1.0/(1.0 + t/t_the)
    f_g = 1.0 - np.exp( -t_thg*t_thg/t/t)
    return edot_g*f_g + edot_b*f_b

def f_alpha(t,t_tha):
    return np.power(1.0+ t/t_tha, -1.5)
   
def f_fiss(t,t_thf):
    return 1.0/(1.0+t/t_thf)

def thermalf(t,t_tha,t_thg,t_thf,t_the,rn):
    if rn==0:
       return f_alpha(t,t_tha)
    elif rn==6:
       return f_beta(t,t_the,t_thg)
    elif rn in [3,4,5]:
       return f_fiss(t,t_thf)
    else:
       return 1
###############
def heatonlyband(plt,fig, ax1,fission, datelist, notelist, betafiles, masslist, labellist, bandcolor):
    event = []
    for date in datelist:
        for note in notelist:
            for betaf in betafiles:
                for massmodel in range(0, len(masslist)):
                    for fissionfile in fission:
                        eventf = masslist[massmodel] + \
                            betaf+date+fissionfile + get_barrier(masslist[massmodel])+note
                        event.append(eventf)
                        dirout = dirt+"output/"+date+"/"+betaf+"/"
    x = []
    y = []
    heatcumtemp = totalHeat(dirout+eventf+"cum_heating_rates.txt")
    for ts in sorted(heatcumtemp.data.keys()):
        t = heatcumtemp.data[ts].tsec
        if t > mintime and t < maxtime:
            x.append(heatcumtemp.data[ts].tsec/daytosec)

    print len(event)
    heatcum = [0]*len(event)
    yevent = [0]*len(event)
    revent = []
    j = 0
    for i in range(len(event)):
        #print i, event[i]
        try:
            statinfo = os.stat(dirout+event[i]+'cum_heating_rates.txt')
            if statinfo.st_size > 1:
                heatcum[j] = totalHeat(dirout+event[i]+'cum_heating_rates.txt')
                yevent[j] = []
                revent.append(event[i])
                j += 1
        except:
            continue
    ymin = []
    ymax = []
    sumh = [0.0]*len(revent)

    for t in x:
        sumh = [0.0]*len(revent)
        for i in range(len(revent)):
            ts = heatcum[i].find_m_closest_to_t(t*daytosec)
            for rn in range(7):  # heatcum[i].data[ts].H.keys():
                sumh[i] += heatcum[i].data[ts].H[rn]
            yevent[i].append(sumh[i])
        y.append(np.mean(sumh))
        ymin.append(min(sumh))
        ymax.append(max(sumh))

    plt.loglog(x, ymin,bandcolor,linestyle='-.',alpha=0.8)
    plt.loglog(x, ymax,bandcolor,linestyle='-',alpha=0.8)
    plt.fill_between(x, ymin, ymax, facecolor=bandcolor,label=labellist, alpha=0.08)
    plt.ylim(min(ymin), max(ymax))
    plt.xlim(mintime/daytosec, maxtime/daytosec)
    return plt
    ########################################
    ########################################
def heatbandeff(fission, datelist, notelist, betafiles, masslist, labellist, bandcolor='#d9d9d9'):

    event = []
    for date in datelist:
        for note in notelist:
            for betaf in betafiles:
                for massmodel in range(0, len(masslist)):
                    for fissionfile in fission:
                      
                        eventf = masslist[massmodel] + \
                            betaf+date+fissionfile + get_barrier(masslist[massmodel])+note
                        event.append(eventf)
                        dirout = dirt+"output/"+date+"/"+betaf+"/"
    x = []
    y = []
    heatcumtemp = totalHeat(dirout+eventf+"cum_heating_rates.txt")
    for ts in sorted(heatcumtemp.data.keys()):
        t = heatcumtemp.data[ts].tsec
        if t > mintime and t < maxtime:
            x.append(heatcumtemp.data[ts].tsec/daytosec)

    print len(event)
    heatcum = [0]*len(event)
    yevent = [0]*len(event)
    revent = []
    ymineff = []
    ymaxeff = []
    yeventeff = [0]*len(event)
    j = 0
    for i in range(len(event)):
        print i, event[i]
        try:
            statinfo = os.stat(dirout+event[i]+'cum_heating_rates.txt')
            if statinfo.st_size > 1:
                heatcum[j] = totalHeat(dirout+event[i]+'cum_heating_rates.txt')
                yevent[j] = []
                yeventeff[j] = []
                revent.append(event[i])
                j += 1
        except:
            continue
    ymin = []
    ymax = []
    sumh = [0.0]*len(revent)


    print len(heatcum)
    print len(revent)
    for t in x:
        sumh = [0.0]*len(revent)
        sumheff = [0.0]*len(revent)

        for i in range(len(revent)):
            ts = heatcum[i].find_m_closest_to_t(t*daytosec)
            for rn in range(7):  # heatcum[i].data[ts].H.keys():
                factor = thermalf(t,t_tha,t_thg,t_thf,t_the,rn)

                sumh[i] += heatcum[i].data[ts].H[rn]*meject
                sumheff[i] += heatcum[i].data[ts].H[rn]*factor*meject
            yevent[i].append(sumh[i])
            yeventeff[i].append(sumheff[i])
        y.append(np.mean(sumh))
        ymin.append(min(sumh))
        ymax.append(max(sumh))
        ymineff.append(min(sumheff))
        ymaxeff.append(max(sumheff))

    for i, ef in enumerate(event):
        if len(notelist)>1:
            mylabel = yelabel(notelist, notelist[i])
            colori = getcolor(notelist[i])
        elif len(masslist)>1:    
            mylabel = masslabel(masslist,masslist[i])    
            colori = getcolor(masslist[i])
        eventf = str(revent[i])
        plt.loglog(x, yeventeff[i], '-', label=mylabel,
                   color=colori, alpha=0.86)


    plt.ylim(min(ymineff), max(ymaxeff))
    plt.xlim(mintime/daytosec, maxtime/daytosec)
    plt.fill_between(x, ymineff, ymaxeff, facecolor=bandcolor,
                     label=labellist, alpha=0.18)

    return plt
    ########################################
    ########################################

    # ###############
def heatandband(addone,fission, datelist, notelist, betafiles, masslist, labellist, bandcolor='#d9d9d9'):
    exdir = addone[0]
    en = len(addone) -1
    exevent = addone[1:]
    event = []
    for date in datelist:
        for note in notelist:
            for betaf in betafiles:
                for massmodel in range(0, len(masslist)):
                    for fissionfile in fission:
                      
                        eventf = masslist[massmodel] + \
                            betaf+date+fissionfile + get_barrier(masslist[massmodel])+note
                        event.append(eventf)
                        dirout = dirt+"output/"+date+"/"+betaf+"/"
    x = []
    y = []
    heatcumtemp = totalHeat(dirout+eventf+"cum_heating_rates.txt")
    for ts in sorted(heatcumtemp.data.keys()):
        t = heatcumtemp.data[ts].tsec
        if t > mintime and t < maxtime:
            x.append(heatcumtemp.data[ts].tsec/daytosec)

    print len(event)
    heatcum = [0]*len(event)
    yevent = [0]*len(event)
    revent = []
    j = 0
    for i in range(len(event)):
        print i, event[i]
        try:
            statinfo = os.stat(dirout+event[i]+'cum_heating_rates.txt')
            if statinfo.st_size > 1:
                heatcum[j] = totalHeat(dirout+event[i]+'cum_heating_rates.txt')
                yevent[j] = []
                revent.append(event[i])
                j += 1
        except:
            continue
    ymin = []
    ymax = []
    sumh = [0.0]*len(revent)

    if len(exevent)>0:
        for i,exevent in enumerate(exevent):
            print exevent
            revent.append(exevent)
            heatcum.append(totalHeat(exdir+exevent+'cum_heating_rates.txt'))
            yevent.append([])
            sumh.append([0.0])
            print "extra\t",exevent
    print len(heatcum)
    print len(revent)
    for t in x:
        sumh = [0.0]*len(revent)
        for i in range(len(revent)):
            ts = heatcum[i].find_m_closest_to_t(t*daytosec)
            for rn in range(7):  # heatcum[i].data[ts].H.keys():
                sumh[i] += heatcum[i].data[ts].H[rn]*meject
            yevent[i].append(sumh[i])
        y.append(np.mean(sumh))
        ymin.append(min(sumh))
        ymax.append(max(sumh))

    for i, ef in enumerate(event):
        if len(notelist)>1:
            mylabel = yelabel(notelist, notelist[i])
            colori = getcolor(notelist[i])
        elif len(masslist)>1:    
            mylabel = masslabel(masslist,masslist[i])    
            colori = getcolor(masslist[i])
        eventf = str(revent[i])
        plt.loglog(x, yevent[i], '-', label=mylabel,
                   color=colori, alpha=0.86)
    ########################################
    if en>1:
        for i,exe in enumerate(exevent):
            print i,exe
            plt.loglog(x, yevent[i+len(event)-1], '-', label=exe[:2], alpha=0.86)
    elif en==1:
        plt.loglog(x, yevent[len(event)], '-', label='TFmkt', alpha=0.86)

    plt.ylim(min(ymin), max(ymax))
    plt.xlim(mintime/daytosec, maxtime/daytosec)
    plt.fill_between(x, ymin, ymax, facecolor=bandcolor,
                     label=labellist, alpha=0.18)

    return plt
    ########################################
    ########################################

###################################################################################################################
def heatandbandmix(event,dirout):
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

    return x, yeventeff
###################################################################################################################
def heatandbandevent11(eventmark, fig, ax1, fission, datelist, notelist, betafiles, masslist, labellist, linecolor, bandcolor='#d9d9d9'):
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
                sumh[i] += heatcum[i].data[ts].H[rn]
                sumheff[i] += heatcum[i].data[ts].H[rn]*factor
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
    ###################################################################################################################
def heatandbandevent(eventmark,fission, datelist, notelist, betafiles, masslist):
    event, dirout = get_events(
        datelist, notelist, betafiles, masslist, fission, dirt)
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
                sumh[i] += heatcum[i].data[ts].H[rn]
                sumheff[i] += heatcum[i].data[ts].H[rn]*factor
            yevent[i].append(sumh[i])
            yeventeff[i].append(sumheff[i])
        ymin.append(min(sumh))
        ymax.append(max(sumh))
        ymineff.append(min(sumheff))
        ymaxeff.append(max(sumheff))
        y.append(np.mean(sumh))

    return revent, x, yeventeff, ymineff, ymaxeff
#################################################################################################