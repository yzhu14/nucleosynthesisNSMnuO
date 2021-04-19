from hnunu.network import abA
from hnunu.prismBasic import *
dirt = '../../data/'
solarA = abA(dirt+"yva_solar_sne08.dat")
solarZ = abA("../../data/"+"solarabz.txt")
date = datelist[0]
minA = 120
maxA = 249
daytosec = 86400
ycutoff = 1e-30

fitmix = [[4.03150890e-16, 2.82572247e-14, 1.82624236e+03, 2.42513800e-14,4.85836134e+02, 2.74045833e+02, 1.45539315e-18, 1.16426023e+02],
[1.04512508e+02, 1.37117780e-18, 1.18310697e-18, 1.94291911e-16, 7.09143652e+02, 2.08542758e-20, 8.92243015e+02, 5.69904325e+02]]

def mfsumscale(prismab, solar, ma, mi):
    init = False
    sumprism = []
    sumsolar = []
    solarmf = []
    az = []
    prismy = []
    prismmf = dict()
    for a in sorted(prismab.data.keys()):
        if a < maxA + 1:
            try:
                yap = float(solar.data[a].ya*a)
            except:
                yap = ycutoff
            sumsolar.append(yap)
            for z, a in sorted(prismab.data[a].y.keys()):
                try:
                    sumprism.append(float(prismab.data[a].y[z, a])*a)
                except:
                    sumprism.append(0.)
                if z < ma+1 and z > mi-1:
                    try:
                        prismmf[z] = prismmf[z]+prismab.data[a].y[z, a]*a
                    except:
                        prismmf[z] = dict()
                        prismmf[z] = prismab.data[a].y[z, a]*a
    for z in range(110):
        az.append(z)
        try:
            prismmf[z] = prismmf[z]/sum(sumprism)
        except:
            prismmf[z] = dict()
            prismmf[z] = ycutoff/sum(sumprism)
        prismy.append(prismmf[z])
    return az, prismy, prismmf
#####################################################

def solarscale(notelist, fission, betafiles, masslist, solar, fitma, fitmi):
    # when plotting raw ab from prism
    # scale solar with sum of a range
    # for example, For ye21, A188-198
    sumsolar = 0.0
    sumprism = 0.0
    for note in notelist:
        for fissionfile in fission:
            for betaf in betafiles:
                for massmodel in masslist:
                    eventf = massmodel+betaf+date+fissionfile+get_barrier(massmodel) + note
                    prismy = abA(dirt+"output/"+date+"/" +
                                 betaf+"/"+"abA"+eventf)
                    for a in sorted(prismy.data.keys()):
                        if fitmi <= a <= fitma:
                            sumprism += prismy.data[a].ya
    sumprism /= len(notelist)*len(fission)*len(betafiles)*len(masslist)
    for a in sorted(solar.data.keys()):
        if fitmi <= a <= fitma:
            sumsolar += solar.data[a].ya
    return sumsolar/sumprism


def sumscalenew(notelist, fission, betafiles, masslist, prismab, solar, fitma, fitmi):
    solarfactor = solarscale(notelist, fission, betafiles,
                             masslist, solar, fitma, fitmi)
    print(solarfactor)
    #solarfactor = 1
    init = False
    sumprism = []
    sumsolar = []
    solary = []
    aa = []
    prismy = []
    prisma = dict()
    for a in range(300):
        if a < maxA + 1:
            try:
                yas = float(solar.data[a].ya)
            except:
                yas = ycutoff
            try:
                yap = float(prismab.data[a].ya)
            except:
                yap = ycutoff
            prismy.append(float(yap))
            aa.append(a)
            solary.append(float(yas)/solarfactor)
        # if a < fitma+1 and a > fitmi-1:
        #     sumsolar.append(float(solar.data[a].ya))
        #     sumprism.append(float(prismab.data[a].ya))
    prismy = np.array(prismy)  # *(sum(sumsolar)/sum(sumprism))
    return aa, prismy, solary
