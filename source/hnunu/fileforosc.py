__author__ = 'yzhu14'
#March 2017
from albinoBasics import *
########################################################################
def get_tracer_id(tracers):
    if tracers=="xy70":
        fileID =["14746", "14782","15364","15377","15378","15379","15388","15400",
                 "15407","15412","15414","17090","17096","17704","17708","17721",
                 "17779","17781","17783","17786","17789","17792","22290"]
        dir="/Users/yzhu14/Research/hnunu/tests/data/neuCap/xy70/forNueosc"
    elif tracers=="dirk":
        fileID=['57221',"61778",'66394','78323','79049',"80224"]
        dir="/Users/yzhu14/Research/hnunu/tests/data/neuCap/dirksample/forNueosc"
    elif tracers=="paperOne":
        fileID=['3']
        dir='/Users/yzhu14/Research/hnunu/data/tracers/paperOne'
    else:
        print("Wrong Dir!")
    return fileID,dir
#######################################################################

tracers="paperOne"
fileID,dir=get_tracer_id(tracers)
#dir='/Users/yzhu14/Research/hnunu/tests/data/tracers/testcos/'
newElist=[16]
start ="0"
whichFileP = 'potential'
whichFileN = 'density'
dir="/Users/yzhu14/Research/hnunu/tests/data/psma/test1/inputfornupotential/17090/"
idrange=443

########################################################################
for nNew in newElist:
    print(nNew,"energy bins")
    #for id in fileID:
    for id in range(0,idrange+1):
        print (id)
        ntestpointsN = sum(1 for line in open(dir +'/' + whichFileN + str(id) + '.txt'))
        ntestpointsP = sum(1 for line in open(dir +'/' + whichFileP + str(id) + '.txt'))
        ntestpoints =min(ntestpointsN,ntestpointsP)
        print(ntestpoints)
        # I just happen t know that there are seven comment lines.
        print(dir +'/' + str(nNew) + 'NE/'  + str(id) + 'intpCheckffcub' + whichFileP+note + '.txt')
        Morebins(dir,id,nNew,whichFileP,whichFileN,note,ntestpoints)
        input(dir,id,nNew,whichFileP,whichFileN,note,tracers,nCommentLines)
        #ve(id,nNew,whichFileP,note)
        #neucap(dir,id,nNew,whichFileP,whichFileN,note)
