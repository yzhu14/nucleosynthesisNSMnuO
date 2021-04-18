from hnunu.albinoBasics import *

tracers = "dirk"
fileID, dirt = get_tracer_id(tracers)
dirt_neufield = '../inoutput/neuField/'
dirt_neuInput = '../inoutput/sf2input/'
newElist = [16]
start = "0"
whichFileP = 'v_potential'
whichFileN = 'v_density'

########################################################################
for nNew in newElist:
    print("Generating neutrino files with: ",nNew, "energy bins")
    for id in fileID:
        print("tracer ID:",id)
        ntestpointsN = sum(1 for line in open(
            dirt_neufield + whichFileN + str(id) + '.txt'))
        ntestpointsP = sum(1 for line in open(
            dirt_neufield + whichFileP + str(id) + '.txt'))
        ntestpoints = min(ntestpointsN, ntestpointsP)
        print(ntestpoints)
        # I just happen to know that there are seven comment lines.
        print(dirt_neufield + '/' + str(nNew) + 'NE/' + str(id) +
              'intpCheckffcub' + whichFileP+note + '.txt')
        Morebins(dirt_neufield, id, nNew, whichFileP, whichFileN, note, ntestpoints)
        input(dirt_neufield, id, nNew, whichFileP,
              whichFileN, note, tracers, nCommentLines)
        # ve(id,nNew,whichFileP,note)
        # neucap(dir,id,nNew,whichFileP,whichFileN,note)
