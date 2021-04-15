from hnunu.tracerForprism import *
############  MAIN #################################
tracers = "dirk"
fileID, dirIN = get_tracer_id(tracers)

dirOut = "../inoutput/prism/thermoInput/"
dirEx = "../inoutput/prism/extend/"
for id in fileID:
    inFileD = open(dirIN + "tracer_" + str(id) + ".txt", "r")
    ntestpoints = sum(1 for line in open(dirIN + "tracer_" + str(id) + ".txt"))
    # print(dirIN+"tracer_"+str(id)+'.txt',ntestpoints)
    line = inFileD.readlines()
    prep_thermodynamic_input(line,ntestpoints, id, dirOut)
    time = [0.0] * (ntestpoints - commentline)
    R = [0.0] * (ntestpoints - commentline)
############  MAIN #################################
    extend_tracer(R,time,line,id,ntestpoints,dirEx,dirOut)
