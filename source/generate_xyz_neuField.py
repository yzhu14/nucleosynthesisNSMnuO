
import numpy as np
from hnunu.albinoBasics import get_tracer_id
from hnunu.nu_potential_tools import generate_tracer_nupotential



if __name__ == "__main__":

    diroutput = "../inoutput/neuField/"
    commentline = 2
    stepAdjust = 50
    dirkID, dirtracer= get_tracer_id('dirk')
    ###########
    # pickSample(stepAdjust)
    print(dirtracer)
    outfilelist = open(dirtracer + "tracerList.txt", "w")
    deleteTracer = np.zeros((10000, 1))
    zmax = 500
    for id in dirkID:
        nline = generate_tracer_nupotential(id, zmax,dirtracer,diroutput)
        outfilelist.write(str(id) + "\t" + str(zmax) + "\t" + str(nline) + "\n")
    outfilelist.close()
