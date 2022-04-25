import numpy
import matplotlib.pyplot as plt
import scipy.stats
import scipy.special
import scipy.interpolate
import scipy.optimize
import pandas
import os
import re
import copy




def stringToFloatArray (theString):
    preNumbers = re.split(',|\[|\]| ', theString)
    return numpy.array ([float(x) for x in preNumbers if x != ''])


def realSizes (cs, theFactor):
    return theFactor * 0.5/numpy.exp (cs)


def readChimpOutput (filename):
    ifs = open (filename, 'r')
    rescaleFactor = 0
    optTrajectory = []
    currentEstep = None
    currentMsteps = []
    lls = []
    for line in ifs:
        line = line.strip()
        if (line.startswith ('Internal Rescale factor')):
            rescaleFactor = float(line.split()[-1])
        elif (line.startswith ('Epoch of population size change are')):
            sizeChangeTimes = stringToFloatArray (line.split(':')[-1])
        elif (line.startswith ('c_rate weights are')):
            assert (rescaleFactor > 0)
            # see if we have to save stuff
            if (len(currentMsteps) > 0):
                # should have been an E step before this
                assert (currentEstep is not None)
                optTrajectory.append ((currentEstep, copy.deepcopy(currentMsteps)))
                currentMsteps = []
            currentEstep = realSizes (stringToFloatArray (line.split(':')[-1]), rescaleFactor)
        elif (line.startswith ('#PARAMS')):
            currentMsteps.append (realSizes (stringToFloatArray (line.split(':')[-1]), rescaleFactor))
        elif (line.startswith ('after')):
            lls.append (stringToFloatArray (line.split(':')[-1]))
    # append the final E step
    optTrajectory.append ((currentEstep, None))

    return (sizeChangeTimes, optTrajectory, lls)


def prepareForPlotting (sizeChangeTimes, sizes):
    plotTimes = numpy.concatenate (([1e-8], sizeChangeTimes, [1e8]))
    plotSizes = numpy.concatenate (([sizes[0]], sizes))
    return (plotTimes, plotSizes)


def plotGradientTrajectory (times, truthData, steps, cMapName='cool', truthColor='black', xLim = (1e1,1e6), yLim=(1e3,1e7)):
    # now some trajectory
    theseColors = getColors (len(steps), cMapName)
    for (idx, theseSizes) in enumerate(steps):
        (plotTimes, plotSizes) = prepareForPlotting (times, theseSizes)
        plt.step (plotTimes, plotSizes, where='pre', color=theseColors[idx])

    # truth
    if (truthData is not None):
        plt.plot (truthData['t (gens)'], truthData['N(t)'], color=truthColor)

    # some formatting
    plt.xlim (xLim)
    plt.ylim (yLim)
    plt.xscale ('log')
    plt.yscale ('log')

    plt.show()


def getColors (numColors, cMapName):
    return plt.get_cmap(cMapName)(numpy.linspace(0,1,numColors))


# load some output files
preprefix = "00p_"
postprefix = "TH_n210_"
chrom = "chr1"

# load all pops and ls
pops = ["FIN", "JPT", "LWK"]
ls = [-5]
outData = {}
for thisL in ls:
    for pop in pops:
        chimpOutputFile = f"../../analysis_CHIMP/out/{preprefix}l{thisL}_{postprefix}unphased_{chrom}_{pop}/unphased_{chrom}_{pop}/unphased_{chrom}_{pop}.out"
        (sizeChangeTimes, optTrajectory, lls) = readChimpOutput (chimpOutputFile)
        outData[(thisL, pop)] = (sizeChangeTimes, optTrajectory, lls)

# make a nice plot
# from Matthew Hahn
yearsPerGen = 26.9
realL = -5
pops = ["FIN", "JPT", "LWK"]
for pop in pops:
    # get lasst estimate for current pop
    finalSizes = outData[(realL, pop)][1][-1][0]

    # and plot them
    (plotTimes, plotSizes) = prepareForPlotting (sizeChangeTimes, finalSizes)
    plt.step (yearsPerGen * plotTimes, plotSizes, where='pre')

# polish plot
plt.xscale ('log')
plt.yscale ('log')
plt.xlim (yearsPerGen * 1e2, yearsPerGen * 1e6)
plt.xlabel ("Years Before Present")
plt.ylim (1e3, 1e5)
plt.ylabel ("Effective Population Size")
plt.legend (pops, loc='lower left')
plt.title ("Population size history inferred from Chr 1 (1000G)")

# save it
plt.savefig (f"1000G_l{realL}.pdf")
