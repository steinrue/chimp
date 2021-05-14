import subprocess
import pandas
import matplotlib.pyplot as plt
import numpy
import time

# mostly based on information from
# https://github.com/stschiff/msmc2
# https://github.com/stschiff/msmc
# https://github.com/stschiff/msmc/blob/master/guide.md
# https://github.com/stschiff/msmc-tools
# Li & Durbin (2011)


# EXTERNAL TOOLS
# available through conda: conda install -c bioconda bcftools
bcftoolsBinary = "bcftools"

# available through conda: conda install -c bioconda tabix
bgzipBinary = "bgzip"

# can be activated on the CRI cluster using
# module load dmd/2.072.1
# module load msmc2/2.1.2
msmc2Binary = "msmc2"


# EXTERNAL PYTHON SCRIPTS
# converts 0/1 nucleotide encoding from msprime to proper A/T
nucleotideConversionScript = "jobs_handler/py_scripts/convertVcf.py"

# produces output for msmc2
# downloaded from https://github.com/stschiff/msmc-tools
multihetScript = "jobs_handler/py_scripts/generate_multihetsep.py"


def msmcTimeToGener (time, mu):
    # according to msmc manual
    return time/mu


def msmcLambdaToSize (daLambda, mu):
    # according to msmc manual
    return (1/daLambda)/(2*mu)


def prepareInput (input_vcf, msmc2_input, tmp_dir, inds_to_extract):
    s_time = time.clock_gettime(time.CLOCK_REALTIME)

    print ("[CONVERT_VCF_TO_MSMC2]")

    print (input_vcf)
    print (msmc2_input)
    print ("temp files stored at -> " + tmp_dir)
    tmpVcfBasename = f"{tmp_dir}tmp"

    # convert alleles in vcf from 0/1 to A/T
    tmpVcf = f"{tmpVcfBasename}_complete.vcf"
    process = subprocess.Popen (f"python {nucleotideConversionScript} {input_vcf} > {tmpVcf}", stdout=subprocess.PIPE, shell=True)
    # make it so
    process.communicate()

    # get names of individuals in file
    process = subprocess.Popen (f"{bcftoolsBinary} query -l {tmpVcf}", stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    indNames = [x.decode('ASCII') for x in proc_stdout.split()]
    print (indNames)
    #assert (len(indNames) == len(inds_to_extract))
    print (inds_to_extract)

    # prepare a file for each individual as required by msmc2
    # IF ONLY FOR SUBSET, SPECIFY HERE
    tmpFiles = []
    commands = []
    for i in range(len(indNames)):

        # should we include this one?
        if (not inds_to_extract[i]):
            continue

        # get the name
        ind = indNames[i]

        # vcf file for this one
        tmpFiles.append(f"{tmpVcfBasename}_{ind}.vcf")
        
        # extract genotypes of this individual
        cmd = f"{bcftoolsBinary} view -s {ind} {tmpVcf} > {tmpFiles[-1]}"
        # and append the command to list to be executed
        commands.append (cmd)
        
        # bgzip individual vcfs
        # BGZIP IS NEEDED FOR THIS
        cmd = f"{bgzipBinary} -f {tmpFiles[-1]}"
        commands.append (cmd)

    # collect all commands and execute them
    cmd = f"; ".join(commands)
    print (cmd)
    process = subprocess.Popen (cmd, stdout=subprocess.PIPE, shell=True)
    process.communicate()

    # and put all the individual files together
    gzList = [f"{x}.gz" for x in tmpFiles]
    cmd = f"python3 {multihetScript} {' '.join(gzList)} > {msmc2_input}"
    print (cmd)
    process = subprocess.Popen (cmd, stdout=subprocess.PIPE, shell=True)
    process.communicate()

    # clean up
    toDelete = [tmpVcf] + gzList
    cmd = "; ".join([f'rm {x}' for x in toDelete])
    print (cmd)
    process = subprocess.Popen (cmd, stdout=subprocess.PIPE, shell=True)
    process.communicate()

    print ("[DONE]")
    f_time = time.clock_gettime(time.CLOCK_REALTIME)
    print (f"{(f_time - s_time)/60.0} minutes")

def wattersonsEstimator (segSites, sampleSize):
    return segSites/sum(1/numpy.arange(1,sampleSize))


def getMsmcIntervals (nrIntervals, theta, pairs):

    # according to Li & Durbin (2011), the formula is
    # (0.1 * numpy.exp (i/n * numpy.log (1 + 10 * T_MAX)) - 0.1) * theta
    # with
    # T_MAX = 15 ( * 2N_E generations)
    # n = # intervals
    # i = index of current interval
    # theta = population scaled mutation rate
    T_MAX = 15

    # i/n
    # fractions = numpy.linspace(0, 1, nrIntervals+1)[:-1]
    fractions = numpy.linspace(0, 1, nrIntervals+1)

    #leftBounds = (0.1 * numpy.exp (fractions * numpy.log (1 + 10 * T_MAX)) - 0.1) * theta # old
    leftBounds = (0.1/pairs) * (numpy.exp (fractions * numpy.log (1 + pairs * 10 * T_MAX)) - 1.) * theta 

    rightBounds = numpy.concatenate ((leftBounds[1:], numpy.array([float("inf")])))

    return (leftBounds, rightBounds)


def getMsmcTimeSegmentPattern (timesInGener, maxIntervals, theta, mu, pairs):
    
    # try out a bunch of intervals and see who fits best
    bestFit = float("inf")
    bestNrIntervals = 0
    for nrIntervals in range(1,maxIntervals+1):

        (currLeftBounds, currRightBounds) = getMsmcIntervals (nrIntervals, theta, pairs)
        currLeftBounds = msmcTimeToGener (currLeftBounds, mu)

        # what's the fit
        mySum = 0
        for time in timesInGener:
            thisDistance = numpy.min (numpy.absolute(time - currLeftBounds))
            # squaresum for now
            mySum += thisDistance * thisDistance

        # does this one fit better?
        if (mySum < bestFit):
            bestFit = mySum
            bestNrIntervals = nrIntervals

    # show the best fit
    (leftBounds, rightBounds) = getMsmcIntervals (bestNrIntervals, theta, pairs)
    leftBounds = msmcTimeToGener (leftBounds, mu)

    # print (leftBounds)

    # now that we have the best one, prepare the right pattern string
    # first, whom to group
    positions = []
    for time in timesInGener:
        # what are the distances
        distances = numpy.absolute (leftBounds- time)
        # where is the minimum
        positions.append (numpy.argmin(distances))
    # positions = numpy.searchsorted (leftBounds, timesInGener)

    # put beginninf and end
    positions = numpy.concatenate (([0], positions, [bestNrIntervals]))

    # print (positions)
    # remove duplicates (means that we lose precision, but what to do ...)
    positions = sorted(list(set(positions)))
    # print (positions)

    # what sizes are the groups?
    stretches = numpy.diff(positions)

    # make an msmc pattern string
    patternString = "1*" + "+1*".join([str(x) for x in stretches])
    return patternString


def combPair(num_haps):
    out = "0-1"
    for i in range(num_haps):
        for j in range(i+1,num_haps):
            if(i==0 and j==1):
                out = "0-1"
            else:
                out += f",{i}-{j}"
    return out

def indPair(num_haps):
    out = ""
    for i in range(num_haps-1):
        if i%2 == 0:
            out += f"{i}-{i+1}"
        else:
            if(i != num_haps-2):
                out += ","
    return out


def prepareAuxiliaryFiles():
    pass


def runMsmc2 (msmc2_input, msmc2_output_base, nrIterations, inds_to_analyze, nrIntervals, changePoints, maxHMMstates, mu, theta, extra_options):
    s_time = time.clock_gettime(time.CLOCK_REALTIME)
    print ("[RUN_MSMC2]")

    ind_pairs = inds_to_analyze.count(",") + 1.0
    
    # if no nr given, we need changePoints
    if (nrIntervals == 0):
        assert (len(changePoints) > 0)
        assert (maxHMMstates > 0)
    # if -1, we use default of msmc and provide no string
    if (nrIntervals == -1):
        assert (maxHMMstates == 0)
    # if changepoints given, the nr for equidistant should be zero
    if (len(changePoints) > 3):
        assert (nrIntervals < 1)


    # now get the appropriate time segment string
    timeSegmentString = ""
    if (nrIntervals == 0):
        # get some time segement pattern string to run msmsc2 from reverse engineering changepoints
        timeSegmentString = " -p "+ getMsmcTimeSegmentPattern (changePoints, maxHMMstates, theta, mu, ind_pairs)
    elif (nrIntervals == -1):
        # default, no option provided here (can provide an explicit string using extra_options param
        pass
    else:
        # just equidistant in log-space according to msmc2 default
        timeSegmentString = " -p "+ f"{nrIntervals}*1"
   
    print ("time segment option : " + timeSegmentString)

    
    # msmc2 -i 16 -I 0-1 -p 2*4 -o pwcsaw_results pwcsaw.msmc
    msmc2Cmd = msmc2Binary + \
        " -i " + str(nrIterations) + \
        " -I " + inds_to_analyze + \
        timeSegmentString + \
        " " + extra_options + " " \
        " -o " + msmc2_output_base + " " + msmc2_input 

    # show it
    print (msmc2Cmd)

    # run it
    process = subprocess.Popen (msmc2Cmd, stdout=subprocess.PIPE, shell=True)
    print (process.communicate())

    f_time = time.clock_gettime(time.CLOCK_REALTIME)
    print (f"{(f_time - s_time)/60.0} minutes")
    print ("[DONE]")


def createHistoryFile (msmc2_output, history_file, mu, minGener, maxGener, resolution):

    print ("[CREATE_REAL_OUTPUT]")

    # read in the output    
    df = pandas.read_csv (msmc2_output, sep="\t")

    # put in understandable units
    generRightBoundary = msmcTimeToGener (numpy.array(df["right_time_boundary"]), mu)
    popSize = msmcLambdaToSize (numpy.array(df["lambda"]), mu)

    # make a grid
    times = numpy.exp(numpy.linspace(numpy.log(minGener), numpy.log(maxGener), resolution))

    # now get the sizes
    sizes = numpy.zeros (len(times))

    # where does the grid fall in terms of the right boundaies
    sizeIdx = numpy.searchsorted (generRightBoundary, times)

    # get the corresponding sizes
    sizesOnGrid = popSize[sizeIdx]

    # and write it to file
    ofs = open (history_file, "w")
    ofs.write ("x, y \n")
    for i in range(len(times)):
        ofs.write (f"{times[i]},{sizesOnGrid[i]}\n")
    ofs.close()

    print ("[DONE]")


