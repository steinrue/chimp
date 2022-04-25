import subprocess
import pandas
# import matplotlib.pyplot as plt
import numpy
import os
numpy.random.seed (4711)
import time

# mostly based on information from
# https://myersgroup.github.io/relate/


# EXTERNAL TOOLS
# relate uses 'Rscript' for plotting
# needs libraries ggplot2, gridextra
# can be obtained with conda using
# conda install -c r r
# conda install -c r r-ggplot2
# conda install -c r r-gridextra


# BINARIES
# can be downloaded at # https://myersgroup.github.io/relate/
#macVlin = "relate_v1.1.2_MacOSX" #mac
lin_bin2 = "relate_v1.1.2_x86_64_static" #linux
lin_bin3 = "relate_v1.1.3_x86_64_static"
relateBinary = "./jobs_handler/py_scripts/binaries/"+lin_bin3+"/bin/Relate"
# There was some bug in Relate in bin3, that's why we have this weird versioning
relateFileFormatBinary = "./jobs_handler/py_scripts/binaries/"+lin_bin2+"/bin/RelateFileFormats"
relatePopSizeScript = "./jobs_handler/py_scripts/binaries/"+lin_bin3+"/scripts/EstimatePopulationSize/EstimatePopulationSize.sh"


def prepareInput (input_vcf_name, input_filtered_vcf, input_haps_file, input_sample_file, num_inds, ind_labs=None):
    s_time = time.clock_gettime(time.CLOCK_REALTIME)
    print ("[CONVERT_VCF_TO_RELATE]")

    vcfFilter = f"vcftools --vcf {input_vcf_name}.vcf --recode --out {input_filtered_vcf} "

    if ind_labs==None:
        for ind in range(num_inds):
            vcfFilter = vcfFilter + f"--indv tsk_{ind} "
    else:
        for ind in ind_labs:
            vcfFilter = vcfFilter + f"--indv {ind} "


    vcfFilter = vcfFilter + "\n"
    vcfFilter = vcfFilter + f"mv {input_filtered_vcf}.recode.vcf {input_filtered_vcf} \n"

    # show it and run
    print(vcfFilter)
    process = subprocess.Popen (vcfFilter, stderr=subprocess.PIPE, shell=True)
    print (process.communicate())
    

    # convert vcf to relate input
    # get basename from vcf, because relate doesn't want extensions
    base_vcf_name = os.path.splitext(input_filtered_vcf)[0]

    

    # binaries/relate_v1.1.2_x86_64_static/bin/RelateFileFormats --mode ConvertFromVcf --haps exp_dataset1.haps --sample exp_dataset1.sample --input ../exp_dataset1
    vcfToRelateCmd = relateFileFormatBinary + \
        " --mode ConvertFromVcf " + \
        " --haps " + input_haps_file + \
        " --sample " + input_sample_file + \
        " --input " + base_vcf_name

    # show it
    print (vcfToRelateCmd)

    # run the command
    process = subprocess.Popen (vcfToRelateCmd, stderr=subprocess.PIPE, shell=True)
    print (process.communicate())

    # relate is algergic to SNPs being listed more than once, so we have to remove duplicates
    # reading everything at once is maybe not the most efficient, but ok for now
    relateInput = pandas.read_csv (input_haps_file, header=None, sep=" ")
    # the second row has the SNP numbers
    # we drop all but the first occurence
    relateInput.drop_duplicates (subset=[2], keep='first', inplace=True)
    # and write with duplicates removed
    ofs = open (input_haps_file, "w")
    relateInput.to_csv (ofs, header=None, sep=" ", index=False)
    ofs.close ()

    f_time = time.clock_gettime(time.CLOCK_REALTIME)
    print (f"{(f_time - s_time)/60.0} minutes")
    print ("[DONE]")


def prepareAuxiliaryFiles (recoProb, numLoci, input_genetic_map_file, input_sample_file, input_poplabels_file):
    
    # we need a genetic map
    # so make a file for a unifom genetic map

    print ("[FAKE_GENETIC_MAP]")

    # (pos, cM/MB, cM)
    # (p, r, rdist)
    # r[i] = (rdist[i+1] - rdist[i])/(p[i+1] - p[i]) * 1e6
    ofs = open (input_genetic_map_file, "w")
    # recoProb * 100 because it is CENTI morgan
    # not sure whether the formating might cause problems for certain parameters
    ofs.write (f"0\t{recoProb*100*1e6}\t0\n")
    ofs.write (f"{int(numLoci):d}\t{recoProb*100*1e6}\t{numLoci*recoProb*100}\n")
    ofs.close ()

    print ("[DONE_FAKING_GENETIC_MAP]")

    print ("[PRODUCE_FILE_WITH_POPLABELS]")

    # first get the labels of the individuals in the vcf
    sampleFile = pandas.read_csv (input_sample_file, header=None, sep="\t")
    # first two rows are bogus, labels are in first column
    indLabels = sampleFile.iloc[2:,0]

    # now write the poplabels file
    ofs = open (input_poplabels_file, "w")
    # header
    ofs.write ("sample population group sex\n")
    # write a line for each sample
    for ind in indLabels:
        ofs.write (f"{ind} POP1 GRP1 NA\n")
    ofs.close ()

    print ("[DONE]")


def runRelate (input_haps_file, input_sample_file, input_genetic_map_file, input_poplabels_file, mu, effectiveSize, popSizeIterations, yearsPerGen, binParameters, treeOutputPrefix, tree_log_file, sizeOutputPrefix, size_log_file):
    s_time = time.clock_gettime(time.CLOCK_REALTIME)
    print ("[RUN_RELATE_TREES]")

    # first estimate local trees using relate
    # ./binaries/relate_v1.1.2_x86_64_static/bin/Relate --mode All -m 1.25e-8 -N 20000 --haps exp_dataset1.haps --sample exp_dataset1.sample --map exp_dataset1.map --seed 4711 --output result_exp_dataset1
    # --memory can be used to control memory usage
    effpop = ""
    if (effectiveSize>0):
        effpop = " -N " + f"{int(effectiveSize):d}"

    relateTreeCmd = relateBinary + \
        " --mode All " + \
        " -m " + f"{mu:.4e}" + \
        effpop + \
        " --haps " + input_haps_file + \
        " --sample " + input_sample_file + \
        " --map " + input_genetic_map_file + \
        " --seed " + f"{numpy.random.randint (0,9999999)}" + \
        " --memory " + "10" + \
        " --output " + treeOutputPrefix + \
        f" 2>&1 | tee -a {tree_log_file}"

    # show it
    print (relateTreeCmd)

    # run the command
    process = subprocess.Popen (relateTreeCmd, stderr=subprocess.PIPE, shell=True)
    process.communicate()

    # in addition to the log-file, this produces the files treeOutputPrefix.mut/.anc which are the input for size inference
    f_time = time.clock_gettime(time.CLOCK_REALTIME)
    print (f"{(f_time - s_time)/60.0} minutes")
    print ("[DONE]")

    print ("[RUN_RELATE_POPSIZE]")
    s_time = time.clock_gettime(time.CLOCK_REALTIME)

    # and then estimate the population size
    # ./binaries/relate_v1.1.2_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh --input result_exp_dataset1 -m 1.25e-8 --poplabels exp_dataset1.poplabels --seed 4711 --output result_size_exp_dataset1 --num_iter 2
    # --threads (could be used for parallelisation)
    # --bins
    # --years_per_gen
    bin_string = ""
    if binParameters[2] > 0:
        bin_string = " --bins " + ",".join([str(x) for x in binParameters])

    relateSizeCmd = relatePopSizeScript + \
        " --input " + treeOutputPrefix + \
        " -m " + f"{mu:.4e}" + \
        " --poplabels " + input_poplabels_file + \
        " --seed " + f"{numpy.random.randint (0,9999999)}" + \
        " --num_iter " + f"{popSizeIterations:d}" + \
        " --years_per_gen " + f"{yearsPerGen:d}" + \
        bin_string + \
        " --output " + sizeOutputPrefix + \
        f" 2>&1 | tee -a {size_log_file}"

    # show it
    print (relateSizeCmd)

    # run the command
    process = subprocess.Popen (relateSizeCmd, stderr=subprocess.PIPE, shell=True)
    process.communicate()


    f_time = time.clock_gettime(time.CLOCK_REALTIME)
    print (f"{(f_time - s_time)/60.0} minutes")
    print ("[DONE]")


def createHistoryFile (output_coalescent_rate_file, output_history_file, minGen, maxGen, resolution):

    print ("[CREATE_REAL_OUTPUT]")

    # read in the output    
    ifs = open(output_coalescent_rate_file)

    # skip first line
    ifs.readline()
    # second line has generation times
    generTimes = ifs.readline()
    # third line has coalescent rates
    coalRates = ifs.readline()
    # should be everything
    ifs.close ()

    # now get the numbers from the lines
    generTimes = numpy.array([float(x) for x in generTimes.split()])
    # skip the first two, cause they just indices
    coalRates = numpy.array([float(x) for x in coalRates.split()][2:])

    # convert the rates to (diploid) sizes
    # 0.5/rate according to the documentation (seems to work)
    popSizes = 0.5/coalRates

    # make a grid
    times = numpy.exp(numpy.linspace(numpy.log(minGen), numpy.log(maxGen), resolution))

    # now get the sizes
    sizes = numpy.zeros (len(times))

    # where does the grid fall in terms of the right boundaries
    # use the right boundaries
    sizeIdx = numpy.searchsorted (generTimes[1:], times)

    # get the corresponding sizes
    sizesOnGrid = popSizes[sizeIdx]

    # and write it to file
    ofs = open (output_history_file, "w")
    ofs.write ("x, y\n")
    for i in range(len(times)):
        ofs.write (f"{times[i]},{sizesOnGrid[i]}\n")
    ofs.close()

    print ("[DONE]")

