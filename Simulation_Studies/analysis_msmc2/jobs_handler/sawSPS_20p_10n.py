from py_scripts.analysis_routines import *
import sys
import os
from pathlib import Path

# RUN WITHOUT ARGUMENTS THIS FILE PRINTS OUT SH FILES FOR THE JOBS
if len(sys.argv) == 1:
    home_path = str(Path.home())
    
    datadir = f"../data/sawSPS/"
    prefix = "20p_10n_"
    dataset = "sawSPS_dataset"


    # make sure you create output directory to place all outfiles
    prefix_d =  "out/" + prefix + dataset
    process = subprocess.Popen("mkdir "+ prefix_d, stderr=subprocess.PIPE, shell=True)
    process.communicate()

    for i in range(16):
        datafile = dataset + str(i+1)
        outDir =  prefix_d + "/" + datafile

        # make sure you create output directory to place all outfiles
        process = subprocess.Popen("mkdir "+ outDir, stderr=subprocess.PIPE, shell=True)
        process.communicate()



        
        i_file = open("job"+str(i+1)+".sh","w+" )

        i_file.write( "#!/bin/bash" + "\n" )
        i_file.write( "#PBS -N MSMC2_job_"+str(i+1) + "\n" )
        i_file.write( "#PBS -S /bin/bash" + "\n" )
        i_file.write( "#PBS -l walltime=40:00:00" + "\n" )
        i_file.write( "#PBS -l nodes=1:ppn=1" + "\n" )
        i_file.write( "#PBS -l mem=16gb" + "\n" )
        #i_file.write( "#PBS -l host=!cri16cn051" + "\n" )
        i_file.write( f"#PBS -o ../{outDir}/{datafile}.out" + "\n" )
        i_file.write( f"#PBS -e ../{outDir}/{datafile}.err" + "\n" )

        i_file.write( "cd $PBS_O_WORKDIR/.." + "\n" )

        # module load all required packages /environs
        i_file.write( "module load gcc/6.2.0 python/3.8.1 \n" )
        i_file.write( "source activate PY_deminf \n" )
        i_file.write( "module load dmd/2.072.1 msmc2/2.1.2 \n" )

        
        # now write the line to run java jar with all arguments
        i_file.write( "python " + "jobs_handler/"+sys.argv[0] + " " + datadir + " " + datafile + " " + outDir + "\n")
        
        i_file.write("\n")
        i_file.close()

        os.system("qsub "+ "job"+str(i+1)+".sh")
        os.system("rm "  + "job"+str(i+1)+".sh" )
        


# HERES THE ACTUAL ANALYSIS SEGMENT IF YOU RUN THIS PROGRAM WITH ARGUMENTS
if len(sys.argv) > 1:

    data_folder = sys.argv[1] 
    input_vcf = sys.argv[2]
    output_dir = sys.argv[3] + "/"

    popsize_out_base = "deminf_results"

    
    ########################
    # MSMC PARAMETERS
    #########################

    # EM iterations (default 20)
    nrIterations = 20

    # individuals to consider in preparing msmc input file
    ##haplotypes to analyze in msmc2 analysis
    ## combPair(n) does all comb of first n haps, ##noPair(n) does non overlapping pairs of first n haps
    inds_to_extract_from_vcf = 5*[True] + 95*[False]
    inds_to_analyze = combPair(10) # all pairs among the first 10

    # Loci along genome
    numLoci = 200e6

    # popgen parameters
    mu = 1.25e-8

    ## Specify PWC Demography Breakpoints
    # give a number if equidistant, otherwise 0
    nrIntervals = 0
    # points where the population size changes
    # if nrIntervals = -1 or > 0, then these are just used to find bounds during plotting at the very end
    changePoints = [ 40, 59, 86, 127, 186, 273, 400, 587, 862, 1265, 1857, 2725, 4000, 5871, 8618, 12649, 18566, 27252, 40000]
    # maximal number of hmm states when trying to match the changepoints best
    maxHMMstates = 64

    extra_options = '-r 1.0 --fixedRecombination'


    #########################################
    # PROCESS THE VCF FILE TO MSMC FILE
    ##########################################

    msmc2_input = f"{output_dir}{input_vcf}.msmc" # msmc file name
    nrSamples = 2*sum(inds_to_extract_from_vcf)

    prepareInput (f"{data_folder}{input_vcf}.vcf", msmc2_input, output_dir, inds_to_extract_from_vcf) # process the data



    #########################################
    # RUN MSMC2 ANALYSIS
    ##########################################
    # the theta used in msmc2 comes from some sort of wattersons estimator?
    # a proxy for the number of segreating sites is the number of lines in the msmc2_inout file
    numSegSites = sum(1 for line in open(msmc2_input))
    # for some reason /2
    theta = wattersonsEstimator(numSegSites, nrSamples)/numLoci/2

    msmc2_output_base = f"{output_dir}{input_vcf}"
    runMsmc2 (msmc2_input, msmc2_output_base, nrIterations, inds_to_analyze, nrIntervals, changePoints, maxHMMstates, mu, theta, extra_options)


    ##################################################
    # PLOT FUNCTION FROM OUTPUT
    ##################################################
    final_file = msmc2_output_base + ".final.txt" # output file from having run msmc2
    minGen = changePoints[0] / 2.
    maxGen = changePoints[-1] *2.
    resolution = 2000

    createHistoryFile (final_file, f"{output_dir}{popsize_out_base}.csv", mu, minGen, maxGen, resolution)

    print("DONE")


    #####
    ## clean extra files
    process = subprocess.Popen("rm "+ msmc2_input, stderr=subprocess.PIPE, shell=True)
    process.communicate()
