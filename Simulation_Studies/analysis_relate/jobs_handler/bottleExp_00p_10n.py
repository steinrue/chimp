import subprocess
import pandas
# import matplotlib.pyplot as plt
from py_scripts.analysis_routines import *
import numpy
import os
import sys
from pathlib import Path

numpy.random.seed (4711)


if len(sys.argv) == 1:
    home_path = str(Path.home())

    datadir = f"{home_path}/labshare/gupadhya_folder/sims_final/bottleExp/"
    prefix = "00p_10n_"
    dataset = "bottleExp_dataset"

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
        i_file.write( "#PBS -N RELATE_job_"+str(i+1) + "\n" )
        i_file.write( "#PBS -S /bin/bash" + "\n" )
        i_file.write( "#PBS -l walltime=36:00:00" + "\n" )
        i_file.write( "#PBS -l nodes=1:ppn=1" + "\n" )
        i_file.write( "#PBS -l mem=16gb" + "\n" )
        i_file.write( f"#PBS -o ../{outDir}/{datafile}.out" + "\n" )
        i_file.write( f"#PBS -e ../{outDir}/{datafile}.err" + "\n" )

        i_file.write( "cd $PBS_O_WORKDIR/.." + "\n" )

        # module load all required packages /environs
        i_file.write( "source activate PY_deminf \n" )
        i_file.write( "module load gcc/6.2.0 zlib/1.2.11 perl/5.24.0 samtools/1.10  vcftools/0.1.16\n" )
      
        # now write the line to run java jar with all arguments
        i_file.write( "python " + "jobs_handler/"+sys.argv[0] + " " + datadir + " " + datafile + " " + outDir + "\n")
        
        i_file.write("\n")
        i_file.close()

        os.system("qsub "+ "job"+str(i+1)+".sh")
        os.system("rm "  + "job"+str(i+1)+".sh" )
        

if len(sys.argv) > 1:

    data_folder = sys.argv[1] 
    input_vcf = sys.argv[2]
    output_dir = sys.argv[3]+"/" # output directory where all output files are placed

    popsize_out_base = "deminf_results"



    ########################
    # RELATE PARAMETERS
    #########################

    # popgen parameters
    mu = 1.25e-8
    recoProb = 1.25e-8
    numLoci = 200e6
    num_inds = 5 #number of INDIVIDUALS, or haplotypes / 2
    
    # effective size ( scaling factor )
    effectiveSize = 11000

    # num iterations when estimating the popsize, default is 10
    popSizeIterations = 10

    
	# how many years per generation ( default: 28 )
	# better, use 1 and treat result in gens
    yearsPerGen = 1

    # parameters for PWC demography divisions, not sure of default
	# number are given in years (or generations, if yearsPerGen = 1)
	# format : [lower, upper, stepsize]
	# changepoints for piecewise history: [0, 10^(lower), 10^(lower+stepsize), 10^(lower+2*stepsize), ..., 10^(upper), infty)
    binParameters = [2.30103,4.30103,-1] # To use relate's default, specify -1 for third index ( [lb,ub,-1] ) then it ignores lb and rb and does default everything



    

    ####################################
    # PROCESS VCF TO PREPARE INPUT FILES
    #####################################

    # relate input file
    input_filtered_vcf = f"{output_dir}{input_vcf}_filter.vcf"
    input_haps_file = f"{output_dir}{input_vcf}.haps"
    input_sample_file = f"{output_dir}{input_vcf}.sample"
    input_poplabels_file = f"{output_dir}{input_vcf}.poplabels"
    input_genetic_map_file = f"{output_dir}{input_vcf}.map"

    # prepare the proper input for relate from the vcf
    prepareInput (data_folder + input_vcf, input_filtered_vcf, input_haps_file, input_sample_file, num_inds)

	# other files to prepare
    prepareAuxiliaryFiles (recoProb, numLoci, input_genetic_map_file, input_sample_file, input_poplabels_file)




    ####################################
    # RUN ANALYSIS
    #####################################

    
	# relate output files (or better: prefixes)
    temp_tag = os.path.basename(os.path.split(sys.argv[3])[0])  +"_"+ os.path.basename(sys.argv[3])

    treeOutputPrefix = f"result_tree_{temp_tag}"
    tree_log_file = treeOutputPrefix + ".log"
    sizeOutputPrefix = f"result_size_{temp_tag}"
    size_log_file = sizeOutputPrefix + ".log"

	
	# do the analysis
    runRelate (input_haps_file, input_sample_file, input_genetic_map_file, input_poplabels_file, mu, effectiveSize, popSizeIterations, yearsPerGen, binParameters, treeOutputPrefix, tree_log_file, sizeOutputPrefix, size_log_file)



    
    ##########################################
    # EXTRACT APPROPRIATE HISTORY AND WRITE CSV
    ##########################################
    minGen = 10**binParameters[0] / 2.
    maxGen = 10**binParameters[1] * 2.
    resolution = 2000

    # this is the final output from relate that we can use to create a history
    output_coalescent_rate_file = sizeOutputPrefix + ".coal"
    output_history_file = output_dir + popsize_out_base+ ".csv"

    createHistoryFile (output_coalescent_rate_file, output_history_file, minGen, maxGen, resolution)



    
    ##########################
    ## move all output files to the appropriate directory
    #################################
    
    mv_command = f"mv result_*{temp_tag}.* {output_dir}"
    os.system(mv_command)
    mv_command = f"mv result_*{temp_tag}_* {output_dir}"
    os.system(mv_command)
    rm_command = f"rm {input_haps_file} {input_poplabels_file} {input_sample_file} {input_genetic_map_file} {input_filtered_vcf}"
    os.system(rm_command)
    rm_command = f"rm {output_dir}*.anc {output_dir}*.mut {output_dir}*.dist {output_dir}*.bin"
    os.system(rm_command)
