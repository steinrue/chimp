import subprocess
import pandas
# import matplotlib.pyplot as plt
import numpy
import os
import sys
from pathlib import Path

# RUN WITHOUT ARGUMENTS THIS FILE PRINTS OUT SH FILES FOR THE JOBS
if len(sys.argv) == 1:
    home_path = str(Path.home())
    
    datadir = f"{home_path}/labshare/gupadhya_folder/sims_final/exp/"
    prefix = "3p_10n_TL_"
    dataset = "exp_dataset"



    # make sure you create output directory to place all outfiles
    prefix_d =  "out/" + prefix + dataset
    process = subprocess.Popen("mkdir "+ prefix_d, stderr=subprocess.PIPE, shell=True)
    process.communicate()

    for i in range(16):
        datafile = dataset + str(i+1)
        outDir =  prefix_d + "/"  + datafile

        # make sure you create output directory to place all outfiles
        process = subprocess.Popen("mkdir "+ outDir, stderr=subprocess.PIPE, shell=True)
        process.communicate()


        i_file = open("job"+str(i+1)+".sh","w+" )

        i_file.write( "#!/bin/bash" + "\n" )
        i_file.write( "#PBS -N CHIMP_job_"+str(i+1) + "\n" )
        i_file.write( "#PBS -S /bin/bash" + "\n" )
        i_file.write( "#PBS -l walltime=60:00:00" + "\n" )
        i_file.write( "#PBS -l nodes=1:ppn=1" + "\n" )
        i_file.write( "#PBS -l mem=16gb" + "\n" )
        i_file.write( f"#PBS -o ../{outDir}/{datafile}.out" + "\n" )
        i_file.write( f"#PBS -e ../{outDir}/{datafile}.err" + "\n" )

        i_file.write( "cd $PBS_O_WORKDIR/.." + "\n" )

      
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


    chimp_jar = "jobs_handler/CHIMP/Chimp.jar"
    CHIMP_command = f"java -Xmx16G -jar {chimp_jar} "

    
    ########################
    # CHIMP PARAMETERS
    #########################
    
    ## demography specifications/basic model parameters
    CHIMP_command = CHIMP_command + " --rec_rate=.0000000125 "
    CHIMP_command = CHIMP_command + " --mut_rate=.0000000125 "
    CHIMP_command = CHIMP_command + " --base_n=10 "
    CHIMP_command = CHIMP_command + " --n_groups=1 "

    CHIMP_command = CHIMP_command + " --t_bounds='1000,2000' "
    CHIMP_command = CHIMP_command + " --dof=1 "
    #CHIMP_command = CHIMP_command + " --spline "

    
    ## tree length options
    CHIMP_command = CHIMP_command + " --tree_length "
    CHIMP_command = CHIMP_command + " --pde_res='1000,250' "

    ## EM options
    #CHIMP_command = CHIMP_command + " --em_cap=50 "
    #CHIMP_command = CHIMP_command + " --s_type=0 "
    #CHIMP_command = CHIMP_command + " --simplex_scale=.1 "
    #CHIMP_command = CHIMP_command + " --ll_converge=.02 "
    #CHIMP_command = CHIMP_command + " --m_evals=35 "


    #CHIMP_command = CHIMP_command + " --metalocus=500 "
    #CHIMP_command = CHIMP_command + " --reg_lambdas='0,0,0,0' "
    #CHIMP_command = CHIMP_command + " --ps_scale=1000 "


    ## deal with files
    CHIMP_command = CHIMP_command + f" --data_dir={data_folder} "
    CHIMP_command = CHIMP_command + f" --vcf_list={input_vcf}.vcf "
    CHIMP_command = CHIMP_command + " --ref_list=../sim_ref.fasta "
    CHIMP_command = CHIMP_command + " --anc_list=../sim_anc.fasta "
    CHIMP_command = CHIMP_command + " --chr_l=200000000 "
    CHIMP_command = CHIMP_command + f" --out_file={output_dir}{popsize_out_base} "

    #CHIMP_command = CHIMP_command + " --print_all "






    
    
    ##########################
    ## run CHIMP
    #################################
    
    process = subprocess.Popen ( CHIMP_command, stderr=subprocess.PIPE, shell=True)
    print (process.communicate())

