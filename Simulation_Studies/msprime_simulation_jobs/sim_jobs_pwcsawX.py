import subprocess
import os
import sys
from pathlib import Path
import msprime
import time
import pandas as pd
import numpy as np

# RUN WITHOUT ARGUMENTS THIS FILE PRINTS OUT SH FILES FOR THE JOBS



if len(sys.argv) == 1:
    home_path = str(Path.home())
    
    #######################################3

    #datadir = "../sim_jobs"
    datadir = f"{home_path}/labshare/gupadhya_folder/sims_final/pwcsawX"
    dataset = "pwcsawX_dataset"
    dem_type="pwcsawX" # "exp" or "bottle" or "bottleExp" or "pwcsaw" or "saw"

    replicates = 16
    n_samples = 200
    rec_rate=.0000000125
    mut_rate=.0000000125
    genome_length= 200000000
    ts_switch = False  # this is switch for also outputing _TS.csv tree sequence file

    np.random.seed(111)

    #######################################3


    
    for i in range(replicates):
        datafile = dataset + str(i+1)
        

        i_file = open("job"+str(i+1)+".sh","w+" )

        i_file.write( "#!/bin/bash" + "\n" )
        i_file.write( f"#PBS -N {dem_type}_sim_{i+1} \n" )
        i_file.write( "#PBS -S /bin/bash" + "\n" )
        i_file.write( "#PBS -l walltime=24:00:00" + "\n" )
        i_file.write( "#PBS -l nodes=1:ppn=1" + "\n" )
        i_file.write( "#PBS -l mem=16gb" + "\n" )
        i_file.write( f"#PBS -o {datadir}/{datafile}.out" + "\n" )
        i_file.write( f"#PBS -e {datadir}/{datafile}.err" + "\n" )
        i_file.write( "cd $PBS_O_WORKDIR" + "\n" )

        # setup appropriate environ
        i_file.write( "source activate PY_deminf \n" )

        
        # now write the line to run java jar with all arguments
        ###################################

        ####################################3
        r_seed= np.random.randint(10000000) # based on seed set at top of this program
        sim_command = f"simulateVCFs('{datadir}','{datafile}',1,{n_samples},{rec_rate},{mut_rate},{genome_length},'{dem_type}',{r_seed},{ts_switch})"
        sim_command = f"python3 -c \"from simVCF_demography import simulateVCFs; {sim_command} \" \n "

        #sim_command = sim_command + f"rm {datadir}/{datafile}.out \n"
        #sim_command = sim_command + f"rm {datadir}/{datafile}.err \n"

        if(i==1):
            sim_command = sim_command + f"mv {datadir}/{datafile}_truth.csv {datadir}/{dataset}_truth.csv \n"
        else:
            sim_command = sim_command + f"rm {datadir}/{datafile}_truth.csv \n"

        
        i_file.write( sim_command + "\n")
        
        i_file.write("\n")
        i_file.close()

        os.system("qsub "+ "job"+str(i+1)+".sh")
        os.system("rm "  + "job"+str(i+1)+".sh" )
        


  

