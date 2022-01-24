import subprocess
import os
import sys
from pathlib import Path
import numpy as np


if len(sys.argv) == 1:

    home_path = str(Path.home())

    ####
    # THIS REQUIRES PWCSAWX DATA TO ALREADY HAVE BEEN SIMULATED.
    # THIS SCRIPT WILL REPROCESS THE PWCSAWX DATA AND CREATE PSEUDOHAPLOID DATA FROM IT
    
    input = f'../data/pwcsawX/pwcsawX_dataset'
    output = f'../data/pwcsawX_pseudo'
    os.system(f"mkdir {output}")
    output = output + '/pwcsawX_pseudo_dataset'
    replicates = 16

    '''
    for i in range(replicates):
        command = f"python pseudo_hap_converter.py {input}{i+1}.vcf {ouput}{i+1}.vcf"
        print(command)
        os.system(command)
    '''

    for i in range(replicates):
        
        i_file = open("job"+str(i+1)+".sh","w+" )

        i_file.write( "#!/bin/bash" + "\n" )
        i_file.write( f"#PBS -N sim_{i+1} \n" )
        i_file.write( "#PBS -S /bin/bash" + "\n" )
        i_file.write( "#PBS -l walltime=00:59:00" + "\n" )
        i_file.write( "#PBS -l nodes=1:ppn=1" + "\n" )
        i_file.write( "#PBS -l mem=16gb" + "\n" )
        i_file.write( f"#PBS -o {output}{i+1}.out" + "\n" )
        i_file.write( f"#PBS -e {output}{i+1}.err" + "\n" )
        i_file.write( "cd $PBS_O_WORKDIR" + "\n" )

        # setup appropriate environ
        i_file.write( "source activate PY_deminf \n" )

        
        # now write the line to run java jar with all arguments
        ###################################

        command = f"python pseudo_hap_converter.py {input}{i+1}.vcf {output}{i+1}.vcf \n"
        command = command + f"rm {output}{i+1}.err \n"

        
        i_file.write( command + "\n")
        i_file.write("\n")
        i_file.close()

        os.system("qsub "+ "job"+str(i+1)+".sh")
        os.system("rm "  + "job"+str(i+1)+".sh" )
        

