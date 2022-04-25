Analysis of real data

### ANALYZING DATA:

First the Chimp.jar file and appropriate binaries for Relate will need to be downloaded and placed in the appropriate directories ("analysis_CHIMP/job_handlers/CHIMP/\*" and "analysis_relate/job_handlers/py_scripts/binaries/\*"). Msmc2 will need to be installed and accessible from command line. 

Then in the analysis_[method]/job_handlers/ directory, each python script will need to be run. The python scripts are names in the format [demography_type]\_[#parameters for inference]\_[#samples analyzed] with some additional modifiers. Each of these scripts writes 16 .sh files each of which will run 1 analysis of 1 dataset. It then submits these through torque, and deletes the job files. The output is written to the analysis_[method]/out directories.


### PLOTTING DATA:

In order to plot data once all the analyses have been run:

Then the plotting_tools/out/script_cons_[method].sh scripts need to be run. These will create consolidated results files by consolidating the results from across the 16 replicates for each of the individual analyses.

Finally, running plot_script.py will create the plots using the consolidated data, and write them to the Fplots directory. 





