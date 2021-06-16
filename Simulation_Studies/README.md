The code in this folder is used in producing the plots presented in the CHIMP paper.
The general sequence to follow in reproducing these plots is as follows. 

- Simulate data
- Analyze data using CHIMP, MSMC2, and Relate 
- Plot results

For Simulating and analyzing data, we ran our code on a computing cluster that uses Torque to manage jobs/submissions, so our code is written to make use of that system. However it can be adjusted for other systems. We recommend using some sort of computing system that allows the various analyses to be run in parallel since it can otherwise take too long. We also had a conda environment with msprime, numpy, scipy, pandas and a few other packages installed, and had msmc2 installed. We also note that each of the scripts should be run from the directory in which it is located.

### SIMULATING DATA:

In the msprime_simulation_jobs folder, running the python scripts titled "sim_jobs_\*" will simulate the various datasets we used in our studies. These scripts both create the jobs as \*.sh files, and immediately submit them using torque, before deleting the \*.sh files. This can be modified as needed.

Simulating the data will also create a file of ordered pairs that represents the true demography that data was simulated under, and this is titled [demography]_truth.csv
Data will be written to "/data/" folder, where directories will be created for the various demographies, and each will have 16 replicates of data.

### ANALYZING DATA:

First the Chimp.jar file and appropriate binaries for Relate will need to be downloaded and placed in the appropriate directories ("analysis_CHIMP/job_handlers/CHIMP/\*" and "analysis_relate/job_handlers/py_scripts/binaries/\*"). Msmc2 will need to be installed and accessible from command line. 

Then in the analysis_[method]/job_handlers/ directory, each python script will need to be run. The python scripts are names in the format [demography_type]\_[#parameters for inference]\_[#samples analyzed] with some additional modifiers. Each of these scripts writes 16 .sh files each of which will run 1 analysis of 1 dataset. It then submits these through torque, and deletes the job files. The output is written to the analysis_[method]/out directories.


### PLOTTING DATA:

In order to plot data once all the analyses have been run:

First the appropriate [demography]_truth.csv files should be moved to plotting_tools/out/truth_trajectories/\* . These were created when data was simulated, and will be used to show the "true" demography in plots when comparing the results of inference. 

Then the plotting_tools/out/script_cons_[method].sh scripts need to be run. These will create consolidated results files by consolidating the results from across the 16 replicates for each of the individual analyses.

Finally, running plot_script.py will create the plots using the consolidated data, and write them to the Fplots directory. 





