# Analysis of the 1000 Genomes Data

This directory contains scripts to analyze Chromosome 1 of the `LWK`, the `JPT`, and the `FIN` populations from the 1000 Genomes dataset.

### ANALYZING DATA:

First the Chimp.jar file needs to be downloaded and placed in the appropriate directory "analysis_CHIMP/job_handlers/CHIMP/\*"

Then in the analysis_CHIMP/job_handlers/ directory, the  python script `prepare1000Gdata.py` needs to be run. This produces commands to download the vcf files and for preprocessing it into unphased vcfs ready to be analyzed. Then, the python scripts 
- `FIN_chr1_00p_l-5_TH_n2510.py`
- `JPT_chr1_00p_l-5_TH_n2510.py`
- `LWK_chr1_00p_l-5_TH_n2510.py`
need to be executed to analyze the respective data. Each script checks for the availability of the reference fasta for Chromosome 1 and a fasta file containing the ancestral alleles. If these files are not present, it prints instructions where to download them. The scripts then submit these jobs through torque, and deletes the job files. The output is written to the analysis_CHIMP/out directory.


### PLOTTING DATA:

In order to plot the results, the script `plot1000G.py` in the directory `plotting_tools/out/` has to be executed. This produces the output file `1000G_l-5.pdf` displaying the inferred size histories for the threee population groups.





