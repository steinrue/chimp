# CHMM Inference

Here we provide an implementation of CHIMP in `java`, a method that infers the population size history of a population from whole genome sequencing data. Background on the method and its theoretical models can be found in this [manuscript](https://doi.org/10.1101/2021.05.22.445274).

The method can be run with the `Chimp.jar` file provided for download in this repository (the relevant libraries are packaged into this jar-file). 

The subdirectory `Simulation_Studies` contains `python` scripts that were used to conduct the simulation studies in the [manuscript](https://doi.org/10.1101/2021.05.22.445274).


# Input and Output

In order to run CHIMP, you need to make sure you have three things for your data:
-  .vcf files containing data for SNPs. Non-bialellic SNPs (e.g. tri-allellic, structural variation, ...) are filtered out by Chimp.
-  .fasta files for the reference sequence. This specifies the alleles at the sites that are not listed in the .vcf ('N' for missing sites).
-  .fasta files for the ancestral sequence. If you do not have this, the reference files can be used here, though the inference may be affected.

CHIMP will analyze the VCFs provided, selecting a specified number of haplotypes from left to right (it will treat diploid individual data as two separate haplotypes) and ignoring additional haplotypes if the file contains more than specified. If this method for selecting the haplotypes might cause bias for your data, we recommend reordering (e.g. randomly permuting) the haplotypes before running CHIMP. 

Any variant listed in the VCF that is not a biallelic SNP for the samples to be analyzed (e.g. tri-allellic, structural variation, ...) will be treated as a missing site, and positions not listed in the VCF are treated as non-segregating. The reference FASTA file is used to specify tracts of missing data, and the ancestral FASTA file helps distinguish ancestral from derived alleles. 

CHIMP will infer a piece-wise constant population size history with 20 epochs. The boundaries between these epochs are distributed exponentially (uniformly on a log scale) between Ne/50 and 5*Ne. Ne is the effective population size computed from Waterson's Estimator across the data, and is also the value at which the population size in each epoch is initialized.

Upon completion, CHIMP will output two files. The PARAM file contains the results of the inference, listed as a single ordered pair for each epoch of the history. Each ordered pair contains [epoch's lower bound in generations, population size during epoch]. All times are given in generations before present, and population sizes are given in number of diploid individuals. The CSV file contains a dense set of ordered pairs (with many pairs for each epoch) that can be used for plotting. 

# Necessary Parameters and a Minimal Example

When running CHIMP, the following options are necessary.By default, CHIMP will use the TMRCA as the hidden state.
```
--rec_rate <recombination_rate>
        Recombination rate, specified per generation per nucleotide (e.g.
        .0000000125 for human data).

--mut_rate <mutation_rate>
        Mutation rate along a lineage, specified per generation per nucleotide
        (e.g. .0000000125 for human data).
        
--base_n <base_n_samples>
        Number of haplotypes considered in each parallel CHMM model (The total
        number of haplotypes analyzed is <base_n> * <n_groups> under our
        composite likelihood scheme, default for <n_groups> is 1). We
        recommend at most choosing --base_n=10 (though for TMRCA the method is
        still tractable up to 30) and adjusting <n_groups> accordingly to tackle
        larger sample sizes.
        
--vcf_list [vcf_1,vcf_2,...,vcf_N]
	List of VCF files for each chromosome to perform inference on. Will
        analyze the first (<base_n> * <n_groups>, default for <n_groups> is
        one) haplotypes in each file.
        
--ref_list [ref_1,ref_2,...,ref_N] 
        List of reference files for each chromosome. Indexing and length of list
        should match [vcf_list] and [anc_list].
        
--anc_list [anc_1,anc_2,...,anc_N]
        List of ancestral files for each chromosome. Can specify same files as
        [ref_list] if you don't have these (though accuracy may be affected if
        reference alleles and ancestral alleles are not the same).

--out_file <output_file>
        Prefix for output files CHIMP will produce with the inference results.
        CSV contains ordered pairs of [time, population size], and PARAM file
        contains, for each epoch, [epoch's lower bound (time), inferred
        population size]. Times are given in generations before present, and
        population sizes are in number of diploid individuals.
```

Thus, a minimal command line example is:
```
java -jar Chimp.jar --vcf_list=chrom1.vcf --ref_list=ref1.fasta
      --anc_list=anc1.fasta --out_file=results --rec_rate=.0000000125
      --mut_rate=.0000000125 --base_n=10
```
This will output the files results.csv and results.param which contain the inferred population size histories. The inference is performed using the default parameters for the model. 


# Advanced Options for Customized Performance

Further customization of the model and the inference procedure can be achieved with the following options. 


```
--n_groups <sub_sample_groups>
        Number of subordinate CHMMs (each considering <base_n> haplotypes)
        analyzed in composite analysis. Total number of haplotypes analyzed is 
        <base_n> * <n_groups>. (default: 1)
        
--t_bounds [min_time,max_time]
        This is used to specify epochs for population size history. The epoch
        boundaries will be spaced exponentially between <min_time> and
        <max_time>. Defaults to [Ne/50, 5*Ne] where Ne is computed from
        Waterson's estimator.

--dof <degrees_of_freedom>
        Number of epochs of population size history between [t_bounds] (with 2
        additional epochs added, one above and one below [t_bounds]). One
        parameter is inferred for each epoch during EM. Note that dof=X will
        yield X+2 epochs. (default: 18)
        
--psh0 <initial_size_guess>
        Initial guess for population size. Default is to use Waterson's estimator
        for N_effective, computed from data (and based on the specified mutation
        rate). Specifying this option will also use <psh0> instead of Waterson's
        estimate in computing default for [t_bounds] and the partitioning CHMM
        states.

--psh0_xs [epoch_bound_1,epoch_bound_2,...,epoch_bound_N]
        Specify custom list of epoch boundaries (overrides [t_bounds], and
        <dof>) to use for population size history model.
        
--psh0_ys [population_epoch1,population_epoch2,...,population_epochN+1]
        Specify population sizes (must also have specified [psh0_xs]) to
        initialize EM. This option overrides psh0. [psh0_ys] must be a vector
        with size 1 greater than that of [psh0_xs].
        
--n_states <num_chmm_states>
        Number of discrete states for CHMM (discrete intervals into which TMRCA
        or L falls). States are partitioned according to equal probability for
        each state under a constant population size prior (Waterson's estimate).
        (default: 50)

--metalocus_size <metalocus_size>
        Number of bases grouped into each metalocus to speed up computations of
        the E-Step. Selecting metalocus_size=1 will activate a locus skipping
        algorithm to aid computational efficiency. (default: 500)

--em_cap <max_em_steps>
        Maximum number of EM steps performed. Inference will end sooner if the
        convergence criterion is met. (default: 50)
        
--m_evals <m_step_cap>
        Maximum number of transition/emission matrix evaluations allowed during
        maximization step. Default scales with number of parameters (2 * <dof> +
        10) since higher dimensional spaces will need more function evaluations
        to optimize effectively. M-step will terminate sooner if optimization
        converges.
        
--ll_converge <em_convergence_criterion>
        Absolute threshold used to determine convergence of EM. If the posterior
        log-likelihood improves by less than this amount after any EM step, we
        consider the method converged. (default: .02)

-l, --tree_length
        Use this to specify use of total branch length (L) instead of tree
        height (TMRCA) as representation of the CHMMs hidden state.
        
--output_steps
	      This option will produce new output files after each EM step, titled
        "EMstep_X" where X is the EM step number. The files will be placed in a
        new folder entitled <out_file>.
        
--pseudo
	      Use this switch to specify that the data is pseudohaploid data. CHIMP 
	      will assume each allele for each variant in the VCF is randomly selected 
	      from either haploid of a diploid individual (independently for each 
	      variant/individual). Emission probabilities will be adjusted accordingly.
	
```

We now show some usage examples using these options. 

 ```
java -jar Chimp.jar --vcf_list=chrom1.vcf,chrom2.vcf --ref_list=ref1.fasta,ref2.fasta
        --anc_list=anc1.fasta,anc2.fasta --out_file=results --rec_rate=.0000000125
        --mut_rate=.0000000125 --base_n=10 --n_groups=15 --tree_length
```

The above example analyzes two chromosomes together (chrom1, chrom2), analyzes 150 haplotypes (in 15 groups of 10), and uses total branch length (L) instead of TMRCA as the hidden state.

 ```
java -jar Chimp.jar --vcf_list=chrom1.vcf --ref_list=ref1.fasta --anc_list=anc1.fasta
        --out_file=results --rec_rate=.0000000125 --mut_rate=.0000000125 --base_n=10
        --n_groups=20 --t_bounds=500,4000 --dof=3
```

The above example analyzes a single chromosome using 200 haplotypes and the TMRCA model for the CHMM. The population size history has 5 epochs partitioned by the times [500,1000,2000,4000]. Note that this is achieved by specifying 3 degrees of freedom (that will be exponentially spaced) between the bounds 500 and 4000. Alternatively we could have specified `--psh0_xs=500,1000,2000,4000` instead of `t_bounds` and `dof` to achieve the same result. 

 ```
java -jar Chimp.jar --vcf_list=chrom1.vcf --ref_list=ref1.fasta --anc_list=anc1.fasta
        --out_file=results --rec_rate=.0000000125 --mut_rate=.0000000125 --base_n=10
        --psh0_xs=1000,2000,3000 --psh0_ys=40000,20000,10000,5000 --ouput_steps
```
In the above example (for 10 haplotypes) we have not only customized the epoch partitions (to be [1000,2000,3000] which would have been unachievable with `t_bounds` and `dof` due to the non-exponential sequence), but we have also specified a prior guess of population sizes (40k,20k,10k,5k) in the epochs. This population history is where the EM will be initialized, instead of assuming a size of N_e in each epoch (computed from Waterson's estimate). The final option will also cause the method to output the intermediate population size histories as well after each EM step. This can be useful for analyzing convergence behavior or overfitting behavior of the method.


 ```
java -jar Chimp.jar --vcf_list=chrom1.vcf --ref_list=ref1.fasta --anc_list=anc1.fasta
        --out_file=results --rec_rate=.0000000125 --mut_rate=.0000000125 --base_n=10
        --n_states=25 --metalocus_size=1000 --em_cap=20 --m_evals=35 --ll_converge=.5
```

In the above example, we have shown the usage of several options that will make CHIMP complete the inference faster, though with less accuracy. By reducing the number of `n_states` we move further away from the SMC (the "true" model), but gain speed for the E-step. Similarly by increasing `metalocus_size` we are assuming that larger segments of the genome are described by the same tree, increasing E-step speed but reducing accuracy of the model. We have also reduced the number of EM steps, the probability matrix computations per M-step, and widened the convergence threshold -- all of which will cause CHIMP to finish quicker, but possibly before it can near the true likelihood maximum. 


# Additional Options (Most users should not need to change any of these options)

The following options may be used.  

```
--partitions [probability_1,probability_2,...,probability_N]
        List of marginal probabilities to partition states of the CHMM. Under a
        constant population prior (from Waterson's estimate) the marginal
        probability of the TMRCA or L falling into each state will be given by
        [partitions]. The probabilities are indexed so that the first
        probability corresponds to the state with the smallest TMRCA (or L) and
        the last corresponds to the largest. [partitions] should sum to 1. This
        option overrides n_states.

--ps_scale <population_rescale_factor>
        Internal Population Rescale Factor. Modify this if ODE stepsize error
        thrown. (default: 1000)
        
--max_tract <max_tract_size>
        Maximum Tract Length of monomorphic tracts during locus
        skipping-algorithm. The locus-skipping algorithm is only used if
        <metalocus_size>=1. Reduce <max_tract> to avoid underflow errors during
        the E step. (default: 1000)

--pde_res [emission_pde_resolution,transition_pde_resolution]
        Density of PDE solver's grid over the range of relevant times/lengths
        when using tree-length as hidden state. Higher density increases
        accuracy but also run-time. First parameter is for the emission
        probability solver (2D), and second is for transition probability solver
        (3D). This field is ignored for TMRCA. (default: 1000,250)

--simplex_scale <simplex_init_scale>
        Size (edge length) of simplex used in M step. Since the search is in log
        space, 1.0 corresponds to an e-fold increase/decrease of the population
        size parameter. (default: .1)

--fix_ss
        Use this switch to use same size of simplex for all EM steps (default is
        to slightly shrink the simplex on each successive step).
        
--disc_storage <temporary_file_prefix>
        Use option to store intermediate data structures temporarily on the disc
        (default stores them in RAM). Using this option may decrease RAM usage,
        but may increase run-time since structures will be stored to and read
        from disc.
       
--chr_l [chromosomal_length1,chromosomal_length2,...,chromosomal_lengthN] 
        List of chromosome lengths for each VCF. Will analyze positions
        [1,length] for each file. Indices of [chr_l] should match those of
        [vcf_list]. If not provided, the reference file will determine the
        contiguous chromosomal tract that is analyzed.
```
      
The following options have not been extensively tested and we recommend caution when using them.

```
--reg_lambdas [regularization1,regularization2,regularization3,regularization4]
        Regularizing coefficients used in likelihood optimization. Coefficients
        of [d0l2, d1l1, d1l2, d2l2] respectively. These coefficients are
        multiplied by [d0l2,d1l1,d1l2,d2l2] respectively and are subtracted from
        the objective function. d0l2 corresponds to the squared deviation from a
        constant population (Waterson?s estimate). d1l1 corresponds to the
        absolute value difference of the population size change between adjacent
        epochs, while d1l2 corresponds to the squared difference of the same.
        d2l2 is an analog for the square of the second derivative (computed for
        discrete epochs instead), and upweighting it will penalize deviations
        from linear behavior. The default is no regularization. (default:
        0,0,0,0)
        
--spline
        This option is not fully tested. It specifies that the population size
        history will be modeled by a cubic spline function (instead of
        piece-wise constant). The PARAM file will not be generated.  This can be
        used in conjunction with [t_bounds] and <dof>, which will specify the
        bounds of the spline and the number of nodes respectively. This should
        not be used with [psh0_xs] or [psh0_ys].

--lin_t
        This will distribute the epoch partitions uniformly on a linear scale
        rather than on a log scale.
```

In implementing this project, we have used the following publicly available libraries:
- [HTSJDK - VCF tool] (https://github.com/samtools/htsjdk)
- [JSAP] (http://www.martiansoftware.com/jsap/)
- [Apache Commons] (https://commons.apache.org/)
