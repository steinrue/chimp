import msprime
import time
import pandas as pd
import numpy as np
import itertools
import os
    
def main():
    print("python main function")
    out_dir="truth_plots"
    out_file_tag="test_data"
    replicates= 0
    n_samples = 10
    rec_rate=.0000000125
    mut_rate=.0000000125
    genome_length= 200000
    dem_type="bottleExp"
    r_seed= 250
    
    simulateVCFs(out_dir,out_file_tag, replicates, n_samples, rec_rate, mut_rate, genome_length, dem_type, r_seed, False)



def simulateVCFs(out_dir, out_file_tag, replicates, n_samples, rec_rate, mut_rate, genome_length, dem_type, r_seed, treeSeq_switch):


    # start timer to time script
    start = time.time()

    # name of output files
    truth_f = out_dir +"/"+ out_file_tag +"_truth.csv" # csv where truth trajectory is stored
    sim_file_prefix = out_dir+"/"+ out_file_tag # filename for output vcfs

    # seed random number generator
    np.random.seed(r_seed)


    ###
    print(f"random seed is {r_seed}")
    print(f"num samples is {n_samples}")
    print(f"rec rate is {rec_rate}")
    print(f"mut rate is {mut_rate}")
    print(f"{replicates} replicates of {genome_length} bases")
    print(f"demography type: {dem_type}")

    #################
    # define population size history for each of the demography types
    ################

    if(dem_type == "bottle"):
        start_p = [10000,5000,10000] # pop sizes
        growth_rates = [0,0, 0 ] # growth rates
        dts =  [1000,2000]  # epoch bounds

    elif(dem_type == "exp"):
        start_p = [40000,20000,10000]
        growth_rates = [0,0, 0 ]
        dts =  [1000,2000]

    elif(dem_type == "saw"):
        start_p = [50000, 50000, 5000, 50000, 5000, 50000, 5000]
        growth_rates = [0., 0.01931,-0.00485047, 0.00121838, -0.000306044, 0.0000768747, 0.]
        dts =  [40.0, 159.243, 633.957, 2523.829, 10047.546, 40000.0]

    elif(dem_type == "sawSPS"):
        start_p = [14312, 14312, 1431,14312, 1431, 14312,  1431]
        growth_rates = [0., 0.023025,-0.005756,0.0014391,-0.00035977, .0000899448, 0.]
        dts =  [33.333, 133.33, 533.33, 2133.33, 8533.33, 34133.31]
        
    elif(dem_type == "pwcsaw"):
        start_p = [50000, 50000, 15811, 5000, 15811, 50000, 15811, 5000, 15811,50000,15811,5000,5000]
        growth_rates =[0,0,0,0,0,0,0,0,0,0,0,0,0]
        dts = [28.3178, 56.5015, 112.735, 224.937, 448.807, 895.488, 1786.73, 3565., 7113.12, 14192.5, 28317.8, 56501.5]

    elif(dem_type == "pwcsawX"):
        start_p = [50000, 50000, 15811, 5000, 15811, 50000, 15811, 5000, 15811,50000,15811,5000,15811,50000,15811,5000]#,15811]
        growth_rates =[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]#,0]
        dts = [28.3178, 56.5015, 112.735, 224.937, 448.807, 895.488, 1786.73, 3565., 7113.12, 14192.5, 28317.8, 56501.5, 112735.2, 224936.0, 448806.0]#,895484.2]

    elif(dem_type == "bottleExp"):
        start_p = [24365,2000,10000]
        growth_rates =[.0025,0,0]
        dts = [1000,4000]

    else:
        print("Demography type not recognized")
        exit()
    
    # initialize population
    p0  = msprime.PopulationConfiguration(sample_size = n_samples,  initial_size= start_p[0], growth_rate = growth_rates[0])


    # create list of all demographic change events
    dem_events = []
    # add all size change events
    for i in range(len(dts)):
        dem_events.append(msprime.PopulationParametersChange(dts[i],start_p[i+1],growth_rate = growth_rates[i+1] ))
    

    # print to terminal the parameters for the population size for every epoch to double check
    dd = msprime.DemographyDebugger(  population_configurations=[p0], migration_matrix=[[0]],  demographic_events=dem_events)
    dd.print_history()



    # define psh function
    def psh(time):
        ii = 0
        tt = time
        for te in dts:
            if(tt < te): break
            ii = ii + 1
        
        if(ii > 0):
            tt = time - dts[ii - 1]
        
        return start_p[ii] * np.exp(-(tt)*growth_rates[ii] )


    # Write ordered pairs of truth function so we can plot them
    x_vals = np.asarray(range(int(dts[0]/40),int(6*dts[-1])))
    y_vals = np.copy(x_vals)
    for i in range(len(x_vals)):
        y_vals[i] = psh(x_vals[i])
    ops = np.matrix([x_vals,y_vals]).getT()
    df = pd.DataFrame(ops, columns = ['t (gens)', 'N(t)']) 

    f = open(truth_f,"w+")
    f.write(df.to_csv(index = False))
    f.close()




    # simulate and write to vcf files
    for nn in range(replicates):

        file_rep = sim_file_prefix
        if replicates > 1:
            file_rep = file_rep + str(nn+1)
        
        tree_sequence = msprime.simulate(length= genome_length, population_configurations= [p0], demographic_events= dem_events, recombination_rate= rec_rate, mutation_rate= mut_rate ,  random_seed=np.random.randint(10000000))
        with open(file_rep +".rvcf", "w") as vcf_file:
            tree_sequence.write_vcf(vcf_file,2)


        
        # store the relevant stats for trees in sequence to a csv
        tree_heights = []
        tree_lengths = []
        tree_spans = []
        for tree in tree_sequence.trees():
            # collect stats for each tree 
            tree_heights.append(  tree.get_time(tree.get_root())   )
            tree_lengths.append(   tree.get_total_branch_length()   )
            tree_spans.append(   tree.get_length()   )
            
        # write the tree sequence stats to file
        out_df = pd.DataFrame(data={'tree_spans':tree_spans, 'tree_heights':tree_heights, 'tree_lengths':tree_lengths})
        if(treeSeq_switch):
            out_df.to_csv(file_rep + "_TS.csv", index=False)
        
    
        print(f"Completed {file_rep} !")

        ancestral_derived_Fix(file_rep)

        
    
    print("done")
    print(time.time() - start , " seconds to complete")



def ancestral_derived_Fix(rawvcf):

    outvcf = open(rawvcf+ ".vcf","w+")
    
    with open(rawvcf+".rvcf") as infile:
        for line in infile:
            if line.find('#') == 0:
                outvcf.write(line)
            else:
                vars = line.split()
                
                ##
                vars[3] = 'A'
                vars[4] = 'T'
                ##

                out_l = ""
                for v in vars:
                    out_l = out_l + v + "\t"
                outvcf.write(out_l + "\n")
                
    outvcf.close()
    os.remove(rawvcf+'.rvcf')

    
if __name__ == "__main__":
    main()
