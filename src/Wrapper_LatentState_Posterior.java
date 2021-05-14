import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.SimpleJSAP;
import com.martiansoftware.jsap.Switch;

import java.awt.*;

public class Wrapper_LatentState_Posterior {

	public static void main(String[] args) throws IOException, JSAPException {
		// TODO Auto-generated method stub
		/*
		//String folder = "/Users/gautam/Desktop/AdaptiveSweeps/Simulations/genicSelectionSweeps/posterior_inference_genic/";
		//String[] target_files = new String[] {"sim_154_sample1"};
		
		String folder = "/Users/gautam/Desktop/10-19_AdaptiveSweeps/RealData/Chrom2/";
		
		String target_file = "chr2_gbr";
		String ref_file = "chr2_ref";
		String anc_file = "human_ancestor_2";
		
		
		
		Wrapper_LatentState_Posterior test = new Wrapper_LatentState_Posterior(folder, target_file, ref_file, anc_file, target_file+"_TEST_");
		//test.InferPosterior();
		
		test.obtain_mode_distribution(folder +target_file+"_ml-10k_TL" + ".csv" , folder + "mode_out_TL.csv");
		*/
		
		Wrapper_LatentState_Posterior test = new Wrapper_LatentState_Posterior();
		test.general_wrap(args);
		
		System.out.println("DONE!!!"); Toolkit.getDefaultToolkit().beep();
	}
	
	Helper helper = new Helper();
	

	Random random = new Random();
	int id_earmark = random.nextInt(100000000);  // id number for temp files created by this EM_HMM object
	
	
	

	void InferPosterior() throws IOException {
		System.out.println("method deprecated");
		/*
		// parameters
		double rec_rate =  0.00000001;
		double mut_rate = 0.00000001;
		int n_samples = 11;
		int subsets = 9;
		
		int chrom_size = 200000000; // (220M for chr2, 170M for chr6)
		int metalocus_size = 10000;
		
		// initial value for Pop History (note this is the number of DIPLOID individuals)
		double psh0 = 20000;
		
		
		
		// Ne = 20k, pdf bins 1%,1%,3%,5%,10%...tree height then tree length, TH_N = 33, TL_N = 10
		//double[] dlengths = new double[] {21734., 24505.1, 29632.7, 35123.8, 43767., 51361., 59079.2, 67281.1, 77391.7, 89021.6, 106578., 134139., 162208., 198121., 227894.} ;
		double[] dlengths = new double[] { 73617.3, 83837.4, 101380., 120167., 145311., 167148., 188459., 208279., 232497., 262140., 298533., 357409., 415249., 492197., 543962.} ;

		
		
		double sfactor = 1000;
		
		/////////////////////////////////////////////////////////////////////////////////////////////
		
		// scale all demographic params and prepare to get transition and emission probs
		RealVector epochs = new ArrayRealVector(new double[] {10000});
		psh0 = sfactor/(2.* psh0);
		epochs.mapMultiplyToSelf(1./sfactor);
		FunctionUnivariate cr = new Function_PWC(epochs.toArray(), psh0);

		
		
		double theta = 2 * sfactor * mut_rate;
		double rho = 2 * sfactor * rec_rate;
		double[] scaled_lengths = new ArrayRealVector(dlengths).mapDivide(sfactor).toArray(); // these are in coalescent time = 2*N_ref generations
		
		
		
		
		// obtain the transition and emission probabilities given demographic parameters
		//TH_1P_Emissions ePcomp = new TH_1P_Emissions(cr, theta, n_samples, scaled_lengths);
		//TH_1P_Transitions tPcomp = new TH_1P_Transitions(cr, rho, n_samples, scaled_lengths);
		
		TL_1P_Emissions ePcomp = new TL_1P_Emissions(cr, theta, n_samples, scaled_lengths, 1000); System.out.println(ePcomp.getMargPDF());
		TL_1P_Transitions tPcomp = new TL_1P_Transitions(cr, rho, n_samples, scaled_lengths, new int[] {120,360});


		
		RealMatrix tP = tPcomp.transitionProbs();
		RealMatrix eP = ePcomp.emissionProbs();
		RealVector sD = ePcomp.getMargPDF();

		System.out.println("choice of time discretization is : ");
		System.out.println(Arrays.toString(dlengths));
		System.out.println("marginal distribution is : ");
		System.out.println(sD);
		
		

		/////////////////////////////////////////////
		// first compute posterior for all the sub samplings of the chrom, and write them to the temp files
		ArrayList<String> tempfiles = new ArrayList<String>();
		
		Expectation_FBskip posterior_computer = new Expectation_FBskip(tP,eP,sD);
		
		for(int i = 0 ; i < subsets ; i++) {
			int sample_lb = i * n_samples; int sample_rb = sample_lb + n_samples - 1;
			
			VcfReader v_file = new VcfReader(target_vcf , reference_file, ancestral_file, chrom_size, new int[] {sample_lb,sample_rb});
			posterior_computer.loadDataFromFile(v_file);
			
			String temp_file_name = "/Users/gautam/Desktop/scratch/"+ id_earmark+ "_temp_posterior_"+i ;
			posterior_computer.computePosteriorDistros(temp_file_name, dlengths,metalocus_size);
			
			tempfiles.add(temp_file_name);
		}
		
		// now combine all the temp files into one single final output file, and delete the temp files
		average_posteriors(tempfiles, output_file);
		*/
	}


	void InferPosterior(JSAPResult arg_box) throws IOException {
		// target file info
		String folder = arg_box.getString("data_file_directory");
		String target_vcf = folder + arg_box.getString("vcf_file");
		String ancestral_file = folder + arg_box.getString("anc_file");
		String reference_file = folder + arg_box.getString("ref_file");
		int chrom_size = arg_box.getInt("chromosomal_length");
		
		String out_file = arg_box.getString("output_file");

		
		// parameters
		double rec_rate =  arg_box.getDouble("recombination_rate"); System.out.println("Recombination Rate = " + rec_rate);
		double mut_rate = arg_box.getDouble("mutation_rate");  System.out.println("Mutation Rate = " + mut_rate);
		int n_samples = arg_box.getInt("base_n_samples");
		int num_groups_of_samples = arg_box.getInt("sub_sample_groups"); System.out.println("Samples Handled = " + num_groups_of_samples + " groups of " + n_samples + " samples");
		
		int metalocus_size = arg_box.getInt("metalocus_size"); System.out.println("MetaLocus Size is: "+metalocus_size+" bases");
		
		// initial value for Pop History (note this is the number of DIPLOID individuals)
		double psh = arg_box.getDouble("pop_size"); System.out.println("Initialize Population Size at size " + psh);
		
		
		/////// setup partitions 
		RealVector partitions = new ArrayRealVector(0);
		if(arg_box.getBoolean("tree_length_switch")) { partitions = helper.getTreeLengthPartitions(arg_box.getDoubleArray("partition_vector"), n_samples);	}
		else { 	partitions = helper.getTreeHeightPartitions(arg_box.getDoubleArray("partition_vector"), n_samples); }
		partitions.mapMultiplyToSelf(psh * 2);
		
		
		
		double sfactor = arg_box.getDouble("population_rescale_factor");   System.out.println("Internal Rescale factor = "+sfactor);

		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n\n");
		
		
		
		/////////////////////////////////////////////////////////////////////////////////////////////
		
		// scale all demographic params and prepare to get transition and emission probs

		psh = sfactor/(2.* psh);
		FunctionUnivariate cr = new Function_PWC(new double[] {1.0/psh}, psh); // CURRENTLY ONLY HANDLE CONSTANT PSH SPECIFICATION 

		double theta = 2 * sfactor * mut_rate;
		double rho = 2 * sfactor * rec_rate;
		double[] scaled_partitions = new ArrayRealVector(partitions).mapDivide(sfactor).toArray(); // these are in coalescent time = 2*N_ref generations
		
		
		
		
		
		// OBTAIN TRANSITION AND EMISSION PROBABILITIES FOR TREE LENGTH AND HEIGHT

		EmissionsComputer ePcomp = null;
		TransitionsComputer tPcomp = null; 
		
		if(arg_box.getBoolean("tree_length_switch"))
		{
			int[] pde_res = arg_box.getIntArray("pde_resolution");
			ePcomp = new TL_1P_Emissions(cr, theta, n_samples, scaled_partitions, pde_res[1]);
			tPcomp = new TL_1P_Transitions(cr, rho, n_samples, scaled_partitions, pde_res[0]);
		}
		else {
			ePcomp = new TH_1P_Emissions(cr, theta, n_samples, scaled_partitions);
			tPcomp = new TH_1P_Transitions(cr, rho, n_samples, scaled_partitions);
		}
		
		
		

		
		RealMatrix tP = tPcomp.transitionProbs();
		RealMatrix eP = ePcomp.emissionProbs();
		RealVector sD = ePcomp.getMargPDF();

		System.out.println("choice of latent variable partitioning is : ");
		System.out.println(Arrays.toString(partitions.toArray()));
		System.out.println("marginal distribution is : ");
		System.out.println(sD);
		
		

		/////////////////////////////////////////////
		// first compute posterior for all the sub samplings of the chrom, and write them to the temp files
		ArrayList<String> tempfiles = new ArrayList<String>();
		
		Expectation_FBskip posterior_computer = new Expectation_FBskip(tP,eP,sD);
		
		for(int i = 0 ; i < num_groups_of_samples ; i++) {
			int sample_lb = i * n_samples; int sample_rb = sample_lb + n_samples - 1;
			
			VcfReader v_file = new VcfReader(target_vcf , reference_file, ancestral_file, chrom_size, new int[] {sample_lb,sample_rb});
			posterior_computer.loadDataFromStream(v_file.getReducedStream(),1);
			
			String temp_file_name = arg_box.getString("scratch_folder")+ id_earmark+ "_temp_posterior_"+i ;
			posterior_computer.computePosteriorDistros(temp_file_name, partitions.toArray(),metalocus_size);
			tempfiles.add(temp_file_name);
		}
		
		// now combine all the temp files into one single final output file, and delete the temp files
		average_posteriors(tempfiles, out_file);
		
		
		// obtain the modes if asked for
		if(arg_box.contains("mode_file")) {  this.obtain_mode_distribution(out_file, arg_box.getString("mode_file")) ; }
	}



	void general_wrap(String[] args) throws IOException, JSAPException {
		SimpleJSAP jsap = new SimpleJSAP("PSH_Inference",
				"Performs demographic inference of population size history under single population model",
				new Parameter[] {
						
						/////////////////////////// GENERAL PARAMS ABOUT DEMOGRAPHY AND SAMPLING /////////////////////////
	
						// RECOMBINATION RATE: --rec_rate
						new FlaggedOption("recombination_rate", JSAP.DOUBLE_PARSER, ".00000001" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "rec_rate",
								"Per base per generation recombination rate."),
						
						// MUTATION RATE: --mut_rate
						new FlaggedOption("mutation_rate", JSAP.DOUBLE_PARSER, ".00000001" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "mut_rate",
								"Per base per generation mutation rate."),
						
						// SAMPLES IN SINGLE TREE: --base_n
						new FlaggedOption("base_n_samples", JSAP.INTEGER_PARSER, "10" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "base_n",
								"Number of Samples that the hidden state inherently handles (equivalently number of emissions)."),
						
						// NUMBER OF GROUPS OF SUBSAMPLES: --n_groups
						new FlaggedOption("sub_sample_groups", JSAP.INTEGER_PARSER, "1" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "n_groups",
								"Number of sub sample groups of base_n samples each that constitute effective sample size."),
						
						// INITIAL POPULATION SIZE HISTORY GUESS: --psh
						new FlaggedOption("pop_size", JSAP.DOUBLE_PARSER, "10000" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "psh",
								"Constant population size."),

						// HIDDEN STATE PROBABILITY PARTITIONS: --partitions
						new FlaggedOption("partition_vector", JSAP.DOUBLE_PARSER, ".05,.1,.1,.1,.1,.1,.1,.1,.1,.1,.05" ,
								JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "partitions",
								"Probability partitions to use for discrete binning of tree statistic into states.").setList(true).setListSeparator(','),
						
						// INTERNAL POPULATION RESCALE FACTOR: --ps_scale
						new FlaggedOption("population_rescale_factor", JSAP.DOUBLE_PARSER, "1000" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ps_scale",
								"Internal Population Rescale Factor (modify only if ODE stepsize error thrown)."),
						
						// # OF LOCI GROUPED IN POSTERIOR FILE : --ml_size
						new FlaggedOption("metalocus_size", JSAP.INTEGER_PARSER, "10000" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ml_size",
								"Number of bases that are grouped into meta-loci that compose the entries of the posterior file."),
						

						
						
						
						
						////////////////////////// PARAMETERS FOR USING TREE LENGTH  /////////////////////
						
						// SWITCH TO USE TOTAL TREE LENGTH INSTEAD OF TMRCA
						new Switch("tree_length_switch", 'l', "tree_length", 
								"Use this switch to specify use of tree length instead of tree height."),
						
						// TREE LENGTH PDE RESOLUTION: --pde_res
						new FlaggedOption("pde_resolution", JSAP.INTEGER_PARSER, "300,1000" ,
								JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "pde_res",
								"Internal PDE solver resolution for discretization (specify discrete num of points for TL transition/emision pde solvers).").setList(true).setListSeparator(','),
						
						
						
						
						
						////////////////////////////// PARAMETERS REGARDING INPUT FILES ETC  ////////////////////////////////
						
						
						// SCRATCH FOLDER FOR DATASTREAM FILES: --scratch 
						new FlaggedOption("scratch_folder", JSAP.STRING_PARSER, "", JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "scratch",
								"Name of folder where you want program to write intermediate files and delete them after use. "),
						
						// DIRECTORY WHERE DATA FILES ARE STORED: --data_dir
						new FlaggedOption("data_file_directory", JSAP.STRING_PARSER, "" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "data_dir",
								"Directory where data files are stored. (include trailing \"/\")"),
						
						// INPUT VCF DATA FILE: --vcf
						new FlaggedOption("vcf_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "vcf",
								"List of VCF files for each chromosome to perform inferrence on. Assume sufficient samples/groups in each file."),
						
						// CHROMOSOMAL LENGTH FOR VCF: --chr_l
						new FlaggedOption("chromosomal_length", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "chr_l",
								"List of chromosome lengths for each vcf."),
						
						// ANCESTRAL DATA FILE: --anc
						new FlaggedOption("anc_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "anc",
								"List of ancestral files for each chromosome. List should synch with vcfs."),
						
						// REFERENCE DATA FILE: --ref_list
						new FlaggedOption("ref_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ref",
								"List of reference files for each chromosome. List should synch with vcfs."),
						
						// OUTPUT FILE: --out_file
						new FlaggedOption("output_file", JSAP.STRING_PARSER, "default_posterior_file" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "out_file",
								"Output file name where result of inference is stored."),
						
						// OUTPUT FILE FOR MODE OF DISTRIBUTION: --mode_file
						new FlaggedOption("mode_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "mode_file",
								"Output file name for file where we store mode of distribution across genome. Not required"),
				}
		
		);
	
		
		
		JSAPResult arguments = jsap.parse(args); // obtain the holder which has parsed and knows all the arguments
		
		// Print Usage options if we run into issues, or help is tagged
		if (!arguments.success()) {
			
	        System.err.println();
	        System.err.println("Usage:  ");
	        System.err.println("                "+ jsap.getUsage());
	        System.err.println();
            System.err.println(jsap.getHelp());
	        System.err.println();

	        
			if(!arguments.success()) {
				for (java.util.Iterator errs = arguments.getErrorMessageIterator();
	                    errs.hasNext();) {
	                System.err.println("Error: " + errs.next());
	            }
		        System.err.println();
			}
			System.exit(1);

	    }
		
		
		//modify default values contingently
		
		
		this.InferPosterior(arguments);
		
	}
	
	
	//////////////////////////////////////////////////////////////////////////
	// Helper Functions ///
	///////////////////////
	
	void average_posteriors(ArrayList<String> t_files, String o_file) {
		
		RealMatrix running_posterior = new Array2DRowRealMatrix(  helper.read_doubleMat_FromCSV(t_files.get(0))  );
		helper.delete_csvFromDisc(t_files.get(0));
		
		// aggregate all the info from each file in a single matrix where we keep adding the values to
		for(int i = 1 ; i < t_files.size() ; i++) {
			RealMatrix next_mat = new Array2DRowRealMatrix( helper.read_doubleMat_FromCSV(t_files.get(i)) );
			running_posterior = running_posterior.add( next_mat) ;
			helper.delete_csvFromDisc(t_files.get(i));
		}
		
		// divide by number of files to get averaged vals. 
		running_posterior = running_posterior.scalarMultiply(1./ t_files.size());
		
		
		// write to the outfile
	
		helper.write_doubleMat_ToCSV(o_file, running_posterior.getData());
		
	}

	void obtain_mode_distribution(String inFile, String outFile) {
		BufferedReader br = null; FileWriter scribe = null;
        String line = "";
        String csvSplitBy = ",";

        try {

            br = new BufferedReader(new FileReader(inFile+ ".csv") );
    		scribe = new FileWriter(outFile + ".csv");
    		
    		String[] head_line = br.readLine().split(csvSplitBy); // This skips first line
    		//if(Double.parseDouble(head_line[0]) < 0) {System.out.print("Input Header Error");System.exit(0);}
    		
    		scribe.append("position, max_bin_num, lb_max_bin" + "\n");
    		while ((line = br.readLine()) != null) {
 
            	
                // use comma as separator
                String[] row_strings = line.split(csvSplitBy);
                if(row_strings.length == 1) {
                	continue;}
                
                double[] row_entries = new double[row_strings.length];
                for(int i = 0 ; i < row_strings.length ; i++) {row_entries[i] = Double.parseDouble(row_strings[i]);}
                
               
                
                // process and output
                
                int max_i = 1; double max_prob = row_entries[max_i];
                for(int i = 1 ; i < row_entries.length;i++) {
                	if (row_entries[i] >= max_prob) {
                		max_i = i; max_prob = row_entries[max_i];
                	}
                }
                
                scribe.append( row_entries[0] + ", " );
                scribe.append( max_i + ", ");
                
                String lower_bound_of_bin = "0";
                if(max_i == 1) {} else {lower_bound_of_bin = head_line[max_i];}
                
                scribe.append(lower_bound_of_bin + "\n");

            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (br != null) {
                try {
                    br.close(); scribe.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        
        System.out.println("done");
        
	}
	
	




}
