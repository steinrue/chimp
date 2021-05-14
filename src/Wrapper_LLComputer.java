
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.SimpleJSAP;
import com.martiansoftware.jsap.Switch;

public class Wrapper_LLComputer {
	
	Helper helper = new Helper();
	
	public static void main(String[] args) throws IOException, JSAPException {
		Wrapper_LLComputer test = new Wrapper_LLComputer(); 
		test.general_wrap(args);
	}

	
	
	/////////////////////////////////////////////////////////////////////////////////////////////

	void infer_PSH_1P(JSAPResult arg_box) throws IOException {
		Helper helper = new Helper(); 
		
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");
		double rec_rate =  arg_box.getDouble("recombination_rate"); System.out.println("Recombination Rate = " + rec_rate);
		double mut_rate = arg_box.getDouble("mutation_rate"); System.out.println("Mutation Rate = " + mut_rate);
		int n_samples = arg_box.getInt("base_n_samples"); 
		int num_groups_of_samples = arg_box.getInt("sub_sample_groups"); System.out.println("Samples Handled = " + num_groups_of_samples + " groups of " + n_samples + " samples");
		double[] lambdas = arg_box.getDoubleArray("regularizing_parameters"); System.out.println("Regularizing Parameters = " + Arrays.toString(lambdas));
		
		// initial value for psh
		double psh0 = arg_box.getDouble("initial_size_guess"); System.out.println("Initialize Population Size at size " + psh0);
		double sfactor = arg_box.getDouble("population_rescale_factor"); System.out.println("Internal Rescale factor = "+sfactor);
		
		
		
		/////// setup partitions 
		RealVector partitions = new ArrayRealVector(0);
		if(arg_box.getBoolean("tree_length_switch")) { partitions = helper.getTreeLengthPartitions(arg_box.getDoubleArray("partition_vector"), n_samples);	}
		else { 	partitions = helper.getTreeHeightPartitions(arg_box.getDoubleArray("partition_vector"), n_samples); }
		partitions.mapMultiplyToSelf(psh0 * 2);
		
		
		if(arg_box.contains("scratch_folder")) {	System.out.println("Intermediate Data Structures will be stored in '" +  arg_box.getString("scratch_folder")+ "' directory" );}
		else { System.out.println("Intermediate Data Structures will be held in RAM");}

		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n\n");
	
		
		
		
		

		/////////////////////////////////////////////////////////////////////////////////////////////
		
		// Scale all demographic parameters 
		psh0 = sfactor / (2. * psh0);
		double theta = 2 * sfactor * mut_rate;
		double rho = 2 * sfactor * rec_rate;
		double[] scaled_partitions = partitions.mapDivide(sfactor).toArray(); // these are in units of coalescent time = 2*N_ref generations
		
		
		
		// CREATE PSH STRUCTURE WITH SCALED PARAMETERS ACCORDING TO SPECIFIED MODEL
		Function_PWC crtt = new Function_PWC(arg_box.getDoubleArray("psh_ints"), arg_box.getDoubleArray("psh_vals"),!arg_box.getBoolean("linear_time"));
		crtt.reciprocalSelf();
		FunctionUnivariate cr = crtt;	
		cr.multiplyScalarSelf(sfactor/2.0);
		cr.dilateXSelf(1./sfactor);
		
		/////////////////////////////////////////////////////////////////////////////////////
		// SETUP LISTS OF FILE+SUBSAMPLING COMBOS TO FEED TO HMMLEARNER
		
		String data_folder = arg_box.getString("data_file_directory") + "/";
		String[] i_vl = arg_box.getStringArray("vcf_file_list");
		String[] i_al = arg_box.getStringArray("anc_file_list");
		String[] i_rl = arg_box.getStringArray("ref_file_list");
		int[] i_bl = arg_box.getIntArray("chromosomal_length_list");
		
		//
		if(i_vl.length != i_al.length || i_vl.length != i_rl.length || i_vl.length != i_bl.length ) { System.err.print("Input lists are not same length."); System.exit(1);}
		
		
		// get sample intervals ready (adjacent non overlapping intervals)
		int[][] samp_sub_intervals = new int[num_groups_of_samples][2] ;
		for(int i = 0 ; i < samp_sub_intervals.length ; i++) {samp_sub_intervals[i][0] = n_samples * i + 1; samp_sub_intervals[i][1] = n_samples * (i+1) ;}
		
		
		
		
		// get file names ready
		String[] chrom_files = new String[num_groups_of_samples * i_vl.length];
		String[] ref_files = new String[num_groups_of_samples * i_vl.length];
		String[] anc_files = new String[num_groups_of_samples * i_vl.length];
		int[] chrom_lengths = new int[num_groups_of_samples * i_vl.length];
		int[][] file_subsamplings = new int[num_groups_of_samples * i_vl.length][2]; 
		int count_j = 0;
		
		for(int i = 0 ; i < i_vl.length ; i++) {					// cycle through all chromosome vcf files
			
			for(int ii = 0 ; ii < num_groups_of_samples ; ii++) {	// cycle through all appropriate subsampling intervals
				chrom_files[count_j] = data_folder + i_vl[i];
				ref_files[count_j] = data_folder + i_rl[i];
				anc_files[count_j] = data_folder + i_al[i];
				chrom_lengths[count_j] = i_bl[i];
				file_subsamplings[count_j] = samp_sub_intervals[ii].clone();
				count_j++;
			}
			
		}
				
		
		
		
		/////////////////////////////////////////////////////////////////////////////////////////
		// SETTING UP HMM LEARNER CLASS, WITH THE APPROPRIATE DATA FILES AND PARAMETERS
		
		Inferer learner = null;
		// for tree height
		if(!arg_box.getBoolean("tree_length_switch")) {
			learner = new Inferer(cr, rho, theta, n_samples, scaled_partitions, "TH_1P", 0,null);
		}
		// for tree length
		else {
			int[] pde_rez = new int[] {};
			if(arg_box.contains("pde_resolution")) {   pde_rez = arg_box.getIntArray("pde_resolution");
				if(pde_rez.length != 2) {System.err.print(" Pde resolution must have 2 values.");System.exit(1);}
			}
			learner = new Inferer(cr, rho, theta, n_samples, scaled_partitions, "TL_1P", 0,pde_rez);
		}
		learner.update_reg_weights(lambdas); System.out.println("\n");
		
		// load data into files or memory depending on if scratch directory is specified or not
		if(arg_box.contains("scratch_folder")) {
			learner.loadVCFsToFiles(chrom_files, ref_files, anc_files, chrom_lengths, file_subsamplings, arg_box.getString("scratch_folder") + "/",1000);
		}
		else {	learner.loadVCFsToStreams(chrom_files, ref_files, anc_files, chrom_lengths, file_subsamplings,1000);	}
		System.out.println("\n");
		
		/////////////////////////////////////////////////////////////////////////////////////////		
		// COMPUTE LL, PRINT, AND PLOT INPUT FUNCTION
		

		System.out.print("State partitions are :  ");
		System.out.println(partitions);
		System.out.print("IG distro across time bins is:  ");
		System.out.println(learner.get_SD());
		System.out.println("c_rate weights are : "+ Arrays.toString(cr.getFunctionParameters()) );
		System.out.println();
		
		
		
		////////////////////
		double c_ll = learner.get_posteriorLL(); // current parameters posteriorLL
		System.out.println("LL is: " + Double.toString(c_ll));
		

		// print the function
		cr.updateFromParams(learner.get_current_params());
		cr.dilateXSelf(sfactor); cr.multiplyScalarSelf(2./sfactor);
		cr.printReciprocalToFile(arg_box.getString("output_file"), 1000);


		

		learner.clear_all_csvs();
		
		System.out.println("Done!!");
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
						
						// REGULARIZATION COEFFICIENT LAMBDA: --reg_lambda
						new FlaggedOption("regularizing_parameter", JSAP.DOUBLE_PARSER, "0" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "reg_lambda",
								"Regularizing coefficient used in likelihood optimization."),
						
						// INITIAL POPULATION SIZE HISTORY GUESS: --psh0
						new FlaggedOption("initial_size_guess", JSAP.DOUBLE_PARSER, "10000" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "psh0",
								"Initial guess for population size."),
						
						// TARGET PSH TO GET LL: --psh_v
						new FlaggedOption("psh_vals", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "psh_v",
								"Population sizes.").setList(true).setListSeparator(','),
						
						// TARGET PSH TO GET LL: --psh_x
						new FlaggedOption("psh_ints", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "psh_x",
								"PSH epoch bounds.").setList(true).setListSeparator(','),
						
						
						// INTERNAL POPULATION RESCALE FACTOR: --ps_scale
						new FlaggedOption("population_rescale_factor", JSAP.DOUBLE_PARSER, "1000" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ps_scale",
								"Internal Population Rescale Factor (modify if ODE stepsize error thrown)."),
						
						// HIDDEN STATE PROBABILITY PARTITIONS: --partitions
						new FlaggedOption("partition_vector", JSAP.DOUBLE_PARSER, ".05,.1,.1,.1,.1,.1,.1,.1,.1,.1,.05" ,
								JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "partitions",
								"Probability partitions to use for discrete binning of tree statistic into states.").setList(true).setListSeparator(','),
						
						
						

						////////////////////////// PARAMETERS FOR USING TREE LENGTH  /////////////////////
						
						// SWITCH TO USE TOTAL TREE LENGTH INSTEAD OF TMRCA
						new Switch("tree_length_switch", 'l', "tree_length", 
								"Use this switch to specify use of tree length instead of tree height."),
						
						// TREE LENGTH PDE RESOLUTION: --pde_res
						new FlaggedOption("pde_resolution", JSAP.INTEGER_PARSER, "250,1000" ,
								JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "pde_res",
								"Internal PDE solver resolution for discretization (specify discrete num of points for TL transition/emision pde solvers).").setList(true).setListSeparator(','),
						
						
						
						
						
						////////////////////////////// PARAMETERS REGARDING INPUT FILES ETC  ////////////////////////////////
						
						
						// SCRATCH FOLDER FOR DATASTREAM FILES: --scratch 
						new FlaggedOption("scratch_folder", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "scratch",
								"Only use option if processed vcf's should be temporarily stored to existing scratch directory instead of RAM."),
						
						// DIRECTORY WHERE DATA FILES ARE STORED: --data_dir
						new FlaggedOption("data_file_directory", JSAP.STRING_PARSER, "" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "data_dir",
								"Directory where data files are stored. (include trailing \"/\")"),
						
						// INPUT VCF DATA FILE LIST: --vcf_list
						new FlaggedOption("vcf_file_list", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "vcf_list",
								"List of VCF files for each chromosome to perform inferrence on. Assume sufficient samples/groups in each file.").setList(true).setListSeparator(','),
						
						// CHROMOSOMAL LENGTHS FOR VCFS: --chr_l
						new FlaggedOption("chromosomal_length_list", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "chr_l",
								"List of chromosome lengths for each vcf.").setList(true).setListSeparator(','),
						
						// ANCESTRAL DATA FILE LIST: --anc_list
						new FlaggedOption("anc_file_list", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "anc_list",
								"List of ancestral files for each chromosome. List should synch with vcfs.").setList(true).setListSeparator(','),
						
						// REFERENCE DATA FILE LIST: --ref_list
						new FlaggedOption("ref_file_list", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ref_list",
								"List of reference files for each chromosome. List should synch with vcfs.").setList(true).setListSeparator(','),
						
						// OUTPUT FILE: --out_file
						new FlaggedOption("output_file", JSAP.STRING_PARSER, "default_output_file" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "out_file",
								"Output file name where result of inference is stored."),
						
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
		
		
		
		infer_PSH_1P(arguments);
		//tf(arguments);
	}
	
	void tf(JSAPResult arg_box) {
		System.out.println(arg_box.getDouble("recombination_rate"));
	}
	
	//////////////////////////////////////////

	
	private class Inferer extends EM_HMM {

		Inferer(FunctionUnivariate cr, double rr, double mr, int ns, double[] ts, String genealogy_rep, int mLsize,
				int[] additional_params) {
			super(cr, rr, mr, ns, ts, genealogy_rep, mLsize, additional_params);
		}

		
		
		
		@Override
		protected double[] params_From_Hmm(HMmodel hmm) {
			double[][][] dem_params = hmm.getDemoParams();
			
			return dem_params[0][0];
		}

		@Override
		protected void params_To_HMM(double[] params, HMmodel hmm) {
			// TODO Auto-generated method stub
			double[][][] dem_params = hmm.getDemoParams();
			
			dem_params[0][0] = params; // sets the transition objects coalescence rate parameters
			dem_params[1][0] = params; // sets the emission object's coalescence rate parameters
			
			hmm.updateDemography(dem_params[0], dem_params[1]);
			
		}
		
	}

	/////////////////////////////////////////////////////////////////////////////////////////////



	
	
	
}
