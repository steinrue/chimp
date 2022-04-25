
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.SimpleJSAP;
import com.martiansoftware.jsap.Switch;




public class Wrapper_LL {
	
	Helper helper = new Helper();
	
	public static void main(String[] args) throws IOException, JSAPException {
		
		Wrapper_LL test = new Wrapper_LL(); 

		test.infer_PSH_1P(test.wrap_arg_parser(args,"PSH Demographic Inference Method."));
	
	}

	
	
	/////////////////////////////////////////////////////////////////////////////////////////////

	void infer_PSH_1P(JSAPResult[] arg_box) throws IOException {
		
		
		////////////////////////////////////////////////////////
		// BASIC PARAMETERS FOR ANALYSIS ///////////
		////////////////////////////////////////////
		
		
		// RECOMBINATION AND MUTATION RATES
		double rec_rate =  arg_box[0].getDouble("recombination_rate");
		double mut_rate =  arg_box[0].getDouble("mutation_rate"); 
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");
		System.out.println("Recombination Rate = " + rec_rate);
		System.out.println("Mutation Rate = " + mut_rate);
		
		
		// DATA FILE LISTS AND CHROMOSOME LENGTHS
		String data_folder = arg_box[1].contains("data_file_directory")? arg_box[1].getString("data_file_directory") + "/" : ""; // only include if given, otherwise just search for files here
		String[] i_vl = arg_box[0].getStringArray("vcf_file_list");
		String[] i_al = arg_box[0].getStringArray("anc_file_list");
		String[] i_rl = arg_box[0].getStringArray("ref_file_list");
		int[] i_bl = new int[i_vl.length]; for (int i =0 ; i < i_bl.length ; i++) {i_bl[i] = 0; }
		if(arg_box[1].contains("chromosomal_length_list")) {i_bl = arg_box[1].getIntArray("chromosomal_length_list");}
		//check
		if(i_vl.length != i_al.length || i_vl.length != i_rl.length || i_vl.length != i_bl.length ) { System.err.print("Input lists are not same length."); System.exit(1);}
		for(int i = 0 ; i < i_vl.length;i++) {
			i_vl[i] = data_folder + i_vl[i];
			i_al[i] = data_folder + i_al[i];
			i_rl[i] = data_folder + i_rl[i];
		}
				
			
		
		// ARRANGE SAMPLES FOR ANALYSIS. create the int[][] sub_hap_indices	
		int[][] sub_hap_indices = null;

		if(arg_box[0].contains("subgroup_file")) { 
			String sub_group_file = arg_box[0].getString("subgroup_file");  System.out.println("Ignore <base_n> and <n_groups>. Samples Handled according to: "+ sub_group_file + "."); 
			sub_hap_indices = helper.hap_ind_file_read(sub_group_file);
		}
		else if (arg_box[0].contains("base_n_samples")) { 
			int[] n_samples = arg_box[0].getIntArray("base_n_samples");
			int[] num_groups_of_samples = null;
			if(arg_box[0].contains("sub_sample_groups")) { num_groups_of_samples = arg_box[0].getIntArray("sub_sample_groups");
				if(n_samples.length != num_groups_of_samples.length) {System.out.println("<base_n> and <n_groups> sizes should match.");System.exit(0);}
				System.out.print("Analyze non-overlapping groups of samples. Each of "+Arrays.toString(num_groups_of_samples)+" groups of "+Arrays.toString(n_samples)+" samples." );
			}
			else { int max_samples = helper.get_max_samples_from_VCFs(i_vl); // read num samples listed in the vcfs (they should match)				
				System.out.println("Analyze all samples in VCF in non-overlapping groups of n_s, for each n_s value in "+ Arrays.toString(n_samples));
				num_groups_of_samples = n_samples.clone();
				for(int i = 0 ; i < n_samples.length;i++) { num_groups_of_samples[i] = max_samples / n_samples[i];}
			}
			sub_hap_indices = helper.default_hap_groups(n_samples, num_groups_of_samples);

		}
		else {System.out.println("Must specify either <base_n> or <subgroup_file>");System.exit(0);}
		
		
		// PSEUDOHAPLOID, AND SIZE CHECK
		if(arg_box[0].getBoolean("pseudo_haploid")) {System.out.println("Treat data as pseudo-haploid data.");}

		// check if number of haps used in tree computation is greater than 10 for L, 30 for TH
		int max_ns = 0; for(int i = 0 ; i < sub_hap_indices.length;i++) {if(sub_hap_indices[i].length > max_ns) {max_ns = sub_hap_indices[i].length;}}
		
		if(max_ns * (arg_box[0].getBoolean("pseudo_haploid") ? 2 : 1 )  >  (arg_box[0].getBoolean("tree_length_switch")? 10 : 30) ) {
			System.out.println("WARNING: n_s value higher than recommended (either in <base_n> or <subgroup_file>). Efficiency of method may be severely decreased.");
		}
		
		
		// REGULARIZING PARAMETERS AND INTERNAL RESCALE FACTOR
		double[] lambdas = arg_box[1].getDoubleArray("regularizing_parameters"); 
		double sfactor = arg_box[1].getDouble("population_rescale_factor"); 
		System.out.println("Regularizing Parameters = " + Arrays.toString(lambdas));
		System.out.println("Internal Rescale factor = "+sfactor);
		
		
		
		// METALOCUS SIZE
		int mL_size = arg_box[0].getInt("metalocus_size");
		int max_tract_size = arg_box[1].getInt("max_tract_size");
		if(mL_size < 1) {System.out.println("metalocus_size must be >0"); System.exit(0);}
		else if (mL_size >1) {System.out.println("Using MetaLocus FB implementation: \n    mL_size = "+ mL_size );}
		else {System.out.println("mL_size = 1 ; Using Locus Skipping Algorithm, (skip up to " + max_tract_size + " loci).");}
		
	
			
		// INITIAL PSH GUESS
		double psh0 = arg_box[0].contains("initial_size_guess") ? arg_box[0].getDouble("initial_size_guess") : -1 ;
		System.out.println( psh0<0? "Use Waterson's Estimate for Initial PSH size" : "Initial PSH size = "+ psh0  );


		
		// TMRCA SWITCH
		if(arg_box[0].getBoolean("tree_length_switch")) {System.out.println("Perform inference using TREE LENGTH for state representation.");}
		else { System.out.println("Perform inference using TMRCA for state representation."); }
		
		
		
		// OTHER SWITCHES
		if(arg_box[1].getBoolean("SFS_switch")) {System.out.println("Perform Inference Using SFS (no linkage information).");}
		if(arg_box[1].getBoolean("binary_emission")) {System.out.println("Use binary classification (\"seg\" or \"non-seg\" as emission) ");}
		else {System.out.println("Use number of derived alleles as emission.");}
		
		
		
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");

		///////////////////

		
		
		////////////////////////////////////////////////////////
		// READ VCF FILES, STORE STREAMS ///////////
		////////////////////////////////////////////
		
		///////////////////
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");
		System.out.println("Read VCF Files Into Streams");
		if(arg_box[1].getBoolean("temporary_file_folder")) {	System.out.println("Intermediate Data Structures will be written to disc in temp folder.");}
		else { System.out.println("Intermediate Data Structures will be held in RAM");}
		///////////////////

		
		// for each stream get file names and sample indices ready
		ArrayList<String> chrom_files = new ArrayList<String>(); 
		ArrayList<String> ref_files = new ArrayList<String>(); 
		ArrayList<String> anc_files = new ArrayList<String>(); 
		ArrayList<Integer> chrom_lengths = new ArrayList<Integer>(); 
		ArrayList<int[]> file_subsamplings = new ArrayList<int[]>(); 
		
		// switch for intermediate psh files after every step
		if(arg_box[0].getBoolean("print_intermediate_psh")) {
			File dir = new File(arg_box[0].getString("output_file"));
			dir.mkdir();
		}
		
		// for each stream get appropriate information
		for(int i = 0 ; i < i_vl.length ; i++) {// cycle through all chromosome vcf files
			for(int ii = 0 ; ii < sub_hap_indices.length ; ii++) {	// cycle through all appropriate subsampling groups
				chrom_files.add(i_vl[i]);
				ref_files.add(i_rl[i]);
				anc_files.add(i_al[i]);
				chrom_lengths.add(i_bl[i]);
				file_subsamplings.add(sub_hap_indices[ii].clone());
			}	
		}
		
		// create prefix for temp files stored on disk if that's what's specified.
		String temp_streams_prefix = null;
		if(arg_box[1].getBoolean("temporary_file_folder")) {
			temp_streams_prefix = arg_box[0].getString("output_file") + "_TEMP";
        	Path temppath = Paths.get(temp_streams_prefix);
        	if (!Files.exists(temppath)) { 	Files.createDirectories(temppath); }
        	temp_streams_prefix = temp_streams_prefix + "/temp";
		}
		
		
		Helper.dataStreamsHandler dataStreams = helper.new dataStreamsHandler(chrom_files, ref_files, anc_files,
				chrom_lengths, file_subsamplings,
				max_tract_size, mL_size , arg_box[1].getBoolean("binary_emission"), temp_streams_prefix);
		
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");

		
		
		
		////////////////////////////////////////////////////////////////////
		// GET WATERSONS ESTIMATE, /////////////////
		// GET INITIAL PSH SIZE, ///////////////
		// CREATE STATE PARTITIONS ///////////
		////////////////////////////////////////////

		if(psh0 < 0) {
			System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");
			System.out.print("Compute Average Waterson's Estimate Across All Subgroups: ");
			System.out.println(dataStreams.get_Watterson());
			psh0 = 	dataStreams.get_Watterson()/(4 * mut_rate);
			System.out.println("Initialize Population Size = " + psh0);
		
			System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n\n");
		}

		
		
		//////////////////////////////////////////////////////////////////////////
		// SCALE ALL DEMOGRAPHIC PARAMETERS             /////////////////
		///////////////////////////////////////////////////////////////	
		
		
		// Scale all demographic parameters 
		
		double cr0 = sfactor / (2. * psh0);
		double theta = 2 * sfactor * mut_rate;
		double rho = 2 * sfactor * rec_rate;
		
		//if(arg_box[0].contains("initial_psh_ys") && ( arg_box[1].getBoolean("splines")|| arg_box[0].contains("time_bounds") || arg_box[0].contains("degrees_of_freedom"))) {System.err.println("time_bounds + dof specified in addition to initial psh shape.");System.exit(0);}		
	
		
		
		
		/////////////////////////////////////////////
		// INITIALIZE PARTITIONS /////////////
		
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n\n");

		// HIDDEN STATES, create partition vector
		double[] part_vec_temp = null;

		// HOLDER FOR PARTITION VALUES FOR EACH N VAL
		AbstractMap<Integer,double[]> scaled_partition_map= new HashMap<Integer,double[]>();
		int[] n_s_set = dataStreams.get_ns_list();

		// DONT ALLOW BOTH VALS AND PROBS TO BE SPECIFIED
		if(arg_box[1].contains("partition_vals") && arg_box[1].contains("partition_probs")) {
			System.out.println("Can only specify one of [partitions] or [partition_vals]");System.exit(0);
		}
		
		// IF PARTITION VALS PROVIDED, SCALE AND WRITE VALS INTO HOLDERS
		if(arg_box[1].contains("partition_vals")) { 
			part_vec_temp = arg_box[1].getDoubleArray("partition_vals"); // is array of partition values
			System.out.println("CHMM states partitioned by "+Arrays.toString(part_vec_temp));
			
			for(int i = 0; i < n_s_set.length ; i++) {		
				// index in the map is given by t_n_s, the number of haplotypes in each group, and the computation of probs is done with 
				// model_n, which is 2 * t_n_s for pesudo_haploids
				int t_n_s = n_s_set[i]; 
				
				// directly write the partition vals 
				scaled_partition_map.put(t_n_s, new ArrayRealVector(part_vec_temp).mapDivide(sfactor).toArray());
			}
		}
		
		// IF PARTITION PROBS PROVIDED, TREAT AS QUANTILES
		else if(arg_box[1].contains("partition_probs")) { 
			part_vec_temp = arg_box[1].getDoubleArray("partition_probs");  // is array of prob quantiles
			System.out.println("CHMM states partitioned by probs: "+Arrays.toString(arg_box[1].getDoubleArray("partition_probs")));
			System.out.println("State partitions are :  ");
			
			for(int i = 0; i < n_s_set.length ; i++) {
				// index in the map is given by t_n_s, the number of haplotypes in each group, and the computation of probs is done with 
				// model_n, which is 2 * t_n_s for pesudo_haploids
				int t_n_s = n_s_set[i]; int model_n = t_n_s * (arg_box[0].getBoolean("pseudo_haploid") ? 2 : 1);
					
				// create partition vector and compute
				RealVector t_part;
				t_part = arg_box[0].getBoolean("tree_length_switch")? helper.getTreeLengthPartitions(part_vec_temp, model_n):helper.getTreeHeightPartitions(part_vec_temp, model_n) ;
				t_part.mapMultiplyToSelf(psh0 * 2);
				
				System.out.println("n_s = "+t_n_s+" :--> "+ t_part); // print state partitions in generation units
				
				// rescale and store
				scaled_partition_map.put(t_n_s, t_part.mapDivide(sfactor).toArray());
			}
		}
		
		
		
		else { 
			
			// THIS ID DEFAULT IF NOTHING IS PROVIDED
			// EQUAL PROB QUANTILE DEFAULT - SAME FOR TH and TL
			/*
			int n_states = arg_box[0].getInt("num_chmm_states");
			part_vec_temp = new ArrayRealVector(n_states, 1./n_states).toArray();
			System.out.println( "Number of CHMM states (equipartitioned) = "+n_states );
			System.out.println("State partitions are :  ");
			
			for(int i = 0; i < n_s_set.length ; i++) {
				// index in the map is given by t_n_s, the number of haplotypes in each group, and the computation of probs is done with 
				// model_n, which is 2 * t_n_s for pesudo_haploids
				int t_n_s = n_s_set[i]; int model_n = t_n_s * (arg_box[0].getBoolean("pseudo_haploid") ? 2 : 1);
					
				// create partition vector and compute
				RealVector t_part;
				t_part = arg_box[0].getBoolean("tree_length_switch")? helper.getTreeLengthPartitions(part_vec_temp, model_n):helper.getTreeHeightPartitions(part_vec_temp, model_n) ;
				t_part.mapMultiplyToSelf(psh0 * 2);
				
				System.out.println("n_s = "+t_n_s+" :--> "+ t_part); // print state partitions in generation units
				
				// rescale and store
				scaled_partition_map.put(t_n_s, t_part.mapDivide(sfactor).toArray());
			}
			*/
			
			
			
			
			// NEW DEFAULT - DIFFERENT FOR TH and TL
			// CREATE LOG UNIFORM STATE BOUNDARIES, ENDPOINTS SPECIFIED BY TH and TL DIFFERENTLY
			int n_states = arg_box[0].getInt("num_chmm_states");
			System.out.println("State partitions are :  ");
			
			for(int i = 0; i < n_s_set.length ; i++) {
				// index in the map is given by t_n_s, the number of haplotypes in each group, and the computation of probs is done with 
				// model_n, which is 2 * t_n_s for pesudo_haploids
				int t_n_s = n_s_set[i]; int model_n = t_n_s * (arg_box[0].getBoolean("pseudo_haploid") ? 2 : 1);

				// get partition end-points for tree height - either from time_bounds, the default, or ps0_xs
				double[] part_ends = new double[] {psh0 / 50. , psh0 * 20. * 1.2}; //default
				if(arg_box[0].contains("time_bounds")) { part_ends =  arg_box[0].getDoubleArray("time_bounds");}		
				if(arg_box[0].contains("initial_psh_xs")) { double[] temptt = arg_box[0].getDoubleArray("initial_psh_xs"); part_ends = new double[] {temptt[0], temptt[temptt.length-1]} ;  }
				
				//part_ends = helper.getTreeHeightPartitions(new double[] {.001,.998,.001}, model_n).toArray(); part_ends = new ArrayRealVector(part_ends).mapMultiplyToSelf(psh0 * 2).toArray();
				part_ends = new double[] {100,500000};
				
				// partition boundaries for tree length
				if(arg_box[0].getBoolean("tree_length_switch")) {
					part_ends = helper.getTreeLengthPartitions(new double[] {.001,.998,.001}, model_n).toArray(); new ArrayRealVector(part_ends).mapMultiplyToSelf(psh0 * 2).toArray();
				}

				// log-uniform discretization of state bounds between ends
				RealVector t_part = new ArrayRealVector(helper.generate_log_discretization(part_ends[0], part_ends[1], n_states));
				System.out.println("n_s = "+t_n_s+" :--> "+ t_part); // print state partitions in generation units
				
				// rescale and store
				scaled_partition_map.put(t_n_s, t_part.mapDivide(sfactor).toArray());
			}

		}
		
			
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");

		
		
		
		///////////////////////////////////////////////////////////////////////////////////////
		// INITIALIZE PSH STRUCTURE WITH SCALED PARAMS /////////////////
		///////////////////////////////////////////////////////////////////////////////////////

		
		////////////////
		RealVector time_bounds;		
		if(arg_box[0].contains("time_bounds")){
			time_bounds = new ArrayRealVector(arg_box[0].getDoubleArray("time_bounds")); 
			if(time_bounds.getDimension()!=2) {System.err.println("Should have two time_bounds elements.");System.exit(1);}
			System.out.println("Min/Max Time Bounds Provided: "+time_bounds);
		}
		else {
			//time_bounds = new ArrayRealVector(new double[] {psh0 / 50. , psh0 * 5.});
			time_bounds = new ArrayRealVector(new double[] {psh0 / 50. , psh0 * 20.});
			
			System.out.println("Min/Max Time Bounds Inferred Based on Initial Population Size: "+ time_bounds);
		}
		int num_dof = arg_box[0].getInt("degrees_of_freedom");
		System.out.println("Number of DOF for PSH parametrization is: " + num_dof );

		////////////////
		
		FunctionUnivariate cr = null;

		if(arg_box[0].contains("initial_psh_xs") ) {
			double[] psh_xs = arg_box[0].getDoubleArray("initial_psh_xs");
			double[] psh_ys = arg_box[0].contains("initial_psh_ys")?  arg_box[0].getDoubleArray("initial_psh_ys"): new ArrayRealVector(arg_box[0].getDoubleArray("initial_psh_xs").length + 1, psh0).toArray() ;
			num_dof = psh_xs.length-1;
			
			System.out.println("PSH epochs provided; Overriding time_bounds with: "+ Arrays.toString(psh_xs));
			if(arg_box[0].contains("initial_psh_ys")) {
				System.out.println("Initial history Overridden with: " + Arrays.toString(psh_ys));
			}
			System.out.println("Number of DOF overriden with: " + num_dof);
			
			Function_PWC crtt = new Function_PWC(psh_xs, psh_ys, !arg_box[1].getBoolean("linear_time"));
			crtt.reciprocalSelf();
			crtt.multiplyScalarSelf(sfactor / 2.0);
			cr = crtt;
		}
		
		
		else if(arg_box[1].getBoolean("splines")) { // create spline functions
			if(arg_box[1].getBoolean("linear_time")) {
				cr = new Function_B3S(time_bounds.getEntry(0), time_bounds.getEntry(1), num_dof, cr0);
			}
			else { // this is the default
				cr = new Function_B3S_LX(time_bounds.getEntry(0), time_bounds.getEntry(1), num_dof, cr0);
			}
		}
		
		else {
			cr = new Function_PWC(time_bounds.getEntry(0),time_bounds.getEntry(1), num_dof, cr0, !arg_box[1].getBoolean("linear_time"));
		}
			
		
		System.out.println("Epoch divisions are: "+Arrays.toString(cr.getDiscontinuities()));
		
		cr.dilateXSelf(1./sfactor);
				
		
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n\n");

		
		
		/////////////////////////////////////////////////////////////////////////////////////////
		// SETTING UP HMM LEARNER CLASS, WITH THE APPROPRIATE DATA FILES AND PARAMETERS /////////
		/////////////////////////////////////////////////////////////////////////////////////////

		Inferer learner = null;

		// dont allow both SFS and TL options
		if(arg_box[1].getBoolean("SFS_switch") && arg_box[0].getBoolean("tree_length_switch")) {
			System.out.println("Cannot Specify both TL and SFS options."); System.exit(0);
		}
		
		
		// for tree height
		if(!arg_box[0].getBoolean("tree_length_switch")) {
			if(arg_box[1].getBoolean("SFS_switch")) {rho = -rho;} // TH computer will handle this appropriately.
			learner = new Inferer(cr, rho, theta, n_s_set, scaled_partition_map, "TH_1P", mL_size, null, arg_box[1].getBoolean("binary_emission"), arg_box[0].getBoolean("pseudo_haploid"));
		}
		// for tree length
		else {
			int[] pde_rez = new int[] {};
			if(arg_box[1].contains("pde_resolution")) {   pde_rez = arg_box[1].getIntArray("pde_resolution");
				if(pde_rez.length != 2 && pde_rez.length != 3) {System.err.print(" Pde resolution must have 2 or 3 values.");System.exit(1);}
			}
			learner = new Inferer(cr, rho, theta, n_s_set, scaled_partition_map, "TL_1P", mL_size, pde_rez, arg_box[1].getBoolean("binary_emission"), arg_box[0].getBoolean("pseudo_haploid"));
		}
		learner.update_reg_weights(lambdas); System.out.println("\n");
		
		
		// load DataStreams into learner
		learner.loadData(dataStreams);
		System.out.println("\n");
		
		
		
		
		
		
		/////////////////////////////////////////////////

		/////////////////////////////////////////////////////////////////////////////////////////		
		// COMPUTE LL, PRINT, AND PLOT INPUT FUNCTION
		

		System.out.println("Transition Matrix:  ");
		System.out.println(learner.hmm_model.getTProbs(n_s_set[0]));
		
		System.out.println("Emission Matrix:  ");
		System.out.println(learner.hmm_model.getEProbs(n_s_set[0]));
		System.out.println("\n \n \n \n ");
		
		
		
		System.out.print("IG distro across time bins is:  ");
		learner.print_marginals();
		System.out.println("c_rate weights are : "+ Arrays.toString(cr.getFunctionParameters()) );
		System.out.println("over the epochs : "+ Arrays.toString(cr.getDiscontinuities()) );
		System.out.println();
		
		
		
		////////////////////
		double c_ll = learner.get_posteriorLL(); // current parameters posteriorLL
		System.out.println("LL is: " + Double.toString(c_ll));
		System.out.println("Initial Function's Regularization Quantities are: "+ learner.hmm_model.getRegularization(theta / dataStreams.get_Watterson()));
		
		
		
		// print the function
		cr.updateFromParams(learner.get_current_params());
		cr.dilateXSelf(sfactor); cr.multiplyScalarSelf(2./sfactor);
		cr.printReciprocalToFile(arg_box[0].getString("output_file"), 1000);

		

		learner.clear_all_csvs();
		
		System.out.println("Done!!");
	}

	
	///////////////////////////////////////////

	// Args for dem inference
	JSAPResult[] wrap_arg_parser(String[] args, String help_description) throws IOException, JSAPException {
		
		// VISIBLE PARAMS //
		////////////////////
		
		SimpleJSAP jsap = new SimpleJSAP("Demographic Inference Method \n",
				help_description,
				new Parameter[] {
						
						
						//////////
						///////// REQUIRED PARAMETERS ///////////////////
						//////////

						
						
						// DEMOGRAPHIC PARAMS
						//////////
						
						// RECOMBINATION RATE: --rec_rate
						new FlaggedOption("recombination_rate", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "rec_rate",
								"Recombination rate, specified per generation per nucleotide (recommend .0000000125 for human data)."),
						
						// MUTATION RATE: --mut_rate
						new FlaggedOption("mutation_rate", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "mut_rate",
								"Mutation rate along a lineage, specified per generation per nucleotide (recommend .0000000125 for human data)."),
						
						// SAMPLES IN SINGLE TREE: --base_n
						new FlaggedOption("base_n_samples", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "base_n",
								"List of n_s values. CHIMP will analyze all the samples in the vcf as non-overlapping subgroups of n_s samples for all n_s values. These parallel CHMM runs are aggregated in the composite likelihood framework. We recommend using an n_s value of 10 at most (though for TMRCA, the method is still tractable up to 30). To analyze a subset of all the samples, use either <n_groups> or <hap_groups>.").setList(true).setListSeparator(','),
						
						
						// FILE PARAMS
						//////////


						// INPUT VCF DATA FILE LIST: --vcf_list
						new FlaggedOption("vcf_file_list", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "vcf_list",
								"List of VCF files for each chromosome to include in analysis. The samples across files should agree.").setList(true).setListSeparator(','),
							
						// REFERENCE DATA FILE LIST: --ref_list
						new FlaggedOption("ref_file_list", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ref_list",
								"List of reference files for each chromosome. Indexing and length of list should match [vcf_list] and [anc_list].").setList(true).setListSeparator(','),
						
						// ANCESTRAL DATA FILE LIST: --anc_list
						new FlaggedOption("anc_file_list", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "anc_list",
								"List of ancestral files for each chromosome. Can specify same files as [ref_list] if you don't have these (though accuracy may be affected if reference alleles and ancestral alleles are not the same).").setList(true).setListSeparator(','),
						
						// OUTPUT FILE: --out_file
						new FlaggedOption("output_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "out_file",
								"Prefix for output files CHIMP will produce with the inference results. CSV contains ordered pairs of [time, population size], and PARAM file contains, for each epoch, [epoch's lower bound (time), inferred population size]. Times are given in generations from present, and population sizes are in number of diploid individuals."),
						
						
						
						//////////
						///////// OPTIONAL VISIBLE PARAMETERS ///////////////////
						//////////
						

						// NUMBER OF GROUPS OF SUBSAMPLES: --n_groups
						new FlaggedOption("sub_sample_groups", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "n_groups",
								"A list of number of groups that correspond to the n_s values in <base_n>. For the i'th entry of [base_n], <n_groups>[i] non-overlapping groups of samples will be included in the analysis, indexed in the order they appear in the VCF.").setList(true).setListSeparator(','),
						
						// FILE FOR SPECIFYING THE HAPLOTYPE SUBGROUPS --hap_groups
						new FlaggedOption("subgroup_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "hap_groups",
								"Use this file to specify the exact haplotype subgroups used in the composite likelihood framework. Each row contains list of haplotype indices (positive integers, in order following VCF). Both <base_n> and <ngroups> get overwritten if this is specified."),
						
						// LEFT AND RIGHT BOUND OF TIME IN CONSIDERED HISTORY: --t_bounds
						new FlaggedOption("time_bounds", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "t_bounds",
								"This is used to specify epochs for population size history. The epoch boundaries will be spaced exponentially between <min_time> and <max_time>. Defaults to [Ne/50, 5*Ne] where Ne is computed from Waterson's estimator.").setList(true).setListSeparator(','),
						
						// NUMBER OF INDEPENDENT PARAMETERS TO USE IN MODELING PSH: --dof
						new FlaggedOption("degrees_of_freedom", JSAP.INTEGER_PARSER, "18" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "dof",
								"Number of epochs of population size history between [t_bounds] (with 2 additional epochs added at boundaries). One parameter is inferred for each epoch during EM. Note that dof=X will yield X+2 epochs."),
						
						// INITIAL POPULATION SIZE HISTORY GUESS: --psh0
						new FlaggedOption("initial_size_guess", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "psh0",
								"Initial guess for population size. Default is to use Waterson's estimate for N_effective, computed from data (and based on the specified mutation rate). Specifying this option will also use <psh0> instead of Waterson's estimate in computing default for [t_bounds] and the partitioning CHMM states."),
							
						// INITIAL POPULATION SIZE HISTORY DOMAIN: --psh0_xs
						new FlaggedOption("initial_psh_xs", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "psh0_xs",
								"Specify custom list of epoch boundaries (overrides [t_bounds], and <dof>) to use for population size history model.").setList(true).setListSeparator(','),
							
						// INITIAL POPULATION SIZE HISTORY VALS GUESS: --psh0_ys
						new FlaggedOption("initial_psh_ys", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "psh0_ys",
								"Specify population sizes (must also have specified [psh0_xs]) to initialize EM. This option overrides psh0. [psh0_ys] must be a vector with size 1 greater than that of [psh0_xs].").setList(true).setListSeparator(','),
						
						// NUMBER OF STATES IN CHMM: --n_states
						new FlaggedOption("num_chmm_states", JSAP.INTEGER_PARSER, "50" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "n_states",
								"Number of discrete states for CHMM (discrete intervals into which TMRCA or L falls). States are partitioned according to equal probability for each state under a constant population size prior (Waterson's estimate)."),
						
						// METALOCUS IMPLEMENTATION, BEST USED WITH TREELENGTH --metalocus_size
						new FlaggedOption("metalocus_size", JSAP.INTEGER_PARSER, "500" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "metalocus_size",
								"Number of bases grouped into each metalocus. Selecting metalocus_size=1 will activate a locus skipping algorithm to aid computational efficiency."),
						
						// MAX NUMBER OF EM STEPS: --em_cap
						new FlaggedOption("max_em_steps", JSAP.INTEGER_PARSER, "50" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "em_cap",
								"Maximum number of EM steps performed. Inference will end sooner if the convergence criterion is met."),
						
						// MAX NUMBER OF EVALS IN M STEP: --m_evals
						new FlaggedOption("m_step_cap", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "m_evals",
								"Maximum number of transition/emission matrix evaluations allowed during maximization step. Default scales with number of parameters (2 * <dof> + 10) since higher dimensional spaces will need more function evaluations to optimize effectively. M-step will terminate sooner if optimization converges."),
						
						// CRITERION FOR EM CONVERGENCE: --ll_converge
						new FlaggedOption("em_convergence_criterion", JSAP.DOUBLE_PARSER, ".02" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ll_converge",
								"Absolute threshold used to determine convergence of EM. If the posterior log-likelihood improves by less than this amount after any EM step, we consider the method converged."),
						
						// SWITCH TO USE TOTAL TREE LENGTH INSTEAD OF TMRCA: --tree_length
						new Switch("tree_length_switch", 'l', "tree_length", 
								"Use this to specify use of total branch length (L) instead of tree height (TMRCA) as representation of the CHMMs hidden state."),
					
						// PRINT ALL INTERMEDIATE STEPS: --output_steps
						new Switch("print_intermediate_psh", JSAP.NO_SHORTFLAG, "output_steps",
								"This option will produce new output files after each EM step, titled \"EMstep_X\" where X is the EM step number. The files will be placed in a new folder entitled <out_file>."),
						
						// SWITCH FOR PSEDUOHAPLOID EMISSION MODEL --pseudo
						new Switch("pseudo_haploid", JSAP.NO_SHORTFLAG, "pseudo", 
								"Use this switch to specify that the data is pseudohaploid data. CHIMP will assume each allele for each variant in the VCF is randomly selected from either haploid of a diploid individual (independently for each variant/individual). Emission probabilities will be adjusted accordingly."),
					
						//
						// SWITCH FOR PRINTING HELP FOR HIDDEN PARAMETERS --help_hidden
						new Switch("helpHidden", JSAP.NO_SHORTFLAG, "help_hidden", 
								"Print help descriptions for all the hidden options."),
						
					
				}
		
		
		);
		
		/////////////
		//////////// INVISIBLE PARAMETERS
		/////////////
		
		SimpleJSAP jsap_hidden = new SimpleJSAP("",
				"",
				new Parameter[] {
						
						///////////// SOMETIMES USEFUL

						// HIDDEN STATE PROBABILITY PARTITIONS: --partitions
						new FlaggedOption("partition_vector", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT,
								JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "partitions",
								"List of marginal probabilities to partition states of the CHMM. Under a constant population prior (from Waterson's estimate) the marginal probability of the TMRCA or L falling into each state will be given by [partitions]. The probabilities are indexed so that the first probability corresponds to the state with the smallest TMRCA (or L) and the last corresponds to the largest. [partitions] should sum to 1. This option overrides n_states.").setList(true).setListSeparator(','),
						
						// INTERNAL POPULATION RESCALE FACTOR: --ps_scale
						new FlaggedOption("population_rescale_factor", JSAP.DOUBLE_PARSER, "1000" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ps_scale",
								"Internal Population Rescale Factor. Modify this if ODE stepsize error thrown."),
						
						// MAXIMUM TRACT SIZE, REDUCE THIS TO AVOID UNDERFLOW IN FB ALGORITHM
						new FlaggedOption("max_tract_size", JSAP.INTEGER_PARSER, "1000" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "max_tract",
								"Maximum Tract Length of monomorphic tracts during locus skipping-algorithm. The locus-skipping algorithm is only used if <metalocus_size>=1. Reduce <max_tract> to avoid underflow errors during the E step."),
						
						// TREE LENGTH PDE RESOLUTION: --pde_res
						new FlaggedOption("pde_resolution", JSAP.INTEGER_PARSER, "1000,250" ,
								JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "pde_res",
								"Density of PDE solver's grid over the range of relevant times/lengths when using tree-length as hidden state. Higher density increases accuracy but also run-time. First parameter is for the emission probability solver (2D), and second is for transition probability solver (3D). This field is ignored for TMRCA.").setList(true).setListSeparator(','),
						// note can have one or two numbers for the TRANSITION rez (one means evenly spaced, 2 means even for each of the T/2, T/n epochs)
											
						// SIMPLEX SCALE PARAMETER: --simplex_scale
						new FlaggedOption("simplex_init_scale", JSAP.DOUBLE_PARSER, ".1" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "simplex_scale",
								"Size (edge length) of simplex used in M step. Since the search is in log space, 1.0 corresponds to an e-fold increase/decrease of the population size parameter."),
						
						// SWITCH FOR FIXING SIMPLEX SCALE OVER ITERATIONS: --fix_ss
						new Switch("fixed_simplex_scale", JSAP.NO_SHORTFLAG, "fix_ss", 
								"Use this switch to use same size of simplex for all EM steps (default is to slightly shrink the simplex on each successive step)."),
						
						// SCRATCH FOLDER FOR DATASTREAM FILES: --disc_storage 
						new Switch("temporary_file_folder", JSAP.NO_SHORTFLAG, "disc_storage",
								"Use option to store intermediate data structures temporarily on the disc (default stores them in RAM). Using this option may decrease RAM usage, but may increase run-time since structures will be stored to and read from disc."),
						
						// CHROMOSOMAL LENGTHS FOR VCFS: --chr_l
						new FlaggedOption("chromosomal_length_list", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "chr_l",
								"List of chromosome lengths for each VCF. Will analyze positions [1,length] for each file. Indices of [chr_l] should match those of [vcf_list]. If not provided, the reference file will determine the contiguous chromosomal tract that is analyzed.").setList(true).setListSeparator(','),
						
						
						///////////// NOT EXTENSIVELY TESTED
						
						// REGULARIZATION COEFFICIENT LAMBDA: --reg_lambdas
						new FlaggedOption("regularizing_parameters", JSAP.DOUBLE_PARSER, "0,0,0,0" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "reg_lambdas",
								"Regularizing coefficients used in likelihood optimization. Coefficients of [d0l2, d1l1, d1l2, d2l2] respectively. These coefficients are multiplied by [d0l2,d1l1,d1l2,d2l2] respectively and are subtracted from the objective function. d0l2 corresponds to the squared deviation from a constant population (Waterson’s estimate). d1l1 corresponds to the absolute value difference of the population size change between adjacent epochs, while d1l2 corresponds to the squared difference of the same. d2l2 is an analog for the square of the second derivative (computed for discrete epochs instead), and upweighting it will penalize deviations from linear behavior. The default is no regularization.").setList(true).setListSeparator(','),
						
						// SWITCH FOR REQUESTING SPLINE PSH MODEL: -s --spline
						new Switch("splines", JSAP.NO_SHORTFLAG, "spline", 
								"This option is not fully tested. It specifies that the population size history will be modeled by a cubic spline function (instead of piece-wise constant). The PARAM file will not be generated.  This can be used in conjunction with [t_bounds] and <dof>, which will specify the bounds of the spline and the number of nodes respectively. This should not be used with [psh0_xs] or [psh0_ys]."),
						
						// SWITCH FOR LINEAR TIME MODEL: --lin_t
						new Switch("linear_time", JSAP.NO_SHORTFLAG, "lin_t", 
								"This will distribute the epoch partitions uniformly on a linear scale rather than on a log scale."),
						
						///////////// SHOULD BE FULLY HIDDEN FROM USER
						
						// SIMPLEX TYPE: --s_type
						new FlaggedOption("simplex_type", JSAP.INTEGER_PARSER, "0" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "s_type",
								"ID for simplex type that EM_HMM uses at each em step. Default corresponds to an equilateral simplex centered on the current value, adjusted to include current value if it is better than the worst vertex."),
						
						// DIRECTORY WHERE DATA FILES ARE STORED: --data_dir
						new FlaggedOption("data_file_directory", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "data_dir",
								"Directory where data files are stored (include trailing \"/\"). Prepends this to all input files."),
						
						// SWITCH TO USE SFS FOR INFERENCE
						new Switch("SFS_switch", JSAP.NO_SHORTFLAG, "sfs", 
								"Use this switch to specify use of Site Frequency Spectrum (will not use linkage info). Computes marginal from TMRCA transition probability matrix. This is a highly inefficient implementation of SFS based demographic inference."),
						
						// SWITCH FOR BINARY EMISSION MODEL --binary_emission
						new Switch("binary_emission", JSAP.NO_SHORTFLAG, "binary_emission", 
								"Use this switch to specify model with only 2 emission possibilities, segregating or non-segregating."),
					
						
											
						//////////////////
						// IN TESTING
						/////////////////
						
	
				}
		
		);
		
		// disable errors since parsing will produce errors for each parser (of the args meant for the other parser)
		PrintStream bu_errS = System.err;
		System.setErr(new PrintStream(new OutputStream() { public void write(int b) { }}));
		
		// Parse with the two parsers
		JSAPResult arguments = jsap.parse(args); // obtain the holder which has parsed and knows all the arguments
		if(arguments.getBoolean("help")) {System.out.println(jsap.getHelp());}
		if(arguments.getBoolean("helpHidden")) { System.out.println(jsap_hidden.getHelp());}
		jsap_hidden.setHelp("");
		JSAPResult hidden_args = jsap_hidden.parse(args); // obtain the holder which has parsed and knows all the arguments
		
		
		// Restore Error Stream To the Backed-up one
		System.setErr(bu_errS);


		// EXIT AND PRINT ERRORS COMMON FOR BOTH PARSERS
		//Iterate over all errors from first parser
		for (@SuppressWarnings("rawtypes")
		java.util.Iterator errs = arguments.getErrorMessageIterator();
                errs.hasNext();) {
				String err1 = (String) errs.next();

				//Iterate over all errors from second parser
				for (@SuppressWarnings("rawtypes")
				java.util.Iterator errs_hd = hidden_args.getErrorMessageIterator();
		                errs_hd.hasNext();) {
					String err2 = (String) errs_hd.next();
					
					// If first and second parser share an error, stop and print it. Otherwise one of the parsers covers the argument, so do nothing
					if(err1.equals(err2)) {
						System.err.println("Error: " + err1);
						System.exit(0);
					}
					
				}
        }

		// MAKE SURE FIRST PARSER HAS ENOUGH INFO TO RUN
		boolean fail_Switch = false;
		java.util.Iterator badparams = arguments.getBadParameterIDIterator();
		while (badparams.hasNext()) {
			String p_name = (String) badparams.next();
			if (p_name != null){
				System.err.println("Error: Issue with <" + p_name +">");
				fail_Switch=true;
			}
		}
		if(fail_Switch) {System.err.println("Issue Parsing."); System.exit(0);}
		
		
		
		//modify default values contingently
		return new JSAPResult[] {arguments, hidden_args};
	}

	///////////////////////////////////////////
	
	
	
	private class Inferer extends EM_HMM {

		Inferer(FunctionUnivariate cr, double rr, double mr, int[] ns, AbstractMap<Integer,double[]> partitions, String genealogy_rep, int mLsize,
				int[] additional_params, boolean binary_emission, boolean pseudohap) {
			super(cr, rr, mr, ns, partitions, genealogy_rep, mLsize, additional_params, binary_emission, pseudohap);
		}

		
		
		
		@Override
		protected double[] params_From_Hmm(HMmodel hmm) {
			double[][] dem_params = hmm.getDemoParams();
			
			return dem_params[0];
		}

		@Override
		protected void params_To_HMM(double[] params, HMmodel hmm) {
			// TODO Auto-generated method stub
			double[][] dem_params = hmm.getDemoParams();
			
			dem_params[0] = params; // sets the transition objects coalescence rate parameters
			
			hmm.updateDemography(dem_params);
			
		}
		
	}

	/////////////////////////////////////////////////////////////////////////////////////////////

	
	
	
	
}
