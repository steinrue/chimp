import java.io.File;
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

public class Wrapper_Inference {
	
	Helper helper = new Helper();
	
	public static void main(String[] args) throws IOException, JSAPException {
		Wrapper_Inference test= new Wrapper_Inference(); 

		test.infer_PSH_1P(test.helper.wrap_arg_parser(args,"PSH Demographic Inference Method."));
		
		/*
		String dir = "/Users/gautam/Desktop/test/";
		String[] arggg = new String[] { 
				"--rec_rate=.0000000125",  "--mut_rate=.0000000125"//, "--ps_scale=1000"
				,"--base_n=10", "--n_groups=20"
				//,"--psh0=8000"
				,"--partitions=.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02"				
				//,"--partitions=.1,.1,.1,.1,.1,.1,.1,.1,.1,.1"
				,"--t_bounds=30,60000"//, "--dof=1"
				//,"--tree_length",  "--pde_res=1000,250"
				,"--em_cap=50"//,"--m_evals=100"
				, "--s_type=0", "--simplex_scale=.1", "--ll_converge=.02"
				,"--metalocus=500"
				//,"--binary_emission"
				,"--data_dir="+dir
				,"--vcf_list=bottleExp_dataset2.vcf"
				,"--chr_l=200000000", "--ref_list=sim_ref.fasta", "--anc_list=sim_anc.fasta"
				,"--out_file="+dir+"deminf"
				//,"--print_all
			
				
		} ;
		test.infer_PSH_1P(test.helper.wrap_arg_parser(arggg,"PSH Demographic Inference Method."));
		
		*/
				
		
	}

	
	
	/////////////////////////////////////////////////////////////////////////////////////////////

	void infer_PSH_1P(JSAPResult arg_box) throws IOException {
		
		
		////////////////////////////////////////////////////////
		// BASIC PARAMETERS FOR ANALYSIS ///////////
		////////////////////////////////////////////
		
		double rec_rate =  arg_box.getDouble("recombination_rate");
		double mut_rate =  arg_box.getDouble("mutation_rate"); 
		int n_samples = arg_box.getInt("base_n_samples");
		int num_groups_of_samples = arg_box.getInt("sub_sample_groups"); 
		double[] lambdas = arg_box.getDoubleArray("regularizing_parameters"); 
		
		
		int n_states = arg_box.getInt("num_chmm_states"); // if this is overriden by partition list, then there's no need to update this definition since its not used.
		double psh0 = arg_box.contains("initial_size_guess") ? arg_box.getDouble("initial_size_guess") : -1 ;
		double sfactor = arg_box.getDouble("population_rescale_factor"); 
		int max_tract_size = arg_box.getInt("max_tract_size") ;
		int mL_size = arg_box.getInt("metalocus_size");
		
		///////////////////
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");
		System.out.println("Recombination Rate = " + rec_rate);
		System.out.println("Mutation Rate = " + mut_rate);
		System.out.println("Samples Handled = " + num_groups_of_samples + " groups of " + n_samples + " samples");
		
		if(arg_box.getBoolean("tree_length_switch")? n_samples > 10 : n_samples > 30 ) {
			System.out.println("WARNING: base_n higher than recommended. Efficiency of method may be severely decreased.");
		}
		
		System.out.println("Regularizing Parameters = " + Arrays.toString(lambdas));
		System.out.println("Internal Rescale factor = "+sfactor);
		
		if(mL_size < 1) {System.out.println("metalocus_size must be >0"); System.exit(0);}
		else if (mL_size >1) {System.out.println("Using MetaLocus FB implementation: \n    mL_size = "+ mL_size );}
		else {System.out.println("mL_size = 1 ; Using Locus Skipping Algorithm, (skip up to " + max_tract_size + " loci).");}
		
		System.out.println( psh0<0? "Use Waterson's Estimate for Initial PSH size" : "Initial PSH size = "+ psh0  );
		System.out.println(!arg_box.contains("partition_vector")? "Number of CHMM states (equipartitioned) = "+n_states : "CHMM states partitioned by "+Arrays.toString(arg_box.getDoubleArray("partition_vector")));
		
		////////////
		
		if(arg_box.getBoolean("tree_length_switch")) {System.out.println("Perform inference using TREE LENGTH for state representation.");}
		else { System.out.println("Perform inference using TMRCA for state representation."); }
		
		if(arg_box.getBoolean("SFS_switch")) {System.out.println("Perform Inference Using SFS (no linkage information).");}
		if(arg_box.getBoolean("binary_emission")) {System.out.println("Use binary classification (\"seg\" or \"non-seg\" as emission) ");}
		else {System.out.println("Use number of derived alleles as emission.");}
		if(arg_box.getBoolean("pseudo_haploid")) {System.out.println("Treat data as pseudo-haploid data.");}
		
		
		
		
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");

		///////////////////

		
		////////////////////////////////////////////////////////
		// READ VCF FILES, STORE STREAMS ///////////
		////////////////////////////////////////////
		
		///////////////////
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");
		System.out.println("Read VCF Files Into Streams");
		if(arg_box.contains("temporary_file_prefix")) {	System.out.println("Intermediate Data Structures will be written to disc with prefix '" +  arg_box.getString("temporary_file_prefix")+ "'." );}
		else { System.out.println("Intermediate Data Structures will be held in RAM");}
		///////////////////

		
		// setup lists of vcf files and subsampling parameters etc for each stream.
		
		String data_folder = arg_box.contains("data_file_directory")? arg_box.getString("data_file_directory") + "/" : ""; // only include if given, otherwise just search for files here
		String[] i_vl = arg_box.getStringArray("vcf_file_list");
		String[] i_al = arg_box.getStringArray("anc_file_list");
		String[] i_rl = arg_box.getStringArray("ref_file_list");
		
		// list of lengths, initialize to zero if not specified (will then read them in from ref files)
		int[] i_bl = new int[i_vl.length]; for (int i =0 ; i < i_bl.length ; i++) {i_bl[i] = 0; }
		if(arg_box.contains("chromosomal_length_list")) {i_bl = arg_box.getIntArray("chromosomal_length_list");}
		
		// check
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
		
		if(arg_box.getBoolean("print_intermediate_psh")) {
			File dir = new File(arg_box.getString("output_file"));
			dir.mkdir();
		}
		
		int count_j = 0;
		for(int i = 0 ; i < i_vl.length ; i++) {// cycle through all chromosome vcf files
			
			for(int ii = 0 ; ii < num_groups_of_samples ; ii++) {	// cycle through all appropriate subsampling intervals
				chrom_files[count_j] = data_folder + i_vl[i];
				ref_files[count_j] = data_folder + i_rl[i];
				anc_files[count_j] = data_folder + i_al[i];
				chrom_lengths[count_j] = i_bl[i];
				file_subsamplings[count_j] = samp_sub_intervals[ii].clone();
				count_j++;
			}
			
		}
		
		String temp_streams_prefix = null;
		if(arg_box.contains("temporary_file_prefix")) {temp_streams_prefix = arg_box.getString("temporary_file_prefix");}
		
		Helper.dataStreamsHandler dataStreams = helper.new dataStreamsHandler(chrom_files, ref_files, anc_files,
				chrom_lengths, file_subsamplings,
				max_tract_size, mL_size ,n_samples, arg_box.getBoolean("binary_emission"), temp_streams_prefix);
		
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");

		
		
		
		
		
		
		
		
		
		
		
		////////////////////////////////////////////////////////////////////
		// GET WATERSONS ESTIMATE, /////////////////
		// GET INITIAL PSH SIZE, ///////////////
		// CREATE STATE PARTITIONS ///////////
		////////////////////////////////////////////

		if(psh0 < 0) {
			System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");
			System.out.println("Compute Waterson's Estimate");
			System.out.println(dataStreams.getFSS());
			psh0 = 	(.5)/helper.estimate_watersons_CR(dataStreams.getFSS(), n_samples,  2*mut_rate);
			System.out.println("Initialize Population Size = " + psh0);
		
			System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n\n");
		}



		
		//////
		// Initialize uniform partitions if unspecified, or follow vectors if specified
		RealVector partitions;
		double[] part_vec_temp = new ArrayRealVector(n_states, 1./n_states).toArray();
		if(arg_box.contains("partition_vector")) {part_vec_temp = arg_box.getDoubleArray("partition_vector");}
		
		
		partitions = arg_box.getBoolean("tree_length_switch")? helper.getTreeLengthPartitions(part_vec_temp, n_samples):helper.getTreeHeightPartitions(part_vec_temp, n_samples) ;
		partitions.mapMultiplyToSelf(psh0 * 2);
		
		
		
	
		
		
		
		

		//////////////////////////////////////////////////////////////////////////
		// SCALE ALL DEMOGRAPHIC PARAMETERS             /////////////////
		// INITIALIZE PSH STRUCTURE WITH SCALED PARAMS /////////////////
		///////////////////////////////////////////////////////////////	
		
		
		// Scale all demographic parameters 
		
		double cr0 = sfactor / (2. * psh0);
		double theta = 2 * sfactor * mut_rate;
		double rho = 2 * sfactor * rec_rate;
		double[] scaled_partitions = partitions.mapDivide(sfactor).toArray(); // these are in units of coalescent time = 2*N_ref generations
		
		//if(arg_box.contains("initial_psh_ys") && ( arg_box.getBoolean("splines")|| arg_box.contains("time_bounds") || arg_box.contains("degrees_of_freedom"))) {System.err.println("time_bounds + dof specified in addition to initial psh shape.");System.exit(0);}		
	
		
		////////////////
		RealVector time_bounds;		
		if(arg_box.contains("time_bounds")){
			time_bounds = new ArrayRealVector(arg_box.getDoubleArray("time_bounds")); 
			if(time_bounds.getDimension()!=2) {System.err.println("Should have two time_bounds elements.");System.exit(1);}
			System.out.println("Min/Max Time Bounds Provided: "+time_bounds);
		}
		else {
			time_bounds = new ArrayRealVector(new double[] {psh0 / 50. , psh0 * 5.});
			System.out.println("Min/Max Time Bounds Inferred Based on Initial Population Size: "+ time_bounds);
		}
		int num_dof = arg_box.getInt("degrees_of_freedom");
		System.out.println("Number of DOF for PSH parametrization is: " + num_dof );

		////////////////
		
		FunctionUnivariate cr = null;

		if(arg_box.contains("initial_psh_xs") ) {
			double[] psh_xs = arg_box.getDoubleArray("initial_psh_xs");
			double[] psh_ys = arg_box.contains("initial_psh_ys")?  arg_box.getDoubleArray("initial_psh_ys"): new ArrayRealVector(arg_box.getDoubleArray("initial_psh_xs").length + 1, psh0).toArray() ;
			num_dof = psh_xs.length-1;
			
			System.out.println("PSH epochs provided; Overriding time_bounds with: "+ Arrays.toString(psh_xs));
			if(arg_box.contains("initial_psh_ys")) {
				System.out.println("Initial history Overridden with: " + Arrays.toString(psh_ys));
			}
			System.out.println("Number of DOF overriden with: " + num_dof);
			
			Function_PWC crtt = new Function_PWC(psh_xs, psh_ys, !arg_box.getBoolean("linear_time"));
			crtt.reciprocalSelf();
			crtt.multiplyScalarSelf(sfactor / 2.0);
			cr = crtt;
		}
		
		
		else if(arg_box.getBoolean("splines")) { // create spline functions
			if(arg_box.getBoolean("linear_time")) {
				cr = new Function_B3S(time_bounds.getEntry(0), time_bounds.getEntry(1), num_dof, cr0);
			}
			else { // this is the default
				cr = new Function_B3S_LX(time_bounds.getEntry(0), time_bounds.getEntry(1), num_dof, cr0);
			}
		}
		
		else {
			cr = new Function_PWC(time_bounds.getEntry(0),time_bounds.getEntry(1), num_dof, cr0, !arg_box.getBoolean("linear_time"));
		}
			
		
		System.out.println("Epoch divisions are: "+Arrays.toString(cr.getDiscontinuities()));
		
		cr.dilateXSelf(1./sfactor);
				
		
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n\n");

		
		
		/////////////////////////////////////////////////////////////////////////////////////////
		// SETTING UP HMM LEARNER CLASS, WITH THE APPROPRIATE DATA FILES AND PARAMETERS /////////
		/////////////////////////////////////////////////////////////////////////////////////////

		Inferer learner = null;

		// dont allow both SFS and TL options
		if(arg_box.getBoolean("SFS_switch") && arg_box.getBoolean("tree_length_switch")) {
			System.out.println("Cannot Specify both TL and SFS options."); System.exit(0);
		}
		
		
		// for tree height
		if(!arg_box.getBoolean("tree_length_switch")) {
			if(arg_box.getBoolean("SFS_switch")) {rho = -rho;} // TH computer will handle this appropriately.
			learner = new Inferer(cr, rho, theta, n_samples, scaled_partitions, "TH_1P", mL_size, null, arg_box.getBoolean("binary_emission"), arg_box.getBoolean("pseudo_haploid"));
		}
		// for tree length
		else {
			int[] pde_rez = new int[] {};
			if(arg_box.contains("pde_resolution")) {   pde_rez = arg_box.getIntArray("pde_resolution");
				if(pde_rez.length != 2 && pde_rez.length != 3) {System.err.print(" Pde resolution must have 2 or 3 values.");System.exit(1);}
			}
			learner = new Inferer(cr, rho, theta, n_samples, scaled_partitions, "TL_1P", mL_size, pde_rez, arg_box.getBoolean("binary_emission"), arg_box.getBoolean("pseudo_haploid"));
		}
		learner.update_reg_weights(lambdas); System.out.println("\n");
		
		
		// load DataStreams into learner
		learner.loadData(dataStreams);
		System.out.println("\n");
		
		/////////////////////////////////////////////////////////////////////////////////////////		
		// ACTUALLY DO LEARNING
		

		System.out.print("State partitions are :  ");
		System.out.println(partitions);
		System.out.print("IG distro across time bins is:  ");
		System.out.println(learner.get_SD());
		System.out.println("c_rate weights are : "+ Arrays.toString(cr.getFunctionParameters()) );
		System.out.println("over the epochs : "+ Arrays.toString(cr.getDiscontinuities()) );
		System.out.println();
		
		
		
		////////////////////
		double c_ll = learner.get_posteriorLL(); // current parameters posteriorLL
		System.out.println("STARTING LL is: " + Double.toString(c_ll));
		System.out.println("Initial Function's Regularization Quantities are: "+ learner.hmm_model.getRegularization(helper.estimate_watersons_CR(dataStreams.getFSS(), n_samples, theta)) );
		
		for(int i = 0 ; i < arg_box.getInt("max_em_steps") ; i++) {
			// if m_step_cap specified, use it, otherwise use heuristic to determine how many
			int m_evals = (arg_box.contains("m_step_cap")) ?  arg_box.getInt("m_step_cap") : 2 * num_dof + 10 ; // heuristic for num of m-evals per iteration.
			System.out.println("Max number of Function Evaluations in M step is: " + m_evals);


			// simplex scale diminishes each iteration unless option selected for fixed scale
			double simplex_scale = (arg_box.getBoolean("fixed_simplex_scale")) ? arg_box.getDouble("simplex_init_scale") : helper.emToSimpScale(arg_box.getDouble("simplex_init_scale"), i);

			
			// DO 1 EM step
			learner.em_single_update(m_evals, simplex_scale , arg_box.getInt("simplex_type"));


			/////// retrieve psh
			FunctionUnivariate temp_cr = cr.copy();
			temp_cr.updateFromParams(learner.get_current_params());

			System.out.println();
			System.out.println("c_rate weights are : "+ Arrays.toString(temp_cr.getFunctionParameters()) );
			System.out.println("\n --------------------------------");
			
			// print the function post EM step
			temp_cr.dilateXSelf(sfactor); temp_cr.multiplyScalarSelf(2./sfactor);
			temp_cr.printReciprocalToFile(arg_box.getString("output_file"), 100);
			if(arg_box.getBoolean("print_intermediate_psh")) {
				temp_cr.printReciprocalToFile(arg_box.getString("output_file") +"/EMstep"+ "_"+i, 400);
			}

			
			// obtain LL after EM
			double temp = learner.get_posteriorLL();
			System.out.println("after " + i + " iteration, ll is: " + Double.toString(temp)  );
			System.out.println("\n -------------------------------- \n -------------------------------- \n -------------------------------- ");
			
			///////////////////////////////
			///////////////////////////////
			// check convergence criterion
			
			if (Math.abs(temp - c_ll) < arg_box.getDouble("em_convergence_criterion")) {
				System.out.println("\n\n\n CONVERGENCE COMPLETE! \n\n\n");
				i = arg_box.getInt("max_em_steps");
			}
			c_ll = temp;
			
		}
		

		learner.clear_all_csvs();
		
		System.out.println("Done!!");
	}

	//////////////////////////////////////////

	
	private class Inferer extends EM_HMM {

		Inferer(FunctionUnivariate cr, double rr, double mr, int ns, double[] ts, String genealogy_rep, int mLsize,
				int[] additional_params, boolean binary_emission, boolean pseudohap) {
			super(cr, rr, mr, ns, ts, genealogy_rep, mLsize, additional_params, binary_emission, pseudohap);
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
