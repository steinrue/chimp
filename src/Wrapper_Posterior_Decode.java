import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

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

public class Wrapper_Posterior_Decode {
	
	Helper helper = new Helper();
	
	public static void main(String[] args) throws IOException, JSAPException {

		Wrapper_Posterior_Decode test = new Wrapper_Posterior_Decode();
		test.computePosterior(test.wrap_posterior_args(args, "Posterior Decoding Method"));
		
		/*
		String dir = "/Users/gautam/Desktop/PD_diag/";
		String[] arggg = new String[] { 
				"--rec_rate=.0000000125",  "--mut_rate=.0000000125"
				,"--base_n=10", "--n_groups=18"
				,"--psh0=10000"
				,"--chr_l=200000000"
				//,"--tree_length",  "--pde_res=1000,250"
				,"--n_states=64"
				,"--psh_params="+dir+"GBR_DI_results.param"
				,"--vcf_file="+dir+"GBR_Chr2.vcf"
				,"--ref_file="+dir+"chr2.ref", "--anc_file="+dir+"homo_sapiens_ancestor_2.fa"
				,"--out_file="+dir+"PD_out"
			
		} ;
		test.computePosterior(test.wrap_posterior_args(arggg,"Posterior Decoding Method."));
		*/
		
	}


	
	
	/////////////////////////////////////////////////////////////////////////////////////////////

	void computePosterior(JSAPResult[] arg_box) throws IOException {
		
		
		////////////////////////////////////////////////////////
		// BASIC PARAMETERS FOR ANALYSIS ///////////
		////////////////////////////////////////////
		
		double rec_rate =  arg_box[0].getDouble("recombination_rate");
		double mut_rate =  arg_box[0].getDouble("mutation_rate"); 
		int n_samples = arg_box[0].getInt("base_n_samples");
		int model_n = arg_box[0].getBoolean("pseudo_haploid") ? 2 * n_samples : n_samples; // for pseudohaploid, trees are for 2n haplotypes
		int num_groups_of_samples = arg_box[0].getInt("sub_sample_groups");  // overwritten if sub_group_file is provided
		String sub_group_file = null; 

		
		int n_states = arg_box[0].getInt("num_chmm_states"); // if this is overriden by partition list, then there's no need to update this definition since its not used.
		double psh0 = arg_box[0].contains("psh_const") ? arg_box[0].getDouble("psh_const") : -1 ;
		double sfactor = arg_box[1].getDouble("population_rescale_factor"); 
		int max_tract_size = arg_box[1].getInt("max_tract_size") ;
		int mL_size = arg_box[0].getInt("metalocus_size");
		
		///////////////////
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");
		System.out.println("Recombination Rate = " + rec_rate);
		System.out.println("Mutation Rate = " + mut_rate);
		
		//check if sub_groups are specified, otherwise we'll use default
		if(arg_box[0].contains("subgroup_file")) { sub_group_file = arg_box[0].getString("subgroup_file"); num_groups_of_samples=-1;  System.out.println("Samples Handled according to "+ sub_group_file);  }
		else{ System.out.println("Samples Handled = " + num_groups_of_samples + " groups of " + n_samples + " samples"); }
				
		System.out.println("Samples Handled = " + num_groups_of_samples + " groups of " + n_samples + " samples");
		if(arg_box[0].getBoolean("pseudo_haploid")) {System.out.println("Treat data as pseudo-haploid data.");}

		if(arg_box[0].getBoolean("tree_length_switch")? model_n > 10 : model_n > 30 ) {
			System.out.println("WARNING: base_n higher than recommended. Efficiency of method may be severely decreased.");
		}
		
		System.out.println("Internal Rescale factor = "+sfactor);
		
		if(mL_size < 1) {System.out.println("metalocus_size must be >0"); System.exit(0);}
		else if (mL_size >1) {System.out.println("Using MetaLocus FB implementation: \n    mL_size = "+ mL_size );}
		else {System.out.println("mL_size = 1 ; Using Locus Skipping Algorithm, (skip up to " + max_tract_size + " loci).");}
		
		System.out.println( psh0<0? "Use Waterson's Estimate for Initial PSH size" : "Initial PSH size = "+ psh0  );
		System.out.println(!arg_box[1].contains("partition_vector")? "Number of CHMM states (equipartitioned) = "+n_states : "CHMM states partitioned by "+Arrays.toString(arg_box[1].getDoubleArray("partition_vector")));
		
		////////////
		
		if(arg_box[0].getBoolean("tree_length_switch")) {System.out.println("Perform inference using TREE LENGTH for state representation.");}
		else { System.out.println("Perform inference using TMRCA for state representation."); }
		
		
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");

		///////////////////

		
		////////////////////////////////////////////////////////
		// READ VCF FILES, STORE STREAMS ///////////
		////////////////////////////////////////////
		
		///////////////////
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n");
		System.out.println("Read VCF Files Into Streams");
		if(arg_box[1].getBoolean("temporary_file_folder")) { System.out.println("Intermediate Data Structures will be written to disc in temp folder.");}
		else { System.out.println("Intermediate Data Structures will be held in RAM");}
		///////////////////

		
		// setup lists of vcf files and subsampling parameters etc for each stream.
		
		String i_v = arg_box[0].getString("vcf_file");
		String i_a = arg_box[0].getString("anc_file");
		String i_r = arg_box[0].getString("ref_file");

		// chromosome length, initialize to zero if not specified (will then read them in from ref files)
		int i_bl = arg_box[1].contains("chromosome_length") ? arg_box[1].getInt("chromosome_length") : 0 ;
		

		// get sample intervals ready 
		int[][] sub_hap_indices = null ;
				
		if (sub_group_file != null) { // Read from the param file
					sub_hap_indices = helper.hap_ind_file_read(sub_group_file);
					if(n_samples != sub_hap_indices[0].length) {System.err.println("base_n should equal subgroup sizes in <hap_groups>"); System.exit(0);}; 
		}
		else { 		// Default (adjacent non overlapping intervals)
			sub_hap_indices = new int[num_groups_of_samples][n_samples] ;
			for(int i = 0 ; i < sub_hap_indices.length ; i++) {
						for(int j = 0; j < n_samples ; j++) { sub_hap_indices[i][j] = n_samples * i + j + 1;  }
			}
		}
		
		
		// get file names ready
		ArrayList<String> chrom_files = new ArrayList<String>(); 
		ArrayList<String> ref_files = new ArrayList<String>(); 
		ArrayList<String> anc_files = new ArrayList<String>(); 
		ArrayList<Integer> chrom_lengths = new ArrayList<Integer>(); 
		ArrayList<int[]> file_subsamplings = new ArrayList<int[]>(); 
		
			
		// get sample intervals ready (adjacent non overlapping intervals) and the filenames as well
		for(int i = 0 ; i < sub_hap_indices.length ; i++) {
			chrom_files.add(i_v);
			ref_files.add(i_r);
			anc_files.add(i_a);
			chrom_lengths.add(i_bl);
			file_subsamplings.add(sub_hap_indices[i]);
		}
		
		
		// SETUP TEMP FOLDER ALWAYS FOR SUBSAMPLED POSTERIORS
		String temp_folder = arg_box[0].getString("output_file") + "_TEMP";
		Path temppath = Paths.get(temp_folder);
    	if (!Files.exists(temppath)) { 	Files.createDirectories(temppath); }
		
    	// for temp streams check if we're storing on disc or in RAM
    	String temp_streams_prefix = arg_box[1].getBoolean("temporary_file_folder") ? temp_folder + "/temp_stream_" : null;

    	// setup stream handler
    	Helper.dataStreamsHandler dataStreams = helper.new dataStreamsHandler(chrom_files, ref_files, anc_files,
				chrom_lengths, file_subsamplings,
				max_tract_size, mL_size , false, temp_streams_prefix );
    	
		
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



		
		//////
		// Initialize uniform partitions if unspecified, or follow vectors if specified
		RealVector partitions;
		double[] part_vec_temp = new ArrayRealVector(n_states, 1./n_states).toArray();
		if(arg_box[1].contains("partition_vector")) {part_vec_temp = arg_box[1].getDoubleArray("partition_vector");}
		
		
		partitions = arg_box[0].getBoolean("tree_length_switch")? helper.getTreeLengthPartitions(part_vec_temp, n_samples):helper.getTreeHeightPartitions(part_vec_temp, n_samples) ;
		partitions.mapMultiplyToSelf(psh0 * 2);
		
		
		//////////////////////////////////////////////////////////////////////////
		// OBTAIN PSH INPUT
		// SCALE ALL DEMOGRAPHIC PARAMETERS             /////////////////
		// INITIALIZE PSH STRUCTURE WITH SCALED PARAMS /////////////////
		///////////////////////////////////////////////////////////////	
		
		
		// Scale all demographic parameters 
		
		double theta = 2 * sfactor * mut_rate;
		double rho = 2 * sfactor * rec_rate;
		double[] scaled_partitions = partitions.mapDivide(sfactor).toArray(); // these are in units of coalescent time = 2*N_ref generations
		
	
		////////////////
		
		// DEFAULT is to use constant pop based on Waterson's Estimate
		Function_PWC cr = new Function_PWC(new double[] {psh0}, psh0, true);
 

		// OVERWRITE WITH PSH FROM FILE IF PROVIDED IN .PARAM FILE
		if(arg_box[0].contains("psh_file")) {
			System.out.println("Reading population size history from: "+ arg_box[0].getString("psh_file")  );
			cr = new Function_PWC(arg_box[0].getString("psh_file"), true);
		}
		
		// USE PARAMS IF PROVIDED AS [psh_xs] AND [psh_ys]
		else if (arg_box[0].contains("psh_xs") && arg_box[0].contains("psh_ys")) {
			double[] psh_xs = arg_box[0].getDoubleArray("psh_xs");
			double[] psh_ys = arg_box[0].getDoubleArray("psh_ys");

			cr = new Function_PWC(psh_xs, psh_ys, true);
		}
		
		
		cr.reciprocalSelf();
		cr.multiplyScalarSelf(sfactor / 2.0);		
		cr.dilateXSelf(1./sfactor);
				
		
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n\n");


		
		///////////////////////////////////////////////////////////////	
		// OBTAIN TRANSITION AND EMISSION PROBABILITIES FOR TREE LENGTH AND HEIGHT


		
		EmissionsComputer ePcomp = null;
		TransitionsComputer tPcomp = null; 
		
		if(arg_box[0].getBoolean("tree_length_switch")) // TL impelementation
		// note that model_n is for total haplotypes in model, double base_n if they are pseudohaps
		{
			int[] pde_rez = arg_box[1].getIntArray("pde_resolution") ;
			
			// constructor for trans prob based on how many pde_rez args, 2 or 3
			if(pde_rez.length == 2) { 	tPcomp = new TL_1P_Transitions(cr, rho, model_n, scaled_partitions, pde_rez[1]);	}
			else if(pde_rez.length == 3) { tPcomp = new TL_1P_Transitions(cr, rho, model_n, scaled_partitions, new int[] {pde_rez[1],pde_rez[2]});}
			else {	System.err.print(" Pde resolution must have 2 or 3 values.");System.exit(1);}

			ePcomp = new TL_1P_Emissions(cr, theta, model_n, scaled_partitions, pde_rez[0], arg_box[0].getBoolean("pseudo_haploid"));
			
		}
		else { // TH implementation
			ePcomp = new TH_1P_Emissions(cr, theta, model_n, scaled_partitions, arg_box[0].getBoolean("pseudo_haploid"));
			tPcomp = new TH_1P_Transitions(cr, rho, model_n, scaled_partitions);
		}
		
		
		RealMatrix tP = tPcomp.transitionProbs().power(mL_size);
		RealMatrix eP = ePcomp.emissionProbs(false);
		RealVector sD = ePcomp.getMargPDF();
		

		System.out.println("Marginal Distribution Over States Is : ");
		System.out.println(sD);
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n\n");

		Expectation_FBskip posterior_computer = new Expectation_FBskip(tP,eP,sD);

		/////////////////////////////////////////////
	
		// first compute posterior for all the sub samplings of the chrom, and write them to the temp files
		ArrayList<String> tempfiles = new ArrayList<String>();
		
	
	
	    // Compute each subsamplings posterior and write to file
	    for(int i = 0 ; i < dataStreams.num_streams ; i++) {

			// load datastream
			int[][] dstream = dataStreams.getStream(i);
			posterior_computer.loadDataFromStream(dstream, mL_size);

			
			String temp_file_name = temp_folder + "/temp_posterior_"+i ;
			
			// this uses ML = mL_size, but if mL_size = 1, then FB is done site-wise, but the
			// output file uses windows of 500 to save space.
			posterior_computer.computePosteriorDistros(temp_file_name, partitions.toArray(),500);
			tempfiles.add(temp_file_name);
			
			// print LL for this stream, and the samples involved
			System.out.print("Posterior computed for samples: "+Arrays.toString(sub_hap_indices[i]) + " :" );
			System.out.println("LL = "+posterior_computer.get_LL());
			
		}
		System.out.println("\n --------------------------------------------------------------------------------------------------------------- \n\n\n");

	    String out_file = arg_box[0].getString("output_file");
	    
		// now combine all the temp files into one single final output file, and delete the temp files
		average_posteriors(tempfiles, out_file);
		
		
		// obtain the modes if asked for
		if(arg_box[0].getBoolean("mode_file")) {  this.obtain_mode_distribution(out_file) ; }
		
		
		dataStreams.deleteStreamFiles();
		Files.delete(temppath);
		System.out.println("Done!!");
	}
	
	
	//////////////////////////////////////////

	// helper methods to process posterior files. 
	
	// averages posteriors across t_files, and then flips columns so position is increasing
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
		
		
		// flip elements so that the positions are correctly ordered
		int l_ind = 1; int u_ind = running_posterior.getColumnDimension()-1; //traverse from both sides, swapping
		RealVector t_vec;
		while(l_ind < u_ind){
			t_vec = running_posterior.getColumnVector(l_ind);
			running_posterior.setColumnVector(l_ind, running_posterior.getColumnVector(u_ind));
			running_posterior.setColumnVector(u_ind, t_vec);
			l_ind++; u_ind--;
		}
		
		
		
		// write to the outfile
		helper.write_doubleMat_ToCSV(o_file, running_posterior.getData());
		
	}

	// creates a new file with info about the mode from the input file.
	void obtain_mode_distribution(String inFile) {
		// takes the posterior in inFile and writes it to an output called <infile>_MODE.csv
		
		BufferedReader br = null; FileWriter scribe = null;
        String line = "";
        String csvSplitBy = ",";

        try {

            br = new BufferedReader(new FileReader(inFile+ ".csv") );
    		scribe = new FileWriter(inFile + "_MODE.csv");
    		
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
	
	
	
	////////////////////
	
	// Args parser for main arguments and hidden arguments
	JSAPResult[] wrap_posterior_args(String[] args, String help_description) throws IOException, JSAPException {
		
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
						new FlaggedOption("base_n_samples", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "base_n",
								"Number of haplotypes considered in each parallel CHMM model (The total number of haplotypes analyzed is <base_n> * <n_groups> under our composite likelihood scheme). We recommend at most choosing --base_n=10 (though for TMRCA the method is still tractable up to 30) and adjusting <n_groups> accordingly to tackle larger sample sizes."),
						

						// FILE PARAMS
						//////////


						// INPUT VCF DATA FILE LIST: --vcf_file
						new FlaggedOption("vcf_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "vcf_file",
								"VCF file for posterior decoding. We assume the VCF is for a single chromosome (or a single continuous segment). Will analyze the first (<base_n> * <n_groups>) haplotypes in the file."),
							
						// REFERENCE DATA FILE LIST: --ref_file
						new FlaggedOption("ref_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ref_file",
								"Reference file for the chromosome."),
						
						// ANCESTRAL DATA FILE LIST: --anc_file
						new FlaggedOption("anc_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "anc_file",
								"Ancestral file for each chromosome. Can specify same file as <ref_file> if you don't have this (though accuracy may be affected if reference alleles and ancestral alleles are not the same)."),
						
						// OUTPUT FILE: --out_file
						new FlaggedOption("output_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "out_file",
								"Prefix for output files containing posterior info. "),
						
						
						
						//////////
						///////// OPTIONAL VISIBLE PARAMETERS ///////////////////
						//////////
						

						// NUMBER OF GROUPS OF SUBSAMPLES: --n_groups
						new FlaggedOption("sub_sample_groups", JSAP.INTEGER_PARSER, "1" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "n_groups",
								"Number of subordinate CHMMs (each considering <base_n> haplotypes) analyzed in composite analysis. Total number of haplotypes analyzed is  <base_n> * <n_groups>. The posterior likelihood file contains the average of the posteriors for each subsampled tree. Such a scheme will partially attenuate spurious dips/spikes in the tree height/length."),
						
						// INITIAL POPULATION SIZE HISTORY GUESS: --psh0
						new FlaggedOption("psh_const", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "psh0",
								"Specify constant population size. Default is to use Waterson's estimate for N_effective, computed from data (and based on the specified mutation rate). Specifying this option will use <psh0> instead of Waterson's estimate in partitioning CHMM states. This value also is assumed to specify the constant population size demography unless <psh_params> or [psh_xs]+[psh_ys] is specified."),
							
						// OUTPUT FILE: --psh_params
						new FlaggedOption("psh_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "psh_params",
								"Provide population size history using .param file produced from CHIMP's demographic inference functionality. This file should contain piecewise constant population size parameters. Providing this option overrides parameters specified using [psh_xs] and [psh_ys]."),
						
						// INITIAL POPULATION SIZE HISTORY DOMAIN: --psh_xs
						new FlaggedOption("psh_xs", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "psh_xs",
								"Use [psh_xs] and [psh_ys] together to specify population size history. [psh_xs] is a list of boundaries between epochs.").setList(true).setListSeparator(','),
							
						// INITIAL POPULATION SIZE HISTORY VALS GUESS: --psh_ys
						new FlaggedOption("psh_ys", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "psh_ys",
								"Use [psh_xs] and [psh_ys] together to specify population size history. [psh_ys] is a list of population sizes. Length must be 1 greater than that of [psh_xs]).").setList(true).setListSeparator(','),
						
						// NUMBER OF STATES IN CHMM: --n_states
						new FlaggedOption("num_chmm_states", JSAP.INTEGER_PARSER, "50" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "n_states",
								"Number of discrete states for CHMM (discrete intervals into which TMRCA or L falls). States are partitioned according to equal probability for each state under a constant population size prior (Waterson's estimate)."),
						
						// METALOCUS IMPLEMENTATION, BEST USED WITH TREELENGTH --metalocus_size
						new FlaggedOption("metalocus_size", JSAP.INTEGER_PARSER, "500" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "metalocus_size",
								"Number of bases grouped into each metalocus. Selecting metalocus_size=1 will activate a locus skipping algorithm to aid computational efficiency."),
						

						// SWITCH TO ALSO OUTPUT MODE OF POSTERIOR: --post_mode
						new Switch("mode_file", JSAP.NO_SHORTFLAG, "post_mode", 
								"Use this to include a file <output_file>_MODE.csv that contains the mode of the posterior."),
						
						// SWITCH TO USE TOTAL TREE LENGTH INSTEAD OF TMRCA: --tree_length
						new Switch("tree_length_switch", 'l', "tree_length", 
								"Use this to specify use of total branch length (L) instead of tree height (TMRCA) as representation of the CHMMs hidden state."),
					
						// SWITCH FOR PSEDUOHAPLOID EMISSION MODEL --pseudo
						new Switch("pseudo_haploid", JSAP.NO_SHORTFLAG, "pseudo", 
								"Use this switch to specify that the data is pseudohaploid data. CHIMP will assume each allele for each variant in the VCF is randomly selected from either haploid of a diploid individual (independently for each variant/individual). The tree height/length in the posterior is computed for all haplotypes of the parents, or 2*<base_n> haplotypes."),
					
						// SWITCH FOR PRINTING HELP FOR HIDDEN PARAMETERS --help_hidden
						new Switch("helpHidden", JSAP.NO_SHORTFLAG, "help_hidden", 
								"Print help descriptions for all the hidden options."),
						
						// FILE FOR SPECIFYING THE HAPLOTYPE SUBGROUPS --hap_groups
						new FlaggedOption("subgroup_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "hap_groups",
								"Use this file to specify the haplotype subgroups used in the composite likelihood framework. Each row contains list of haplotype indices (positive integers, in order following VCF). Row length must match <base_n>. <ngroups> gets overwritten"),
						
										
				}
		
		
		);
		
		/////////////
		//////////// INVISIBLE PARAMETERS
		/////////////
		
		SimpleJSAP jsap_hidden = new SimpleJSAP("",
				"",
				new Parameter[] {
						
						
						/////////////
						//////////// INVISIBLE PARAMETERS
						/////////////

						
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
											
						
						
						// SCRATCH FOLDER FOR DATASTREAM FILES: --disc_storage 
						new Switch("temporary_file_folder", JSAP.NO_SHORTFLAG, "disc_storage",
								"Use option to store intermediate data structures temporarily on the disc (default stores them in RAM). Using this option may decrease RAM usage, but may increase run-time since structures will be stored to and read from disc."),
						
						// CHROMOSOMAL LENGTHS FOR VCFS: --chr_l
						new FlaggedOption("chromosome_length", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "chr_l",
								"Chromosome length for VCF. Will analyze positions [1,length]. If not provided, the reference file will determine the contiguous chromosomal tract that is analyzed."),
						

	
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
				System.err.println("Error: Issue with <" + p_name +">. ");
				fail_Switch=true;
			}
		}
		if(fail_Switch) {System.err.println("Issue Parsing. "); System.exit(0);}
		
		
		
		//modify default values contingently
		return new JSAPResult[] {arguments, hidden_args};
	}

	
	
	
	
}
