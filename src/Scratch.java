import java.util.Arrays;

import org.apache.commons.math3.linear.*;


public class Scratch {
	
	/*
	 * THIS CLASS IS EXTRANEOUS. 
	 * WAS JUST USED AS A PLACE TO WRITE A SCRIPT TO COMPUTE THE TRANSITION AND EMISSION PROBABILITIES
	 * FROM THE TREE SEQUENCE OUTPUT FILES FROM MSPRIME SIMULATION, TO COMPARE AGAINST OUR METHODS'
	 * SEEMED TO MATCH UP FAIRLY WELL. 
	 * 
	 */
	
	
	
	

	String folder = "/Users/gautam/Desktop/scratch/data_cons/";
	String[] files = new String[] {"con_dataset1", "con_dataset2", "con_dataset3", "con_dataset4", "con_dataset5","con_dataset6", "con_dataset7", "con_dataset8", "con_dataset9", "con_dataset10","con_dataset11", "con_dataset12", "con_dataset13", "con_dataset14", "con_dataset15","con_dataset16", "con_dataset17", "con_dataset18", "con_dataset19", "con_dataset20","con_dataset21", "con_dataset22", "con_dataset23", "con_dataset24", "con_dataset25","con_dataset26", "con_dataset27", "con_dataset28", "con_dataset29", "con_dataset30"};
	//String[] files = new String[] {"bot_dataset1", "bot_dataset2", "bot_dataset3", "bot_dataset4", "bot_dataset5","bot_dataset6", "bot_dataset7", "bot_dataset8", "bot_dataset9", "bot_dataset10","bot_dataset11", "bot_dataset12", "bot_dataset13", "bot_dataset14", "bot_dataset15","bot_dataset16", "bot_dataset17", "bot_dataset18", "bot_dataset19", "bot_dataset20","bot_dataset21", "bot_dataset22", "bot_dataset23", "bot_dataset24", "bot_dataset25","bot_dataset26", "bot_dataset27", "bot_dataset28", "bot_dataset29", "bot_dataset30"};

	
	//double[] partitions = new double[] {18936, 26578, 35648, 49998}; // for const 10K size
	double[] partitions = new double[] {72380, 93420, 115900, 148400};
	
	int n_samples = 10;
	int base_pairs = 200000000; // per vcf file
	
	
	public static void main(String[] args){
		Scratch ss  = new Scratch();
		
		int tloth = 1; // 0 for tree height, 1 for tree length
		
		ss.ProbDists(tloth);
		ss.analyticDists();

		System.out.println("DONE!");
		
	}
	
	
	public void ProbDists(int tloth){
		Helper helper = new Helper();
		
		RealMatrix transEvents = new Array2DRowRealMatrix(partitions.length+1 , partitions.length+1);
		RealMatrix emitEvents = new Array2DRowRealMatrix(partitions.length+1 , n_samples);
		RealVector margEvents = new ArrayRealVector(partitions.length + 1);
		
		
		for(int ff = 0 ; ff < files.length ; ff++) {
			String file_pre = folder + files[ff];
			
			double[][] tseq = helper.read_doubleMat_FromCSV(file_pre+"_TS");
			VcfReader vread = new VcfReader(file_pre+".vcf", folder + "sim_ref.fasta", folder + "sim_anc.fasta", base_pairs, n_samples);
			int[][] datastream = vread.getReducedStream();
			
			// have tseq and datastream ready to do tree analysis
			
			// get current tree and mmtract, with index and length
			double[] tree = new double[] {0, tseq[2][0]-1}; 
			double[] mmtract = new double[] {0, datastream[1][0]};
			
			int state = (int) helper.placeInPartition(tseq[tloth][0], partitions);
			int prev_state = state;
			int n_da = datastream[0][0]; //number of derived alleles at head of current mms tract
			
			int genome_pos = 1 ;
			
			while(genome_pos < base_pairs + 1) {
				
				// if tracts need to be reset, then reset
				
				if(tree[1] < 0) {
					double temp = tree[1];
					tree[0] = tree[0] + 1 ;
					tree[1] = tseq[2][(int) tree[0]] + temp;
					state = (int) helper.placeInPartition(tseq[tloth][(int) tree[0]], partitions);
				}
				
				if(mmtract[1] < 0.5) { // another way of asking if the int is equal to zero
					mmtract[0] = mmtract[0] + 1;
					mmtract[1] = datastream[1][(int) mmtract[0]];
					n_da = datastream[0][(int) mmtract[0]];
				}
				
				

				// handle transition of tree
				transEvents.setEntry(prev_state, state, 1 + transEvents.getEntry(prev_state, state));
				
				// handle emission of tree
				emitEvents.setEntry(state, n_da, 1 + emitEvents.getEntry(state, n_da));
				
				// handle marginal state
				margEvents.setEntry(state, margEvents.getEntry(state) + 1);
				
				genome_pos++;
				prev_state = state;
				if(n_da > 0) {n_da = 0;} // set to zero after we pass the head of the tract
				
				// shorten tseq and mm tracts
				tree[1] = tree[1] - 1 ;
				mmtract[1] = mmtract[1] - 1 ;
				
			}
			
		}
		
		transEvents = transEvents.scalarMultiply(1. / helper.totalEntriesMat(transEvents));
		emitEvents = emitEvents.scalarMultiply(1. / helper.totalEntriesMat(emitEvents));
		margEvents.mapDivideToSelf(margEvents.getL1Norm());
		
		
		System.out.println("Transition Probabilities are: \n");
		System.out.println(transEvents);
		System.out.println("\n \n \n");
		
		System.out.println("Emission Probabilities are: \n");
		System.out.println(emitEvents);
		System.out.println("\n \n \n");
		
		System.out.println("Marginal Probabilities are: \n");
		System.out.println(margEvents);
		System.out.println("\n \n \n");
		
		
		
	}
	
	public void analyticDists() {
		Helper helper = new Helper();
		
		double rec_rate =  0.00000001;
		double mut_rate = 0.00000001;
		
		// initial value for psh
		double psh0 = 10000;
		RealVector dlengths = new ArrayRealVector(partitions);

		
		double sfactor = 1000;
		
		double [] xs = new double[] {2000,4000};
		double[] ys = new double[] {10000,10000,10000};

		/////////////////////////////////////////////////////////////////////////////////////////////
		
		// scale all demographic params and get ready to be input into EM_HMM

		psh0 = sfactor / (2. * psh0);
		Function_PWC cr = new Function_PWC(xs, ys, false);
		cr.reciprocalSelf();
		cr.multiplyScalarSelf(sfactor / 2.0);
		cr.dilateXSelf(1./sfactor);
	
		double theta = 2 * sfactor * mut_rate;
		double rho = 2 * sfactor * rec_rate;
		double[] scaled_lengths = new ArrayRealVector(dlengths).mapDivide(sfactor).toArray(); // these are in coalescent time = 2*N_ref generations
		
		
		TL_1P_Transitions tp_comp = new TL_1P_Transitions(cr, rho, n_samples, scaled_lengths, new int[] {200,200});
		TL_1P_Emissions ep_comp = new TL_1P_Emissions(cr,theta, n_samples, scaled_lengths, 1000);
		
		System.out.println("*****************************************************");
		System.out.println("***************************************************** \n \n");

		System.out.println("Analytic Results: \n");
		
		System.out.println("Transition PDF: \n");
		System.out.println(tp_comp.getPDF());
		
		System.out.println("\n \n Emission PDF: \n");
		System.out.println(ep_comp.getPDF());
		
		System.out.println("\n \n Marg PDF: \n");
		System.out.println(ep_comp.getMargPDF());
	
		
	}
	
	
}
