import java.util.Arrays;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class TH_1P_Emissions extends EmissionsComputer {

	RealMatrix combinatoricMatrix;
	
	TH_1P_Emissions(FunctionUnivariate cr, double rr, int ns, double[] dts, boolean pseudo_hap) {
		super(cr, rr, ns, dts, pseudo_hap);
		// TODO Auto-generated constructor stub
		
		num_emissions = n_samples;
		combinatoricMatrix = computeMEMatrix(n_samples);
	}

	
	
	
	
	
	//////////////////////////////////////////////////////////////////////////
	// INHERITED ABSTRACT METHODS TO IMPLEMENT
	//////////////////////////////////////////////////////////////////////////////
	
	
	@Override
	void updateCDF() {


//		d_times[d_times.length-1] = 584.7962407069716;
//		System.out.println(Arrays.toString(d_times));
		
		
		double[][] res_prob = helper.odeSolsAtTimes(this, this.getInitialResProbs(), d_times, t_step_precision, true, this.getAbsorbingStateIndices()) ;
		
		// hold onto residence probabilities only of absorbing states (where we have found tmrca)
		
		RealMatrix res_prob_absorbing = new Array2DRowRealMatrix(n_samples, num_t_bins + 1);
		
		res_prob_absorbing.setRow(0, res_prob[stateToInt(new int[] {1,0})]);
		for(int i = 1 ; i < res_prob_absorbing.getRowDimension() ; i++) {
			// note since 1,1 is dummy state, row 1 corresp to (1,2) ; row 2 to (1,3) etc //
			res_prob_absorbing.setRow(i, res_prob[stateToInt(new int[] {1,i+1})]);
		}

		res_prob_absorbing = res_prob_absorbing.transpose();
		
		// the final CDF is obtained by multiplying these res probs by the combinatoric emission factors, and summing over the absorbing states
		currentCDF = res_prob_absorbing.multiply(combinatoricMatrix);
		upToDate = true;
		
	}

	
	@Override
	int[] ancestralSim() {
		
		// setup the tree_lengths to return, at each of 2 loci
		int seg_sites = 0;
				
				
		// initialize the starting representation of state, where all lineages subtend one modern sample
		int[] state = new int[n_samples];  // how many samples subtended by extant lineages. the extant lineages are first "lineages" elements of state
		Arrays.fill(state, 1);
		
		int lineages = n_samples;          // number of ancestral lineages extant
		double sim_time = 0.;			// simulation time
		
		boolean mutated = false;			// marker to keep track if one mutation occurred or not
		Random ind_gen = new Random();
		
		
		while(lineages > 1) {
			
			// compute rate of events occurring (coalescence or mutation)
			FunctionUnivariate event_rate = coalescence_rate.copy();
			event_rate.multiplyScalarSelf(helper.nChoose2(lineages));
			if(!mutated) { event_rate.addScalarSelf( lineages * mutation_rate/2. );}  // divide by 2 by definition of m_rate
			
			//////// compute time of next event and set that as the simulation time
			sim_time = event_rate.inverseIntegrate0(- Math.log(  1. -  Math.random()  ) +  event_rate.integrate0(sim_time)  );
			////////////////////
			
			//// now determine which event occurred and update state parameters accordingly.
			
			// mutation can occur only if it has not occurred yet, and otherwise with some probability based on rates at sim_time
			if( !mutated && Math.random() < (lineages * mutation_rate*0.5 / event_rate.evaluate(sim_time))) {
				seg_sites = state[ ind_gen.nextInt(lineages) ]; // write samples subtended into holder (output)
				mutated = true;
			}
			
			// in all other cases, it will be a coalescence event.
			else {
				
				// update state 
				int c_l = ind_gen.nextInt(lineages - 1);
				state[c_l] = state[c_l] + state[lineages - 1];    // rightmost extant lineage feeds into another random lineage					
				
				// update times/lineages
				lineages = lineages - 1;  // coalescence will decrease number of lineages by 1.
			}
			
		}
		
		int t_index = helper.placeInPartition(sim_time, d_times);
		
		return new int[] {t_index, seg_sites}; 
		
	}

	

	///////////////////////////////////////////
	//access 
	@Override
	void updateDemoFromPs(double[][] params) {

		coalescence_rate.updateFromParams(params[0]);
		mutation_rate = params[2][0];
		
		upToDate = false;
	}

	@Override
	double[][] getPsFromDemo() {
		double[][] out = new double[3][];
		
		out[0] = coalescence_rate.getFunctionParameters();
		out[2] = new double[] {mutation_rate};
		
		return out;
	}

	
	
	/////////////////////////////////////////////////////////////////
	// HELPER FUNCTIONS //////////
	/////////////////////////////////////////////////////////////////
		
	private RealMatrix computeMEMatrix(int n){
		// n is number of samples
		// computes a probability matrix. first index is number of extant lineages at time of mutation. 
		// second index is output (number of differentiating alleles), ie number of samples subtended by mutation
		// values of mem hold probability of observing an output given that mutation occurred with so many lineages
		
		Helper.ChooseComputer chooser = helper.new ChooseComputer(n_samples);
		
		RealMatrix mem = new Array2DRowRealMatrix(n,n);
		
		// for s=1 is dummy state (not possible), and s=0 always emits 0
		mem.setEntry(0, 0, 1.);
		
		// starts with ks = 2, but this corresponds to row 1 since 1,1 is ignored as dummy, and row 1 is for ks 2 instead
		for(int s = 2 ; s < n+1 ; s++) {
			RealVector sEmits = new ArrayRealVector(n); // use this to store the probabilities for this row of mem
			
			for(int i = 1 ; i < n ; i++) {
				// if mutated element size of vector is i, then the number of vectors formed by picking other element sizes is 
				// given by choose(n-i-1 , s-2) (ball and partition method).
				sEmits.setEntry(i, chooser.get(n-i-1, s-2)); //*************
			}
			// note our choose function returns 0 if i is too large (ie if there are not enough possible partition spots)
			
			sEmits.mapMultiplyToSelf(1.0 / sEmits.getL1Norm()); // normalize over all possible vector configs
			mem.setRowVector(s-1, sEmits);  // write into appropriate row of mem (ks - 1 is row)
		}
		return mem;
	}
	
	
	
	
	
	
	
	
	
	//////////////////////////////////////////////////////////
	// function to compare to sims to evaluate how well we're doing
	///////////////////////////////////////////////////////////
	
	public void compareSim(int n) {
		//RealMatrix cm = this.getCDF();
		//RealMatrix sm = this.simulatedCDF(n);
		
		RealMatrix cm = this.getPDF();
		RealMatrix sm = this.cdfToPDF(this.simulatedCDF(n));
		
		String barrier = "*************************************************************************************************";
		
		System.out.println( "\n" +barrier + "\n" + barrier);
		System.out.println("Computed Matrix is:");
		System.out.println(cm);
		
		System.out.println( "\n" +barrier + "\n" + barrier);
		System.out.println("Simulated Matrix (n = " + n + ") is:");
		System.out.println(sm);
		
		System.out.println( "\n" +barrier + "\n" + barrier);
		System.out.println("StdErr Comparison is:");
		System.out.println(helper.stdErrorMatComp(sm, cm, n));
		
		System.out.println( "\n" +barrier + "\n" + barrier);

	}
	
	
	
//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

	public static void main(String[] args) {
		Helper helper = new Helper();
		
		double rec_rate = 0.0000000125;
		double mut_rate = 0.0000000125;
		int n_samples = 5;
		
		// initial value for psh
		double psh0 = 17000; // to define states
				
		/////// //////////////////////////////////////
		//double[] pp = new double[] {.001,.001,.002,.002,.002,.002,.002,.002,.002,.002,.002,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.1};
		//double[] pp = new double[] {.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02,.02};
		double[] pp = new double[] {.1,.1,.1,.1,.1,.1,.1,.1,.1,.1};

		
		RealVector dlengths = helper.getTreeHeightPartitions(pp, n_samples);
		dlengths.mapMultiplyToSelf(psh0 * 2);
					
				
		double sfactor = 1000;
				
		double [] xs = helper.generate_log_discretization(56.501,56501, 12);
		double[] ys = new double[] {17000,17000,17000,17000,17000,17000,17000,17000,17000,17000,17000,17000}; 
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
		
		
		TH_1P_Emissions ep_comp = new TH_1P_Emissions(cr,theta, n_samples, scaled_lengths, false); 
		//System.out.println(ep_comp.pseudo_hap_filter);
		//System.out.println(helper.mat_rowsums(ep_comp.pseudo_hap_filter));
		System.out.println(ep_comp.emissionProbs(false));
		System.out.println(helper.mat_rowsums(ep_comp.emissionProbs(false)));

		
		/*	
		System.out.println("********** \n \n STATE PARTITIONS:\n");
		System.out.println(dlengths);
		System.out.println("\n \n **********");
		System.out.println("************************");
		System.out.println("MARGINAL PDF:");
		System.out.println(ep_comp.getMargPDF());
		System.out.println("************************");
		System.out.println("SFS:");
		System.out.println(ep_comp.getSFS());
		System.out.println("\n \n NUMERIC JOINT PDF:");
		System.out.println(ep_comp.getPDF());
		System.out.println("\n \n NUMERIC EMISSION PROBS:");
		System.out.println(helper.pdfToRowProbMatrix(ep_comp.getPDF(), 1));
		System.out.println("\n \n \n ************************");
		*/
		//System.out.println("SIMULATION COMPARISON (10^8 RUNS) ");
		//ep_comp.compareSim(100000000);
		
		System.out.println("************************");

		System.out.println("************************");

	}

}
