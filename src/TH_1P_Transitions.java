import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class TH_1P_Transitions extends TransitionsComputer {

	
	
	TH_1P_Transitions(FunctionUnivariate cr, double rr, int ns, double[] dts) {
		super(cr, rr, ns, dts);
		// TODO Auto-generated constructor stub
	}

	
	
	//////////////////////////////////////////////////////////////////////////
	// INHERITED ABSTRACT METHODS TO IMPLEMENT
	//////////////////////////////////////////////////////////////////////////////
	
	
	@Override
	void updateCDF() {
		
		double[][] res_prob_soi = helper.odeSolsAtTimes(this, getInitialResProbs(), d_times, t_step_precision, false, null);
		
		// first obtain the residence_probability at all of the discrete times for states of interest
		int[] soi = new int[] {stateToInt(new int[] {1,0,0,0}), stateToInt(new int[] {1,0,0,1}), stateToInt(new int[] {1,0,1,1})}; // states of interest
				
		// use the residence_probs to compute the joint CDF (discretized), ie prob that Ta < t1 and Tb < t2
		double[][] disc_CDF = new double[num_t_bins+1][num_t_bins+1]; // + 1 to include infinity bin
		
		for(int i1 = 0 ; i1 < num_t_bins ; i1++) {
			for(int i2 = i1 ; i2 < num_t_bins ; i2++) {
				
				// structure of for loops ensures that t1 <= t2
				double prob; // joint probability of Ta < t1 and Tb < t2
				double fcp; // this is probability of the final coalescence event (in t2-t1 time)
							
				fcp = singlePairCoalescence(d_times[i1], d_times[i2]);
				prob = res_prob_soi[soi[0]][i1] + res_prob_soi[soi[1]][i1] + res_prob_soi[soi[2]][i1] *fcp;
						
				disc_CDF[i1][i2] = prob;
				disc_CDF[i2][i1] = prob;
			}
					
			// compute jointCDF in the last row and column (where one of the variables is < infinity, so really this isn't a joint cdf)
			// ie, sum up residence probs for states where coalescence at A has occured before t1 
			// (note  1010, 0111, 0110, 0101, 0100 are not possible)
			double t_a_cdf = res_prob_soi[soi[0]][i1] + res_prob_soi[soi[1]][i1] + res_prob_soi[soi[2]][i1] ;
			disc_CDF[i1][num_t_bins] = t_a_cdf;
			disc_CDF[num_t_bins][i1] = t_a_cdf;
					
					
		}
				
		disc_CDF[num_t_bins][num_t_bins] = 1; // final entry, that both mrca times are less than infinity
		
		currentCDF = new Array2DRowRealMatrix(disc_CDF);
		upToDate = true;

	}

	
	@Override
	double[] ancestralSim() {
		
		// specifically for recombination restricted to 1
		
		// setup the tree_heights to return, at each of 2 loci
		double th1 = 0.; 
		double th2 = 0.;
		boolean[] th_found = new boolean[] {false, false};
						
		
		// initialize starting state (n,0,0,0) for Markov process, the start time (0) and absorbing states
		int[] state = new int[] {n_samples,0,0,0};
		double sim_time = 0. ;
				
					 	
		// begin markov process
		while(  !(th_found[0] && th_found[1])   ) {
					
			double[] rates = new double[5];
					
			rates[0] = helper.nChoose2(state[0]); // ab - ab coalescence
			rates[1] = state[0] * state[1]; // ab - a coalescence 
			rates[2] = state[0] * state[2]; // ab - b coalescence
			rates[3] = state[1] * state[2]; // a - b coalescence
			if(state[3] == 0) { rates[4] = state[0]; } // ab recombination
					
			// compute the total leaving rate from state as a function of time
			FunctionUnivariate total_leaving_rate = coalescence_rate.copy();
			total_leaving_rate.multiplyScalarSelf(rates[0] + rates[1] + rates[2] + rates[3] );
			total_leaving_rate.addScalarSelf(rates[4] * recombination_rate/2.);
							
			// randomly pick from appropriate part of distribution and obtain time of next event
							
			/////
			sim_time = total_leaving_rate.inverseIntegrate0(- Math.log(  1. -  Math.random()  ) +  total_leaving_rate.integrate0(sim_time)  );
			////////
							
			double cr = coalescence_rate.evaluate(sim_time);
			RealVector v_rates = new ArrayRealVector(rates);
			RealVector rs = new ArrayRealVector(new double[] {cr, cr, cr, cr, recombination_rate/2.});
			v_rates = v_rates.ebeMultiply(rs);
			v_rates.mapDivideToSelf(  v_rates.getL1Norm()  );
					
			int event = helper.roll(v_rates);
					
					
					
			// handle each event
			if(event == 0) { 
				state[0] = state[0] - 1 ; }
					
			if(event == 1) { 
				state[1] = state[1] - 1 ; }
					
			if(event == 2) { 
				state[2] = state[2] - 1 ;}
					
			if(event == 3) { 
				state[0] = state[0] + 1 ; 
				state[1] = state[1] - 1; 
				state[2] = state[2] - 1;}
					
			if(event == 4) { 
				state[0] = state[0] - 1;
				state[1] = state[1] + 1;
				state[2] = state[2] + 1;
				state[3] = state[3] + 1;}
					
					
					
			// terminate loop and save TL values as appropriate
			if(!th_found[0] && (state[0]+state[1] == 1) ) {
				th_found[0] = true;
				th1 = sim_time;
			}
			if(!th_found[1] && (state[0]+state[2] == 1) ) {
				th_found[1] = true;
				th2 = sim_time;
			}
							
							
		}
						
		return new double[] {th1 , th2}; 
	}

	
	double[] getInitialResProbs() {
		
		double[] out = new double[num_states];
		out[stateToInt(new int[] {n_samples,0,0,0})] = 1.;
		return out;
	}

	
	///////////////////////////////////////////
	//access 
	@Override
	void updateDemoFromPs(double[][] params) {

		coalescence_rate.updateFromParams(params[0]);
		recombination_rate = params[1][0];
		
		upToDate = false;
	}

	@Override
	double[][] getPsFromDemo() {
		double[][] out = new double[3][];
		
		out[0] = coalescence_rate.getFunctionParameters();
		out[1] = new double[] {recombination_rate};
		
		return out;
	}

	
	
	
	/////////////////////////////////////////////////////////////////
	// HELPER FUNCTIONS //////////
	/////////////////////////////////////////////////////////////////
	
	// THIS COMPUTES PROBABILITY OF COALESCENCE WITHIN T2-T1 ELAPSING FROM START AT T1, UNDER RATE CRATE
	private double singlePairCoalescence(double t1, double t2) {
		if(t1 > t2) {
			System.out.println("interval bounds for time not ascending in order");
			System.exit(0);
		}
		return  1 - Math.exp(-( coalescence_rate.integrate0(t2) - coalescence_rate.integrate0(t1)));
	}



	

	//////////////////////////////////////////////////////////
	// function to compare to sims to evaluate how well we're doing
	///////////////////////////////////////////////////////////
	
	private void compareSim(int n) {
		
		RealMatrix cm = this.getPDF();
		RealMatrix sm = this.cdfToPDF(this.simulatedCDF(n));
		
		//RealMatrix cm = this.getCDF();
		//RealMatrix sm = this.simulatedCDF(n);
		
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
		int n_samples = 10;
		
		// initial value for psh
		double psh0 = 11559; // to define states
		
		/////// //////////////////////////////////////
		double[] pp = new ArrayRealVector(10,.1).toArray();
		

		RealVector dlengths = helper.getTreeHeightPartitions(pp, n_samples);
		dlengths.mapMultiplyToSelf(psh0 * 2);
			
		
		double sfactor = 1000;
		
		double [] xs = new double[] {57.0, 113.653, 226.614, 451.849, 900.95, 1796.417, 3581.90, 7142.002, 14240.539, 28394.411, 56616.017, 112887.476, 225087.931, 448806.0};
		//double[] ys = new double[] {17000,17000,17000,17000,17000,17000,17000,17000,17000,17000,17000,17000}; 
		double[] ys = new double[] {50000, 15811, 5000, 15811, 50000, 15811, 5000, 15811,50000,15811,5000,15811,50000,15811,5000}; 

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		// scale all demographic params and get ready to be input into EM_HMM

		psh0 = sfactor / (2. * psh0);
		Function_PWC cr = new Function_PWC(xs, ys, false);
		cr.reciprocalSelf();
		cr.multiplyScalarSelf(sfactor / 2.0);
		cr.dilateXSelf(1./sfactor);
	
		double theta = 2 * sfactor * mut_rate;
		double rho = 2 * sfactor * rec_rate; //rho = -rho; // setting it less than zero forces SFS scheme
		double[] scaled_lengths = new ArrayRealVector(dlengths).mapDivide(sfactor).toArray(); // these are in coalescent time = 2*N_ref generations
		
		
		TH_1P_Transitions tp_comp = new TH_1P_Transitions(cr,rho, n_samples, scaled_lengths); 

		
		
		
		System.out.println("********** \n \n STATE PARTITIONS:\n");
		System.out.println(dlengths);
		System.out.println("\n \n **********");
		System.out.println("************************");
		System.out.println("NUMERIC JOINT CDF:");
		System.out.println(tp_comp.getCDF());
		System.out.println("************************");
		System.out.println("\n \n NUMERIC JOINT PDF:");
		System.out.println(tp_comp.getPDF());
		System.out.println("************************");
		System.out.println("\n \n NUMERIC TRANSITION PROBS:");
		System.out.println(tp_comp.transitionProbs().power(1));
		System.out.println("\n \n \n ************************");
		//System.out.println("SIMULATION COMPARISON (10^8 RUNS) ");
		//tp_comp.compareSim(100000000);
		
		System.out.println("************************");

		System.out.println("************************");
		RealMatrix trap = tp_comp.transitionProbs();
		System.out.println(trap.power(500));

	}



	

	
	
	
	
}
