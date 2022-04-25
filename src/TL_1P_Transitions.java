import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;


public class TL_1P_Transitions extends TransitionsComputer {
	
	// Switch to use Linear Interpolations Instead of Spline Interp
	boolean switch_linInt = true;
	
	// coalescence rate function discontinuity points
	double[] cr_discontinuities;
	
	// Parameters of class that change with population size history, and record of the most recently computed values. 
	// use these to store the final stage of the most recent evaluation. can easily compute discrete CDF matrix from this.
	Function_Spline f1;
	Function_BiSpline f1001_vals;
	Function_Spline f1000_vals;
	RealVector lambda_grid;
	double[][] ode_grid_sol;

	
	// v vectors defining the characteristic directions for each state
	int[] cVecA;
	int[] cVecB;
	
	// params regarding the internal discretization
	int num_disc_times;
	RealVector t_grid;	
	
	
	// CONSTRUCTOR AND HELPING FUNCTIONS
	TL_1P_Transitions(FunctionUnivariate cr, double rr, int ns, double[] dts, int[] gsize) {
		super(cr, rr, ns, dts);
		

		// populate the characteristic vector directions
		buildVchars();
		
		createTimeGrid_uniformInsert(gsize[0], gsize[1]);

	}
	
	TL_1P_Transitions(FunctionUnivariate cr, double rr, int ns, double[] dts, int gsize) {
		super(cr, rr, ns, dts);
		
		// split so that equal number of divisions in each regime.
		int[] gss = new int[2];
		gss[0] = (int) gsize / 2;
		gss[1] = gsize - gss[0];
		
		// populate the characteristic vector directions
		buildVchars();
		
		createTimeGrid_uniformInsert(gss[0], gss[1]);
	}


	private void buildVchars() {
		cVecA = new int[num_states];
		cVecB = new int[num_states];
		
		for(int s = 0 ; s < num_states ; s++) {
			int[] state = intToState(s);
			int aLins = state[0] + state[1];
			int bLins = state[0] + state[2];
					
			if( aLins > 1) {cVecA[s] = aLins ;}
			if( bLins > 1) {cVecB[s] = bLins ;}

			
		}
		
	}
	

	private void createTimeGrid_uniformInsert(int num_a, int num_b) {
		
		num_disc_times = num_a + num_b + 1;
		double[] ts = new double[num_disc_times];
		
		double bin_delta = (d_times[d_times.length-1]/n_samples) /num_a;
		
		// FIRST DEAL WITH THE GRID IN THE DENSE REGION
		for(int grid_pos = 1 ; grid_pos < num_a + 1; grid_pos++) {
				ts[grid_pos] = bin_delta * grid_pos;			
		}
		
		// NOW DEAL WITH GRID IN THE SUPPORT REGION
		double lb = ts[num_a];
		double interval = d_times[num_t_bins-1]/2. - lb;
		double delta = interval / num_b;
		
		for(int b = 1 ; b < num_b + 1 ; b++) {
			int grid_pos = num_a + b;
			ts[grid_pos] = lb + delta * (b);
		}
		
		
		
		// have the uniform grid now. but we want to incorporate the d_times into this
		// insert d_times values in, replacing the values they are nearest to.
		
		for(int i = 0 ; i < num_t_bins ; i++) {
			double dtime = d_times[i] / n_samples;
			
			double prev_distance = 10000000.;
			for(int t = 0 ; t < num_disc_times; t++) {
				double distance = Math.abs(ts[t] - dtime);
				if(distance > prev_distance) {
					int target_i = Math.max(t - 1, 0);
					ts[target_i] = dtime;
					t = num_disc_times;
				}
				prev_distance = distance;
		
			}
			
		}
		

		t_grid = new ArrayRealVector(ts);
	}
		
	///////////////////////////////////////////////////////////////////////////
	// INHERITED ABSTRACT METHODS. MUST IMPLEMENT
	////////////////////////////////////////////////////////////////////////////
	
	
	@Override
	void updateCDF() {
		// only go through trouble of recomputing everything and solving pde if we need to, otherwise just pull values. 
		if (!upToDate) { computeFinalFMats(); }
				
		double[][] cdf = new double[num_t_bins + 1][num_t_bins + 1];
				
		for(int i = 0 ; i < num_t_bins ; i++) {
			for(int j = 0 ; j < num_t_bins ; j++) {
				cdf[i][j] = evaluate_F100X(d_times[i], d_times[j]);
			}
			double bv = f1.evaluate_0(d_times[i]);
			cdf[i][num_t_bins] = bv;
			cdf[num_t_bins][i] = bv;
		}
				
		cdf[num_t_bins][num_t_bins] = 1. ;
				
		this.currentCDF = new Array2DRowRealMatrix(cdf) ;
		adjustCDF();

		upToDate = true;
		
	}

	void adjustCDF() {
		// check if you need to correct CDF due to entries > 1
		boolean needs_correction = false; // matrix probabilities are ill behaved
		boolean uncorrectable = false; // probs are ill behaved in an uncorrectable way
		for(int i = num_t_bins; i > -1 ; i--) {
			for(int j = num_t_bins; j > -1 ; j-- ) {
				double t_ent = currentCDF.getEntry(i, j);
				if(t_ent > 1.) { needs_correction = true; }
				if(t_ent < -.00000001 || t_ent > 2) { uncorrectable = true; }
			}
		}
		if(uncorrectable) {System.err.println("ISSUES COMPUTING TRANSITION PROB, CANNOT CORRECT"); }
		
		// if no correction needed, exit
		if(!needs_correction) {return;}
	
		// correct CDF
		
		RealMatrix t_pdf = this.cdfToPDF(currentCDF);
		RealMatrix t_pmat = helper.pdfToRowProbMatrix(t_pdf, 0);
		
		// find probability of staying in same state for 2 next to last states
		
		double stay_prob1 = t_pmat.getEntry(num_t_bins - 2, num_t_bins - 2);
		double stay_prob2 = t_pmat.getEntry(num_t_bins - 1, num_t_bins - 1);
		double stay_prob_corner = 2*stay_prob2 - stay_prob1;
				
		// project same prob into last state, and find corner val of pdf
		RealVector t_vec = t_pdf.getColumnVector(num_t_bins).getSubVector(0, num_t_bins);
		double corner_val = t_vec.getL1Norm() * (stay_prob_corner) / (1. - stay_prob_corner);

		// find what the corner val in cdf should be based on this.
		corner_val = corner_val + 2*currentCDF.getEntry(num_t_bins, num_t_bins - 1) - currentCDF.getEntry(num_t_bins-1, num_t_bins-1);
		
		System.out.println("\n Warning, numerical drift, cdf adjustment needed: " + corner_val);
		
		currentCDF.setEntry(num_t_bins, num_t_bins, corner_val);
		currentCDF = currentCDF.scalarMultiply(1. / corner_val);
		
	}
	
	
	
	@Override
	double[] ancestralSim() {

		// specifically for recombination restricted to 1
		
		// setup the tree_lengths to return, at each of 2 loci
		double tl1 = 0.; 
		double tl2 = 0.;
		boolean[] tl_found = new boolean[] {false, false};
						
		// intermediate holders, time of last coalescence
		double tlc1 = 0;
		double tlc2 = 0;
				
		// initialize starting state (n,0,0,0) for Markov process, the start time (0) and absorbing states
		int[] state = new int[] {n_samples,0,0,0};
		double sim_time = 0. ;
				
						 	
		// begin markov process
		while(  !(tl_found[0] && tl_found[1])   ) {
					
			double[] rates =new double[5];
					
			rates[0] = helper.nChoose2(state[0]); // ab - ab coalescence
			rates[1] = state[0] * state[1]; // ab - a coalescence 
			rates[2] = state[0] * state[2]; // ab - b coalescence
			rates[3] = state[1] * state[2]; // a - b coalescence
			if(state[3] == 0) { rates[4] = state[0]; } // ab recombination
			//rates[4] = state[0]; // this line makes it a multiple recombination model		
			
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
				tl1 = tl1 + (state[0] + state[1]) * (sim_time - tlc1) ;
				tlc1 = sim_time;
				tl2 = tl2 + (state[0] + state[2]) * (sim_time - tlc2) ;
				tlc2 = sim_time;
				state[0] = state[0] - 1 ; }
					
			if(event == 1) { 
				tl1 = tl1 + (state[0] + state[1]) * (sim_time - tlc1) ;
				tlc1 = sim_time;
				state[1] = state[1] - 1 ; }
			
			if(event == 2) { 
				tl2 = tl2 + (state[0] + state[2]) * (sim_time - tlc2) ;
				tlc2 = sim_time;
				state[2] = state[2] - 1 ; }
					
			if(event == 3) { state[0] = state[0] + 1 ; 
				state[1] = state[1] - 1; 
				state[2] = state[2] - 1; }
					
			if(event == 4) { state[0] = state[0] - 1;
				state[1] = state[1] + 1;
				state[2] = state[2] + 1;
				state[3] = state[3] + 1; }
					
					
					
			// terminate loop and save TL values as appropriate
			if(!tl_found[0] && (state[0]+state[1] == 1) ) {
				tl_found[0] = true;
				//tl1 = sim_time;
			}
			if(!tl_found[1] && (state[0]+state[2] == 1) ) {
				tl_found[1] = true;
				//tl2 = sim_time;
			}
							
							
		}
						
		return new double[] {tl1 , tl2}; 		
		
	}

	@Override
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

	
	
	
	
	
	
	
	
	
	//////////////////////////////////////
	// COMPARE TO SIMS /////
	////////////////////
	

	void compareSim(int n) {
		
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
	
	
	
///////////////////////////////////////////////////////////////////////////
	
	/// PRIVATE METHODS FOR OBTAINING CDF. MEAT OF THIS CLASS ////
	
	// compute sum of f_1000(x,y) and f_1001(x,y) given that f1001 and f1000 are already stored and up to date
	private double evaluate_F100X(double x, double y) {
		double out = 0;
		out = out + f1000_vals.evaluate_0(Math.min(x, y));
		out = out + f1001_vals.value0(x, y);
		return out;
	}
	
	private void createLambdaGrid() {
		double[] lambdas = new double[num_disc_times];
		
		for(int grid_pos = 1 ; grid_pos < num_disc_times ; grid_pos++) {
			lambdas[grid_pos] = coalescence_rate.integrate0(t_grid.getEntry(grid_pos));
		}

		lambda_grid = new ArrayRealVector(lambdas);
	}
	
	private void f_obtainOrder(ArrayList<Integer> f_key, ArrayList< F_s > f_probability_functions) {
		
		// this is the portion to try optimal ordering to reduce heap usage
		int s;
		
		// first state
		
		for(int i = n_samples - 1 ; i > 0 ; i-- ) {
			s = stateToInt( new int[] {i,0,0,0} );
			f_probability_functions.add(new F_s(s));
			f_key.add(s);
			
			s = stateToInt( new int[] {i,1,1,1} );
			f_probability_functions.add(new F_s(s));
			f_key.add(s);
			
			s = stateToInt( new int[] {i+1,0,0,1} );
			f_probability_functions.add(new F_s(s));
			f_key.add(s);
			
			s = stateToInt( new int[] {i,0,1,1} );
			f_probability_functions.add(new F_s(s));
			f_key.add(s);
			
		}
		
		s = stateToInt( new int[] {1,0,0,1} );
		f_probability_functions.add(new F_s(s));
		f_key.add(s);
		
		int singular_index = f_key.indexOf(stateToInt(new int[] {n_samples - 1 , 1, 1, 1}));
		f_probability_functions.set(singular_index, null);
		f_probability_functions.remove( singular_index);
		f_key.remove(singular_index);
		
		
	}
	
	
	
	private void computeFinalFMats() {
		
		// get now rather than having to make a copy every time we integrate on a char line
		cr_discontinuities = coalescence_rate.getDiscontinuities();
				
		// make sure the lambda grid is up to date (since it also depends on coalescence rates).
		createLambdaGrid();
		
		// integrate the ODE to obtain the residence probabilities at the times of interest, for all states
		ode_grid_sol = helper.odeSolsAtTimes(this, getInitialResProbs(), t_grid.toArray(), t_step_precision, false, null);
		
		//////////////
		
		
		// initialize the f function array containing an f object for each state
		ArrayList< F_s > f_probability_functions = new ArrayList< F_s > (num_states);
		ArrayList<Integer> f_key = new ArrayList<Integer> (num_states);
		
		// now fill all the state objects in order which we should cycle through list in some ordering
		f_obtainOrder(f_key, f_probability_functions);
		
		
		// get references for final state objects, easier to access them this way
		F_s f1000 = f_probability_functions.get(f_key.indexOf( stateToInt(new int[] {1,0,0,0}) ));
		F_s f1001 = f_probability_functions.get(f_key.indexOf( stateToInt(new int[] {1,0,0,1}) ));
		F_s f1011 = f_probability_functions.get(f_key.indexOf( stateToInt(new int[] {1,0,1,1}) ));
		
		F_s fprob; // dummy reference variable
		
		boolean computed_something_new = true;
		while(! (f1000.computed && f1001.computed)  ) {
			// make sure infinite loop doesn't go forever, check if previous pass through while loop was productive
			if(!computed_something_new) {
				System.out.println("While loop iteration did not compute more states");
				System.exit(0);
			}
			
			computed_something_new = false;
			
			
			// compute whatever we are able to 
			for(int i = 0 ; i < f_probability_functions.size() ; i++) {
				fprob = f_probability_functions.get(i);
				fprob.computeIfAble();
				
				if(fprob.computed) { // if we successfully computed
					
					// take values, and transform into every other coordinate frame immediately downstream
					fprob.shiftCoordFeed(f_probability_functions, f_key);
					
					// now remove from the Fs and the keys, mark this iteration as productive
					computed_something_new = true; // this pass of while loop was productive
					f_key.remove(i);
					f_probability_functions.set(i,null);
					f_probability_functions.remove(i);
					i--; // correct i for the removed element. note f_p_f size should automatically decrease now.
				}
			}	
			
		}
	
		/////
	
		RealVector xys = t_grid.mapMultiply(n_samples);
		
		;
		f1001_vals = new Function_BiSpline(xys, xys, f1001.F_s_vals[ num_disc_times - 1 ], switch_linInt);
		f1000_vals = new Function_Spline(xys.toArray(), f1000.F_s_vals[num_disc_times - 1].getRowVector(num_disc_times - 1).toArray(), -1., switch_linInt);


		////////////////////////////////////////
		RealVector f1_ft = f1000.F1_vals;
		f1_ft = f1_ft.add(f1001.F1_vals);
		f1_ft = f1_ft.add(f1011.F1_vals);
		
		f1 = new Function_Spline(xys.toArray(), f1_ft.toArray(),-1., switch_linInt);
		//////////////////////////////////////////////////////
		
	}
	

	private class F_s{
		// the coordinations are tied to the global t_grid and the characteristics for state s.
		
		int state;	// integer identifying which state this is
		RealMatrix[] F_s_vals;
		
		ArrayList<Integer> necessaryStates = new ArrayList<Integer>(num_states);
		ArrayList<Integer> downStreamStates = new ArrayList<Integer>(num_states);
		
		ArrayList<Integer> obtainedStates; 
		ArrayList<RealMatrix[]> obtainedFvals;
		
		RealVector F1_vals;
		
		boolean computed;
		boolean absorbing;
		boolean singleXlin;
		boolean computable1d;	 // a boolean to mark that the a and b loci have collected tree length at the same rate up till this point (X000, and X111 states)	
		boolean symmetric; // if not symmetric, this is a x01y format state, see f_obtainOrder
		
		int vA;
		int vB;
		
		F_s(int s) {
			
			state = s;
			absorbing = (s == stateToInt(new int[] {1,0,0,0}) || s == stateToInt(new int[] {1,0,0,1}));
			singleXlin = (absorbing || s == stateToInt(new int[] {1,0,1,1}) ) ;
			
			int[] s_xxxx = intToState(state);
			symmetric = (s_xxxx[1] == s_xxxx[2]);
			computable1d = (s_xxxx[3]==0);
			computable1d = computable1d || (s_xxxx[1]==1 && s_xxxx[2]==1 );
	
			// populate necessary states with the indices of states needed to compute this F matrix, do same for downStream states
			for(int i = 0 ; i < num_states ; i++) {
				int[] sp = intToState(i);
				if(recoRates.getEntry(i, state) + coalWeights.getEntry(i, state) > 0.) {
					if(cVecA[i] == n_samples && cVecB[i] == n_samples) { } // don't inherit from singularly defined states
					else {	necessaryStates.add(i); }
				}
				if(recoRates.getEntry(state, i) + coalWeights.getEntry(state, i) > 0) {
					if(sp[1]==1 && sp[2]==0) {} // don't include X10Y states in downstream states
					else {	downStreamStates.add(i); }
				}
			}
			
			// initialize list of obtained F_s
			obtainedStates = new ArrayList<Integer>(necessaryStates.size());
			obtainedFvals = new ArrayList<RealMatrix[]>(necessaryStates.size());
			
			// starts in default false position for computed
			computed = false;
			
			// initialize vA and vB 
			vA = cVecA[state];
			vB = cVecB[state];
			
		}
		
		void computeIfAble() {
			
			// make sure not already computed
			if(computed) {return;} 
			
			// go through necessary states and make sure it has what it needs, otherwise terminate.
			for(int i = 0 ; i < necessaryStates.size() ; i++) {
				if( !obtainedStates.contains(necessaryStates.get(i)) ) { return; }
			}
			
			//System.out.println("computing: "+state+": " + Arrays.toString(intToState(state)));
			System.out.print(".");

			/////////////////////////
			/////  MEAT OF METHOD ///////////////////////////////////////////////////////////////////////////
			/////////////////////////
			
			F_s_vals = create0s();
		
			// GET VALUES OF BOUNDARY ODE ON T_GRID
			double[] boundary_corner = ode_grid_sol[state];
					
			// GET RATES FLOWING INTO THIS STATE ON TIME GRID FOR THE RELEVANT STATES
			RealMatrix inRates = getInfluxRates();
			
			
			// COMPUTE EXP[ -H2 ] VALUES ON GRID (START FROM T = 0)
			RealVector exp_h2_s = getEH2grid();   // note that this is really ex( - h2 ) even though we have named it as eh2
			

			// COMPUTE MARGINAL BOUNDARY ON THE X,Y = NT BOUNDARIES
			computeOuterBoundaries(boundary_corner, exp_h2_s, inRates);
			// note the inner boundaries, where x,y = va*t,vb*t are already all set to zero
							
			// COMPUTE VALUES ALONG EACH INTERIOR CHARACTERISTIC
			computeInteriorCharacteristics(exp_h2_s , inRates);
		
			// IF THIS IS F1 STATE, THEN COMPUTE CONTRIBUTION
			contributeF1(boundary_corner, exp_h2_s, inRates);
			////////////////////////////////////////////////////////////////////////////////////////////////////
			
			// mark this object as being computed (clears memory)
			markComputed();
			
			
		}
		

		///////////////////////////////////////////////////
		// AUXILLIARY METHODS ////////////////////
		////////////////////////////////

		private RealMatrix getInfluxRates() {
			// this computes a matrix of rates. each row is a different valid incoming state, each column is a different point on time grid. 
			// the values are the value by which residence prob should be multiplied to find the actual influx
			
			if( necessaryStates.size() == 0 ) {return null; }
			
			RealMatrix influxRates = new Array2DRowRealMatrix(obtainedStates.size(), num_disc_times);
			RealVector influxRR = new ArrayRealVector(obtainedStates.size());
			RealVector influxCW = new ArrayRealVector(obtainedStates.size());

			for(int u = 0 ; u < obtainedStates.size() ; u++) { 
				int inState = obtainedStates.get(u);
				influxRR.setEntry(u, recoRates.getEntry(inState, state));
				influxCW.setEntry(u, coalWeights.getEntry(inState, state));
			}
						
			for(int t = 0 ; t < num_disc_times ; t++) {
				double crate = coalescence_rate.evaluate(t_grid.getEntry(t));
				influxRates.setColumnVector(t, influxCW.mapMultiply(crate).add(influxRR));
			}
			
			return influxRates;
		}
		
		private	RealVector getEH2grid() {
			RealVector exp_h2_s = new ArrayRealVector( new double[num_disc_times] );
			double exp_h2 = 1;
			exp_h2_s.setEntry(0, exp_h2);
						
			for(int t = 1 ; t < num_disc_times ; t++) {
				double h2s = - coalWeights.getEntry(state, state) * (lambda_grid.getEntry(t) - lambda_grid.getEntry(t-1)  )  ; // coalescence part
				h2s = h2s - recoRates.getEntry(state, state) * (  t_grid.getEntry(t) - t_grid.getEntry(t-1)  );
							
				exp_h2 = exp_h2 * Math.exp(- h2s);
				exp_h2_s.setEntry(t, exp_h2 );
			}
			return exp_h2_s;
		}
		
		private RealVector integrateOnChar( RealVector g2h2 , double iVal , int[] char_ij, String var) {
			
			double node0x = -1.;
			int i=-1;
			int j=-1;
			int ti0 = char_ij[0];
			
			// CASEWORK TO CREATE CORRECT NODE ENTRY, BASED ON TYPE OF CHARACTERISTIC WE ARE INTEGRATING
			if(char_ij.length == 3) { // CHARACTERISTIC ON INTERIOR
				i = char_ij[1];
				j = char_ij[2];
				
				double[] coords = this.getCoords(ti0, i, j);
				double x0 = coords[1];
				double y0 = coords[2];
				double t0 = coords[0];
			
				
			 // CHECK AT WHICH TIME YOU HIT THE UPSTREAM GRID BOUNDARY (0 BOUND) AS YOU FOLLOW CHARACTERISTIC.
				
				double rho = 0; // how far you follow characteristic till the values are all 0 (you hit grid boundary of upstream state)
				for(int s = 0 ; s < this.obtainedStates.size() ; s++) {
					int vp_a = cVecA[obtainedStates.get(s)];
					int vp_b = cVecB[obtainedStates.get(s)];
						
					double rho_x = d_times[num_t_bins-1];
					double rho_y = d_times[num_t_bins-1];
						
					// if vp = v then you are parallel to boundary, so pick max t value as the node point)
					if((vp_a - vA) != 0) {rho_x = (x0 - t0 * vp_a) / (vp_a - vA);  } // compute rho along each dimension
					if((vp_b - vB) != 0) {rho_y = (y0 - t0 * vp_b) / (vp_b - vB);  }

					rho = Math.max(rho, Math.min(rho_x, rho_y)); // pick which dimension constrains more (min), then make sure it is non-negative (max)
					// min function can yield negative values especially for absorbing state, so the max will just ensure rho is set to 0 so its well defined
				}
					
				node0x = rho + t0;
				
			}
			
			else if(char_ij.length == 1) { // CHARACTERISTIC IN ONE DIMENSION, ALONG SURFACE OF THE RECTANGULAR CONE
				double[] coords = this.getCoords(ti0, ti0, ti0);
				double x0 = coords[1];
				double y0 = coords[2];
				double t0 = coords[0];
				
				
				double rho = 0;
				
				for(int s = 0 ; s < this.obtainedStates.size(); s++) {
					int vp; 
					int vV;
					double rho_temp = d_times[num_t_bins-1];
					double z0;
							
					if(var == "X") {  vp = cVecA[obtainedStates.get(s)]; vV = vA; z0 = x0;}
					else if(var == "Y") {  vp = cVecB[obtainedStates.get(s)]; vV = vB; z0 = y0;}
					else {vp = 0; vV = 0; System.out.println("Invalid Variable Type"); System.exit(0); z0 = 0;}

					if((vp - vV) != 0) {rho_temp = (z0 - t0 * vp) / (vp - vV);  }

					rho = Math.max(rho, rho_temp);
					
				}
			
				
				node0x = rho + t0;

			}
			
			else {ti0 = 0;  System.out.println("invalid char defined");  System.exit(0);;}
		
			
		
			// t0 is index of start time, g2h2 is the vector to integrate. note the g2h2 length should be set appropriately. iVal is starting value.
			Function_Spline g2_eh2 = new Function_Spline( t_grid.getSubVector(ti0, num_disc_times - ti0).toArray() , g2h2.toArray() , node0x, switch_linInt);						
			
			///////////////////////////////////////////////////////////////////////////
			return g2_eh2.integrate(iVal, cr_discontinuities);
			
		}
		
		private void computeOuterBoundaries(double[] boundary_corner, RealVector exp_h2_s , RealMatrix inRates) {
			/// THIS METHOD FILLS IN THE OUTER BOUNDARY SURFACES ALONG THE X = NT AND Y = NT SURFACES
			
			// write in the value of the last characteristic we will go along
			F_s_vals[num_disc_times - 1].setEntry(num_disc_times - 1, num_disc_times - 1, boundary_corner[num_disc_times - 1]);
													

			// proceed characteristic by characteristic
			for(int c = 1 ; c < num_disc_times - 1 ; c++) {
							
				// obtain the exp( - h2 ) vector by taking a subvector from the full grid version and scaling appropriately
				RealVector eh2 = exp_h2_s.getSubVector(c, num_disc_times - c);
				eh2.mapDivideToSelf(eh2.getEntry(0));
				
				// obtain g*e integrand
				RealVector g2X = new ArrayRealVector(num_disc_times - c);
				RealVector g2Y = new ArrayRealVector(num_disc_times - c);
							
				for(int u = 0 ; u < obtainedStates.size() ; u++) { // need to obtain contribution form each input state separately
					RealMatrix[] fspm = obtainedFvals.get(u);
						
				
					
					RealVector fsprimeX = new ArrayRealVector(num_disc_times - c);
					RealVector fsprimeY = new ArrayRealVector(num_disc_times - c);
					for(int t = c ; t < num_disc_times ; t++) {
						fsprimeX.setEntry(t - c, fspm[t].getEntry(c, t) ); 
						fsprimeY.setEntry(t - c, fspm[t].getEntry(t, c) ); 
					}
						
					
					// Add two lines to ensure that the boundary is not included in the integration if the function is not well defined at the boundary. 
					// ie if we get influx at the corner, but not along the rest of the characteristic, then we do not want to include the corner contribution under an integral that spans more than the corner.
					if(cVecA[obtainedStates.get(u)] == n_samples) {fsprimeX.setEntry(0, 0);}
					if(cVecB[obtainedStates.get(u)] == n_samples) {fsprimeY.setEntry(0, 0);}
						
					
					RealVector q = inRates.getRowVector(u).getSubVector(c,  num_disc_times - c);
					g2X = g2X.add(fsprimeX.ebeMultiply(q));
					g2Y = g2Y.add(fsprimeY.ebeMultiply(q));
					
				}
				
				// get vector containing the integrand of the integral. 
				RealVector g2_eh2X = g2X.ebeDivide(eh2) ;
				RealVector g2_eh2Y = g2Y.ebeDivide(eh2) ;

				
				
				// integrate along characteristic, and multiply by exp( - h2 )
				RealVector char_f_valsX = integrateOnChar(g2_eh2X, boundary_corner[c] , new int[] {c}, "X").ebeMultiply(eh2);
				RealVector char_f_valsY = integrateOnChar(g2_eh2Y, boundary_corner[c] , new int[] {c}, "Y").ebeMultiply(eh2);

				
				// write values into the F_s_values structure
				for(int t = c ; t < num_disc_times ; t++) {
					F_s_vals[t].setEntry(c, t, char_f_valsX.getEntry(t-c));
					F_s_vals[t].setEntry(t, c, char_f_valsY.getEntry(t-c));
				}
							
							
			}
		}

		private void computeInteriorCharacteristics( RealVector exp_h2_s, RealMatrix inRates) {
			
			// use shortcut for unrecombined states.
			
			if(computable1d) { 
				//if(!absorbing) {return;}
				project_1d_Interior(); 
				return;
			}
			
			// otherwise compute normally
			
			for(int i = 1 ; i < num_disc_times - 1 ; i++) {       // WE ONLY WANT TO COMPUTE THE INTERIOR CHARACTERISTICS, IE NOT ON THE BOUNDARY SURFACES. 
				int jlb = 1;
				if(symmetric) { jlb = i; } // in symmetric case, we can do half of the characteristics
				for(int j = jlb ; j < num_disc_times - 1 ; j++) {
					
					
					int t0 = Math.max(i, j); // this gives time index of the first point on our grid for this characteristic. this point should be filled in already
					double iVal = F_s_vals[t0].getEntry(i, j);
			
					// obtain the exp( - h2 ) vector by taking a subvector from the full grid version and scaling appropriately
					RealVector eh2 = exp_h2_s.getSubVector(t0, num_disc_times - t0);
					eh2.mapDivideToSelf(eh2.getEntry(0));
			
					// obtain g*e integrand
					RealVector g2 = new ArrayRealVector(num_disc_times - t0);
					for(int u = 0 ; u < obtainedStates.size() ; u++) { // need to obtain contribution form each input state separately
						RealMatrix[] fspm = obtainedFvals.get(u);
							
						RealVector fsprime = new ArrayRealVector(num_disc_times - t0);
						
						// Ensure that the boundary is not included in the integration if the function is exclusively defined on a boundary rather than a grid 
						// ie: if we get influx at the surface, but not along the rest of the characteristic, then we do not want to include the surface contribution under an integral that spans the interior
						boolean bc_only = (cVecA[obtainedStates.get(u)] == n_samples && i >= j)  ||  (cVecB[obtainedStates.get(u)] == n_samples && i <= j) ; // true if its boundary confined only on the start point of fsprime
						
						if(bc_only) {fsprime.setEntry(0, 0);} 
						else {
							for(int t = t0 ; t < num_disc_times ; t++) {
								fsprime.setEntry(t - t0, fspm[t].getEntry(i, j) ); // populate fsprime influx residences if well defined.
							}
						}
						
						
						RealVector q = inRates.getRowVector(u).getSubVector(t0,  num_disc_times - t0);
						g2 = g2.add(fsprime.ebeMultiply(q));
						
					}
					
			
					// get vector containing the integrand of the integral. 
					RealVector g2_eh2 = g2.ebeDivide(eh2) ;
			

					
					// integrate along characteristic, and multiply by exp( - h2 )
					RealVector char_f_vals;
					char_f_vals = integrateOnChar(g2_eh2, iVal , new int[] {t0,i,j}, "").ebeMultiply(eh2) ;  

					
					// write values into the F_s_values structure
					for(int t = t0 ; t < num_disc_times ; t++) {
						F_s_vals[t].setEntry(i, j, char_f_vals.getEntry(t-t0));
						if(symmetric) {F_s_vals[t].setEntry(j, i, char_f_vals.getEntry(t-t0));}
					}
			
				}
			}
				
		}
	
		private void project_1d_Interior() {
			
			// fill in the interior of an unrecombined state simply by projecting appropriately from boundaries. 
			for(int t = 1 ; t < num_disc_times; t++) {
				RealMatrix f_slice = F_s_vals[t];
				for(int i = 1 ; i < t + 1 ; i++) {
					double limiting_val = f_slice.getEntry(i, t);
					for(int j = i ; j < t + 1 ; j++) {
						f_slice.setEntry(i, j, limiting_val);
						f_slice.setEntry(j, i, limiting_val);
					}
				}
			}
		}
		
		void contributeF1(double[] boundary_corner, RealVector exp_h2_s , RealMatrix inRates) {
			
			if(singleXlin) {
				
				F1_vals = new ArrayRealVector(num_disc_times);
				
				for(int c = 1 ; c < num_disc_times - 1 ; c++) {
					
					// obtain the values of F(x/2 , x , xn/2) which is the contribution to F1(x/2,x)
					
					double[] F_vert = new double[num_disc_times - c ];
					
					for(int t = c ; t < num_disc_times ; t++) {
						F_vert[t - c] =  F_s_vals[t].getEntry(c, t);
					}
					
					// what time is the x/2 point for this characteristic that we want value at
					double target_t = t_grid.getEntry(c) * n_samples / 2;
					
					// use spline if the target time is within our time grid.
					if(target_t <= t_grid.getEntry(num_disc_times - 1)) {
						
						// now define the spline_func we'll use to interpolate the value at the in between cut-off time point 
						Function_Spline spl = new Function_Spline(t_grid.getSubVector(c, num_disc_times - c).toArray(), F_vert, -1., switch_linInt);
						
						// infer value and write to holder vector
						F1_vals.setEntry(c, spl.evaluate(target_t));
					}
					
					// if target time is beyond grid, 
					else {
						F1_vals.setEntry(c, F_vert[F_vert.length - 1]);
					}
					
				}
				
			}
			
			
			
		}
	
		
		/////////////////////////////////////////////////
		
		private RealMatrix[] create0s() {
			RealMatrix[] temp = new RealMatrix[num_disc_times];
			for(int i = 0 ; i < num_disc_times ; i++) {
				temp[i] = new Array2DRowRealMatrix(i+1, i+1);
			}
			return temp;
		}
		
		private void markComputed() {
			
			// delete references to needless objects, only hang on to final set of values.
			necessaryStates = null;
			obtainedStates = null;
			obtainedFvals = null;
			
			computed = true;
			
		}
		
		private RealMatrix[] xySwitch(RealMatrix[] m) {
			RealMatrix[] out = new RealMatrix[m.length] ;
			for(int i = 0 ; i < m.length ; i++) {
				out[i] = m[i].transpose();
			}
			return out;
		}
		
		
		
		////////////////////////////////////////////////////
		
		void shiftCoordFeed( ArrayList< F_s > f_prob_functions, ArrayList<Integer> f_key) { 
						
			// initialize the holding matrices for the values, all set to 0 //
			ArrayList<RealMatrix[]> shiftedMatrices = new ArrayList<RealMatrix[]>(downStreamStates.size());
			for(int d = 0 ; d < downStreamStates.size() ; d++) {
				shiftedMatrices.add( create0s() );
			}
			
			// HANDLE CASE WHERE GRID IS FULLY UNDEFINED EXCEPT FOR AT A POINT AT EACH TIME
			if( vA == n_samples && vB == n_samples ) {
				// in this case there is 1 point per time for grid, at outer corner, but that gets overwritten by ODE, so just pass zeros
				return ;
			}
			
			// setup the scaled grid from which we take subvectors and add a term later to obtain xp and yp grids
			RealVector x0h = t_grid.mapMultiply(n_samples - vA);
			RealVector y0h = t_grid.mapMultiply(n_samples - vB); 
			ArrayList<RealVector> xph = new ArrayList<RealVector>(downStreamStates.size());
			ArrayList<RealVector> yph = new ArrayList<RealVector>(downStreamStates.size());
						
			// deal with the t = 0 matrix for each downstream state, and add the appropriate scaled grid (xph, yph) into the holder
			for(int d = 0 ; d < downStreamStates.size() ; d++) {
				shiftedMatrices.get(d)[0] = new Array2DRowRealMatrix(1,1); //t = 0 matrix
				xph.add(t_grid.mapMultiply(n_samples - cVecA[downStreamStates.get(d)]));
				yph.add(t_grid.mapMultiply(n_samples - cVecB[downStreamStates.get(d)]));
			}
		
			
			// iterate through each time step, performing the splineInference and Interpolation at each time
			for(int t = 1 ; t < num_disc_times ; t++ ) {
				double ct = t_grid.getEntry(t); // obtain the current time
				
				RealVector xhs = x0h.getSubVector(0, t+1);
				RealVector yhs = y0h.getSubVector(0, t+1);
						
				RealVector x0_grid = xhs.mapAdd(ct * vA);  	// the discrete x values along the characteristics for this state's grid
				RealVector y0_grid = yhs.mapAdd(ct * vB);  	// "" for y values instead
				
				Function_Spline one_interpSlice = null; // 1d version if needed 
				Function_BiSpline interpSlice = null;
				
				
				
				
				// learn the function with spline interpolation
				if(vA == n_samples) { // partially defined grid. one dimension's coordinates are degenerate
					one_interpSlice = new Function_Spline(y0_grid.toArray(), F_s_vals[t].getRow(t), -1., switch_linInt);
				}
				
				else if(vB == n_samples || computable1d) { // partially defined in other direction. 
					one_interpSlice = new Function_Spline(x0_grid.toArray(), F_s_vals[t].getColumn(t), -1., switch_linInt);
				}
				
				
				else {// well defined grid
					interpSlice = new Function_BiSpline( x0_grid, y0_grid, F_s_vals[t], switch_linInt);

				}
				// Now the function learning is complete. now map this info onto the downstream grids
				
				
				
				
				for(int d = 0 ; d < downStreamStates.size() ; d++) {
					int target_state = downStreamStates.get(d); // actual state_int of the target
					int vAt = cVecA[target_state];  // target state vA char value
					int vBt = cVecB[target_state];	// "" for B locus
					

					// if target grid is identical, can ignore the interpolation and just copy over the matrix
					if( vAt == vA && vBt ==vB) {
						shiftedMatrices.get(d)[t] = F_s_vals[t].copy();
					}
					
					// interpolation procedure from the inferred function (interpSlice)
					else {
						RealMatrix target_matrix = shiftedMatrices.get(d)[t] ; // initialize the holder matrix
						
						// obtain the xp and yp grid of x,y vals for the target state
						RealVector xp_grid = xph.get(d).getSubVector(0, t+1).mapAdd(ct * vAt);
						RealVector yp_grid = yph.get(d).getSubVector(0, t+1).mapAdd(ct * vBt);
						
						
						// interpolate at the correct coordinate value, and write into target_matrix
						for(int i = 0 ; i < t+1 ; i++) {
							
							if(vA == n_samples ) {
								target_matrix.setEntry(t, i, one_interpSlice.evaluate_0(yp_grid.getEntry(i)) );
							}
							else if(vB == n_samples) {
								target_matrix.setEntry(i, t, one_interpSlice.evaluate_0(xp_grid.getEntry(i)));
							}
							
							else if(computable1d){
								for(int j = 0 ; j < t+1; j++) {
									double val = one_interpSlice.evaluate_0(Math.min(xp_grid.getEntry(i), yp_grid.getEntry(j)));
									target_matrix.setEntry(i, j, val);
								}
							}
							
							else {// this is a normal grid 
								for(int j = 0 ; j < t+1; j++) {
									target_matrix.setEntry(i, j, interpSlice.value0(  xp_grid.getEntry(i), yp_grid.getEntry(j) )   );
								}
							}
								
						}
					
						
					}
					
					
				}
			
				
			}
			
			
			// HAVE SO FAR COMPUTED MATRICES TO PASS DOWNSTREAM. 
			// NEED TO PASS THEM, AND ALSO PASS DOWNSTREAM FOR THE SYMMETRIC CASE
			
			
			// PASS VALUES
			for(int j = 0 ; j < downStreamStates.size() ; j++) {
				int target_si = downStreamStates.get(j);
				F_s copy_target = f_prob_functions.get(f_key.indexOf(target_si)) ;
				
				copy_target.obtainState(state, shiftedMatrices.get(j));
			}
			
			if(!symmetric) { // in this case this state is X01Y, and we need to replicate passing on the X10Y state
				int[] si = intToState(state);
				si[1] = 1 - si[1];
				si[2] = 1 - si[2];
				int mirror_state = stateToInt(si);
				
				for(int j = 0 ; j < downStreamStates.size() ; j++ ) {
					int target_si = downStreamStates.get(j);
					F_s copy_target = f_prob_functions.get(f_key.indexOf(target_si));
					
					if(copy_target.necessaryStates.contains(mirror_state)) { // DOUBLE CHECK MIRRORED STATE IS APPROPRIATE FOR THIS DOWNSTREAM STATE
						copy_target.obtainState(mirror_state, xySwitch(shiftedMatrices.get(j)));
					}
				}
				
				
				
			}
			
			
			
			
			
			
		}
		
		void obtainState(int s, RealMatrix[] shiftedVals) {
			if(!necessaryStates.contains(s)) {
				System.out.println("state " + state + " was incorrectly given shifted-values from state " + s );
				return;
			}
			
				
			obtainedStates.add(s);
			obtainedFvals.add(shiftedVals);
		}
		
		double[] getCoords(int t, int i, int j){
			
			double time = t_grid.getEntry(t);
			
			double x = t_grid.getEntry(i) * (n_samples - vA) + time * vA;
		
			double y = t_grid.getEntry(j) * (n_samples - vB) + time * vB;
			
			double[] coords = new double[] {time, x, y};
			
			return coords;
		}
		
	}
	
	
	
	
	
	
	
//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


	public static void main(String[] args) {
		
		
		Helper helper = new Helper();
		
		double rec_rate = 0.0000000125;
		double mut_rate = 0.0000000125;
		int n_samples = 10;
		
		// initial value for psh
		double psh0 = 10000; // to define states
		
		///////  .02,.05,.15,.25... for n = 33 , N = 15000, for tree height
		double[] pp = new double[] {.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05};
		RealVector dlengths = helper.getTreeLengthPartitions(pp, n_samples);
		dlengths.mapMultiplyToSelf(psh0 * 2);

		
		double sfactor = 1000;
		
		double [] xs = new double[] {500,1000,2000,3000,4500,6750,10000,15000,22500};
		double[] ys = new double[] {8711.760278276202, 7261.460809309039, 10492.929636574328, 8343.054523899298, 8476.626532956956, 6918.35376842471, 8600.189111402131, 16931.04993072053, 8574.763346698017, 10290.756873186861};
		
		
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
		
		
		TL_1P_Transitions tp_comp = new TL_1P_Transitions(cr,rho, n_samples, scaled_lengths, new int[] {180,180}); 

			
		System.out.println("********** \n \n STATE PARTITIONS:\n");
		System.out.println(dlengths);
		System.out.println("\n \n **********");
		System.out.println("************************");
		System.out.println("NUMERIC JOINT PDF:");
		System.out.println(tp_comp.getPDF());
		System.out.println("\n \n NUMERIC TRANSITION PROBS:");
		System.out.println(helper.pdfToRowProbMatrix(tp_comp.getPDF(), 0));
		System.out.println("\n \n \n ************************");
		//System.out.println("SIMULATION COMPARISON (10^8 RUNS) ");
		//tp_comp.compareSim(100000000);
		
		System.out.println("************************");

		System.out.println("************************");


	}
		
	
	 
}
