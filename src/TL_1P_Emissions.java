import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;



public class TL_1P_Emissions  extends EmissionsComputer{

	
	
	// concerning time grid/ internal discretization
	int num_disc_times;
	RealVector t_grid;
	
	// concerning characteristic val for each state
	int[] vChar;
	
	// use this to compute combinatoric factors
	CFactor cfacs;

	
	// things to be recomputed every time we update parameters. the final objects we need are the lists of absorbing states' Fvals
	
	RealVector lambda_grid;
	double[][] ode_grid_sol;

	double[] cr_discontinuities;
	
	// quantities for absorbing states, want to compute as meat of this method
	ArrayList< Function_Spline > absorbingFvals;
	ArrayList< Integer > absorbingKey; 

	
	
	// CONSTRUCTOR AND HELPING FUNCTIONS

	
	TL_1P_Emissions(FunctionUnivariate cr, double mr, int ns, double[] dts, int gsize, boolean pseudo_hap){
		super(cr, mr, ns, dts, pseudo_hap);
		
		
		buildVchars();
		createTimeGrid(gsize);
		createLambdaGrid();
		
		cfacs = new CFactor();
		
	}
	
	private void createTimeGrid(int gs) {
		
		num_disc_times = gs + 1;
		double[] ts = new double[num_disc_times];
		
		// uniformly spaced grid
		double maxt = d_times[num_t_bins - 1]/2. ;
		double dt = maxt / gs ;
		
		// write values into ts and lambdas
		for(int ti = 0 ; ti < num_disc_times; ti++) {
			ts[ti] = dt * ti ;
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
	
	
	private void buildVchars() {
		
		vChar = new int[num_states];
		for(int s = 0 ; s < num_states ; s++) {
			
			int[] state = intToState(s);
			int lins = state[0];	
			if( lins > 1) {vChar[s] = lins ;}
			
		}
		
	}


///////////////////////////////////////////////////////////////////////////
// INHERITED ABSTRACT METHODS. MUST IMPLEMENT
////////////////////////////////////////////////////////////////////////////
	
	
	
	@Override
	void updateCDF() {
		double[][] cdf = new double[num_t_bins + 1][n_samples];
		
		// run method to populate the f functions if final F vals are not up to date
		if(!upToDate) { computeFinalFMats(); }
		
		// use the absorbing state F vals to compute
		
		for(int i = 0 ; i < num_t_bins ; i++) {
			double x = d_times[i];
			cdf[i][0] = absorbingFvals.get(absorbingKey.indexOf(0)).evaluate(x);
			
			for(int ss = 1 ; ss < n_samples; ss++) { // for each emission entry in the row
				double contribution = 0 ; 
				
				for(int ks = 2 ; ks < n_samples + 1 ; ks++) {// sum contributions from each absorbing state
					double fp = absorbingFvals.get(absorbingKey.indexOf(ks)).evaluate(x);
					contribution = contribution + fp * cfacs.get(ks, ss);
				}
				
				cdf[i][ss] = contribution;
				
			}
			
		}
		
	
		
		// now populate last row of CDF (just distro of ss)
		// note last column (j = num_disc_times) of ode_grid_sols is inifnity vals
		
		cdf[num_t_bins][0] = ode_grid_sol[stateToInt(new int[] {1,0})][num_disc_times]; 
		for(int ss = 1 ; ss < n_samples; ss++) { // for each emission entry in the row
			double contribution = 0 ; 
			
			for(int ks = 2 ; ks < n_samples + 1 ; ks++) {// sum contributions from each absorbing state
				double fp = ode_grid_sol[stateToInt(new int[] {1,ks})][num_disc_times];
				contribution = contribution + fp * cfacs.get(ks, ss);
			}
			
			cdf[num_t_bins][ss] = contribution;
			
		}

		
	
		// final length bin should be computed from ode_inf_sol
		currentCDF = new Array2DRowRealMatrix(cdf);
		upToDate = true;
		
	}


	@Override
	int[] ancestralSim() {

		// setup the tree_lengths to return, at each of 2 loci
		double tl = 0.; 
		int seg_sites = 0;
				
				
		// initialize the starting representation of state, where all lineages subtend one modern sample
		int[] state = new int[n_samples];  // how many samples subtended by extant lineages. the extant lineages are first "lineages" elements of state
		Arrays.fill(state, 1);         
		int lineages = n_samples;          // number of ancestral lineages extant
		double sim_time = 0.;			// simulation time
		double lc_time = 0.;
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
				// update tree length
				tl = tl + (sim_time - lc_time) * lineages;
				
				// update state 
				int c_l = ind_gen.nextInt(lineages - 1);
				state[c_l] = state[c_l] + state[lineages - 1];    // rightmost extant lineage feeds into another random lineage					
				
				// update times/lineages
				lineages = lineages - 1;  // coalescence will decrease number of lineages by 1.
				lc_time = sim_time;
			}
			
		}
		
		int l_index = helper.placeInPartition(tl, d_times);
		
		return new int[] {l_index, seg_sites}; 
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

	
	
////////////////////////////////////////////////////////////////////
	// helping methods ///////////
	class CFactor{

				Hashtable<List<Integer>, Double> cTensor;
				Helper.ChooseComputer choose = helper.new ChooseComputer(n_samples);
				
				CFactor(){
					cTensor = new Hashtable<List<Integer>, Double>();
					buildHash();
				}
				
				void buildHash() {
					for(int ks = 2 ; ks < n_samples + 1 ; ks++) {
						for(int ss = 1 ; ss < n_samples  ; ss++) {
							double fac = computeCoef(ks,ss);
							if(fac!=0) {cTensor.put( Arrays.asList(ks, ss), fac);}
						}
					}
				}
				
				double computeCoef( int ks, int s) {
					return  choose.get( n_samples - s - 1, ks - 2) / choose.get1(n_samples - 1, ks - 1);
				}

				double get(int ks, int ss) {
					Double c = cTensor.get(Arrays.asList(ks, ss));
					if(c == null) {c = 0.;}
					return c;
				}
				
				
				
	}
	
	private void computeFinalFMats() {
		
		// get these now rather than having to make a copy every time we integrate on a char line
		cr_discontinuities = coalescence_rate.getDiscontinuities();
				
		// create lambda grid based on the coalescence rate function
		createLambdaGrid();
		
		
		// integrate the ODE to obtain the residence probabilities at the times of interest, for all states
		ode_grid_sol = helper.odeSolsAtTimes(this, getInitialResProbs(), t_grid.toArray(), t_step_precision, true, this.getAbsorbingStateIndices());
		
				
		// initialize the holders to hold absorbing states' final vals
		absorbingFvals = new ArrayList< Function_Spline > (n_samples);
		absorbingKey = new ArrayList< Integer > (n_samples);
		
		
		// initialize the f function array containing an f object for each state
		ArrayList< F_s > f_probability_functions = new ArrayList< F_s > (num_states);
		ArrayList<Integer> f_key = new ArrayList<Integer> (num_states); // values are the k* value of the absorbing state (k = 1 obvi)
		
		// now fill all the state objects in order which we should cycle through list in some ordering
		f_obtainOrder(f_key, f_probability_functions);
		
		
		
	
		F_s fprob; // dummy reference variable
		
		boolean computed_something_new = true;
		while( f_probability_functions.size() > 0 ) {
			
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
					
					// if it was absorbing store the desired values
					if(fprob.absorbing) {
						
						RealVector xg = t_grid.mapMultiply(n_samples);  // x grid
						RealVector f_x2t = fprob.F_s_vals[num_disc_times - 1]; // F values on x = 2t line
						
						absorbingFvals.add(  new Function_Spline(xg.toArray(), f_x2t.toArray(), -1, false)  ); // infer spline and store
						absorbingKey.add(intToState(fprob.state)[1]);
					}
					
					
					// take values, and transform into every other coordinate frame immediately downstream
					fprob.shiftCoordFeed(f_probability_functions, f_key);
					
					// now remove from the Fs and the keys, mark this iteration as productive
					computed_something_new = true; // this pass of while loop was productive
					f_key.remove(i);
					f_probability_functions.set(i, null);
					f_probability_functions.remove(i);
					i--; // correct i for the removed element. note f_p_f size should automatically decrease now.
				}
			}	
			
		}
	
		/////
	
		// we have computed the absorbing states' F_s(x/2, x)
		
		//////////////////////////////////////////////////////
		System.out.print("!");

	}


	private void f_obtainOrder(ArrayList<Integer> f_key, ArrayList< F_s > f_probability_functions) {
		
		// this is the portion to try optimal ordering to reduce heap usage
		int s;
		
		for(int i = n_samples; i > 1 ; i-- ) {
			s = stateToInt( new int[] {i,0} );
			f_probability_functions.add(new F_s(s));
			f_key.add(s);
			
			for(int j = i ; j > 0 ; j--) {
				s = stateToInt( new int[] {j,i} );
				f_probability_functions.add(new F_s(s));
				f_key.add(s);
			}
			
		}
		
		s = stateToInt( new int[] {1,0} );
		f_probability_functions.add(new F_s(s));
		f_key.add(s);
		
	}

	private void createLambdaGrid() {
		double[] lambdas = new double[num_disc_times];

		for(int ti = 0 ; ti < num_disc_times; ti++) {
			lambdas[ti] = coalescence_rate.integrate0(t_grid.getEntry(ti));
		}
		
		lambda_grid = new ArrayRealVector(lambdas);
	}


	private class F_s{
		
		int state;	// integer identifying which state this is
		RealVector[] F_s_vals; // holder for Fvals
		
		int necessaryState = -1;  // there is only ever one state that feeds into a given state
		ArrayList<Integer> downStreamStates = new ArrayList<Integer>(2);
		
		boolean obtainedState; 
		RealVector[] obtainedFvals;
				
		boolean computed;
		boolean absorbing;
		
		int v_char;
		
		F_s(int s) {
			
			// state and type
			state = s;
			absorbing = (intToState(s)[0] == 1);
					
			// populate necessary states with the indices of states needed to compute this F matrix, do same for downStream states
			boolean ncb = false;
			for(int i = 0 ; i < num_states ; i++) {
				if(mutRates.getEntry(i, state) + coalWeights.getEntry(i, state) > 0.) {
					if(ncb) {System.out.println("issue, more than 1 influx state"); System.exit(0);}
					if(vChar[i] != n_samples) {  // don't inherit from singularly defined states (only defined on single line)
						necessaryState = i; ncb = true;
					}
				}
				if(mutRates.getEntry(state, i) + coalWeights.getEntry(state, i) > 0) {
					downStreamStates.add(i);
				}
			}
			
			// initialize list of obtained F_s, if its a state defined singularly, then consider it as ready to be computed
			if( necessaryState == -1 ) {obtainedState = true;}
			else {	obtainedState = false; }
			
			// starts in default false position for computed
			computed = false;
			
			// initialize v
			v_char = vChar[state];
			
		}
		
		void computeIfAble() {
			
			// make sure not already computed
			if(computed) {return;} 
			
			// go through necessary states and make sure it has what it needs, otherwise terminate.
			if(!obtainedState) {return;}
			
			//System.out.println("computing: "+state+": " + Arrays.toString(intToState(state)));
			
			// initialize structure to hold values
			F_s_vals = create0s();
					
			//compute values along each characteristic separately
			computeInteriorCharacteristics();
		
			// clear memory and mark this object as being computed
			markComputed();
			
			
		}
		
	
		
		///////////////////////////////////////////////////
		// AUXILLIARY METHODS ////////////////////
		////////////////////////////////

		private RealVector getInfluxRates() {
			
			if( necessaryState == -1 ) {return null; }
			
			RealVector inRate = new ArrayRealVector( num_disc_times);
	
			double influxMR = mutRates.getEntry(necessaryState, state);
			double influxCW = coalWeights.getEntry(necessaryState, state);
						
			for(int t = 0 ; t < num_disc_times ; t++) {
				double crate = coalescence_rate.evaluate(t_grid.getEntry(t));
				inRate.setEntry(t, influxCW * crate + influxMR);
			}
			
			return inRate;
		}
		
		private	RealVector getEH2grid() {
			RealVector exp_h2_s = new ArrayRealVector( new double[num_disc_times] );
			double exp_h2 = 1;
			exp_h2_s.setEntry(0, exp_h2);
						
			for(int t = 1 ; t < num_disc_times ; t++) {
				double h2s = - coalWeights.getEntry(state, state) * (lambda_grid.getEntry(t) - lambda_grid.getEntry(t-1)  )  ; // coalescence part
				h2s = h2s - mutRates.getEntry(state, state) * (  t_grid.getEntry(t) - t_grid.getEntry(t-1)  );
							
				exp_h2 = exp_h2 * Math.exp(- h2s);
				exp_h2_s.setEntry(t, exp_h2 );
			}
			return exp_h2_s;
		}
		
		private RealVector integrateOnChar( RealVector g2h2 , double iVal , int ti0) {
			
			double node0x = -1.;
			
			if(necessaryState == -1) {}
			
			else {
				double[] coords = this.getCoords(ti0, ti0);
				double x0 = coords[1];
				double t0 = coords[0];
					
				int vp = vChar[necessaryState];   
				if((vp - v_char) != 0) {node0x = (x0 - t0 * vp) / (vp - v_char) + t0;  } //if not node0x is -1, dummy val
				
			}
			
			// t0 is index of start time, g2h2 is the vector to integrate. note the g2h2 length should be set appropriately. iVal is starting value.
			Function_Spline g2_eh2 = new Function_Spline( t_grid.getSubVector(ti0, num_disc_times - ti0).toArray() , g2h2.toArray() , node0x, false);						
			return g2_eh2.integrate(iVal, cr_discontinuities);
			
			
			
		}
		
		private void computeInteriorCharacteristics( ) {
			
			double[] boundary_vals = ode_grid_sol[state];      // vals along x = nt boundary from ode
			RealVector exp_h2_grid = getEH2grid(); 				// full grid of exp_h2 vals
			RealVector inRates = getInfluxRates();				// full grid of influx rates from upstream state
			
			// fill final corner entry of the last time slice, since this is just a value picked from the boundary vals, no integration
			F_s_vals[num_disc_times-1].setEntry(num_disc_times-1, boundary_vals[num_disc_times-1]);
			
			// find values on each characteristic one at a time
			for(int t0 = 1 ; t0 < num_disc_times - 1 ; t0++) {       // WE ONLY WANT TO COMPUTE THE INTERIOR CHARACTERISTICS, IE NOT ON THE BOUNDARY SURFACES. 
					
					// t0 is index of first time point along the char (x = nt boundary)
					double iVal = boundary_vals[t0];
			
					///////////////////////////////////////////////////////////////////
					// Obtain g_s vector on char, first obtain f_s' on characteristic
					
					RealVector g2 = new ArrayRealVector(num_disc_times - t0);
					
					if(necessaryState != -1) {// if there was no state necessary, then g2 is already initialized to 0s
						
						RealVector fsp_char = new ArrayRealVector(num_disc_times - t0); // F_s' along char
						for(int t = t0 ; t < num_disc_times ; t++) {
							fsp_char.setEntry(t - t0, obtainedFvals[t].getEntry(t0) ); // populate fsprime influx residences if well defined.
						}
						
						// obtain the appropriate rate vector (subvector of inRates)
						RealVector q = inRates.getSubVector(t0,  num_disc_times - t0);
						
						// multiply to get g
						g2 = fsp_char.ebeMultiply(q); 
					}
					
					
					////////////////////////////////////////////////////
					
					// obtain the exp( - h2 ) vector by taking a subvector from the full grid version and scaling appropriately
					RealVector eh2 = exp_h2_grid.getSubVector(t0, num_disc_times - t0);
					eh2.mapDivideToSelf(eh2.getEntry(0));
		
					
					// get vector containing the integrand of the integral. 
					RealVector g2_eh2 = g2.ebeDivide(eh2) ;
			
					/////////////////////////////////////////////////////
					
					// integrate along characteristic, and multiply by exp( - h2 )
					RealVector char_f_vals;
					
					char_f_vals = integrateOnChar(g2_eh2, iVal , t0).ebeMultiply(eh2) ; 

					
					// write values into the F_s_values structure
					for(int t = t0 ; t < num_disc_times ; t++) {
						F_s_vals[t].setEntry(t0, char_f_vals.getEntry(t-t0));
					}
			
				
			}
			
		}
		
		/////////////////////////////////////////////////
		
		private RealVector[] create0s() {
			RealVector[] temp = new RealVector[num_disc_times];
			for(int i = 0 ; i < num_disc_times ; i++) {
				temp[i] = new ArrayRealVector(i+1);
			}
			return temp;
		}
		
		
		private void markComputed() {

			// delete references to needless objects, only hang on to final set of values.
			obtainedFvals = null;
			
			computed = true;
			
		}
		
		
		////////////////////////////////////////////////////
		
		void shiftCoordFeed( ArrayList< F_s > f_prob_functions, ArrayList<Integer> f_key) { 
			
			// initialize the holding matrices for the values, all set to 0 //
			ArrayList<RealVector[]> shiftedVals = new ArrayList<RealVector[]>(downStreamStates.size());
			for(int d = 0 ; d < downStreamStates.size() ; d++) {
				shiftedVals.add( create0s() );
			}
			
			// HANDLE CASE WHERE GRID IS FULLY UNDEFINED EXCEPT FOR AT A POINT AT EACH TIME
			if( v_char == n_samples ) {
				// in this case there is 1 point per time for grid, at outer corner, but that gets overwritten by ODE, so just pass zeros
				return;
			}
			
			// setup the scaled grid from which we take subvectors and add a term later to obtain xp and yp grids
			RealVector x0h = t_grid.mapMultiply(n_samples - v_char);
			ArrayList<RealVector> xph = new ArrayList<RealVector>(downStreamStates.size());
						
			// deal with the t = 0 matrix for each downstream state, and add the appropriate scaled grid (xph, yph) into the holder
			for(int d = 0 ; d < downStreamStates.size() ; d++) {
				shiftedVals.get(d)[0] = new ArrayRealVector(1); //t = 0 matrix
				xph.add(t_grid.mapMultiply(n_samples - vChar[downStreamStates.get(d)]));
			}
		
			
			// iterate through each time step, performing the splineInference and Interpolation at each time
			for(int t = 1 ; t < num_disc_times ; t++ ) {
				double ct = t_grid.getEntry(t); // obtain the current time
				
				RealVector xhs = x0h.getSubVector(0, t+1);
				RealVector x0_grid = xhs.mapAdd(ct * v_char);  	// the discrete x values along the characteristics for this state's grid
				Function_Spline interpSlice = new Function_Spline( x0_grid.toArray(), F_s_vals[t].toArray(), -1, false); // 1d spline
				
				// Now the function learning is complete. now map this info onto the downstream grids
				
				for(int d = 0 ; d < downStreamStates.size() ; d++) {
					int target_state = downStreamStates.get(d); // actual state_int of the target
					int vT = vChar[target_state];  // target state vA char value
					

					// if target grid is identical, can ignore the interpolation and just copy over the matrix
					if( v_char == vT) {
						shiftedVals.get(d)[t] = F_s_vals[t].copy();
					}
					
					// interpolation procedure from the inferred function (interpSlice)
					else {
						RealVector target_vec = shiftedVals.get(d)[t] ; // initialize the holder matrix
						
						// obtain the xp and yp grid of x,y vals for the target state
						RealVector xp_grid = xph.get(d).getSubVector(0, t+1).mapAdd(ct * vT);
						
						// interpolate at the correct coordinate value, and write into target_matrix
						for(int i = 0 ; i < t+1 ; i++) {
						
								target_vec.setEntry(i, interpSlice.evaluate_0(xp_grid.getEntry(i))   );
								
						}
					
						
					}
					
					
				}
				
				

				
			}
			
			
			// HAVE SO FAR COMPUTED MATRICES TO PASS DOWNSTREAM. 
			// NEED TO PASS THEM, AND ALSO PASS DOWNSTREAM FOR THE SYMMETRIC CASE
			
			
			// PASS VALUES
			for(int d = 0 ; d < downStreamStates.size() ; d++) {
				int target_si = downStreamStates.get(d);
				F_s copy_target = f_prob_functions.get(f_key.indexOf(target_si)) ;
				
				copy_target.obtainState(state, shiftedVals.get(d));
			}
				
			
		}
		
		void obtainState(int s, RealVector[] shiftedVals) {
			if(necessaryState != s) {
				System.out.println("state " + state + " was incorrectly given shifted-values from state " + s );
				return;
			}
			
				
			obtainedState = true;
			obtainedFvals = shiftedVals;
		}
		
		double[] getCoords(int t, int i){
			
			double time = t_grid.getEntry(t);
			
			double x = t_grid.getEntry(i) * (n_samples - v_char) + time * v_char;
		
			double[] coords = new double[] {time, x};
			
			return coords;
		}

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
	
	
	
	
	
	
//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

	
	

	
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
		
		
		TL_1P_Emissions ep_comp = new TL_1P_Emissions(cr,theta, n_samples, scaled_lengths, 1000, false);

		
		System.out.println("********** \n \n STATE PARTITIONS:\n");
		System.out.println(dlengths);
		System.out.println("\n \n **********");
		System.out.println("************************");
		System.out.println("MARGINAL PDF:");
		System.out.println(ep_comp.getMargPDF());
		System.out.println("\n \n NUMERIC JOINT PDF:");
		System.out.println(ep_comp.getPDF());
		System.out.println("\n \n NUMERIC EMISSION PROBS:");
		System.out.println(helper.pdfToRowProbMatrix(ep_comp.getPDF(), 1));
		System.out.println("\n \n \n ************************");
		//System.out.println("SIMULATION COMPARISON (10^8 RUNS) ");
		//ep_comp.compareSim(100000000);
		
		System.out.println("************************");

		System.out.println("************************");

		

	}

	

	
}
