import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;



public abstract class EmissionsComputer implements FirstOrderDifferentialEquations {
	
	Helper helper = new Helper();
	
	// FIELDS
	
	// demographic parameters
	FunctionUnivariate coalescence_rate;							// coalescence rate
	double mutation_rate;											// mutation rate
			
	// sample and discretization params
	int n_samples;
	double[] d_times; 
	int num_t_bins;
	int num_emissions;

	
	// ode state parameters
	private Hashtable<Integer,List<Integer>> state_map;				// ode state space, maps ints to ordered sets
	private Hashtable<List<Integer>, Integer> rev_state_map;		// maps ordered sets to ints
	int num_states;
	
	// pseudohaploid filter 
	boolean pseudo_hap = false;
	RealMatrix pseudo_hap_filter;     // only gets filled if pseudohap option is selected
	
	// integrator params
	final double t_step_precision = 1e-11;
	
	// values we hold as we compute things
	protected RealMatrix coalWeights;
	protected RealMatrix mutRates;
	
	boolean upToDate = false;
	RealMatrix currentCDF;
	
	
	EmissionsComputer(FunctionUnivariate cr, double mr, int ns, double[] dts, boolean pseudo){
		coalescence_rate = cr.copy();
		mutation_rate = mr;
		n_samples = ns;  
		
		// process pseudo hap status info
		pseudo_hap = pseudo;
		// note: if pseudo_hap is true, n_samples should be multiple of 2 (twice number of pseudohaps)
		if(pseudo_hap && (n_samples % 2 != 0)) {
			System.err.print("number of haplotypes in model not consistent with pseudo_hap switch");System.exit(0);}
		if(pseudo_hap) {this.build_pseudo_hap_filter();}		
	
		
		d_times = dts.clone();
		num_t_bins = d_times.length;
		
		buildStates(n_samples); // the hashmaps and the num_states fields are now filled
		updateCMRates();
		

		
	}
	
	
	
	
	// NEED TO IMPLEMENT IN EACH INSTANCE
	
	abstract void updateCDF();
	
	abstract int[] ancestralSim();
	
		
	
	
	/////////////////////////////////////////////////////////////////////////
	// ODE DEFAULT METHODS
	//////////////////////////////////////////////////////////////////////////
	
	// RETURNS NUMBER OF STATES IN STATE SPACE
	public int getDimension() {
		return num_states;
	}
		
	// COMPUTES DERIVATIVES AND WRITES THEM TO Y_DOT. REQUIRED FOR FODifEQ IMPLEMENTING CLASS
	public void computeDerivatives(double t, double[] y, double[] y_dot)
			throws MaxCountExceededException, DimensionMismatchException {
			
		RealVector temp = new ArrayRealVector(y);
		RealMatrix transition_m = mutRates.add( coalWeights.scalarMultiply( coalescence_rate.evaluate(t) ));
				
		temp = transition_m.preMultiply(temp);
		for(int i = 0 ; i < y_dot.length ; i++) {
			y_dot[i] = temp.getEntry(i);
		}
			
	}
	


	
	
	
	
	//////////////////////////////////////////////////
	
	// DEFAULT 1 POP CRRATE UPDATE AND STATE BUILDING METHODS

	private void updateCMRates() {
		
		coalWeights = new Array2DRowRealMatrix(num_states, num_states);
		mutRates = new Array2DRowRealMatrix(num_states, num_states);

		int[] iS;
		int[] fS;
		
		for(int i = 0 ; i < num_states ; i++) {
			for(int j = 0 ; j < num_states ; j++) {
				iS = intToState(i);
				fS = intToState(j);
				
				if(iS[0] > 1) {  // if iS[0] < 2 it is non-existent or absorbing state, so away rates are all muted	
					
					// coalescence events decrease index0 by 1, leaves index1 alone. 
					if(   (iS[1] == fS[1])  &&  (iS[0] == 1 + fS[0])   ) {
						coalWeights.setEntry(i, j, helper.nChoose2(iS[0]));
					}
					
					// mutation event requires index1 = 0, and writes index0 to both indices of final state
					if(   (iS[1] == 0)   &&   (iS[0] == fS[0])   &&   (fS[0] == fS[1])) {
						mutRates.setEntry(i, j, iS[0]);
					}	
				}

			}
		}
		
		// now need to set diagonals so there is conservation of probability
		RealVector t1 = new ArrayRealVector(num_states, -1.);
				
		coalWeights = coalWeights.add( new DiagonalMatrix( coalWeights.operate(t1).toArray() ) ) ;
		mutRates = mutRates.add( new DiagonalMatrix( mutRates.operate(t1).toArray() ) ) ;
				
		mutRates = mutRates.scalarMultiply(mutation_rate / 2.0); // this is how reco rate is defined
	
	}

	
	private void buildStates(int ns) {
		
		num_states = 0;
		state_map = new Hashtable<Integer, List<Integer>>();
		rev_state_map = new Hashtable<List<Integer>, Integer>();
		
		// go through all possible values for the indices for each state
		for(int k = 1; k < n_samples + 1 ; k ++) {
			for(int ks = 0 ; ks < n_samples + 1 ; ks++) {

					boolean valid_state = true;
						
					// mark state as invalid if it meets conditions
					if(ks < k && ks > 0) {valid_state = false;}
					if(ks == 1 && k == 1) {valid_state = false;}
						
					//
					if(valid_state) {
						state_map.put(num_states , Arrays.asList(k,ks));
						rev_state_map.put( Arrays.asList(k,ks), num_states );
						num_states++;
					}
						
				
			}
		}
		
	
	}
	
	int[] getAbsorbingStateIndices(){
		int[] out = new int[n_samples];
		for(int ks = 1 ; ks < n_samples ; ks++) {
			out[ks] = stateToInt(new int[] {1,ks+1});
		}
		out[0] = stateToInt(new int[] {1,0});
		return out;
	}
	
	private void build_pseudo_hap_filter() {
		
		int model_hap_n = n_samples;       // number of total haps in the model
		int pseudo_hap_n = n_samples / 2 ;  // number of pseudohaploids
		

		pseudo_hap_filter = new Array2DRowRealMatrix(model_hap_n, pseudo_hap_n + 1);
		
		Helper.ChooseComputer chooser = helper.new ChooseComputer(model_hap_n);
		
		// now fill entries of filter matrix. entry i,j is
		// prob that there are j emissions among the n/2 pseudo_haps given that 
		// there are i emissions among the model_haps
		
		for(int i = 0 ; i < pseudo_hap_filter.getRowDimension(); i++) {
			for(int j = 0 ; j < pseudo_hap_filter.getColumnDimension(); j++) {
				
				// entries are hypergeometric 
				double cfactor = 1;
				cfactor = cfactor * chooser.get(model_hap_n - i, pseudo_hap_n - j);
				cfactor = cfactor * chooser.get(i, j);
				cfactor = cfactor / chooser.get1(model_hap_n, pseudo_hap_n);
				
				/*
				cfactor = 1;
				cfactor = cfactor * chooser.get(pseudo_hap_n, j);
				cfactor = cfactor * chooser.get(pseudo_hap_n, i-j);
				cfactor = cfactor / chooser.get1(model_hap_n, i);
				*/
				
				pseudo_hap_filter.setEntry(i, j, cfactor);
			}
		}
		
		pseudo_hap_filter.setColumnVector(0, pseudo_hap_filter.getColumnVector(pseudo_hap_n).add(pseudo_hap_filter.getColumnVector(0)));
		pseudo_hap_filter = pseudo_hap_filter.getSubMatrix(0, pseudo_hap_filter.getRowDimension()-1, 0, pseudo_hap_filter.getColumnDimension()-2);
		
		
		
		// DO FINAL DIMENSION CHECK
		
		if(pseudo_hap_filter.getRowDimension() != 2 * pseudo_hap_filter.getColumnDimension()) {
			System.err.println("PseudoHaploid Issue");
			System.exit(0);
		}

	}
	
	/////////////////////////////////////////////////////////////////////////
	
	// DEFAULT METHODS FOR OBTAINING PROB MATRICES ONCE WE HAVE CDF MATRIX
	
	//////////////////////////////////////////////////////////////////////////
	
	private RealMatrix getCDF() {
		
		if(upToDate) {}
		else {
			double stime = System.nanoTime();

			updateCDF();
			
			System.out.println("^^^ Emission Probs computed: " + Double.toString((System.nanoTime() - stime)/1000000000.) + " seconds");

		}
		return currentCDF.copy();
	}
	
	RealMatrix getPDF() {	
		return cdfToPDF( getCDF() );
	}
	
	RealVector getSFS() {
		RealVector out;
		return helper.mat_colsums(this.getPDF());
	}
	
	
	RealMatrix emissionProbs(boolean binary_emissions) {
		
		// obtain transition probabilities from PDFs.
		RealMatrix out = helper.pdfToRowProbMatrix( getPDF(), 1);
		
		// filter back to pseudo_hap data space if its pseudohaploid data
		if(pseudo_hap) {out = out.multiply(pseudo_hap_filter);}
		
		// if emissions are binary, then use filter to binary probs 
		// (group all columns but first into a single, second column)
		if(binary_emissions) {out = helper.adjustBinaryEmission(out);}
				
		
		return out;
		
	}
	
	RealVector getMargCDF() {
		return helper.mat_rowsums(getCDF());
	}
	
	RealVector getMargPDF(){
		return helper.mat_rowsums(getPDF());
	}
	
	
	RealMatrix cdfToPDF(RealMatrix t_cdf) {
		
		RealMatrix emitCDF = t_cdf;
		RealMatrix emitPDF = new Array2DRowRealMatrix(emitCDF.getRowDimension(), emitCDF.getColumnDimension()); 

		emitPDF.setEntry(0, 0, emitCDF.getEntry(0, 0));
		
		for(int j = 0 ; j < n_samples ; j++) {
			emitPDF.setEntry(0, j, emitCDF.getEntry(0, j));
			for(int i = 1 ; i < num_t_bins + 1 ; i++) {
				emitPDF.setEntry(i, j, emitCDF.getEntry(i, j) - emitCDF.getEntry(i-1, j));
			}
		}

		// now we have the PDFs
		fix(emitPDF, .0000000001);// fix entries that are negative due to numerics. (or arbitrarily close to zero)
				
		return emitPDF;
		
	}
	
	
	///////// private methods /////////////
	
	// MAKE SURE PDF VALUES ARE NON-NEGATIVE, AND NORMALIZE
	private void  fix(RealMatrix m, double epsilon) {
		for(int i = 0 ; i < m.getRowDimension() ; i++) {
			for(int j = 0 ; j < m.getColumnDimension() ; j++) {
				if(m.getEntry(i, j) <=  epsilon) {m.setEntry(i, j, 0.);}
			}
		}
	}
	
	
	/////////////////////////////////////////////////////////////////////////
	// UPDATING DEMOGRAPHY FROM PARAMETERS
	/////////////////////////////////////////////

	abstract void updateDemoFromPs(double[][] params);
	
	abstract double[][] getPsFromDemo();
	
	

	
	/////////////////////////////////////////////////////////////////////////
	
	// SIMULATION METHODS

	//////////////////////////////////////////////////////////////////////////
	
	RealMatrix simulatedCDF(int trials) {
		
		RealMatrix eprobs = new Array2DRowRealMatrix(num_t_bins + 1, n_samples);
		
		int[] temp;
		int t_bin;
		int segs;
		
		for(int t = 0 ; t < trials ; t++ ) {
			temp = ancestralSim();
			
			t_bin = temp[0];
			segs = temp[1];
			
			// following line used for PDF
			//eprobs.setEntry(t_bin, segs, 1 + eprobs.getEntry(t_bin,segs));
			
			// following loops used for CDF
			for(int i1 = t_bin ; i1 < num_t_bins + 1 ; i1++ ) {
					eprobs.setEntry(i1, segs, 1 + eprobs.getEntry(i1,segs));
			}
			
		}
			
		// At this stage eprobs is cdf in length, exact discrete val for seg sites
		return eprobs.scalarMultiply(1./trials);
		
	}

	
	
	
	
	/////////////////////////////////////////////////////////////////////////
	
	// ODE STATE PRINTING, MANIPULATING

	//////////////////////////////////////////////////////////////////////////
	

	
	int stateToInt(int[] state){
		return rev_state_map.get( Arrays.stream( state).boxed().collect( Collectors.toList() ) );
		//return rev_state_map.get(Arrays.asList(state[0],state[1],state[2],state[3]));
	}
		
	int[] intToState(int state){
		List<Integer> temp= state_map.get(state);
		return temp.stream().mapToInt(i->i).toArray();
	}

	void printStateTransitions() {

	}
	
	void printStateKey() {
		
		for(int i = 0 ; i < num_states ; i++) {
			int[] s = intToState(i);
			System.out.print(i + "  corresponds to state   (");
			for(int nn = 0 ; nn < s.length ; nn++) {
				System.out.print( s[nn]+",");
			}
			System.out.println(")");
			
			if (stateToInt(intToState(i)) != i) {
				System.out.println("issue going back and forth from state to int");
				System.exit(0);
			}
		}
	}
	
	double[] getInitialResProbs() {
		// get initialization for this ancestral process in terms of this state vector
		
		double[] out = new double[num_states];
		out[stateToInt(new int[] {n_samples,0})] = 1.;
		return out;
	
	}
	
	
	
	
	
	
}
