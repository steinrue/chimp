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


public abstract class TransitionsComputer implements FirstOrderDifferentialEquations {
	
	Helper helper = new Helper();
	
	// FIELDS
	
	// demographic parameters
	FunctionUnivariate coalescence_rate;							// coalescence rate
	double recombination_rate;										// recombination rate
	boolean sfs_scheme = false;
		
	
	// sample and discretization params
	int n_samples;
	double[] d_times; 
	int num_t_bins;
	
	// ode state parameters
	private Hashtable<Integer,List<Integer>> state_map;				// ode state space, maps ints to ordered sets
	private Hashtable<List<Integer>, Integer> rev_state_map;		// maps ordered sets to ints
	int num_states;
	
	// integrator params
	final double t_step_precision = 1e-11;
	
	// values we hold as we compute things
	protected RealMatrix coalWeights;
	protected RealMatrix recoRates;
	
	boolean upToDate = false;
	RealMatrix currentCDF;
	
	
	TransitionsComputer(FunctionUnivariate cr, double rr, int ns, double[] dts){
		coalescence_rate = cr.copy();
		sfs_scheme = (rr < 0);
		recombination_rate = Math.abs(rr);
		
		
		n_samples = ns;
		d_times = dts.clone();
		num_t_bins = d_times.length;
		
		buildStates(n_samples); // the hashmaps and the num_states fields are now filled
		updateCRRates();
		
		
		
	}
	
	
	
	
	// NEED TO IMPLEMENT IN EACH INSTANCE
	
	abstract void updateCDF();
	
	abstract double[] ancestralSim();
	
	abstract double[] getInitialResProbs();
	
	
	
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
			RealMatrix transition_m = recoRates.add( coalWeights.scalarMultiply( coalescence_rate.evaluate(t) ));
				
			temp = transition_m.preMultiply(temp);
			for(int i = 0 ; i < y_dot.length ; i++) {
				y_dot[i] = temp.getEntry(i);
			}
			
	}
	

	
	
	
	//////////////////////////////////////////////////
	
	// DEFAULT 1 POP CRRATE UPDATE AND STATE BUILDING METHODS
	private void updateCRRates() {
		coalWeights = new Array2DRowRealMatrix(num_states, num_states);		
		recoRates = new Array2DRowRealMatrix(num_states, num_states);		

		int[] iS;
		int[] fS;
		
		for(int i = 0 ; i < num_states ; i++) {
			iS = intToState(i);
			for(int j = 0 ; j < num_states ; j++) {
				fS = intToState(j);
				double cw = 0;
				double rw = 0;
			///////////////	
			// absorbing states
				if((iS[0]==1)   &&   (iS[1]==0)   && (iS[2]==0)) {}
			
			// coalescence possibilities	
				else if(  (fS[0]==iS[0]-1)  &&  (fS[1]==iS[1])  &&  (fS[2]==iS[2])  &&  (fS[3]==iS[3]) ) {
					cw =  helper.nChoose2(iS[0]);
				}
				
				else if( (fS[0]==iS[0])  &&  (fS[1]==iS[1]-1)  &&  (fS[2]==iS[2])  &&  (fS[3]==iS[3])  ) {
					cw =   helper.nChoose2(iS[1]) + iS[0]*iS[1] ;
				}
				
				else if( (fS[0]==iS[0])  &&  (fS[1]==iS[1])  &&  (fS[2]==iS[2]-1)  &&  (fS[3]==iS[3]) ) {
					cw=   helper.nChoose2(iS[2]) + iS[0]*iS[2] ;
				}
				
				else if(  (fS[0]==iS[0]+1)  &&  (fS[1]==iS[1]-1)  &&  (fS[2]==iS[2]-1)  &&  (fS[3]==iS[3])  ) {
					cw =  iS[1]*iS[2] ;
				}
			
			// recombination possibility	
				else if(  (fS[0]==iS[0]-1)  &&  (fS[1]==iS[1]+1)  &&  (fS[2]==iS[2]+1)  &&  (fS[3]==iS[3]+1)  ) {
					rw = iS[0];
				}

				coalWeights.setEntry(i, j, cw);   // for coalescence
				recoRates.setEntry(i, j, rw);   // these are really weights not rates (for reco)
			}
		}
	
		// now need to set diagonals so there is conservation of probability
		RealVector t1 = new ArrayRealVector(num_states, -1.);
		
		coalWeights = coalWeights.add( new DiagonalMatrix( coalWeights.operate(t1).toArray() ) ) ;
		recoRates = recoRates.add( new DiagonalMatrix( recoRates.operate(t1).toArray() ) ) ;
		
		recoRates = recoRates.scalarMultiply(recombination_rate / 2.0); // this is how reco rate is defined
	}

	private void buildStates(int ns) {
		num_states = 0;
		state_map = new Hashtable<Integer, List<Integer>>();
		rev_state_map = new Hashtable<List<Integer>, Integer>();
		
		// go through all possible values for the indices for each state
		for(int kab = 0; kab < ns + 1 ; kab ++) {
			for(int recos = 0 ; recos < 2 ; recos++) {
				for(int ka = 0 ; ka < 2 ; ka++) {
					for(int kb = 0 ; kb < 2 ; kb++) {
						boolean valid_state = true;
						
						// mark state as invalid if it meets conditions
						if(ka > recos || kb > recos) {valid_state = false;}
						if(ka + kab < 1 || kb+kab < 1 ) {valid_state = false;}
						if(ka + kab > ns || kb+kab > ns) {valid_state=false;}
						
						// this one is empirically invalid for max 1 recombination
						if(kab == 0){valid_state = false;}
						
						
						//
						if(valid_state) {
							state_map.put(num_states , Arrays.asList(kab,ka,kb,recos));
							rev_state_map.put( Arrays.asList(kab,ka,kb,recos), num_states );
							num_states++;
						}
						
					}
				}
			}
		}
	}
	
	
	
	/////////////////////////////////////////////////////////////////////////
	
	// DEFAULT METHODS FOR OBTAINING PROB MATRICES ONCE WE HAVE CDF MATRIX
	
	//////////////////////////////////////////////////////////////////////////
	
	RealMatrix getCDF() {
		
		if(upToDate) {}
		else {
			double stime = System.nanoTime();
			
			updateCDF();
		
			System.out.println(">>> Transition Probs computed: " + Double.toString((System.nanoTime() - stime)/1000000000.) + " seconds");

		}
		return currentCDF.copy();
	}
	
	RealMatrix getPDF() {
		return cdfToPDF( getCDF() );
	}

	RealMatrix transitionProbs() {
		
		// obtain transition probabilities from PDFs.
		RealMatrix out = getPDF();
		
		// this is only needed if we want SFS model where transition probability always reverts to marginal (infite recombination limit)
		if(sfs_scheme) {
			RealVector marg = helper.mat_rowsums(out);
			for(int i = 0 ; i < out.getRowDimension(); i++) {
				out.setRowVector(i, marg);
			}
		}
		
		
		
		else {
			out = helper.pdfToRowProbMatrix( out, 0);
			for(int i = 0; i < out.getRowDimension(); i++) {
				if(out.getEntry(i, i) == 1.0) {
					System.out.println("-- Warning, locked state in CHMM transition --");
					i = out.getRowDimension();
				}
			}
		}
		
		
		return out;
		
	}

	RealMatrix cdfToPDF(RealMatrix t_cdf) {
		
		RealMatrix jointCDF = t_cdf;
		RealMatrix jointPDF = new Array2DRowRealMatrix(jointCDF.getRowDimension(), jointCDF.getColumnDimension()); 

		jointPDF.setEntry(0, 0, jointCDF.getEntry(0, 0));
		for(int i = 1; i < jointCDF.getRowDimension() ; i++ ) {
			jointPDF.setEntry(0, i, jointCDF.getEntry(0, i) - jointCDF.getEntry(0, i-1));
			jointPDF.setEntry(i, 0, jointCDF.getEntry(i, 0) - jointCDF.getEntry(i-1, 0));

			for(int j = 1; j < jointCDF.getColumnDimension() ; j++) {
				// draw the rectangular regions chart, and you see that PDFs are obtained from CDFs as follows.
				jointPDF.setEntry(i, j, jointCDF.getEntry(i, j) - jointCDF.getEntry(i, j-1) - jointCDF.getEntry(i-1, j) + jointCDF.getEntry(i-1, j-1) ); 
			}
		}

		// now we have the PDFs
		fix(jointPDF, .0000000001);// fix entries that are negative due to numerics. (or arbitrarily close to zero)
				
		return jointPDF;
		
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
		
		RealMatrix tprobs = new Array2DRowRealMatrix(num_t_bins + 1, num_t_bins+1);
	
		double[] temp;
		int i;
		int j;
		
		for(int t = 0 ; t < trials ; t++ ) {
			temp = ancestralSim();
			
			i = helper.placeInPartition(temp[0], d_times); 
			j = helper.placeInPartition(temp[1], d_times);
			
			// following line for PDF
			//tprobs.setEntry(i, j, 1 + tprobs.getEntry(i,j));
			
			// following loops used for CDF
			for(int i1 = i ; i1 < num_t_bins + 1 ; i1++ ) {
				for(int j1 = j ; j1 < num_t_bins + 1 ; j1++) {
					tprobs.setEntry(i1, j1, 1 + tprobs.getEntry(i1,j1));
				}
			}
			
			
		}
		
		
		
		// At this stage tprobs is really just the joint pdf (scaled by trials) for t1mrca and t2mrca
		
		return tprobs.scalarMultiply(1./trials);
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
		for(int i = 0 ; i < num_states;i++) {
			for(int j = 0 ; j < num_states; j++) {
				double cw = coalWeights.getEntry(i, j);
				double rw = recoRates.getEntry(i, j) * 2 / recombination_rate;
				
				if(cw+rw == 0.) {}
				else {
					System.out.print(Arrays.toString(intToState(i)) + " flow to --->");
					System.out.print(Arrays.toString(intToState(j)) + " is ::");
					System.out.println("c"+cw+ " - r"+ rw);
				}
				
			}
		}
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
	
	
	
	
	
	
	
	
}
