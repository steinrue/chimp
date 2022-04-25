import java.io.FileWriter;
import java.io.IOException;
import java.util.AbstractMap;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MathUnsupportedOperationException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.LineSearch;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;


public abstract class EM_HMM {
	
	Helper helper = new Helper();

	
	// Number of Samples
	int[] num_samples; // list of all n_s values in the data we will be looking at. for pseudohaps this is # of pseudohaps, NOT 2*n (ie, not model_n)
	
	// HMMs used. 
	protected HMmodel hmm_model;  // current state of the EM process
	private HMmodel calc_hmm;   // a dummy model used for maximization computations
	
	// regularization coefficients
	private RealVector regularization_coefficients = new ArrayRealVector(4);
	
	
	//counter tracking number of EM steps
	int em_iter; 
	
	//optimizer
	FlexibleOptimizer optim; 
	
	// data Handler
	Helper.dataStreamsHandler dataStreams = null;
	
	//// size of metalocus
	int metaLocus_size; 
	

	
	// Output holders for FB algorithm aggregated over all the chromosomes/subsampling iterations. 
	// expectations for each n_s model are different entry in map (each entry is array of {transExp, emitExp, startExp} )
	boolean fb_output_upToDate;
	AbstractMap<Integer,RealMatrix[]> current_expectation_matrices = new HashMap<Integer,RealMatrix[]>();
	double current_ll = -1;
	
	// fields for printing log file //
	Random random = new Random();
	int log_earmark = random.nextInt(100000000);
	boolean print_log = false;
	FileWriter scribe = null;
	
	
	
	/////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////
	
	
	EM_HMM(FunctionUnivariate cr, double rr, double mr, int[] ns, AbstractMap<Integer,double[]> partitions, String genealogy_rep, int mL_size,int[] additional_params, boolean binary_emission, boolean pseudo){
		num_samples = ns;
		metaLocus_size = mL_size;
		
		
		//initialize hmm models with the appropriate parameters
		hmm_model = new HMmodel(cr, rr, mr, ns, partitions, genealogy_rep, additional_params, binary_emission, pseudo); 
		calc_hmm = new HMmodel(cr, rr, mr, ns, partitions, genealogy_rep, additional_params, binary_emission, pseudo);
		
		// check switch for printing optimization log, and create file if needed.
		if(print_log) {
			try { scribe = new FileWriter( "log_" + log_earmark + ".csv");  }
			catch (IOException e) { e.printStackTrace(); }
		}
		
	}


	///////////////////////////////////////
	// ABSTRACT METHODS, MUST IMPLEMENT ///
	///////////////////////////////////////
	
	// for these methods, do not convert to log values. we do that later. just go from demo to 
	// actual inferencing values. 
	
	abstract protected double[] params_From_Hmm(HMmodel hmm);
	
	abstract protected void params_To_HMM(double[] params, HMmodel hmm);
	
	
	///////////////////////////////////////////////////////////////////////////////
	///// METHODS WE CALL FROM OTHER PROGRAMS/FILES /////
	///////////////////////////////////////////////////////////////////////////////
	
	void print_marginals() {
		
		for(int i = 0 ; i < num_samples.length; i++) {
			int t_ns = num_samples[i];
			System.out.println("Marginal for n_s = " + t_ns+ " -> "+ hmm_model.getSS(t_ns));
		}		
	}

	double[] get_current_params() {
		return params_From_Hmm(hmm_model);
	}
	
	public void update_reg_weights(double[] lambdas) {
		regularization_coefficients = new ArrayRealVector(lambdas);
	}

	
	
	/////// use these to load list of vcf's either to streams in RAM or to files in specified folder
	
	void loadData( Helper.dataStreamsHandler dssHandler) {
		
		dataStreams = dssHandler;
		this.fb_clear_output();
		em_iter = 0;
	}
	
	
	// Perform E step and M step and update. 
	void em_single_update(int m_step_eval_limit, double simplex_span, int simplex_type) {
		
		System.out.println("\n ----------------- Begin EM --------------------");
		em_iter++;
		// E- STEP ////////////////////////////////////////
		///////////////////////
		
		// perform forward backwards only if we havent already
		if(!fb_output_upToDate) {this.fb_overAllData();}
		// gets expectations, and updates  current_expectation_matrices
		
				
		// M- STEP ////////////////////////////////////////
		///////////////////////
		double stime = System.nanoTime();

		
		// Prepare the optimization parameters. The objective function (to be maxed) and initial guess
		
		Q_Function q_f =  new Q_Function(current_expectation_matrices, calc_hmm);
		double[] initial_guess = params_From_Hmm(hmm_model);
		
		// READY OPTIMIZER
		optim = new FlexibleOptimizer(simplex_type, simplex_span, m_step_eval_limit, q_f);

		// PERFORM OPTIMIZATION 
		double[] updated_vals = optim.maximize(initial_guess, em_iter);
		
		
		// UPDATE HMM using the params obtained from optimization ( exponentiate if needed)
		params_To_HMM(updated_vals, hmm_model); 	// update smc appropriately from the optimized params
		this.fb_clear_output(); // resetting params in hmm_model means fb algorithm outputs are no longer up to date
		
		System.out.println("\n \n**************************\n \n M step time: " + Double.toString((System.nanoTime() - stime)/1000000000.) +" seconds \n\n ***************************** \n\n");
		System.gc();
	}

	
	double get_posteriorLL() {
		
		if (!fb_output_upToDate) { 
			this.fb_overAllData();
		}

		return current_ll;

	}

	
	// clear all temp csv files from the disc
	void clear_all_csvs() {
			
			// close the log file //
			if(print_log) {
				try { 	scribe.flush(); scribe.close(); } catch (IOException e) { e.printStackTrace(); 	}
				if(!print_log) {helper.delete_csvFromDisc("log_" + log_earmark+ ".csv"); }
			}
			
			// clear the temporary stream files if they exist
			dataStreams.deleteStreamFiles();
			
			
			
	}
	
	// clears the forward backwards algorithm information (most recent expectation matrices and the most recent LL)
	private void fb_clear_output() {
		fb_output_upToDate = false;
		current_expectation_matrices.clear();
		current_ll = -1;
	}

	// performs forward backwards using the Current model_hmm values and the data loaded in currently //
	// does this from RAM or Files depending on streams_in_mem boolean
	// after fb, writes the ouput expectation matrices and ll to the appropriate holders
	
	protected void fb_overAllData() {
		
		// make sure data is loaded appropriately
		if( dataStreams == null ) {System.out.println("No Data Loaded for EM_HMM"); System.exit(0);}
				
		double stime = System.nanoTime();

		// clear all values, so expectations+likelihoods are counted from scratch
		current_expectation_matrices.clear();		
		double totalLL = 0 ;
		
		
		for(int f = 0 ; f < dataStreams.num_streams; f++) {
			
			// get the stream and its associated n_s value
			int n_s = dataStreams.get_ns(f);
			int[][] dstream = dataStreams.getStream(f);
			
			// if n_s is a new value, initialize it.
			if(!current_expectation_matrices.containsKey(n_s)) {
				int[] mat_dims = new int[] {hmm_model.getTProbs(n_s).getRowDimension(),hmm_model.getEProbs(n_s).getColumnDimension()};
				
				RealMatrix tE = new Array2DRowRealMatrix(mat_dims[0],mat_dims[0]);  
				RealMatrix eE = new Array2DRowRealMatrix(mat_dims[0],mat_dims[1]); 	
				RealMatrix sE = new Array2DRowRealMatrix(1,mat_dims[0]);
				
				current_expectation_matrices.put(n_s, new RealMatrix[] {tE,eE,sE});
			}
			
			// compute expectations for this stream only
			ExpectationsComputer exp_computer = new Expectation_FBskip(hmm_model.getTProbs(n_s), hmm_model.getEProbs(n_s), hmm_model.getSS(n_s));
			exp_computer.loadDataFromStream(dstream, metaLocus_size);
			exp_computer.computeExps();
			
			// add the expectations to the total for the appropriate n_s
			RealMatrix[] temp = current_expectation_matrices.get(n_s);
			temp[0] = temp[0].add(exp_computer.getExpT());
			temp[1] = temp[1].add(exp_computer.getExpE());
			temp[2].setRowVector(0, temp[2].getRowVector(0).add(exp_computer.getExpS()));
			
			current_expectation_matrices.put(n_s, temp);
			
			totalLL = totalLL + exp_computer.get_LL();
			
		}
		
		current_ll = totalLL;
		
		fb_output_upToDate = true;
		
		System.out.println("\n \n**************************\n \n E step time: " + Double.toString((System.nanoTime() - stime)/1000000000.) +" seconds \n\n ***************************** \n\n");

	}

	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	class HMmodel {
		
		private boolean binary_emission;
		private boolean pseudohap;
				
		// we need a separate computer for each ns value, ns is key in map
		private AbstractMap<Integer,TransitionsComputer> tPcomp;
		private AbstractMap<Integer,EmissionsComputer> ePcomp;
		
		// we compute one of each type of matrix for each value of ns. ns is the key in following maps
		private AbstractMap<Integer,RealMatrix> transProbs; // trans prob per metalocus (TP per locus raised to mlsize power)
		private AbstractMap<Integer,RealMatrix> emitProbs;
		private AbstractMap<Integer,RealVector> ssDistro;
		
				
		private double[][] demo_ps; //current demographic parameters
		
		HMmodel(FunctionUnivariate cr, double rr, double mr, int[] ns, AbstractMap<Integer,double[]> ts, String genealogy_rep, int[] additional_params, boolean bin_em, boolean pseudo){
			
			binary_emission = bin_em;
			pseudohap = pseudo;
			
			int[] keys = num_samples;
			
			// initialize all the map holders
			tPcomp = new HashMap<Integer, TransitionsComputer>();
			ePcomp = new HashMap<Integer, EmissionsComputer>();
			transProbs = new HashMap<Integer, RealMatrix>();
			emitProbs = new HashMap<Integer, RealMatrix>();
			ssDistro = new HashMap<Integer, RealVector>();

			
			// go through each key, create tP, eP computers for each
			for(int i = 0 ; i < keys.length ; i++) {
				int n_s = keys[i]; // num of haplotypes in tree for model
				int model_n = n_s * (pseudo ? 2 : 1);
				
//				if(genealogy_rep.equals("TH_1P")) {
				// can we hack this in here?
				if(genealogy_rep.equals("TH_1P") || (model_n == 2)) {
//					if(additional_params != null) {System.out.print("Invalid additional params"); System.exit(0);}
					if (genealogy_rep.equals("TH_1P") && (additional_params != null)) {
						System.out.print("Invalid additional params"); System.exit(-1);
					}
					tPcomp.put(n_s, new TH_1P_Transitions(cr, rr, model_n, ts.get(n_s)));
					ePcomp.put(n_s, new TH_1P_Emissions(cr, mr, model_n, ts.get(n_s), pseudohap));
				}
				
				else if (genealogy_rep.equals("TL_1P")) {
					int[] ap = additional_params;

					if(additional_params.length == 2) { tPcomp.put(n_s, new TL_1P_Transitions(cr, rr, model_n, ts.get(n_s), ap[1])); }
					else if(additional_params.length == 3) { tPcomp.put(n_s, new TL_1P_Transitions(cr, rr, model_n, ts.get(n_s), new int[] {ap[1],ap[2]})); }
					else {System.out.print("Invalid additional params"); System.exit(0);}
					
					ePcomp.put(n_s, new TL_1P_Emissions(cr, mr, model_n, ts.get(n_s), ap[0], pseudohap));

				}
				
				else {System.out.println("Need a genealogy representation specified"); System.exit(0);}
				
				
				int mls = metaLocus_size;
				transProbs.put(n_s, tPcomp.get(n_s).transitionProbs().power(mls)); // metalocus exponentiation
				emitProbs.put(n_s, ePcomp.get(n_s).emissionProbs(binary_emission));
				ssDistro.put(n_s, ePcomp.get(n_s).getMargPDF()); // obtain the initial distro from pdf

			}
			
			
			// initialize demography
			demo_ps = tPcomp.get(num_samples[0]).getPsFromDemo(); // pulls the vals that were fed to the trans computer. missing mutation still
			demo_ps[2] = ePcomp.get(num_samples[0]).getPsFromDemo()[2]; // pulls mutation rate from emission computer
			
		}
		
		
		////////////////////////////////////////////////////////
		// ACCESS METHODS 
		
		// update the split time. pick new split time, compute emission and transition probs and write it in 
		void updateDemography(double[][] dem_params) {
			
			// params is list of params. 
			// index 0 is the params needed to update the coalescence_rates in each
			// index 1 is the rr, and index 2 is the mr
			
			demo_ps = dem_params.clone();
						
			int[] keys = num_samples;
			
			for(int i = 0 ; i < keys.length; i++) {
				int key = keys[i];
				tPcomp.get(key).updateDemoFromPs(demo_ps);
				ePcomp.get(key).updateDemoFromPs(demo_ps);
				System.out.println("-----------------------");
				
				int mls = metaLocus_size;
				transProbs.put(key, tPcomp.get(key).transitionProbs().power(mls));
				emitProbs.put(key, ePcomp.get(key).emissionProbs(binary_emission));	
				ssDistro.put(key, ePcomp.get(key).getMargPDF()) ;
			}

			System.out.println("-----------------------");
			System.out.println("-----------------------");
			
		



			

			
		}
		
		double[][] getDemoParams(){
			return demo_ps;
		}
		
		RealMatrix getTProbs(int ns) {
			
			return transProbs.get(ns).copy();
		}
		
		RealMatrix getEProbs(int ns) {
			return emitProbs.get(ns).copy();
		}
		
		RealVector getSS(int ns) {
			return ssDistro.get(ns).copy();
		}
		
		RealVector getRegularization(double estimated_cr) {
			// array of regularizing quantities
			// [d0l2, d1l1, d1l2, d2l2]
			
			RealVector out = new ArrayRealVector(4);
			FunctionUnivariate t_crate = tPcomp.get(tPcomp.keySet().toArray()[0]).coalescence_rate;
			
			out.setEntry(0, t_crate.reg_d0l2(estimated_cr));
			out.setEntry(1, t_crate.reg_d1l1());
			out.setEntry(2, t_crate.reg_d1l2());
			out.setEntry(3, t_crate.reg_d2l2());

			return out;
		}
		
	}
	
	// maximization function class, constructed using the expectation matrices
	public class Q_Function implements MultivariateFunction{
		
		PointValuePair best_eval = new PointValuePair(new double[] {}, Double.NEGATIVE_INFINITY);		
		
		AbstractMap<Integer,RealMatrix[]> tesX0; // all the transitions, emissions, starts, expectations for all the different n_s values used

		private HMmodel hmC;
				
		private double e_CR; // estimated effective coalescent rate from Wattersons estimator
		
		Q_Function(AbstractMap<Integer,RealMatrix[]> x0, HMmodel hm){
				tesX0 = x0;

				hmC = hm;
					
				e_CR = hm.getDemoParams()[2][0] / dataStreams.get_Watterson();				
		}
				
		public double value(double[] log_args) {
					System.gc();

			
					params_To_HMM(log_args, hmC);
					
					System.out.println("#PARAMS: " + Arrays.toString(log_args));

					double total_ll = 0;
					
					for(int key = 0 ; key < num_samples.length ; key++) {
						int n_s = num_samples[key]; // iterate over each n_s value
						
						// these are hmm probs for the given n_s value
						RealMatrix tP = hmC.getTProbs(n_s);
						RealMatrix eP = hmC.getEProbs(n_s);
						RealVector sP = hmC.getSS(n_s);
						
						
						RealMatrix[] expectations = tesX0.get(n_s);
						RealMatrix tX = expectations[0];
						RealMatrix eX = expectations[1];
						RealVector sX = expectations[2].getRowVector(0);
						
						
						double ll_ns = 0;
						
						// add the transition contributions
						for(int i = 0 ; i < tP.getRowDimension(); i++) {
							for(int j = 0 ; j < tP.getColumnDimension() ; j++) {
								if(tX.getEntry(i, j) > 0. && tP.getEntry(i, j) > 0.) {ll_ns = ll_ns + Math.log(tP.getEntry(i, j))*(tX.getEntry(i, j));}
							}
						}
						
						// add the emission contributions.
						for(int i = 0 ; i < eP.getRowDimension(); i++) {
							for(int j = 0 ; j < eP.getColumnDimension() ; j++) {
								if(eX.getEntry(i, j) > 0. && eP.getEntry(i, j) > 0.) { ll_ns = ll_ns + Math.log(eP.getEntry(i, j))*(eX.getEntry(i, j)); }
							}
						}
						
						// add initial distro contribution
						for(int i = 0 ; i < sP.getDimension() ; i++) {
							if(sX.getEntry(i) > 0. && sP.getEntry(i) > 0.) { ll_ns = ll_ns + Math.log(sP.getEntry(i))*(sX.getEntry(i));}
						}
						
						// ll_ns is the contribution of this n_s's trans, emit, start matrices to the log likelihood
						
						total_ll = total_ll + ll_ns;
					}
					
					// total_ll is ll aggregate over all three matrices for each n_s value
					 
					
					if(total_ll == 0.0) {
						System.out.println("-- Warning, Q-function forced value --");
						total_ll = Double.NEGATIVE_INFINITY;
					}
					
					// take care of regularization contribution
					RealVector reg_qs = hmC.getRegularization(e_CR);
					
					// Scale reg by total genetic length across all streams.
					reg_qs.mapMultiplyToSelf(dataStreams.getGenLen());
					
					System.out.println("(Q_eval:"+total_ll +" ; Reg:"+ reg_qs+")");
					
					if(print_log) { try{scribe.append(total_ll + ", ");
						double[] tarr = helper.ebeExp(log_args);
						for(int ii = 0 ; ii < tarr.length;ii++) {
							scribe.append(tarr[ii] + ", ");
						} scribe.append("\n");
					}catch (IOException e){;}}
					
				
					double 	qval = total_ll - regularization_coefficients.dotProduct(reg_qs);
					
					System.out.println("(Q_eval (adjusted):"+qval+")");

					// keep record of the best evaluation for this Q-function
					if(qval > best_eval.getValue()) { best_eval = new PointValuePair(log_args, qval);}
					
					return qval;
		}
		
	}
	
	
	private class ApacheNelderMeadSimplex_custom extends NelderMeadSimplex{
		
		double simplex_initialization_scale;
		int simplex_type;
		MultivariateFunction qfunc_sbuild;
		
		public ApacheNelderMeadSimplex_custom( int s_type, double s_step, int space_dim, MultivariateFunction q){
			super(space_dim);
			simplex_type = s_type;
			simplex_initialization_scale = s_step;
			qfunc_sbuild = q; // this is the objective function used only for building the equilateral+-hybrid simplex
		}
		

		public void build(final double[] startPoint) {
			super.build(startPoint);
			if (this.getDimension() != startPoint.length) {
				     throw new DimensionMismatchException(this.getDimension(), startPoint.length);	        
			}
			
			if(simplex_type == 0 ) {build_EquilateralAdjust(startPoint);} // returns initial guess, also reinitializes n_m_s appropriately
			else if(simplex_type == 1) {build_Equilateral( startPoint);} //""
			else if(simplex_type == 2) {build_EquiHybrid( startPoint);} //""
			else if(simplex_type == 3) {build_default(startPoint);} //""
			else if(simplex_type == 4) {build_AltDefault(startPoint);} //""
			// currently type5 is used by Powell (look in EM step to see this)
			else if(simplex_type == 6) {build_EquilateralRandom( startPoint);} //""

		
			
			else {System.err.println("Invalid Simplex Type Provided"); System.exit(0);} // no valid option selected
		
		}
		
		private void build_default(final double[] startPoint) {
			System.out.println("Build default (apache default) simplex");
			int dimension = this.getDimension();
			PointValuePair[] simplex;
			
			simplex = new PointValuePair[dimension + 1];
			simplex[0] = new PointValuePair(startPoint, Double.NaN);
			
			
			// Set remaining vertices.
			for (int i = 0; i < dimension; i++) {
			      final double[] vertexI = new double[dimension];
			      for (int k = 0; k < dimension; k++) {
			    	  if(k <= i) { vertexI[k] = startPoint[k] + simplex_initialization_scale; }
			    	  else {vertexI[k] = startPoint[k]; }
			      }
			      
			      simplex[i + 1] = new PointValuePair(vertexI, Double.NaN);
			}
			
			this.setPoints(simplex);
		}
		
		private void build_AltDefault(final double[] startPoint) {
			System.out.println("Build alternating default simplex");

			
			int dimension = this.getDimension();
			PointValuePair[] simplex;
			
			simplex = new PointValuePair[dimension + 1];
			simplex[0] = new PointValuePair(startPoint, Double.NaN);
			
			
			// Set remaining vertices.
			for (int i = 0; i < dimension; i++) {
			      final double[] vertexI = new double[dimension];
			      for (int k = 0; k < dimension; k++) {
			    	  if(k <= i) { vertexI[k] = startPoint[k] + simplex_initialization_scale * Math.pow(-1, k); }
			    	  else {vertexI[k] = startPoint[k]; }
			      }
			      
			      simplex[i + 1] = new PointValuePair(vertexI, Double.NaN);
			}
			this.setPoints(simplex);
		}
		
		private void build_Equilateral(final double[] startPoint) {
			System.out.println("Build equilateral simplex");

			
			// s is scale of size (of start simplex), it is edge length / sqrt(2)
			int d = this.getDimension();
			double s = simplex_initialization_scale;
			
			// CREATE THE START SIMPLEX (here is equilateral method, first point located at origin, rest all in positive n-drant)
			
			double[][] s_simp = new double[d + 1][d];
			
			// first point is appropriate same displacement along all dimensions
			for(int ii = 0 ; ii < s_simp[0].length; ii++) { s_simp[0][ii] = (1 - Math.sqrt(1 + d)) * s/d; }
			// rest points are each just step of size s along single dimension
			for(int i = 0 ; i < d; i++) { 	s_simp[i+1][i] = s_simp[i+1][i] + s; }
			
			RealMatrix start_simplex = new Array2DRowRealMatrix(s_simp);

			// REPOSITION SO CENTROID IS OVERLAID ON ORIGIN
			
			// now find the centroid
			int num_points = d+1;
			RealVector centroid = new ArrayRealVector(num_points,1.0);
			centroid = start_simplex.preMultiply(centroid); centroid.mapDivideToSelf(num_points);
			
			// reposition to match centroid to origin
			for(int i = 0 ; i < start_simplex.getRowDimension() ; i++) {
				start_simplex.setRowVector(i, start_simplex.getRowVector(i).subtract(centroid));
			}
			
			// initial guess is the startPoint 
			RealVector ig = new ArrayRealVector(startPoint);
			
		
			// now using this start simplex with the ig initial guess should place the centroid origin on this start point
			PointValuePair[] simplex = new PointValuePair[d + 1];
			for(int i = 0 ; i < simplex.length ; i++) {
				double[] vertexI = ig.add(start_simplex.getRowVector(i)).toArray();
				simplex[i] = new PointValuePair(vertexI, Double.NaN);
			}
			
			this.setPoints(simplex);
		}
		
		private void build_EquiHybrid(final double[] startPoint) {
			System.out.println("Build hybrid of equilateral and negative-equilateral simplex");

			
			// s is scale of size (of start simplex), it is edge length / sqrt(2)
			int d = this.getDimension();
			double s = simplex_initialization_scale;

			// CREATE POSITIVE EQUILATERAL SIMPLEX
			double[][] s_simp = new double[d + 1][d];
			
			// first point is appropriate same displacement along all dimensions
			for(int ii = 0 ; ii < s_simp[0].length; ii++) { s_simp[0][ii] = (1 - Math.sqrt(1 + d)) * s/d; }
			// rest points are each just step of size s along single dimension
			for(int i = 0 ; i < d; i++) { 	s_simp[i+1][i] = s_simp[i+1][i] + s; }
			
			RealMatrix start_simplex = new Array2DRowRealMatrix(s_simp);
			
			// REPOSITION SO CENTROID IS OVERLAID ON ORIGIN
			
			// now find the centroid
			int num_points = d+1;
			RealVector centroid = new ArrayRealVector(num_points,1.0);
			centroid = start_simplex.preMultiply(centroid); centroid.mapDivideToSelf(num_points);
						
			// reposition to match centroid to origin
			for(int i = 0 ; i < start_simplex.getRowDimension() ; i++) {
				start_simplex.setRowVector(i, start_simplex.getRowVector(i).subtract(centroid));
			}
			
			//NOW HAVE POSITIVE EQUILATERAL SIMPLEX
			////////////////////////////////////////////////////////////
			// initial guess is the startPoint 
			RealVector ig = new ArrayRealVector(startPoint);
			
			///////////////////////////////////////////////////////////
			
			// add positive and negative equilateral simplex to a list to keep track
			PointValuePair[] points = new PointValuePair[2 * start_simplex.getRowDimension()];
			for(int i = 0 ; i < start_simplex.getRowDimension() ; i++) {
				double[] p;
				
				// add positive simplex version of vertex
				p = ig.add(start_simplex.getRowVector(i)).toArray();
				points[i] = new PointValuePair(p, qfunc_sbuild.value(p));
				// add negative simplex version of vertex
				p = ig.subtract(start_simplex.getRowVector(i)).toArray();
				points[i + start_simplex.getRowDimension()] = new PointValuePair(p, qfunc_sbuild.value(p));
			}
			
			//////////////////////////////////////////////////////////////
			// setup comparator and arrange LL in order from best (highest) to worst (lowest)
			final Comparator<PointValuePair> pvp_comp  = new Comparator<PointValuePair>() {
			            public int compare(final PointValuePair o1, final PointValuePair o2) {
			                final double v1 = o1.getValue();
			                final double v2 = o2.getValue();
			                return Double.compare(v2, v1);
			            }
			 };
			 // rearrange points in proper order
			 Arrays.sort(points, pvp_comp);
			 
			//////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////

			 // only keep the best points
			 PointValuePair[] simplex = new PointValuePair[d+1];
			 for(int i = 0 ; i < simplex.length ; i++) {
				 simplex[i] = points[i];
			 }
			 
			 this.setPoints(simplex);
			 
			 
		}
		
		private void build_EquilateralAdjust(final double[] startPoint) {
			System.out.println("Build equilateral simplex, and replace worst point by origin if its better");

			
			// s is scale of size (of start simplex), it is edge length / sqrt(2)
			int d = this.getDimension();
			double s = simplex_initialization_scale;

			// CREATE POSITIVE EQUILATERAL SIMPLEX
			double[][] s_simp = new double[d + 1][d];
			
			// first point is appropriate same displacement along all dimensions
			for(int ii = 0 ; ii < s_simp[0].length; ii++) { s_simp[0][ii] = (1 - Math.sqrt(1 + d)) * s/d; }
			// rest points are each just step of size s along single dimension
			for(int i = 0 ; i < d; i++) { 	s_simp[i+1][i] = s_simp[i+1][i] + s; }
			
			RealMatrix start_simplex = new Array2DRowRealMatrix(s_simp);
			
			// REPOSITION SO CENTROID IS OVERLAID ON ORIGIN
			
			// now find the centroid
			int num_points = d+1;
			RealVector centroid = new ArrayRealVector(num_points,1.0);
			centroid = start_simplex.preMultiply(centroid); centroid.mapDivideToSelf(num_points);
						
			// reposition to match centroid to origin
			for(int i = 0 ; i < start_simplex.getRowDimension() ; i++) {
				start_simplex.setRowVector(i, start_simplex.getRowVector(i).subtract(centroid));
			}
			
			//NOW HAVE POSITIVE EQUILATERAL SIMPLEX
			////////////////////////////////////////////////////////////
			// initial guess is the startPoint 
			RealVector ig = new ArrayRealVector(startPoint);
			
			///////////////////////////////////////////////////////////
			
			// add simplex vertices to a list to keep track
			PointValuePair[] points = new PointValuePair[start_simplex.getRowDimension() + 1];
			for(int i = 0 ; i < start_simplex.getRowDimension() ; i++) {
				double[] p;
				
				// add vertex to list
				p = ig.add(start_simplex.getRowVector(i)).toArray();
				points[i] = new PointValuePair(p, qfunc_sbuild.value(p));
				
			}
			// add the origin as well as a final point to compare to
			double[] p = ig.toArray();
			points[start_simplex.getRowDimension()] = new PointValuePair(p, qfunc_sbuild.value(p));
			
			//////////////////////////////////////////////////////////////
			// setup comparator and arrange LL in order from best (highest) to worst (lowest)
			final Comparator<PointValuePair> pvp_comp  = new Comparator<PointValuePair>() {
			            public int compare(final PointValuePair o1, final PointValuePair o2) {
			                final double v1 = o1.getValue();
			                final double v2 = o2.getValue();
			                return Double.compare(v2, v1);
			            }
			 };
			 // rearrange points in proper order
			 Arrays.sort(points, pvp_comp);
			 
			//////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////

			 // only keep the best points
			 PointValuePair[] simplex = new PointValuePair[d+1];
			 for(int i = 0 ; i < simplex.length ; i++) {
				 simplex[i] = points[i];
			 }
			 
			 this.setPoints(simplex);
			 
			 
		}
		
		private void build_EquilateralRandom(final double[] startPoint) {
						System.out.println("Build equilateral simplex, but rotate it randomly");

			
						// s is scale of size (of start simplex), it is edge length / sqrt(2)
						int d = this.getDimension();
						double s = simplex_initialization_scale;

						// CREATE POSITIVE EQUILATERAL SIMPLEX
						double[][] s_simp = new double[d + 1][d];
						
						// first point is appropriate same displacement along all dimensions
						for(int ii = 0 ; ii < s_simp[0].length; ii++) { s_simp[0][ii] = (1 - Math.sqrt(1 + d)) * s/d; }
						// rest points are each just step of size s along single dimension
						for(int i = 0 ; i < d; i++) { 	s_simp[i+1][i] = s_simp[i+1][i] + s; }
						
						RealMatrix start_simplex = new Array2DRowRealMatrix(s_simp);
						
						// REPOSITION SO CENTROID IS OVERLAID ON ORIGIN
						
						// now find the centroid
						int num_points = d+1;
						RealVector centroid = new ArrayRealVector(num_points,1.0);
						centroid = start_simplex.preMultiply(centroid); centroid.mapDivideToSelf(num_points);
									
						// reposition to match centroid to origin
						for(int i = 0 ; i < start_simplex.getRowDimension() ; i++) {
							start_simplex.setRowVector(i, start_simplex.getRowVector(i).subtract(centroid));
						}
						
						//NOW HAVE POSITIVE EQUILATERAL SIMPLEX
						////////////////////////////////////////////////////////////
						
						//////// ROTATE BY RANDOM ROTATION, USE UNIFORM HAAR MEASURE, 
						// TAKE RANDOM SQUARE MATRIX, WITH INDEPENDENT IDENTICAL GAUSSIAN ENTRIES, 
						// TAKE QR DECOMP AND USE Q MATRIX AS RANDOM ROTATION TO APPLY TO ALL VALS
						////
						Random rand_gen = new Random();
						RealMatrix gauss_mat = new Array2DRowRealMatrix(d,d);
						for(int rr = 0;  rr < gauss_mat.getRowDimension() ; rr++) {
							for(int cc = 0 ; cc < gauss_mat.getColumnDimension(); cc++) {
								gauss_mat.setEntry(rr, cc, rand_gen.nextGaussian());
							}
						}
						
						// obtain random rotation matrix
						QRDecomposition qr = new QRDecomposition(gauss_mat);
						RealMatrix random_rotation = qr.getQ();
						
						// rotate simplex vertices
						start_simplex = start_simplex.transpose();
						start_simplex = random_rotation.multiply(start_simplex);
						start_simplex = start_simplex.transpose();
						
						
						
						///////////////////////////////////////////////////////////////
						// initial guess is the startPoint 
						RealVector ig = new ArrayRealVector(startPoint);
						
						///////////////////////////////////////////////////////////
						
						// add simplex vertices to a list to keep track
						PointValuePair[] points = new PointValuePair[start_simplex.getRowDimension() + 1];
						for(int i = 0 ; i < start_simplex.getRowDimension() ; i++) {
							double[] p;
							
							// add vertex to list
							p = ig.add(start_simplex.getRowVector(i)).toArray();
							points[i] = new PointValuePair(p, qfunc_sbuild.value(p));
							
						}
						// add the origin as a final point to include it in the comparison as well
						double[] p = ig.toArray();
						points[start_simplex.getRowDimension()] = new PointValuePair(p, qfunc_sbuild.value(p));
						
						//////////////////////////////////////////////////////////////
						// setup comparator and arrange LL in order from best (highest) to worst (lowest)
						final Comparator<PointValuePair> pvp_comp  = new Comparator<PointValuePair>() {
						            public int compare(final PointValuePair o1, final PointValuePair o2) {
						                final double v1 = o1.getValue();
						                final double v2 = o2.getValue();
						                return Double.compare(v2, v1);
						            }
						 };
						 // rearrange points in proper order
						 Arrays.sort(points, pvp_comp);
						 
						//////////////////////////////////////////////////////////////
						//////////////////////////////////////////////////////////////

						 // only keep the best points
						 PointValuePair[] simplex = new PointValuePair[d+1];
						 for(int i = 0 ; i < simplex.length ; i++) {
							 simplex[i] = points[i];
						 }
						 
						 this.setPoints(simplex);
						 
			
		}



	}

	private class FlexibleOptimizer {
		
		int optimizer_type; // type of optimizer/simplex
		double optimizer_scale; // the scale at which the optimizer is initialized (related to step size)
		int max_m_evals; // max function evals allowed in optimizing 
		MultivariateFunction q_function; // q-function
		
		
		public FlexibleOptimizer(int o_type, double s_size, int m_evals, MultivariateFunction q) {
			optimizer_type = o_type;
			optimizer_scale = s_size;
			max_m_evals = m_evals;
			q_function = q;
		}
		
		public double[] maximize(double[] i_guess, int em_step) {
		
			// use apache nelder mead
			//double[] out = apache_NM(i_guess);	
			
			// use Mike Hutt Nelder Mead
			double[] out = MH_NM(i_guess);
			
			
			return out;
			
		}
		
		
		private double[] MH_NM(final double[] startPoint) {
			
			// define wrapper function of q-function since MHUTT_NM only minimizes, and we need to maximize
			final class wrap_func implements MultivariateFunction
			{	
				@Override
				public double value(double[] x){
					return -q_function.value(x);
				}
			}		
			
			int nm_iters = 1000; // number of iterations for NM
			double eps = .0000001; // precision for convergence criterion
			
			MHutt_NMSimplex nm_simplex = new MHutt_NMSimplex(nm_iters); // max_iter in NM opt.
			PointValuePair pvp = nm_simplex.minimize(max_m_evals, startPoint, eps, optimizer_scale, new wrap_func(), optimizer_type);
			return pvp.getPoint();
		}
		
		private double[] apache_NM( final double[] startPoint){
			

			SimplexOptimizer opt = new SimplexOptimizer(.0000001, .0000001);
			NelderMeadSimplex n_m_s = new ApacheNelderMeadSimplex_custom(optimizer_type, optimizer_scale, startPoint.length, q_function ); // start_config, rho, khi , gamma, sigma

			// for apache optimizer need:
			ObjectiveFunction llFunction = new ObjectiveFunction( q_function );
			InitialGuess iG = new InitialGuess(startPoint);

			
			PointValuePair pvp;
			try { pvp = opt.optimize(new MaxIter(10000), new MaxEval(max_m_evals) , n_m_s , GoalType.MAXIMIZE, iG ,llFunction);;}
			catch(Exception TooManyEvaluationsException) { 	pvp = n_m_s.getPoint(0); }
			
			return pvp.getPoint();
			
		}
		
	}
	
	
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//
	public static void main(String[] args) {
		
		Random rr = new Random();
		for(int i = 0 ; i < 10 ; i++) {System.out.println(rr.nextGaussian());}
		
	}

	
}

