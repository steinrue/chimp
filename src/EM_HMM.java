import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;
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
	int num_samples;
	
	// HMMs used. 
	protected HMmodel hmm_model;  // current state of the EM process
	private HMmodel calc_hmm;   // a dummy model used for maximization computations
	
	
	// Simplex parameters
	private SimplexOptimizer opt;
	private NelderMeadSimplex n_m_s;
	private RealVector regularization_coefficients = new ArrayRealVector(4);
	
	
	// data Handler
	Helper.dataStreamsHandler dataStreams = null;
	
	//// size of metalocus
	int metaLocus_size; 
	

	
	// Output holders for FB algorithm aggregated over all the chromosomes/subsampling iterations
	boolean fb_output_upToDate;
	RealMatrix[] current_expectation_matrices = null;
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
	
	
	EM_HMM(FunctionUnivariate cr, double rr, double mr, int ns, double[] ts, String genealogy_rep, int mL_size,int[] additional_params, boolean binary_emission, boolean pseudo){
		num_samples = ns;
		metaLocus_size = mL_size;
		
		
		//initialize hmm models with the appropriate parameters
		hmm_model = new HMmodel(cr, rr, mr, ns, ts, genealogy_rep, additional_params, binary_emission, pseudo); 
		calc_hmm = new HMmodel(cr, rr, mr, ns, ts, genealogy_rep, additional_params, binary_emission, pseudo);

		// initialize optimizer with basic accuracy params
		opt = new SimplexOptimizer(.0000001, .0000001); 
		
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
	
	RealVector get_SD() {
		return hmm_model.getSS();
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
		
	}
	
	
	// Perform E step and M step and update. 
	void em_single_update(int m_step_eval_limit, double simplex_span, int simplex_type) {
		
		// redirect to Powell Optimizer right off the bat if simplex type specifies to
		if(simplex_type == 5) {
			em_single_update_Powell( m_step_eval_limit, simplex_span);
			return;
		}
		
		
		double stime = System.nanoTime();
		System.out.println("\n ----------------- Begin EM --------------------");
		
		// E- STEP ////////////////////////////////////////
		///////////////////////
		// perform forward backwards only if we havent already
		if(!fb_output_upToDate) {this.fb_overAllData();}

		// get expectations 
		
		RealMatrix transExp = this.current_expectation_matrices[0];  
		RealMatrix emitExp = this.current_expectation_matrices[1]; ;	
		RealVector startExp = this.current_expectation_matrices[2].getColumnVector(0);
		
		
		// M- STEP ////////////////////////////////////////
		///////////////////////
		
		// Prepare the optimization parameters. The objective function (to be maxed), initial guess, & n_m_s
		ObjectiveFunction llFunction = new ObjectiveFunction(  new Q_Function(transExp, emitExp, startExp, dataStreams.getFSS(), calc_hmm) );
		// next either pick simplex initialization point. 
		InitialGuess iG = new InitialGuess(params_From_Hmm(hmm_model));
		// initialize simplex
		n_m_s = new nelderMeadSimplex_custom(simplex_type, simplex_span, iG.getInitialGuess().length, llFunction.getObjectiveFunction() ); // start_config, rho, khi , gamma, sigma

		// PERFORM OPTIMIZATION (most of params for this are written into opt already in the initialization call
		PointValuePair pvp;
		try { pvp = opt.optimize(new MaxIter(10000), new MaxEval(m_step_eval_limit) , n_m_s , GoalType.MAXIMIZE, iG ,llFunction);;}
		catch(Exception TooManyEvaluationsException) { 	pvp = n_m_s.getPoint(0); }
		
		// UPDATE HMM using the params obtained from optimization ( exponentiate if needed)
		params_To_HMM(pvp.getPoint(), hmm_model); 	// update smc appropriately from the optimized params
		this.fb_clear_output(); // resetting params in hmm_model means fb algorithm outputs are no longer up to date
		
		System.out.println("\n \n**************************\n \n EM step time: " + Double.toString((System.nanoTime() - stime)/1000000000.) +"\n\n ***************************** \n\n");

		System.gc();
	}

	void em_single_update_Powell(int m_step_eval_limit, double simplex_span) {
		
		double stime = System.nanoTime();
		System.out.println("\n ----------------- Begin EM --------------------");
		
		// E- STEP ////////////////////////////////////////
		///////////////////////
		// perform forward backwards only if we havent already
		if(!fb_output_upToDate) {this.fb_overAllData();}

		// get expectations 
		
		RealMatrix transExp = this.current_expectation_matrices[0];  
		RealMatrix emitExp = this.current_expectation_matrices[1]; ;	
		RealVector startExp = this.current_expectation_matrices[2].getColumnVector(0);
		
		
		// M- STEP ////////////////////////////////////////
		///////////////////////
		
		// Prepare the optimization parameters. The objective function (to be maxed), initial guess, & n_m_s
		Q_Function qref = new Q_Function(transExp, emitExp, startExp, dataStreams.getFSS(), calc_hmm);
		ObjectiveFunction llFunction = new ObjectiveFunction(  qref );
		// next either pick simplex initialization point. 
		InitialGuess iG = new InitialGuess(params_From_Hmm(hmm_model));
		// initialize simplex
		//n_m_s = new nelderMeadSimplex_custom(simplex_type, simplex_span, iG.getInitialGuess().length, llFunction.getObjectiveFunction() ); // start_config, rho, khi , gamma, sigma
		
		
		// PERFORM OPTIMIZATION (most of params for this are written into opt already in the initialization call
		PointValuePair pvp; 
		PowellOptimizer_custom optPP = new PowellOptimizer_custom(.0000001, .0000001); optPP.step_size_specify(simplex_span);
		try { pvp = optPP.optimize(new MaxIter(10000), new MaxEval(m_step_eval_limit) , GoalType.MAXIMIZE, iG ,llFunction);;}
		catch(Exception TooManyEvaluationsException) { 	pvp = qref.best_eval; } // need to obtain best value from the record in the Q-function itself
		
		// UPDATE HMM using the params obtained from optimization ( exponentiate if needed)
		params_To_HMM(pvp.getPoint(), hmm_model); 	// update smc appropriately from the optimized params
		this.fb_clear_output(); // resetting params in hmm_model means fb algorithm outputs are no longer up to date
		
		System.out.println("\n \n**************************\n \n EM step time: " + Double.toString((System.nanoTime() - stime)/1000000000.) +"\n\n ***************************** \n\n");

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
		current_expectation_matrices = null;
		current_ll = -1;	
	}

	// performs forward backwards using the Current model_hmm values and the data loaded in currently //
	// does this from RAM or Files depending on streams_in_mem boolean
	// after fb writes the ouput expectation matrices and ll to the appropriate holders
	
	protected void fb_overAllData() {
		
		// make sure data is loaded appropriately
		if( dataStreams == null ) {System.out.println("No Data Loaded for EM_HMM"); System.exit(0);}
				
		int num_streams = dataStreams.num_streams;
		
			
		double totalLL = 0 ;
		RealMatrix transExp = new Array2DRowRealMatrix(hmm_model.getTProbs().getRowDimension(),hmm_model.getTProbs().getColumnDimension());  
		RealMatrix emitExp = new Array2DRowRealMatrix(hmm_model.getEProbs().getRowDimension(), hmm_model.getEProbs().getColumnDimension()); ;	
		RealVector startExp = new ArrayRealVector(hmm_model.num_states);

		
		ExpectationsComputer exp_computer = new Expectation_FBskip(hmm_model.getTProbs(), hmm_model.getEProbs(), hmm_model.getSS());
		
		for(int f = 0 ; f < num_streams; f++) {
			
			int[][] dstream = dataStreams.getStream(f);
			
			exp_computer.loadDataFromStream(dstream, metaLocus_size);
			exp_computer.computeExps();

			// add on the expected number of transitions, emissions, and residences for each file
			transExp = exp_computer.getExpT().add(transExp);  
			emitExp = exp_computer.getExpE().add(emitExp);	
			startExp = exp_computer.getExpS().add(startExp);
			
			totalLL = totalLL + exp_computer.get_LL();
			
		}
		
		current_expectation_matrices = new RealMatrix[] {transExp, emitExp, new Array2DRowRealMatrix(startExp.toArray())};
		current_ll = totalLL;
		
		fb_output_upToDate = true;
	}

	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	class HMmodel {
		private boolean binary_emission;
		private boolean pseudohap;
		private TransitionsComputer tPcomp;
		private EmissionsComputer ePcomp;
		
		private RealMatrix transProbs; // trans prob per metalocus (TP per locus raised to mlsize power)
		private RealMatrix emitProbs;
		private RealVector ssDistro;
		
		public int num_states;
		public int num_emissions;
				
		HMmodel(FunctionUnivariate cr, double rr, double mr, int n, double[] ts, String genealogy_rep, int[] additional_params, boolean bin_em, boolean pseudo){
			
			binary_emission = bin_em;
			pseudohap = pseudo;
			int model_n = pseudohap ? 2 * n : n ;
			
			if(genealogy_rep.equals("TH_1P")) {
				if(additional_params != null) {System.out.print("Invalid additional params"); System.exit(0);}
				tPcomp = new TH_1P_Transitions(cr, rr, model_n, ts);
				ePcomp = new TH_1P_Emissions(cr, mr, model_n, ts, pseudohap);
			}
			
			if(genealogy_rep.equals("TL_1P")) {
				int[] ap = additional_params;

				if(additional_params.length == 2) { tPcomp = new TL_1P_Transitions(cr, rr, model_n, ts, ap[1]); }
				else if(additional_params.length == 3) { tPcomp = new TL_1P_Transitions(cr, rr, model_n, ts, new int[] {ap[1],ap[2]}); }
				else {System.out.print("Invalid additional params"); System.exit(0);}
				
				ePcomp = new TL_1P_Emissions(cr, mr, model_n, ts, ap[0], pseudohap);

			}
			
			

			int mls = metaLocus_size;
			transProbs = tPcomp.transitionProbs().power(mls); // metalocus exponentiation
			
			emitProbs = ePcomp.emissionProbs(binary_emission);
			ssDistro = ePcomp.getMargPDF(); // obtain the initial distro from pdf
			
			
			num_states = ssDistro.getDimension();
			num_emissions = emitProbs.getColumnDimension();
			
		}
		
		
		////////////////////////////////////////////////////////
		// ACCESS METHODS 
		
		// update the split time. pick new split time, compute emission and transition probs and write it in 
		void updateDemography(double[][] trans_params, double[][] emit_params) {
			
			// params is list of params. 
			// index 0 is the params needed to update the coalescence_rates in each
			// index 1 is the rr, and index 2 is the mr
			
			tPcomp.updateDemoFromPs(trans_params);
			ePcomp.updateDemoFromPs(emit_params);
			
			System.out.println("-----------------------");
			System.out.println("-----------------------");

		
			int mls = metaLocus_size;
			transProbs = tPcomp.transitionProbs().power(mls);
			
			emitProbs = ePcomp.emissionProbs(binary_emission);	
			ssDistro = ePcomp.getMargPDF();
			

			
		}
		
		double[][][] getDemoParams(){
			double[][][] out = new double[2][][];
			out[0] = tPcomp.getPsFromDemo();
			out[1] = ePcomp.getPsFromDemo();
			
			return out;
		}
		
		RealMatrix getTProbs() {
			
			return transProbs.copy();
		}
		
		RealMatrix getEProbs() {
			return emitProbs.copy();
		}
		
		RealVector getSS() {
			return ssDistro.copy();
		}
		
		RealVector getRegularization(double estimated_cr) {
			// array of regularizing quantities
			// [d0l2, d1l1, d1l2, d2l2]
			
			RealVector out = new ArrayRealVector(4);
			
			out.setEntry(0, tPcomp.coalescence_rate.reg_d0l2(estimated_cr));
			out.setEntry(1, tPcomp.coalescence_rate.reg_d1l1());
			out.setEntry(2, tPcomp.coalescence_rate.reg_d1l2());
			out.setEntry(3, tPcomp.coalescence_rate.reg_d2l2());

			return out;
		}
		
	}
	
	// maximization function class, constructed using the expectation matrices
	public class Q_Function implements MultivariateFunction{
		
		PointValuePair best_eval = new PointValuePair(new double[] {}, Double.NEGATIVE_INFINITY);		
		
		private RealMatrix transX0;
		private RealMatrix emitX0;
		private RealVector startX0;
		private HMmodel hmC;
				
		private double e_CR; // estimated effective coalescent rate from Wattersons estimator
		
		Q_Function(RealMatrix tX, RealMatrix eX, RealVector sX, double eNe, HMmodel hm){
				transX0 = tX;
				emitX0 = eX;
				startX0 = sX;
				hmC = hm;
								
				e_CR = helper.estimate_watersons_CR(dataStreams.getFSS(), num_samples,  hm.getDemoParams()[1][1][0]);
				
		}
				
		public double value(double[] log_args) {
					System.gc();

			
					params_To_HMM(log_args, hmC);
					
					System.out.println("#PARAMS: " + Arrays.toString(log_args));


					RealMatrix tP = hmC.getTProbs();
					RealMatrix eP = hmC.getEProbs();
					RealVector sP = hmC.getSS();

					///// do a dimension check /////
					if(eP.getColumnDimension() != emitX0.getColumnDimension() || eP.getRowDimension() != emitX0.getRowDimension()) {
						System.out.println("expectation-prob matrix dimension mismatch."); System.exit(0);
					}
					if(tP.getColumnDimension() != transX0.getColumnDimension() || tP.getRowDimension() != transX0.getRowDimension()) {
						System.out.println("expectation-prob matrix dimension mismatch."); System.exit(0);
					}
					
					
					
					double a = 0;
					double b = 0;
					double c = 0;
					
					// add the transition contributions
					for(int i = 0 ; i < tP.getRowDimension(); i++) {
						for(int j = 0 ; j < tP.getColumnDimension() ; j++) {
							if(transX0.getEntry(i, j) > 0. && tP.getEntry(i, j) > 0.) { a = a + Math.log(tP.getEntry(i, j))*(transX0.getEntry(i, j));}
						}
					}
					// add the emission contributions.
					for(int i = 0 ; i < eP.getRowDimension(); i++) {
						for(int j = 0 ; j < eP.getColumnDimension() ; j++) {
							if(emitX0.getEntry(i, j) > 0. && eP.getEntry(i, j) > 0.) { b = b + Math.log(eP.getEntry(i, j))*(emitX0.getEntry(i, j)); }
						}
					}
					// add initial distro contribution
					for(int i = 0 ; i < sP.getDimension() ; i++) {
						if(startX0.getEntry(i) > 0. && sP.getEntry(i) > 0.) { c = c + Math.log(sP.getEntry(i))*(startX0.getEntry(i));}
					}
					double tt = a+b+c;
					
					
					if(tt == 0.0) {
						System.out.println("-- Warning, Q-function forced value --");
						tt = Double.NEGATIVE_INFINITY;
					}
					RealVector reg_qs = hmC.getRegularization(e_CR);
					System.out.println("(Q_eval:"+tt +" ; Reg:"+ reg_qs+")");
					
					
					if(print_log) { try{scribe.append(tt + ", ");
						double[] tarr = helper.ebeExp(log_args);
						for(int ii = 0 ; ii < tarr.length;ii++) {
							scribe.append(tarr[ii] + ", ");
						} scribe.append("\n");
					}catch (IOException e){;}}
					
					
					
					
					double 	qval = tt - regularization_coefficients.dotProduct(reg_qs);
					
					// keep record of the best evaluation for this Q-function
					if(qval > best_eval.getValue()) { best_eval = new PointValuePair(log_args, qval);}
					
					return qval;
		}
					

	}
	
	
	
	private class nelderMeadSimplex_custom extends NelderMeadSimplex{
		
		double simplex_initialization_scale;
		int simplex_type;
		MultivariateFunction qfunc_sbuild;
		
		public nelderMeadSimplex_custom( int s_type, double s_step, int space_dim, MultivariateFunction q){
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

		
	public class PowellOptimizer_custom  extends MultivariateOptimizer {
	    
		// Code was copy-pasted from apache library, but then tweaked to allow for the simplex size to be specified
		// this tweak is in definition of step_size field, and is used in doOptimize() method
		
		
	    //Minimum relative tolerance.
	    private final double MIN_RELATIVE_TOLERANCE = 2 * FastMath.ulp(1d);
	    
	    // Relative threshold.
	    private final double relativeThreshold;
	    
	    // Absolute threshold.
	    private final double absoluteThreshold;
	    
	    // Line search.
	    private final LineSearch line;
	
	    // initial step size = 1.
	    private double step_size = Double.NaN;
	    
	    
	    public PowellOptimizer_custom(double rel, double abs, ConvergenceChecker<PointValuePair> checker) {
	        this(rel, abs, FastMath.sqrt(rel), FastMath.sqrt(abs), checker);
	    }
	
	    public PowellOptimizer_custom(double rel, double abs, double lineRel, double lineAbs,ConvergenceChecker<PointValuePair> checker) {
	        super(checker);
	
	        if (rel < MIN_RELATIVE_TOLERANCE) {
	            throw new NumberIsTooSmallException(rel, MIN_RELATIVE_TOLERANCE, true);
	        }
	        if (abs <= 0) {
	            throw new NotStrictlyPositiveException(abs);
	        }
	        relativeThreshold = rel;
	        absoluteThreshold = abs;
	
	        // Create the line search optimizer.
	        line = new LineSearch(this, lineRel, lineAbs, 1d);
	       
	    }
	
	    public PowellOptimizer_custom(double rel, double abs) {
	        this(rel, abs, null);
	    }
	
	    public PowellOptimizer_custom(double rel, double abs, double lineRel, double lineAbs) {
	        this(rel, abs, lineRel, lineAbs, null);
	    }
	
	    /** {@inheritDoc} */
	    @Override
	    protected PointValuePair doOptimize() {
	        //if(true) {return doOptimizeOriginal();}
	    	
	    	checkParameters();
	
	        final GoalType goal = getGoalType();
	        final double[] guess = getStartPoint();
	        final int n = guess.length;
	
	        final double[][] direc = new double[n][n];
	        for (int i = 0; i < n; i++) {
	            direc[i][i] = step_size;
	        }
	
	        final ConvergenceChecker<PointValuePair> checker
	            = getConvergenceChecker();
	
	        double[] x = guess;
	        double fVal = computeObjectiveValue(x);
	        double[] x1 = x.clone();
	        while (true) {
	            incrementIterationCount();
	
	            double fX = fVal;
	            double fX2 = 0;
	            double delta = 0;
	            int bigInd = 0;
	            double alphaMin = 0;
	            
	            for (int i = 0; i < n; i++) {
	            	System.out.println("\n \n \n Unidirectional search in D"+i+" direction \n ");
	            	
	                final double[] d = MathArrays.copyOf(direc[i]);
	
	                fX2 = fVal;
	
	                final UnivariatePointValuePair optimum = line.search(x, d);
	                fVal = optimum.getValue();
	                alphaMin = optimum.getPoint();
	                final double[][] result = newPointAndDirection(x, d, alphaMin);
	                x = result[0];
	
	                if (FastMath.abs(fX2 - fVal) > delta) {
	                    delta = fX2 - fVal;
	                    bigInd = i;
	                }
	            }
	
	            // Default convergence check.
	            boolean stop = 2 * FastMath.abs(fX - fVal) <= 
	                (relativeThreshold * (FastMath.abs(fX) + FastMath.abs(fVal)) +
	                 absoluteThreshold);   
	
	            final PointValuePair previous = new PointValuePair(x1, fX);
	            final PointValuePair current = new PointValuePair(x, fVal);
	            if (!stop && checker != null) { // User-defined stopping criteria.
	                stop = checker.converged(getIterations(), previous, current);
	            }
	            if (stop) {
	                if (goal == GoalType.MINIMIZE) {
	                    return (fVal < fX) ? current : previous;
	                } else {
	                    return (fVal > fX) ? current : previous;
	                }
	            } 
	
	            final double[] d = new double[n];
	            final double[] x2 = new double[n];
	            for (int i = 0; i < n; i++) {
	                d[i] = x[i] - x1[i];
	                x2[i] = 2 * x[i] - x1[i];
	            }
	
	            x1 = x.clone();
	            fX2 = computeObjectiveValue(x2);
	
	            
	            ///////// now optimize in aggregate direction, and replace biggest biggest contributor with aggregate.

		        System.out.println("\n \n \n Unidirectional search in D* direction \n ");
	            final UnivariatePointValuePair optimum = line.search(x, d);
	            fVal = optimum.getValue();
	            alphaMin = optimum.getPoint();
	            final double[][] result = newPointAndDirection(x, d, alphaMin);
	            x = result[0];
	
	            final int lastInd = n - 1; // 
	            direc[bigInd] = direc[lastInd];
	            direc[lastInd] = result[1];
	
	        }
	    }
	
	    protected PointValuePair doOptimizeOriginal() {
	        // this is the original implementation with the minor modification to adjust initial step sizes. 
	    	
	    	checkParameters();
	
	        final GoalType goal = getGoalType();
	        final double[] guess = getStartPoint();
	        final int n = guess.length;
	
	        final double[][] direc = new double[n][n];
	        for (int i = 0; i < n; i++) {
	            direc[i][i] = step_size;
	        }
	
	        final ConvergenceChecker<PointValuePair> checker
	            = getConvergenceChecker();
	
	        double[] x = guess;
	        double fVal = computeObjectiveValue(x);
	        double[] x1 = x.clone();
	        while (true) {
	            incrementIterationCount();
	
	            double fX = fVal;
	            double fX2 = 0;
	            double delta = 0;
	            int bigInd = 0;
	            double alphaMin = 0;
	            
	            for (int i = 0; i < n; i++) {
	            	System.out.println("\n \n \n Unidirectional search in D"+i+" direction \n ");
	            	
	                final double[] d = MathArrays.copyOf(direc[i]);
	
	                fX2 = fVal;
	
	                final UnivariatePointValuePair optimum = line.search(x, d);
	                fVal = optimum.getValue();
	                alphaMin = optimum.getPoint();
	                final double[][] result = newPointAndDirection(x, d, alphaMin);
	                x = result[0];
	
	                if ((fX2 - fVal) > delta) {/////////*************
	                    delta = fX2 - fVal;
	                    bigInd = i;
	                }/////////*************
	            }
	
	            // Default convergence check.
	            boolean stop = 2 * (fX - fVal) <=      /////////*************
	                (relativeThreshold * (FastMath.abs(fX) + FastMath.abs(fVal)) +
	                 absoluteThreshold);   /////////*************
	
	            final PointValuePair previous = new PointValuePair(x1, fX);
	            final PointValuePair current = new PointValuePair(x, fVal);
	            if (!stop && checker != null) { // User-defined stopping criteria.
	                stop = checker.converged(getIterations(), previous, current);
	            }
	            if (stop) {/////////*************
	                if (goal == GoalType.MINIMIZE) {
	                    return (fVal < fX) ? current : previous;
	                } else {
	                    return (fVal > fX) ? current : previous;
	                }
	            } /////////*************
	
	            final double[] d = new double[n];
	            final double[] x2 = new double[n];
	            for (int i = 0; i < n; i++) {
	                d[i] = x[i] - x1[i];
	                x2[i] = 2 * x[i] - x1[i];
	            }
	
	            x1 = x.clone();
	            fX2 = computeObjectiveValue(x2);
	
	            if (fX > fX2) {/////////*************
	                double t = 2 * (fX + fX2 - 2 * fVal);
	                double temp = fX - fVal - delta;
	                t *= temp * temp;
	                temp = fX - fX2;
	                t -= delta * temp * temp;
	        /////////*************
	                if (t < 0.0) {
		            	System.out.println("Unidirectional search in D"+n+" direction");
	                    final UnivariatePointValuePair optimum = line.search(x, d);
	                    fVal = optimum.getValue();
	                    alphaMin = optimum.getPoint();
	                    final double[][] result = newPointAndDirection(x, d, alphaMin);
	                    x = result[0];
	
	                    final int lastInd = n - 1;
	                    direc[bigInd] = direc[lastInd];
	                    direc[lastInd] = result[1];
	                }
	            }
	        }
	    }
	
	    
	    protected void step_size_specify(double new_ss) {
	    	step_size = new_ss;
	    }
	    
	    
	    private double[][] newPointAndDirection(double[] p,
	                                            double[] d,
	                                            double optimum) {
	        final int n = p.length;
	        final double[] nP = new double[n];
	        final double[] nD = new double[n];
	        for (int i = 0; i < n; i++) {
	            nD[i] = d[i] * optimum;
	            nP[i] = p[i] + nD[i];
	        }
	
	        final double[][] result = new double[2][];
	        result[0] = nP;
	        result[1] = nD;
	
	        return result;
	    }
	
	    
	    private void checkParameters() {
	        if (getLowerBound() != null ||
	            getUpperBound() != null) {
	            throw new MathUnsupportedOperationException(LocalizedFormats.CONSTRAINT);
	        }
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

