import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.martiansoftware.jsap.JSAPException;


public class Expectation_FBskip extends ExpectationsComputer {
	Helper helper = new Helper();

	
	// data loaded in is stored here
	private int[][] obs;
	private int mL_size; // size of metalocus (in bases)
	private int n_loci;
	private int chrom_length;
	
	// combined emission/observation information for use in mL implementation
	private RealVector[] mL_emissions;
	
	// matrix decomposition info for monomorphic seg site type
	private SingleSiteMovementMatrix mono_nss;
	private SingleSiteMovementMatrix mono_mis;

	// choose computer needed for metaLocus emissions
		
	
//////////////////////////////////////////////////////////////////////////////

	
	
	public static void main(String[] args) {
	}

		
	
//////////////////////////////////////////////////////////////////////////////

	// CONSTRUCTOR AND PRIVATE CLASS FOR DECOMPOSITION INFO ////////////

	Expectation_FBskip(RealMatrix tp, RealMatrix ep, RealVector sd) {
		super(tp, ep, sd);
		

	}
	
	private class SingleSiteMovementMatrix{
			
			public RealMatrix m_b;			// The B matrix (diagonal from emission column). W = T.B
			public RealMatrix m_mat;    	// the matrix to be exponentiated, := W
			public RealMatrix m_p;		 	// diagonalizING matrix for w. W = P * Diag[e_vals] * P^-1
			public RealMatrix m_pi;			// p^-1 or p_inverse
			public RealMatrix m_diag;		// diagonal matrix of eigenvalues
			public RealVector m_d_log_rel; 	// vector of log of eigenvalues in diag, relative to some sigma_bar
			public double m_d_log_scaling; 	// median of log(vector of eigenvalues of w)
			
			public RealMatrix m_mat_inv;    // inverse of main matrix W
			
			public RealMatrix m_eigen_diff; 	// matrix of pairwise eigenvalue differences.
			public RealMatrix m_log_eigenrow;	// each entry has the log of the eigenvalue (relative) corresponding to the row num
			
			public int type; // type of monomorphic site. 1 means non segregating
			
			SingleSiteMovementMatrix(RealVector b, int t){
				type = t;
				m_b = new DiagonalMatrix(b.toArray());
				m_mat = tM0.multiply(m_b);
				m_mat_inv = MatrixUtils.inverse(m_mat);
				
				EigenDecomposition eigen = new EigenDecomposition( m_mat );

				// save the eigendecomposition into usable form
				m_p = eigen.getV();
				m_pi = MatrixUtils.inverse(m_p);
				m_diag = eigen.getD();
				
				// set up the log rep of the diagonal matrix, using the median log value
				// and the deviations from this to represent each eigenvalue
				m_d_log_rel = new ArrayRealVector(m_diag.getRowDimension());
				for(int i = 0 ; i < m_d_log_rel.getDimension(); i++) {m_d_log_rel.setEntry(i, Math.log(m_diag.getEntry(i, i)));}
				m_d_log_scaling = ( m_d_log_rel.getMaxValue() + m_d_log_rel.getMinValue() ) / 2.;
				m_d_log_rel.mapSubtractToSelf(m_d_log_scaling);

				build_eigen_diff();
			}
			
			void build_eigen_diff(){
				// first create m_log_eigenrow
				m_log_eigenrow = new Array2DRowRealMatrix(m_d_log_rel.getDimension(), m_d_log_rel.getDimension());
				for(int r = 0 ; r < m_d_log_rel.getDimension() ; r++) {
					m_log_eigenrow.setRow(r, m_d_log_rel.toArray());
				}
				
				// m_eigen_diff. set the diagonal entries to 1, since we divide by this matrix elementwise later, and dont want to have 0 vals
				m_eigen_diff = diag_diff_operation(1);
				//for(int i = 0 ; i < m_eigen_diff.getRowDimension(); i++) {m_eigen_diff.setEntry(i, i, 1);}
				
			}
		
			RealMatrix diag_diff_operation(int pow) {
				// this will return the unity ordered matrix X_ab = D_aa ^ pow - D_bb ^ pow
				
				RealMatrix temp;
				
				temp = ebe_exp(m_log_eigenrow.scalarMultiply(pow));
				temp = temp.subtract(temp.transpose());
				
				return temp;
				
			}
		
			RealMatrix ebe_exp(RealMatrix m) {
				
				RealMatrix temp = new Array2DRowRealMatrix(m.getRowDimension(), m.getColumnDimension());
				for(int i = 0 ; i < m.getRowDimension(); i++) {
					for(int j = 0 ; j < m.getColumnDimension(); j++) {
						temp.setEntry(i, j, Math.exp(m.getEntry(i, j)));
					}
				}
				
				return temp;	
			}
		
			RealMatrix obtain_Q_unity(int tract_length) {
			
				// compute all entries of q, this is not correct for diagonal entries
				RealMatrix q = this.diag_diff_operation(tract_length - 1);
				q = helper.ebeDivide(q, this.m_eigen_diff);
				
				//compute diagonal entries
				RealVector temp = m_d_log_rel.mapMultiply(tract_length -2);
				for(int i = 0 ; i < temp.getDimension(); i++) {
					q.setEntry(i, i, Math.exp(temp.getEntry(i)) * (tract_length - 1) );
				}
			
				// now q is the true matrix scaled by exp(sigma_bar * (tract_length - 2) )
				return q;
			}
		
		}
		
	
	
	
	//////////////////////////////////
	// ABSTRACT METHOD NEED TO IMPLEMENT///
	//////////////////////////////////
	
	//@Override
	void loadDataFromStream(int[][] observed_data, int mLs) {
		obs = observed_data;
		n_loci = obs[0].length; // number of tracts or mLs depending on mLsize
		mL_size = mLs;
		mL_emissions = null;
		
		// get the actual length of observation in sites
		
		// do this if we are not using mL version, and doing locus-skipping implementation
		if(mL_size == 1) { // mL_size=1 defaults to locus skipping
			chrom_length = 0 ;  
			for(int i = 0 ; i < n_loci ; i++) { chrom_length = chrom_length + obs[1][i]; }
		}
		
	
		// relevant things if we are using mL
		else if(mL_size > 1) {
			chrom_length = n_loci * mL_size;
			mL_emissions = new RealVector[n_loci];
			for(int i = 0 ; i < n_loci ; i++) {mL_emissions[i] = mL_emitP(i);}
		}
		
		else {System.err.print("Metalocus Size Invalid"); System.exit(0);}
		
	}
	
	// compute the expectations by performing forwards backwards algorithm
	@Override
	void computeExps() {
		// reroute to new method if need to use metalocus framework
		if(mL_size > 1) {computeExps_mL(); return;} // mL_size=1 defaults to locus skipping
		
		// get matrix parts for monomorphic tracts of non segregating sites. only do this if doing locus skipping method, otherwise don't need ot
		mono_nss = new SingleSiteMovementMatrix(eM0.getColumnVector(0), 1);
		mono_mis = new SingleSiteMovementMatrix(new ArrayRealVector( tM0.getColumnDimension(), 1.0) ,2);
		///
		
		
		double stime = System.nanoTime(); // just to time this method
		
		initializeExp();
		
		RealVector[] salphas = new RealVector[n_loci];
		RealVector[] salphas_left_adjacent = new RealVector[n_loci];
		double[] scale_cs = new double[n_loci];
		double[] norm_unity_factors = new double[n_loci];
		
		// FORWARDS ALGORITHM
		forwardsFill(salphas, salphas_left_adjacent, scale_cs, norm_unity_factors);
		
		// HERE WE PROCEED WITH BACKWARDS ALGORITHM 
		// compute sbeta at the last locus
		
		RealVector sbeta = new ArrayRealVector(sM0.getDimension(), 1.);
		
		// Contribution to gamma and xi and expectation matrices
		buildExpectation_SingleSite(salphas[n_loci-1], sbeta, salphas_left_adjacent[n_loci-1], scale_cs[n_loci-1], obs[0][n_loci-1]);
		
		
		for(int l = n_loci - 2 ; l > -1 ; l--) {
			
			// compute beta at last spot in the tract (from the leftmost point of following tract)
			sbeta = tM0.operate(emissionColumn(obs[0][l+1]).ebeMultiply(sbeta));
			sbeta.mapDivideToSelf(scale_cs[l+1]);
			
			// compute beta at first position in the tract (left most ss)
			int tract_length = obs[1][l];
			
			SingleSiteMovementMatrix mono = null;
			int type = obs[2][l];
			if(type == 1) {mono = mono_nss;}
			else if(type == 2) {mono = mono_mis;}
			else {System.out.println("Invalid SS Tract Type"); System.exit(0);}
			
			
			if(tract_length == 1) {}// nothing to do
			else {
				
				// Contribution to gamma and xi and exp for monomorphic sites if skipping
				// ****************************
				buildExpectation_MonoTract(salphas[l],sbeta, mono, tract_length);
				
				double[] exp_sig = mono.m_d_log_rel.mapMultiply(tract_length - 1).toArray();
				for(int i = 0; i < exp_sig.length ; i++) {exp_sig[i] = Math.exp(exp_sig[i]);}
				
				// now set beta =  P * D^# * P^-1 * beta
				sbeta = mono.m_p.multiply(new DiagonalMatrix(exp_sig)).operate(mono.m_pi.operate(sbeta));
				
				// renormalize by appropriate factor from forwards alg
				sbeta.mapDivideToSelf(norm_unity_factors[l+1]);
				
			}
			
			// now sbeta should sync with leftmost (ss) spot in tract
			// Contribution to gamma and xi and exp for the leftmost (ss) spot in the tract
			buildExpectation_SingleSite(salphas[l], sbeta, salphas_left_adjacent[l], scale_cs[l], obs[0][l]);
			
		}
		
		// find gamma @ locus 0 to compute starting distribution
		buildExpectation_SD(salphas[0], sbeta);
		
		
		System.out.println("\n <--> FB algorithm: " + Double.toString((System.nanoTime() - stime)/1000000000.) + " seconds\n");
		System.out.println("LL="+posteriorLL+"\n");

		return;
	}

///////////////////////////////////////////////////////////////////////////
	
	// These are FB algorithm implemented by skipping along monomorphic tracts
	
	private void forwardsFill(RealVector[] salphas, RealVector[] salphas_l, double[] scale_cs, double[] norm_unity_factors) {
		
		// salphas is holder of RealVectors at each segregating site locus (start of each tract)
		// salphas_l is holder of alpha RealVectors at the sites left adjacent to the tract starts.
		// scale_cs is holder for the c values for the ss loci 
		// norm_unity_factors are contributions to product of cs for monomorphic tract to the left that
		// does not include the exp( w_d_log_scaling* [tract_length-1] ) contribution.
		
		//FIND INITIAL VALUES AT START OF GENOME SEQUENCE
		
		// total posterior ll counter (cumulative c, or log(product of all c's)
		double cumulative_ll_c = 0;
		
		// first compute scale_c(0) 
		RealVector temp0 = emissionColumn(obs[0][0]).ebeMultiply(sM0);
		
		double c = temp0.getL1Norm();
		temp0.mapDivideToSelf(c); // temp0 is rescaled

		// write into initial positions of salphas and c, and update cumulative log likelihood
		salphas[0] = temp0;
		salphas_l[0] = new ArrayRealVector(temp0.getDimension());
		scale_cs[0] = c;
		cumulative_ll_c = cumulative_ll_c + Math.log(c);
		
		
		for(int l = 1 ; l < n_loci;l++) {
			
			int tract_length = obs[1][l-1];
			
			SingleSiteMovementMatrix mono = null;
			int type = obs[2][l-1];
			if(type == 1) {mono = mono_nss;}
			else if(type == 2) {mono = mono_mis;}
			else {System.out.println("Invalid SS Tract Type"); System.exit(0);}
						
			if(tract_length == 1) {}  //  temp0 = salpha at position before next seg site
			else {
				

				// compute the order 1 scaled diagonal matrix raised to the appropriate power to get to l-1 locus
				double[] exp_sig = mono.m_d_log_rel.mapMultiply(tract_length - 1).toArray();
				for(int i = 0; i < exp_sig.length ; i++) {exp_sig[i] = Math.exp(exp_sig[i]);}
				
				// now set temp0 = temp0 * P * D^# * P^-1  (where D is renormalized diagonal of order 1)
				// this neglects overall factor of exp(# * w_d_log_scaling)
				temp0 = mono.m_p.multiply(new DiagonalMatrix(exp_sig)).multiply(mono.m_pi).preMultiply(temp0);
				
				// renormalize (not strictly necessary, but doesnt hurt to keep things order 1)
				c = temp0.getL1Norm();
				temp0.mapDivideToSelf(c);
				
				// add contribution (from previous_locus + 1 till l-1) of log(product of c's)
				cumulative_ll_c = cumulative_ll_c + Math.log(c) + mono.m_d_log_scaling * (tract_length - 1) ;
				norm_unity_factors[l] = c;
				
			}
			
			
			
			// write in salpha for previous site (left of current)
			salphas_l[l] = temp0; 
			
			// temp0 and c are currently set for locus l-1. now extend to l. 
			temp0 = emissionColumn(obs[0][l]).ebeMultiply(tM0.preMultiply(temp0));
			c = temp0.getL1Norm();
			temp0.mapDivideToSelf(c);
			
			// write to holder for locus l, and add contribution to log likelihood counter
			salphas[l] = temp0;
			scale_cs[l] = c;
			cumulative_ll_c = cumulative_ll_c + Math.log(c);
			
		}
		
		posteriorLL = cumulative_ll_c;
		
	}


	private void buildExpectation_SingleSite(RealVector salpha, RealVector sbeta, RealVector salpha_l, double cc, int nss_obs) {
		// input at locus salpha, sbeta, alpha for site to the left, scale_c value, and observed number of ss

		RealVector gamma = salpha.ebeMultiply(sbeta);
		
		RealMatrix xi = new Array2DRowRealMatrix(sM0.getDimension(), sM0.getDimension());
		RealVector temp = new ArrayRealVector(sM0.getDimension());
		for(int r = 0 ; r < sM0.getDimension();r++) {
			temp = tM0.getRowVector(r).ebeMultiply(emissionColumn(nss_obs).ebeMultiply(sbeta));
			temp.mapMultiplyToSelf(salpha_l.getEntry(r)/cc);
			xi.setRowVector(r, temp);
		}
		
		// now we have gamma and xi, contribute to the expectation holders
		expTr = expTr.add(xi);
		if(nss_obs != -1) {expEm.setColumnVector(nss_obs, expEm.getColumnVector(nss_obs).add(gamma));}
		
	}
	

	private void buildExpectation_MonoTract(RealVector salpha_lb, RealVector sbeta_rb, SingleSiteMovementMatrix w, int tract_length) {
		
		// compute the u matrix (to be hadamard multiplied by the Q matrix)
		RealMatrix u = w.m_pi.multiply(sbeta_rb.outerProduct(salpha_lb)).multiply(w.m_p);
		
		// compute q matrix (unity order scaled)
		RealMatrix q = w.obtain_Q_unity(tract_length);
		
		RealMatrix temp0 = helper.ebeMultiply(u,q); 
		temp0 = w.m_p.multiply(temp0); // now this is P.(U had Q)
		
		/////////////////////////////////
		// COMPUTE CONTRIBUTIONS TO EMISSION EXPECTATIONS
		////////////////////////////////
		
		RealMatrix temp1 = temp0.multiply(w.m_diag);   // this is P.(U had Q).D
		
		// now multiply by P^-1 but only retain diagonal entries in sum_gamma
		RealVector gamma_sum =  new ArrayRealVector(sM0.getDimension());
		for(int i = 0 ; i < gamma_sum.getDimension(); i++) {
			gamma_sum.setEntry(i, temp1.getRowVector(i).dotProduct(w.m_pi.getColumnVector(i)));
		}
		gamma_sum.mapMultiplyToSelf((tract_length - 1) / gamma_sum.getL1Norm() );		// rescale so l1 norm is number of monomorphic summant sites
		
		// contribute to expEm if its nss tract only, contribute 0 ss vals
		if(w.type==1) {	expEm.setColumnVector(0, gamma_sum.add(expEm.getColumnVector(0))); }
		
		/////////////////////////////////
		// COMPUTE CONTRIBUTIONS TO TRANSITION EXPECTATIONS
		////////////////////////////////
		
		RealMatrix temp2 = temp0.multiply(w.m_pi);
		temp2 = w.m_b.multiply(temp2).transpose();
		RealMatrix xi_sum = helper.ebeMultiply(temp2, tM0);
		xi_sum = xi_sum.scalarMultiply((tract_length - 1)/helper.totalEntriesMat(xi_sum));
		
		// contribute to expTr regardless of tract type, all tracts will have state transitions
		expTr = expTr.add(xi_sum) ;

	}
	
	private void buildExpectation_SD( RealVector salpha, RealVector sbeta) {
		
		RealVector gamma = salpha.ebeMultiply(sbeta);
		gamma.mapMultiplyToSelf(gamma.getL1Norm());
		
		expSd = gamma;
		
	}

	private RealVector emissionColumn(int i) {
		if(i == -1) {return new ArrayRealVector(eM0.getRowDimension(),1.0);}		
		else { return eM0.getColumnVector(i);	}
	}
	
///////////////////////////////////////////////////////////////////////////
	
	// this is implemented to go site by site (no skipping)
	
	// METHODS INVOLVED IN WRITING POSTERIOR DISTROS DOWN
	public void computePosteriorDistros(String output_file, double[] head_bins,  int avg_window_size) throws IOException {
		
		// reroute to new method if need to use metalocus framework
		if(mL_size > 1) {computePosteriorDistros_ML(output_file, head_bins); return;} // mL_size=1 defaults to locus skipping
				
		
		double stime = System.nanoTime(); // just to time this method
				
		RealVector[] salphas = new RealVector[n_loci];
		RealVector[] salphas_left_adjacent = new RealVector[n_loci];
		double[] scale_cs = new double[n_loci];
		double[] norm_unity_factors = new double[n_loci];
		
		// FORWARDS ALGORITHM
		forwardsFill(salphas, salphas_left_adjacent, scale_cs, norm_unity_factors);
		
		
		// setup the output files, and metaloci recorders for the full and reduced printers
		FileWriter full_scribe= new FileWriter(output_file+ ".csv"); full_scribe.append(getHeader(head_bins)+" \n");
		//FileWriter red_scribe= new FileWriter(output_file+ "_reduced.csv"); red_scribe.append(getHeader(head_bins)+" \n");
		metalocusRecorder full_rec = new metalocusRecorder(full_scribe, avg_window_size);
		//metalocusRecorder red_rec = new metalocusRecorder(red_scribe, 1);
		
		
		
		// HERE WE PROCEED WITH BACKWARDS ALGORITHM 
		
		// sbeta and salpha at the last locus
		RealVector sbeta = new ArrayRealVector(sM0.getDimension(), 1.);
		RealVector salpha = salphas[n_loci - 1]; 
		int position = chrom_length;
		
		// record posteriors at the last site
		RealVector gamma = salpha.ebeMultiply(sbeta); gamma.mapDivideToSelf(gamma.getL1Norm()) ;
		full_rec.push(position, gamma);
		//red_rec.push(position, gamma);

		
		for(int tract = n_loci - 2 ; tract > -1 ; tract-- ) { // iterate through each tract in reverse order
			
			// compute beta then alpha at last spot in the tract (from the leftmost point of following tract)
			position--; if(((chrom_length - position) % 100000) == 0) {System.out.println("Computed posterior for "+(chrom_length - position)+" bases");}
			sbeta = tM0.operate(emissionColumn(obs[0][tract+1]).ebeMultiply(sbeta));
			sbeta.mapDivideToSelf(sbeta.getL1Norm());
			salpha = salphas_left_adjacent[tract+1];
			
			// record to full
			gamma = salpha.ebeMultiply(sbeta); gamma.mapDivideToSelf(gamma.getL1Norm()) ;
			full_rec.push(position, gamma);

			int tract_length = obs[1][tract]; // length of current tract
			int tract_type = obs[2][tract];
			
			SingleSiteMovementMatrix mono = null;
			if(tract_type == 1) {mono = mono_nss;}
			else if(tract_type == 2) {mono = mono_mis;}
			else {System.out.println("Invalid SS Tract Type"); System.exit(0);}
			
			
			// walk backwards along all the monomorphic steps
			for(int l = tract_length - 1 ; l > 0 ; l-- ) {
				
				position--; if(((chrom_length - position) % 100000) == 0) {System.out.println("Computed posterior for "+(chrom_length - position)+" bases");}
				
				// compute alpha and beta at new position. note for alpha we are reversing the forwards process
				sbeta = mono.m_mat.operate(sbeta);
				salpha = mono.m_mat_inv.preMultiply(salpha);
				
				// rescale to keep things in line
				sbeta.mapDivideToSelf(sbeta.getL1Norm());
				salpha.mapDivideToSelf(salpha.getL1Norm());
				
				// record to full
				gamma = salpha.ebeMultiply(sbeta); gamma.mapDivideToSelf(gamma.getL1Norm()) ;
				full_rec.push(position, gamma);
				
			}
			
			
			// now salpha and sbeta are the leftmost position of the tract, now want to print to reduced file
			//red_rec.push(position, gamma);

			
			
		}
		
		// make sure the recorders print whatever is left on the stack
		full_rec.print_current_stack();
		//red_rec.print_current_stack();
		
		
		// close files, print time, and terminate method
		full_scribe.flush(); full_scribe.close();
		//red_scribe.flush(); red_scribe.close();
		System.out.println("\n <--> FB algorithm for posteriors: " + Double.toString((System.nanoTime() - stime)/1000000000.) + " seconds\n");
		System.out.println("LL="+posteriorLL+"\n");
		return;
	}
		
	private String getHeader(double[] h_bins) {
		String header = "";
		header = header + "-1,0";//position ";
		/*
		for(int i = 0 ; i < this.sM0.getDimension(); i++) {
			header = header + ", " + "bin " + Integer.toString(i+1) ;
		}
		*/
		
		for(int i = 0 ; i < h_bins.length ; i++) {
			int val = (int) h_bins[i];
			header = header + ",  " + val;
		}
		

		return header;
	}
	
	// this is not tied to mL implementation of data stream
	private class metalocusRecorder{
		FileWriter printer; // file to which we print the locus and distribution
		int max_size; 		// prespecified size of the metaloci
		
		// variables that are changed as we push more values here
		int current_size;	// the number of sites that have been input so far, prints once it reaches max
		RealVector current_total_posterior; // sum of posterior for all sites that have been pushed onto this object so far
		
		int left_position=-1; // stores the position of the leftmost site in the current stack of sites
		int right_position=-1; // stores the position of the rightmost site in the current stack of sites
		
		
		
		metalocusRecorder(FileWriter fw, int loci_size){
			printer = fw;
			max_size = loci_size;
			
			current_size = 0;
			current_total_posterior = new ArrayRealVector(sM0.getDimension());
		}
		
		void push(int pos, RealVector local_posterior) throws IOException{
			if(current_size == 0) {left_position = pos; right_position = pos;}
			else {  
				if(pos < left_position) {left_position = pos;}
				if(pos > right_position) {right_position = pos;}
			}
			
			current_total_posterior = current_total_posterior.add(local_posterior);
			current_size++;
			
			if(current_size == max_size) {
				print_current_stack();
			}
			
			
			
		}
		
		void print_current_stack() throws IOException {
			if(current_size == 0) {return;}
		
			int avg_position = (right_position - left_position + 1) / 2 + left_position;
			
			
			RealVector avg_posterior = current_total_posterior.mapDivide(current_size);
			current_total_posterior = new ArrayRealVector(sM0.getDimension());
			
			// reset variables
			current_size = 0;
			current_total_posterior = new ArrayRealVector(sM0.getDimension());
			
			// print as a line in the file, using the avg position, and the avg posterior over stack
			printer.append(avg_position+"");
			for(int i = 0 ; i < avg_posterior.getDimension(); i++) {
				printer.append(", " + avg_posterior.getEntry(i));
			}
			printer.append(" \n ");
			
			
		}
		
		
	}

///////////////////////////////////////////////////////////////////////////

	

////////////////////////////////////////////////////////////
// METALOCUS IMPLEMENTATION - This is the default that gets used
////////////////////////////////////////////////////////////

	// 
	
	// HAVE ALREADY CONFIRMED THAT ml_size=1 MATCHES RESULTS OF LOCUS SKIPPING, so now ml_size=1 defaults to locus skipping
	void computeExps_mL() {

		double stime = System.nanoTime(); // just to time this method

		initializeExp();
		
		// using obs, write in expectations for expTr and expEm (defined above)
		// note this is written so that we do not change the references for input matrices. we want to write to them 
		// in place, keeping their current references.
		
		// a temporary holder for trans expectations
		RealMatrix teX = expTr.copy();
				
		
		// initialize holders and fill salpha and sbeta vals
		RealVector[] salpha = new RealVector[n_loci];
		RealVector[] sbeta = new RealVector[n_loci];
		double[] scale_lnc = new double[n_loci];
				
		fbAB_sf_mL_fill(salpha, sbeta, scale_lnc);
				
				
		// compute and store the log_likelihoods			
		double likelihood = 0;
		for(int t = 0 ; t < n_loci ; t++) {
			likelihood += scale_lnc[t];
		}
		posteriorLL = likelihood;
		

		
		// expTr should store expected number of transitions from hidden state i to j , given the observed data
		// expEm should store expected number of emissions from hidden state i to observed d , given the observed data
		
		
		// compute gammas and add emission expectations
		for(int t = 0 ; t < n_loci ; t++) {
			RealVector gamma = salpha[t].ebeMultiply(sbeta[t]);
			
			// expected emissions. several emissions should be observed from given metalocus
			for(int j = 0 ; j < obs.length ; j++) {
				expEm.setColumnVector(j,expEm.getColumnVector(j).add( gamma.mapMultiply(obs[j][t]))    );				
			}
		}
		
		// compute xis and add transition expectations
		for(int t = 0 ; t < n_loci - 1 ; t++) {
					
			RealMatrix xi = new Array2DRowRealMatrix(sM0.getDimension(), sM0.getDimension());
					
			// log-probs of states emitting observed sequence at next mL, regularized using scale_c
			RealVector temp = mL_emissions[t+1].copy();
			temp.mapSubtractToSelf(scale_lnc[t+1]);
			temp = helper.ebeExp(temp);
			temp = temp.ebeMultiply(sbeta[t+1]);
					
			// compute xi(t) row by row. 
			for(int r = 0 ; r < xi.getRowDimension() ; r++) {
				RealVector temp2 = tM0.getRowVector(r).mapMultiply(salpha[t].getEntry(r)) ;
				xi.setRowVector(r, temp.ebeMultiply(temp2));
			}
					
			// add on expected transitions
			// for t = n_loci-1 this matrix is defined as zeros, so no harm
			teX = teX.add(xi); 
		}

		// in place fill expSd and expTr
		expSd.setSubVector(0, salpha[0].ebeMultiply(sbeta[0])); 
		expTr.setSubMatrix(teX.getSubMatrix(0, teX.getRowDimension()-1, 0, teX.getColumnDimension()-1).getData(),0, 0);

		
		System.out.print("\n <--> FB algorithm: " + Double.toString((System.nanoTime() - stime)/1000000000.) + " seconds\n");
		System.out.println("LL="+posteriorLL+"\n");
	}
	
	// forwards backwards scaled quantities salphas, sbetas and scale_cs
	
	private void fbAB_sf_mL_fill(RealVector[] salpha, RealVector[] sbeta, double[] scale_lnc ){
				
		// FORWARDS ALGORITHM
		
		RealVector state_dist_cum = sM0; // prob of being in state given data forwards up to this point

		//  iteratively compute successive elements of salpha and scale_c
		RealVector temp;
		
		for(int t = 0; t < n_loci ; t++) {
			
			// log-probs of states emitting observed sequence at this mL
			RealVector bmL = mL_emissions[t].copy();
			
			// state_dist_cum is prob of being in various states given iterative data
			RealVector l_res_prob = helper.ebeLog(state_dist_cum);
			
			//ebe multiplication of emission and residence vectors
			temp = bmL.add(l_res_prob);
			
			// Set states with low res-prob to 0 resprob, otherwise
			// they mess numerics in next steps
			for(int iii = 0 ; iii < temp.getDimension(); iii++) {
				if(temp.getEntry(iii) < -675.){ temp.setEntry(iii, Double.NEGATIVE_INFINITY); }
			}
				
			// pull out common factor, scaling things reasonably for norm computation
			double lnc = getMinDefVal(temp); 
			temp.mapSubtractToSelf(lnc); 
			
			// compute norm and rescale
			temp = helper.ebeExp(temp);
			double t_norm = temp.getL1Norm(); 
			temp.mapDivideToSelf(t_norm); 		
			
			// compute scale factor c (as log) from norm and pulled factor
			lnc = lnc + Math.log(t_norm); 	
				
			// write values into holders
			scale_lnc[t] = lnc;
			salpha[t] = temp;
			
			// state distro for next site
			state_dist_cum = tM0.preMultiply(temp);
		}
				
		
		// BACKWARDS ALGORITHM
				
		// compute sbeta at the last time
		sbeta[n_loci - 1] = new ArrayRealVector(sM0.getDimension() , 1.0) ;
				
		// iteratively compute sbeta backwards
		for(int t = n_loci - 2 ; t >= 0 ; t--) {
			
			// log-probs of states emitting observed sequence at next mL
			RealVector bmL = mL_emissions[t+1].copy();
			
			// bring bmL to regular scale by dividing by scale factor, then bring into regular space
			bmL.mapSubtractToSelf(scale_lnc[t+1]);
			bmL = helper.ebeExp(bmL);
			
			// compute sbeta here, note we have already divided by scale factor
			sbeta[t]= tM0.operate(bmL.ebeMultiply(sbeta[t+1])) ;
		}
				
		return;
	}
	
	
	private RealVector mL_emitP(int mL_obs) {
		// mL_obs is the index in obs[x][mL_obs] corresp to column for emission
		
		// initialize
		RealVector out = new ArrayRealVector(eM0.getRowDimension()); 
		
		// obtain the number of each emission observed in the full metaLocus
		if(eM0.getColumnDimension() != obs.length) {System.err.print("number of emissions mismatch"); System.exit(0);}
		
		int[] mL_emission = new int[obs.length];
		for(int i = 0 ; i < mL_emission.length; i++) {
			mL_emission[i] = obs[i][mL_obs];
		}
		
		
		// 
		//int mL_sites_left = mL_size; // counter of how many sites in metaLocus still have to be counted
		
		for(int j = 0 ; j < mL_emission.length; j++ ) {
			
			
			// multiply by the log( probability ^ # observances )
			for(int s = 0 ; s < out.getDimension(); s++) {		// iterate over states
				double increment = 0;
				
				// if zero instances, having em0 entry as 0  means log(0) * 0 yields NaN
				// so instead just do nothing, (treat it as prob 1).
				// if not zero instances, then do normally, and if eM0 entry is 0
				// final val here will be -Inf
				if(mL_emission[j] != 0) {
					increment +=  mL_emission[j] * Math.log(eM0.getEntry(s, j));
				}
				out.setEntry(s, out.getEntry(s) + increment);
			}
			
			// readjust the sites left in mL
			//mL_sites_left = mL_sites_left - mL_emission[j];
		}
		
		
		// each index is log prob of that index's state emitting mL_obs for metalocus 
		return out;
	}

	
	private double getMinDefVal(RealVector in) {
		// return the minimum value that is well defined (not negative infinity)
		double out = Double.POSITIVE_INFINITY;
		for(int i = 0 ; i < in.getDimension() ; i++) {
			if(Double.isFinite(in.getEntry(i))) {
				out = Math.min(out, in.getEntry(i));
			}
		}
		return out;
	}
	
	
	public void computePosteriorDistros_ML(String output_file, double[] head_bins) throws IOException {
		

		double stime = System.nanoTime(); // just to time this method

					
		// initialize holders and fill salpha and sbeta vals
		RealVector[] salpha = new RealVector[n_loci];
		RealVector[] sbeta = new RealVector[n_loci];
		double[] scale_lnc = new double[n_loci];
				
		fbAB_sf_mL_fill(salpha, sbeta, scale_lnc);
				
		// compute and store the log_likelihoods			
		double likelihood = 0;
		for(int t = 0 ; t < n_loci ; t++) {
			likelihood += scale_lnc[t];
		}
		posteriorLL = likelihood;
		
		
		
		
		// setup File Writers
		FileWriter ofile= new FileWriter(output_file+ ".csv"); ofile.append(getHeader(head_bins)+" \n");
		mL_scribe scribe = new mL_scribe(ofile, mL_size);
		
		

		
		// compute gammas and add emission expectations
		for(int t = n_loci - 1 ; t >= 0 ; t--) {
			RealVector gamma = salpha[t].ebeMultiply(sbeta[t]);
			
			// print gamma for this locus
			scribe.print_mL_post(t, gamma);
	
		}
		
		
		
		
		ofile.flush(); ofile.close();
		System.out.println("\n <--> FB algorithm: " + Double.toString((System.nanoTime() - stime)/1000000000.) + " seconds\n");
		System.out.println("LL="+posteriorLL+"\n");
		
	}

	
	private class mL_scribe{
		FileWriter printer;
		int mlsize;
		
		mL_scribe(FileWriter fw, int ml_s){
			printer = fw;
			mlsize = ml_s;
		}
		
		void print_mL_post(int mL_num, RealVector post_dist) throws IOException {
			int middle_base_pos =  (int) ( (mL_num + 0.5) * mlsize );
			
			// print as a line in the file, using the avg position, and the avg posterior over stack
			printer.append(middle_base_pos+"");
			for(int i = 0 ; i < post_dist.getDimension(); i++) {
				printer.append(", " + post_dist.getEntry(i));
			}
			printer.append(" \n ");

		}
	}
	
	
	
	////////////////////////////////////////////////////////////
	
	
	/* HERE WE HAVE THE Expectation_FBfull class, that is deprecated, but would compute forward backwards without skipping loci
	 * public class Expectation_FBfull extends ExpectationsComputer{
	
	
	// data loaded in is stored here
	private int[] obs;
	private int n_loci;
	
	Expectation_FBfull(RealMatrix tp, RealMatrix ep, RealVector sd) {
		super(tp, ep, sd);
		// TODO Auto-generated constructor stub
	}

	///////////////////////////////////////////////////////////////////////////
	/// Abstract methods need to be implemented/ ////////
	/////////////////////////////////////////////////////////////////////
	@Override
	void loadDataFromFile(VcfReader v_file) {
		System.out.println("Expectation_FBfull class is DEPRECATED");
		
		//obs = v_file.getFullStream();
		//n_loci = obs.length;
		
	}

	@Override
	void computeExps() {
		
		initializeExp();
		
		// using obs, write in expectations for expTr and expEm (defined above)
		// note this is written so that we do not change the references for input matrices. we want to write to them 
		// in place, keeping their current references.
		
		RealMatrix teX = expTr.copy();
				
		RealMatrix[][] gx = gammaXiCompute_sf();
		RealMatrix[] gamma = gx[0];
		RealMatrix[] xi = gx[1];
		
		// expTr should store expected number of transitions from hidden state i to j , given the observed data
		// expEm should store expected number of emissions from hidden state i to observed d , given the observed data
		
		
		for(int t = 0 ; t < n_loci - 1 ; t++) {
			teX = teX.add(xi[t]);
			expEm.setColumnMatrix(obs[t], expEm.getColumnMatrix(obs[t]).add(gamma[t].transpose()));
		}
		expEm.setColumnMatrix(obs[n_loci-1], expEm.getColumnMatrix(obs[n_loci-1]).add(gamma[n_loci-1].transpose()));
		
		expSd.setSubVector(0, gamma[0].getRowVector(0)); 
		expTr.setSubMatrix(teX.getSubMatrix(0, teX.getRowDimension()-1, 0, teX.getColumnDimension()-1).getData(),0, 0);

	}

	
	///////////////////////////////////////////////////////////////////////////

	
	// COMPUTE AND RETURN GAMMA AND XI (SCALE FACTOR IMPLEMENTATION)
	private RealMatrix[][] gammaXiCompute_sf(){
			
			RealMatrix[][] gammaXi = new RealMatrix[2][n_loci];
			RealMatrix[] gamma = new RealMatrix[n_loci];
			RealMatrix[] xi = new RealMatrix[n_loci];
			
			RealVector[][] ab = fbAlphaBeta_sf();
			RealVector[] salpha = ab[0];
			RealVector[] sbeta = ab[1];
			RealVector[] scale_c = ab[2];
			
			// compute and store the log_likelihood			
			double likelihood = 0;
			for(int t = 0 ; t < n_loci ; t++) {
				likelihood += Math.log(scale_c[t].getEntry(0));
			}
			posteriorLL = likelihood;
			
		
			
			// compute gamma
			RealVector temp = new ArrayRealVector();
			for(int t = 0 ; t < n_loci ; t++) {
				temp = salpha[t].ebeMultiply(sbeta[t]);
				gamma[t] = new Array2DRowRealMatrix(temp.toArray()).transpose();
			}
			
			// compute xi
			RealMatrix temp1 = new Array2DRowRealMatrix(sM0.getDimension(), sM0.getDimension());  // only reason this doesn't need to go inside loop over time is because "ScalarMultiply" method below returns a new matrix and doesnt operate in place
			for(int t = 0 ; t < n_loci - 1 ; t++) {
				// compute xi(t) row by row. 
				for(int r = 0 ; r < temp1.getRowDimension() ; r++) {
					temp = tM0.getRowVector(r).ebeMultiply(eM0.getColumnVector(obs[t+1]).ebeMultiply(sbeta[t+1]));
					temp.mapMultiplyToSelf(salpha[t].getEntry(r));
					temp1.setRowVector(r, temp);
				}
				xi[t] = temp1.scalarMultiply(1.0 / scale_c[t+1].getEntry(0)) ;
			}
			
			// assign gamma and xi to holder and return
			gammaXi[0]=gamma;
			gammaXi[1]= xi;
			return gammaXi;
			
	}
	
	// COMPUTE FORWARDS AND BACKWARDS SCALED PDFS (include scale factors)
	private RealVector[][] fbAlphaBeta_sf(){
			
			RealVector[][] abc = new RealVector[3][];
			RealVector[] salpha = new RealVector[n_loci];
			RealVector[] sbeta = new RealVector[n_loci];
			RealVector[] scale_c = new RealVector[n_loci];
			
			// first compute scale_c(0) 
			RealVector temp0 = eM0.getColumnVector(obs[0]).ebeMultiply(sM0);
			scale_c[0] = new ArrayRealVector().append(temp0.getL1Norm());
			
			// next compute salpha(0). you can use temp0 for this
			salpha[0] = temp0.mapDivideToSelf(scale_c[0].getEntry(0));
			
			//  iteratively compute successive elements of salpha and scale_c
			double temp1 = 0;
			for(int t = 1; t < n_loci ; t++) {
				temp0 = eM0.getColumnVector(obs[t]).ebeMultiply(tM0.preMultiply(salpha[t-1]));
				temp1 = temp0.getL1Norm();
				temp0.mapDivideToSelf(temp1);
				
				scale_c[t] = new ArrayRealVector().append(temp1);
				salpha[t] = temp0;
			}
			
			
			// compute sbeta at the last time
			double[] temp2 = new double[sM0.getDimension()];
			Arrays.fill(temp2, 1.0);
			sbeta[n_loci - 1] = new ArrayRealVector(temp2) ;
			
			// iteratively compute sbeta backwards
			for(int t = n_loci - 2 ; t >= 0 ; t--) {
				sbeta[t]= tM0.operate(eM0.getColumnVector(obs[t+1]).ebeMultiply(sbeta[t+1])).mapDivideToSelf(scale_c[t+1].getEntry(0))  ;
			}
			
			abc[0] = salpha;
			abc[1] = sbeta;
			abc[2] = scale_c;
			
			return abc;
	}
		
			
	
}
	 * */



}