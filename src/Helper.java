import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.vcf.VCFFileReader ;

public class Helper {

	public static void main(String[] args) throws IOException {
		Helper help = new Helper();
		
//		String a = "e";
//		String b = "ba";
		System.out.println(Arrays.toString(help.generate_log_discretization(10,100,4)));
//
//		System.out.println(-13 % 10);
		
		System.out.println (Arrays.toString ((help.getPsmcPartition (8, 50, 0.1 * 1d/45)).toArray()));
	}

	
	int nChoose2(int n) {
		if(n < 2) {return 0;}
		else{return n * (n-1) / 2 ;}
	}
	

	// CLASS TO EFFICIENTLY LOOKUP N CHOOSE K DIRECTLY FROM PASCALS TRIANGLE
	class ChooseComputer{
		int[][] choose_table;
		int max_n;
		
		ChooseComputer(int i){
			max_n = i;
			build_choose_table();
		}
		
		void build_choose_table(){
			choose_table = new int[max_n + 1][max_n + 1];
			choose_table[0][0] = 1;
			for(int i = 1 ; i < max_n + 1 ; i++) {
				choose_table[i][0] = 1 ;
				for(int j = 1 ; j < i+1 ; j++) {
					choose_table[i][j] = choose_table[i-1][j-1] + choose_table[i-1][j];
				}
			}
		}
		
		double get(int n, int k){
			
			if(n==k) {
				return 1;
			}
			
			if(k > n ||  k < 0) {
				return 0;
			}
			
			
			
			else {
				return choose_table[n][k];
			}
			
		}
		
		double get1(int n, int k){
		
			if(k > n ||  k < 0) {
				return 1;
			}
			
			else {
				return choose_table[n][k];
			}
			
		}
		
	}
	
	
	// convergence criteria for when we integrate ode to infinity to obtain distro across absorbing states at infinity
	class Converge implements EventHandler {
			
			double precision;
			int[] aStates; // absorbing states
			
			Converge(double p, int[] aSs){
				precision = p;
				aStates = aSs;
			}
			
			@Override
			public Action eventOccurred(double arg0, double[] arg1, boolean arg2) {
				return EventHandler.Action.STOP;
			}

			@Override
			public double g(double arg0, double[] arg1) {

				double arp = 0; //absorbing residence probability
				for(int i = 0 ; i < aStates.length ; i++) {
					arp = arp + arg1[aStates[i]];
				}
				
				return (1 - precision) - arp;
			}

			@Override
			public void init(double arg0, double[] arg1, double arg2) {
				
			}

			@Override
			public void resetState(double arg0, double[] arg1) {
				
			}
			
	}
		
	double[][] odeSolsAtTimes(FirstOrderDifferentialEquations ode, double[] initial_rps, double[] times, double t_prec, boolean till_inf , int[] abs_indices) {
		
		/*
		 * return double[][] where first index is state, second is the discrete times contained in times at which residence prob is evaluated
		 */
		
		int n_cols = times.length+1;
		if(!till_inf) {n_cols--;}
		
		double[][] out_res_probs = new double[ode.getDimension()][n_cols];
		
		
		FirstOrderIntegrator integrator = new DormandPrince853Integrator(t_prec, 100.0, 1.0e-12, 1.0e-12); // adaptive integrator. better

		double[] residence_probs = initial_rps;
		
		double time_counter = 0.;				
		
		
//		/// time some stuff
//		double stime = System.nanoTime();

		
		for(int t = 0 ; t < times.length; t++) {
			
			if( (times[t] - time_counter) <= t_prec) {} // do nothing if time step is too small
			else {
				integrator.integrate(ode, time_counter, residence_probs, times[t], residence_probs);
				time_counter = times[t];
			}
			

			
//			System.out.println("current time: " + Double.toString((System.nanoTime() - stime)/1000000000.) + " seconds");
			
			
			
			
			// store value of residence_probs at time time_counter for state1 and state2
			for(int s = 0 ; s < ode.getDimension() ; s++) {
				out_res_probs[s][t] = residence_probs[s];
			}	
			
		}
		
		
		if(till_inf) {
			
			Converge converge_condition = new Converge(.0000000001, abs_indices);

			// if we are already converged at the last time, we just take the current state again
			// otherwise, we need to run until enough convergence
			if (converge_condition.g(time_counter, residence_probs) > 0) {
				// not yet converged, need to go some more
				integrator.addEventHandler(converge_condition, 100., .01 , 1000000);
				integrator.integrate(ode, time_counter, residence_probs, 100000, residence_probs); 
			}
			
			// save whatever the right distribution for infitnity (converged) is
			for(int s = 0 ; s < ode.getDimension() ; s++) {
				out_res_probs[s][times.length] = residence_probs[s];
			}	
			
		}
		
		

//		System.out.println("current time: " + Double.toString((System.nanoTime() - stime)/1000000000.) + " seconds");
//
//		for (int timeIdx=n_cols-3; timeIdx<n_cols; timeIdx++) {
//			double[] print_array = new double[ode.getDimension()];
//			for (int stateIdx=0; stateIdx<print_array.length; stateIdx++) {
//				print_array[stateIdx] = out_res_probs[stateIdx][timeIdx];	
//			}
//			System.out.println(Arrays.toString(print_array));
//		}
		

		return out_res_probs;
		
	}
	
	
	
	////////////////////////////////////////
	// reading and writing data-streams into temp files
	///////////////////////
	
	void write_dataStreamToCSV(String temp_file, int[][] d_stream){
		double stime = System.nanoTime();

		
		FileWriter scribe;
		try {
			scribe = new FileWriter( temp_file+ ".csv");
			int n_sites = d_stream[0].length;
			int n_params = d_stream.length;
				

			for(int s = 0 ; s < n_sites ; s++) {
				for(int p = 0 ; p < n_params ; p++) {
					scribe.append(d_stream[p][s] + ", ");
				}
				scribe.append("\n");
			}
			
			scribe.flush();
			scribe.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.print("Failed to write datastream to file \n");
			System.exit(0);
		}
		
		
		System.out.println(Double.toString((System.nanoTime() - stime)/1000000000.) +" sec to write datastream " + temp_file + " to temp csv");

	}
	
	int[][] read_dataStreamFromCSV(String temp_file) {
		int[][] out = new int[][] {};
		double stime = System.nanoTime();

		try {
			// initialize scanners
			Scanner scanner = new Scanner(new File(temp_file+".csv"));
			Scanner line_scan;

			// figure out how many columns this has, so we can create that many holder arraylists.
			String first_line = scanner.nextLine();
			int num_columns = 0;
			line_scan = new Scanner(first_line); line_scan.useDelimiter("\\s*,\\s*");
			while(line_scan.hasNext()) { line_scan.nextInt(); num_columns++; }
			
			// close and reset scanners
			line_scan.close(); scanner.close();
			scanner = new Scanner(new File(temp_file+".csv")); // reset the buffer/scanner to start of file
			
			// Fill holder structures with the data from the file going line by line
			ArrayList<ArrayList<Integer>> data = new ArrayList<ArrayList<Integer>>(num_columns);
			for(int i = 0 ; i < num_columns ; i++) {data.add(new ArrayList<Integer>());}
			while(scanner.hasNextLine()) {
				// write in one line at a time
				line_scan = new Scanner(scanner.nextLine());  line_scan.useDelimiter("\\s*,\\s*");
				// process the line
				for(int col = 0 ; col < num_columns; col++) {
					data.get(col).add(line_scan.nextInt());
				}
				line_scan.close();
			}
			
			// copy from holder structures into the out object
			out = new int[num_columns][data.get(0).size()];
			for(int i = 0 ; i < num_columns ; i++) { 
				ArrayList<Integer> tt_array = data.get(i);
				for(int j = 0 ; j < tt_array.size() ; j++) {
					out[i][j] = tt_array.get(j);
				}
			}
			
			line_scan.close();
			scanner.close();
					
		} 
		
		catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.print("Could not read dataStream from the file"); System.exit(0);
		}

		
		System.out.println(Double.toString((System.nanoTime() - stime)/1000000000.) +" sec to read datastream " + temp_file + " from temp csv");
		
		return out;
	}
	
	void delete_csvFromDisc(String temp_file) {
		File file = new File(temp_file + ".csv");
		file.delete();
	}
	
	
	// does exact same thing as the methods above, just handles doubles instead of ints
	
	double[][] read_doubleMat_FromCSV(String temp_file){
		double[][] out = new double[][] {};
		double stime = System.nanoTime();

		try {
			// initialize scanners
			Scanner scanner = new Scanner(new File(temp_file+".csv"));
			Scanner line_scan;

			// figure out how many columns this has, so we can create that many holder arraylists.
			int num_columns = 0;
			int header_size = 0;
			line_scan = new Scanner(scanner.nextLine()); line_scan.useDelimiter("\\s*,\\s*");
			
			// first loop finds size of header, second finds number of columns
			while(!line_scan.hasNextDouble()) {line_scan = new Scanner(scanner.nextLine()); line_scan.useDelimiter("\\s*,\\s*"); header_size++;}
			while(line_scan.hasNext()) { line_scan.next(); num_columns++; }
			
			// close and reset scanners
			line_scan.close(); scanner.close();
			scanner = new Scanner(new File(temp_file+".csv")); // reset the buffer/scanner to start of file
			for(int i = 0 ; i < header_size ; i++) {scanner.nextLine();} // fast-forward through header
			
			
			// Fill holder structures with the data from the file going line by line
			ArrayList<ArrayList<Double>> data = new ArrayList<ArrayList<Double>>(num_columns);
			for(int i = 0 ; i < num_columns ; i++) {data.add(new ArrayList<Double>());}
			while(scanner.hasNextLine()) {
				// write in one line at a time
				String nl = scanner.nextLine().trim();
				if(nl.length() != 0) {
					line_scan = new Scanner(nl);  line_scan.useDelimiter("\\s*,\\s*");
					// process the line
					for(int col = 0 ; col < num_columns; col++) {
						double next_element = line_scan.nextDouble(); 
						data.get(col).add(next_element);
					}
					line_scan.close();
				} else {}
			}
			
			// copy from holder structures into the out object
			out = new double[num_columns][data.get(0).size()];
			for(int i = 0 ; i < num_columns ; i++) { 
				ArrayList<Double> tt_array = data.get(i);
				for(int j = 0 ; j < tt_array.size() ; j++) {
					out[i][j] = tt_array.get(j);
				}
			}
			
			line_scan.close();
			scanner.close();
					
		} 
		
		catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.print("Could not read dataStream from the file"); System.exit(0);
		}

		
		System.out.println(Double.toString((System.nanoTime() - stime)/1000000000.) +" sec to read datastream " + temp_file + " from temp csv");
		
		return out;
	}
	
	void write_doubleMat_ToCSV(String temp_file, double[][] d_stream) {
		double stime = System.nanoTime();

		
		FileWriter scribe;
		try {
			scribe = new FileWriter( temp_file+ ".csv");
			int n_sites = d_stream[0].length;
			int n_params = d_stream.length;
				
			for(int s = 0 ; s < n_sites ; s++) {
				scribe.append(d_stream[0][s]+"");
				for(int p = 1 ; p < n_params ; p++) {
					scribe.append( ", " + d_stream[p][s]);
				}
				scribe.append("\n");
			}
			
			scribe.flush();
			scribe.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.print("Failed to write datastream to file");
		}
		
		
		System.out.println(Double.toString((System.nanoTime() - stime)/1000000000.) +" sec to write datastream " + temp_file + " to temp csv");
		
	}

	////////////////////////////////////////

	
	RealVector mat_rowsums(RealMatrix mat) {
		
		double[] ones = new double[mat.getColumnDimension()];
		Arrays.fill(ones, 1.);
		
		return mat.operate(new ArrayRealVector(ones));
	}
	
	RealVector mat_colsums(RealMatrix mat) {
		
		double[] ones = new double[mat.getRowDimension()];
		Arrays.fill(ones, 1.);
		return mat.preMultiply(new ArrayRealVector(ones));
		
	}

	double[] generate_log_discretization(double lb ,  double rb , int n ){
		double[] times = new double[n-1];
		times[0] = lb;
		times[n-2] = rb;
		
		double a = Math.log(rb/lb);
		for(int i = 1 ; i < n-2 ; i++) {
			times[i] = Math.exp(i * a/(n-2))*lb;
		}
		
		return times;
	}

	double[] generate_uniform_discretization(double lb, double rb, int n){
		double[] times = new double[n-1];
		for(int i = 0 ; i < n-1 ; i++) {
			times[i] = lb +  i * (rb - lb) / (n-2);
		}
		return times;
	}
	
	
	int placeInPartition(double t, double[] dts) {
		
		if(t<0.) {
			System.out.println("length invalid");
			System.exit(0);
		}
		for(int i = 0 ; i < dts.length; i++) {
			if(t < dts[i]) {
				return i;
			}
		}
		return dts.length;
		
	}
	
	RealMatrix pdfToRowProbMatrix(RealMatrix pdf, int type) {
		// type determines how we handle 0 prob rows. 
		// type = 0 means transition matrix
		// type = 1 means emission matrix
		
		RealMatrix out = pdf.copy();
		
		
		///////////////////////////
		//// USE MATTHIAS'S UNDERFLOW PREVENTION SCHEME
		// ENDS UP CAUSING NON-INVERTIBLE MATRIX ERRORS DOWN THE ROAD FOR SOME CASES
		double eps = 1e-12;
		for(int i = 0 ; i < out.getRowDimension(); i++) {
			for(int j = 0 ; j < out.getColumnDimension(); j++ ) {
				if (out.getEntry(i, j) < eps) {
					out.setEntry(i, j, eps);
				}
			}
		}
		out = out.scalarMultiply(1./this.totalEntriesMat(out));
		//////////////////////////
		
		

		
		// make sure pdf is a proper pdf
		double tot = this.totalEntriesMat(out);
		if(Math.abs(tot - 1.0 ) > 0.00000001) {
			System.out.println("Warning: PDF Matrix sum = " + tot); 
		}
		
		boolean underflow = false;

		for(int i = out.getRowDimension() - 1 ; i > -1 ; i--) {
			RealVector r = out.getRowVector(i);
			double r_sum = r.getL1Norm();
			
			if(r_sum == 0 ) { // how to handle 0 prob rows, where we never see events
				
				underflow = true;
				
				/*
				//THIS SEGMENT DOES NOT FIX UNDERFLOW PROBLEM BECAUSE RESULTING PROB MATRIX IS SINGULAR AND WONT INVERT. NEED TO USE THE FIX FURTHER ON
				// following scheme is patch fix for case where the first few rows have underflow (not for if the last few have it)
				// set all rows under this one to be last informative one. (mainly used for very recent states where the chmm might never spend time)
				r = out.getRowVector(i + 1);
				for(int ii = i ; ii > -1; ii--) {
					out.setRowVector(ii, r.copy());
				}
				*/		

				// PATCH FIX FOR IF THERE IS UNDERFLOW. ASSUME THAT IT DOES NOT TRANSITION FROM THIS STATE, AND EMITS ONLY ZEROS FROM THIS STATE
				if(type == 0) {	out.setEntry(i, i, 1.);}
				else if (type == 1) {out.setEntry(i, 0, 1.);}
				else {System.err.print("Invalid prob_matrix type"); System.exit(0);}
			}
			else {
				r.mapDivideToSelf(r_sum); 
				out.setRowVector(i, r);
			}
		}
		if(underflow) {	System.err.println("ERROR: UNDERFLOW IN STATES"); }

		return out;
	}
	
	double totalEntriesMat(RealMatrix m) {
		RealVector t1 = new ArrayRealVector(m.getColumnDimension(), 1.);
		RealVector t2 = new ArrayRealVector(m.getRowDimension(), 1.);
		
		return m.operate(t1).dotProduct(t2);
	}

	int roll (RealVector row) {
		if ( Math.abs( row.getL1Norm() - 1.0) > .000000001) {
			System.out.print("this isnt a probability row vector (doesnt sum to 1), sums to ");
			System.out.println(row.getL1Norm());
			System.exit(0);
		}		
		int selection = 0;
		double frac = Math.random();
		while(frac > row.getEntry(selection)) {
			frac = frac - row.getEntry(selection);
			selection++;
		}
		return selection;	
	}

	
	double[] ebeExp(double[] inp) {
		double[] t = new double[inp.length];
		for(int i = 0 ; i < inp.length ; i++) {
			t[i] = Math.exp(inp[i]);
		}
		return t;
	}
	
	RealVector ebeExp(RealVector in) {
		RealVector out = in.copy();
		for(int i = 0 ; i < out.getDimension(); i++) {
			out.setEntry(i, Math.exp(in.getEntry(i)));
		}
		return out;
	}
	
	double[] ebeLog(double[] in) {
		double[] t = new double[in.length];
		for(int i = 0 ; i < t.length ; i++) {
			t[i] = Math.log(in[i]);
		}
		return t;
	}

	RealVector ebeLog(RealVector in) {
		RealVector out = in.copy();
		for(int i = 0 ; i < out.getDimension(); i++) {
			out.setEntry(i, Math.log(in.getEntry(i)));
		}
		return out;
	}
	
	RealMatrix ebeMultiply(RealMatrix a , RealMatrix b) {
		RealMatrix out = new Array2DRowRealMatrix(a.getRowDimension(), a.getColumnDimension());
		for(int i = 0 ; i < out.getRowDimension(); i++) {
			for(int j = 0 ; j < out.getColumnDimension(); j++) {
				out.setEntry(i, j, a.getEntry(i, j) * b.getEntry(i, j));
			}
		}
		return out;
	}

	
	RealMatrix ebeDivide(RealMatrix a , RealMatrix b) {
		RealMatrix out = new Array2DRowRealMatrix(a.getRowDimension(), a.getColumnDimension());
		for(int i = 0 ; i < out.getRowDimension(); i++) {
			for(int j = 0 ; j < out.getColumnDimension(); j++) {
				out.setEntry(i, j, a.getEntry(i, j) / b.getEntry(i, j));
			}
		}
		return out;
	}

	
	RealMatrix stdErrorMatComp(RealMatrix a , RealMatrix b, int trials) {
		return ebeDiv(a.subtract(b), compStdErr(a,trials));
	}
	
	RealVector stdErrorVecComp(RealVector a , RealVector b , int trials) {
		return a.subtract(b).ebeDivide(compStdErr(a, trials));
	}
	
	private RealMatrix compStdErr(RealMatrix m, int n) {
		RealMatrix stderr = m.copy();
		for(int i = 0 ; i < stderr.getRowDimension(); i++) {
			for(int j = 0 ; j < stderr.getColumnDimension(); j++) {
				double err = stderr.getEntry(i, j);
				err = Math.sqrt( Math.pow(1-err, 2)*err + Math.pow(err, 2)*(1-err) );
				stderr.setEntry(i, j, err/Math.sqrt(n));
			}
		}
		
		return stderr;
	}
	
	private RealVector compStdErr(RealVector m, int n) {
		RealVector stderr = m.copy();
		for(int i = 0 ; i < stderr.getDimension(); i++) {
				double err = stderr.getEntry(i);
				err = Math.sqrt( Math.pow(1-err, 2)*err + Math.pow(err, 2)*(1-err) );
				stderr.setEntry(i, err/Math.sqrt(n));
		}
		
		return stderr;
	}
	
	private RealMatrix ebeDiv(RealMatrix n , RealMatrix d) {
		if(n.getColumnDimension() != d.getColumnDimension()) {System.exit(0);}
		if(n.getRowDimension() != d.getRowDimension()) {System.exit(0);}

		
		RealMatrix q = n.copy();
		for(int i = 0 ; i < q.getRowDimension() ; i++) {
			for(int j = 0 ; j < q.getColumnDimension() ; j++) {
				q.setEntry(i, j, n.getEntry(i, j) / d.getEntry(i, j));
			}
		}
		return q;
	}
	
	
	double[] partitionProbDistributor(double min_bin_prob, int num_states, boolean break_last) {
		
		// log uniform distro
		double lb = Math.log(min_bin_prob);
		double rb = 0.;
		
		double[] probs = new double[num_states];
		probs[0] = lb;
		for(int i = 1 ; i < num_states - 1 ; i++) {
			probs[i] = probs[i-1] + (rb - lb)/(num_states-1) ;
		}
		probs = this.ebeExp(probs);
		
		
		// Convert CDF to PDF 
		for(int i = probs.length-1 ; i > 0 ; i--) {
			probs[i] = probs[i] - probs[i-1];
		}
		double[] out_p = probs;
		
		// break up the last pdf bin into pieces the size of the next to last
		if(break_last) {
			double chopped_size = .15; //probs[probs.length - 2]; 
			int extra_bins = (int) ( probs[probs.length - 1] / chopped_size - 1E-10); // the 1E-10 ensures that there is not a tiny trailing bin if it just barely clears the floor
			
			out_p = new double[probs.length + extra_bins];
			for(int i = 0 ; i < probs.length; i++) {out_p[i] = probs[i];}
			for(int i = probs.length - 1 ; i < out_p.length - 1; i++) {
				out_p[i] = chopped_size;
			}
			out_p[out_p.length-1] = 1.0 - (new ArrayRealVector(out_p).getL1Norm());
		}
		
		
		
		return out_p;
	}
	
	// a psmc/msmc2 inspired partition in coalescent time, akin to the positive part of an exponential function - 1
	// psmc-default: t_max=15, recency_factor=0.1
	// msmc2 divides recency_factor by the number of pairs
	RealVector getPsmcPartition (int nr_intervals, double t_max, double recency_factor) {
		
		
		// we omit the first one, because it would be zero
		double[] ts = new double[nr_intervals];
		
		for (int i=0; i<ts.length; i++) {
			// reindex to omit first one
			double x = i+1;
			
			// get the fractions between 0 and 1
			x /=nr_intervals;

			// get the exponent
			x *= Math.log (1 + t_max/recency_factor);
			
			// and exponentiate
			ts[i] = recency_factor*Math.exp(x) - recency_factor;
		}

		return new ArrayRealVector (ts);
	}
	
	RealVector getTreeHeightPartitions(double[] bin_probs, int samples) {
		// for constant population, computes the tree height/length partitions according to bin_probabilities 
		// in units of coalescent time
		
		// CHECK FOR VALID PARAMS
		double ss = 0;
		for(int i = 0 ; i < bin_probs.length; i++) {
			ss = ss + bin_probs[i];
			if(bin_probs[i] < 0) {System.err.print("Distribution Probabilities must be positive"); System.exit(0);}
		}
		if(Math.abs(ss - 1) > .0000000001 ) {System.err.print("Distribution Probabilities must sum to one."); System.exit(0);}
		
		
		int resolution = 10000; // number of internal points within 1 coalescent time unit
		// X DOMAIN DISCRETE VALS
		double[] x_vals = new double[resolution * 10 + 1];
		for(int i = 0 ; i < 10 * resolution + 1 ; i++) { x_vals[i] = 1.0 * i / resolution;	}
		
		
	
		// FILL THE Y PDF VALS
		double[] pdf_vals = new double[x_vals.length]; 
		
		// need to compute coefficients. do it once now to avoid repeating
		double[] c_nk = new double[samples+1];
		for(int k = 2 ; k < samples + 1 ; k++) {
			double val = Math.log(2*k - 1);
			for(int kk = 1 ; kk < k+1 ; kk++) {
				val = val + Math.log(1 - (k - 1.) / (samples + kk - 1.) );
			}
			c_nk[k] = Math.exp(val);
		}
		
		
		// for all x points compute pdf by summing over convolution distro of exponentials
		for(int i = 0 ; i < pdf_vals.length ; i++) { // given this x val compute pdf
			double expx = Math.exp(x_vals[i]);
			
			double pdft = 0;
			for(int k = 2 ; k < samples + 1 ; k++) {
				double rate = this.nChoose2(k);
				double summand = Math.pow(-1., k) * rate * Math.pow(expx, -rate) * c_nk[k];
				pdft = pdft + summand;
			}
			
			pdf_vals[i] = pdft;
		}
		
		
		// compute integral (cdf) by trapezoidal method
		double[] cdf_vals = new double[x_vals.length]; 
		double integral_counter = 0;
		for(int i = 1 ; i < cdf_vals.length ; i++) { // given this x val compute pdf
			double summand = (pdf_vals[i] + pdf_vals[i-1])/(2.0 * resolution);
			integral_counter = integral_counter + summand;
			cdf_vals[i] = integral_counter; 	
		}
	
		
		// FIND PARTITION LOCATIONS
		double[] partition_poss = new double[bin_probs.length - 1] ;
		integral_counter = 0;
		for(int p = 0 ; p < partition_poss.length ; p++) {
			integral_counter = integral_counter + bin_probs[p];
			int i = 0;
			while(cdf_vals[i]<=integral_counter) {i++;}
			partition_poss[p] = x_vals[i];
		}
		
		
		
		// RETURN PARTITIONS
		return new ArrayRealVector(partition_poss);
		
		
	}
	
	RealVector getTreeLengthPartitions(double[] bin_probs, int samples) {
		// for constant population, computes the tree height/length partitions according to bin_probabilities 
		// in units of coalescent time
		
		// CHECK FOR VALID PARAMS
		double ss = 0;
		for(int i = 0 ; i < bin_probs.length; i++) {
			ss = ss + bin_probs[i];
			if(bin_probs[i] < 0) {System.err.print("Distribution Probabilities must be positive"); System.exit(0);}
		}
		if(Math.abs(ss - 1) > .0000000001 ) {System.err.print("Distribution Probabilities must sum to one."); System.exit(0);}
		
		
		int resolution = 1000; // number of internal points within 1 coalescent time unit
		
		// X DOMAIN DISCRETE VALS
		// estimate range based on empiric heuristic
		int distro_peak_estimate = (int) (2.4 * Math.log(samples)/Math.log(3) - 0.9);
		int x_lb =  Math.max(distro_peak_estimate - 5, 0);
		int x_rb =  distro_peak_estimate + 16;
		int dom_span = x_rb - x_lb;
		
		
		double[] x_vals = new double[resolution * dom_span + 1];
		for(int i = 0 ; i < x_vals.length ; i++) { x_vals[i] = x_lb + 1.0 * i / resolution;	}
		x_vals[0] = 0; // set this for completeness
		
	
		// FILL THE Y PDF VALS
		double[] pdf_vals = new double[x_vals.length]; 
	
		
		
		// for all x points compute pdf 
		for(int i = 0 ; i < pdf_vals.length ; i++) { // given this x val compute pdf
			double expx = Math.exp(- x_vals[i] / 2.);
			
			double pdfv = (samples - 1)/2.;
			pdfv = pdfv * Math.pow(1 - expx, samples - 2) * expx ;
			pdf_vals[i] = pdfv;
		}
		
		
		// compute integral (cdf) by trapezoidal method
		double[] cdf_vals = new double[x_vals.length]; 
		double integral_counter = 0;
		for(int i = 1 ; i < cdf_vals.length ; i++) { // given this x val compute pdf
			double summand = (pdf_vals[i] + pdf_vals[i-1])/(2.0 * resolution);
			integral_counter = integral_counter + summand;
			cdf_vals[i] = integral_counter; 	
		}
	
		
		// FIND PARTITION LOCATIONS
		double[] partition_poss = new double[bin_probs.length - 1] ;
		integral_counter = 0;
		for(int p = 0 ; p < partition_poss.length ; p++) {
			integral_counter = integral_counter + bin_probs[p];
			int i = 0;
			while(cdf_vals[i]<=integral_counter) {i++;}
			partition_poss[p] = x_vals[i];
		}
		
		
		
		// RETURN PARTITIONS
		return new ArrayRealVector(partition_poss);
		
		
	}
	
	
	class dataStreamsHandler {
		
		private Random random = new Random();
	
		
		// switch indicating whether streams are stored in memory (faster, more RAM) or in files (slower, less RAM)
		private boolean streams_in_mem = false; 
	
		int num_streams;
		
		// datastream holder
		private int[][][] chromosome_streams = null; // array of streams if streams_in_mem is true

		// array of filenames with datastreams for each chromosome/subsampling
		private String[] data_stream_files = null;
		
		/// Average Watterson's estimate for theta across all streams
		private double average_Wattersons_theta;
		
		// total length of non-missing genome across all chromosomes and all subsamplings (ie total length of n-sample genome acted over as independent)
		private double total_genetic_length;
		
		private int[][] haps_indices;
		
		dataStreamsHandler(ArrayList<String> file_names, ArrayList<String> reference_names, ArrayList<String> ancestral_names,  ArrayList<Integer> ls, ArrayList<int[]> hap_index_lists, int max_tract_size, int mL_size, boolean binary_emission, String temp_prefix){
			// filenames, ref_names, anc_names, ls, sample_intervals, 
			// each are array of equal length, and define the corresponding field for a datastream that we use for analysis
			
			// max tract size and mL_size and n_samples are also important to define structure of data streams
			
			// target_folder is null if we are using memory, otherwise it is name of temporary file prefix
		
			
			
			num_streams = file_names.size();
			if(num_streams == 0) {System.out.println("No Files Provided"); System.exit(0);}
			

			// We compute Watersons estimator for each datastream
			double wattersons_sum = 0;
			
			// check whether to use memory or files for the streams
			if(temp_prefix == null) { streams_in_mem = true;}
			
			// initialize stream structures
			if(streams_in_mem) {chromosome_streams = new int[num_streams][][] ;}
			else { 	data_stream_files = new String[num_streams];  }
			
			//initialize the haps indices list for above streams
			haps_indices = new int[num_streams][];
			
			// iterate through each stream and store appropriately. also scan to compute waterson's estimate
			for(int s = 0 ; s < num_streams ; s++) {
				// make sure that the sample_index_intervals are consistent with the num-samples
				int[] subgroup_indices = hap_index_lists.get(s);
				
				VcfReader temp = new VcfReader(file_names.get(s), reference_names.get(s), ancestral_names.get(s), ls.get(s), subgroup_indices);
				temp.update_max_tract_size(max_tract_size);
				int[][] stream = temp.getStream(mL_size);
				
				if(binary_emission) {stream = adjustBinaryEmission(stream, mL_size > 1);}
				
				// store stream into appropriate structure
				if(streams_in_mem) {
					chromosome_streams[s] = stream;
				}
				else {
					data_stream_files[s] = temp_prefix + "_datastream_" + s;
					write_dataStreamToCSV(data_stream_files[s], stream);
				}
				haps_indices[s] = subgroup_indices;
				
				// scan stream, count total number seg sites and total sites processed
				int[] ss_frac = getFracSS(mL_size, stream, binary_emission);
				// compute watterson's theta, and add to running total
				wattersons_sum = wattersons_sum + estimate_watersons_theta(1.0*ss_frac[0]/ss_frac[1], subgroup_indices.length)*ss_frac[1];
				total_genetic_length = total_genetic_length + ss_frac[1];
			}
			
			average_Wattersons_theta = wattersons_sum / total_genetic_length;
		}
		
		
		protected double estimate_watersons_theta(double f_ss, int n_s) {
			// find wattersons theta, given fraction of segregating sites computed already
			// just divides f_ss by harmonic sum
			
			double harmonic = 0;
			for(int i = 1 ; i < n_s ; i++) {
				harmonic = harmonic + 1.0/i ;
			}
					
			return f_ss/harmonic;
		}
		
		
		double get_Watterson() { // get AVG Watterson's theta across all streams
			return average_Wattersons_theta;
		}
		
		double getGenLen() {
			return total_genetic_length;
		}
		
		int[][] getStream(int s){
			if(s >= num_streams) {
				System.out.println("Trying to access invalid stream.");
				System.exit(0);
			}
			
			if(streams_in_mem) {
				return chromosome_streams[s];
			}
			else {
				return read_dataStreamFromCSV(data_stream_files[s]);
			}
		}
		
		int get_ns(int s) {
			// return n_s for the s'th stream
			return haps_indices[s].length;
		}
		
		void deleteStreamFiles() {
			if(!streams_in_mem) {
				for(int f = 0 ; f < num_streams ; f++) {
					delete_csvFromDisc(data_stream_files[f]);
				}
			}

		}
		
		
		int[] get_ns_list() { // returns set of unique n_s values as specified by the hap_index_lists
			ArrayList<Integer> ns_list = new ArrayList<Integer>();
			for(int i = 0 ; i < haps_indices.length; i++) {
				if (!ns_list.contains(new Integer(haps_indices[i].length))) {
					ns_list.add(new Integer(haps_indices[i].length));
				}
			}
			// copy to int[] instead of arrayvector
			int[] out = new int[ns_list.size()];
			for(int i = 0 ; i < out.length; i++) {
				out[i] = ns_list.get(i);
			}
			
			return out;
		}
		
		
		
	}
	
	int[][] adjustBinaryEmission(int[][] stream, boolean ml){
		// This adjusts DATASTREAMS to be appropriate for a binary emission implementation 
		// where each site on genome emits "segregating (0)" or "non_segregating (1)"
		// NB original stream object is changed IN PLACE for non-MetaLocus implementation 
		// and for ML the shape is modified and copied to new structure
		int[][] out = null;
		
		
		// first handle stream in non-metalocus implementation (non-seg-site skipping stream)
		if(!ml) {
			out = stream;
			for(int i = 0 ; i < stream[0].length ; i++) {
				if(out[0][i]>=1) {out[0][i]=1;}
			}
		}
		
		// now handle metalocus case
		else {
			out = new int[2][stream[0].length] ;
			
			
			
			for(int l = 0 ; l < out[0].length ; l++) {
				out[0][l] = stream[0][l];
				
				int num_ss = 0;
				for(int ss = 1 ; ss < stream.length; ss++) {
					num_ss = num_ss + stream[ss][l];
				}
				
				// aggregate num of seg sites gets written to the second row
				out[1][l] = num_ss;
			}
		}
		
		
		return out;
		
	}

	RealMatrix adjustBinaryEmission(RealMatrix em_prob) {
		// This adjusts Emission Probability Matrix to be appropriate for a binary emission implementation 
		// where each site on genome emits "segregating (0)" or "non_segregating (1)"
		// note that the Matrix shape is altered, so a new structure is returned
		
		RealMatrix out = new Array2DRowRealMatrix(em_prob.getRowDimension(), 2);
		for(int i = 0 ; i < out.getRowDimension(); i++) {
			double prob_nss = em_prob.getEntry(i, 0);
			
			out.setEntry(i, 0, prob_nss);
			out.setEntry(i, 1, 1 - prob_nss);
		}
		
		return out;
	}
	
	int[] getFracSS(int mls, int[][] stream, boolean binary_emission) {
		// input is metaLocus_size, number of samples, and datastream
		int num_tracts = stream[0].length;		
		int nss = 0; // number of segregating sites in stream
		int tsg = 0; // total sites spanned by stream
		
		if(mls == 1) {
			for(int i = 0 ; i < num_tracts ; i++) {
				if(stream[2][i] == 1) {
					if(stream[0][i] > 0) {nss++;}
					tsg = tsg + stream[1][i];
				}
			}
		}
		
		else {
			if(binary_emission && stream.length != 2 ) {
				System.out.println("Stream structure shape mismatch"); System.exit(0);
			}
			
			for(int i = 0 ; i < num_tracts ; i++) {
				for(int j = 1 ; j < stream.length ; j++ ) {
					int inc  = stream[j][i];
					nss = nss + inc;
					tsg = tsg + inc;
				}
				tsg = tsg + stream[0][i];
			}
		}
		
		return new int[] {nss,tsg} ;
	}
	

	double emToSimpScale(double simp_scale_start, double em_iteration) {
		double out = (simp_scale_start * 1.5) * (Math.pow(0.95, em_iteration));
		return out;
	}
	

	int[][] hap_ind_file_read(String file_name) {
		
		int[][] out = null;
		try {
			// initialize scanner
			Scanner scanner = new Scanner(new File(file_name));
			Scanner line_scan;
			
			// Fill holder structures with the data from the file going line by line
			ArrayList<int[]> holder = new ArrayList<int[]>();
						
			while(scanner.hasNextLine()) {
				// write in one line at a time
				line_scan = new Scanner(scanner.nextLine());  line_scan.useDelimiter("\\s*,\\s*");
				ArrayList<Integer> subgroup = new ArrayList<Integer>();
				
				// process the line
				while(line_scan.hasNextInt()) {
					subgroup.add(line_scan.nextInt());
				}
				
				// see if it looks good
				if (subgroup.size() < 2) {
					System.out.println("Haplotype group to analyze has size less then two. Maybe file-formatting incorrect. Has to be one list of comma-separated integers on each line.");
					System.exit(-1);
				}
				if (Collections.min(subgroup) <= 0) {
					System.out.println ("Minimum index in haplotype group is 0 or less. Start indexing at 1.");
					System.exit(-1);
				}

				// convert to int[]
				int[] subgroup_array = new int[subgroup.size()];
				for(int i = 0 ; i < subgroup_array.length ; i++) {subgroup_array[i] = subgroup.get(i);}
				
				// write to holder
				line_scan.close();
				holder.add(subgroup_array);
				
			}
			
			// copy from holder structures into the out object
			out = new int[holder.size()][];
			for(int i = 0 ; i < out.length; i++) {out[i]= holder.get(i);}
			
			scanner.close();
		} 
		
		catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.print("Could not read dataStream from the file"); System.exit(0);
		}
		return out;		

	}

	int get_max_samples_from_VCFs(String[] vcfs) {
		int out = -1;
		for(int i = 0 ;i < vcfs.length ; i++) {
			VCFFileReader v_read = new VCFFileReader(new File(vcfs[i]),false);
			
			// get num of haps in vcf
			//int v_samp = v_read.iterator().next().getSampleNames().size();
			
			GenotypesContext temp = v_read.iterator().next().getGenotypes();
			int v_samp = temp.getMaxPloidy(0)*temp.getSampleNames().size();

			// record this if its the first, otherwise make sure it agrees with other files'
			if(out > -1) {
				if( v_samp != out) { System.out.println("Error: VCFs have different numbers of samples.");System.exit(0); }
			}
			else { out = v_samp; }
			v_read.close();
		}
		
		return out;
	}
	
	int[][] default_hap_groups(int[] n_samps, int[] n_groups){
		
		ArrayList<int[]> out = new ArrayList<int[]>();
		
		for(int i_ns = 0 ; i_ns < n_samps.length; i_ns++) { // go through each ordered pair of an element n_s and the number of groups associated
			int ns = n_samps[i_ns];
			int groups = n_groups[i_ns];
			
			// for each ordered pair create the appropriate number of groups
			for(int i = 0 ; i < groups ; i++) {
				int[] temp = new int[ns];
				
				// for each group, use non-overlapping n_s samples 
				for(int j = 0; j < ns ; j++) {
					temp[j] = ns * i + j + 1;  
				}
				out.add(temp);
			}		
		}
		
		// list -> array and return
		int[][] haps_groups = new int[out.size()][];
		for(int i = 0 ; i < haps_groups.length;i++) { haps_groups[i]=out.get(i);}
		return haps_groups;
	
	}
	
	
}
