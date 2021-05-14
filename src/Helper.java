import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.InputMismatchException;
import java.util.Random;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.SimpleJSAP;
import com.martiansoftware.jsap.Switch;


import java.util.Scanner;

public class Helper {

	public static void main(String[] args) throws IOException {

		
		
		Helper helper = new Helper();
		
		RealMatrix a = new Array2DRowRealMatrix(2,3);
		int b = 21;
		System.out.println(a);
		
		ChooseComputer chooser = helper.new ChooseComputer(10);
		System.out.println(chooser.get1(10, 10));
		
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
		
		
		
		for(int t = 0 ; t < times.length; t++) {
			
			if( (times[t] - time_counter) <= t_prec) {} // do nothing if time step is too small
			else {
				integrator.integrate(ode, time_counter, residence_probs, times[t], residence_probs);
				time_counter = times[t];
			}
			
			
			// store value of residence_probs at time time_counter for state1 and state2
			for(int s = 0 ; s < ode.getDimension() ; s++) {
				out_res_probs[s][t] = residence_probs[s];
			}	
			
		}
		
		
		if(till_inf) {
			
			Converge converge_condition = new Converge(.0000000001, abs_indices);
			integrator.addEventHandler(converge_condition, 100., .01 , 1000000);
			integrator.integrate(ode, time_counter, residence_probs, 100000, residence_probs); 
			
			for(int s = 0 ; s < ode.getDimension() ; s++) {
				out_res_probs[s][times.length] = residence_probs[s];
			}	
			
		}
		
		

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
		
		// make sure pdf is a proper pdf
		double tot = this.totalEntriesMat(out);
		if(Math.abs(tot - 1.0 ) > 0.0000001) {
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
			else {r.mapDivideToSelf(r_sum); out.setRowVector(i, r);}
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
		
		/// fraction that are segregating site across all loaded data (for Watersons Estimate)
		private double frac_ss;
		
		dataStreamsHandler(String[] file_names, String[] reference_names, String[] ancestral_names,  int[] ls, int[][] sample_index_intervals, int max_tract_size, int mL_size, int n_samples, boolean binary_emission, String temp_prefix){
			// filenames, ref_names, anc_names, ls, sample_intervals, 
			// each are array of equal length, and define the corresponding field for a datastream that we use for analysis
			
			// max tract size and mL_size and n_samples are also important to define structure of data streams
			
			// target_folder is null if we are using memory, otherwise it is name of temporary file prefix
		
			
			if(file_names.length == 0) {System.out.println("No Files Provided"); System.exit(0);}
			
			
			// Elements to compute Watersons estimator for each datastream
			// expectations computer class here is dummy, just used since the vcf's are read through that class
			double num_ss = 0; // number of segregating sites
			double num_sg = 0; // number of non-missing sites processed
			
			// check whether to use memory or files for the streams
			if(temp_prefix == null) { streams_in_mem = true;}
			
			// initialize stream structures
			num_streams = file_names.length;
			if(streams_in_mem) {chromosome_streams = new int[num_streams][][] ;}
			else { 	data_stream_files = new String[num_streams];  }
			
			
			// iterate through each stream and store appropriately. also scan to compute waterson's estimate
			for(int s = 0 ; s < num_streams ; s++) {
				// make sure that the sample_index_intervals are consistent with the num-samples
				int[] interval = sample_index_intervals[s];
				if(interval[1] - interval[0] + 1 != n_samples) {System.out.print("Inconsistency in specified interval and num samples"); System.exit(0);}
				
				VcfReader temp = new VcfReader(file_names[s], reference_names[s], ancestral_names[s], ls[s], interval);
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
				
				// scan stream, count total number seg sites and total sites processed
				int[] ss_frac = getFracSS(mL_size, n_samples, stream, binary_emission);
				num_ss = num_ss + ss_frac[0]; // total number of segregating sites
				num_sg = num_sg + ss_frac[1]; // total number of non missing sites processed
			}
			
			frac_ss = num_ss / num_sg ;

		}
		
		
		
		
		double getFSS() {
			return frac_ss;
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
		
		void deleteStreamFiles() {
			if(!streams_in_mem) {
				for(int f = 0 ; f < num_streams ; f++) {
					delete_csvFromDisc(data_stream_files[f]);
				}
			}

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
	
	int[] getFracSS(int mls, int samps, int[][] stream, boolean binary_emission) {
		// input is metaLocus_size, number of samples, and datastream
		int num_tracts = stream[0].length;		
		int nss = 0; // number of segregating sites of genome
		int tsg = 0; // total sites of genome
		
		if(mls == 1) {
			for(int i = 0 ; i < num_tracts ; i++) {
				if(stream[2][i] == 1) {
					if(stream[0][i] > 0) {nss++;}
					tsg = tsg + stream[1][i];
				}
			}
		}
		
		else {
			if((binary_emission && stream.length != 2) || (!binary_emission && stream.length != samps) ) {
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
	
	protected double estimate_watersons_CR(double f_ss, int n_s, double theta) {
		// find wattersons estimate for pop size, given fraction of segregating sites computed already

		// estimates coalescent rate based on fraction of segregating sites, number of samples, and scaled mutation rate THETA
		// if theta is rescaled, then coalescent rate is appropriately rescaled. 
		
		double harmonic = 0;
		for(int i = 1 ; i < n_s ; i++) {
			harmonic = harmonic + 1.0/i ;
		}
		
		double x1 =  f_ss/harmonic;
		
		return theta / x1;
	}
	
	double emToSimpScale(double simp_scale_start, double em_iteration) {
		double out = (simp_scale_start * 1.5) * (Math.pow(0.95, em_iteration));
		return out;
	}
	
	JSAPResult wrap_arg_parser(String[] args, String help_description) throws IOException, JSAPException {
		SimpleJSAP jsap = new SimpleJSAP("CHMM Scheme",
				help_description,
				new Parameter[] {
						
						
						//////////
						///////// REQUIRED PARAMETERS ///////////////////
						//////////

						
						
						// DEMOGRAPHIC PARAMS
						//////////
						
						// RECOMBINATION RATE: --rec_rate
						new FlaggedOption("recombination_rate", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "rec_rate",
								"Recombination rate, specified per generation per nucleotide (recommend .0000000125 for human data)."),
						
						// MUTATION RATE: --mut_rate
						new FlaggedOption("mutation_rate", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "mut_rate",
								"Mutation rate along a lineage, specified per generation per nucleotide (recommend .0000000125 for human data)."),
						
						// SAMPLES IN SINGLE TREE: --base_n
						new FlaggedOption("base_n_samples", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "base_n",
								"Number of haplotypes considered in each parallel CHMM model (The total number of haplotypes analyzed is <base_n> * <n_groups> under our composite likelihood scheme). We recommend at most choosing --base_n=10 (though for TMRCA the method is still tractable up to 30) and adjusting <n_groups> accordingly to tackle larger sample sizes."),
						
						
						// FILE PARAMS
						//////////


						// INPUT VCF DATA FILE LIST: --vcf_list
						new FlaggedOption("vcf_file_list", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "vcf_list",
								"List of VCF files for each chromosome to perform inference on. Will analyze the first (<base_n> * <n_groups>) haplotypes in each file.").setList(true).setListSeparator(','),
							
						// REFERENCE DATA FILE LIST: --ref_list
						new FlaggedOption("ref_file_list", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ref_list",
								"List of reference files for each chromosome. Indexing and length of list should match [vcf_list] and [anc_list].").setList(true).setListSeparator(','),
						
						// ANCESTRAL DATA FILE LIST: --anc_list
						new FlaggedOption("anc_file_list", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "anc_list",
								"List of ancestral files for each chromosome. Can specify same files as [ref_list] if you don't have these (though accuracy may be affected if reference alleles and ancestral alleles are not the same).").setList(true).setListSeparator(','),
						
						// OUTPUT FILE: --out_file
						new FlaggedOption("output_file", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "out_file",
								"Prefix for output files CHIMP will produce with the inference results. CSV contains ordered pairs of [time, population size], and PARAM file contains, for each epoch, [epoch's lower bound (time), inferred population size]. Times are given in generations from present, and population sizes are in number of diploid individuals."),
						
						
						
						//////////
						///////// OPTIONAL VISIBLE PARAMETERS ///////////////////
						//////////
						

						// NUMBER OF GROUPS OF SUBSAMPLES: --n_groups
						new FlaggedOption("sub_sample_groups", JSAP.INTEGER_PARSER, "1" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "n_groups",
								"Number of subordinate CHMMs (each considering <base_n> haplotypes) analyzed in composite analysis. Total number of haplotypes analyzed is  <base_n> * <n_groups>."),
						
						// LEFT AND RIGHT BOUND OF TIME IN CONSIDERED HISTORY: --t_bounds
						new FlaggedOption("time_bounds", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "t_bounds",
								"This is used to specify epochs for population size history. The epoch boundaries will be spaced exponentially between <min_time> and <max_time>. Defaults to [Ne/50, 5*Ne] where Ne is computed from Waterson's estimator.").setList(true).setListSeparator(','),
						
						// NUMBER OF INDEPENDENT PARAMETERS TO USE IN MODELING PSH: --dof
						new FlaggedOption("degrees_of_freedom", JSAP.INTEGER_PARSER, "18" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "dof",
								"Number of epochs of population size history between [t_bounds] (with 2 additional epochs added at boundaries). One parameter is inferred for each epoch during EM. Note that dof=X will yield X+2 epochs."),
						
						// INITIAL POPULATION SIZE HISTORY GUESS: --psh0
						new FlaggedOption("initial_size_guess", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "psh0",
								"Initial guess for population size. Default is to use Waterson's estimate for N_effective, computed from data (and based on the specified mutation rate). Specifying this option will also use <psh0> instead of Waterson's estimate in computing default for [t_bounds] and the partitioning CHMM states."),
							
						// INITIAL POPULATION SIZE HISTORY DOMAIN: --psh0_xs
						new FlaggedOption("initial_psh_xs", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "psh0_xs",
								"Specify custom list of epoch boundaries (overrides [t_bounds], and <dof>) to use for population size history model.").setList(true).setListSeparator(','),
							
						// INITIAL POPULATION SIZE HISTORY VALS GUESS: --psh0_ys
						new FlaggedOption("initial_psh_ys", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "psh0_ys",
								"Specify population sizes (must also have specified [psh0_xs]) to initialize EM. This option overrides psh0. [psh0_ys] must be a vector with size 1 greater than that of [psh0_xs].").setList(true).setListSeparator(','),
						
						// NUMBER OF STATES IN CHMM: --n_states
						new FlaggedOption("num_chmm_states", JSAP.INTEGER_PARSER, "50" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "n_states",
								"Number of discrete states for CHMM (discrete intervals into which TMRCA or L falls). States are partitioned according to equal probability for each state under a constant population size prior (Waterson's estimate)."),
						
						// METALOCUS IMPLEMENTATION, BEST USED WITH TREELENGTH --metalocus_size
						new FlaggedOption("metalocus_size", JSAP.INTEGER_PARSER, "500" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "metalocus_size",
								"Number of bases grouped into each metalocus. Selecting metalocus_size=1 will activate a locus skipping algorithm to aid computational efficiency."),
						
						// MAX NUMBER OF EM STEPS: --em_cap
						new FlaggedOption("max_em_steps", JSAP.INTEGER_PARSER, "50" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "em_cap",
								"Maximum number of EM steps performed. Inference will end sooner if the convergence criterion is met."),
						
						// MAX NUMBER OF EVALS IN M STEP: --m_evals
						new FlaggedOption("m_step_cap", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "m_evals",
								"Maximum number of transition/emission matrix evaluations allowed during maximization step. Default scales with number of parameters (2 * <dof> + 10) since higher dimensional spaces will need more function evaluations to optimize effectively. M-step will terminate sooner if optimization converges."),
						
						// CRITERION FOR EM CONVERGENCE: --ll_converge
						new FlaggedOption("em_convergence_criterion", JSAP.DOUBLE_PARSER, ".02" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ll_converge",
								"Absolute threshold used to determine convergence of EM. If the posterior log-likelihood improves by less than this amount after any EM step, we consider the method converged."),
						
						// SWITCH TO USE TOTAL TREE LENGTH INSTEAD OF TMRCA: --tree_length
						new Switch("tree_length_switch", 'l', "tree_length", 
								"Use this to specify use of total branch length (L) instead of tree height (TMRCA) as representation of the CHMMs hidden state."),
					
						// PRINT ALL INTERMEDIATE STEPS: --output_steps
						new Switch("print_intermediate_psh", JSAP.NO_SHORTFLAG, "output_steps",
								"This option will produce new output files after each EM step, titled \"EMstep_X\" where X is the EM step number. The files will be placed in a new folder entitled <out_file>."),
						
						// SWITCH FOR PSEDUOHAPLOID EMISSION MODEL --pseudo
						new Switch("pseudo_haploid", JSAP.NO_SHORTFLAG, "pseudo", 
								"Use this switch to specify that the data is pseudohaploid data. CHIMP will assume each allele for each variant in the VCF is randomly selected from either haploid of a diploid individual (independently for each variant/individual). Emission probabilities will be adjusted accordingly."),
					
						
						
						
						
						
						
						
						/////////////
						//////////// INVISIBLE PARAMETERS
						/////////////

						
						///////////// SOMETIMES USEFUL

						// HIDDEN STATE PROBABILITY PARTITIONS: --partitions
						new FlaggedOption("partition_vector", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT,
								JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "partitions",
								"List of marginal probabilities to partition states of the CHMM. Under a constant population prior (from Waterson's estimate) the marginal probability of the TMRCA or L falling into each state will be given by [partitions]. The probabilities are indexed so that the first probability corresponds to the state with the smallest TMRCA (or L) and the last corresponds to the largest. [partitions] should sum to 1. This option overrides n_states.").setList(true).setListSeparator(','),
						
						// INTERNAL POPULATION RESCALE FACTOR: --ps_scale
						new FlaggedOption("population_rescale_factor", JSAP.DOUBLE_PARSER, "1000" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "ps_scale",
								"Internal Population Rescale Factor. Modify this if ODE stepsize error thrown."),
						
						// MAXIMUM TRACT SIZE, REDUCE THIS TO AVOID UNDERFLOW IN FB ALGORITHM
						new FlaggedOption("max_tract_size", JSAP.INTEGER_PARSER, "1000" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "max_tract",
								"Maximum Tract Length of monomorphic tracts during locus skipping-algorithm. The locus-skipping algorithm is only used if <metalocus_size>=1. Reduce <max_tract> to avoid underflow errors during the E step."),
						
						// TREE LENGTH PDE RESOLUTION: --pde_res
						new FlaggedOption("pde_resolution", JSAP.INTEGER_PARSER, "1000,250" ,
								JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "pde_res",
								"Density of PDE solver's grid over the range of relevant times/lengths when using tree-length as hidden state. Higher density increases accuracy but also run-time. First parameter is for the emission probability solver (2D), and second is for transition probability solver (3D). This field is ignored for TMRCA.").setList(true).setListSeparator(','),
						// note can have one or two numbers for the TRANSITION rez (one means evenly spaced, 2 means even for each of the T/2, T/n epochs)
											
						// SIMPLEX SCALE PARAMETER: --simplex_scale
						new FlaggedOption("simplex_init_scale", JSAP.DOUBLE_PARSER, ".1" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "simplex_scale",
								"Size (edge length) of simplex used in M step. Since the search is in log space, 1.0 corresponds to an e-fold increase/decrease of the population size parameter."),
						
						// SWITCH FOR FIXING SIMPLEX SCALE OVER ITERATIONS: --fix_ss
						new Switch("fixed_simplex_scale", JSAP.NO_SHORTFLAG, "fix_ss", 
								"Use this switch to use same size of simplex for all EM steps (default is to slightly shrink the simplex on each successive step)."),
						
						// SCRATCH FOLDER FOR DATASTREAM FILES: --disc_storage 
						new FlaggedOption("temporary_file_prefix", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "disc_storage",
								"Use option to store intermediate data structures temporarily on the disc (default stores them in RAM). Using this option may decrease RAM usage, but may increase run-time since structures will be stored to and read from disc."),
						
						// CHROMOSOMAL LENGTHS FOR VCFS: --chr_l
						new FlaggedOption("chromosomal_length_list", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "chr_l",
								"List of chromosome lengths for each VCF. Will analyze positions [1,length] for each file. Indices of [chr_l] should match those of [vcf_list]. If not provided, the reference file will determine the contiguous chromosomal tract that is analyzed.").setList(true).setListSeparator(','),
						
						
						///////////// NOT EXTENSIVELY TESTED
						
						// REGULARIZATION COEFFICIENT LAMBDA: --reg_lambdas
						new FlaggedOption("regularizing_parameters", JSAP.DOUBLE_PARSER, "0,0,0,0" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "reg_lambdas",
								"Regularizing coefficients used in likelihood optimization. Coefficients of [d0l2, d1l1, d1l2, d2l2] respectively. These coefficients are multiplied by [d0l2,d1l1,d1l2,d2l2] respectively and are subtracted from the objective function. d0l2 corresponds to the squared deviation from a constant population (Waterson’s estimate). d1l1 corresponds to the absolute value difference of the population size change between adjacent epochs, while d1l2 corresponds to the squared difference of the same. d2l2 is an analog for the square of the second derivative (computed for discrete epochs instead), and upweighting it will penalize deviations from linear behavior. The default is no regularization.").setList(true).setListSeparator(','),
						
						// SWITCH FOR REQUESTING SPLINE PSH MODEL: -s --spline
						new Switch("splines", JSAP.NO_SHORTFLAG, "spline", 
								"This option is not fully tested. It specifies that the population size history will be modeled by a cubic spline function (instead of piece-wise constant). The PARAM file will not be generated.  This can be used in conjunction with [t_bounds] and <dof>, which will specify the bounds of the spline and the number of nodes respectively. This should not be used with [psh0_xs] or [psh0_ys]."),
						
						// SWITCH FOR LINEAR TIME MODEL: --lin_t
						new Switch("linear_time", JSAP.NO_SHORTFLAG, "lin_t", 
								"This will distribute the epoch partitions uniformly on a linear scale rather than on a log scale."),
						
						///////////// SHOULD BE FULLY HIDDEN FROM USER
						
						// SIMPLEX TYPE: --s_type
						new FlaggedOption("simplex_type", JSAP.INTEGER_PARSER, "0" , JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "s_type",
								"ID for simplex type that EM_HMM uses at each em step. Default corresponds to an equilateral simplex centered on the current value, adjusted to include current value if it is better than the worst vertex."),
						
						// DIRECTORY WHERE DATA FILES ARE STORED: --data_dir
						new FlaggedOption("data_file_directory", JSAP.STRING_PARSER, JSAP.NO_DEFAULT , JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "data_dir",
								"Directory where data files are stored (include trailing \"/\"). Prepends this to all input files."),
						
						// SWITCH TO USE SFS FOR INFERENCE
						new Switch("SFS_switch", JSAP.NO_SHORTFLAG, "sfs", 
								"Use this switch to specify use of Site Frequency Spectrum (will not use linkage info). Computes marginal from TMRCA transition probability matrix. This is a highly inefficient implementation of SFS based demographic inference."),
						
						// SWITCH FOR BINARY EMISSION MODEL --binary_emission
						new Switch("binary_emission", JSAP.NO_SHORTFLAG, "binary_emission", 
								"Use this switch to specify model with only 2 emission possibilities, segregating or non-segregating."),
					
						
					
						
						//////////////////
						// IN TESTING
						/////////////////
						
						
						
						
						
						
				}
		
		);
	
		
		
		JSAPResult arguments = jsap.parse(args); // obtain the holder which has parsed and knows all the arguments
		
		// Print Usage options if we run into issues, or help is tagged
		if (!arguments.success()) {
			
	        System.err.println();
	        System.err.println("Usage:  ");
	        System.err.println("                "+ jsap.getUsage());
	        System.err.println();
            System.err.println(jsap.getHelp());
	        System.err.println();

	        
			if(!arguments.success()) {
				for (java.util.Iterator errs = arguments.getErrorMessageIterator();
	                    errs.hasNext();) {
	                System.err.println("Error: " + errs.next());
	            }
		        System.err.println();
			}
			System.exit(1);

	    }
		
		
		//modify default values contingently
		
		
		
		return arguments;
	}
	
	
	
}
