import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class Function_B3S implements FunctionUnivariate{ 
	Helper helper = new Helper();
	
	// these are bounds of the real x space the entire function is defined on, before transforming to cardinal knot space
	double x_low_bound; 
	double x_up_bound;	
	int x_divisions;
	
	// used for transforming coordinates, these are lns of the x bounds
	//private double lnx_lb;
	//private double lnx_rb;
	
	double[] b_weights;  		// each is coefficient of b-spline spanning (index-3) to (index+1)
	private double[] b_knots;
	RealMatrix c_3splines;		// each row is coef of x'^3, x'^2, x', 1 for poly spanning (i-3) -> (i-2)
	private RealMatrix b_to_c;  // multiply by (c_j, c_j-1, c_j-2, c_j-3) to obtain coefficients of x'^3, x'^2, x', 1 in spline segment from j to j+1
	
	///////////
	final double bw_scale = 10.; // hard coded scaling, empirically find that 
	//this is rough scale to go between the inference parameters b_weights and pwc-coalescence-rates
	
	
	// for aid in computing the intgrals, integral checkpoints
	RealVector integral_checks;
	
	
	Function_B3S(double x_lb, double x_ub, int x_div, double y_c){
		x_low_bound = x_lb;
		x_up_bound = x_ub;
		x_divisions = x_div;
		
		
		initialize_spline_skeleton(y_c);
		updateCfromB();
		
		generate_integral_checks();
	}
	
	// set up basic elements of the spline structure, the knots and b weights, etc.
	void initialize_spline_skeleton( double constant){
		
		// b_to_c are the hard-coded coefficients to go from b-spline weights to cubic polynomial coefficients
		// hard coded for cardinal b-splines, (uniform integer spacing for knots)
		double[][] temp = new double[][] {{1./6.,-1./2.,1./2.,-1./6.},{0,1./2.,-1.,1./2.},{0,1./2.,0,-1./2.},{0,1./6.,2./3.,1./6.}};
		b_to_c = new Array2DRowRealMatrix(temp);
		
		b_knots = new double[x_divisions + 7];
		for(int i = 0 ; i < b_knots.length ; i++) {b_knots[i] = -3 + i;}
		
		b_weights = new double[x_divisions + 3];
		for(int i = 0 ; i < b_weights.length ; i++) {b_weights[i] = constant;}
		
	}
	
	// use the b-weights to find the cubic spline polynomials, or vis versa
	void updateCfromB() {
	
		// cubic poly for each 1-interval from -3 to n+3
		c_3splines = new Array2DRowRealMatrix(x_divisions + 6,4);
		
		// for each interval, pick the matching weight and trailing 3 weight from b_weights
		// and use to compute the polynomial coefficients
		for(int i = 0 ; i < c_3splines.getRowDimension() ; i++) {
			
			RealVector bs_coef = new ArrayRealVector(4);
			for(int j = 0; j < 4 ; j++) {
				try { bs_coef.setEntry(j, b_weights[i - j]);}
				catch(Exception obe){}
			}
			
			// now bs_coefficients should have the appropriate vector of (c_i, c_i-1, c_i-2, c_i-3)
			c_3splines.setRowVector(i, b_to_c.operate(bs_coef));
		}
		
		
		
		
	}
	
	void updateBfromC() {

		RealVector c2b = MatrixUtils.inverse(b_to_c).getRowVector(0);
		
		for(int i = 0 ; i < b_weights.length ; i++) {
			b_weights[i] = c2b.dotProduct(c_3splines.getRowVector(i));
		}
		
	}
	
	// perform integral up to each knot point (in x-raw space) and store values
	void generate_integral_checks(){
		integral_checks = new ArrayRealVector(x_divisions + 1);
		
		// integral before knot0 is just rectangular area
		integral_checks.setEntry(0, c_3splines.getEntry(3, 3) * x_low_bound);
		
		// now for each knot from 1 till x_divisions, compute integral
		for(int i = 1 ; i < integral_checks.getDimension() ; i++) { // for each knot
			
			int j = i - 1;
			RealVector relevant_poly = c_3splines.getRowVector(j + 3); // poly from i-1 to i
			
			double a = relevant_poly.getEntry(0); // cubic coef
			double b = relevant_poly.getEntry(1); // quad coef
			double c = relevant_poly.getEntry(2); // lin coef
			double d = relevant_poly.getEntry(3); // const coef
			double m1 = x_divisions / (x_up_bound - x_low_bound) ;
			
			RealVector t_poly = new ArrayRealVector(4);
			t_poly.setEntry(0, a/4);
			t_poly.setEntry(1, b/3);
			t_poly.setEntry(2, c/2  );
			t_poly.setEntry(3, d );
			
			double int_chunk;
			int_chunk = t_poly.dotProduct(new ArrayRealVector(4,1.0));
			int_chunk = int_chunk  / m1 ;

			
			integral_checks.setEntry(i, int_chunk + integral_checks.getEntry(j));
		}
		
		
		
	}
	
	//////////////////////////////////////////////////////////////////
	
	//transform raw X coordinate to cardinal x space
	double transform_RX2X(double rawx) {		
		return (x_divisions)*(rawx - x_low_bound)/(x_up_bound - x_low_bound);
	}
	
	double transform_X2RX(double cardx) {
		return x_low_bound + (cardx / x_divisions) * (x_up_bound - x_low_bound);
	}
	
	
	//////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////

	
	@Override
	public double evaluate(double raw_X) {
		
		double x = transform_RX2X(raw_X);
		double out_value;
			
		if(x <= 0) {out_value = c_3splines.getEntry(3, 3);}
		else if(raw_X >= x_up_bound) {out_value = c_3splines.getEntry(x_divisions+3, 3);}
		else {
		
			int bin = (int) Math.floor(x);
			double xp = x - bin;
		
			RealVector poly_coefficients = c_3splines.getRowVector(bin + 3);
			RealVector xpowers = new ArrayRealVector(new double[] {Math.pow(xp, 3), Math.pow(xp, 2), xp, 1});
		
			out_value = poly_coefficients.dotProduct(xpowers);
		}
		
		//return out_value;
		return Math.max(out_value, 0.0000000001);
	}

	@Override
	public double integrate0(double raw_x) {
		if(raw_x < -0.0000000001) {System.out.println("integrand < 0 : Error"); System.exit(0);}
		
		double output;
		
		// convert to appropriate xbar. by convention xbar is 0 < xb <= 1, for use in generate_integral_checks function
		double xbar = transform_RX2X(raw_x);
		int j = ((int) Math.ceil(xbar)) - 1;    // the knot value below x, indicating the j value to use. 
		xbar = xbar - j;    // final xbar val
		
		
		// handle integration differently for each of these cases, first two cases use rectangular area, last uses the poly
		// furst condition defined by x_raw instead of j since j is ill defined for s_raw = 0
		if(raw_x <= x_low_bound || j < 0) {output = c_3splines.getEntry(3, 3) * raw_x;}
		else if(j >= x_divisions) { 
			output = this.integral_checks.getEntry(x_divisions);
			output = output + c_3splines.getEntry(x_divisions+3, 3) * (raw_x - x_up_bound);
		}
		else {
			RealVector relevant_poly = c_3splines.getRowVector(j + 3);
			double a = relevant_poly.getEntry(0); // cubic coef
			double b = relevant_poly.getEntry(1); // quad coef
			double c = relevant_poly.getEntry(2); // lin coef
			double d = relevant_poly.getEntry(3); // const coef
			double m1 = x_divisions / (x_up_bound - x_low_bound) ;
			
			RealVector t_poly = new ArrayRealVector(4);
			t_poly.setEntry(0, a/4);
			t_poly.setEntry(1, b/3);
			t_poly.setEntry(2, c/2  );
			t_poly.setEntry(3, d );
			
			output = t_poly.dotProduct(new ArrayRealVector(new double[] {Math.pow(xbar, 4), Math.pow(xbar, 3),Math.pow(xbar, 2), xbar}));
			output = output  / m1 ;
						
			// add the previous chunks from integral_checkpoints
			output = output + integral_checks.getEntry(j);
		}
		
		
		

		return output;
	}

	@Override
	public double inverseIntegrate0(double in_val) {
		
		System.out.println("LogBSpline class incompatible with inverseIntegrate ability");	
		System.exit(0);
		
		return 0;
	}

	///////////////////////////////
	// Regularizing methods

	@Override
	public double reg_d1l1() {

		double accumulated_integral= 0;
		
		//sum contribution of each polynomial between lb, ub
		for(int p = 3 ; p < x_divisions + 3 ; p++) {
			
			RealVector poly_coefficients = c_3splines.getRowVector(p);
			double a = poly_coefficients.getEntry(0); // cubic coef
			double b = poly_coefficients.getEntry(1); // quad coef
			double c = poly_coefficients.getEntry(2); // lin coef

			double contribution = 0;
		
			// derivative is 3a*x^2  + 2b*x + c
			
			// derivative is zero at
			double zero_1 = (-b - Math.sqrt(b*b-3*a*c))/(3*a);
			double zero_2 = (-b + Math.sqrt(b*b-3*a*c))/(3*a);

			// integrate abs val of derivative from 0 to 1, find points where it crosses zero
			double[] i_xs = new double[3];
			if(zero_1 > 0 && zero_1 < 1) {i_xs[0] = zero_1 ;}
			if(zero_2 > 0 && zero_2 < 1) {i_xs[1] = zero_2 ;}
			i_xs[2] = 1.0;
			
			
			double seg_lb = 0;
			for(int i = 0 ; i < 3; i++) {
				double seg_ub = i_xs[i];
				if(seg_ub > seg_lb) { // only do if seg_ub is not 0, otherwise move to next x interval mark
					// determine if function is pos or neg here
					double sign = (seg_ub + seg_lb) / 2.;
					if(3*a*sign*sign + 2*b*sign + c > 0.) {sign = 1;}
					else {sign = -1;}
					
					RealVector evec = new ArrayRealVector(4,1.);
					for(int ii = 0 ; ii < 4 ; ii++) {
						double tt = Math.pow(seg_ub, 3-ii) - Math.pow(seg_lb, 3-ii);
						evec.setEntry(ii, tt);
					}
					
					
					// evaluate integral over this segment of derivative 
					// value of derivative does not cross zero over this segment (might reach 0 at a bound)
					double seg_integral = 0;
					seg_integral = poly_coefficients.dotProduct(evec);
					
					
					contribution = contribution + sign * seg_integral;
					
				}
			}

			accumulated_integral = accumulated_integral + contribution;
		}
		
		return accumulated_integral;
	}

	@Override
	public double reg_d1l2() {
		
		double accumulated_integral= 0;
		
		//sum contribution of each polynomial between lb, ub
		for(int p = 3 ; p < x_divisions + 3 ; p++) {
			
			RealVector poly_coefficients = c_3splines.getRowVector(p);
			
			double a = poly_coefficients.getEntry(0); // cubic coef
			double b = poly_coefficients.getEntry(1); // quad coef
			double c = poly_coefficients.getEntry(2); // lin coef
			
			double contribution = 0;
			contribution = contribution + 9*(a*a)/5 ;
			contribution = contribution + 3*(a*b);
			contribution = contribution + 4*(b*b)/3;
			contribution = contribution + 2*c*(a+b);
			contribution = contribution + (c*c);


			accumulated_integral = accumulated_integral + contribution;
		}
		
		return accumulated_integral;
	}
	
	@Override
	public double reg_d2l2() {
		
		double accumulated_integral= 0;
		
		//sum contribution of each polynomial between lb, ub
		for(int p = 3 ; p < x_divisions + 3 ; p++) {
			
			RealVector poly_coefficients = c_3splines.getRowVector(p);
			
			double a = poly_coefficients.getEntry(0); // cubic coef
			double b = poly_coefficients.getEntry(1); // quad coef
			
			double contribution = (12 * a * a + 12 * a * b + 4 * b * b);
			
			accumulated_integral = accumulated_integral + contribution;
		}
		
		
		return accumulated_integral;
	}

	@Override
	public double reg_d0l2(double ref) {
		double accumulated_integral= 0;
		
		//sum contribution of each polynomial between lb, ub
		for(int p = 3 ; p < x_divisions + 3 ; p++) {
			
			RealVector poly_coefficients = c_3splines.getRowVector(p);
			
			double a = poly_coefficients.getEntry(0); // cubic coef
			double b = poly_coefficients.getEntry(1); // quad coef
			double c = poly_coefficients.getEntry(2); // lin coef
			double d = poly_coefficients.getEntry(3) - ref; // const coef

			// ^^ we want l2 norm of function - constant (specified by final value in -infinity time)
			// thats why d is defined by subtracting off that term.
			
			double contribution = 0;
			contribution = contribution + (a*a)/7 ;
			contribution = contribution + (a*b)/3;
			contribution = contribution + (2*a*c + b*b)/5;
			contribution = contribution + (a*d + b*c)/2;
			contribution = contribution + (2*b*d + c*c)/3;
			contribution = contribution + (c*d);
			contribution = contribution + (d*d);

			accumulated_integral = accumulated_integral + contribution;
		}
		
		return accumulated_integral;
	}

	///////////////////////////////
	///////////////////////////////

	
	
	
	@Override
	public void updateFromParams(double[] params) {

		if(params.length != b_weights.length) {System.out.println("wrong num of updating params"); System.exit(0);}
		
		
		b_weights = params.clone();
		for(int i = 0 ; i < b_weights.length; i++) {b_weights[i] = b_weights[i] / bw_scale;} // scale to make it better size for nelder meade
		
		
		// fully update the internal spline polys
		updateCfromB();
		generate_integral_checks();
	}

	@Override
	public double[] getFunctionParameters() {
		double[] out = b_weights.clone();
		for(int i = 0 ; i < out.length; i++) {out[i] = out[i] * bw_scale;}
		return out;		
	}

	@Override
	public FunctionUnivariate copy() {
		Function_B3S out = new Function_B3S(x_low_bound, x_up_bound, x_divisions, 1.);
		out.updateFromParams(this.getFunctionParameters());
		
		return out;
	}

	@Override
	public void multiplyScalarSelf(double s) {

		for(int i = 0; i < b_weights.length ; i++) {b_weights[i] = s * b_weights[i];}
		c_3splines = c_3splines.scalarMultiply(s);
		
		integral_checks.mapMultiplyToSelf(s);
	}

	@Override
	public void addScalarSelf(double s) {
		
		c_3splines.setColumnVector(3, c_3splines.getColumnVector(3).mapAdd(s));
		this.updateBfromC();
		
		generate_integral_checks();
		
	}

	@Override
	public void dilateXSelf(double s) {
		x_low_bound = x_low_bound * s;
		x_up_bound = x_up_bound * s;
		
		integral_checks.mapMultiplyToSelf(s);
	}

	@Override
	public double[] getDiscontinuities() {
		// returns points discontinuous in first derivative and higher (but still technically continuous)
		return new double[] {x_low_bound, x_up_bound};
	}

	public void printToFile(String out_file, int density) {
			
			FileWriter scribe;
			try {
				
				scribe = new FileWriter( out_file+ ".csv");
				scribe.append("x, y \n");
				//scribe.append("B3S params:" + Arrays.toString(getFunctionParameters()) + " : function \n");

				
				for(int i = 0 ; i < b_knots.length - 1 ; i++) {
					double dx = (b_knots[i+1] - b_knots[i])/density;
					double x_val = b_knots[i]; // value uniformly spaced in cardinal space
					
					while(x_val < b_knots[i+1]) {
						double raw_x = transform_X2RX(x_val); //transform to raw X domain space
						if(raw_x > 0) {
							scribe.append(raw_x + ", " + evaluate(raw_x) + "\n");
						}
						x_val = x_val + dx;

					}
				}
				
	
				scribe.flush();
				scribe.close();
				
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.out.print("Failed to write function to file");
			}
			
			

	}
	
	@Override
	public void printReciprocalToFile(String out_file, int density) {
		
		FileWriter scribe;
		try {
			
			scribe = new FileWriter( out_file+ ".csv");
			scribe.append("x, y \n");
			//scribe.append("B3S params:" + Arrays.toString(getFunctionParameters()) + " : reciprocal function \n");

			
			for(int i = 0 ; i < b_knots.length - 1 ; i++) {
				double dx = (b_knots[i+1] - b_knots[i])/density;
				double x_val = b_knots[i]; // value uniformly spaced in cardinal space
				
				while(x_val < b_knots[i+1]) {
					double raw_x = transform_X2RX(x_val); //transform to raw X domain space
					if(raw_x > 0) {
						double y_val = 1.0/evaluate(raw_x);
						scribe.append(raw_x + ", " + y_val + "\n");
					}
					x_val = x_val + dx;


				}
			}
			

			scribe.flush();
			scribe.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.print("Failed to write function to file");
		}
		
		
	}
	
	public void printIntegralToFile(String out_file, int density) {
		
		FileWriter scribe;
		try {
			
			scribe = new FileWriter( out_file+ ".csv");
			scribe.append("x, INT(x) \n");

			for(int i = 0 ; i < b_knots.length - 1 ; i++) {
				double dx = (b_knots[i+1] - b_knots[i])/density;
				double x_val = b_knots[i]; // value uniformly spaced in cardinal space
				
				while(x_val < b_knots[i+1]) {
					double raw_x = transform_X2RX(x_val); //transform to raw X domain space
					if(raw_x > 0) {
						scribe.append(raw_x + ", " + integrate0(raw_x) + "\n");
					}
					x_val = x_val + dx;
				}
			}
			

			scribe.flush();
			scribe.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.print("Failed to write function to file");
		}
		
		

}

	//////////////////////////////////////////////////////////////////
	
	
	
	public static void main(String[] args) {
		
		Helper helper = new Helper();
				
		
		Function_B3S fun = new Function_B3S(1,100,8,1);
				
		System.out.println(Arrays.toString(fun.getFunctionParameters()));

		
		RealVector pars = new ArrayRealVector(new double[] {6,6,6,6,7,6,6,5,6,6,6});
		//RealVector pars = new ArrayRealVector(new double[] {1,2,3,3,5,6,7,8,10,10,10});

		pars.mapMultiplyToSelf(100);
		fun.updateFromParams(pars.toArray());

		
		System.out.println(Arrays.toString(fun.getDiscontinuities()));
		System.out.println(Arrays.toString(fun.getFunctionParameters()));

		System.out.println(fun.reg_d1l1());
		
		fun.printToFile("/Users/gautam/Desktop/scratch/fun_spline", 10000);
		
		System.out.println("Done!!!");
	}


}
