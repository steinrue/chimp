
import java.util.ArrayList;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.linear.*;

public class Function_Spline {
	Helper helper = new Helper();
	
	int definingArraySize; // num discrete vals in the array used to fit spline. later 0nodes mean we may store fewer points than this.
	
	double[] x_vals; //discrete x vals 
	double[] y_vals; // discrete y vals
	int num_vals;
	
	double domain_lb;
	double domain_ub;
	
	PolynomialSplineFunction spline_function;
	
	boolean linear_approx; // use switch for doing linear interp instead of spline interp

	
	Function_Spline(RealVector xs, RealVector ys, double node0x, boolean lin_interp) {
		this(xs.toArray(), ys.toArray(), node0x, lin_interp);
	}
	
	Function_Spline(double[] xs, double[] ys, double node0x, boolean lin_interp){
		
		linear_approx = lin_interp;
		
		// check to make sure xs and ys are valid
		if(xs.length != ys.length) {
			System.out.println("Invalid x values for spline_function");
			System.exit(0);
		}
		
		for(int i = 0 ; i < xs.length-1 ; i++) {
			if(xs[i] >= xs[i+1]) {
				System.out.println("Invalid x values for spline_function");
				System.exit(0);
			}
		}
		
		num_vals = xs.length;
		definingArraySize = num_vals;
		domain_lb = xs[0];
		domain_ub = xs[num_vals-1];
		
		
		// find where node0x fits in, and adjust accordingly
		double x_prec = 1E-12;
		if(node0x > domain_lb + x_prec && node0x < domain_ub) {
			int node_index = helper.placeInPartition(node0x, xs);
			if( node0x - xs[node_index-1] < 1E-12) {node_index--;} // if node0x falls on grid need to adjust by one.
			
			num_vals = node_index + 1;
			x_vals = new double[num_vals];
			y_vals = new double[num_vals];
			domain_ub = node0x;
			
			for(int i = 0 ; i < num_vals - 1 ; i++) {
				x_vals[i] = xs[i];
				y_vals[i] = ys[i];
			}
			x_vals[node_index] = node0x;
			y_vals[node_index] = 0;
		}
		
		else { // node0 is irrelevant
			x_vals = xs.clone();
			y_vals = ys.clone();
			
		}
		
		
		if(num_vals == 2 || linear_approx) {
			linear_approx = true;
		}
		

		else {
			SplineInterpolator sp_int = new SplineInterpolator();
			spline_function = sp_int.interpolate(x_vals,y_vals);
		}
		
	
		
	}

	public double evaluate(double x) {
		
		if(linear_approx) {
			
			for(int i = 1 ; i < num_vals ; i++) {
				if(x_vals[i] >= x) {
					return (x - x_vals[i-1])*(y_vals[i] - y_vals[i-1])/(x_vals[i] - x_vals[i-1])+ y_vals[i-1] ;
				}
			}
			
			return Double.NaN;
		}
		
		else {
			return spline_function.value(x);
		}
	}
	
	public double evaluate_0( double x) {
		if(x < domain_lb) {return 0;}
		
		if(x > domain_ub) {return y_vals[num_vals - 1];}
		
		return evaluate(x);
	}
	
	
	
	// this returns an array of values that are the integral from lb till that point in x_vals,
	// also allows specification of discontinuities and 0-nodes
	public RealVector integrate(double iVal, double[] dcs) {
		
		// first process the relevant discontinuities and nodes.
		ArrayList<Double> discontinuities = new ArrayList<Double> (dcs.length);
		ArrayList<Integer> disc_bins = new ArrayList<Integer> (dcs.length);
		
		
		for(int i = 0 ; i < dcs.length ; i++) {
			double dis = dcs[i];
			if(dis >= domain_lb && dis <= domain_ub) {
				discontinuities.add(dis);      
				disc_bins.add(helper.placeInPartition(dis, x_vals));
			}
		}
		
		
		
		// now we have list of relevant discontinuities and nodes and the bins where they have effect
		
		// note if discontinuity is on a x_domain value, we assume it is associated with the interval to the left of the value
		// this convention is inherited from how function_domain values are operated, and assumption we integrate forward in time
		
		// begin integration process
		double integral_counter = iVal;
		
		RealVector integral_vals = new ArrayRealVector(definingArraySize);
		integral_vals.setEntry(0, integral_counter);
		
		PolynomialFunction[] polys = null;
		if(!linear_approx) { polys = this.spline_function.getPolynomials(); }
		
		for(int t= 1; t < definingArraySize; t++) {
			double y = 0; 
			
			if(t > num_vals - 1) {} // if we are looking at index above where we have hit 0, y remains 0 as there is no more contribution to integral
			
			
			else if(disc_bins.contains(t)) {
				// two contributions to y, yL and yR from either side till nearest discrete domain point
				double disco = discontinuities.get(disc_bins.indexOf(t));

				
				double yL;
				double deltaL = disco - x_vals[t-1];
				if(t == 1) { yL = deltaL * y_vals[t-1]; }
				else {
					double x1 = x_vals[t-2];
					double y1 = y_vals[t-2];
					double x2 = x_vals[t-1];
					double y2 = y_vals[t-1];
					double dvL = y2 + deltaL * ((y2 - y1)/(x2 - x1)); 
							
					yL = (dvL + y2) * deltaL / 2.;
				}
				
				double yR;
				double deltaR = x_vals[t] - disco ;
				if(t == num_vals - 1) { yR = deltaR * y_vals[t]; }
				else {
					double x1 = x_vals[t];
					double y1 = y_vals[t];
					double x2 = x_vals[t+1];
					double y2 = y_vals[t+1];
					double dvR = y1 - deltaR * ((y2 - y1)/(x2 - x1));
					
					yR = (dvR + y1) * deltaR / 2.;
				}
				
				y = yL + yR;
				
			}
			
			else if(linear_approx || disc_bins.contains(t-1) || disc_bins.contains(t-2) || disc_bins.contains(t+1) || disc_bins.contains(t+2) || disc_bins.contains(t+3) || disc_bins.contains(t-3)) {
				double delta = x_vals[t] - x_vals[t-1];
				y = delta * ( y_vals[t] + y_vals[t-1] )/2.;
			}
			
			else {
				double delta = x_vals[t] - x_vals[t-1];

				double[] polyCs = polys[t-1].getCoefficients(); // these coeff describe f(x) where x = 0 corresp to times[c + t - 1]
				int degree = polys[t-1].degree();
				
				//
				if(degree >= 3) { y = delta * polyCs[3] / 4. ; }
				if(degree >= 2) { y = delta * ( y + polyCs[2] / 3.); }
				if(degree >= 1) { y = delta * ( y + polyCs[1] / 2.); }
				if(degree >= 0) { y = delta * ( y + polyCs[0]); }
				// now y is integral over this interval
			}
			
			
			integral_counter = integral_counter + y;
			integral_vals.setEntry(t, integral_counter);
			
			
			
		}
		
		
		
		
		
		return integral_vals;
	}

	
	public static void main(String[] args) {
		System.out.println(4E-2 + 3E-4);		
		
		/*
		double[] xs = new double[] {0,1,2,3,4,5,6,7,8,9,10};
		double[] ys = new double[] {22,13,6,1,-2,-3,-2,1,6,13,22};

		spline_function test = new spline_function(xs, ys);
		System.out.println(test.evaluate(10));
		*/
		double[] xs = new double[] {0,1,2,3};
		double[] ys = new double[] {0,2,5,7};

		Function_Spline test = new Function_Spline(xs, ys, -1, false);
		//System.out.println(test.evaluate(4.5));
		System.out.println(test.integrate(0, new double[] {3.1}));
		
	}

}
