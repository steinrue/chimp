import org.apache.commons.math3.analysis.interpolation.PiecewiseBicubicSplineInterpolatingFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class Function_BiSpline {
	
	double[] x_vals; //discrete x vals 
	double[] y_vals; // discrete y vals
	
	double x_lb;
	double x_ub;
	double y_lb;
	double y_ub;
	
	int x_len;
	int y_len;
	
	//double margin = .0000000001;
	
	double[][] f_vals;
	
	PiecewiseBicubicSplineInterpolatingFunction func; 
	
	boolean linear_approx;
	
	
	
	///////////////////////////////
	// CONSTRUCTORS //
	/////////////////////////////////
	Function_BiSpline(RealVector xs, RealVector ys, RealMatrix fs, boolean lin_interp){
		linear_approx = lin_interp;
		
		//this(xs.toArray(), ys.toArray(), fs.getData());
		int xlen = xs.getDimension();
		int ylen = ys.getDimension();
		
		int vxlen = fs.getRowDimension();
		int vylen = fs.getColumnDimension();
		
		
		// check to make sure xs are in ascending order.
		if(xlen <= 1 || ylen <=1 ) {
			System.out.println("need more than 1 discrete point along dimension");
			System.exit(0);
		}
					
		if( ! (xlen == vxlen && ylen == vylen)  ) {
			System.out.println("invalid coordinates to define grid");
			System.exit(0);
		}
		for(int i = 0 ; i < xlen-1 ; i++) {
			if(xs.getEntry(i) >= xs.getEntry(i+1)) {
				System.out.println("Invalid x values for spline_function");
				System.exit(0);
			}
		}
		for(int i = 0 ; i < ylen-1 ; i++) {
			if(ys.getEntry(i) >= ys.getEntry(i+1)) {
				System.out.println("Invalid y values for spline_function");
				System.exit(0);
			}
		}
						
		x_vals = xs.toArray();
		y_vals = ys.toArray();
		x_len = xlen;
		y_len = ylen;
						
		x_lb = x_vals[0];
		x_ub = x_vals[x_len - 1];
		y_lb = y_vals[0];
		y_ub = y_vals[y_len - 1];
						
		f_vals = fs.getData();
		
				
		if(xlen < 5 || ylen < 5 || linear_approx) {linear_approx = true;}
		else { func =  new PiecewiseBicubicSplineInterpolatingFunction(x_vals, y_vals, f_vals); }
						
		
	}

	Function_BiSpline(double[] xs, double[] ys, double[][] fs, boolean lin_interp){
			linear_approx = lin_interp;

		// check to make sure xs are in ascending order.
			if(xs.length <= 1 || ys.length <=1 ) {
				System.out.println("need more than 1 discrete point along dimension");
				System.exit(0);
			}
			
			if( ! (xs.length == fs.length && ys.length == fs[0].length)  ) {
				System.out.println("invalid coordinates to define grid");
				System.exit(0);
			}
			for(int i = 0 ; i < xs.length-1 ; i++) {
				if(xs[i] >= xs[i+1]) {
					System.out.println("Invalid x values for spline_function");
					System.exit(0);
				}
			}
			for(int i = 0 ; i < ys.length-1 ; i++) {
				if(ys[i] >= ys[i+1]) {
					System.out.println("Invalid y values for spline_function");
					System.exit(0);
				}
			}
				
			x_vals = xs.clone();
			y_vals = ys.clone();
			x_len = x_vals.length;
			y_len = y_vals.length;
				
			x_lb = x_vals[0];
			x_ub = x_vals[x_len - 1];
			y_lb = y_vals[0];
			y_ub = y_vals[y_len - 1];
				
			f_vals = fs.clone();
			for(int i = 0 ; i < fs.length ; i++) {
				f_vals[i] = fs[i].clone();
			}
				
				
			if(x_len < 5 || y_len < 5) {linear_approx = true;}
			else { func =  new PiecewiseBicubicSplineInterpolatingFunction(x_vals, y_vals, f_vals); }
				
		
	}
	
	
	////////////////////////////////////////////
	// OBTAIN VALUE, GIVEN THE BOUNDARY STRUCTURE
	///////////////////////////////////////////////
	
	private double value(double x, double y) {
		
		if(linear_approx) {
	
			// find indices of which farmpatch the coordinate is in, and the bounding coordinates
			int xbox = timeToBin(x, x_vals);  
			int ybox = timeToBin(y, y_vals);
			
			double x0 = x_vals[xbox - 1];
			double x1 = x_vals[xbox];
			double y0 = y_vals[ybox - 1];
			double y1 = y_vals[ybox];
			
			// scale coordinates to be on unit square
			double x_scale = (x - x0)/(x1 - x0);
			double y_scale = (y - y0)/(y1 - y0);
			
			// find the 4 bounding vals of the farmpatch to do a linear interpolation.
			double v_x0y0 = f_vals[xbox-1][ybox-1];
			double v_x0y1 = f_vals[xbox-1][ybox];
			double v_x1y0 = f_vals[xbox][ybox-1];
			double v_x1y1 = f_vals[xbox][ybox];

			// bilinear interpolation pulled from wikipedias page on this for unit square
			double vxy = v_x0y0 * (1 - x_scale) * (1 - y_scale);
			vxy = vxy + v_x0y1 * (1 - x_scale) * y_scale ;
			vxy = vxy + v_x1y0 * x_scale * (1 - y_scale) ;
			vxy = vxy + v_x1y1 * x_scale * y_scale;
			
			return vxy;
			
			
		}
		
	
		else {
			return func.value(x, y);
		}
	}
	
	double value0(double x, double y) {
		
		// this is set up to return values that are structured the way we see in pde solver
		
		if(x < x_lb || y < y_lb) {
			return 0;
		}
		
		//if(x < x_lb) {x = x_lb;}
		if(x > x_ub) {x = x_ub;}
		//if(y < y_lb) {y = y_lb;}
		if(y > y_ub) {y = y_ub;}
		
		return value(x,y);
		
	}

	private int timeToBin(double t, double[] dts) {
		
		for(int i = 1 ; i < dts.length; i++) {
			if(t <= dts[i]) {
				return i;
			}
		}
		return -1; // this happens if t is below lowest dts val or above highest
	}
	
	
	
	//////////////////////////////////////////////////
	
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		//RealVector xs = new ArrayRealVector( new double[] {2,3,4,5,6,7} );
		//RealVector ys = xs.copy();
		//RealMatrix vals = new Array2DRowRealMatrix( new double[][] {{0,2,10,2,12,4},{5,10,4,1,7,4} , {1,6,3,3,2,1},{1,6,3,3,2,1},{5,10,4,1,7,4},{0,2,10,2,12,4}}  );
		
		RealVector xs = new ArrayRealVector( new double[] {2,3,4,5,} );
		RealVector ys = xs.copy();
		RealMatrix vals = new Array2DRowRealMatrix( new double[][] {{0,2,10,2},{5,10,4,1} , {1,6,3,3},{1,6,3,3}}  );
				
		Function_BiSpline test = new Function_BiSpline(xs.toArray(),ys.toArray(),vals.getData(), false);
		
		System.out.println(test.value0(2.5,2.5));
		System.out.println(test.timeToBin(.5,  new double[] {1,2,3}));
		
	}

}
