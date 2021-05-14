import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class Function_PWC implements FunctionUnivariate {
	Helper helper = new Helper();

	double[] domain;
	double[] function; 

	int discrete_fbins;  // number of discrete finite bins. do not include the last bin that goes to infinity
	boolean lx_scale = false;
	
	Function_PWC(double[] dom , double[] vals, boolean logx){
		boolean valid = true;
		
		for(int i = 0 ; i < dom.length - 1 ;  i++) {
			if(dom[i] > dom[i+1]) {valid = false;}
		}
		
		if(dom.length != vals.length - 1) {valid = false;}
		
		
		if(valid) {
			function = vals.clone();
			domain = dom.clone();
			discrete_fbins = domain.length;
		}
		else {
			System.out.println("Invalid attempt to initialize Function_PWExp structure");
			System.exit(0);
		}
		
		lx_scale = logx;
		
	}
	
	Function_PWC(double[] dom , double constant_val, boolean logx){
		domain = dom.clone();
		function = new double[dom.length+1];
		for(int i = 0 ;i < function.length;i++) {function[i] = constant_val;}
		
		discrete_fbins = domain.length;
		lx_scale = logx;

	}

	Function_PWC(double x_lb, double x_rb, int num_div, double const_val, boolean logx){
		if(num_div < 1) {System.out.println("Need 0 divs for this constructor"); System.exit(0);}
		
		// make the domain divisions log uniform if specified
		if(logx) {	domain = helper.generate_log_discretization(x_lb, x_rb, num_div + 2);  lx_scale = true;}
		else { domain = helper.generate_uniform_discretization(x_lb, x_rb, num_div + 2);		}
		
		function = new double[domain.length+1];
		for(int i = 0 ;i < function.length;i++) {function[i] = const_val;}
		
		discrete_fbins = domain.length;
	}
	
	

	// CORE METHODS FOR FunctionUnivariate INTERFACE
	
	@Override
	public double evaluate(double x) {

		for(int i = 0 ; i < discrete_fbins; i++) {
			if(domain[i] >= x) {
				return function[i];
			}
		}
		
		return function[discrete_fbins];
		
	}
	
	@Override
	public double integrate0(double x) {

		double sum = 0;
		
		int i = 0;
		double left = 0;   // lower bound of bin in question
		
		while(i < discrete_fbins && domain[i] <= x ) {
			sum = sum + (domain[i] - left) * function[i];
			left = domain[i] ;
			i++;
		}
		
		sum = sum + (x - left) * function[i] ;  // add the partial contribution of the incomplete bin
		
		return sum;
	}

	@Override
	public double inverseIntegrate0(double in_val) {
		double sum = in_val ; 
		int i = 0 ; 
		double gap;
		boolean cont = true;
		
		while(i < discrete_fbins && cont) {
			if(i == 0) {   gap  = domain[0];	}
			else{ gap = domain[i] - domain[i-1];}
			
			double newsum = sum - gap * function[i] ; 
			if(newsum < 0) { cont = false;}
			else {
				sum = newsum ;
				i++;
			}	
		}
		
		double lb = 0;
		if(i > 0) {lb = domain[i-1];}
		return lb + sum /function[i]  ;
	}

	//////////////////////////
	// REGULARIZING QUANTITIES

	@Override
	public double reg_d1l1() {
		
		if(discrete_fbins < 2) {return 0;}

		// actually compute the finite difference scheme regularity 
		
		double accumulated_integral = 0;
		double sb_width = 1;
		for(int bin = 1 ; bin < discrete_fbins ; bin++) {
			
			sb_width = domain[bin] - domain[bin - 1] ;
			if(lx_scale) {sb_width = Math.log(domain[bin]) - Math.log(domain[bin - 1]);}
			
			double der = function[bin] - function[bin-1];
			accumulated_integral = accumulated_integral + Math.abs(der);
		}
		
		// now do it for the last finite difference, but for bin width use the second bounded bin ie current sb_width val
		double der = function[discrete_fbins] - function[discrete_fbins-1];
		accumulated_integral = accumulated_integral + Math.abs(der) ;
				
		return accumulated_integral;
	}

	@Override
	public double reg_d1l2() {
		
		if(discrete_fbins < 2) {return 0;}

		// actually compute the finite difference scheme regularity 
		
		double accumulated_integral = 0;
		double sb_width = 1;
		for(int bin = 1 ; bin < discrete_fbins ; bin++) {
			
			sb_width = domain[bin] - domain[bin - 1] ;
			if(lx_scale) {sb_width = Math.log(domain[bin]) - Math.log(domain[bin - 1]);}
			
			double der = function[bin] - function[bin-1];
			
			accumulated_integral = accumulated_integral + Math.pow(der,2) / sb_width;
		}
		
		// now do it for the last finite difference, but for bin width use the second bounded bin ie current sb_width val
		double der = function[discrete_fbins] - function[discrete_fbins-1];
		accumulated_integral = accumulated_integral + Math.pow(der,2) / sb_width;
				
		return accumulated_integral;
	}

	@Override
	public double reg_d2l2() {
		if(discrete_fbins < 2) {return 0;}

		// actually compute the finite difference scheme regularity 
	
		double accumulated_integral = 0;
		for(int bin = 1 ; bin < discrete_fbins ; bin++) {
			double sb_width = domain[bin] - domain[bin - 1] ;
			if(lx_scale) {sb_width = Math.log(domain[bin]) - Math.log(domain[bin - 1]);}
			
			double sec_der = function[bin + 1] - 2 * function[bin] + function[bin-1];
			
			double contribution = Math.pow(sec_der,2) / Math.pow(sb_width,3);
			
			accumulated_integral = accumulated_integral + contribution;
		}
				
		return accumulated_integral;
	}

	@Override
	public double reg_d0l2(double ref) {
		
		if(discrete_fbins < 1) {return 0;}

		// actually compute the finite difference scheme regularity 
		
		
		double accumulated_integral = 0;
		for(int bin = 1 ; bin < discrete_fbins ; bin++) {
			double sb_width; 
			if(lx_scale) { sb_width = Math.log(domain[bin]) - Math.log(domain[bin - 1]); }
			else { sb_width = domain[bin] - domain[bin - 1]; }			
			
			double diff = function[bin] - ref;
			double contribution = Math.pow(diff,2) * sb_width;
			
			accumulated_integral = accumulated_integral + contribution;
		}
				
		return accumulated_integral;
	}


	//////////////////////////
	//////////////////////////

	
	
	
	public void printToFile(String out_file, int density) {
		
		FileWriter scribe;
		try {
			
			scribe = new FileWriter( out_file+ ".csv");
			scribe.append("x, y \n");
			//scribe.append("PWC params:" + Arrays.toString(getFunctionParameters()) + " : function \n");
			
			double dom_lb = 0;
			double dx = 0;
			
			for(int i = 0 ; i < domain.length + 1 ; i++) {
				
				if(i!= domain.length) {	dx = (domain[i] - dom_lb)/density; }
				else {}
				
				for(int j = 0 ; j < density ; j++) {
					double x_val = dom_lb + j * dx ;
					scribe.append(x_val + ", " + evaluate(x_val) + "\n");	
				}
				
				
				if(i!= domain.length) {dom_lb = domain[i];}
				
			}
			

			scribe.flush();
			scribe.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.print("Failed to write function to file");
		}
		
		

}

	public void printReciprocalToFile(String out_file, int density) {
		printParamsToFile(out_file);
		
		FileWriter scribe;
		try {
			
			scribe = new FileWriter( out_file+ ".csv");
			scribe.append("x, y \n");

			
			double dom_lb = 0;
			double dx = 0;
			
			for(int i = 0 ; i < domain.length + 1 ; i++) {
				
				if(i!= domain.length) {	dx = (domain[i] - dom_lb)/density; }
				else {}
				
				for(int j = 0 ; j < density ; j++) {
					double x_val = dom_lb + j * dx ;
					double y_val = 1.0/evaluate(x_val);
					scribe.append(x_val + ", " + y_val + "\n");	
				}
				
				
				if(i!= domain.length) {dom_lb = domain[i];}
				
			}
			

			scribe.flush();
			scribe.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.print("Failed to write function to file");
		}
		
		
	}

	public void printParamsToFile(String out_file) {

		FileWriter scribe;
		try {
			
			scribe = new FileWriter( out_file+ ".param");
			scribe.append("epoch_left_bounds, demographic_parameter \n");

			
			scribe.append("0.0" + ", " + (1.0/function[0]) + "\n");
			for(int i = 1 ; i < function.length ; i++) {
				scribe.append(domain[i-1] + ", " + (1.0/function[i]) + "\n");
			}
			
			scribe.flush();
			scribe.close();
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.out.print("Failed to write function to file");
		}
		
		
	}
	
	
	// OBTAIN AND UPDATE RELEVANT PARAMETERS OF FUNCTION
	
	// the parameters that come in and go out could be negative, so we use the exp to 
	// convert between the log parameter space and the values written into function
	@Override
	public void updateFromParams(double[] params) {
		if(params.length != function.length) {
			System.out.println("Updating PWExp with incompatible params"); System.exit(0);
		}
		
		/*
		 * use this for non-log space optimization
		double[] temp = params.clone();
		for(int i = 0 ; i < temp.length ; i++) {if(temp[i] <= 0) {temp[i] = 0.00000001;}}
		function = params.clone(); 
		*/
		
		function = helper.ebeExp(params);
	}

	@Override
	public double[] getFunctionParameters() {
		
		// USE this line for non log space optimization
		// return function.clone(); 
		
		return helper.ebeLog(function);
	}

	public FunctionUnivariate copy() {
		Function_PWC out = new Function_PWC(domain, function, lx_scale);
		return out;
	}

	@Override
	public double[] getDiscontinuities() {
		return domain.clone();
	}





	@Override
	public void multiplyScalarSelf(double s) {
		
		for(int i = 0 ; i < function.length; i++) {
			function[i] = function[i] * s;
		}
		
	}



	@Override
	public void addScalarSelf(double s) {
		
		for(int i = 0 ; i < function.length; i++) {
			function[i] = function[i] + s;
		}
		
	}
	
	public void reciprocalSelf() {
		for(int i = 0 ; i < function.length; i++) {
			function[i] = 1.0/function[i];
		}
	}

	public void dilateXSelf(double s) {
		
		for(int i = 0 ; i < domain.length; i++) {
			domain[i] = s * domain[i];
		}
		
	}

	

	public static void main(String[] args) {

		Helper helper = new Helper();
		
		FunctionUnivariate fun = new Function_PWC(4,40,8,10000,false);
		//fun.updateFromParams(helper.ebeLog(new double[] {100,200,500,200,200,200,700,200,200,200}));
		//fun.updateFromParams(helper.ebeLog(new double[] {100,200,300,300,500,600,700,800,900,1000}));
		
		System.out.println(Arrays.toString(fun.getDiscontinuities()));
		System.out.println(Arrays.toString(helper.ebeExp(fun.getFunctionParameters())));

		System.out.println(fun.reg_d2l2());

		fun.printReciprocalToFile("/Users/gautam/Desktop/scratch/fun_pwc", 10000);
		System.out.println("done");
		
		
	}






	
	
}
