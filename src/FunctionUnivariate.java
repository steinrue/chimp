import java.io.FileWriter;
import java.io.IOException;

public interface FunctionUnivariate {
	
	// This is the interface for handling functions, namely demographies. 
	// following methods must be defined
	
	//////////////////////////
	// CORE METHODS
	/////////////////////////
	
	double evaluate(double x);
	
	double integrate0(double x);
	
	double inverseIntegrate0(double in_val); 
	//^^^ only needed for SIMULATING transition, emission matrices
	
	// gives positive quantity that reflects how "normal" function is (0=normal)
	
	double reg_d1l1();
	
	double reg_d1l2();
	
	double reg_d2l2();
	
	double reg_d0l2(double ref);
		
	// go from parametrizing numbers to function, assume domain is fixed
	void updateFromParams(double[] params);
	
	double[] getFunctionParameters();
	
	
	
	// create copy
	FunctionUnivariate copy();
	
	
	// perform basic manipulations 
	void multiplyScalarSelf(double s);
	
	void addScalarSelf(double s);
		
	void dilateXSelf(double s);
	
	
	
	
	// give x points where we expect discontinuities or edges
	double[] getDiscontinuities();
	
	
	
	// produce files showing the function we have. main way that end user will 
	// interact with current state of function
	
	public void printToFile(String out_file, int density);
	
	public void printReciprocalToFile(String out_file, int density);
	

}
