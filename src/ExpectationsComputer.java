import java.io.File;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;


public abstract class ExpectationsComputer {
	
	Helper helper = new Helper();
	
	
	protected RealMatrix tM0;
	protected RealMatrix eM0;
	protected RealVector sM0;
	
	// compute following when we run this method
	
	protected RealMatrix expTr;
	protected RealMatrix expEm;
	protected RealVector expSd;
	
	protected double posteriorLL;
	

	ExpectationsComputer(RealMatrix tp, RealMatrix ep, RealVector sd){
		
	
		tM0 = tp.copy();
		eM0 = ep.copy();
		sM0 = sd.copy();

			
	}
	
	
	//////////////////////////////////
	
	abstract void loadDataFromStream(int[][] observed_data, int mL_size);
	
	abstract void computeExps();
		
	///////////////////////////////////////////////
	
	protected void initializeExp() {
		
		expTr = new Array2DRowRealMatrix(tM0.getRowDimension(), tM0.getColumnDimension());
		expEm = new Array2DRowRealMatrix(eM0.getRowDimension(), eM0.getColumnDimension());
		expSd = new ArrayRealVector(sM0.getDimension());
		posteriorLL = 0;
	
	}
	

	//////////////////////////////////////////////////////////////////////////////////
	// Access Methods ////////
	//////////////////////////

	
	
	public RealMatrix getExpT() {
		return expTr;
	}
	
	
	public RealMatrix getExpE() {
		return expEm;
	}


	public RealVector getExpS() {
		return expSd;
	}

	
	public double get_LL() {
		return posteriorLL;
	}
	
	
}
