import java.io.*;
import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.PointValuePair;

public class MHutt_NMSimplex {

	//////////////////////////////////////////////////
	
	public static void main(String[] args) throws IOException
	{
		//double start[] = {-1.7,1.5};
		//int dim = 2;

		double start[] = {5.,-3.0,-2.3};
		int dim = 3;

		double eps = 1.0e-4;
		double scale = 1;

		double t_Scale = 3.0;
		final class test_func implements MultivariateFunction
		{	
			@Override
			public double value(double[] x){

				// TODO Auto-generated method stub
				//return t_Scale*(100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])+(1.0-x[0])*(1.0-x[0]));
				return -Math.exp(-((x[0]-t_Scale)*(x[0]- t_Scale) + x[1]*x[1] + x[2]*x[2]) / 1.);
			}
		}		
		
		MHutt_NMSimplex simplex = new MHutt_NMSimplex(1000);
		
		PointValuePair pvp = simplex.minimize(400, start, eps, scale, new test_func(),0);
		
		System.out.println(Arrays.toString(pvp.getPoint()) + ": " + pvp.getValue());
	}
	
	//////////////////////////////////////////////////
	static int MAX_IT;      /* maximum number of iterations */

	
	
	
	public MHutt_NMSimplex(int max_iter) {
		// currently we don't use simplex_type info
		MAX_IT = max_iter;
	}
	
	
	public PointValuePair minimize(int max_m_evals, double i_guess[], double EPSILON, double scale, MultivariateFunction q, int simplex_type) {
		
		int n = i_guess.length; // dimensionality of space
		final double MAX_EV = max_m_evals; //  max number of evals (may go few over)    // can simplify by using combining this with k (k is the tracker for this that the original code already uses)

		
		double ALPHA  = 1.0;       /* reflection coefficient */
		double BETA   = 0.5;       /* contraction coefficient */
		double GAMMA  = 2.0;       /* expansion coefficient */
		double SIGMA  = 0.5;       /* shrink coefficient */

		// if using adaptive transformations to do better in high-d space as per scipy-nelder-mead
		//GAMMA = 1. + 2./n ;
		//BETA = 0.75 - 1./(2. * n) ;
		//SIGMA = 1. - 1./n ;
		

		
		Constraint c = new Constraint(); // artificial constraint on the coordinates
	    MultivariateFunction objf = q; // objective function
	    
	    double[] start = i_guess.clone(); // start position
	    // EPSILON is the convergence precision
	    // scale is the scale of initial simplex
	    
	    
	    double v[][] = new double[n+1][n];
	    double f[]   = new double[n+1];
	    double vr[]  = new double[n];
	    double ve[]  = new double[n];
	    double vc[]  = new double[n];
	    double vm[]  = new double[n];
	    double fr;      /* value of function at reflection point */
	  	double fe;      /* value of function at expansion point */
	  	double fc;      /* value of function at contraction point */
	    double pn,qn;   /* values used to create initial simplex */
	    double fsum,favg,s,cent;
	    int vs;         /* vertex with smallest value */
	  	int vh;         /* vertex with next smallest value */
	  	int vg;         /* vertex with largest value */
	    int i,j,m,row = 0;
	    int k;   	      /* track the number of function evaluations */
	    int itr;	      /* track the number of iterations */

	    
	    
	    
	    
	    ///////////////////////////////////////////////
	    // INITIALIZE SIMPLEX
	    ///////////////////////////////////
	    // create the initial simplex 
	    
	    if(simplex_type == 0) { // NELDER MEAD SCIPY-OPTIMIZE INITIALIZATION
		    for (i=0;i<n;i++) {
		  		v[0][i] = start[i];
		  	}
		    // for next points
		    for(i = 0 ; i < n; i++) {
		    	v[i+1] = start.clone();
		    	if(v[i+1][i] != 0) { v[i+1][i] = v[i+1][i] * (1.05); }
		    	else { v[i+1][i] = 0.00025 ; }
		    }
	    }
	    
	    else if (simplex_type ==1) { // MIKE-HUTT DEFAULT
	    	
	    	// assume one of the vertices is 0,0 
		    
		    //System.out.format("Starting from : %f\n",start[0]);

		    pn = scale*(Math.sqrt(n+1)-1+n)/(n*Math.sqrt(2));
		    qn = scale*(Math.sqrt(n+1)-1)/(n*Math.sqrt(2));
		    
		    for (i=0;i<n;i++) {
		  		v[0][i] = start[i];
		  	}
		    
		    //System.out.format("pn:%f qn:%f\n",pn,qn);

		    for (i=1;i<=n;i++) {
		      for (j=0;j<n;j++) {
		        if (i-1 == j) {
		          v[i][j] = pn + start[j];
		        }
		        else {
		          v[i][j] = qn + start[j];
		        }
		      }
		      c.getConstrainedValues(v[i],n);
		    }
	    }
	    
	    
	    else {	System.out.println("INVALID SIMPLEX ID."); System.exit(0); } // invalid key
	    
	    //////////////////////////////////////////////////////
	    //////////////////////////////////////////////////////
	    //////////////////////////////////////////////////////
	    
	    
	    
	    
	    /* find the initial function values */
	    for (j=0;j<=n;j++) {
	      f[j] = objf.value(v[j]);
	    }

	    k = n+1;
	    
	    /* print out the initial values */
	    /*
	  	System.out.println("Initial Values");
	  	for (j=0;j<=n;j++) {
	  	  for (i=0;i<n;i++) {
	  		  System.out.format("%f %f\n",v[j][i],f[j]);
			  }
	  	}
	  	*/
	    
	  	/* begin the main loop of the minimization */
	    boolean maxed_eval = false;
	    int num_evals = 0;
	  	for (itr=1;itr<=MAX_IT && !maxed_eval;itr++) {     
	  		/* find the index of the largest value */
	  		vg=0;
	  		for (j=0;j<=n;j++) {
	  			if (f[j] > f[vg]) {
	  				vg = j;
	  			}
	  		}

	  		/* find the index of the smallest value */
	  		vs=0;
	  		for (j=0;j<=n;j++) {
	  			if (f[j] < f[vs]) {
	  				vs = j;
	  			}
	  		}

	  		/* find the index of the second largest value */
	  		vh=vs;
	  		for (j=0;j<=n;j++) {
	  			if (f[j] > f[vh] && f[j] < f[vg]) {
	  				vh = j;
	  			}
	  		}

	  		/* calculate the centroid */
	  		for (j=0;j<=n-1;j++) {
	  			cent=0.0;
	  			for (m=0;m<=n;m++) {
	  				if (m!=vg) {
	  					cent += v[m][j];
	  				}
	  			}
	  			vm[j] = cent/n;
	  		}

	  		/* reflect vg to new vertex vr */
	  		for (j=0;j<=n-1;j++) {
	  			/*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
	  			vr[j] = vm[j]+ALPHA*(vm[j]-v[vg][j]);
	  		}
	      
	      c.getConstrainedValues(vr,n);
	  		fr = objf.value(vr); num_evals++; maxed_eval = num_evals > MAX_EV;
	  		k++;

	  		if (fr < f[vh] && fr >= f[vs]) {
	  			for (j=0;j<=n-1;j++) {
	  				v[vg][j] = vr[j];
	  			}
	  			f[vg] = fr;
	  		}

	  		/* investigate a step further in this direction */
	  		if ( fr <  f[vs]) {
	  			for (j=0;j<=n-1;j++) {
	  				/*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
	  				ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
	  			}

	        c.getConstrainedValues(ve,n);
	  			fe = objf.value(ve); num_evals++; maxed_eval = num_evals > MAX_EV;
	  			k++;

	  			/* by making fe < fr as opposed to fe < f[vs], 			   
	  			   Rosenbrocks function takes 63 iterations as opposed 
	  			   to 64 when using double variables. */

	  			if (fe < fr) {
	  				for (j=0;j<=n-1;j++) {
	  					v[vg][j] = ve[j];
	  				}
	  				f[vg] = fe;
	  			}
	  			else {
	  				for (j=0;j<=n-1;j++) {
	  					v[vg][j] = vr[j];
	  				}
	  				f[vg] = fr;
	  			}
	  		}

	  		/* check to see if a contraction is necessary */
	  		if (fr >= f[vh]) {
	  			if (fr < f[vg] && fr >= f[vh]) {
	  				/* perform outside contraction */
	  				for (j=0;j<=n-1;j++) {
	  					/*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
	  					vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
	  				}

	          c.getConstrainedValues(vc,n);
	  				fc = objf.value(vc); num_evals++; maxed_eval = num_evals > MAX_EV;
	  				k++;
	  			}
	  			else {
	  				/* perform inside contraction */
	  				for (j=0;j<=n-1;j++) {
	  					/*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
	  					vc[j] = vm[j]-BETA*(vm[j]-v[vg][j]);
	  				}

	          c.getConstrainedValues(vc,n);
	  				fc = objf.value(vc); num_evals++; maxed_eval = num_evals > MAX_EV;
	  				k++;
	  			}


	  			if (fc < f[vg]) {
	  				for (j=0;j<=n-1;j++) {
	  					v[vg][j] = vc[j];
	  				}
	  				f[vg] = fc;
	  			}
	  			/* at this point the contraction is not successful,
	  			   we must halve the distance from vs to all the 
	  			   vertices of the simplex and then continue.
	  			   10/31/97 - modified to account for ALL vertices. 
	  			*/
	  			else {
	  				for (row=0;row<=n;row++) {
	  					if (row != vs) {
	  						for (j=0;j<=n-1;j++) {
	  							v[row][j] = v[vs][j]+(v[row][j]-v[vs][j])*SIGMA;
	  						}
	  					}
	  				}

	          c.getConstrainedValues(v[vg],n);
	  				f[vg] = objf.value(v[vg]); num_evals++; maxed_eval = num_evals > MAX_EV;
	  				k++;

	          c.getConstrainedValues(v[vh],n);
	  				f[vh] = objf.value(v[vh]); num_evals++; maxed_eval = num_evals > MAX_EV;
	  				k++;
	  			}
	  		}

	  		/* print out the value at each iteration */
	  		/*
	  		System.out.format("Iteration %d\n",itr);
	    	for (j=0;j<=n;j++) {
	    	  for (i=0;i<n;i++) {
	    		  System.out.format("%f %f\n",v[j][i],f[j]);
	  		  }
	    	}
	  		*/

	  		/* test for convergence */
	  		fsum = 0.0;
	  		for (j=0;j<=n;j++) {
	  			fsum += f[j];
	  		}
	  		favg = fsum/(n+1);
	  		s = 0.0;
	  		for (j=0;j<=n;j++) {
	  			s += Math.pow((f[j]-favg),2.0)/(n);
	  		}
	  		s = Math.sqrt(s);
	  		if (s < EPSILON) break;
	  	}
	  	/* end main loop of the minimization */

	  	/* find the index of the smallest value */
	  	vs=0;
	  	for (j=0;j<=n;j++) {
	  		if (f[j] < f[vs]) {
	  			vs = j;
	  		}
	  	}

	  	/*
	  	System.out.format("The minimum was found at\n"); 
	  	for (j=0;j<n;j++) {
	  		//System.out.format("%e\n",v[vs][j]);
	  		start[j] = v[vs][j];
	  	}
	  	*/
	  	
	  	k++;
	  	System.out.format("%d Function Evaluations\n",k);
	  	System.out.format("%d Iterations through program\n",itr);

	  	
	  	PointValuePair final_pvp = new PointValuePair(v[vs].clone(), f[vs]);
	  	return final_pvp;
	}

	////////////////////////////////////////////////////
	
	
	public class Constraint{
		 
		double round(double num, int precision){ // precision is digit to round all values to. 
			double rnum;
			int tnum;

			rnum = num*Math.pow(10,precision);
			tnum = (int)(rnum < 0 ? rnum-0.5 : rnum + 0.5);
			rnum = tnum/Math.pow(10,precision);

			return rnum;
		}
		
		public void getConstrainedValues(double x[], int n) {

		}
	}
	

		
}
