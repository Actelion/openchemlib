package com.actelion.research.chem.alignment3d.transformation;



/**
 * Rotation derivatives for rotations expressed by means of the exponential map
 * http://www.cs.cmu.edu/~spiff/exp-map/
 * Grassi,98
 * --> we take the transpose because of different convention!
 */
 

public class RotationDerivatives {
	private double[] v;
    private double theta;
    private double cosp;
    private double sinp = Math.sin(.5*theta);
    Quaternion q;
    
    public RotationDerivatives(double[] v) {
    	this.v = v;
        theta = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        cosp = Math.cos(.5*theta);
        sinp = Math.sin(.5*theta);
    	ExponentialMap em = new ExponentialMap(v[0],v[1],v[2]);
        q = em.toQuaternion();
    }

	
	/* -----------------------------------------------------------------
	 * 'Partial_R_Partial_Vi' Given a quaternion 'q' computed from the
	 * current 2 or 3 degree of freedom EM vector 'v', and the partial
	 * derivative of the quaternion with respect to the i'th element of
	 * 'v' in 'dqdvi' (computed using 'Partial_Q_Partial_3V' or
	 * 'Partial_Q_Partial_2V'), compute and store in 'dRdvi' the i'th
	 * partial derivative of the rotation matrix 'R' with respect to the
	 * i'th element of 'v'.
	 * -----------------------------------------------------------------*/
	void dRdvi(Quaternion dqdvi, double[][] dRdvi){
	    double[] prod = new double[9];

	 	 /* This efficient formulation is arrived at by writing out the
	     * entire chain rule product dRdq * dqdv in terms of 'q' and 
	     * noticing that all the entries are formed from sums of just
	     * nine products of 'q' and 'dqdv' */
	    prod[0] = -4*q.getQ1()*dqdvi.getQ1();
	    prod[1] = -4*q.getQ2()*dqdvi.getQ2();
	    prod[2] = -4*q.getQ3()*dqdvi.getQ3();
	    prod[3] = 2*(q.getQ2()*dqdvi.getQ1() + q.getQ1()*dqdvi.getQ2());
	    prod[4] = 2*(q.getQ0()*dqdvi.getQ3() + q.getQ3()*dqdvi.getQ0());
	    prod[5] = 2*(q.getQ3()*dqdvi.getQ1() + q.getQ1()*dqdvi.getQ3());
	    prod[6] = 2*(q.getQ0()*dqdvi.getQ2() + q.getQ2()*dqdvi.getQ0());
	    prod[7] = 2*(q.getQ3()*dqdvi.getQ2() + q.getQ2()*dqdvi.getQ3());
	    prod[8] = 2*(q.getQ0()*dqdvi.getQ1() + q.getQ1()*dqdvi.getQ0());

	    /* first row, followed by second and third */
	    dRdvi[0][0] = prod[1] + prod[2];
	    dRdvi[1][0] = prod[3] - prod[4];
	    dRdvi[2][0] = prod[5] + prod[6];

	    dRdvi[0][1] = prod[3] + prod[4];
	    dRdvi[1][1] = prod[0] + prod[2];
	    dRdvi[2][1] = prod[7] - prod[8];

	    dRdvi[0][2] = prod[5] - prod[6];
	    dRdvi[1][2] = prod[7] + prod[8];
	    dRdvi[2][2] = prod[0] + prod[1];

	    /* the 4th row and column are all zero */

	} 
	
	public void dRdv(int i, double[][] dRdvi){
	    Quaternion dqdvi;
	    ExponentialMap em = new ExponentialMap(v[0],v[1],v[2]);
	    q = em.toQuaternion();
	    dqdvi = dqdvi(i);
	    dRdvi(dqdvi,dRdvi);

	}


	
	
	


	
	/**
	 * derivatives of the quaternions with respect to the exponential map representation (vector p): 
	 * described in: F. Sebastian Grassia (1998) Practical Parameterization of Rotations Using the
	 *	Exponential Map, Journal of Graphics Tools, 3:3, 29-48, DOI: 10.1080/10867651.1998.10487493
	 * @param p0
	 * @param p1
	 * @param p2
	 * @return
	 */
	private Quaternion dqdvi(int i)
	{
	    double w;
	    double[] xyz = new double[3];
	    
	    assert(i>=0 && i<3);

	    /* This is an efficient implementation of the derivatives given
	     * in Appendix A of the paper with common subexpressions factored out */
	    if (theta < ExponentialMap.EPSILON){
			final int i2 = (i+1)%3; 
			final int i3 = (i+2)%3;
			double Tsinc = 0.5 - theta*theta/48.0;
			double vTerm = v[i] * (theta*theta/40.0 - 1.0) / 24.0;
			
			w = -.5*v[i]*Tsinc;
			xyz[i]  = v[i]* vTerm + Tsinc;
			xyz[i2] = v[i2]*vTerm;
			xyz[i3] = v[i3]*vTerm;
	    }
	    else{
			final int i2 = (i+1)%3; 
			final int i3 = (i+2)%3;
			final double  ang = 1.0/theta;
			final double ang2 = ang*ang*v[i]; 
		    final double sang = sinp*ang;
			final double  cterm = ang2*(.5*cosp - sang);
			
			xyz[i]  = cterm*v[i] + sang;
			xyz[i2] = cterm*v[i2];
			xyz[i3] = cterm*v[i3];
			w = -.5*v[i]*sang;
	    }
	    return new Quaternion(w,xyz[0],xyz[1],xyz[2]);
	}

		
		
	

	

}
