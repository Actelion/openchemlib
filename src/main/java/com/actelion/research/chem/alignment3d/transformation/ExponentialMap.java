package com.actelion.research.chem.alignment3d.transformation;

import com.actelion.research.chem.Coordinates;

/**
 * described in:
 * F. Sebastian Grassia (1998) Practical Parameterization of Rotations Using the
 * Exponential Map, Journal of Graphics Tools, 3:3, 29-48, DOI: 10.1080/10867651.1998.10487493
 * @author wahljo1
 *
 */

public class ExponentialMap {

	public static final double EPSILON = 1E-6;
	public static final double CUTOFF_ANGLE = Math.PI;
	private Coordinates p;
	
	public ExponentialMap(Quaternion q) {
		calc(q);
		checkParameterization();
	}
	
	public ExponentialMap(Coordinates p) {
		this.p = p;
		checkParameterization();
	}
	
	public ExponentialMap(double p1, double p2, double p3) {
		this.p = new Coordinates(p1,p2,p3);
	}
	
	public Coordinates getP() {
		return p;
	}


	/**
	 * https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation#Unit_quaternions: accessed May 2021
	 * @param q
	 */
	private void calc(Quaternion q) {
		Coordinates axis = new Coordinates(0,0,0);
		Coordinates qv = new Coordinates(q.getQ1(),q.getQ2(),q.getQ3());
		double qw = q.getQ0();
		double qvLength = qv.getLength();
		double angle = 2*Math.atan2(qvLength, qw);
		if(angle>EPSILON)
			axis = qv.scaleC(1.0/Math.sin(angle/2.0));
		axis.scale(angle);
		this.p = axis;
	}
	

	public Quaternion toQuaternion() {

	    double   cosp, sinp, theta;
	    theta = p.getLength();
	    
	    cosp = Math.cos(.5*theta);
	    sinp = Math.sin(.5*theta);

	    double scalar = cosp;
	    Coordinates v;
	    if (theta < EPSILON)
	      v = p.scaleC(0.5-theta*theta/48.0);

	    else
	      v = p.scaleC(sinp/theta);


	    return new Quaternion(scalar,v.x,v.y,v.z);
}
	
	
	/* -----------------------------------------------------------------
	 * 'Check_Parameterization' To escape the vanishing derivatives at
	 * shells of 2PI rotations, we reparameterize to a rotation of (2PI -
	 * theta) about the opposite axis when we get too close to 2PI
	 * -----------------------------------------------------------------*/
	private void checkParameterization(){
	    double theta = p.getLength();

	    if (theta > CUTOFF_ANGLE){
			double scl = theta;
			if (theta > 2*Math.PI){	/* first get theta into range 0..2PI */
			    theta = theta%(2*Math.PI);
			    scl = theta/scl;
			    p.scale(scl);
	
			}
			if (theta > CUTOFF_ANGLE){
			    scl = theta;
			    theta = 2*Math.PI - theta;
			    scl = 1.0 - 2*Math.PI/scl;
			    p.scale(scl);

			}
	    }

	}
	
	
	
	public static void main(String[] args) {
		double c = 0.707106781;

		double[][] transforms2 = {{1,0,0,0},{0,1,0,0},{0,0,1,0},
				{0,0,0,1},
				{c,c,0,0},
				{c,0,c,0},
				{c,0,0,c},
				{-0.5,0.5,0.5,-0.5},
				{0.5,-0.5,0.5,-0.5},
				{0.5,0.5,0.5,-0.5},
				{0.5,-0.5,-0.5,-0.5},
				{0.5,0.5,-0.5,-0.5}
				};
		for(double[] transform : transforms2) {
			ExponentialMap eMap = new ExponentialMap(new Quaternion(transform[0],transform[1],
					transform[2],transform[3]));
			System.out.println(eMap.toQuaternion());
		}
	}
	

}
