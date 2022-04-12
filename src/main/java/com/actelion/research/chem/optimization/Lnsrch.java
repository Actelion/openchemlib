package com.actelion.research.chem.optimization;

import com.actelion.research.chem.phesa.EvaluableOverlap;
import java.util.Arrays;





/**
 * taken from DD_chem3d
 *  
 */
class Lnsrch{

	
	/**
	 * Minimize the function according to the line search algorithm. (ie find lamda so that newpos = pos+lambda*dir minimizes the energy)
	 * This function expects the initial FGValue, the direction of search and the initial move
	 * 
	 *  http://www.mcs.anl.gov/~anitescu/CLASSES/2012/LECTURES/S310-2012-lect4.pdf
	 *  
	 * @param f0
	 * @param function
	 * @param grad
	 * @param dir
	 * @param fMove
	 * @return
	 */
	public static final Object[] minimizeEnergyAroundDirection(final Evaluable function, double f0, final double[] grad, final double[] dir, final double fMove) {
		final double CAPPA = .9;
		final double STPMIN = 1e-6;
		final double STPMAX = .1;
		double fA=0, fB=0, fC=0;
		double slopeA=0, slopeB=0, slopeC=0;
		double cube = 0;
		//Compute length of Gradient
		final double len = grad.length;
		final double sNorm = OptimizerLBFGS.getNorm(dir);

		//Normalize the search vector and find the projected gradient
		double slope = 0;
		for(int i=0; i<len; i++) {
			dir[i] /= sNorm;
			slope += dir[i]*grad[i];
		}
			
		double step = Math.min(2.0*Math.abs(fMove/slope), sNorm);
		if(step>STPMAX) step = STPMAX;
		else if(step<STPMIN) step = STPMIN;
		
		double lambda = 0;
		double[] initial = function.getState();
		double initialF = f0;
		double[] v = new double[initial.length];
		try {
			for(int reSearch=0; reSearch<2; reSearch++) {
				fB = f0;
				slopeB = slope;
	
				//Quadratic interpolation: f(r) = a.r^2 + b.r + c  ( a = (gA-gB)/(2A-2B), b = gA - (gA-gB)/(A-B)) 
				//
				//    A
				//     \          /B
				//       -       /
				//         -- --
				//
				//    0            lambda
				for(int counter = 0; counter<20; counter++) {
					fA = fB;
					slopeA = slopeB;
					
					//evaluate at lambda+step
					lambda += step;
					move(function, dir, lambda, initial,v);
					fB = function.getFGValue(grad);
					slopeB = 0; for(int i=0; i<len; i++) slopeB += dir[i]*grad[i];

					if(Math.abs(slopeB/slope)<=CAPPA && fB<=fA) { //success
						return new Object[]{fB, grad, Boolean.TRUE};
					}				
					
					//go to cubic interpolation if gradient changes sign or function increases
					if(fB>=fA || slopeB*slopeA<0) {
						break;
					}
					
					//Adapt step		
					if(slopeB>slopeA) {
						double parab = (fA - fB) / (slopeB - slopeA);
						if(parab>2*step) step = 2 * step;
						else if(parab<2*step) step = step / 2;
						else step = parab;
					} else {
						step*=2;
					}
					if(step>STPMAX) step = STPMAX;
					else if(step<STPMIN) step = STPMIN;
				} // end-while

				fC = fB;
				//Cubic interpolation (http://www.mathworks.com/access/helpdesk_r13/help/toolbox/optim/tutori5b.html)

				//Cubic Interpolation has failed, reset to best current point

				double fL, sgL;
				if(fA<=fB && fA<=fC) { //A is min
					fL = fA;
					sgL = slopeA;
					lambda += cube-step;					
				} else if(fB<=fA && fB<=fC) { //B is min
					fL = fB;
					sgL = slopeB;
					lambda += cube;					
				} else /*fC<=fA && fC<=fB*/ { //C is min
					fL = fC;
					sgL = slopeC;
				}

				//try to restart from best point with smaller stepsize
				if(fL>f0) {
					move(function, dir, lambda, initial,v);
					f0 = function.getFGValue(grad);
					//System.err.println("ERR fL>f0");
					return new Object[]{f0, grad, Boolean.FALSE};
				}
				f0 = fL;

				if(sgL>0) {
					lambda = 0;
					for(int i=0; i<dir.length; i++) dir[i] = -dir[i];
					slope -= -sgL;
				}  else {
					slope = sgL;					
				}
				step = Math.max( STPMIN, Math.max(cube, step-cube) / 10);				
			}
			
			//Already restarted, return best current point
			move(function, dir, lambda, initial,v);
			f0 = function.getFGValue(grad);
			if(f0>initialF) {
				move(function, dir, 0, initial,v);
				f0 = initialF;
			}
			return new Object[]{f0, grad, Boolean.FALSE};

		} catch(Exception e) {
			e.printStackTrace();
			function.setState(initial);
			f0 = function.getFGValue(grad);
			return new Object[]{new Double(f0), grad, Boolean.FALSE};				
		}
	}	

	private final static void move(Evaluable eval, double[] dir, double lambda, double[] transformOld, double[] transform) {
		for(int i=0;i<transform.length;i++){
			transform[i] = transformOld[i] + lambda*dir[i];

		}

		eval.setState(transform);


}
	
}
