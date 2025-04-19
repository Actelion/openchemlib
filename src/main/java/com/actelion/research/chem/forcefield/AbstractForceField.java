package com.actelion.research.chem.forcefield;

import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;

public abstract class AbstractForceField implements ForceField {
	public static final double FUNCTOL = 1e-4;
    public static final double MOVETOL = 1e-7;
    public static final double EPS = 3e-8;
    public static final double TOLX = 4.0*EPS;
    public static final double MAXSTEP = 100.0;
    
	protected ArrayList<ForceFieldChangeListener> listeners = new ArrayList<ForceFieldChangeListener>();
    protected StereoMolecule mMol;
	protected final int mDim; ;
	protected double[] mPos;
	protected double[] mNewpos;
	protected double[] mGrad;
	protected int[] mFixedAtoms;
	protected double mTotalEnergy;
	protected long mTimeInterval; //time interval for
	protected volatile boolean mIsInterrupted;

	public AbstractForceField(StereoMolecule mol) {
		int implicitHydrogens = 0;
		for(int at=0;at<mol.getAtoms();at++) {
			implicitHydrogens += mol.getImplicitHydrogens(at);
		}
		if(implicitHydrogens>0) {
			throw new IllegalArgumentException("molecule needs explicit hydrogen atoms for force field calculations");
		}
		mMol = mol;
		mDim = 3*mol.getAllAtoms();
        mGrad = new double[mDim];
        mPos = new double[mDim];
        mNewpos = new double[mDim];
        mIsInterrupted = false;
        mTimeInterval = 20;

        // get the atom positions to be placed in the pos array.
        for (int i=0; i<mol.getAllAtoms(); i++) {
            mPos[3*i    ] = mol.getAtomX(i);
            mPos[3*i + 1] = mol.getAtomY(i);
            mPos[3*i + 2] = mol.getAtomZ(i);
        }
	}

	/**
	 * Return the variance across all atoms in a molecule
	 * for the specified coordinate.
	 *  @param listener
	 *  @return variance for the specified coordinate.
	 */
	@Override
	public void addListener(ForceFieldChangeListener listener) {
		listeners.add(listener);
	}
	
	public void addGradient(double[] grad) {
		assert grad.length==mGrad.length;
		updateGradient();
		for(int i=0;i<mGrad.length;i++) {
			grad[i] += mGrad[i];
		}
	}
	
	public void getState(double[] pos) {
		assert pos.length==mPos.length;
		for(int i=0;i<mPos.length;i++) {
			pos[i] = mPos[i];
		}
	}

	public void setState(double[] pos) {
		assert pos.length==mPos.length;
		for(int i=0;i<mPos.length;i++) {
			mPos[i] = pos[i];
		}
	}
	
	
	public double coordVariance(int c) {
        double m = 0.0;
        double s = 0.0;
        int k = 1;
        for (int i=0; i<mMol.getAllAtoms(); i++) {
            double v;
            switch (c) {
                case 0:
                    v = mMol.getAtomX(i);
                    break;
                case 1:
                    v = mMol.getAtomY(i);
                    break;
                default:
                    v = mMol.getAtomZ(i);
                    break;
            }
            double tm = m;
            m += (v - tm) / (double)k;
            s += (v - tm) * (v - m);
            k++;
        }
        return (k > 1 ? s / (double)(k - 1) : 0.0);
    }

    @Override
    public int minimise() {
        return minimise(4000, 1e-4, 1e-6);
    }
    
	@Override
	public void setFixedAtoms(int[] fixedAtoms) {
	   mFixedAtoms = fixedAtoms;
   }
   
 
	public void zeroGradient() {
		if (mFixedAtoms!=null) {
			for (int i:mFixedAtoms) {
				mGrad[3*i] = 0.0;
				mGrad[3*i+1] = 0.0;
				mGrad[3*i+2] = 0.0;
	   		}
		}
	}
	
    /**
     * Minimise the current molecule.
     *  @param maxIts The maximum number of iterations to run for.
     *  @param gradTol The gradient tolerance.
     *  @param funcTol The energy tolerance.
     *  @return Return code, 0 on success.
     */
    public int minimise(int maxIts, double gradTol, double funcTol) {
        int res = 1;

        for (int i=0; i<mMol.getAllAtoms(); i++) {
    	    mPos[3*i    ] = mMol.getAtomX(i); //+ delta[0];
    	    mPos[3*i + 1] = mMol.getAtomY(i); //+ delta[1];
    	    mPos[3*i + 2] = mMol.getAtomZ(i); //+ delta[2];
        }
        res = run_minimiser(maxIts, gradTol, funcTol);

        if (res == 0) {
            for (int i=0; i<mMol.getAllAtoms(); i++) {
                mMol.setAtomX(i, mPos[3*i  ]);
                mMol.setAtomY(i, mPos[3*i+1]);
                mMol.setAtomZ(i, mPos[3*i+2]);
            }
        }
   	    for(ForceFieldChangeListener listener: listeners) {
		    listener.stateChanged();
	    }
        return res;
    }

	public int run_minimiser(int maxIts, double gradTol, double funcTol) {
		double sum,maxStep,fp;
		mGrad = new double[mDim];
		double[] dGrad = new double[mDim];
		double[] hessDGrad = new double[mDim];
		double[] newPos = new double[mDim];
		double[] xi = new double[mDim];
		double[] invHessian = new double[mDim*mDim];
		for (int i=0; i<mDim; i++)
		         newPos[i] = mPos[i];
		// evaluate the function and gradient in our current position:
		fp = getTotalEnergy(mPos);
		updateGradient();
		zeroGradient();
		sum = 0.0;
		//memset(invHessian,0,dim*dim*sizeof(double));
		for (int i=0; i<mDim; i++) {
		     // initialize the inverse hessian to be identity
		     invHessian[i*mDim+i] = 1.0;
		     // the first line dir is -grad:
		     xi[i] = -mGrad[i];
		     sum += mPos[i]*mPos[i];
		}
		
		     // pick a max step size:
	     maxStep = MAXSTEP * Math.max(Math.sqrt(sum), mDim);
	     long timePassed;
	     long t0 = System.currentTimeMillis();
	     for (int iter=1; iter<=maxIts && !mIsInterrupted; iter++) {
	         // do the line search:
	         int status = linearSearch(mPos,fp,xi,newPos,maxStep);
	         if (status < 0)
	             return 2;
	
	         // save the function value for the next search:
	         fp = mTotalEnergy;
	
	         // set the direction of this line and save the gradient:
	         double test=0.0;
	         for (int i=0; i<mDim; i++) {
	             xi[i] = newPos[i] - mPos[i];
	             mPos[i] = newPos[i];
	             double temp = Math.abs(xi[i])/Math.max(Math.abs(mPos[i]),1.0);
	             if (temp > test)
	                 test = temp;
	             dGrad[i] = mGrad[i];
	         }
	
	         if (test < TOLX) {
	             return 0;
	         }
	
	         // update the gradient:
	         double gradScale = updateGradient();
	         zeroGradient();
	         // is the gradient converged?
	         test = 0.0;
	         double term = Math.max(mTotalEnergy*gradScale, 1.0);
	         for (int i=0; i<mDim; i++) {
	             double tmp = Math.abs(mGrad[i])*Math.max(Math.abs(mPos[i]), 1.0);
	             test = Math.max(test, tmp);
	             dGrad[i] = mGrad[i] - dGrad[i];
	         }
	
	         test /= term;
	
	         if (test < gradTol) {
	             return 0;
	         }
	
	         // compute hessian*dGrad:
	         double fac = 0, fae = 0, sumDGrad = 0, sumXi = 0;
	         for(int i=0; i<mDim; i++) {
	             int itab = i*mDim;
	             hessDGrad[i] = 0.0;
	
	             for (int j=0; j<mDim; j++)
	                 hessDGrad[i] += invHessian[itab+j] * dGrad[j];
	
	             fac += dGrad[i] * xi[i];
	             fae += dGrad[i] * hessDGrad[i];
	             sumDGrad += dGrad[i] * dGrad[i];
	             sumXi += xi[i] * xi[i];
	         }
	
	         if (fac > Math.sqrt(EPS*sumDGrad*sumXi)) {
	             fac = 1.0/fac;
	             double fad = 1.0/fae;
	             for (int i=0; i<mDim; i++)
	                 dGrad[i] = fac*xi[i] - fad*hessDGrad[i];
	
	             for (int i=0; i<mDim; i++) {
	                 int itab = i*mDim;
	                 for (int j=i; j<mDim; j++) {
	                     invHessian[itab+j] += fac*xi[i]*xi[j]
	                             - fad*hessDGrad[i]*hessDGrad[j]
	                             + fae*dGrad[i]*dGrad[j];
	                     invHessian[j*mDim+i] = invHessian[itab+j];
	                 }
	             }
	         }
	
	         // generate the next direction to move:
	         for (int i=0; i<mDim; i++) {
	             int itab = i*mDim;
	             xi[i] = 0.0;
	             for (int j=0; j<mDim; j++)
	                 xi[i] -= invHessian[itab+j]*mGrad[j];
	         }
	         long t1 = System.currentTimeMillis();
	         timePassed = t1-t0;
	         if(timePassed>=mTimeInterval) {
	        	 for(ForceFieldChangeListener listener: listeners) {
	        		 listener.stateChanged();
	        	 }
	        	 t0=t1;
	         }
	     }
	     return 1;
	 }

 /**
  *
  */
	 private int linearSearch(double[] oldPt,
	     double oldVal,
	     double[] dir,
	     double[] newPt,
	     double maxStep) {
	     final int MAX_ITER_LINEAR_SEARCH = 1000;
	     int ret = -1;
	     double [] tmpPt = new double[mDim];
	     double sum = 0.0, slope = 0.0, test = 0.0, lambda = 0.0;
	     double lambda2 = 0.0, lambdaMin = 0.0, tmpLambda = 0.0, val2 = 0.0;
	
	     // get the length of the direction vector:
	     sum = 0.0;
	     for (int i=0; i<mDim; i++)
	         sum += dir[i]*dir[i];
	     sum = Math.sqrt(sum);
	
	     // rescale if we're trying to move too far:
	     if (sum > maxStep)
	         for (int i=0; i<mDim; i++)
	             dir[i] *= maxStep/sum;
	
	     // make sure our direction has at least some component along
	     // -grad
	     slope = 0.0;
	     for (int i=0; i<mDim; i++)
	         slope += dir[i]*mGrad[i];
	
	     if (slope >= 0.0)
	         return ret;
	
	     test = 0.0;
	     for (int i=0; i<mDim; i++) {
	         double temp = Math.abs(dir[i])/Math.max(Math.abs(oldPt[i]),1.0);
	         if (temp > test)
	             test=temp;
	     }
	
	     lambdaMin = MOVETOL/test;
	     lambda = 1.0;
	     int it = 0;
	     while (it < MAX_ITER_LINEAR_SEARCH) {
	         if (lambda < lambdaMin) {
	             // the position change is too small.
	             ret = 1;
	             break;
	         }
	
	         for(int i=0; i<mDim; i++)
	             newPt[i]=oldPt[i]+lambda*dir[i];
	         mTotalEnergy = getTotalEnergy(newPt);
	
	         // we're converged on the function:
	         if (mTotalEnergy-oldVal <= FUNCTOL*lambda*slope)
	             return 0;
	
	         // if we made it this far, we need to backtrack:
	         // it's the first step:
	         if (it == 0)
	             tmpLambda = -slope / (2.0*(mTotalEnergy - oldVal - slope));
	         else {
	             double rhs1 = mTotalEnergy - oldVal - lambda*slope;
	             double rhs2 = val2 - oldVal - lambda2*slope;
	             double a = (rhs1/(lambda*lambda) - rhs2/(lambda2*lambda2))
	                     /(lambda-lambda2);
	             double b = (-lambda2*rhs1/(lambda*lambda)
	                     + lambda*rhs2/(lambda2*lambda2))/(lambda-lambda2);
	             if (a == 0.0)
	                 tmpLambda = -slope/(2.0*b);
	             else {
	                 double disc = b*b-3*a*slope;
	                 if (disc < 0.0)
	                     tmpLambda = 0.5*lambda;
	                 else if (b <= 0.0)
	                     tmpLambda = (-b + Math.sqrt(disc))/(3.0*a);
	                 else
	                     tmpLambda = -slope/(b + Math.sqrt(disc));
	             }
	
	             if (tmpLambda > 0.5*lambda)
	                 tmpLambda = 0.5*lambda;
	         }
	
	         lambda2 = lambda;
	         val2 = mTotalEnergy;
	         lambda = Math.max(tmpLambda, 0.1*lambda);
	         ++it;
	     }
	     // nothing was done
	     for(int i=0; i<mDim; i++)
	         newPt[i]=oldPt[i];
	     return ret;
	}

	@Override
	public void interrupt() {
		 mIsInterrupted = true;
	 }
}
