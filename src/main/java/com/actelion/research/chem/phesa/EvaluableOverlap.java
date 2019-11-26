package com.actelion.research.chem.phesa;

import java.util.ArrayList;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.phesa.pharmacophore.PPGaussian;



/**
 * @author J.Wahl, February 2018
 * describes the overlap function, the state (relative orientation, translation of the two molecules)
 * returns the objective function and the gradient of the overlap with respect to translation and rotation
 * accessed by the optimization algorithm
 */


public class EvaluableOverlap implements Evaluable  {

	private static final int PENALTY = 80; 

	private PheSAAlignment shapeAlign;
	private double[] transform;
    private double [][] qDersAt;
    private double [][] rDersAt;
    private double [][] sDersAt;
    private double [][] uDersAt;
    private double [][] qDersPP;
    private double [][] rDersPP;
    private double [][] sDersPP;
    private double [][] uDersPP;
    private double[][] results;
    private Coordinates[] fitAtGaussModCoords;
    private Coordinates[] fitPPGaussModCoords;
    private Coordinates[] fitPPDirectionalityMod;

    
    public EvaluableOverlap(PheSAAlignment shapeAlign, double[] transform) {
		this.shapeAlign = shapeAlign; 
		this.transform = transform;
	    this.fitAtGaussModCoords = new Coordinates[shapeAlign.getMolGauss().getAtomicGaussians().size()];
	    this.fitPPGaussModCoords = new Coordinates[shapeAlign.getMolGauss().getPPGaussians().size()];
	    this.qDersAt = new double[fitAtGaussModCoords.length][3];
	    this.rDersAt = new double[fitAtGaussModCoords.length][3];
	    this.sDersAt = new double[fitAtGaussModCoords.length][3];
	    this.uDersAt = new double[fitAtGaussModCoords.length][3];
	    this.qDersPP = new double[fitPPGaussModCoords.length][3];
	    this.rDersPP = new double[fitPPGaussModCoords.length][3];
	    this.sDersPP = new double[fitPPGaussModCoords.length][3];
	    this.uDersPP = new double[fitPPGaussModCoords.length][3];
	    this.fitPPDirectionalityMod = new Coordinates[shapeAlign.getMolGauss().getPPGaussians().size()];
		this.results = new double[shapeAlign.getRefMolGauss().getAtomicGaussians().size()][shapeAlign.getRefMolGauss().getAtomicGaussians().size()];
		
	}
	
	public EvaluableOverlap(EvaluableOverlap e) {
		this.shapeAlign = e.shapeAlign;
		this.transform = e.transform;	
		this.qDersAt = e.qDersAt;
		this.rDersAt = e.rDersAt;
		this.sDersAt = e.sDersAt;
		this.uDersAt = e.uDersAt;
		this.qDersPP = e.qDersPP;
		this.rDersPP = e.rDersPP;
		this.sDersPP = e.sDersPP;
		this.uDersPP = e.uDersPP;
		this.fitAtGaussModCoords = e.fitAtGaussModCoords;
		this.fitPPGaussModCoords = e.fitPPGaussModCoords;
		this.results = e.results;

	}
	
	@Override
	public void setState(double[] transform){
		this.transform=transform;
	}
	
	public double[] getState(double[] v){
		for(int i=0;i<this.transform.length;i++) {
			v[i] = transform[i];
			
		}
		return v;
	}
	
	public double[] getState() {
		return this.getState(new double[transform.length]);
	}
	
	

	
	public PheSAAlignment getAlignment() {
		return this.shapeAlign;
	}
	


	@Override
	public double getFGValue(double[] grad) {
		MolecularVolume refMolGauss = shapeAlign.getRefMolGauss();
		MolecularVolume fitMolGauss = shapeAlign.getMolGauss();
		double value = 0.0;
		double[] atomGrad = new double[grad.length];
		value += this.getFGValueOverlap(atomGrad,refMolGauss.getAtomicGaussians(),fitMolGauss.getAtomicGaussians(), refMolGauss.getExclusionGaussians(),
						qDersAt,rDersAt,sDersAt,uDersAt,fitAtGaussModCoords);
			
		
		double[] ppGrad = new double[grad.length];
		value += this.getFGValueOverlapPP(ppGrad,refMolGauss.getPPGaussians(),fitMolGauss.getPPGaussians(),
						qDersPP,rDersPP,sDersPP,uDersPP,fitPPGaussModCoords,fitPPDirectionalityMod);

		for(int i=0;i<grad.length;i++) 
			grad[i] = atomGrad[i]+ ppGrad[i];
				

		return value;
		
		
	}
	
	
	private void getQuatGradient(double[][] q0Ders, double[][] q1Ders, double[][] q2Ders, double[][] q3Ders, ArrayList<? extends Gaussian3D> refMolGauss,ArrayList<? extends Gaussian3D> fitMolGauss,Coordinates[] fitModCoords,
			double q0, double q1, double q2, double q3, double invnorm2) {

		    /**
		     * we first calculate the partial derivatives with respect to the four elements of the quaternion q,r,s,u
		     * the final gradient has 7 elements, the first four elements are the gradients for the quaternion (rotation),
		     * the last three elements are for the translation
		     */
	
			int i=0;
		    for(Coordinates fitCenterModCoord:fitModCoords){
		        double xk=fitCenterModCoord.x;
		        double yk=fitCenterModCoord.y;
		        double zk=fitCenterModCoord.z;   
		        
		        double dxdq0 =  invnorm2*2.0*(q0*xk - q3*yk + q2*zk);
		        double dydq0 =  invnorm2*2.0*(q3*xk + q0*yk - q1*zk);
		        double dzdq0 =  invnorm2*2.0*(-q2*xk + q1*yk + q0*zk);
		        
		        double dxdq1 =  invnorm2*2.0*(q1*xk + q2*yk + q3*zk);
		        double dydq1 =  invnorm2*2.0*(q2*xk - q1*yk - q0*zk);
		        double dzdq1 =  invnorm2*2.0*(q3*xk + q0*yk - q1*zk);
		        
		        double dxdq2 =  invnorm2*2.0*(-q2*xk + q1*yk + q0*zk);
		        double dydq2 =  invnorm2*2.0*(q1*xk + q2*yk + q3*zk);
		        double dzdq2 =  invnorm2*2.0*(-q0*xk + q3*yk - q2*zk);
		        
		        double dxdq3 =  invnorm2*2.0*(-q3*xk - q0*yk + q1*zk);
		        double dydq3 =  invnorm2*2.0*(q0*xk - q3*yk + q2*zk);
		        double dzdq3 =  invnorm2*2.0*(q1*xk + q2*yk + q3*zk);

		        
		        q0Ders[i][0] = dxdq0;
		        q0Ders[i][1] = dydq0;
		        q0Ders[i][2] = dzdq0;
		        q1Ders[i][0] = dxdq1;
		        q1Ders[i][1] = dydq1;
		        q1Ders[i][2] = dzdq1;
		        q2Ders[i][0] = dxdq2;
		        q2Ders[i][1] = dydq2;
		        q2Ders[i][2] = dzdq2;
		        q3Ders[i][0] = dxdq3;
		        q3Ders[i][1] = dydq3;
		        q3Ders[i][2] = dzdq3;
		        i+=1;
		    }
		
	}
	


	/**
	 * calculates the gradient of the overlap function with respect to the three components of translation (dx,dy,dz)
	 * and the quaternion with elements q,r,s,u composing the rotation q,r,s,u
     * derivatives are described in: Griewank, Markey and Evans, The Journal of Chemical Physics, 71, 3449, 1979
     * the intersection volume of two Atomic Gaussians is given by (Grant, Gallardo and Pickup, Journal of Computational Chemistry,16,1653,1996 
     * 
     *                                pi           3/2                 alpha_i*alpha_j*Rij(T)**2
     * equation 1: Vij = p_i*p_j*(---------------)        * exp( - -------------------- )      
     *                            alpha_i + alpha_j                     alpha_i + alpha_j 
     * 
     * Rij is the distance between the two atomic Gaussians and depends on the transformation T (orientation, translation) of the molecule to be fitted
     * we can use the chain rule:
     * dVij/dT = dVij/dRij * dRij/dT 
     *  
     * dVij/dRij = 2*(alpha_i*alpha_j)/(alpha_i + alpha_j) * Rij * Vij
     * 
     * Rij = T(j) - i  the transformation T consists of a rotation R and a translation t, the rotation is described by a rotation matrix R(q) that can be expressed 
     * by means of a quaternion q(q,r,s,u)
     * Rij=(Xi-R(q)*Xj-t)   Xi are the coordinates of Atomic Gaussian i  
     * dRij/dt = -1
     * for the rotational part, we need the derivation of the rotational Matrix R(q) with respect to q,r,s and u  
     * dRij/dq = -dR(q)/dq *Xj  and accordingly for r,s,u
     * dR(q)/dq yields a 3x3 matrix, multiplied by a vector with 3 elements (Xj), therefore every derivative of the rotation with respect to the elements of the quaternion has three elements:
     * dxdq,dydq,dzdq usw. 
     * the total derivative is then the sum over all pairs of overlapping Atomic Gaussians
     * 
     * to force the quaternions into unity, a penalty term is added for deviation from unity
	 * @param grad 
	 */
	
	private double getFGValueOverlap(double[] grad,ArrayList<AtomicGaussian> refMolGauss,ArrayList<AtomicGaussian> fitMolGauss,ArrayList<ExclusionGaussian> exclusionGaussians,
			double[][] qDers,double[][] rDers,double[][] sDers,double[][] uDers, Coordinates[] fitGaussModCoords) {
		double q=transform[0];
	    double r=transform[1];
	    double s=transform[2];
	    double u=transform[3];
	    Quaternion quat = new Quaternion(q,r,s,u);
	    double norm2 = quat.normSquared();
	    double norm = Math.sqrt(norm2);
	    double invnorm2 = 1.0/norm2;
	    double invnorm = 1/norm;
	    
	    //System.out.println(Arrays.toString(transform));

	    /**
	     * we first calculate the partial derivatives with respect to the four elements of the quaternion q,r,s,u
	     * the final gradient has 7 elements, the first four elements are the gradients for the quaternion (rotation),
	     * the last three elements are for the translation
	     */
	    
	    double[][] rotMatrix = quat.getRotMatrix().getArray();


	    for(int k=0;k<fitMolGauss.size();k++) {
	    	fitGaussModCoords[k]=  fitMolGauss.get(k).getRotatedCenter(rotMatrix, invnorm2, new double[] {transform[4],transform[5],transform[6]});    //we operate on the transformed coordinates of the molecule to be fitted

	    }

	    this.getQuatGradient(qDers, rDers, sDers, uDers, refMolGauss, fitMolGauss, fitGaussModCoords,q,r,s,u,invnorm2);


		/**
		 * derivative of ShapeOverlap with respect to the four elements of the quaternion and three elements of translation
		 * 
		 */
	    
	    double totalOverlap = 0.0;
	    Coordinates fitCenterModCoord;
		for(int i=0; i<refMolGauss.size();i++){
			Gaussian3D refAt = refMolGauss.get(i);
			for(int j=0; j<fitMolGauss.size();j++){
				double atomOverlap = 0.0;
				Gaussian3D fitAt = fitMolGauss.get(j);
				fitCenterModCoord = fitGaussModCoords[j];
				double alphaSum = refAt.getWidth() + fitAt.getWidth();

				double dx = refAt.getCenter().x-fitCenterModCoord.x;
				double dy = refAt.getCenter().y-fitCenterModCoord.y;
				double dz = refAt.getCenter().z-fitCenterModCoord.z;
				double Rij2 = dx*dx + dy*dy + dz*dz;
				if(Rij2>=Gaussian3D.DIST_CUTOFF) 
					continue;
				atomOverlap = refAt.getHeight()*fitAt.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refAt.getWidth() * fitAt.getWidth()* Rij2)/alphaSum) *
							QuickMathCalculator.getInstance().getPrefactor(refAt.getAtomicNo(),fitAt.getAtomicNo());
					
				
				if (atomOverlap>0.0) {
					totalOverlap += atomOverlap;
					double gradientPrefactor = atomOverlap*-2*refAt.getWidth()*fitAt.getWidth()/(refAt.getWidth()+fitAt.getWidth());
					double qder = qDers[j][0]*dx+qDers[j][1]*dy+qDers[j][2]*dz; 
					double rder = rDers[j][0]*dx+rDers[j][1]*dy+rDers[j][2]*dz; 
					double sder = sDers[j][0]*dx+sDers[j][1]*dy+sDers[j][2]*dz; 
					double uder = uDers[j][0]*dx+uDers[j][1]*dy+uDers[j][2]*dz; 

					grad[0] += gradientPrefactor*qder;
					grad[1] += gradientPrefactor*rder;
					grad[2] += gradientPrefactor*sder;
					grad[3] += gradientPrefactor*uder;
					grad[4] += gradientPrefactor*dx;
					grad[5] += gradientPrefactor*dy;
					grad[6] += gradientPrefactor*dz;
				    }



				}
		}
			for(int k=0; k<exclusionGaussians.size();k++){
				ExclusionGaussian refEx = exclusionGaussians.get(k);
				for(int j=0; j<fitMolGauss.size();j++){
					double atomOverlap = 0.0;
					Gaussian3D fitAt = fitMolGauss.get(j);
					fitCenterModCoord = fitGaussModCoords[j];
					double alphaSum = refEx.getWidth() + fitAt.getWidth();

					double dx = refEx.getCenter().x-fitCenterModCoord.x;
					double dy = refEx.getCenter().y-fitCenterModCoord.y;
					double dz = refEx.getCenter().z-fitCenterModCoord.z;
					double Rij2 = dx*dx + dy*dy + dz*dz;

					if(Rij2>=Gaussian3D.DIST_CUTOFF) 
						continue;
					atomOverlap = -1*refEx.getHeight()*fitAt.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refEx.getWidth() * fitAt.getWidth()* Rij2)/alphaSum) *
								QuickMathCalculator.getInstance().getPrefactor(refEx.getAtomicNo(),fitAt.getAtomicNo());
						
					
					if (atomOverlap>0.0) {
						totalOverlap += atomOverlap;
						double gradientPrefactor = atomOverlap*-2*refEx.getWidth()*fitAt.getWidth()/(refEx.getWidth()+fitAt.getWidth());
						double qder = qDers[j][0]*dx+qDers[j][1]*dy+qDers[j][2]*dz; 
						double rder = rDers[j][0]*dx+rDers[j][1]*dy+rDers[j][2]*dz; 
						double sder = sDers[j][0]*dx+sDers[j][1]*dy+sDers[j][2]*dz; 
						double uder = uDers[j][0]*dx+uDers[j][1]*dy+uDers[j][2]*dz; 

						grad[0] += gradientPrefactor*qder;
						grad[1] += gradientPrefactor*rder;
						grad[2] += gradientPrefactor*sder;
						grad[3] += gradientPrefactor*uder;
						grad[4] += gradientPrefactor*dx;
						grad[5] += gradientPrefactor*dy;
						grad[6] += gradientPrefactor*dz;
					    }


					}
				
			}
		grad[0] += PENALTY*(1-invnorm)*this.transform[0]; //penalty term to force quaternion into unity
		grad[1] += PENALTY*(1-invnorm)*this.transform[1];
		grad[2] += PENALTY*(1-invnorm)*this.transform[2];
		grad[3] += PENALTY*(1-invnorm)*this.transform[3];




		return (-1.0*totalOverlap+0.5*PENALTY*(norm-1)*(norm-1)); //the negative overlap is returned as the objective, since we minimize the objective in the optimization algorithm

	
	}
	    
	    
	   private double getFGValueOverlapPP(double[] grad, ArrayList<PPGaussian> refMolGauss,ArrayList<PPGaussian> fitMolGauss, double[][] qDers,double[][] rDers,double[][] sDers,double[][] uDers, Coordinates[] fitGaussModCoords, Coordinates[] fitPPDirectionalityMod) {
			double q=transform[0];
		    double r=transform[1];
		    double s=transform[2];
		    double u=transform[3];
		    Quaternion quat = new Quaternion(q,r,s,u);
		    double norm2 = quat.normSquared();
		    double norm = Math.sqrt(norm2);
		    double invnorm2 = 1.0/norm2;
		    double invnorm = 1/norm;
		    
		    double[][] rotMatrix = quat.getRotMatrix().getArray();


		    for(int k=0;k<fitMolGauss.size();k++) {
		    	fitGaussModCoords[k]=  fitMolGauss.get(k).getRotatedCenter(rotMatrix, invnorm2, new double[] {transform[4],transform[5],transform[6]});    //we operate on the transformed coordinates of the molecule to be fitted
		    	fitPPDirectionalityMod[k] = fitMolGauss.get(k).getRotatedDirectionality(rotMatrix, invnorm2);
		    }
		    this.getQuatGradient(qDers, rDers, sDers, uDers, refMolGauss, fitMolGauss, fitGaussModCoords,q,r,s,u,invnorm2);


			/**
			 * derivative of ShapeOverlap with respect to the four elements of the quaternion and three elements of translation
			 * 
			 */
		    
		    double totalOverlap = 0.0;
		    Coordinates fitCenterModCoord;
			for(int i=0; i<refMolGauss.size();i++){
				PPGaussian refAt = refMolGauss.get(i);
				for(int j=0; j<fitMolGauss.size();j++){
					PPGaussian fitAt = fitMolGauss.get(j);
					Coordinates fitPPDirectionalityVector = fitPPDirectionalityMod[j];
					double atomOverlap = 0.0;
					fitCenterModCoord = fitGaussModCoords[j];
					double alphaSum = refAt.getWidth() + fitAt.getWidth();
					double xi = refAt.getCenter().x;
					double yi = refAt.getCenter().y;
					double zi = refAt.getCenter().z;
					double dx = refAt.getCenter().x-fitCenterModCoord.x;
					double dy = refAt.getCenter().y-fitCenterModCoord.y;
					double dz = refAt.getCenter().z-fitCenterModCoord.z;
					double Rij2 = dx*dx + dy*dy + dz*dz;

					if(Rij2>=Gaussian3D.DIST_CUTOFF) {
						continue;
					}
					atomOverlap = refAt.getWeight()*fitAt.getWeight()*refAt.getHeight()*fitAt.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refAt.getWidth() * fitAt.getWidth()* Rij2)/alphaSum) *
							QuickMathCalculator.getInstance().getPrefactor(refAt.getAtomicNo(),fitAt.getAtomicNo());
					if (atomOverlap>0.0) {
						double sim = refAt.getSimilarity(fitAt, fitPPDirectionalityVector);
						atomOverlap *= sim;
						totalOverlap += atomOverlap;
						double gradientPrefactor = atomOverlap*-2*refAt.getWidth()*fitAt.getWidth()/(refAt.getWidth()+fitAt.getWidth());
						double qder = qDers[j][0]*dx+qDers[j][1]*dy+qDers[j][2]*dz; 
						double rder = rDers[j][0]*dx+rDers[j][1]*dy+rDers[j][2]*dz; 
						double sder = sDers[j][0]*dx+sDers[j][1]*dy+sDers[j][2]*dz; 
						double uder = uDers[j][0]*dx+uDers[j][1]*dy+uDers[j][2]*dz; 

					    grad[0] += sim*gradientPrefactor*qder+atomOverlap*(qDers[j][0]*xi+qDers[j][1]*yi+qDers[j][2]*zi)/3.0;
					    grad[1] += sim*gradientPrefactor*rder+atomOverlap*(rDers[j][0]*xi+rDers[j][1]*yi+rDers[j][2]*zi)/3.0;
					    grad[2] += sim*gradientPrefactor*sder+atomOverlap*(sDers[j][0]*xi+sDers[j][1]*yi+sDers[j][2]*zi)/3.0;
					    grad[3] += sim*gradientPrefactor*uder+atomOverlap*(uDers[j][0]*xi+uDers[j][1]*yi+uDers[j][2]*zi)/3.0;
					    grad[4] += sim*gradientPrefactor*dx+atomOverlap*xi/3.0;
					    grad[5] += sim*gradientPrefactor*dy+atomOverlap*yi/3.0;
					    grad[6] += sim*gradientPrefactor*dz+atomOverlap*zi/3.0;
								    	
					    }
					}
				}

			grad[0] += PENALTY*(1-invnorm)*this.transform[0]; //penalty term to force quaternion into unity
			grad[1] += PENALTY*(1-invnorm)*this.transform[1];
			grad[2] += PENALTY*(1-invnorm)*this.transform[2];
			grad[3] += PENALTY*(1-invnorm)*this.transform[3];


			return (-1.0*totalOverlap+0.5*PENALTY*(norm-1)*(norm-1)); //the negative overlap is returned as the objective, since we minimize the objective in the optimization algorithm
			
		
		}


	public EvaluableOverlap clone() {
		return new EvaluableOverlap(this);
	}
		
}
