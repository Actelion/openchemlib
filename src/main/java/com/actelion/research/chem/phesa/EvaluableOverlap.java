package com.actelion.research.chem.phesa;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.alignment3d.transformation.ExponentialMap;
import com.actelion.research.chem.alignment3d.transformation.Quaternion;
import com.actelion.research.chem.alignment3d.transformation.RotationDerivatives;
import com.actelion.research.chem.optimization.Evaluable;
import com.actelion.research.chem.phesa.pharmacophore.pp.PPGaussian;



/**
 * @author J.Wahl, February 2018
 * describes the overlap function, the state (relative orientation, translation of the two molecules)
 * returns the objective function and the gradient of the overlap with respect to translation and rotation
 * accessed by the optimization algorithm
 */


public class EvaluableOverlap implements Evaluable  {

	private static final int PENALTY = 80; 
	
	private double ppWeight;

	private PheSAAlignment shapeAlign;
	private double[] transform;
	private Coordinates[] cachedCoords; //coords before rotation + translation is applied, but COM at coordinate origin: Xi-p  
	private Coordinates[] cachedCoordsPP;
	private Coordinates origCOM;
    private double [][] dv0At;
    private double [][] dv1At;
    private double [][] dv2At;
    private double [][] dv0PP;
    private double [][] dv1PP;
    private double [][] dv2PP;
    private double[][] results;
    private Coordinates[] fitAtGaussModCoords;
    private Coordinates[] fitPPGaussModCoords;
    private Coordinates[] fitPPDirectionalityMod;

    
    public EvaluableOverlap(PheSAAlignment shapeAlign, double[] transform) {
    	this(shapeAlign, transform, 0.5);
    }
    
    public EvaluableOverlap(PheSAAlignment shapeAlign, double[] transform, double ppWeight) {
    	this.ppWeight = ppWeight;
		this.shapeAlign = shapeAlign; 
		this.transform = transform;
	    this.fitAtGaussModCoords = new Coordinates[shapeAlign.getMolGauss().getAtomicGaussians().size()];
	    this.fitPPGaussModCoords = new Coordinates[shapeAlign.getMolGauss().getPPGaussians().size()];
	    this.dv0At = new double[fitAtGaussModCoords.length][3];
	    this.dv1At = new double[fitAtGaussModCoords.length][3];
	    this.dv2At = new double[fitAtGaussModCoords.length][3];
	    this.dv0PP = new double[fitPPGaussModCoords.length][3];
	    this.dv1PP = new double[fitPPGaussModCoords.length][3];
	    this.dv2PP = new double[fitPPGaussModCoords.length][3];
	    this.fitPPDirectionalityMod = new Coordinates[shapeAlign.getMolGauss().getPPGaussians().size()];
		this.results = new double[shapeAlign.getRefMolGauss().getAtomicGaussians().size()][shapeAlign.getRefMolGauss().getAtomicGaussians().size()];
		cachedCoords = new Coordinates[shapeAlign.getMolGauss().getAtomicGaussians().size()];
		cachedCoordsPP = new Coordinates[shapeAlign.getMolGauss().getPPGaussians().size()];
		origCOM = new Coordinates();
		for(int i=0;i<shapeAlign.getMolGauss().getAtomicGaussians().size();i++) {
			cachedCoords[i] = shapeAlign.getMolGauss().getAtomicGaussians().get(i).center;
		}
		for(int i=0;i<shapeAlign.getMolGauss().getPPGaussians().size();i++) {
			cachedCoordsPP[i] = shapeAlign.getMolGauss().getPPGaussians().get(i).center;
		}
		for(Coordinates coords : cachedCoords){
			origCOM.add(coords);
		}
		origCOM.scale(1.0/cachedCoords.length);
		for(Coordinates coords : cachedCoords) {
			coords.sub(origCOM);
		}
		for(Coordinates coords : cachedCoordsPP) {
			coords.sub(origCOM);
		}
		

		
	}
	
	public EvaluableOverlap(EvaluableOverlap e) {
		this.shapeAlign = e.shapeAlign;
		this.transform = e.transform;	
	    this.dv0At = e.dv0At;
	    this.dv1At = e.dv1At;
	    this.dv2At = e.dv2At;
	    this.dv0PP = e.dv0PP;
	    this.dv0PP = e.dv1PP;
	    this.dv0PP = e.dv2PP;
		this.fitAtGaussModCoords = e.fitAtGaussModCoords;
		this.fitPPGaussModCoords = e.fitPPGaussModCoords;
		this.results = e.results;

	}
	
	private void getTransformedCoordinates(Coordinates[] modCoords,List<? extends Gaussian3D> fitMolGauss) {
		 ExponentialMap em = new ExponentialMap(transform[0],transform[1],transform[2]);
		 Quaternion q = em.toQuaternion();
		 double[][] rotMatrix = q.getRotMatrix().getArray();
		 for(int k=0;k<fitMolGauss.size();k++) {
		    Coordinates center = new Coordinates(fitMolGauss.get(k).center);
		    center.sub(origCOM);
		    center.rotate(rotMatrix);
		    center.add(origCOM);
		    center.x += transform[3];
		    center.y += transform[4];
		    center.z += transform[5]; 
		    modCoords[k] = center;
		  }
	}
	
	@Override
	public void setState(double[] transform){
		assert this.transform.length==transform.length;
		for(int i=0;i<transform.length;i++) {
			this.transform[i] = transform[i];
		}
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
		ShapeVolume refMolGauss = shapeAlign.getRefMolGauss();
		ShapeVolume fitMolGauss = shapeAlign.getMolGauss();
		double value = 0.0;
		double[] atomGrad = new double[grad.length];
		List<VolumeGaussian> volumeGaussians = new ArrayList<>();
		if(refMolGauss instanceof MolecularVolume)
			volumeGaussians = ((MolecularVolume)refMolGauss).getVolumeGaussians();
		value += (1.0-ppWeight)*this.getFGValueOverlap(atomGrad,refMolGauss.getAtomicGaussians(),fitMolGauss.getAtomicGaussians(),volumeGaussians,
						dv0At,dv1At,dv2At,fitAtGaussModCoords);
			
		
		double[] ppGrad = new double[grad.length];
		value += ppWeight*this.getFGValueOverlapPP(ppGrad,refMolGauss.getPPGaussians(),fitMolGauss.getPPGaussians(),
						dv0PP,dv1PP,dv2PP,fitPPGaussModCoords,fitPPDirectionalityMod);

		for(int i=0;i<grad.length;i++) 
			grad[i] = (1.0-ppWeight)*atomGrad[i]+ ppWeight*ppGrad[i];
				

		return value;
		
		
	}
	
	
	private void getEMapGradient(double[][] dRdv0, double[][] dRdv1, double[][] dRdv2, Coordinates[] cachedCoords) {

		    /**
		     * we first calculate the partial derivatives with respect to the three elements of the exponential map
		     * the final gradient has 7 elements, the first four elements are the gradients for the quaternion (rotation),
		     * the last three elements are for the translation
		     */
		
		double[] v = new double[] {transform[0],transform[1],transform[2]}; //exponential map
		RotationDerivatives rotationDerivatives = new RotationDerivatives(v);
		 for(int a=0;a<cachedCoords.length;a++){
			Coordinates xi = cachedCoords[a];
			double[][] dRdvi = new double[3][3];
			rotationDerivatives.dRdv(0, dRdvi);
			Coordinates dRij_dv0 = xi.rotateC(dRdvi);
			rotationDerivatives.dRdv(1, dRdvi);
			Coordinates dRij_dv1 = xi.rotateC(dRdvi);
			rotationDerivatives.dRdv(2, dRdvi);
			Coordinates dRij_dv2 = xi.rotateC(dRdvi);
	
			dRdv0[a][0] = dRij_dv0.x;
			dRdv0[a][1] = dRij_dv0.y;
			dRdv0[a][2] = dRij_dv0.z;
			
			dRdv1[a][0] = dRij_dv1.x;
			dRdv1[a][1] = dRij_dv1.y;
			dRdv1[a][2] = dRij_dv1.z;
			
			dRdv2[a][0] = dRij_dv2.x;
			dRdv2[a][1] = dRij_dv2.y;
			dRdv2[a][2] = dRij_dv2.z;
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
	
	private double getFGValueOverlap(double[] grad,List<AtomicGaussian> refMolGauss,List<AtomicGaussian> fitMolGauss,List<VolumeGaussian> volGaussians,
			double[][] dRdv0, double[][] dRdv1, double[][] dRdv2, Coordinates[] fitGaussModCoords) {


	    /**
	     * we first calculate the partial derivatives with respect to the four elements of the quaternion q,r,s,u
	     * the final gradient has 7 elements, the first four elements are the gradients for the quaternion (rotation),
	     * the last three elements are for the translation
	     */
	    
	    getTransformedCoordinates(fitGaussModCoords, fitMolGauss);

	    this.getEMapGradient(dRdv0, dRdv1, dRdv2,cachedCoords);


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
					double dv0 = dRdv0[j][0]*dx+dRdv0[j][1]*dy+dRdv0[j][2]*dz; 
					double dv1 = dRdv1[j][0]*dx+dRdv1[j][1]*dy+dRdv1[j][2]*dz; 
					double dv2 = dRdv2[j][0]*dx+dRdv2[j][1]*dy+dRdv2[j][2]*dz; 

					grad[0] += gradientPrefactor*dv0;
					grad[1] += gradientPrefactor*dv1;
					grad[2] += gradientPrefactor*dv2;
					grad[3] += gradientPrefactor*dx;
					grad[4] += gradientPrefactor*dy;
					grad[5] += gradientPrefactor*dz;
				    }



				}
		}
			for(int k=0; k<volGaussians.size();k++){
				VolumeGaussian refVol = volGaussians.get(k);
				for(int j=0; j<fitMolGauss.size();j++){
					double atomOverlap = 0.0;
					Gaussian3D fitAt = fitMolGauss.get(j);
					fitCenterModCoord = fitGaussModCoords[j];
					double alphaSum = refVol.getWidth() + fitAt.getWidth();

					double dx = refVol.getCenter().x-fitCenterModCoord.x;
					double dy = refVol.getCenter().y-fitCenterModCoord.y;
					double dz = refVol.getCenter().z-fitCenterModCoord.z;
					double Rij2 = dx*dx + dy*dy + dz*dz;

					if(Rij2>=Gaussian3D.DIST_CUTOFF) 
						continue;
					atomOverlap = refVol.getRole()*refVol.getHeight()*fitAt.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refVol.getWidth() * fitAt.getWidth()* Rij2)/alphaSum) *
								QuickMathCalculator.getInstance().getPrefactor(refVol.getAtomicNo(),fitAt.getAtomicNo());
					
					
					if (Math.abs(atomOverlap)>0.0) {
						totalOverlap += atomOverlap;
						double gradientPrefactor = atomOverlap*-2*refVol.getWidth()*fitAt.getWidth()/(refVol.getWidth()+fitAt.getWidth());
						double dv0 = dRdv0[j][0]*dx+dRdv0[j][1]*dy+dRdv0[j][2]*dz; 
						double dv1 = dRdv1[j][0]*dx+dRdv1[j][1]*dy+dRdv1[j][2]*dz; 
						double dv2 = dRdv2[j][0]*dx+dRdv2[j][1]*dy+dRdv2[j][2]*dz; 

						grad[0] += gradientPrefactor*dv0;
						grad[1] += gradientPrefactor*dv1;
						grad[2] += gradientPrefactor*dv2;
						grad[3] += gradientPrefactor*dx;
						grad[4] += gradientPrefactor*dy;
						grad[5] += gradientPrefactor*dz;
					    }


					}
				
			}





		return (-1.0*totalOverlap); //the negative overlap is returned as the objective, since we minimize the objective in the optimization algorithm

	
	}
	    
	    
	   private double getFGValueOverlapPP(double[] grad, List<PPGaussian> refMolGauss,List<PPGaussian> fitMolGauss, double[][] dRdv0, double[][] dRdv1, double[][] dRdv2, Coordinates[] fitGaussModCoords, Coordinates[] fitPPDirectionalityMod) {
		   	ExponentialMap eMap = new ExponentialMap(transform[0],transform[1],transform[2]);
		    
		    double[][] rotMatrix = eMap.toQuaternion().getRotMatrix().getArray();


		    for(int k=0;k<fitMolGauss.size();k++) {
		    	fitPPDirectionalityMod[k] = fitMolGauss.get(k).getRotatedDirectionality(rotMatrix, 1.0);
		    }
		    getTransformedCoordinates(fitGaussModCoords,fitMolGauss);

		    this.getEMapGradient(dRdv0, dRdv1, dRdv2,cachedCoordsPP);


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
					
					atomOverlap = refAt.getWeight()*refAt.getHeight()*fitAt.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refAt.getWidth() * fitAt.getWidth()* Rij2)/alphaSum) *
							QuickMathCalculator.getInstance().getPrefactor(refAt.getAtomicNo(),fitAt.getAtomicNo());
					if (atomOverlap>0.0) {
						double sim = refAt.getSimilarity(fitAt, fitPPDirectionalityVector);
						atomOverlap *= sim;
						totalOverlap += atomOverlap;
						double gradientPrefactor = atomOverlap*-2*refAt.getWidth()*fitAt.getWidth()/(refAt.getWidth()+fitAt.getWidth());
						double dv0 = dRdv0[j][0]*dx+dRdv0[j][1]*dy+dRdv0[j][2]*dz; 
						double dv1 = dRdv1[j][0]*dx+dRdv1[j][1]*dy+dRdv1[j][2]*dz; 
						double dv2 = dRdv2[j][0]*dx+dRdv2[j][1]*dy+dRdv2[j][2]*dz;   

					    grad[0] += sim*gradientPrefactor*dv0+atomOverlap*(dRdv0[j][0]*xi+dRdv0[j][1]*yi+dRdv0[j][2]*zi)/3.0;
					    grad[1] += sim*gradientPrefactor*dv1+atomOverlap*(dRdv1[j][0]*xi+dRdv1[j][1]*yi+dRdv1[j][2]*zi)/3.0;
					    grad[2] += sim*gradientPrefactor*dv2+atomOverlap*(dRdv2[j][0]*xi+dRdv2[j][1]*yi+dRdv2[j][2]*zi)/3.0;
					    grad[3] += sim*gradientPrefactor*dx+atomOverlap*xi/3.0;
					    grad[4] += sim*gradientPrefactor*dy+atomOverlap*yi/3.0;
					    grad[5] += sim*gradientPrefactor*dz+atomOverlap*zi/3.0;
								    	
					    }
					}
				}

			return (-1.0*totalOverlap); //the negative overlap is returned as the objective, since we minimize the objective in the optimization algorithm
			
		
		}


	public EvaluableOverlap clone() {
		return new EvaluableOverlap(this);
	}
		
}
