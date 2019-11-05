package com.actelion.research.chem.phesaflex;


import java.util.Arrays;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.phesa.AtomicGaussian;
import com.actelion.research.chem.phesa.Evaluable;
import com.actelion.research.chem.phesa.Gaussian3D;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.QuickMathCalculator;
import com.actelion.research.chem.phesa.PheSAAlignment;
import com.actelion.research.chem.phesa.pharmacophore.PPGaussian;



/**
 * @author J.Wahl, Oktober 2019
 * functionality for optimizing PheSA overlap (Pharmacophore+Shape) allowing for molecular flexibility
 */


public class EvaluableFlexibleOverlap implements Evaluable  {

	private static final double SCALE = 250;
	private static final double DELTA = -0.01;
	private double e0 = 0.0;
	private StereoMolecule fitMol;
	private PheSAAlignment shapeAlign;
	private boolean[] isHydrogen;
	private double[] v; //coordinates of the atoms
	private double[][] precalcPow;
	private double[] precalcExp;
    private double oAA;
    private double oAApp;
    private ForceFieldMMFF94 ff;

	public EvaluableFlexibleOverlap(PheSAAlignment shapeAlign, StereoMolecule fitMol, boolean[] isHydrogen,double[] v) {
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		ff = new ForceFieldMMFF94(fitMol, ForceFieldMMFF94.MMFF94SPLUS);
		this.shapeAlign = shapeAlign; 
		this.fitMol = fitMol;
		this.isHydrogen = isHydrogen;
		this.v = v;
		for(int a=0,i=0;i<fitMol.getAllAtoms();i++) {
			v[a++]=fitMol.getAtomX(i);
			v[a++]=fitMol.getAtomY(i);
			v[a++]=fitMol.getAtomZ(i);
		}

		this.oAA = this.getFGValueShapeSelf(new double[v.length], shapeAlign.getRefMolGauss(),true);
		this.oAApp = this.getFGValueShapeSelfPP(new double[v.length], shapeAlign.getRefMolGauss(),true);
		
	}
	
	public EvaluableFlexibleOverlap(EvaluableFlexibleOverlap e) {
		this.shapeAlign = e.shapeAlign;
		this.fitMol = e.fitMol;
		this.isHydrogen = e.isHydrogen;
		this.v = e.v;	
		this.precalcPow = e.precalcPow;
		this.precalcExp = e.precalcExp;
		this.ff = e.ff;

	}
	
	@Override
	public void setState(double[] v){
		//System.out.println(Arrays.toString(v));
		this.v=v;
		ff.setState(v);
	}
	
	public double[] getState(double[] v){
		for(int i=0;i<this.v.length;i++) {
			v[i] = this.v[i];
			
		}
		return v;
	}
	
	public void setE0(double e0 ) {
		this.e0 = e0;
	}
	
	public double[] getState() {
		return this.getState(new double[v.length]);
	}
	

	
	public PheSAAlignment getAlignment() {
		return this.shapeAlign;
	}
	

	
	@Override
	public double getFGValue(double[] grad) {
		double ePot = 0.0;
		double T = 0.0;
		double[] overlapGrad = new double[grad.length]; 
		double[] energyGrad = new double[grad.length]; 
		double[] selfOverlapGradFit = new double[grad.length];
		double[] overlapGradPP = new double[grad.length]; 
		double[] selfOverlapGradFitPP = new double[grad.length];
		double oBB = this.getFGValueShapeSelf(selfOverlapGradFit, shapeAlign.getMolGauss(),false);
		double oAB = this.getFGValueShape(overlapGrad);
		double oBBpp = this.getFGValueShapeSelfPP(selfOverlapGradFitPP, shapeAlign.getMolGauss(),false);
		double oABpp = this.getFGValueShapePP(overlapGradPP);
		ff.addGradient(energyGrad);
		ePot = ff.getTotalEnergy();
		double[] dOBB = selfOverlapGradFit;
		double[] dOAB = overlapGrad;
		double[] dOBB_dOAB = new double[grad.length];
		double[] dOBBpp = selfOverlapGradFitPP;
		double[] dOABpp = overlapGradPP;
		double[] dOBBpp_dOABpp = new double[grad.length];
		T = 0.5*(oAB/(oBB+oAA-oAB))+0.5*(oABpp/(oBBpp+oAApp-oABpp));
		System.out.println("new round");
		System.out.println(0.5*(oAB/(oBB+oAA-oAB)));
		System.out.println(0.5*(oABpp/(oBBpp+oAApp-oABpp)));
		double value = -SCALE*Math.exp(DELTA*(ePot-e0))*T + (ePot-e0);
		for(int i=0;i<grad.length;i++) {
			dOBB_dOAB[i] = dOBB[i]-dOAB[i];
			dOBBpp_dOABpp[i] = dOBBpp[i]-dOABpp[i];
		}
		double[] dT = new double[grad.length];
		for(int j=0;j<grad.length;j++) {
			dT[j] = dOAB[j]*(1/(oAA+oBB-oAB))-oAB*Math.pow(oAA+oBB-oAB,-2)*dOBB_dOAB[j] + 
					dOABpp[j]*(1/(oAApp+oBBpp-oABpp))-oAB*Math.pow(oAApp+oBBpp-oABpp,-2)*dOBBpp_dOABpp[j];
		}
		
		for(int k=0;k<grad.length;k++) {
			grad[k] = energyGrad[k]-(dT[k]*SCALE*Math.exp(DELTA*(ePot-e0))+T*SCALE*DELTA*Math.exp(DELTA*(ePot-e0))*energyGrad[k]);
		}


		return value;
		
		
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
	
	public double getFGValueShape(double[] grad) {
		
		MolecularVolume molGauss = shapeAlign.getMolGauss();

		MolecularVolume refMolGauss = shapeAlign.getRefMolGauss();

		for(int i=0;i<grad.length;i++) {
			grad[i] = 0;
		}

		/**
		 * derivative of ShapeOverlap with respect to the four elements of the quaternion and three elements of translation
		 * 
		 */ 
	    double totalOverlap = 0.0;
		for(AtomicGaussian refAt:refMolGauss.getAtomicGaussians()){
			double xi = refAt.getCenter().x;
			double yi = refAt.getCenter().y;
			double zi = refAt.getCenter().z;
			for(AtomicGaussian fitAt:molGauss.getAtomicGaussians()){
				int a = fitAt.getAtomId();
				double atomOverlap = 0.0;
				double xj = v[3*a];
				double yj = v[3*a+1];
				double zj = v[3*a+2];
				double dx = xi-xj;
				double dy = yi-yj;
				double dz = zi-zj;
				double Rij2 = dx*dx + dy*dy + dz*dz;
				double alphaSum = refAt.getWidth() + fitAt.getWidth();
				double gradientPrefactor=0.0;
				if(Rij2<Gaussian3D.DIST_CUTOFF) {
					atomOverlap = refAt.getHeight()*fitAt.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refAt.getWidth() * fitAt.getWidth()* Rij2)/alphaSum) *
							QuickMathCalculator.getInstance().getPrefactor(refAt.getAtomicNo(),fitAt.getAtomicNo());
					
					if (atomOverlap>0.0) {
						totalOverlap += atomOverlap;
						gradientPrefactor = atomOverlap*-2*refAt.getWidth()*fitAt.getWidth()/(refAt.getWidth()+fitAt.getWidth());
					}

				}
				grad[3*a] += (2*xj-2*xi)*gradientPrefactor;
				grad[3*a+1] += (2*yj-2*yi)*gradientPrefactor;
				grad[3*a+2] += (2*zj-2*zi)*gradientPrefactor;
				}

		
		}

		return totalOverlap; 
	
	}
	
	public double getFGValueShapePP(double[] grad) {
		
		MolecularVolume molGauss = shapeAlign.getMolGauss();

		MolecularVolume refMolGauss = shapeAlign.getRefMolGauss();

		for(int i=0;i<grad.length;i++) {
			grad[i] = 0;
		}

		/**
		 * derivative of ShapeOverlap with respect to the four elements of the quaternion and three elements of translation
		 * 
		 */ 
	    double totalOverlap = 0.0;
		for(PPGaussian refPP:refMolGauss.getPPGaussians()){
			double xi = refPP.getCenter().x;
			double yi = refPP.getCenter().y;
			double zi = refPP.getCenter().z;
			for(PPGaussian fitPP:molGauss.getPPGaussians()){
				int a = fitPP.getAtomId();
				double atomOverlap = 0.0;
				double xj = v[3*a];
				double yj = v[3*a+1];
				double zj = v[3*a+2];
				double dx = xi-xj;
				double dy = yi-yj;
				double dz = zi-zj;
				double Rij2 = dx*dx + dy*dy + dz*dz;
				double alphaSum = refPP.getWidth() + fitPP.getWidth();
				double gradientPrefactor=0.0;
				if(Rij2<Gaussian3D.DIST_CUTOFF) {
					atomOverlap = refPP.getHeight()*fitPP.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refPP.getWidth() * fitPP.getWidth()* Rij2)/alphaSum) *
							QuickMathCalculator.getInstance().getPrefactor(refPP.getAtomicNo(),fitPP.getAtomicNo());
					
					if (atomOverlap>0.0) {
						double sim = refPP.getSimilarity(fitPP);
						atomOverlap *= sim;
						totalOverlap += atomOverlap;
						gradientPrefactor = atomOverlap*-2*refPP.getWidth()*fitPP.getWidth()/(refPP.getWidth()+fitPP.getWidth());
						grad[3*a] += (2*xj-2*xi)*gradientPrefactor*sim;
						grad[3*a+1] += (2*yj-2*yi)*gradientPrefactor*sim;
						grad[3*a+2] += (2*zj-2*zi)*gradientPrefactor*sim;
						fitPP.getPharmacophorePoint().getDirectionalityDerivativeCartesian(grad, v, fitPP.getPharmacophorePoint().getDirectionality(), sim);					}

				}

				}

		
		}

		return totalOverlap; 
	
	}
	
	public double getFGValueShapeSelf(double[] grad, MolecularVolume molVol,boolean rigid) {
		double xi,yi,zi,xj,yj,zj;
		
		for(int i=0;i<grad.length;i++) {
			grad[i] = 0;
		}

		/**
		 * derivative of ShapeOverlap with respect to Cartesian coordinates
		 */ 
	    double totalOverlap = 0.0;
	    for(AtomicGaussian refAt:molVol.getAtomicGaussians()){
	    	int a = refAt.getAtomId();
			if(rigid) {
				xi = refAt.getCenter().x;
				yi = refAt.getCenter().y;
				zi = refAt.getCenter().z;
			}
			else {
				xi = v[3*a];
				yi = v[3*a+1];
				zi = v[3*a+2];
			}
			for(AtomicGaussian fitAt:molVol.getAtomicGaussians()){
				int b = fitAt.getAtomId();
				double atomOverlap = 0.0;

				if(rigid) {
					xj = fitAt.getCenter().x;
					yj = fitAt.getCenter().y;
					zj = fitAt.getCenter().z;
				}
				else {
					xj = v[3*b];
					yj = v[3*b+1];
					zj = v[3*b+2];
				}
				double dx = xi-xj;
				double dy = yi-yj;
				double dz = zi-zj;
				double Rij2 = dx*dx + dy*dy + dz*dz;
				double alphaSum = refAt.getWidth() + fitAt.getWidth();
				double gradientPrefactor = 0.0;
				
				if(Rij2<Gaussian3D.DIST_CUTOFF) {
					atomOverlap = refAt.getHeight()*fitAt.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refAt.getWidth() * fitAt.getWidth()* Rij2)/alphaSum) *
							QuickMathCalculator.getInstance().getPrefactor(refAt.getAtomicNo(),fitAt.getAtomicNo());
					
					if (atomOverlap>0.0) {
						totalOverlap += atomOverlap;
						gradientPrefactor = atomOverlap*-2*refAt.getWidth()*fitAt.getWidth()/(refAt.getWidth()+fitAt.getWidth());
					}

				}
					grad[3*b] += (2*xj-2*xi)*gradientPrefactor;
					grad[3*b+1] += (2*yj-2*yi)*gradientPrefactor;
					grad[3*b+2] += (2*zj-2*zi)*gradientPrefactor;
				}

			
		}


		return totalOverlap; 
		
	
	}
	
	public double getFGValueShapeSelfPP(double[] grad, MolecularVolume molVol,boolean rigid) {
		double xi,yi,zi,xj,yj,zj;
		
		for(int i=0;i<grad.length;i++) {
			grad[i] = 0;
		}

		/**
		 * derivative of ShapeOverlap with respect to Cartesian coordinates
		 */ 
	    double totalOverlap = 0.0;
	    for(PPGaussian refPP:molVol.getPPGaussians()){
	    	int a = refPP.getAtomId();
			if(rigid) {
				xi = refPP.getCenter().x;
				yi = refPP.getCenter().y;
				zi = refPP.getCenter().z;
			}
			else {
				xi = v[3*a];
				yi = v[3*a+1];
				zi = v[3*a+2];
			}
			for(PPGaussian fitPP:molVol.getPPGaussians()){
				int b = fitPP.getAtomId();
				double atomOverlap = 0.0;

				if(rigid) {
					xj = fitPP.getCenter().x;
					yj = fitPP.getCenter().y;
					zj = fitPP.getCenter().z;
				}
				else {
					xj = v[3*b];
					yj = v[3*b+1];
					zj = v[3*b+2];
				}
				double dx = xi-xj;
				double dy = yi-yj;
				double dz = zi-zj;
				double Rij2 = dx*dx + dy*dy + dz*dz;
				double alphaSum = fitPP.getWidth() + fitPP.getWidth();
				double gradientPrefactor = 0.0;
				
				if(Rij2<Gaussian3D.DIST_CUTOFF) {
					atomOverlap = refPP.getHeight()*fitPP.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refPP.getWidth() * fitPP.getWidth()* Rij2)/alphaSum) *
							QuickMathCalculator.getInstance().getPrefactor(refPP.getAtomicNo(),fitPP.getAtomicNo());
					
					if (atomOverlap>0.0) {
						double sim = refPP.getSimilarity(fitPP);
						atomOverlap *= sim;
						totalOverlap += atomOverlap;
						gradientPrefactor = atomOverlap*-2*refPP.getWidth()*fitPP.getWidth()/(refPP.getWidth()+fitPP.getWidth());
						grad[3*a] += (2*xj-2*xi)*gradientPrefactor*sim;
						grad[3*a+1] += (2*yj-2*yi)*gradientPrefactor*sim;
						grad[3*a+2] += (2*zj-2*zi)*gradientPrefactor*sim;
						fitPP.getPharmacophorePoint().getDirectionalityDerivativeCartesian(grad, v, fitPP.getPharmacophorePoint().getDirectionality(), sim);					}

				}

				}


			
		}


		return totalOverlap; 
		
	
	}
	
	

	@Override
	public EvaluableFlexibleOverlap clone() {
		return new EvaluableFlexibleOverlap(this);
	}
		
}
