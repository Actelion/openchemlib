package com.actelion.research.chem.phesaflex;


import java.util.Arrays;
import java.util.Map;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.optimization.Evaluable;
import com.actelion.research.chem.phesa.AtomicGaussian;
import com.actelion.research.chem.phesa.Gaussian3D;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.QuickMathCalculator;
import com.actelion.research.chem.phesa.VolumeGaussian;
import com.actelion.research.chem.phesa.PheSAAlignment;
import com.actelion.research.chem.phesa.pharmacophore.PPGaussian;



/**
 * @author JW, Oktober 2019
 * functionality for optimizing PheSA overlap (Pharmacophore+Shape) allowing for molecular flexibility
 */


public class EvaluableFlexibleOverlap implements Evaluable  {

	//private static final double SCALE = -250;
	//private static final double DELTA = -0.01;
	private static final double LAMBDA = 0.0625;
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
    private Map<String, Object> ffOptions;
    private double ppWeight;
    
	public EvaluableFlexibleOverlap(PheSAAlignment shapeAlign, StereoMolecule refMol, StereoMolecule fitMol, double ppWeight, boolean[] isHydrogen,double[] v, Map<String, Object> ffOptions) {
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		this.ffOptions = ffOptions;
		ff = new ForceFieldMMFF94(fitMol, ForceFieldMMFF94.MMFF94SPLUS, this.ffOptions);
		this.shapeAlign = shapeAlign;
		this.fitMol = fitMol;
		this.isHydrogen = isHydrogen;
		this.ppWeight = ppWeight;
		this.v = v;
		for(int i=0;i<fitMol.getAllAtoms();i++) {
			v[3*i]=fitMol.getAtomX(i);
			v[3*i+1]=fitMol.getAtomY(i);
			v[3*i+2]=fitMol.getAtomZ(i);
		}

		this.oAA = this.getFGValueShapeSelf(new double[3*refMol.getAllAtoms()], shapeAlign.getRefMolGauss(),true);
		this.oAApp = this.getFGValueSelfPP(new double[3*refMol.getAllAtoms()], shapeAlign.getRefMolGauss(),true);
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
		assert this.v.length==v.length;
		for(int i=0;i<v.length;i++) {
			this.v[i] = v[i];
		}
		ff.setState(v);
		for(int a=0,i=0;i<fitMol.getAllAtoms();i++) {
			fitMol.setAtomX(i,v[a++]);
			fitMol.setAtomY(i,v[a++]);
			fitMol.setAtomZ(i,v[a++]);
		}
		shapeAlign.getMolGauss().update(fitMol);
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
		double oBBpp = this.getFGValueSelfPP(selfOverlapGradFitPP, shapeAlign.getMolGauss(),false);
		double oABpp = this.getFGValuePP(overlapGradPP);
		ff.addGradient(energyGrad);
		ePot = ff.getTotalEnergy();
		double[] dOBB = selfOverlapGradFit;
		double[] dOAB = overlapGrad;
		double[] dOBB_dOAB = new double[grad.length];
		double[] dOBBpp = selfOverlapGradFitPP;
		double[] dOABpp = overlapGradPP;
		double[] dOBBpp_dOABpp = new double[grad.length];
		T = (1.0-ppWeight)*(oAB/(oBB+oAA-oAB))+ppWeight*(oABpp/(oBBpp+oAApp-oABpp));
		//double value = SCALE*Math.exp(DELTA*(ePot-e0))*T + (ePot-e0);
		double strainEnergy = ePot-e0;
		double strainPrefactor = strainEnergy < FlexibleShapeAlignment.ENERGY_CUTOFF ? 0.0 : strainEnergy-FlexibleShapeAlignment.ENERGY_CUTOFF;
		double value = -T + LAMBDA*strainPrefactor*strainPrefactor;

		for(int i=0;i<grad.length;i++) {
			dOBB_dOAB[i] = dOBB[i]-dOAB[i];
			dOBBpp_dOABpp[i] = dOBBpp[i]-dOABpp[i];
		}
		double[] dT = new double[grad.length];
		for(int j=0;j<grad.length;j++) {
			dT[j] = (1.0-ppWeight)*dOAB[j]*(1/(oAA+oBB-oAB))-(1.0-ppWeight)*oAB*Math.pow(oAA+oBB-oAB,-2)*dOBB_dOAB[j] + 
					ppWeight*dOABpp[j]*(1/(oAApp+oBBpp-oABpp))-ppWeight*oAB*Math.pow(oAApp+oBBpp-oABpp,-2)*dOBBpp_dOABpp[j];
		}
		for(int k=0;k<grad.length;k++) {

			grad[k] = -dT[k] + strainPrefactor*2*LAMBDA*energyGrad[k];
		}


		return value;
		
		
	}

	/**
	 * calculates the gradient of the overlap function with respect to the cartesian coordinates of the atoms
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
		
		for(VolumeGaussian refVG:refMolGauss.getVolumeGaussians()){
			double xi = refVG.getCenter().x;
			double yi = refVG.getCenter().y;
			double zi = refVG.getCenter().z;
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
				double alphaSum = refVG.getWidth() + fitAt.getWidth();
				double gradientPrefactor=0.0;
				if(Rij2<Gaussian3D.DIST_CUTOFF) {
					atomOverlap = refVG.getRole()*refVG.getHeight()*fitAt.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refVG.getWidth() * fitAt.getWidth()* Rij2)/alphaSum) *
							QuickMathCalculator.getInstance().getPrefactor(refVG.getAtomicNo(),fitAt.getAtomicNo());
					
					if (Math.abs(atomOverlap)>0.0) {
						totalOverlap += atomOverlap;
						gradientPrefactor = atomOverlap*-2*refVG.getWidth()*fitAt.getWidth()/(refVG.getWidth()+fitAt.getWidth());
					}

				}
				grad[3*a] += (2*xj-2*xi)*gradientPrefactor;
				grad[3*a+1] += (2*yj-2*yi)*gradientPrefactor;
				grad[3*a+2] += (2*zj-2*zi)*gradientPrefactor;
				}

		
		}


		return totalOverlap; 
	
	}
	
	public double getFGValuePP(double[] grad) {
		
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
					atomOverlap = refPP.getWeight()*refPP.getHeight()*fitPP.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refPP.getWidth() * fitPP.getWidth()* Rij2)/alphaSum) *
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
	
	public double getFGValueShapeSelf(double[] grad, MolecularVolume molGauss,boolean rigid) {
		
		for(int i=0;i<grad.length;i++) {
			grad[i] = 0;
		}

		/**
		 * derivative of ShapeOverlap with respect to Cartesian coordinates
		 */ 
	    double totalOverlap = 0.0;
	    for(AtomicGaussian refAt:molGauss.getAtomicGaussians()){
			for(AtomicGaussian  fitAt:molGauss.getAtomicGaussians()){
				totalOverlap+=getGradientContribution(refAt,fitAt,grad,rigid);
			}
			for(VolumeGaussian fitAt:molGauss.getVolumeGaussians()){
				totalOverlap+=getGradientContribution(refAt,fitAt,grad,rigid);
			}
	    }
	    for(AtomicGaussian refAt:molGauss.getAtomicGaussians()){
	    	for(VolumeGaussian fitAt:molGauss.getVolumeGaussians()){
	    		totalOverlap+=getGradientContribution(refAt,fitAt,grad,rigid);
	    	}
	    }

		return totalOverlap; 
		
	
	}
	
	public double getGradientContribution(Gaussian3D refAt, Gaussian3D fitAt, double[] grad,boolean rigid) {
		double xi,yi,zi,xj,yj,zj;
		int b = fitAt.getAtomId();
		double atomOverlap = 0.0;
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
		double d = 1.0;
		if(refAt instanceof VolumeGaussian) 
			d*= ((VolumeGaussian)refAt).getRole();
		if(fitAt instanceof VolumeGaussian) 
			d*= ((VolumeGaussian)fitAt).getRole();
		
		if(Rij2<Gaussian3D.DIST_CUTOFF) {
			atomOverlap = d*refAt.getHeight()*fitAt.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refAt.getWidth() * fitAt.getWidth()* Rij2)/alphaSum) *
					QuickMathCalculator.getInstance().getPrefactor(refAt.getAtomicNo(),fitAt.getAtomicNo());
			
			if (atomOverlap>0.0) {
				gradientPrefactor = atomOverlap*-2*refAt.getWidth()*fitAt.getWidth()/(refAt.getWidth()+fitAt.getWidth());
			}

		}
			grad[3*b] += (2*xj-2*xi)*gradientPrefactor;
			grad[3*b+1] += (2*yj-2*yi)*gradientPrefactor;
			grad[3*b+2] += (2*zj-2*zi)*gradientPrefactor;
			
			return atomOverlap;
		}
	

	
	
	public double getFGValueSelfPP(double[] grad, MolecularVolume molVol,boolean rigid) {
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
					atomOverlap = refPP.getWeight()*refPP.getHeight()*fitPP.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refPP.getWidth() * fitPP.getWidth()* Rij2)/alphaSum) *
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
