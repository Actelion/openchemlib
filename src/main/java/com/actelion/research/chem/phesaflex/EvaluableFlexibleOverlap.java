package com.actelion.research.chem.phesaflex;


import java.util.Map;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer.PheSASetting;
import com.actelion.research.chem.alignment3d.transformation.ExponentialMap;
import com.actelion.research.chem.alignment3d.transformation.Quaternion;
import com.actelion.research.chem.alignment3d.transformation.Rotation;
import com.actelion.research.chem.alignment3d.transformation.RotationDerivatives;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;
import com.actelion.research.chem.alignment3d.transformation.Translation;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.optimization.Evaluable;
import com.actelion.research.chem.phesa.AtomicGaussian;
import com.actelion.research.chem.phesa.Gaussian3D;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.QuickMathCalculator;
import com.actelion.research.chem.phesa.ShapeVolume;
import com.actelion.research.chem.phesa.VolumeGaussian;
import com.actelion.research.chem.phesa.pharmacophore.pp.PPGaussian;
import com.actelion.research.chem.phesa.PheSAAlignment;



/**
 * @author JW, Oktober 2019
 * functionality for optimizing PheSA overlap (Pharmacophore+Shape) allowing for molecular flexibility
 */


public class EvaluableFlexibleOverlap implements Evaluable  {

	//private static final double SCALE = -250;
	//private static final double DELTA = -0.01;
	private static final double LAMBDA = 0.0625;
	private double e0 = 0.0;
	private Conformer fitConf;
	private StereoMolecule refMol;
	private PheSAAlignment shapeAlign;
	private boolean[] isHydrogen;
	private double[] v; // internal coordinates of the atoms
	private double[][] precalcPow;
	private double[] precalcExp;
    private double oAA;
    private double oAApp;
    private ForceFieldMMFF94 ff;
    private Map<String, Object> ffOptions;
    private PheSASetting settings;
	private Coordinates[] origCoords;
	private Coordinates[] cachedCoords; //used for gradient calculation: ligand coordinates with adjusted dihedral angles, but before rotation and translation
	private double[][] dRdvi1;
	private double[][] dRdvi2;
	private double[][] dRdvi3;
	private Coordinates origCOM;
	private BondRotationHelper torsionHelper;
    
	public EvaluableFlexibleOverlap(PheSAAlignment shapeAlign, StereoMolecule refMol, Conformer fitConf, BondRotationHelper torsionHelper,
			PheSASetting settings, boolean[] isHydrogen, Map<String, Object> ffOptions) {
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		this.ffOptions = ffOptions;
		this.shapeAlign = shapeAlign;
		this.fitConf = fitConf;
		this.isHydrogen = isHydrogen;
		this.settings = settings;
		this.refMol = refMol;
		this.torsionHelper = torsionHelper;
		init();
	}
	
	private void init() {
		ff = new ForceFieldMMFF94(fitConf.getMolecule(), ForceFieldMMFF94.MMFF94SPLUS, this.ffOptions);
		setInitialState();
		origCoords = new Coordinates[fitConf.getMolecule().getAllAtoms()];
		cachedCoords = new Coordinates[fitConf.getMolecule().getAllAtoms()];
		origCOM = new Coordinates();
		for(int a=0;a<fitConf.getMolecule().getAllAtoms();a++) {
			origCoords[a] = new Coordinates(fitConf.getMolecule().getCoordinates(a));
			cachedCoords[a] = new Coordinates(fitConf.getMolecule().getCoordinates(a));
			origCOM.add(cachedCoords[a]);
		}
		origCOM.scale(1.0/cachedCoords.length);
		dRdvi1 = new double[3][3];
		dRdvi2 = new double[3][3];
		dRdvi3 = new double[3][3];
		this.oAA = this.getFGValueShapeSelf(new double[3*refMol.getAllAtoms()], shapeAlign.getRefMolGauss(),true);
		this.oAApp = this.getFGValueSelfPP(shapeAlign.getRefMolGauss(),true);
	}
	
	public void setInitialState() {
		int elements = 3+3+torsionHelper.getRotatableBonds().length; //3 translational, 3 rotational, 3 torsion
		v = new double[elements];
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 0.0;
		Quaternion quat = new Quaternion(1.0,0.0,0.0,0.0);
		ExponentialMap emap = new ExponentialMap(quat);
		v[3] = emap.getP().x;
		v[4] = emap.getP().y;
		v[5] = emap.getP().z;
		for(int b=0;b<torsionHelper.getRotatableBonds().length;b++) {
			int[] atoms = torsionHelper.getTorsionAtoms()[b];
			v[6+b] = TorsionDB.calculateTorsionExtended(fitConf, atoms);
		}
	}
	
	private void resetLigCoordinates() {
		for(int a=0;a<fitConf.getMolecule().getAllAtoms();a++) {
			fitConf.setX(a, origCoords[a].x);
			fitConf.setY(a, origCoords[a].y);
			fitConf.setZ(a, origCoords[a].z);
		}
	}
	
	public EvaluableFlexibleOverlap(EvaluableFlexibleOverlap e) {
		this.shapeAlign = e.shapeAlign;
		this.fitConf = e.fitConf;
		this.isHydrogen = e.isHydrogen;
		this.v = e.v;	
		this.precalcPow = e.precalcPow;
		this.precalcExp = e.precalcExp;
		init();
	}
	

	
	@Override
	public void setState(double[] v){
		assert this.v.length==v.length;
		for(int i=0;i<v.length;i++) {
			if(i>5) { //torsions
				if(v[i]>Math.PI) {
					v[i] -= 2*Math.PI;
				}
			}
			this.v[i] = v[i];
	
		}
		updateLigandCoordinates();
		shapeAlign.getMolGauss().update(fitConf);
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
	
	public void updateLigandCoordinates() {
		resetLigCoordinates();
		//1. update dihedral angles
		//2. translate COM of ligand to rotation center
		//3. apply rotation 
		//4. translate back
		//5. apply translation 
		updateDihedralAngles();
		for(int a=0;a<fitConf.getMolecule().getAllAtoms();a++) {
			cachedCoords[a] = new Coordinates(fitConf.getCoordinates(a));
		}
	
		ExponentialMap eMap = new ExponentialMap(v[3],v[4],v[5]);
		Quaternion q = eMap.toQuaternion();
		Translation trans1 = new Translation(origCOM.scaleC(-1.0));
		Translation trans2 = new Translation(origCOM);
		Rotation rot = new Rotation(q.getRotMatrix().getArray());
		Translation t = new Translation(v[0],v[1],v[2]);
		TransformationSequence transformation = new TransformationSequence();
		transformation.addTransformation(trans1);
		transformation.addTransformation(rot);
		transformation.addTransformation(trans2);
		transformation.addTransformation(t);
		transformation.apply(fitConf);
	}
	
	private void updateDihedralAngles() {
		for(int b=0;b<torsionHelper.getRotatableBonds().length;b++) {
			double targetTorsion = v[6+b];
			int[] atoms = torsionHelper.getTorsionAtoms()[b];
			double currentTorsion = TorsionDB.calculateTorsionExtended(fitConf, atoms);
			double deltaTorsion = targetTorsion - currentTorsion;
			torsionHelper.rotateAroundBond(b, deltaTorsion,fitConf,false);
		}
	}
	

	/**
	 * inspired by: 10.1021/acs.jcim.7b00618 (Reflex3D)
	 * objective function:
	 * f = -T + lambda *strain^2            
	 * A is the fixed reference molecule, whereas B is flexible aligned on A3
	 *                                                             
	 * where T is the shape Tanimoto : T = O_AB/(O_BB+O_AA-O_AB): with O_AB being the shape overlap between molecules A and B
	 *                                                             and O_AA and O_BB the self overlaps
	 * next we calculate the derivate of the shape tanimoto with respect to the cartesian coordinates:
	 * dT/dx
	 * we use the product rule: T = u*v
	 * dT/dx = du/dx*v+u*dv/dx
	 * u = O_AB and v = (O_BB+O_AA-O_AB)^-1 
	 * du/dx = dO_AB/dx
	 * for v we use the chain rule:
	 * dv/dx = -(dO_BB/dx-dO_AB/dx)*(O_AA+O_BB-O_AB)^-2        
	 * dT/dx=dO_AB/dx * 1/(O_AA+O_BB-O_AB)-O_AB*(dO_BB/dx-dO_AB/dx)*(O_AA+O_BB-O_AB)^-2 
	 * and the strain term:
	 * s = lambda*strain^2 with strain = e-c (e: potential energy, c: energy cutoff)                                                    
	 * using the chain rule:
	 * ds/dx = 2*lambda*(e-c)*de/dx
	 */
	@Override
	public double getFGValue(double[] grad) {
		// first calculate the gradient with respect to the coordinates
		double[] coordGrad = new double[fitConf.getMolecule().getAllAtoms()*3];
		double ePot = 0.0;
		double T = 0.0;
		double[] overlapGrad = new double[coordGrad.length]; 
		double[] energyGrad = new double[coordGrad.length]; 
		double[] selfOverlapGradFit = new double[coordGrad.length];
		double oBB = this.getFGValueShapeSelf(selfOverlapGradFit, shapeAlign.getMolGauss(),false);
		double oAB = this.getFGValueShape(overlapGrad);
		double oBBpp = this.getFGValueSelfPP(shapeAlign.getMolGauss(),false);
		double oABpp = this.getFGValuePP();
		ff.setState(getCartState());
		ff.addGradient(energyGrad);
		ePot = ff.getTotalEnergy();
		double[] dOBB = selfOverlapGradFit;
		double[] dOAB = overlapGrad;
		double[] dOBB_dOAB = new double[coordGrad.length];
		T = (1.0-settings.getPpWeight())*(oAB/(oBB+oAA-oAB))+settings.getPpWeight()*(oABpp/(oBBpp+oAApp-oABpp));
		//double value = SCALE*Math.exp(DELTA*(ePot-e0))*T + (ePot-e0);
		double strainEnergy = ePot-e0;
		double strainPrefactor = strainEnergy < FlexibleShapeAlignment.ENERGY_CUTOFF ? 0.0 : strainEnergy-FlexibleShapeAlignment.ENERGY_CUTOFF;
		double value = -T + LAMBDA*strainPrefactor*strainPrefactor;
		for(int i=0;i<coordGrad.length;i++) {
			dOBB_dOAB[i] = dOBB[i]-dOAB[i];
		}
		double[] dT = new double[coordGrad.length];
		for(int j=0;j<coordGrad.length;j++) {
			dT[j] = (1.0-settings.getPpWeight())*dOAB[j]*(1/(oAA+oBB-oAB))-(1.0-settings.getPpWeight())*oAB*Math.pow(oAA+oBB-oAB,-2)*dOBB_dOAB[j];// + 
					//ppWeight*dOABpp[j]*(1/(oAApp+oBBpp-oABpp))-ppWeight*oAB*Math.pow(oAApp+oBBpp-oABpp,-2)*dOBBpp_dOABpp[j];
		}
		for(int k=0;k<coordGrad.length;k++) {

			coordGrad[k] = -dT[k] + strainPrefactor*2*LAMBDA*energyGrad[k];
		}
		//to inner coordinates
		//1. with respect to translational DOG
		for(int a=0;a<fitConf.getMolecule().getAllAtoms();a++) {
			grad[0] += coordGrad[3*a]; 
			grad[1] += coordGrad[3*a+1]; 
			grad[2] += coordGrad[3*a+2]; 
		}
		//2. orientational 
		//with respect to vector of exponential mapping p
		// dE/dpj = Tj*vi'*dE/dx
		//vi': atomic position (after adjustment of torsion values)
		double[] p = new double[] {v[3],v[4],v[5]};
		RotationDerivatives transformDerivatives = new RotationDerivatives(p);
		transformDerivatives.dRdv(0, dRdvi1);
		transformDerivatives.dRdv(1, dRdvi2);
		transformDerivatives.dRdv(2, dRdvi3);
		for(int a=0;a<fitConf.getMolecule().getAllAtoms();a++) {
			Coordinates vi = cachedCoords[a];
			Coordinates Tj_vi = vi.rotateC(dRdvi1);
			grad[3] += coordGrad[3*a]*Tj_vi.x+coordGrad[3*a+1]*Tj_vi.y+coordGrad[3*a+2]*Tj_vi.z;
			Tj_vi = vi.rotateC(dRdvi2);
			grad[4] += coordGrad[3*a]*Tj_vi.x+coordGrad[3*a+1]*Tj_vi.y+coordGrad[3*a+2]*Tj_vi.z;
			Tj_vi = vi.rotateC(dRdvi3);
			grad[5] += coordGrad[3*a]*Tj_vi.x+coordGrad[3*a+1]*Tj_vi.y+coordGrad[3*a+2]*Tj_vi.z;
		}
		//3. torsional gradient
		for(int b=0;b<torsionHelper.getRotatableBonds().length;b++) {
			int[] rotatedAtoms = torsionHelper.getSmallerSideAtomLists()[b];
			int j = torsionHelper.getRotationCenters()[b];
			int k = torsionHelper.getTorsionAtoms()[b][1] == j ? torsionHelper.getTorsionAtoms()[b][2] : torsionHelper.getTorsionAtoms()[b][1];
			Coordinates v1 = fitConf.getCoordinates(k).subC(fitConf.getCoordinates(j));

			for(int i : rotatedAtoms) {
				Coordinates v2 = 
						fitConf.getCoordinates(i).subC(fitConf.getCoordinates(j));
				Coordinates dx_dphi = v1.cross(v2);
				grad[6+b] += dx_dphi.x*coordGrad[3*i] + dx_dphi.y*coordGrad[3*i+1] + 
						dx_dphi.z*coordGrad[3*i+2];
			}
			
				//state[5+b+1] = TorsionDB.calculateTorsionExtended(ligConf, atoms);
		}
		return value;
		
		
	}

	/**
	 * calculates the gradient of the overlap function with respect to the cartesian coordinates of the atoms
	 * 
	 * 
	 * the overlap between molecules A and B is expressed as a sum over the atom-atom overlaps Vij
	 *  
     *                     pi           3/2                 alpha_i*alpha_j*Rij^2
     * Vij = p_i*p_j*(---------------)        * exp( - -------------------- )      
     *               alpha_i + alpha_j                     alpha_i + alpha_j 
     * 
     *     
     *              alpha_i*alpha_j*Rij^2
     *     = a*exp(- --------------------------)
     *                   alpha_i + alpha_j 
     * 
     * Rij^2 = (xi-xj)^2                  
     * we therefore need the derivative of Vij with respect to the atomic coordinates of atom j (molecule B)
     * 
     * we use the chain rule: 
     * 
     * d[e^u(x)]/dx = e^u(x)*du(x)/dx
     * 
     * with u(x) = -c*Rij^2   and c=(alpha_i*alpha_j)/(alpha_i+alpha_j)
     * 
     * u(x) = -c*(xi^2-2xi*xj+xj^2)
     * 
     * du(x)/dx=-c*(-2xi+2xj)
     * 
     * dVij = -Vij*(alpha_i*alpha_j)/(alpha_i+alpha_j) * (2xj - 2xi)
     * 
     * 
     * 
	 */
	
	public double getFGValueShape(double[] grad) {
		
		ShapeVolume molGauss = shapeAlign.getMolGauss();

		ShapeVolume refMolGauss = shapeAlign.getRefMolGauss();

		for(int i=0;i<grad.length;i++) {
			grad[i] = 0;
		}

	    double totalOverlap = 0.0;
		for(AtomicGaussian refAt:refMolGauss.getAtomicGaussians()){
			double xi = refAt.getCenter().x;
			double yi = refAt.getCenter().y;
			double zi = refAt.getCenter().z;
			for(AtomicGaussian fitAt:molGauss.getAtomicGaussians()){
				int a = fitAt.getAtomId();
				double atomOverlap = 0.0;
				double xj = fitConf.getX(a);
				double yj = fitConf.getY(a);
				double zj = fitConf.getZ(a);
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
				grad[3*a] += (xj-xi)*gradientPrefactor;
				grad[3*a+1] += (yj-yi)*gradientPrefactor;
				grad[3*a+2] += (zj-zi)*gradientPrefactor;
				}

		
		}
		if(refMolGauss instanceof MolecularVolume) {
			for(VolumeGaussian refVG:((MolecularVolume)refMolGauss).getVolumeGaussians()){
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
				grad[3*a] += (xj-xi)*gradientPrefactor;
				grad[3*a+1] += (yj-yi)*gradientPrefactor;
				grad[3*a+2] += (zj-zi)*gradientPrefactor;
				}

		
		}
		}

		return totalOverlap; 
	
	}
	
	public double getFGValuePP() {

		ShapeVolume molGauss = shapeAlign.getMolGauss();

		ShapeVolume refMolGauss = shapeAlign.getRefMolGauss();

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
				double atomOverlap = 0.0;
				Coordinates fitCenterModCoord = fitPP.getCenter();
				double xj = fitCenterModCoord.x;
				double yj = fitCenterModCoord.y;
				double zj = fitCenterModCoord.z;
				double dx = xi-xj;
				double dy = yi-yj;
				double dz = zi-zj;
				double Rij2 = dx*dx + dy*dy + dz*dz;
				double alphaSum = refPP.getWidth() + fitPP.getWidth();
				if(Rij2<Gaussian3D.DIST_CUTOFF) {
					atomOverlap = refPP.getWeight()*refPP.getHeight()*fitPP.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refPP.getWidth() * fitPP.getWidth()* Rij2)/alphaSum) *
							QuickMathCalculator.getInstance().getPrefactor(refPP.getAtomicNo(),fitPP.getAtomicNo());
					if (atomOverlap>0.0) {
						double sim = refPP.getInteractionSimilarity(fitPP);
						atomOverlap *= sim;
						totalOverlap += atomOverlap;
					}

				}

				}
		
		}

		return totalOverlap; 
	
	}
	
	public double getFGValueShapeSelf(double[] grad, ShapeVolume molGauss,boolean rigid) {
		
		for(int i=0;i<grad.length;i++) {
			grad[i] = 0;
		}

		/**
		 * derivative of ShapeOverlap with respect to Cartesian coordinates
		 */ 
	    double totalOverlap = 0.0;
	    for(AtomicGaussian refAt:molGauss.getAtomicGaussians()){
			for(AtomicGaussian  fitAt:molGauss.getAtomicGaussians()){
				totalOverlap+=getGradientContributionSelf(refAt,fitAt,grad,rigid);
			}
			if(molGauss instanceof MolecularVolume) {
			for(VolumeGaussian fitAt:((MolecularVolume)molGauss).getVolumeGaussians()){
				totalOverlap+=getGradientContributionSelf(refAt,fitAt,grad,rigid);
			}
			}
	    }
		if(molGauss instanceof MolecularVolume) {
	    for(VolumeGaussian refAt: ((MolecularVolume)molGauss).getVolumeGaussians()){
	    	for(VolumeGaussian fitAt : ((MolecularVolume)molGauss).getVolumeGaussians()){
	    		totalOverlap+=getGradientContributionSelf(refAt,fitAt,grad,rigid);
	    	}
	    }
		}

		return totalOverlap; 

		
	
	}
	
	public double getGradientContributionSelf(Gaussian3D refAt, Gaussian3D fitAt, double[] grad,boolean rigid) {
		double xi,yi,zi,xj,yj,zj;
		int b = fitAt.getAtomId();
		double atomOverlap = 0.0;
		int a = refAt.getAtomId();
		if(rigid) { // rigid is always a self overlap of the ref molecule
			xi = refAt.getCenter().x;
			yi = refAt.getCenter().y;
			zi = refAt.getCenter().z;
		}
		else  {
			xi = fitConf.getX(a);
			yi = fitConf.getY(a);
			zi = fitConf.getZ(a);
		}

		if(rigid) {
			xj = fitAt.getCenter().x;
			yj = fitAt.getCenter().y;
			zj = fitAt.getCenter().z;
		}
		else {
			xj = fitConf.getX(b);
			yj = fitConf.getY(b);
			zj = fitConf.getZ(b);
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
	

	
	
	public double getFGValueSelfPP(ShapeVolume molVol,boolean rigid) {
		double xi,yi,zi,xj,yj,zj;

		/**
		 * derivative of ShapeOverlap with respect to Cartesian coordinates
		 */ 
	    double totalOverlap = 0.0;
	    for(PPGaussian refPP:molVol.getPPGaussians()){
			if(rigid) {
				xi = refPP.getCenter().x;
				yi = refPP.getCenter().y;
				zi = refPP.getCenter().z;
			}
			else {
				xi = refPP.getCenter().x;
				yi = refPP.getCenter().y;
				zi = refPP.getCenter().z;
			}
			for(PPGaussian fitPP:molVol.getPPGaussians()){
				double atomOverlap = 0.0;

				if(rigid) {
					xj = fitPP.getCenter().x;
					yj = fitPP.getCenter().y;
					zj = fitPP.getCenter().z;
				}
				else {
					xj = fitPP.getCenter().x;
					yj = fitPP.getCenter().y;
					zj = fitPP.getCenter().z;
				}
				double dx = xi-xj;
				double dy = yi-yj;
				double dz = zi-zj;
				double Rij2 = dx*dx + dy*dy + dz*dz;
				double alphaSum = fitPP.getWidth() + fitPP.getWidth();
				if(Rij2<Gaussian3D.DIST_CUTOFF) {
					atomOverlap = refPP.getWeight()*refPP.getHeight()*fitPP.getHeight()*QuickMathCalculator.getInstance().quickExp(-( refPP.getWidth() * fitPP.getWidth()* Rij2)/alphaSum) *
							QuickMathCalculator.getInstance().getPrefactor(refPP.getAtomicNo(),fitPP.getAtomicNo());
					
					if (atomOverlap>0.0) {
						double sim = refPP.getInteractionSimilarity(fitPP);
						atomOverlap *= sim;
						totalOverlap += atomOverlap;
					}
				}
				}
		}
		return totalOverlap; 
	}
	
	public double[] getCartState(){
		double[] cartState = new double[3*fitConf.getMolecule().getAllAtoms()];
		for(int a=0;a<fitConf.getMolecule().getAllAtoms();a++) {
			cartState[3*a] = fitConf.getCoordinates(a).x;
			cartState[3*a+1] = fitConf.getCoordinates(a).y;
			cartState[3*a+2] = fitConf.getCoordinates(a).z;
			
		}
		return cartState;
	}
	
	public Conformer getFitConf() {
		return fitConf;
	}
	
	
	@Override
	public EvaluableFlexibleOverlap clone() {
		return new EvaluableFlexibleOverlap(this);
	}
		
}
