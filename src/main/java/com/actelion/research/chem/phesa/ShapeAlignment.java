package com.actelion.research.chem.phesa;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.interactionstatistics.InteractionSimilarityTable;
import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.SingularValueDecomposition;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.DoubleStream;
import java.util.ArrayList;

/** 
 * @version: 1.0, February 2018
 * Author: J. Wahl
 * this class provides functionalities to calculate the overlap between two molecules
 * 
*/


public class ShapeAlignment {

	private MolecularVolume refMolGauss;
	private MolecularVolume molGauss;

	
	
	
	public ShapeAlignment(StereoMolecule refMol, StereoMolecule mol) {
		this.refMolGauss= new MolecularVolume(refMol);
		this.molGauss = new MolecularVolume(mol);
	}
	
	public ShapeAlignment(MolecularVolume refMolGauss, MolecularVolume molGauss) {
		this.refMolGauss= refMolGauss;
		this.molGauss = molGauss;
	}
	


	public MolecularVolume getRefMolGauss() {
		return refMolGauss;
	}

	public MolecularVolume getMolGauss() {
		return molGauss;
	}

	/**
	 * Move COM of the molecular volume to the origin of the lab-frame and orient molecules so that their principal moments
	 * of inertia coincide with the 3 axis of the coordinate system
	 * @param mol
	 * @param molVol
	 */
	public static Matrix preProcess(StereoMolecule mol, MolecularVolume molVol) {
		Coordinates COM = molVol.getCOM();
		int nrOfAtoms = mol.getAllAtoms();	

		
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords1 = mol.getCoordinates(i);
			coords1.sub(COM);
		}

		for (AtomicGaussian ag : molVol.getAtomicGaussians()){
			ag.getCenter().sub(COM);  //translate atomicGaussians. Moves center of mass to the origin.
		}

		
		for (PPGaussian pg : molVol.getPPGaussians()){
			pg.getCenter().sub(COM);  //translate atomicGaussians. Moves center of mass to the origin.
		}
		return initialOrientation(mol,molVol);
	}


	/**
	 * 
	 * Calculate singular value decomposition of the covariance matrix of volume-weighted atomic coordinates by
	 * applying a singular value decomposition, we obtain the principal moments of inertia. The molecule is aligned
	 * so that the axis of the lab-frames coincide with the principal moments of intertia.
	 * @return
	 */
	
	public static Matrix initialOrientation(StereoMolecule mol,MolecularVolume molGauss) {
		
		Matrix massMatrix = new Matrix(3,3); 
		double volume = 0.0;
		for (AtomicGaussian ag : molGauss.getAtomicGaussians()){
			volume += ag.getVolume();
			double value = ag.getVolume()*ag.getCenter().x*ag.getCenter().x;
			massMatrix.addToElement(0,0,value);
			value = ag.getVolume()*ag.getCenter().x*ag.getCenter().y;
			massMatrix.addToElement(0,1,value);
			value = ag.getVolume()*ag.getCenter().x*ag.getCenter().z;
			massMatrix.addToElement(0,2,value);
			value = ag.getVolume()*ag.getCenter().y*ag.getCenter().y;
			massMatrix.addToElement(1,1,value);
			value = ag.getVolume()*ag.getCenter().y*ag.getCenter().z;
			massMatrix.addToElement(1,2,value);
			value = ag.getVolume()*ag.getCenter().z*ag.getCenter().z;
			massMatrix.addToElement(2,2,value);	
		}
		massMatrix.set(0,0,massMatrix.get(0,0)/volume);
		massMatrix.set(0,1,massMatrix.get(0,1)/volume);
		massMatrix.set(0,2,massMatrix.get(0,2)/volume);
		massMatrix.set(1,1,massMatrix.get(1,1)/volume);
		massMatrix.set(1,2,massMatrix.get(1,2)/volume);
		massMatrix.set(2,2,massMatrix.get(2,2)/volume);
		massMatrix.set(1,0,massMatrix.get(0,1));
		massMatrix.set(2,0,massMatrix.get(0,2));
		massMatrix.set(2,1,massMatrix.get(1,2));
		SingularValueDecomposition svd = new SingularValueDecomposition(massMatrix.getArray(),null,null);
		Matrix u= new Matrix(svd.getU());
		double det = u.det();
		if(det<0) {
			u.set(0,2,-u.get(0, 2));
			u.set(1,2,-u.get(1, 2));
			u.set(2,2,-u.get(2, 2));
		}
		
		for (AtomicGaussian ag : molGauss.getAtomicGaussians()){
			ag.getCenter().rotate(u.getArray());
		}
		
		for (PPGaussian pg : molGauss.getPPGaussians()){
			pg.getCenter().rotate(u.getArray());

		}

		
		int nrOfAtoms = mol.getAllAtoms();	
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords1 = mol.getCoordinates(i);
			coords1.rotate(u.getArray());
		}
		
		return u;
		
	}
	
	
	/**
	 * .
	 * generate initial orientations of the molecule: 
     * mode1: 4 orientations: initial orientation and 180 degree rotation about each axis
	 * mode2: mode1 and 90 degree rotations about each axis
	 * a transformation vector consists of 7 elements: the first 4 elements form a Quaternion and describe the rotation
	 * the last three elements are the translation vector
	 * @param mode
	 * @return
	 */
	public static double[][] initialTransform(int mode) {
		
		double cos45 = Math.cos(Math.toRadians(45.0));
		double cos90 = Math.cos(Math.toRadians(90.0));
		double cos135 = Math.cos(Math.toRadians(135));
		double sin45 = Math.sin(Math.toRadians(45.0));
		double sin90 = Math.sin(Math.toRadians(90.0));
		double sin135 = Math.sin(Math.toRadians(135));
	
		switch(mode){
		case 1:
			double[][] transforms1 = {{1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,1.0,0.0,0.0,0.0,0.0},
					{0.0,0.0,0.0,1.0,0.0,0.0,0.0}};
			return transforms1;
		case 2:
			double[][] transforms2 = {{1.0,0.0,0.0,0.0,0.0,0.0,0.0},{cos45,sin45,0.0,0.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0,0.0,0.0,0.0},{cos135,sin135,0.0,0.0,0.0,0.0,0.0},
					{cos45,0.0,sin45,0.0,0.0,0.0,0.0},{0.0,0.0,1.0,0.0,0.0,0.0,0.0},{cos135,0.0,sin135,0.0,0.0,0.0,0.0},
					{cos45,0.0,0.0,sin45,0.0,0.0,0.0},{0.0,0.0,0.0,1.0,0.0,0.0,0.0},{cos135,0.0,0.0,sin135,0.0,0.0,0.0}};
			return transforms2;
		
			
		default:
		
			double [][] transform = {{1.0,0.0,0.0,0.0,0.0,0.0,0.0}};
			return transform;
		}
	
			
	}
		
	/**
	 * calculate the Overlap of the two molecular volumes as a function a transform vector that is applied to the query molecule
	 * overlap Volume of two molecular Volumes  formulated as a summed overlap of atomic Gaussians
	 * taken from Grant, Gallardo, Pickup, Journal of Computational Chemistry, 17, 1653-1666, 1996
	 * @param transform
	 * @return
	 */
	

	public double getTotalAtomOverlap(double[] transform){
		Quaternion quat = new Quaternion(transform[0],transform[1],transform[2],transform[3]);
		double Vtot = 0.0;
		double[][] rotMatrix = quat.getRotMatrix().getArray();
		ArrayList<AtomicGaussian> atomicGaussians = molGauss.getAtomicGaussians();
		Coordinates[] fitCenterModCoords = new Coordinates[atomicGaussians.size()];
		double normFactor = 1/(transform[0]*transform[0]+transform[1]*transform[1]+transform[2]*transform[2]+transform[3]*transform[3]);
		for(int k=0;k<atomicGaussians.size();k++) {
    			fitCenterModCoords[k] = atomicGaussians.get(k).getRotatedCenter(rotMatrix, normFactor, new double[] {transform[4], transform[5], transform[6]}); //we operate on the transformed coordinates of the molecule to be fitted


		}

		for(AtomicGaussian refAt:refMolGauss.getAtomicGaussians()){
			int index = 0;
			for(AtomicGaussian fitAt:molGauss.getAtomicGaussians()){
				Vtot += refAt.getQuickVolumeOverlap(fitAt, fitCenterModCoords[index],10.0);
				index+=1;	
			}
		}
		
		return Vtot;

	}
	
	public double getTotalAtomOverlap(){
		double Vtot = 0.0;
		for(AtomicGaussian refAt:refMolGauss.getAtomicGaussians()){
			for(AtomicGaussian fitAt:molGauss.getAtomicGaussians()) {
					Vtot+=refAt.getQuickVolumeOverlap(fitAt,10.0);
			}
						
			}				
		
		return Vtot;

	}
	
	

	
	
	public double getTotalPPOverlap(double[] transform){
		Quaternion quat = new Quaternion(transform[0],transform[1],transform[2],transform[3]);
		double Vtot = 0.0;
		double[][] rotMatrix = quat.getRotMatrix().getArray();
		ArrayList<PPGaussian> ppGaussians = molGauss.getPPGaussians();
		Coordinates[] fitCenterModCoords = new Coordinates[ppGaussians.size()];
		Coordinates[] fitDirectionalityMod = new Coordinates[ppGaussians.size()];
		double normFactor = 1/(transform[0]*transform[0]+transform[1]*transform[1]+transform[2]*transform[2]+transform[3]*transform[3]);
	    for(int k=0;k<ppGaussians.size();k++) {
	    	fitCenterModCoords[k]=  ppGaussians.get(k).getRotatedCenter(rotMatrix, normFactor, new double[] {transform[4],transform[5],transform[6]});    //we operate on the transformed coordinates of the molecule to be fitted
	    	fitDirectionalityMod[k] = ppGaussians.get(k).getRotatedDirectionality(rotMatrix, normFactor);
	    }

		for(PPGaussian refPP:refMolGauss.getPPGaussians()){
			int index = 0;
			for(PPGaussian fitPP:molGauss.getPPGaussians()){
				Vtot+=refPP.getSimilarity(fitPP, fitDirectionalityMod[index])* refPP.getVolumeOverlap(fitPP, fitCenterModCoords[index]);
				index+=1;
			
		}
			
		}
		return Vtot;

	}



	/**
	 * calculate the self-overlap of the base molecule
	 * @return
	 */
	
	
	public double getSelfAtomOverlap(MolecularVolume molGauss){
			double Vtot = 0.0;
			for(AtomicGaussian at:molGauss.getAtomicGaussians()){
				for(AtomicGaussian at2:molGauss.getAtomicGaussians()){
					Vtot += at.getVolumeOverlap(at2);
	
				}
			}
	
			return Vtot;
		}




	public double getSelfPPOverlap(MolecularVolume molGauss){
		double Vtot = 0.0;
		for(PPGaussian pp:molGauss.getPPGaussians()){
			for(PPGaussian pp2:molGauss.getPPGaussians()){
				Vtot+=pp.getSimilarity(pp2)* pp.getVolumeOverlap(pp2);
			}
		}
		return Vtot;
	}
	

		
	
	public double getSelfAtomOverlapRef(){
		
		return getSelfAtomOverlap(refMolGauss);
	}
	
	public double getSelfAtomOverlapFit(){
		
		return getSelfAtomOverlap(molGauss);
	}
	
	public double getSelfPPOverlapRef(){
		
		return getSelfPPOverlap(refMolGauss);
	}
	
	public double getSelfPPOverlapFit(){
		
		return getSelfPPOverlap(molGauss);
	}
		

	
	
	public static void rotateMol(StereoMolecule mol,Quaternion rotor, double[] transl) {

		double normFactor = 1/rotor.normSquared();

		
		int nrOfAtoms = mol.getAllAtoms();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords = mol.getCoordinates(i);
			double[][] m = rotor.getRotMatrix().getArray();
			double x0 = coords.x;
			double y0 = coords.y;
			double z0 = coords.z;
			coords.x = x0*m[0][0]+y0*m[0][1]+z0*m[0][2];
			coords.y = x0*m[1][0]+y0*m[1][1]+z0*m[1][2];
			coords.z = x0*m[2][0]+y0*m[2][1]+z0*m[2][2];
			
			coords.scale(normFactor);
			coords.add(transl[0],transl[1],transl[2]);

		}
		
	}
	
	
	public double[] findAlignment(double[][] transforms) {
		double [] alignment = {1.0,0.0,0.0,0.0,0.0,0.0,0.0};
		double Oaa = getSelfAtomOverlapRef();
		double Obb = getSelfAtomOverlapFit();
		double ppOaa = getSelfPPOverlapRef();
		double ppObb = getSelfPPOverlapFit();
		EvaluableOverlap eval = new EvaluableOverlap(this, new double[7]);
		ShapeOptimizerLBFGS opt = new ShapeOptimizerLBFGS(200,0.001);
		double maxSimilarity = 0.0;
		double ppScaling = 1.0;
		for(double [] transform:transforms) { //iterate over all initial alignments (necessary since optimizer just finds next local minimum, so we need different initial guesses
			eval.setState(transform);
			double[] bestTransform = opt.optimize(eval);
			double atomOverlap = 0.0;
			double ppOverlap = 0.0;
			float similarity = 0.0f;
			atomOverlap = getTotalAtomOverlap(bestTransform);
			ppOverlap = getTotalPPOverlap(bestTransform);
			float atomSimilarity = (float)(atomOverlap/(Oaa+Obb-atomOverlap));
			float ppSimilarity = 0.0f;
			if(getRefMolGauss().getPPGaussians().size()==0 && getMolGauss().getPPGaussians().size()==0 )
				ppSimilarity = 1.0f;
			else ppSimilarity=(float)(ppOverlap/(ppOaa+ppObb-ppOverlap));
			similarity = (1.0f/(1+(float)ppScaling))* (atomSimilarity + (float)ppScaling*ppSimilarity) ;
			if (similarity>maxSimilarity) {
				maxSimilarity = similarity;
				alignment = bestTransform;
			}
		}
			
			return DoubleStream.concat(Arrays.stream(new double[] {maxSimilarity}), Arrays.stream(alignment)).toArray();
		}
		
		
	
	
	public static void rotateMol(StereoMolecule mol, double[] transform) {
		Quaternion rotor = new Quaternion(transform[0],transform[1], transform [2], transform[3]);
		double[] translate =  {transform[4], transform[5], transform[6]};
		ShapeAlignment.rotateMol(mol, rotor, translate);

	}
	
	public double getShapeSimilarityWithoutOptimization() {
		double Oaa = getSelfAtomOverlapRef();
		double Obb = getSelfAtomOverlapFit();
		double Oab = getTotalAtomOverlap();
		return (Oab)/(Oaa+Obb-Oab);
	}
	
}

	
	
	

	

	
	
	
	


