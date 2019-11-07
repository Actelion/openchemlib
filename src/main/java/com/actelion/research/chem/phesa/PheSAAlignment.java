package com.actelion.research.chem.phesa;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.phesa.pharmacophore.PPGaussian;
import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.SingularValueDecomposition;
import java.util.Arrays;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.ArrayList;

/** 
 * @version: 1.0, February 2018
 * Author: J. Wahl
 * this class provides functionalities to calculate the overlap between two molecules
 * 
*/


public class PheSAAlignment {

	private MolecularVolume refMolGauss;
	private MolecularVolume molGauss;
	public enum axis {X,Y,Z};


	
	
	
	public PheSAAlignment(StereoMolecule refMol, StereoMolecule mol) {
		this.refMolGauss = new MolecularVolume(refMol);
		this.molGauss = new MolecularVolume(mol);
	}
	
	public PheSAAlignment(MolecularVolume refMolGauss, MolecularVolume molGauss) {
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
	public static Matrix preProcess(Conformer conf, MolecularVolume molVol) {
		Coordinates COM = molVol.getCOM();
		int nrOfAtoms = conf.getSize();

		
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords1 = conf.getCoordinates(i);
			coords1.sub(COM);
		}

		molVol.translateToCOM(COM);
		

		return createCanonicalOrientation(conf,molVol);
	}
	
	
	
	


	
	public static Matrix createCanonicalOrientation(Conformer conf,MolecularVolume molGauss) {
		Matrix m = PheSAAlignment.getCovarianceMatrix(molGauss);
		SingularValueDecomposition svd = new SingularValueDecomposition(m.getArray(),null,null);
		Matrix u = new Matrix(svd.getU());
		double det = u.det();
		if(det<0) {
			u.set(0,1,-u.get(0, 1));
			u.set(1,1,-u.get(1, 1));
			u.set(2,1,-u.get(2, 1));
		}
		rotateMol(conf,u);
		molGauss.update(conf);
		Matrix rotMat = u;
		
		if(!isCanonicalOrientation(molGauss)) {
			rotateMolAroundAxis180(conf,axis.X);
			molGauss.update(conf);
			if(isCanonicalOrientation(molGauss)) {
				u.set(0,1,-u.get(0, 1));
				u.set(1,1,-u.get(1, 1));
				u.set(2,1,-u.get(2, 1));
				u.set(0,2,-u.get(0, 2));
				u.set(1,2,-u.get(1, 2));
				u.set(2,2,-u.get(2, 2));
				rotMat = u;
			}
			else {
				rotateMolAroundAxis180(conf,axis.X); // rotate back
				molGauss.update(conf);
				rotateMolAroundAxis180(conf,axis.Y);
				molGauss.update(conf);
				if(isCanonicalOrientation(molGauss)) {
					u.set(0,0,-u.get(0, 0));
					u.set(1,0,-u.get(1, 0));
					u.set(2,0,-u.get(2, 0));
					u.set(0,2,-u.get(0, 2));
					u.set(1,2,-u.get(1, 2));
					u.set(2,2,-u.get(2, 2));
					rotMat = u;
				}
				else {
					rotateMolAroundAxis180(conf,axis.Y);
					molGauss.update(conf);
					rotateMolAroundAxis180(conf,axis.Z);
					molGauss.update(conf);
					if(isCanonicalOrientation(molGauss)) {
						u.set(0,0,-u.get(0, 0));
						u.set(1,0,-u.get(1, 0));
						u.set(2,0,-u.get(2, 0));
						u.set(0,1,-u.get(0, 1));
						u.set(1,1,-u.get(1, 1));
						u.set(2,1,-u.get(2, 1));
						rotMat = u;
					}
				}
			}
		}
		
		return rotMat;
	}
	
	private static void rotateMolAroundAxis180(Conformer conf,axis a) {
		if (a == axis.X) {
			IntStream.range(0,conf.getSize()).forEach(i -> {
				Coordinates coords = conf.getCoordinates(i);
				coords.y = -coords.y;
				coords.z = -coords.z;
			});
		}
		else if (a == axis.Y) {
			IntStream.range(0,conf.getSize()).forEach(i -> {
				Coordinates coords = conf.getCoordinates(i);
				coords.x = -coords.x;
				coords.z = -coords.z;
			});
		}
		
		else  {
			IntStream.range(0,conf.getSize()).forEach(i -> {
				Coordinates coords = conf.getCoordinates(i);
				coords.x = -coords.x;
				coords.y = -coords.y;
			});
		}

	}
	
	private static Matrix getCovarianceMatrix(MolecularVolume molGauss) {
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
		
		return massMatrix;
	}
	
	private static boolean isCanonicalOrientation(MolecularVolume molGauss) {
		double xxPos = 0;
		double xxNeg = 0;
		double yyPos = 0;
		double yyNeg = 0;
		int nXPos = 0;
		int nXNeg = 0;
		int nYPos = 0;
		int nYNeg = 0;
		
		for (AtomicGaussian ag : molGauss.getAtomicGaussians()){
			double x = ag.center.x;
			double y = ag.center.y;
						
			if(x>0) {
				xxPos += x*x;
				nXPos++;
			}
			else { 
				xxNeg += x*x;
				nXNeg++;
			}
			
			if(y>0) {
				yyPos += y*y;
				nYPos++;
			}
			else { 
				yyNeg += y*y;
				nYNeg++;
			}

		}
		
		xxPos/=nXPos;
		yyPos/=nYPos;	
		xxNeg/=nXNeg;
		yyNeg/=nYNeg;	
		
		if(xxPos>xxNeg && yyPos>yyNeg)
			return true;
		else
			return false;
		
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
		
		double c = 0.707106781;
	
		switch(mode){
		case 1:
			double[][] transforms1 = {{1.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,1.0,0.0,0.0,0.0,0.0},
					{0.0,0.0,0.0,1.0,0.0,0.0,0.0}};
			return transforms1;
		case 2:
			double[][] transforms2 = {{1,0,0,0,0,0,0},{0,1,0,0,0,0,0},{0,0,1,0,0,0,0},
					{0,0,0,1,0,0,0},
					{c,c,0,0,0,0,0},
					{c,0,c,0,0,0,0},
					{c,0,0,c,0,0,0},
					{-0.5,0.5,0.5,-0.5,0,0,0},
					{0.5,-0.5,0.5,-0.5,0,0,0},
					{0.5,0.5,0.5,-0.5,0,0,0},
					{0.5,-0.5,-0.5,-0.5,0,0,0},
					{0.5,0.5,-0.5,-0.5,0,0,0}
					};
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
				Vtot += refAt.getVolumeOverlap(fitAt, fitCenterModCoords[index],Gaussian3D.DIST_CUTOFF);
				index+=1;	
			}
		}
		
		for(ExclusionGaussian refEx:refMolGauss.getExclusionGaussians()){
			int index = 0;
			for(AtomicGaussian fitAt:molGauss.getAtomicGaussians()){
				Vtot -= refEx.getVolumeOverlap(fitAt, fitCenterModCoords[index],Gaussian3D.DIST_CUTOFF);
				index+=1;	
			}
		}
		if(Vtot<0)
			Vtot = 0.0;
		
		return Vtot;

	}
	
	public double getTotalAtomOverlap(){
		double Vtot = 0.0;
		for(AtomicGaussian refAt:refMolGauss.getAtomicGaussians()){
			for(AtomicGaussian fitAt:molGauss.getAtomicGaussians()) {
					Vtot+=refAt.getVolumeOverlap(fitAt,Gaussian3D.DIST_CUTOFF);
			}
						
			}				
		
		return Vtot;

	}
	
	
	
	
	public double getTotalPPOverlap(double[] transform){
		Quaternion quat = new Quaternion(transform[0],transform[1],transform[2],transform[3]);
		double Vtot = 0.0;
		double correctionFactor = refMolGauss.getPPGaussians().size()/refMolGauss.getPPGaussians().stream().mapToDouble(g -> g.getWeight()).sum();
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
				Vtot+=refPP.getSimilarity(fitPP, fitDirectionalityMod[index])* refPP.getVolumeOverlap(fitPP, fitCenterModCoords[index],10.0)*correctionFactor;
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
		double correctionFactor = molGauss.getPPGaussians().size()/molGauss.getPPGaussians().stream().mapToDouble(g -> g.getWeight()).sum();
		double Vtot = 0.0;
		for(PPGaussian pp:molGauss.getPPGaussians()){
			for(PPGaussian pp2:molGauss.getPPGaussians()){
				Vtot+=pp.getSimilarity(pp2)* pp.getVolumeOverlap(pp2)*correctionFactor;
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
		

	
	
	public static void rotateMol(Conformer conf,Quaternion rotor, double[] transl) {

		double normFactor = 1/rotor.normSquared();

		
		int nrOfAtoms = conf.getSize();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords = conf.getCoordinates(i);
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
	
	public static void rotateMol(StereoMolecule mol,double[][] m) {
		int nrOfAtoms = mol.getAllAtoms();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords = mol.getCoordinates(i);
			double x0 = coords.x;
			double y0 = coords.y;
			double z0 = coords.z;
			coords.x = x0*m[0][0]+y0*m[0][1]+z0*m[0][2];
			coords.y = x0*m[1][0]+y0*m[1][1]+z0*m[1][2];
			coords.z = x0*m[2][0]+y0*m[2][1]+z0*m[2][2];
		}
		
	}
	
	public static void rotateCoords(Coordinates coords,double[][] m) {
		double x0 = coords.x;
		double y0 = coords.y;
		double z0 = coords.z;
		coords.x = x0*m[0][0]+y0*m[0][1]+z0*m[0][2];
		coords.y = x0*m[1][0]+y0*m[1][1]+z0*m[1][2];
		coords.z = x0*m[2][0]+y0*m[2][1]+z0*m[2][2];
		
	}
	
	public static void translateMol(StereoMolecule mol,double[] translate) {
		int nrOfAtoms = mol.getAllAtoms();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords = mol.getCoordinates(i);
			coords.x += translate[0];
			coords.y += translate[1];
			coords.z += translate[2];
		}
		
	}
	
	public static void multiplyMatrix(double[][] r, double[][] s, double[][] rs) {
		rs[0][0] = r[0][0]*s[0][0] + r[0][1]*s[1][0] + r[0][2]*s[2][0];
		rs[0][1] = r[0][0]*s[0][1] + r[0][1]*s[1][1] + r[0][2]*s[2][1];
		rs[0][2] = r[0][0]*s[0][2] + r[0][1]*s[1][2] + r[0][2]*s[2][2];
		
		rs[1][0] = r[1][0]*s[0][0] + r[1][1]*s[1][0] + r[1][2]*s[2][0];
		rs[1][1] = r[1][0]*s[0][1] + r[1][1]*s[1][1] + r[1][2]*s[2][1];
		rs[1][2] = r[1][0]*s[0][2] + r[1][1]*s[1][2] + r[1][2]*s[2][2];
		
		rs[2][0] = r[2][0]*s[0][0] + r[2][1]*s[1][0] + r[2][2]*s[2][0];
		rs[2][1] = r[2][0]*s[0][1] + r[2][1]*s[1][1] + r[2][2]*s[2][1];
		rs[2][2] = r[2][0]*s[0][2] + r[2][1]*s[1][2] + r[2][2]*s[2][2];
	}
	
	public static void multiplyInverseMatrix(double[][] r, double[][] s, double[][] rs) {
		rs[0][0] = r[0][0]*s[0][0] + r[0][1]*s[0][1] + r[0][2]*s[0][2];
		rs[0][1] = r[0][0]*s[1][0] + r[0][1]*s[1][1] + r[0][2]*s[1][2];
		rs[0][2] = r[0][0]*s[2][0] + r[0][1]*s[2][1] + r[0][2]*s[2][2];
		
		rs[1][0] = r[1][0]*s[0][0] + r[1][1]*s[0][1] + r[1][2]*s[0][2];
		rs[1][1] = r[1][0]*s[1][0] + r[1][1]*s[1][1] + r[1][2]*s[1][2];
		rs[1][2] = r[1][0]*s[2][0] + r[1][1]*s[2][1] + r[1][2]*s[2][2];
		
		rs[2][0] = r[2][0]*s[0][0] + r[2][1]*s[0][1] + r[2][2]*s[0][2];
		rs[2][1] = r[2][0]*s[1][0] + r[2][1]*s[1][1] + r[2][2]*s[1][2];
		rs[2][2] = r[2][0]*s[2][0] + r[2][1]*s[2][1] + r[2][2]*s[2][2];
	}
	
	public static void getRotationMatrix(double theta, Coordinates axis, double[][] r) {
		double x = axis.x;
		double y = axis.y;
		double z = axis.z;
		double c = Math.cos(theta);
		double s = Math.sin(theta);
		double t = 1-c;
		r[0][0] = c+x*x*t;
		r[0][1] = x*y*t-z*s;
		r[0][2] = x*z*t+y*s;
		r[1][0] = x*y*t+z*s;
		r[1][1] = c+y*y*t;
		r[1][2] = y*z*t-x*s;
		r[2][0] = z*x*t-y*s;
		r[2][1] = z*y*t+x*s;
		r[2][2] = c+z*z*t;

	}
	
	
	public double[] findAlignment(double[][] transforms) {
		return findAlignment(transforms,true);
	}
	
	public double[] findAlignment(double[][] transforms, boolean optimize) {
		double [] alignment = {1.0,0.0,0.0,0.0,0.0,0.0,0.0};
		double Oaa = getSelfAtomOverlapRef();
		double Obb = getSelfAtomOverlapFit();
		double ppOaa = getSelfPPOverlapRef();
		double ppObb = getSelfPPOverlapFit();
		EvaluableOverlap eval = new EvaluableOverlap(this, new double[7]);
		OptimizerLBFGS opt = new OptimizerLBFGS(200,0.001);
		double maxSimilarity = 0.0;
		double ppScaling = 1.0;
		for(double [] transform:transforms) { //iterate over all initial alignments (necessary since optimizer just finds next local minimum, so we need different initial guesses
			eval.setState(transform);
			double[] bestTransform;
			if(optimize)
				bestTransform = opt.optimize(eval);
			else
				bestTransform = transform;
			double atomOverlap = 0.0;
			double ppOverlap = 0.0;
			float similarity = 0.0f;
			ppOverlap = getTotalPPOverlap(bestTransform);
			float ppSimilarity = 0.0f;
			if(getRefMolGauss().getPPGaussians().size()==0 && getMolGauss().getPPGaussians().size()==0 )
				ppSimilarity = 1.0f;
			else ppSimilarity=(float)(ppOverlap/(ppOaa+ppObb-ppOverlap));
			if(ppSimilarity>1.0) //can happen because of weights
				ppSimilarity = 1.0f;
			atomOverlap = getTotalAtomOverlap(bestTransform);
			float atomSimilarity = (float)(atomOverlap/(Oaa+Obb-atomOverlap));

			similarity = (1.0f/(1+(float)ppScaling))* (atomSimilarity + (float)ppScaling*ppSimilarity) ;

			if (similarity>maxSimilarity) {
				maxSimilarity = similarity;
				alignment = bestTransform;
			}
		}
			return DoubleStream.concat(Arrays.stream(new double[] {maxSimilarity}), Arrays.stream(alignment)).toArray();
		}
		
		
	
	public static void rotateMol(Conformer conf, Matrix rotMat) {
		int nrOfAtoms = conf.getSize();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords1 = conf.getCoordinates(i);
			coords1.rotate(rotMat.getArray());
		}

	}
	
	
	
	public static void rotateMol(StereoMolecule mol, double[] transform) {
		Quaternion rotor = new Quaternion(transform[0],transform[1], transform [2], transform[3]);
		double[] translate =  {transform[4], transform[5], transform[6]};
		PheSAAlignment.rotateMol(mol, rotor, translate);

	}
	
	public static void rotateMol(StereoMolecule mol, Matrix rotMat) {
		int nrOfAtoms = mol.getAllAtoms();
		for (int i=0;i<nrOfAtoms;i++) {
			Coordinates coords1 = mol.getCoordinates(i);
			coords1.rotate(rotMat.getArray());
		}

	}
	
	
	
	public static void rotateMol(Conformer conf, double[] transform) {
		Quaternion rotor = new Quaternion(transform[0],transform[1], transform [2], transform[3]);
		double[] translate =  {transform[4], transform[5], transform[6]};
		PheSAAlignment.rotateMol(conf, rotor, translate);

	}
	
	public double getShapeSimilarityWithoutOptimization() {
		double Oaa = getSelfAtomOverlapRef();
		double Obb = getSelfAtomOverlapFit();
		double Oab = getTotalAtomOverlap();
		return (Oab)/(Oaa+Obb-Oab);
	}
	
}

	
	
	

	

	
	
	
	


