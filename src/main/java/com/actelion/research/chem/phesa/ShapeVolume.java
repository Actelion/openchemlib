package com.actelion.research.chem.phesa;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.SingularValueDecomposition;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.transformation.ExponentialMap;
import com.actelion.research.chem.alignment3d.transformation.Quaternion;
import com.actelion.research.chem.alignment3d.transformation.Rotation;
import com.actelion.research.chem.alignment3d.transformation.Transformation;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;
import com.actelion.research.chem.alignment3d.transformation.Translation;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.phesa.PheSAAlignment.axis;
import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;
import com.actelion.research.chem.phesa.pharmacophore.pp.ExitVectorPoint;
import com.actelion.research.chem.phesa.pharmacophore.pp.PPGaussian;


public class ShapeVolume {
	
	protected List<PPGaussian> ppGaussians;
	protected List<AtomicGaussian> atomicGaussians;
	protected Coordinates com;
	
	public ShapeVolume() {
		ppGaussians = new ArrayList<>();
		atomicGaussians = new ArrayList<>();
		
	}
	
	public ShapeVolume(ShapeVolume original) {
		this.atomicGaussians = new ArrayList<AtomicGaussian>();
		this.ppGaussians = new ArrayList<PPGaussian>();
		for(AtomicGaussian ag : original.getAtomicGaussians()) {
			this.atomicGaussians.add(new AtomicGaussian(ag));
		}
		for(PPGaussian pg : original.getPPGaussians()) {
			this.ppGaussians.add(new PPGaussian(pg));
		}
		
		this.com = new Coordinates(original.com);
		
	}
	
	public void addPharmacophorePoint(PPGaussian ppGaussian) {
		ppGaussians.add(ppGaussian);
	}
	
	public void addAtomVolume(AtomicGaussian atomGaussian) {
		atomicGaussians.add(atomGaussian);
	}
	
	public void update(StereoMolecule mol) {
		updateCoordinates(getAtomicGaussians(),mol.getAtomCoordinates());
		updateCoordinates(getPPGaussians(),mol.getAtomCoordinates());
	}
	
	public void update(Conformer conf) {
		updateCoordinates(getAtomicGaussians(),conf.getCoordinates());
		updateCoordinates(getPPGaussians(),conf.getCoordinates());
	}
	
	public void transform(Transformation transform) {
		transformGaussians(transform);
	}
	
	/**
	 * Move COM of the molecular volume to the origin of the lab-frame and orient molecules so that their principal moments
	 * of inertia coincide with the 3 axis of the coordinate system
	 * @param mol
	 * @param molVol
	 */
	public Rotation preProcess(Conformer conf) {
		Coordinates COM = getCOM();
		if(conf!=null) {
			int nrOfAtoms = conf.getSize();
	
			
			for (int i=0;i<nrOfAtoms;i++) {
				Coordinates coords1 = conf.getCoordinates(i);
				coords1.sub(COM);
			}
			
		}	translateToCOM(COM);
		
		Matrix rotation =  createCanonicalOrientation(conf);
		Rotation rot = new Rotation(rotation.getArray());
		return rot;
	}
	
	public void removeRings() {
		List<Integer> toRemove = new ArrayList<>();
		int i=0;
		for(PPGaussian ppg : ppGaussians) {
			if (ppg.getPharmacophorePoint().getFunctionalityIndex()==PharmacophoreCalculator.AROM_ID)
				toRemove.add(i);
			i++;
		}
		Collections.reverse(toRemove);
		for(int index : toRemove) {
			ppGaussians.remove(index);
		}
	}
	
	public  Matrix createCanonicalOrientation(Conformer conf) {
		Matrix m = getCovarianceMatrix();
		SingularValueDecomposition svd = new SingularValueDecomposition(m.getArray(),null,null);
		Matrix u = new Matrix(svd.getU());
		double det = u.det();
		if(det<0) {
			u.set(0,1,-u.get(0, 1));
			u.set(1,1,-u.get(1, 1));
			u.set(2,1,-u.get(2, 1));
		}
		Rotation rot = new Rotation(u.getArray());
		if(conf!=null) {
			rot.apply(conf);
			update(conf);

		}
		else {
			transformGaussians(rot);
		}
		
		if(!isCanonicalOrientation()) {
			if(conf!=null) {
				PheSAAlignment.rotateMolAroundAxis180(conf,axis.X);
				update(conf);
			}
			else {
				rotate180DegreeAroundAxis(axis.X);
			}
			if(isCanonicalOrientation()) {
				u.set(0,1,-u.get(0, 1));
				u.set(1,1,-u.get(1, 1));
				u.set(2,1,-u.get(2, 1));
				u.set(0,2,-u.get(0, 2));
				u.set(1,2,-u.get(1, 2));
				u.set(2,2,-u.get(2, 2));

			}
			else {
				if(conf!=null) {
					PheSAAlignment.rotateMolAroundAxis180(conf,axis.X); // rotate back
					update(conf);
					PheSAAlignment.rotateMolAroundAxis180(conf,axis.Y);
					update(conf);
				}
				else {
					rotate180DegreeAroundAxis(axis.X);
					rotate180DegreeAroundAxis(axis.Y);
				}
				if(isCanonicalOrientation()) {
					u.set(0,0,-u.get(0, 0));
					u.set(1,0,-u.get(1, 0));
					u.set(2,0,-u.get(2, 0));
					u.set(0,2,-u.get(0, 2));
					u.set(1,2,-u.get(1, 2));
					u.set(2,2,-u.get(2, 2));
				}
				else {
					if(conf!=null) {
						PheSAAlignment.rotateMolAroundAxis180(conf,axis.Y);
						update(conf);
						PheSAAlignment.rotateMolAroundAxis180(conf,axis.Z);
						update(conf);
					}
					else {
						rotate180DegreeAroundAxis(axis.Y);
						rotate180DegreeAroundAxis(axis.Z);
					}
					
					if(isCanonicalOrientation()) {
						u.set(0,0,-u.get(0, 0));
						u.set(1,0,-u.get(1, 0));
						u.set(2,0,-u.get(2, 0));
						u.set(0,1,-u.get(0, 1));
						u.set(1,1,-u.get(1, 1));
						u.set(2,1,-u.get(2, 1));
					}
					}
				}
		}

		return u;
	}
	
	protected boolean isCanonicalOrientation() {
		double xxPos = 0;
		double xxNeg = 0;
		double yyPos = 0;
		double yyNeg = 0;
		int nXPos = 0;
		int nXNeg = 0;
		int nYPos = 0;
		int nYNeg = 0;
		
		for (AtomicGaussian ag : atomicGaussians){
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
	
	public Coordinates getCOM() {
		this.calcCOM();
		return this.com;
	}
	
	public void calcCOM(){ 
		double volume = 0.0;
		double comX = 0.0;
		double comY = 0.0;
		double comZ = 0.0;
		for(AtomicGaussian atGauss : this.atomicGaussians){
			volume += atGauss.getVolume();
			comX += atGauss.getCenter().x*atGauss.getVolume();
			comY += atGauss.getCenter().y*atGauss.getVolume();
			comZ += atGauss.getCenter().z*atGauss.getVolume();
		}

		comX = comX/volume;
		comY = comY/volume;
		comZ = comZ/volume;
		this.com = new Coordinates(comX,comY,comZ);

	}
	
	
	protected void updateCoordinates(List<? extends Gaussian3D> gaussians, Coordinates[] coords) {
		for(Gaussian3D gaussian : gaussians) {
			gaussian.updateCoordinates(coords);
		}
		
	}
	

	protected void transformGaussians(Transformation transform) {
		transformGaussians(atomicGaussians,transform);
		transformGaussians(ppGaussians,transform);
	}
	
	protected static void transformGaussians(List<? extends Gaussian3D> gaussians, Transformation transform) {
		for (Gaussian3D gaussian : gaussians) {
			gaussian.transform(transform);
		}	
	}
	
	protected static void rotateGaussian(List<? extends Gaussian3D> gaussians, double[][] rot) {
		for (Gaussian3D gaussian : gaussians) {
			Coordinates coords1 = new Coordinates(gaussian.getCenter());
			coords1.rotate(rot);
			gaussian.setCenter(coords1);
		}
		
	}
	
	protected static void rotateGaussians180DegreeAroundAxis(List<? extends Gaussian3D> gaussians ,PheSAAlignment.axis a) {
		for (Gaussian3D gaussian : gaussians) {
			Coordinates coords1 = new Coordinates(gaussian.getCenter());
			PheSAAlignment.rotateCoordsAroundAxis180(coords1, a);
			gaussian.setCenter(coords1);
		}
	}
	
	protected void rotate180DegreeAroundAxis(PheSAAlignment.axis a) {
		rotateGaussians180DegreeAroundAxis(atomicGaussians,a);
		rotateGaussians180DegreeAroundAxis(ppGaussians,a);
	}
	

	
	public String encode() {
		StringBuilder sb = new StringBuilder();
		sb.append(Integer.toString(atomicGaussians.size()));
		sb.append("  ");
		atomicGaussians.forEach(e -> {
			sb.append(e.encode());
			sb.append("  ");
		});
		sb.append(Integer.toString(ppGaussians.size()));
		sb.append("  ");
		ppGaussians.forEach(e -> {
			sb.append(e.encode());
			sb.append("  ");
		});
		return sb.toString();
	}
	
	public static ShapeVolume decode(String s) {
		String[] splitString = s.split("  ");
		int nrOfAtomicGaussians = Integer.decode(splitString[0].trim());
		int firstIndex = 1;
		int lastIndex = 1+nrOfAtomicGaussians;
		List<AtomicGaussian> atomicGaussians = new ArrayList<AtomicGaussian>();
		List<PPGaussian> ppGaussians = new ArrayList<PPGaussian>();
		
		for(int i=firstIndex;i<lastIndex;i++) {
			atomicGaussians.add(AtomicGaussian.fromString(splitString[i].trim()));
		}
		int nrOfPPGaussians = Integer.decode(splitString[lastIndex]);
		firstIndex = lastIndex+1;
		lastIndex = firstIndex + nrOfPPGaussians;
		for(int i=firstIndex;i<lastIndex;i++) {
			ppGaussians.add(PPGaussian.fromString(splitString[i],null));
		}
		ShapeVolume receptorVol = new ShapeVolume();
		atomicGaussians.forEach(e -> receptorVol.addAtomVolume(e));
		ppGaussians.forEach(e -> receptorVol.addPharmacophorePoint(e));
		
		return receptorVol;
		
	}
	
	public List<PPGaussian> getExitVectorGaussians() {
		return ppGaussians.stream().filter(e -> (e.getPharmacophorePoint() instanceof ExitVectorPoint)).collect(Collectors.toList());
	}

	public List<PPGaussian> getPPGaussians() {
		return ppGaussians;
	}

	public void setPPGaussians(List<PPGaussian> ppGaussians) {
		this.ppGaussians = ppGaussians;
	}

	public List<AtomicGaussian> getAtomicGaussians() {
		return atomicGaussians;
	}

	public void setAtomicGaussians(List<AtomicGaussian> atomicGaussians) {
		this.atomicGaussians = atomicGaussians;
	}
	
	public Matrix getCovarianceMatrix() {
		Matrix massMatrix = new Matrix(3,3); 
		double volume = 0.0;
		for (AtomicGaussian ag : atomicGaussians){
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
	
	/**
	 * calculate the self-overlap of the base molecule
	 * @return
	 */
	
	
	public double getSelfAtomOverlap(){
			double Vtot = 0.0;
			for(AtomicGaussian at:atomicGaussians){
				for(AtomicGaussian at2:atomicGaussians){
					Vtot += at.getVolumeOverlap(at2);
				}

			}

			return Vtot;
		}



	public double getSelfPPOverlap(){
		double Vtot = 0.0;
		for(PPGaussian pp:ppGaussians){
			for(PPGaussian pp2:ppGaussians){
				Vtot+=pp.getWeight()*pp.getSimilarity(pp2)* pp.getVolumeOverlap(pp2);
			}
		}
		return Vtot;
	}
	
	/**
	 * calculate the Overlap of the two molecular volumes as a function a transform vector that is applied to the query molecule
	 * overlap Volume of two molecular Volumes  formulated as a summed overlap of atomic Gaussians
	 * taken from Grant, Gallardo, Pickup, Journal of Computational Chemistry, 17, 1653-1666, 1996
	 * returns a double[2]: the first double is the total overlap, whereas the second value is the specific 
	 * contribution of additional volume gaussians (inclusion, exclusion)
	 * @param transform
	 * @return
	 */
	

	public double[] getTotalAtomOverlap(double[] transform, ShapeVolume fitVol){
		double[] result = new double[2];
		ExponentialMap eMap = new ExponentialMap(transform[0],transform[1],transform[2]);
		double Vtot = 0.0;
		double Vvol = 0.0;
		Coordinates com = fitVol.getCOM();
		double[][] rotMatrix = eMap.toQuaternion().getRotMatrix().getArray();
		List<AtomicGaussian> fitGaussians = fitVol.atomicGaussians;
		Coordinates[] fitCenterModCoords = new Coordinates[fitGaussians.size()];
		for(int k=0;k<fitGaussians.size();k++) {
				Coordinates center = new Coordinates(fitGaussians.get(k).getCenter());
				center.sub(com);
			    center.rotate(rotMatrix);
			    center.add(com);
			    center.x += transform[3];
			    center.y += transform[4];
			    center.z += transform[5];
			    fitCenterModCoords[k] = center;
		}

		for(AtomicGaussian refAt:atomicGaussians){
			int index = 0;
			for(AtomicGaussian fitAt:fitVol.atomicGaussians){
				Vtot += refAt.getVolumeOverlap(fitAt, fitCenterModCoords[index],Gaussian3D.DIST_CUTOFF);
				index+=1;	
			}
		}

		if(Vtot<0)
			Vtot = 0.0;
		result[0] = Vtot;
		result[1] = Vvol;
		return result;

	}
		
		
	public double getTotalPPOverlap(double[] transform, ShapeVolume fitVol){
		ExponentialMap eMap = new ExponentialMap(transform[0],transform[1],transform[2]);
		double Vtot = 0.0;
		Coordinates com = fitVol.getCOM();
		double[][] rotMatrix = eMap.toQuaternion().getRotMatrix().getArray();
		List<PPGaussian> fitPPGaussians = fitVol.ppGaussians;
		Coordinates[] fitCenterModCoords = new Coordinates[fitPPGaussians.size()];
		Coordinates[] fitDirectionalityMod = new Coordinates[fitPPGaussians.size()];

		for(int k=0;k<fitPPGaussians.size();k++) {
				Coordinates center = new Coordinates(fitPPGaussians.get(k).getCenter());
				center.sub(com);
			    center.rotate(rotMatrix);
			    center.add(com);
			    center.x += transform[3];
			    center.y += transform[4];
			    center.z += transform[5];
			    fitCenterModCoords[k] = center;
			    fitDirectionalityMod[k] = fitPPGaussians.get(k).getRotatedDirectionality(rotMatrix, 1.0);
		}


		for(PPGaussian refPP:ppGaussians){
			int index = 0;
			for(PPGaussian fitPP:fitVol.ppGaussians){
				Vtot+=refPP.getWeight()*refPP.getSimilarity(fitPP, fitDirectionalityMod[index])* refPP.getVolumeOverlap(fitPP, fitCenterModCoords[index],10.0);
				index+=1;
			
		}
			
		}
		return Vtot;

	}
	
	public void translateToCOM(Coordinates com) {


		for (AtomicGaussian ag : getAtomicGaussians()){
			ag.getCenter().sub(com);  //translate atomicGaussians. Moves center of mass to the origin.
		}

		
		for (PPGaussian pg : getPPGaussians()){
			pg.getCenter().sub(com);  //translate atomicGaussians. Moves center of mass to the origin.
		}
		
		calcCOM();
	}
	
	
	
	

}
