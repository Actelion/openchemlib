package com.actelion.research.chem.phesa.pharmacophore.pp;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.PeriodicTable;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.transformation.Transformation;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;
import com.actelion.research.chem.phesa.EncodeFunctions;
import com.actelion.research.chem.phesa.Gaussian3D;
import com.actelion.research.chem.phesa.MolecularVolume;

import java.util.Base64;
import java.util.Base64.Decoder;
import java.util.Base64.Encoder;

/** 
 * @version: 1.0, February 2018
 * Author: J. Wahl
 * basic class to describe Gaussian functions used for the calculation of Molecular Volumes
 * Gaussian functions have a center (3D coordinates), a width and a height 
 * this class provides functionalities for calculating higher order overlaps of Gaussians

*/

public class PPGaussian extends Gaussian3D {
	IPharmacophorePoint pp;

	
	public PPGaussian(int atomicNo,IPharmacophorePoint pp){
		super(pp.getCenterID(),atomicNo,pp.getCenter(), 1.0);
		this.pp = pp;

		
	}
	
	public PPGaussian(PPGaussian original){
		super(original.atomId,original.atomicNo,new Coordinates(original.center),original.weight);
		this.pp = original.pp.copyPharmacophorePoint();
		

	}
	
	private PPGaussian(String encodedGaussian, StereoMolecule mol) {
		decode(encodedGaussian,mol);
	}
	
	public static PPGaussian fromString(String encodedGaussian, StereoMolecule mol) {
		return new PPGaussian(encodedGaussian, mol);
	}
	
	public Coordinates getRotatedDirectionality(double[][] rotMatrix, double scaleFactor) {
		Coordinates directMod = pp.getRotatedDirectionality(rotMatrix,scaleFactor);
		return directMod;
	}
	
	public double getVectorSimilarity(PPGaussian ppGauss2,Coordinates directionalityMod) {
		
		return this.pp.getVectorSimilarity(ppGauss2.getPharmacophorePoint(),directionalityMod);

		
	}

	public double getVectorSimilarity(PPGaussian ppGauss2) {
		return getVectorSimilarity(ppGauss2,ppGauss2.getPharmacophorePoint().getDirectionality());
		
	}
	
	public IPharmacophorePoint getPharmacophorePoint() {
		return pp;
	}
	
	@Override
	public void setAtomId(int atomID) {
		this.pp.setCenterID(atomID);
	}
	
	@Override
	public int getAtomId() {
		return this.pp.getCenterID();
	}

	public double getSimilarity(PPGaussian ppGauss2, Coordinates directionality) {
		double ppSimilarity = 1.0;
		double vectorSim = getVectorSimilarity(ppGauss2,directionality);
		double similarity = (Math.max(0, vectorSim)+2*ppSimilarity)/3.0;
		return similarity;
	}
	
	public double getSimilarity(PPGaussian ppGauss2) {
		
		return getSimilarity(ppGauss2, ppGauss2.getPharmacophorePoint().getDirectionality());
		
	}
	
	public double getInteractionSimilarity(PPGaussian ppGauss2) {

		return pp.getSimilarity(ppGauss2.pp);
	}
	
	@Override 
	
	public void setCenter(Coordinates center) {
		this.center = center;
		this.pp.getCenter().x = center.x;
		this.pp.getCenter().y = center.y;
		this.pp.getCenter().z = center.z;
	}
	


	@Override
	public String encode() { //encodes all information of an atomicGaussian using the Base64 encoder
		Encoder encoder = Base64.getEncoder();
		StringBuilder molVolString = new StringBuilder();
		molVolString.append(Integer.toString(atomicNo));
		molVolString.append(" ");
		molVolString.append(encoder.encodeToString(EncodeFunctions.doubleToByteArray(weight)));
		molVolString.append(" ");
		molVolString.append(pp.encode());

		return molVolString.toString();
	}
	

	public void decode(String string64, StereoMolecule mol)  {
		Decoder decoder = Base64.getDecoder();
		String[] strings = string64.split(" ");
		if(strings.length==1) { // no pharmacophore information encoded
			return;
		}
		atomicNo = Integer.decode(strings[0]);
		weight = EncodeFunctions.byteArrayToDouble(decoder.decode(strings[1].getBytes()));
		StringBuilder sb = new StringBuilder();
		for(int i=2;i<strings.length;i++) {
			sb.append(strings[i]);
			sb.append(" ");
		}
		pp = PharmacophorePointFactory.fromString(sb.toString(), mol);
		center = pp.getCenter();
		alpha = calculateWidth(); //the width of the Gaussian depends on the atomic radius of the atom
		volume = calculateVolume();
		coeff = calculateHeight();
		this.atomId = pp.getCenterID();
	}

	@Override
	public double calculateHeight() {
		return MolecularVolume.p;
	}
	
	@Override 
	public void transform(Transformation transformation) {
		if(!(transformation instanceof TransformationSequence)) {
			this.pp.applyTransformation(transformation);
			this.center = pp.getCenter();
		}
		else {
			TransformationSequence seq = (TransformationSequence) transformation;
			for(Transformation trans : seq.getTransformations()) {
				this.pp.applyTransformation(trans);
				this.center = pp.getCenter();
			}
			}
		}
	
	@Override
	public void updateCoordinates(Coordinates[] coords) {
		pp.updateCoordinates(coords);
		center = pp.getCenter();
	}

	
	@Override
	public void updateAtomIndeces(int[] map) {
		atomId = map[atomId];
		pp.updateAtomIndices(map);
	}

	@Override
	public double calculateWidth() {
		double vdwR = PeriodicTable.getElement(atomicNo).getVDWRadius();
		return MolecularVolume.alpha_pref/(vdwR*vdwR);
	}
}

