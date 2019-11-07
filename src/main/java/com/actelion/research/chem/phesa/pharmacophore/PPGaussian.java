package com.actelion.research.chem.phesa.pharmacophore;

import com.actelion.research.util.EncoderFloatingPointNumbers;

import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.List;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.PeriodicTable;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
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
		super(original.atomId,original.atomicNo,original.center,original.weight);
		this.pp = original.pp.copyPharmacophorePoint();
		

	}
	
	private PPGaussian(String encodedGaussian, StereoMolecule mol) {
		decode(encodedGaussian,mol);
	}
	
	public static PPGaussian fromString(String encodedGaussian, StereoMolecule mol) {
		return new PPGaussian(encodedGaussian, mol);
	}
	
	public Coordinates getRotatedDirectionality(double[][] rotMatrix, double scaleFactor) {
		Coordinates direct = pp.getDirectionality();
		Coordinates directMod = new Coordinates();
		directMod.x = direct.x*rotMatrix[0][0] + direct.y*rotMatrix[0][1] + direct.z*rotMatrix[0][2];
		directMod.y = direct.x*rotMatrix[1][0] + direct.y*rotMatrix[1][1] + direct.z*rotMatrix[1][2];
		directMod.z = direct.x*rotMatrix[2][0] + direct.y*rotMatrix[2][1] + direct.z*rotMatrix[2][2];
		//centerModCoords = this.getCenter().rotateC(rotMatrix); //we operate on the transformed coordinates of the molecule to be fitted
		directMod.scale(scaleFactor); // scale by the invers

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
	

	public double getSimilarity(PPGaussian ppGauss2, Coordinates directionality) {
		double ppSimilarity = getInteractionSimilarity(ppGauss2);
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
	public void updateCoordinates(StereoMolecule mol) {
		pp.updateCoordinates(mol);
		center = pp.getCenter();
	}
	
	@Override
	public void updateCoordinates(Conformer conf) {
		pp.updateCoordinates(conf);
		center = pp.getCenter();
	}
	
	@Override
	public void updateAtomIndeces(int[] map) {
		atomId = map[atomId];
		pp.updateAtomIndeces(map);
	}

	@Override
	public double calculateWidth() {
		double vdwR = PeriodicTable.getElement(atomicNo).getVDWRadius();
		return MolecularVolume.alpha_pref/(vdwR*vdwR);
	}
}

