package com.actelion.research.chem.phesa;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.PeriodicTable;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.util.EncoderFloatingPointNumbers;

public class ExclusionGaussian extends Gaussian3D {
	
	private Coordinates shiftVector;
	private Coordinates referenceVector;
	
	public ExclusionGaussian(int atomId,int atomicNo,Coordinates center, Coordinates shiftVector){
		super(atomId,atomicNo,center.addC(shiftVector),1.0);
		this.shiftVector = shiftVector;
		this.referenceVector = center;
	}
	
	public ExclusionGaussian(ExclusionGaussian original){
		super(original.atomId,original.atomicNo,original.center, original.weight);
		this.shiftVector = new Coordinates(original.shiftVector);
		this.referenceVector = new Coordinates(original.referenceVector);
	}
	
	private ExclusionGaussian(String encodedGaussian) {
		decode(encodedGaussian);
	}
	
	public static ExclusionGaussian fromString(String encodedGaussian) {
		return new ExclusionGaussian(encodedGaussian);
	}


	@Override
	public String encode() { //encodes all information of a atomicGaussian using the Base64 encoder
		double[] coords = new double[] {referenceVector.x,referenceVector.y,referenceVector.z};
		double[] shift = new double[] {shiftVector.x,shiftVector.y, shiftVector.z};
		StringBuilder molVolString = new StringBuilder();
		molVolString.append(Integer.toString(atomicNo));
		molVolString.append(" ");
		molVolString.append(Integer.toString(atomId));
		molVolString.append(" ");
		molVolString.append(EncoderFloatingPointNumbers.encode(coords,13));
		molVolString.append(" ");
		molVolString.append(EncoderFloatingPointNumbers.encode(shift,13));
		return molVolString.toString();
	}
	
	public void decode(String string64)  {
		String[] strings = string64.split(" ");
		atomicNo = Integer.decode(strings[0]);
		atomId = Integer.decode(strings[1]);
		double [] coords = EncoderFloatingPointNumbers.decode(strings[2]);
		double [] shift = EncoderFloatingPointNumbers.decode(strings[3]);
		alpha = calculateWidth(); //the width of the Gaussian depends on the atomic radius of the atom
		volume = calculateVolume();
		coeff = calculateHeight();
		referenceVector = new Coordinates(coords[0],coords[1],coords[2]);
		shiftVector = new Coordinates(shift[0],shift[1],shift[2]);
		center = referenceVector.addC(shiftVector);
		weight = 1.0;
		
	}

	@Override
	public double calculateHeight() {
		return MolecularVolume.p;
	}

	@Override
	public double calculateWidth() {
		double vdwR = PeriodicTable.getElement(atomicNo).getVDWRadius();
		return MolecularVolume.alpha_pref/(vdwR*vdwR);
	}
	
	public void updateCoordinates(StereoMolecule mol) {
		referenceVector = mol.getCoordinates(atomId);
		center = referenceVector.addC(shiftVector);
	}
	
	public void updateCoordinates(Conformer conf) {
		referenceVector = conf.getCoordinates(atomId);
		center = referenceVector.addC(shiftVector);
	}
	
	public void setShiftVector(Coordinates shift) {
		this.shiftVector = shift;
	}
	
	public Coordinates getShiftVector() {
		return this.shiftVector;
	}
	
	public Coordinates getReferenceVector() {
		return this.referenceVector;
	}
	
	public void setReferenceVector(Coordinates referenceVector) {
		this.referenceVector = referenceVector;
	}
	
	
	
	

	
	

}
