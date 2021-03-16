package com.actelion.research.chem.phesa;

import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.PeriodicTable;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.util.EncoderFloatingPointNumbers;

public class VolumeGaussian extends Gaussian3D {
	public static final int INCLUSION = 1;
	public static final int EXCLUSION = -1;
	private Coordinates shiftVector;
	private Coordinates referenceVector;
	private int role; //-1 for exclusion, +1 for inclusion
	
	public VolumeGaussian(int atomId,int atomicNo,Coordinates center, Coordinates shiftVector, int role){
		super(atomId,atomicNo,center.addC(shiftVector),1.0);
		this.shiftVector = shiftVector;
		this.referenceVector = center;
		this.role = role;
	}
	
	public VolumeGaussian(VolumeGaussian original){
		super(original.atomId,original.atomicNo,original.center, original.weight);
		this.shiftVector = new Coordinates(original.shiftVector);
		this.referenceVector = new Coordinates(original.referenceVector);
		this.role = original.role;
	}
	
	private VolumeGaussian(String encodedGaussian, StereoMolecule mol) {
		decode(encodedGaussian, mol);
	}
	
	public static VolumeGaussian fromString(String encodedGaussian, StereoMolecule mol) {
		return new VolumeGaussian(encodedGaussian, mol);
	}


	@Override
	public String encode() { //encodes all information of a atomicGaussian using the Base64 encoder
		double[] shift = new double[] {shiftVector.x,shiftVector.y, shiftVector.z};
		StringBuilder molVolString = new StringBuilder();
		molVolString.append(Integer.toString(atomicNo));
		molVolString.append(" ");
		molVolString.append(Integer.toString(atomId));
		molVolString.append(" ");
		molVolString.append(EncoderFloatingPointNumbers.encode(shift,13));
		molVolString.append(" ");
		molVolString.append(Integer.toString(role));
		return molVolString.toString();
	}
	
	public void decode(String string64, StereoMolecule mol)  {
		String[] strings = string64.split(" ");
		atomicNo = Integer.decode(strings[0]);
		atomId = Integer.decode(strings[1]);
		double [] shift = EncoderFloatingPointNumbers.decode(strings[2]);
		role = Integer.decode(strings[3]);
		alpha = calculateWidth(); //the width of the Gaussian depends on the atomic radius of the atom
		volume = calculateVolume();
		coeff = calculateHeight();
		referenceVector = new Coordinates(mol.getAtomX(atomId), mol.getAtomY(atomId), mol.getAtomZ(atomId));
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
		referenceVector = new Coordinates(mol.getCoordinates(atomId));
		center = referenceVector.addC(shiftVector);
	}
	
	public void updateCoordinates(Conformer conf) {
		referenceVector = new Coordinates(conf.getCoordinates(atomId));
		center = referenceVector.addC(shiftVector);
	}
	
	public void setShiftVector(Coordinates shift) {
		this.shiftVector = shift;
	}
	
	public Coordinates getShiftVector() {
		return this.shiftVector;
	}
	
	public void addShift(Coordinates shift) {
		shiftVector.add(shift);
		center = referenceVector.addC(shiftVector);
	}
	
	public void translateRef(Coordinates trans) {
		referenceVector.add(trans);
		center = referenceVector.addC(shiftVector);
	}
	
	public void rotateShift(Matrix rotMat) {
		//PheSAAlignment.rotateCoords(shiftVector, rotMat.getArray());
		shiftVector.rotate(rotMat.getArray());
		center = referenceVector.addC(shiftVector);
	}
	
	public Coordinates getReferenceVector() {
		return this.referenceVector;
	}
	
	public void setReferenceVector(Coordinates referenceVector) {
		this.referenceVector = referenceVector;
	}
	
	public int getRole() {
		return role;
	}
	
	
	

	
	

}
