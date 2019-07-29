package com.actelion.research.chem.phesa;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.util.EncoderFloatingPointNumbers;

public class DonorPoint implements IPharmacophorePoint {
	int donorAtom;
	int donorHydrogen;
	Coordinates directionality;
	int interactionClass;
	Coordinates center;

	
	public DonorPoint(StereoMolecule mol, int d, int h, int interactionClass) {
		donorAtom = d;
		donorHydrogen = h;
		this.interactionClass = interactionClass;
		updateCoordinates(mol);
	}
	
	private DonorPoint(String ppString, StereoMolecule mol) {
		decode(ppString,mol);
	}
	
	public static DonorPoint fromString(String ppString, StereoMolecule mol) {
		return new DonorPoint(ppString,mol);
	}

	@Override
	public void updateCoordinates(StereoMolecule mol) {
		center = new Coordinates(mol.getAtomX(donorHydrogen),mol.getAtomY(donorHydrogen),mol.getAtomZ(donorHydrogen));
		directionality = mol.getCoordinates(donorHydrogen).subC(mol.getCoordinates(donorAtom));
		directionality.scale(1.0/directionality.getLength());
		
	}
	

	@Override
	public void updateCoordinates(Conformer conf) {
		center = new Coordinates(conf.getX(donorHydrogen),conf.getY(donorHydrogen),conf.getZ(donorHydrogen));
		directionality = conf.getCoordinates(donorHydrogen).subC(conf.getCoordinates(donorAtom));
		directionality.scale(1.0/directionality.getLength());
		
	}

	@Override
	public Coordinates getCenter() {
		return center;
	}

	@Override
	public Coordinates getDirectionality() {
		return directionality;
	}

	@Override
	public String encode() {
		StringBuilder molVolString = new StringBuilder();
		molVolString.append("d");
		molVolString.append(" ");
		molVolString.append(Integer.toString(donorAtom));
		molVolString.append(" ");
		molVolString.append(Integer.toString(donorHydrogen));
		molVolString.append(" ");
		molVolString.append(Integer.toString(interactionClass));
		return molVolString.toString();
	}
	

	private void decode(String ppString, StereoMolecule mol) {
		String[] strings = ppString.split(" ");
		donorAtom = Integer.decode(strings[1]);
		donorHydrogen = Integer.decode(strings[2]);
		interactionClass = Integer.decode(strings[3]);
		updateCoordinates(mol);
	}

	@Override
	public double getSimilarity(IPharmacophorePoint pp) {
		if(pp instanceof DonorPoint) {
			return 1.0;
		}
		return 0.0;
	}

	@Override
	public int getInteractionClass() {
		return interactionClass;
	}



	@Override
	public int getCenterID() {
		return donorHydrogen;
	}

	@Override
	public void setDirectionality(Coordinates directionality) {
		this.directionality = directionality;
		
	}
	

	
	
		
	
}
