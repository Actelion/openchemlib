package com.actelion.research.chem.phesa.pharmacophore.pp;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.interactionstatistics.InteractionSimilarityTable;

public class DonorPoint implements IPharmacophorePoint {
	private int donorAtom;
	private int donorHydrogen;
	private Coordinates directionality;
	private int interactionClass;
	private Coordinates center;

	
	public DonorPoint(StereoMolecule mol, int d, int h, int interactionClass) {
		donorAtom = d;
		donorHydrogen = h;
		this.interactionClass = interactionClass;
		updateCoordinates(mol.getAtomCoordinates());
	}
	
	private DonorPoint(String ppString, StereoMolecule mol) {
		decode(ppString,mol);
	}
	
	public DonorPoint(DonorPoint dP) {
		donorAtom = dP.donorAtom;
		donorHydrogen = dP.donorHydrogen;
		directionality = new Coordinates(dP.directionality);
		interactionClass = dP.interactionClass;
		center = new Coordinates(dP.center);
	}
	
	public static DonorPoint fromString(String ppString, StereoMolecule mol) {
		return new DonorPoint(ppString,mol);
	}

	@Override
	public void updateCoordinates(Coordinates[] coords) {
		center = new Coordinates(coords[donorHydrogen].x,coords[donorHydrogen].y,coords[donorHydrogen].z);
		directionality = coords[donorHydrogen].subC(coords[donorAtom]);
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
		updateCoordinates(mol.getAtomCoordinates());
	}

	@Override
	public double getSimilarity(IPharmacophorePoint pp) {
		if(pp instanceof DonorPoint) {
			return 1.0;
		}
		return 0.0;
	}

	public int getInteractionClass() {
		return interactionClass;
	}



	@Override
	public int getCenterID() {
		return donorHydrogen;
	}
	
	@Override
	public void setCenterID(int centerID) {
		donorHydrogen = centerID;
	}

	@Override
	public void setDirectionality(Coordinates directionality) {
		this.directionality = directionality;
		
	}
	
	@Override
	public void updateAtomIndices(int[] map) {
		donorAtom = map[donorAtom];
		donorHydrogen = map[donorHydrogen];
		
	}

	/**
	 *
	 * @return donorAtom and donorHydrogen
	 */
	@Override
	public int[] getAtomIndices() {
		int [] a = {donorAtom, donorHydrogen};
		return a;
	}

	@Override
	public IPharmacophorePoint copyPharmacophorePoint() {
		// TODO Auto-generated method stub
		return new DonorPoint(this);
	}
	
	@Override
	public Coordinates getRotatedDirectionality(double[][] rotMatrix,double scaleFactor) {
		Coordinates directMod = new Coordinates();
		directMod = directionality.rotateC(rotMatrix);
		directMod.scale(scaleFactor); // scale by the invers
		return directMod;
	}

	@Override
	public void getDirectionalityDerivativeCartesian(double[] grad, double[] v, Coordinates di, double sim) {
		grad[3*donorHydrogen] = sim*di.x/3.0;
		grad[3*donorHydrogen+1] = sim*di.y/3.0;
		grad[3*donorHydrogen+2] = sim*di.z/3.0;
		grad[3*donorAtom] = sim*-di.x/3.0;
		grad[3*donorAtom+1] = sim*-di.y/3.0;
		grad[3*donorAtom+2] = sim*-di.z/3.0;
		
	}
	
	@Override
	public int getFunctionalityIndex() {
		return IPharmacophorePoint.Functionality.DONOR.getIndex();
	}
	

	
	
		
	
}
