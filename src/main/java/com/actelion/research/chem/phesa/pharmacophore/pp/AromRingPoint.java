package com.actelion.research.chem.phesa.pharmacophore.pp;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.util.ArrayUtils;

public class AromRingPoint implements IPharmacophorePoint {
	private int referenceAtom;
	private List<Integer> ringAtoms;
	private Coordinates center;
	private static final Coordinates directionality = new Coordinates(1.0,0.0,0.0);

	/**
	 *
	 * @param mol
	 * @param a is included in ringAtoms
	 * @param ringAtoms
	 */
	public AromRingPoint(StereoMolecule mol, int a, List<Integer> ringAtoms) {
		referenceAtom = a;
		this.ringAtoms = ringAtoms;
		updateCoordinates(mol.getAtomCoordinates());
	}
	
	public AromRingPoint(AromRingPoint aP) {
		referenceAtom = aP.referenceAtom;
		center = new Coordinates(aP.center);
		ringAtoms = new ArrayList<Integer>();
		for(int ringAtom : aP.ringAtoms) {
			ringAtoms.add(ringAtom);
		}
	}

	@Override
	public Coordinates getCenter() {
		return center;
	}

	@Override
	public void updateCoordinates(Coordinates[] coords) {
		Coordinates com = new Coordinates(0,0,0);
		for(int ringAtom : ringAtoms) {
			com.add(coords[ringAtom]);
		}
		com.scale(1.0/( ringAtoms.size()));

		center = com;
	}
	

	

	@Override
	public Coordinates getDirectionality() {
		// TODO Auto-generated method stub
		return directionality;
	}
	
	private AromRingPoint(String ppString, StereoMolecule mol) {
		decode(ppString,mol);
	}
	
	public static AromRingPoint fromString(String ppString, StereoMolecule mol) {
		return new AromRingPoint(ppString,mol);
	}
	

	private void decode(String ppString, StereoMolecule mol) {
		String[] strings = ppString.split(" ");
		referenceAtom = Integer.decode(strings[1]);
		ringAtoms = new ArrayList<Integer>();
		for(int i=2;i<strings.length;i++) {
			ringAtoms.add(Integer.decode(strings[i]));
		}
		updateCoordinates(mol.getAtomCoordinates());
	}

	@Override
	public String encode() {
		StringBuilder molVolString = new StringBuilder();
		molVolString.append("r");
		molVolString.append(" ");
		molVolString.append(Integer.toString(referenceAtom));
		molVolString.append(" ");
		//molVolString.append(Integer.toString(neighbours.size()));
		//molVolString.append(" ");
		for(Integer ringAtom : ringAtoms) {
			molVolString.append(ringAtom);
			molVolString.append(" ");
		}
		return molVolString.toString().trim();
	}

	@Override
	public double getSimilarity(IPharmacophorePoint pp) {
		double result = 0.0;
		if(pp instanceof AromRingPoint) {
			result = 1.0;
		}
		return result;
	}

	@Override
	public int getCenterID() {
		return referenceAtom;
	}
	
	@Override
	public void setCenterID(int centerID) {
		referenceAtom = centerID;
	}

	@Override
	public void setDirectionality(Coordinates directionality) {
		return;

	}
	

	@Override
	public void updateAtomIndices(int[] map) {
		referenceAtom = map[referenceAtom];
		for(int i=0;i<ringAtoms.size();i++) {
			int neighbour = map[ringAtoms.get(i)];
			ringAtoms.set(i, neighbour);
		}
	}

	@Override
	public int[] getAtomIndices() {
		return ArrayUtils.toIntArray(ringAtoms);
	}

	@Override
	public IPharmacophorePoint copyPharmacophorePoint() {
		// TODO Auto-generated method stub
		return new AromRingPoint(this);
	}

	@Override
	public void getDirectionalityDerivativeCartesian(double[] grad, double[] v, Coordinates di, double sim) {
		return; //no directionality 
		
	}
	
	@Override 
	
	public double getVectorSimilarity(IPharmacophorePoint pp2,Coordinates directionalityMod) {
		return 1.0;
	}
		
	@Override
	 public double getVectorSimilarity(IPharmacophorePoint pp2) {
		return 1.0;
	}
	
	@Override
	public int getFunctionalityIndex() {
		return Functionality.AROM_RING.getIndex();
	}

	@Override
	public Coordinates getRotatedDirectionality(double[][] rotMatrix,double scaleFactor) {
		Coordinates directMod = new Coordinates(directionality);
		return directMod;
	}


}
