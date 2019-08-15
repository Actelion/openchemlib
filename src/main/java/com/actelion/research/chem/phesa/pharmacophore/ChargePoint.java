package com.actelion.research.chem.phesa.pharmacophore;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;

public class ChargePoint implements IPharmacophorePoint {
	private int chargeAtom;
	private List<Integer> neighbours;
	private int charge;
	private Coordinates center;
	private Coordinates directionality = new Coordinates(1.0,1.0,1.0);
	
	public ChargePoint(StereoMolecule mol, int a, List<Integer> neighbours, int charge) {
		if(charge!=1 && charge!=-1) 
			throw new IllegalArgumentException("charge should be +1 or -1");
		chargeAtom = a;
		this.neighbours = neighbours;
		this.charge = charge;
		updateCoordinates(mol);
	}

	@Override
	public Coordinates getCenter() {
		return center;
	}

	@Override
	public void updateCoordinates(StereoMolecule mol) {
		Coordinates com = new Coordinates(mol.getCoordinates(chargeAtom));
		if(neighbours!=null) {
			for(int neighbour:neighbours) {
				com.add(mol.getCoordinates(neighbour));
			}
			com.scale(1.0/(neighbours.size()+1));
		}

		center = com;
	}
	

	@Override
	public Coordinates getDirectionality() {
		// TODO Auto-generated method stub
		return directionality;
	}
	
	private ChargePoint(String ppString, StereoMolecule mol) {
		decode(ppString,mol);
	}
	
	public static ChargePoint fromString(String ppString, StereoMolecule mol) {
		return new ChargePoint(ppString,mol);
	}
	

	private void decode(String ppString, StereoMolecule mol) {
		String[] strings = ppString.split(" ");
		chargeAtom = Integer.decode(strings[1]);
		charge = Integer.decode(strings[2]);
		neighbours = new ArrayList<Integer>();
		for(int i=3;i<strings.length;i++) {
			neighbours.add(Integer.decode(strings[i]));
		}
		updateCoordinates(mol);
	}

	@Override
	public String encode() {
		StringBuilder molVolString = new StringBuilder();
		molVolString.append("i");
		molVolString.append(" ");
		molVolString.append(Integer.toString(chargeAtom));
		molVolString.append(" ");
		molVolString.append(Integer.toString(charge));
		molVolString.append(" ");
		//molVolString.append(Integer.toString(neighbours.size()));
		//molVolString.append(" ");
		for(Integer neighbour : neighbours) {
			molVolString.append(neighbour);
			molVolString.append(" ");
		}
		return molVolString.toString().trim();
	}

	@Override
	public double getSimilarity(IPharmacophorePoint pp) {
		double result = 0.0;
		if(pp instanceof ChargePoint) {
			result = charge*((ChargePoint)pp).charge > 0 ? 1.0 : 0.0;
		}
		return result;
	}

	@Override
	public int getCenterID() {
		return chargeAtom;
	}

	@Override
	public void setDirectionality(Coordinates directionality) {
		return;

	}
	
	public int getCharge() {
		return charge;
	}

}
