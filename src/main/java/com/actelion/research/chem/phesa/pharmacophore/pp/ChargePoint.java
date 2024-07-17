package com.actelion.research.chem.phesa.pharmacophore.pp;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.ArrayUtils;

public class ChargePoint implements IPharmacophorePoint {
	private int chargeAtom;
	private List<Integer> neighbours;
	private int charge;
	private Coordinates center;
	private static final Coordinates directionality = new Coordinates(1.0,0.0,0.0);
	
	public ChargePoint(StereoMolecule mol, int a, List<Integer> neighbours, int charge) {
		if(charge==0) 
			throw new IllegalArgumentException("charge should not be 0 ");
		chargeAtom = a;
		this.neighbours = neighbours;
		this.charge = charge;
		updateCoordinates(mol.getAtomCoordinates());
	}
	
	public ChargePoint(ChargePoint cP) {
		chargeAtom = cP.chargeAtom;
		charge = cP.charge;
		center = new Coordinates(cP.center);
		neighbours = new ArrayList<Integer>();
		for(int neighbour : cP.neighbours) {
			neighbours.add(neighbour);
		}
	}

	@Override
	public Coordinates getCenter() {
		return center;
	}

	@Override
	public void updateCoordinates(Coordinates[] coords) {
		Coordinates com = new Coordinates(coords[chargeAtom]);
		if(neighbours!=null) {
			for(int neighbour:neighbours) {
				com.add(coords[neighbour]);
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
	
	@Override
	public Coordinates getRotatedDirectionality(double[][] rotMatrix,double scaleFactor) {
		Coordinates directMod = new Coordinates(directionality);
		return directMod;
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
		updateCoordinates(mol.getAtomCoordinates());
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
	public void setCenterID(int centerID) {
		chargeAtom = centerID;
	}

	@Override
	public void setDirectionality(Coordinates directionality) {
		return;

	}
	
	public int getCharge() {
		return charge;
	}
	
	@Override
	public void updateAtomIndices(int[] map) {
		chargeAtom = map[chargeAtom];

		for(int i=0;i<neighbours.size();i++) {
			int neighbour = map[neighbours.get(i)];
			neighbours.set(i, neighbour);
		}
	}

	@Override
	public int[] getAtomIndices() {
		int [] a = {chargeAtom};
		return a;
	}


	@Override
	public IPharmacophorePoint copyPharmacophorePoint() {
		// TODO Auto-generated method stub
		return new ChargePoint(this);
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
		if(charge<0)
			return IPharmacophorePoint.Functionality.NEG_CHARGE.getIndex();
		else
			return IPharmacophorePoint.Functionality.POS_CHARGE.getIndex();
	}

	public int getChargeAtom() {
		return chargeAtom;
	}

	public void setChargeAtom(int chargeAtom) {
		this.chargeAtom = chargeAtom;
	}

	public List<Integer> getNeighbours() {
		return neighbours;
	}

	public void setNeighbours(List<Integer> neighbours) {
		this.neighbours = neighbours;
	}

}
