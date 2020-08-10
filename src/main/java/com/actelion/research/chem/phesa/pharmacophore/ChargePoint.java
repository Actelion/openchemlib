package com.actelion.research.chem.phesa.pharmacophore;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

public class ChargePoint implements IPharmacophorePoint {
	private int chargeAtom;
	private List<Integer> neighbours;
	private int charge;
	private Coordinates center;
	private Coordinates directionality = new Coordinates(0.0,0.0,0.0);
	
	public ChargePoint(StereoMolecule mol, int a, List<Integer> neighbours, int charge) {
		if(charge!=1 && charge!=-1) 
			throw new IllegalArgumentException("charge should be +1 or -1");
		chargeAtom = a;
		this.neighbours = neighbours;
		this.charge = charge;
		updateCoordinates(mol);
	}
	
	public ChargePoint(ChargePoint cP) {
		chargeAtom = cP.chargeAtom;
		charge = cP.charge;
		directionality = new Coordinates(cP.directionality);
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
	public void updateCoordinates(Conformer conf) {
		Coordinates com = new Coordinates(conf.getCoordinates(chargeAtom));
		if(neighbours!=null) {
			for(int neighbour:neighbours) {
				com.add(conf.getCoordinates(neighbour));
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
	
	@Override
	public void updateAtomIndeces(int[] map) {
		chargeAtom = map[chargeAtom];

		for(int i=0;i<neighbours.size();i++) {
			int neighbour = map[neighbours.get(i)];
			neighbours.set(i, neighbour);
		}

		
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
