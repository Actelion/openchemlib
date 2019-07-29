package com.actelion.research.chem.phesa;

import java.util.ArrayList;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;


public class AcceptorPoint implements IPharmacophorePoint {
	int acceptorAtom;
	ArrayList<Integer> neighbours;
	Coordinates directionality;
	int interactionClass;
	Coordinates center;
	
	public AcceptorPoint(StereoMolecule mol, int a, ArrayList<Integer> neighbours, int interactionClass) {
		acceptorAtom = a;
		this.neighbours = neighbours;
		this.interactionClass = interactionClass;
		updateCoordinates(mol);
	}
	
	private AcceptorPoint(String ppString, StereoMolecule mol) {
		decode(ppString,mol);
	}
	
	public static AcceptorPoint fromString(String ppString, StereoMolecule mol) {
		return new AcceptorPoint(ppString,mol);
	}

	@Override
	public void updateCoordinates(StereoMolecule mol) {
		center = new Coordinates(mol.getAtomX(acceptorAtom),mol.getAtomY(acceptorAtom),mol.getAtomZ(acceptorAtom));
		if(neighbours.size()==1) {
			int aa1 = neighbours.get(0);
			directionality = center.subC(mol.getCoordinates(aa1));}
		else {
			int aa1 = neighbours.get(0);
			int aa2 = neighbours.get(1);
			Coordinates v1 = center.subC(mol.getCoordinates(aa1));
			Coordinates v2 = center.subC(mol.getCoordinates(aa2));
			directionality = v1.addC(v2);
		}
		directionality.scale(1.0/directionality.getLength());
		// TODO Auto-generated method stub
		
	}
	
	@Override
	public void updateCoordinates(Conformer conf) {
		center = new Coordinates(conf.getX(acceptorAtom),conf.getY(acceptorAtom),conf.getZ(acceptorAtom));
		if(neighbours.size()==1) {
			int aa1 = neighbours.get(0);
			directionality = center.subC(conf.getCoordinates(aa1));}
		else {
			int aa1 = neighbours.get(0);
			int aa2 = neighbours.get(1);
			Coordinates v1 = center.subC(conf.getCoordinates(aa1));
			Coordinates v2 = center.subC(conf.getCoordinates(aa2));
			directionality = v1.addC(v2);
		}
		directionality.scale(1.0/directionality.getLength());
		// TODO Auto-generated method stub
		
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
		molVolString.append("a");
		molVolString.append(" ");
		molVolString.append(Integer.toString(acceptorAtom));
		molVolString.append(" ");
		molVolString.append(Integer.toString(interactionClass));
		molVolString.append(" ");
		//molVolString.append(Integer.toString(neighbours.size()));
		//molVolString.append(" ");
		for(Integer neighbour : neighbours) {
			molVolString.append(neighbour);
			molVolString.append(" ");
		}
		return molVolString.toString().trim();
	}
	

	private void decode(String ppString, StereoMolecule mol) {
		String[] strings = ppString.split(" ");
		acceptorAtom = Integer.decode(strings[1]);
		interactionClass = Integer.decode(strings[2]);
		neighbours = new ArrayList<Integer>();
		for(int i=3;i<strings.length;i++) {
			neighbours.add(Integer.decode(strings[i]));
		}
		updateCoordinates(mol);
	}
	
	@Override
	public double getSimilarity(IPharmacophorePoint pp) {
		if(pp instanceof AcceptorPoint) {
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
		return acceptorAtom;
	}
	
	@Override
	public void setDirectionality(Coordinates directionality) {
		this.directionality = directionality;
		
	}


	
		
	
}

