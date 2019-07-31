package com.actelion.research.chem.phesa;

import java.util.ArrayList;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;

public class AcceptorPoint implements IPharmacophorePoint {
	private int acceptorAtom;
	private ArrayList<Integer> neighbours;
	private Coordinates directionality;
	private int interactionClass;
	private Coordinates center;
	private int acceptorID; //necessary to assign different directionalities to two acceptor points in sp2 oxygen
	
	public AcceptorPoint(StereoMolecule mol, int a, ArrayList<Integer> neighbours, int interactionClass) {
		this(mol, a, neighbours, interactionClass, 0);
	}
	
	public AcceptorPoint(StereoMolecule mol, int a, ArrayList<Integer> neighbours, int interactionClass, int acceptorID) {
		acceptorAtom = a;
		this.neighbours = neighbours;
		this.interactionClass = interactionClass;
		this.acceptorID = acceptorID;
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
			directionality = center.subC(mol.getCoordinates(aa1));
		}
			
			
		else if(neighbours.size()==2 && acceptorID!=0) {
			int aa1 = neighbours.get(0);
			Coordinates v1 = center.subC(mol.getCoordinates(aa1));
			int aa2 = neighbours.get(1);
			Coordinates v2 = mol.getCoordinates(aa2).subC(center);
			Coordinates rotAxis = v1.cross(v2).unit();
			double theta = acceptorID == 1 ? 45.0/180.0*Math.PI :  -45.0/180.0*Math.PI;
			directionality = v1.rotate(rotAxis, theta);
		}
		
		
		else if(neighbours.size()==3) {
			int aa1 = neighbours.get(0);
			int aa2 = neighbours.get(1);
			int aa3 = neighbours.get(2);
			Coordinates v1 = center.subC(mol.getCoordinates(aa1)).unit();
			Coordinates v2 = center.subC(mol.getCoordinates(aa2)).unit();
			Coordinates v3 = center.subC(mol.getCoordinates(aa3)).unit();
			directionality = v3.add(v2).add(v1);
		}
		
		
		else {
			int aa1 = neighbours.get(0);
			int aa2 = neighbours.get(1);
			Coordinates v1 = center.subC(mol.getCoordinates(aa1)).unit();
			Coordinates v2 = center.subC(mol.getCoordinates(aa2)).unit();
			directionality = v1.addC(v2);
		}
		directionality.unit();
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
		molVolString.append(Integer.toString(acceptorID));
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
		acceptorID = Integer.decode(strings[3]);
		neighbours = new ArrayList<Integer>();
		for(int i=4;i<strings.length;i++) {
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
	
	public int getAcceptorID() {
		return acceptorID;
	}


	
		
	
}

