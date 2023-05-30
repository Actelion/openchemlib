package com.actelion.research.chem.chemicalspaces.synthon;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.MoleculeStandardizer;
import com.actelion.research.chem.SmilesCreator;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.coords.CoordinateInventor;

public class SynthonReactor {
	
	public static final int CONNECTOR_OFFSET = 92;
	
	public static StereoMolecule react(List<StereoMolecule> synthons) {		
		List<StereoMolecule> buildingBlocks = new ArrayList<StereoMolecule>();
		synthons.stream().forEach(e -> buildingBlocks.add(new StereoMolecule(e)));
		buildingBlocks.forEach(e -> e.ensureHelperArrays(Molecule.cHelperCIP));
		Set<Integer> rgrps = new HashSet<Integer>();
		for(int m=0;m<buildingBlocks.size();m++) {
			StereoMolecule bb = buildingBlocks.get(m);
			new CoordinateInventor().invent(bb);
			for(int a=0;a<bb.getAtoms();a++) {
				int atomNo = bb.getAtomicNo(a);
				if(atomNo>=CONNECTOR_OFFSET) {
					rgrps.add(atomNo-CONNECTOR_OFFSET+1);

				}

				
			}
		}
		StereoMolecule prod = null;
		for(int i : rgrps) {
			for(int j=0;j<buildingBlocks.size();j++) {
				StereoMolecule bb1 = buildingBlocks.get(j);
				if(bb1.getAllAtoms()==0)
					continue;
				boolean reacted = false;
				for(int connAtom1 : findConnectorAtoms(i,bb1)) {
					if(reacted)
						continue;
					if(connAtom1>-1) {
						SynthonConnector sc1 = new SynthonConnector();
						sc1.connAtom = connAtom1;
						sc1.bb = bb1;
						sc1.neighbourAtom = bb1.getConnAtom(connAtom1, 0);
						int b1 = bb1.getBond(connAtom1, sc1.neighbourAtom);
						sc1.bondType = bb1.getBondType(b1);
						sc1.bond = b1;
						sc1.bondOrder = bb1.getBondOrder(b1);
						for(int k=0;k<buildingBlocks.size();k++) {
							StereoMolecule bb2 = buildingBlocks.get(k);
							if(bb2.getAllAtoms()==0)
								continue;
							for(int connAtom2 : findConnectorAtoms(i,bb2)) {
								if(reacted)
									continue;
								if(j==k && connAtom1==connAtom2)
									continue;
								SynthonConnector sc2 = new SynthonConnector();
								sc2.connAtom = connAtom2;
								sc2.bb = bb2;
								sc2.neighbourAtom = bb2.getConnAtom(connAtom2, 0);
								int b2 = bb2.getBond(connAtom2, sc2.neighbourAtom );
								sc2.bondType = bb2.getBondType(b2);
								sc2.bond = b2;
								sc2.bondOrder = bb2.getBondOrder(b2);
								prod = combineSynthons(sc1,sc2);
								reacted = true;
						}
					}
				}
			}
			}
		}
		return prod;
		
	}
	

		
	private static StereoMolecule combineSynthons(SynthonConnector sc1, SynthonConnector sc2) {
		StereoMolecule mol = null;
		SynthonConnector referenceConnector = null;
		SynthonConnector connector = null;
		if(sc2.bondType==Molecule.cBondTypeUp || sc2.bondType==Molecule.cBondTypeDown ) {
			referenceConnector = sc2;
			connector = sc1;
		}
		else {
			referenceConnector = sc1;
			connector = sc2;
		}

		int upDownIdentifier = referenceConnector.bb.getBondAtom(0, referenceConnector.bond) == referenceConnector.connAtom ? -1 : 1; // if up/down if bond is defined from the connected atom, assign 1, otherwise -1
		int u1 = referenceConnector.connAtom;
		int u2 = connector.connAtom;
		int a1 = referenceConnector.neighbourAtom;
		int a2 = connector.neighbourAtom;
		int newConnectorAtom = -1;
		if(referenceConnector.bb == connector.bb) {
			referenceConnector.bb.markAtomForDeletion(referenceConnector.connAtom);
			referenceConnector.bb.markAtomForDeletion(connector.connAtom);
			int[] map = referenceConnector.bb.deleteMarkedAtomsAndBonds();	
			referenceConnector.neighbourAtom = map[referenceConnector.neighbourAtom];
			connector.neighbourAtom = map[connector.neighbourAtom];
			newConnectorAtom = connector.neighbourAtom;
			
		}
		else {
			referenceConnector.bb.markAtomForDeletion(u1);
			connector.bb.markAtomForDeletion(u2);
			alignSynthons(referenceConnector.bb, connector.bb,u1,u2, a1, a2);
	
			int[] map = referenceConnector.bb.deleteMarkedAtomsAndBonds();	
			referenceConnector.neighbourAtom = map[referenceConnector.neighbourAtom];
			map = connector.bb.deleteMarkedAtomsAndBonds();	
			connector.neighbourAtom = map[connector.neighbourAtom];
			map = referenceConnector.bb.addMolecule(connector.bb);
			newConnectorAtom = map[connector.neighbourAtom];
		}
		int bond = -1;
		if(upDownIdentifier==-1) {
			bond = referenceConnector.bb.addBond(newConnectorAtom, referenceConnector.neighbourAtom);
		}
		else {
			bond = referenceConnector.bb.addBond(referenceConnector.neighbourAtom,newConnectorAtom);
		}
		referenceConnector.bb.setBondOrder(bond, referenceConnector.bondOrder);
		referenceConnector.bb.setBondType(bond, referenceConnector.bondType);

		referenceConnector.bb.ensureHelperArrays(Molecule.cHelperCIP);

		new CoordinateInventor().invent(referenceConnector.bb);

		try {
			MoleculeStandardizer.standardize(referenceConnector.bb, 0);
			referenceConnector.bb.ensureHelperArrays(Molecule.cHelperCIP);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		if(referenceConnector.bb != connector.bb) {
			connector.bb.clear();
		}
		else {
			connector.bb = new StereoMolecule();
		}
		mol = referenceConnector.bb;
		return mol;
	}
	
	private static List<Integer> findConnectorAtoms(int rGrp, StereoMolecule bb) {
		List<Integer> atoms = new ArrayList<>();
		for(int a=0;a<bb.getAtoms();a++) {
			int atomNo = bb.getAtomicNo(a);
			if(atomNo>=CONNECTOR_OFFSET) {
				if(atomNo-CONNECTOR_OFFSET+1==rGrp) {
					atoms.add(a);
				}
			}
		}
		return atoms;
	}
	
	/**
	 * align synthons in 2D so that they are in the correct orientation for the formation of the new bond
	 * consider the bonds A1---U1 in Synthon1 and A2---U2 in Synthon2, the two bond vectors should be aligned antiparallel
	 * and A2 should be translated to U1, this is important for the correct assignment of the stereochemistry, which should be
	 * preserved for the synthon reaction
	 * @param s1
	 * @param s2
	 * @param u1
	 * @param u2
	 * @param a1
	 * @param a2
	 */
	public static void alignSynthons(StereoMolecule s1, StereoMolecule s2,int u1, int u2, int a1, int a2) {
		Coordinates v1 = s1.getCoordinates(u1).subC(s1.getCoordinates(a1));
		Coordinates v2 = s2.getCoordinates(a2).subC(s2.getCoordinates(u2));
		v1.unit();
		v2.unit();

		double alpha = Math.acos(v1.dot(v2));

		Coordinates n;
		Coordinates cross = v1.cross(v2);
		if(cross.dist()<0.001) {
			n = new Coordinates(0.0,0.0,1.0);
			alpha = Math.PI;
		}
		else {
			n = cross.unit();
		}
			

		//Rodrigues' rotation formula 

		Coordinates t = s1.getCoordinates(a1).scaleC(-1.0);
		for(int a=0;a<s1.getAtoms();a++) {
			s1.getCoordinates(a).add(t);
		}
		t = s2.getCoordinates(u2).scaleC(-1.0);
		for(int a=0;a<s2.getAtoms();a++) {
			s2.getCoordinates(a).add(t);
		}

		for(int a=0;a<s2.getAtoms();a++) {
			Coordinates newCoords = eulerRodrigues(s2.getCoordinates(a),n,-alpha);
			s2.setAtomX(a,newCoords.x);
			s2.setAtomY(a,newCoords.y);
			s2.setAtomZ(a,newCoords.z);
		}


	}
	
	public static Coordinates eulerRodrigues(Coordinates v, Coordinates k, double theta) {
		Coordinates c1 = v.scaleC(Math.cos(theta));
		Coordinates c2 = k.cross(v).scale(Math.sin(theta));
		Coordinates c3 = k.scaleC(k.dot(v)*(1-Math.cos(theta)));
		Coordinates vNew = c1.addC(c2);
		vNew.add(c3);
		return vNew;
	}
	
	private static class SynthonConnector {
		StereoMolecule bb;
		int bond;
		int bondOrder;
		int connAtom;
		int neighbourAtom;
		int bondType;
		


	}
	

	

}
