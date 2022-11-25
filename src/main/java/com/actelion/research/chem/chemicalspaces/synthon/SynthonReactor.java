package com.actelion.research.chem.chemicalspaces.synthon;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.MoleculeStandardizer;
import com.actelion.research.chem.SmilesCreator;
import com.actelion.research.chem.StereoMolecule;

public class SynthonReactor {
	
	public static final int CONNECTOR_OFFSET = 92;
	
	public static StereoMolecule react(List<StereoMolecule> synthons) {
		int[] bonds = new int[10];
		List<int[]> deletionMap = new ArrayList<int[]>();
		List<int[]> atomMap = new ArrayList<int[]>();
		
		List<StereoMolecule> buildingBlocks = new ArrayList<StereoMolecule>();
		synthons.stream().forEach(e -> buildingBlocks.add(new StereoMolecule(e)));
		buildingBlocks.forEach(e -> e.ensureHelperArrays(Molecule.cHelperCIP));
		Map<Integer,List<int[]>> rgrps = new HashMap<Integer,List<int[]>>();
		
		for(int m=0;m<buildingBlocks.size();m++) {
			StereoMolecule bb = buildingBlocks.get(m);
			for(int a=0;a<bb.getAtoms();a++) {
				int atomNo = bb.getAtomicNo(a);
				if(atomNo>=CONNECTOR_OFFSET) {
					rgrps.putIfAbsent(atomNo-CONNECTOR_OFFSET+1, new ArrayList<int[]>());
					int connAtom = bb.getConnAtom(a, 0);
					int b = bb.getBond(a, connAtom);
					int upDownIdentifier = bb.getBondAtom(0, b) == a ? -1 : 1; // if up/down if bond is defined from the connected atom, assign 1, otherwise -1
					int bondType = bb.getBondType(b);
					rgrps.get(atomNo-CONNECTOR_OFFSET+1).add(new int[]{m,connAtom,a,upDownIdentifier,bondType});
					int bondOrder = bb.getConnBondOrder(a, 0);
					bonds[0] = bondOrder;
					//bb.markAtomForDeletion(a);
	
					
				}

				
			}
		}
		List<Integer> keysToDelete = new ArrayList<>();
		rgrps.forEach((k,v) -> {
			if(v.size()<2)
				keysToDelete.add(k);
		});
		rgrps.keySet().removeAll(keysToDelete);
		rgrps.forEach((k,v) -> {
			buildingBlocks.get(v.get(0)[0]).markAtomForDeletion(v.get(0)[2]);
			buildingBlocks.get(v.get(1)[0]).markAtomForDeletion(v.get(1)[2]);
			int bb1 = v.get(0)[0];
			int bb2 = v.get(1)[0];
			int u1 = v.get(0)[2];
			int u2 = v.get(1)[2];
			int a1 = v.get(0)[1];
			int a2 = v.get(1)[1];
			alignSynthons(buildingBlocks.get(bb1), buildingBlocks.get(bb2), u1, u2, a1, a2);
		});
		
		// remove R groups
		for(int m=0;m<buildingBlocks.size();m++) {
			StereoMolecule bb = buildingBlocks.get(m);
			int[] map = bb.deleteMarkedAtomsAndBonds();
			deletionMap.add(map);
			final int m_ = m;
			rgrps.forEach((k,v) -> {
				v.stream().forEach(e -> {
					if(e[0]==m_)
						e[1] = map[e[1]];
				});
				
			});

			
		}
		StereoMolecule product = new StereoMolecule();
		atomMap.add(product.addMolecule(buildingBlocks.get(0)));
		for(int m=1;m<buildingBlocks.size();m++) {
			int[] map = product.addMolecule(buildingBlocks.get(m));
			atomMap.add(map);
			final int m_ = m;
			rgrps.forEach((k,v) -> {
				v.stream().forEach(e -> {
					if(e[0]==m_) {
						e[1] = map[e[1]];
						
					}
				});
				
			});
		
		}
		rgrps.forEach((k,v) -> {
			int atom1 = v.get(0)[1];
			int atom2 = v.get(1)[1];
			int bondType = v.get(0)[4];
			int upDownIdentifier = v.get(0)[3];
			int bond=-1;
			if(upDownIdentifier==1)
				bond = product.addBond(atom1, atom2);
			else 
				bond = product.addBond(atom2, atom1);
				
			product.setBondOrder(bond, bonds[k-1]);
			product.setBondType(bond, bondType);
		});
		product.ensureHelperArrays(Molecule.cHelperCIP);
		try {
			MoleculeStandardizer.standardize(product, 0);
			product.ensureHelperArrays(Molecule.cHelperCIP);
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		return product;
		
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
	

	

}
