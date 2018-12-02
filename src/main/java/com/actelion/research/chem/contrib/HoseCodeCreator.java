package com.actelion.research.chem.contrib;

import java.util.Arrays;
import java.util.Vector;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

public class HoseCodeCreator {

	public final static int FULL_HOSE_CODE=0;
	public final static int HOSE_CODE_CUT_C_SP3_SP3=1;
	
	final static boolean DEBUG=false;
	
	/** 
	 * This descriptor requires proper up/down bonds, because it encodes stereo parities. 
	 * If a passed molecule is generated from idcode parsing, make sure that coordinates 
	 * and up/down/bonds are available, i.e. that the IDCodeParser was instantiated with 
	 * the respective option. 
	 */ 
	public static String[][] getHoseCodes(StereoMolecule mol, int maxSphereSize, int type) {
		String[][] ids=new String[mol.getAtoms()][maxSphereSize];
		mol.ensureHelperArrays(Molecule.cHelperRings); 

		for (int rootAtom=0; rootAtom<mol.getAtoms(); rootAtom++) { 
			ids[rootAtom]=getHoseCodesForAtom(mol, rootAtom, maxSphereSize, type);
			
		}
		return ids;
	}
	
	private static boolean isCsp3(ExtendedMolecule mol, int atomID) {
		if (mol.getAtomicNo(atomID)!=6) return false;
		if (mol.getAtomCharge(atomID)!=0) return false;
		if ((mol.getImplicitHydrogens(atomID)+mol.getConnAtoms(atomID))!=4) return false;
		return true;
	}
	
	private static String[] getHoseCodesForAtom(StereoMolecule mol, int rootAtom, int maxSphereSize, int type) {
		StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds()); 
		Vector<String> ids=new Vector();
		int min = 0; 
		int max = 0;
		boolean[] atomMask = new boolean[mol.getAtoms()];
		int[] atomList = new int[mol.getAtoms()]; 
		for (int sphere=0; sphere<maxSphereSize && max<mol.getAtoms(); sphere++) { 
			if (max == 0) { 
				atomList[0] = rootAtom; 
				atomMask[rootAtom] = true; 
				max = 1; 
			} 
			else { 
				int newMax = max; 
				for (int i=min; i<max; i++) { 
					int atom = atomList[i]; 
					for (int j=0; j<mol.getConnAtoms(atom); j++) {
						int connAtom = mol.getConnAtom(atom, j); 
						if (DEBUG) System.out.println("---> "+atom+" to "+connAtom);
						if (!atomMask[connAtom]) {
							switch (type) {
							case FULL_HOSE_CODE:
								atomMask[connAtom] = true; 
								atomList[newMax++] = connAtom;
								break;
							case HOSE_CODE_CUT_C_SP3_SP3:
								if ( ! (isCsp3(mol, atom) && isCsp3(mol,connAtom))) {
									if (DEBUG) System.out.println("NO SKIP");
									atomMask[connAtom] = true; 
									atomList[newMax++] = connAtom;
								} else {
									if (DEBUG) System.out.println("SKIP");
								}
								break;
							}
						} 
					} 
				} 
				min = max; 
				max = newMax; 
			} 

			mol.copyMoleculeByAtoms(fragment, atomMask, true, null); 

			// TO GET ONLY THE SKELETON
			/*
			for (int atom=0; atom<fragment.getAllAtoms(); atom++)  {
				fragment.setAtomicNo(atom, 6); 
			}
			*/
			
			ids.add(new Canonizer(fragment, Canonizer.ENCODE_ATOM_CUSTOM_LABELS).getIDCode()); 
		}
		return ids.toArray(new String[ids.size()]);
	}

	public static String[] getHoseCodesFromDiaID(String diastereotopicID, int maxSphereSize, int type) {
		// We need atom coordinates to properly determine stereo features of fragments later
		StereoMolecule molecule= new IDCodeParser(true).getCompactMolecule(diastereotopicID);
		// One of the atom has to be marked !
		int atomID=-1;
		for (int i=0; i<molecule.getAllAtoms(); i++) {
			// we need to find the marked atom
			String atomCustomLabel=molecule.getAtomCustomLabel(i);
			if (atomCustomLabel!=null && atomCustomLabel.endsWith("*")) {
				atomID=i;
				break;
			};
		}
		if (atomID>=0) {
			return HoseCodeCreator.getHoseCodesForAtom(molecule, atomID, maxSphereSize, type);
		}
		return new String[0];
	}
	
	public static void main(String[] args) {
		StereoMolecule molecule= new IDCodeParser(false).getCompactMolecule("deT@@DjU_k``b`@@");
		StereoMolecule otherMolecule = new StereoMolecule(molecule.getAtoms(), molecule.getBonds()); 
		boolean[] atomMask = new boolean[molecule.getAtoms()];
		Arrays.fill(atomMask, true);
		molecule.copyMoleculeByAtoms(otherMolecule, atomMask, true, null);
		System.out.println(new Canonizer(otherMolecule, Canonizer.ENCODE_ATOM_CUSTOM_LABELS).getIDCode()); 		
		
		// String id="deT@`@f\bbbRK]@PT@@";
		String id="fi{qa@DyZkQPSI`cHhhdhdhddhekF\\\\fNXBBjfjjjaXTh@RB@QJh";
		
		String[] hoses= HoseCodeCreator.getHoseCodesFromDiaID(id, 20, FULL_HOSE_CODE);
		for (int i=0; i<hoses.length; i++) {
			System.out.println(hoses[i]);
		}
		
		hoses= HoseCodeCreator.getHoseCodesFromDiaID(id, 8, HOSE_CODE_CUT_C_SP3_SP3);
		for (int i=0; i<hoses.length; i++) {
			System.out.println(hoses[i]);
		}		
	}
}
