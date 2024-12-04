package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;

public class ProteinSynthesizer {
	
	private Molecule3D mProtein;
	private int mTerminalC =-1;
	
	/**
	 * an amino acid is added to the protein structure and a peptide coupling is performed, 
	 * if coupling cannot be performed, the function results false and the residue is not added
	 */
	public boolean addResidue(Molecule3D residue) { 
		boolean coupled = false;
		int toDelete = -1;
		int newTerminalN = -1;
		int newTerminalC = -1;
		for(int atom=0;atom<residue.getAtoms();atom++) {
			if(residue.getAtomicNo(atom)==7 && residue.getAtomCustomLabel(atom)==null)
				newTerminalN = atom;
			if(residue.getAtomicNo(atom)==6 && residue.getAtomCustomLabel(atom)==null)
				newTerminalC = atom;
		}
		if(mProtein ==null) { // first residue
			mProtein = residue;
			mTerminalC = newTerminalC;
			coupled = true;
		}
		else if(newTerminalN>-1 && mTerminalC >-1) { //coupling should be performed
			Coordinates coordsC = mProtein.getCoordinates(mTerminalC);
			Coordinates coordsN = residue.getCoordinates(newTerminalN);
			if(coordsC.distanceSquared(coordsN)<Residue.BOND_CUTOFF_SQ) {
				boolean notFound = true;
				for(int i = 0; i< mProtein.getConnAtoms(mTerminalC) && notFound; i++) {
					int a = mProtein.getConnAtom(mTerminalC, i);
					int b = mProtein.getBond(mTerminalC, a);
					if(mProtein.getAtomicNo(a)==8 && mProtein.getBondOrder(b)==1) {
						notFound=false;
						toDelete = mProtein.getConnAtom(mTerminalC, i);
					}
				}
		
				if(toDelete>=0) {
					mProtein.deleteAtom(toDelete);
					int[] atomMap = mProtein.addMolecule(residue);
					mProtein.addBond(mTerminalC, atomMap[newTerminalN], 1);
					mTerminalC = mProtein.getAllAtoms()-(residue.getAllAtoms()-newTerminalC);
					coupled = true;
				}
			}
		}
		mProtein.ensureHelperArrays(Molecule.cHelperNeighbours);
		return coupled;
	}
	
	public Molecule3D getProtein() {
		return mProtein;
	}
}
