package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;

public class ProteinSynthesizer {
	
	private Molecule3D protein;
	private int terminalC=-1;
	
	public ProteinSynthesizer() {
	}
	
	
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
		if(protein==null) { // first residue
			protein = residue;
			terminalC = newTerminalC;
			coupled = true;
		}
		else if(newTerminalN>-1 && terminalC>-1) { //coupling should be performed
			Coordinates coordsC = protein.getCoordinates(terminalC);
			Coordinates coordsN = residue.getCoordinates(newTerminalN);
			if(coordsC.distanceSquared(coordsN)<Residue.BOND_CUTOFF_SQ) {
				boolean notFound = true;
				for(int i=0;i<protein.getConnAtoms(terminalC) && notFound;i++) {
					int a = protein.getConnAtom(terminalC, i);
					int b = protein.getBond(terminalC, a);
					if(protein.getAtomicNo(a)==8 && protein.getBondOrder(b)==1) {
						notFound=false;
						toDelete = protein.getConnAtom(terminalC, i);
					}
				}
		
				if(toDelete>=0) {
					protein.deleteAtom(toDelete);
					int[] atomMap = protein.addMolecule(residue);
					protein.addBond(terminalC, atomMap[newTerminalN], 1);
					terminalC = protein.getAllAtoms()-(residue.getAllAtoms()-newTerminalC);
					coupled = true;
				}
			}
		}
		protein.ensureHelperArrays(Molecule.cHelperNeighbours);
		return coupled;
	}
	
	public Molecule3D retrieveProtein() {
		return protein;
	}

}
