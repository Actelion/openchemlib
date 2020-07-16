package com.actelion.research.chem.io.pdb.parser;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.io.pdb.converter.AminoAcidsLabeledContainer;
import com.actelion.research.chem.io.pdb.converter.BondsCalculator;

/**
 * ModelGroupAtoms
 * handles a group of atoms, e.g. a Ligand Molecule or an AminoAcid residue
 *
 *
 * <p>Copyright: Idorsia Pharmaceuticals Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Idorsia Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * Created by korffmo1 on 11.04.18.
 */
public class Residue {
	public static double BOND_CUTOFF_SQ = 3.24;
    private List<AtomRecord> records;
    private boolean isTerminal;
    private Molecule3D mol;

    public Residue(List<AtomRecord> records) {
        this.records = records;
        this.isTerminal = isTerminal();
        mol = constructFragment();
    }
        

	public boolean isTerminal() {
		isTerminal = false;
		for(AtomRecord record : records) {
			if(record.isTerminalC())
				isTerminal = true;
		}
		return isTerminal;
	}
    

    
    private Molecule3D constructFragment() {
    	boolean isAA = false;
    	if(AminoAcidsLabeledContainer.INSTANCE.get(getResname())!=null)
    		isAA = true;
    	Molecule3D fragment;
    	if(isAA) 
    		fragment = constructFragmentFromIDCode(getResname());
    	else 
    		fragment = constructFragmentFromGeometry(getResname());
    	fragment.ensureHelperArrays(Molecule.cHelperNeighbours);
    	return fragment;
    }
    
    private Molecule3D constructFragmentFromIDCode(String resname) {
		Molecule3D fragment;
		Map<String,AtomRecord> recordMap = new HashMap<String,AtomRecord>();
		for(AtomRecord record : records) {
			String atomName = record.getAtomName();
			recordMap.put(atomName,record);
		}
		fragment = AminoAcidsLabeledContainer.INSTANCE.get(resname).createResidue(recordMap);
		if(fragment==null) {
			fragment = constructFragmentFromGeometry(resname);
		}
		return fragment;
    }
    
	private Molecule3D constructFragmentFromGeometry(String resname) {
		Molecule3D fragment = new Molecule3D();
		for(AtomRecord record : records) {
			int atomicNo = record.getAtomicNo();
			int atom = fragment.addAtom(atomicNo);
			fragment.setAtomName(atom, record.getAtomName());
			fragment.setAtomAmino(atom, record.getResName());
			fragment.setAtomSequence(atom,record.getSerialId());
			fragment.setResSequence(atom, record.getResNum());
			fragment.setAtomAmino(atom, record.getResName());
			fragment.setAtomChainId(atom, record.getChainID());
			fragment.setAtomX(atom,record.getX());
			fragment.setAtomY(atom,record.getY());
			fragment.setAtomZ(atom,record.getZ());
		}
		try {
			BondsCalculator.createBonds(fragment, true);
			BondsCalculator.calculateBondOrders(fragment);
		} catch (Exception e) {
			System.err.println("Cannot process structure");
		}
		fragment.ensureHelperArrays(Molecule.cHelperCIP);
		return fragment;
	}
	
	public int getResnum() {
		return records.get(0).getResNum();
	}
	
	public Molecule3D getMolecule() {
		return mol;
	}
	
	public String getResname() {
		return records.get(0).getResName();
	}
	
	public String getChainID() {
		return records.get(0).getChainID();
	}
	
	public String getInsertionCode() {
		return records.get(0).getInsertionCode();
	}
	
	
	/*
	 * check if two residues are connected by a peptide bond, this is the residue with the (supposedly) terminal nitrogen,
	 * resC should contain the terminal carbon
	 */
	public boolean areBonded(Residue resN) { 
		AtomRecord recordC = null; 
		AtomRecord recordN = null;
		boolean areBonded = false;
		for(AtomRecord rec : this.records) {
			if (rec.getAtomName().equals("C")) {
				recordN = rec;
				break;
			}
		}
		
		for(AtomRecord rec : resN.records) {
			if (rec.getAtomName().equals("N")) {
				recordC = rec;
				break;
			}
		}
		if(recordN!=null && recordC!=null) {
			double dx = recordN.getX()-recordC.getX();
			double dy = recordN.getY()-recordC.getY();
			double dz = recordN.getZ()-recordC.getZ();
			double distSq = dx*dx + dy*dy + dz*dz;
			if(distSq<BOND_CUTOFF_SQ)
				areBonded=true;
		}
		return areBonded;
		
	}
		
		

}
