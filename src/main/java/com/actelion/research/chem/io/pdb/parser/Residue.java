/*
 * Copyright (c) 1997 - 2016
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the the copyright holder nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Modest v. Korff
 */

package com.actelion.research.chem.io.pdb.parser;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.io.pdb.calc.AminoAcids;
import com.actelion.research.chem.io.pdb.calc.BondOrderCalculator;
import com.actelion.research.chem.io.pdb.calc.BondsCalculator;

import java.util.*;

/**
 * Residue handles a group of atoms, e.g. a Ligand Molecule or an AminoAcid residue
 *
 * Created by korffmo1 on 11.04.18.
 */
public class Residue {
	public static double BOND_CUTOFF_SQ = 3.24;
	private boolean mIsTerminal;
	private final boolean mIsProtein;
	private final boolean mAddBonds,mIncludeAltLocs;
    private Molecule3D mMol;
	private final int mResNum;
	private final String mResName,mChainID,mInsertionCode;
	private final String mFragmentName;
	private Map<String,AtomRecord> mAtomRecordMap;
	private ArrayList<AtomRecord> mAtomList;

	/**
	 * @param addBonds
	 * @param includeAltLocations
	 */
	public Residue(boolean isProtein, boolean addBonds, boolean includeAltLocations, AtomRecord firstAtom) {
		mIsProtein = isProtein;
		mAddBonds = addBonds;
		mIncludeAltLocs = includeAltLocations;

		mResNum = firstAtom.getAuthSeqID();
		mResName = firstAtom.getResName();
		mChainID = firstAtom.getChainID();
		mInsertionCode = firstAtom.getInsertionCode();

		// We use '1' as Model-No.
		mFragmentName = firstAtom.getResName()+"_1_"+firstAtom.getChainID()+"_"+firstAtom.getAuthSeqID();

		if (mIsProtein)
			mAtomRecordMap = new HashMap<>();
		else
			mAtomList = new ArrayList<>();

		addAtom(firstAtom);
	}

	public void addAtom(AtomRecord atom) {
		if (mIsProtein)
			mAtomRecordMap.put(atom.getLabelAtomName(), atom);
		else mAtomList.add(atom);

		if(atom.isTerminalC())
			mIsTerminal = true;
	}

	public void build() {
		mMol = mIsProtein ? AminoAcids.createResidue(getResname(), mAtomRecordMap) : null;
		if (mMol == null)
			mMol = constructFragmentFromGeometry(mIsProtein ? mAtomRecordMap.values() : mAtomList);
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
		mMol.setName(mFragmentName);
	}

	public boolean isTerminal() {
		return mIsTerminal;
	}

    private void constructFragment(List<AtomRecord> records) {
		AtomRecord ar = records.get(0);
		mMol = AminoAcids.createResidue(getResname(), records);
		if (mMol == null)
			mMol = constructFragmentFromGeometry(records);
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
		// We use '1' as Model-No.
		mMol.setName(ar.getResName()+"_1_"+ar.getChainID()+"_"+ar.getAuthSeqID());
    }

	private Molecule3D constructFragmentFromGeometry(Collection<AtomRecord> records) {
		Molecule3D fragment = new Molecule3D();
		String altLoc = null;
		for(AtomRecord record : records) {
			if (!mIncludeAltLocs) {
				if (altLoc == null)
					altLoc = record.getLabelAltID();
				else if (!altLoc.equals(record.getLabelAltID()))
					continue;
			}
			int atomicNo = record.getAtomicNo();
			int atom = fragment.addAtom(atomicNo);
			fragment.setAtomName(atom, record.getLabelAtomName());
			fragment.setAtomAmino(atom, record.getResName());
			fragment.setAtomSequence(atom,record.getSerialId());
			fragment.setResSequence(atom, record.getAuthSeqID());
			fragment.setAtomAmino(atom, record.getResName());
			fragment.setAtomChainId(atom, record.getChainID());
			fragment.setAtomX(atom,record.getX());
			fragment.setAtomY(atom,record.getY());
			fragment.setAtomZ(atom,record.getZ());
			AtomRecord bridgeRecord = record.getCovalentBridgeAtom();
			if (bridgeRecord != null) {
				int bridgeAtom = fragment.addAtom(0);
				fragment.setAtomCustomLabel(bridgeAtom,"]cov");
				fragment.setAtomX(bridgeAtom, (bridgeRecord.getX() + record.getX()) / 2);
				fragment.setAtomY(bridgeAtom, (bridgeRecord.getY() + record.getY()) / 2);
				fragment.setAtomZ(bridgeAtom, (bridgeRecord.getZ() + record.getZ()) / 2);
				fragment.addBond(atom, bridgeAtom, Molecule.cBondTypeSingle);
				fragment.setCovalentLigand(true);
			}
		}

		if (mAddBonds) {
			try {
				BondsCalculator.createBonds(fragment, true,null);
				new BondOrderCalculator(fragment).calculateBondOrders();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		return fragment;
	}
	
	public int getResnum() {
		return mResNum;
	}
	
	public Molecule3D getMolecule() {
		return mMol;
	}
	
	public String getResname() {
		return mResName;
	}
	
	public String getChainID() {
		return mChainID;
	}
	
	public String getInsertionCode() {
		return mInsertionCode;
	}
	
	/*
	 * check if two residues are connected by a peptide bond, this is the residue with the (supposedly) terminal nitrogen,
	 * resC should contain the terminal carbon
	 *
	public boolean areBonded(Residue resN) { 
		AtomRecord recordC = null; 
		AtomRecord recordN = null;
		boolean areBonded = false;
		for(AtomRecord rec : this.mAtomRecordList) {
			if (rec.getAtomName().equals("C")) {
				recordN = rec;
				break;
			}
		}
		
		for(AtomRecord rec : resN.mAtomRecordList) {
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
	}*/
}
