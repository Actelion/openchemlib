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
 */

package com.actelion.research.chem;

import java.util.*;

public class Molecule3D extends StereoMolecule implements Comparable<Molecule3D> {
	public static final int RIGID					= 0x00000001;
	public static final int LIGAND					= 0x00000002;
	public static final int BACKBONE				= 0x00000004;
	public static final int FLAG1					= 0x00000008; //Flag used for different purpose
	public static final int IMPORTANT				= 0x00000010;
	public static final int PREOPTIMIZED			= 0x00000020;
	public static final int ATTACHED_HYDROGEN_COUNT = 0x000003C0; // explicitly known hydrogens without coordinates
	public static final int ATTACHED_HYDROGEN_COUNT_SHIFT = 6;

	public static final int INFO_DESCRIPTION = 0;
	public static final int INFO_ATOMSEQUENCE = 1;
	public static final int INFO_INTERACTION_ATOM_TYPE = 2;
	public static final int INFO_ATOMNAME = 3;
	public static final int INFO_AMINO = 4;
	public static final int INFO_PPP = 5;
	public static final int INFO_CHAINID = 6;
	public static final int INFO_BFACTOR = 7;
	public static final int INFO_RESSEQUENCE = 8;

	private static final int MAX_INFOS = 9;

	//Molecule information
	private int nMovables = -1;
	private final Map<String, Object> auxiliaryInfos = new HashMap<String, Object>();

	//Atom information
	private int[] atomFlags;
	private Object[][] infos;
	private double[] partialCharges;

	//Molecule properties (calculated during the first call at getRings)
	//private List<Integer>[] mAtomToRings = null;


	public Molecule3D(Molecule3D mol) {
		super(mol);

		auxiliaryInfos.putAll(mol.auxiliaryInfos);

		atomFlags = new int[mol.getMaxAtoms()];
		partialCharges = new double[mol.getMaxAtoms()];
		infos = new Object[mol.getMaxAtoms()][MAX_INFOS];

		for(int atom=0; atom<mol.getAllAtoms(); atom++) {
			atomFlags[atom] = mol.atomFlags[atom];
			partialCharges[atom] = mol.partialCharges[atom];
			for (int i=0; i<MAX_INFOS; i++)
				infos[atom][i] = mol.infos[atom][i];
		}
	}
	
	public Molecule3D(StereoMolecule mol) {
		super(mol);
		atomFlags = new int[mol.getMaxAtoms()];
		partialCharges = new double[mol.getMaxAtoms()];
		infos = new Object[mol.getMaxAtoms()][MAX_INFOS];
	}

	public Molecule3D() {
		this(5, 5);
	}

	public Molecule3D(int a, int b) {
		super(a, b);

		setName("Molecule");
		atomFlags   = new int[a];
		partialCharges = new double[a];
		infos = new Object[a][MAX_INFOS];
//		mm2AtomTypes = new int[a];
//		mmffAtomTypes = new int[a];
//		classIds = new int[a];
	}

	@Override
	public int getOccupiedValence(int atom) {
		return super.getOccupiedValence(atom) + getAttachedHydrogenCount(atom);
	}

	public int getAttachedHydrogenCount(int atom) {
		return (atomFlags[atom] & ATTACHED_HYDROGEN_COUNT) >> ATTACHED_HYDROGEN_COUNT_SHIFT;
	}

	/**
	 * Defines a number of attached hydrogens without known 3D-coordinates
	 * @param atom
	 * @param value 0 to 3
	 */
	public void setAttachedHydrogenCount(int atom, int value) {
		atomFlags[atom] |= (value << ATTACHED_HYDROGEN_COUNT_SHIFT);
	}

		@Override
	public String toString() {
		return (getName() !=null ? getName() : "");
	}

	public void clear() {
		super.deleteMolecule();
	}

	public final void setAtomFlag(int atm, int flag, boolean value) {
		nMovables = -1;
		if(value) atomFlags[atm] |= flag;
		else atomFlags[atm] &= ~flag;
	}

	public final boolean isAtomFlag(int atm, int flag) {
		return (atomFlags[atm] & flag) != 0;
	}
	public final int getAtomFlags(int atm) {
		return atomFlags[atm];
	}
	public final void setAtomFlags(int atm, int flags) {
		nMovables = -1;
		atomFlags[atm] = flags;
	}
	public final void setAllAtomFlag(int flag, boolean value) {
		nMovables = -1;
		for (int i = 0; i < getAllAtoms(); i++) setAtomFlag(i, flag, value);
	}

	/**
	 * This method will append a Molecule3D to the end.
	 * @param m2
	 * @return the index dividing the 2 molecules
	 */
	public int fusion(Molecule3D m2) {
		if(m2==this) throw new IllegalArgumentException("Cannot fusion a molecule with itself");
		int index = getAllAtoms();
		int[] oldToNew = new int[m2.getAllAtoms()];
		for(int i=0; i<m2.getAllAtoms(); i++) {
			oldToNew[i] = addAtom(m2, i);
		}
		for(int i=0; i<m2.getAllBonds(); i++) {
			addBond(oldToNew[m2.getBondAtom(0, i)], oldToNew[m2.getBondAtom(1, i)], m2.getBondOrder(i));
		}
		return index;
	}

/*	public int fusionKeepCoordinatesObjects(Molecule3D m2) {
		int i = getAllAtoms();
		int res = fusion(m2);
		System.arraycopy(m2.getCoordinates(), 0, getCoordinates(), i, m2.getAllAtoms());
		return res;

	}*/

	public void setBondOrder(int bond, int order) {
		if (order == 1)
			super.setBondType(bond, Molecule.cBondTypeSingle);
		else if (order == 2)
			super.setBondType(bond, Molecule.cBondTypeDouble);
		else if (order == 3)
			super.setBondType(bond, Molecule.cBondTypeTriple);
		else
			System.out.println("ERROR: Unusual bond order:"+order);
	}

	public final void setAtomDescription(int atm, String s) {
		infos[atm][INFO_DESCRIPTION] = s;
	}

	public final String getAtomDescription(int atm) {
		return (String) infos[atm][INFO_DESCRIPTION];
	}

	public final void setPPP(int atm, int[] a) {
		infos[atm][INFO_PPP] = a;
	}

	public final int[] getPPP(int atm) {
		return (int[]) infos[atm][INFO_PPP];
	}
	public final void setAtomSequence(int atm, int a) {
		infos[atm][INFO_ATOMSEQUENCE] = a;
	}
	
	public void setResSequence(int atm, int a) {
		infos[atm][INFO_RESSEQUENCE] = a;
	}

	public final int getResSequence(int atm) {
		return infos[atm][INFO_RESSEQUENCE]==null?-1 : (Integer) infos[atm][INFO_RESSEQUENCE];
	}
	
	public final int getAtomSequence(int atm) {
		return infos[atm][INFO_ATOMSEQUENCE]==null?-1 : (Integer) infos[atm][INFO_ATOMSEQUENCE];
	}

	public final void setAtomChainId(int atm, String a) {
		infos[atm][INFO_CHAINID] = a;
	}

	public final String getAtomChainId(int atm) {
		return (String) infos[atm][INFO_CHAINID];
	}

	public final void setAtomName(int atm, String a) {
		infos[atm][INFO_ATOMNAME] = a;
	}
	

	public final String getAtomName(int atm) {
		return (String) infos[atm][INFO_ATOMNAME];
	}

	public final void setAtomAmino(int atm, String a) {
		infos[atm][INFO_AMINO] = a;
	}

	public final String getAtomAmino(int atm) {
		return (String) infos[atm][INFO_AMINO];
	}
	

	
	public final double getAtomBfactor(int atm) {
		return (double) infos[atm][INFO_BFACTOR];
	}
	
	public final void setAtomBfactor(int atm, double bfactor) {
		infos[atm][INFO_BFACTOR] = bfactor;
	}

	public final int getBond(int a1, int a2) {
		for(int connBond=0; connBond<getAllConnAtomsPlusMetalBonds(a1); connBond++) {
			if(getConnAtom(a1, connBond)==a2) return getConnBond(a1, connBond);
		}
		return -1;

	}


	////////////////////////////// UTILITIES ////////////////////////////////////////


	public final Coordinates[] getCoordinates() {

		Coordinates [] coords = new Coordinates[getAtoms()];

		for (int i = 0; i < getAtoms(); i++) {
			coords[i]=getCoordinates(i);
		}

		return coords;
	}


	public final void deleteAtoms(List<Integer> atomsToBeDeleted) {
		Collections.sort(atomsToBeDeleted);
		for (int i = atomsToBeDeleted.size()-1; i>=0; i--) {
			deleteAtom(atomsToBeDeleted.get(i));
		}
	}

	public String getName() {
		String name = super.getName();
		return name == null ? "" : name;
	}

	public String getShortName() {
		String name = getName();
		if(name.indexOf(' ')>0) name = name.substring(0, name.indexOf(' '));
		if(name.length()>12) name = name.substring(0, 12);
		return name;
	}

	public Map<String, Object> getAuxiliaryInfos() {
		return auxiliaryInfos;
	}
	public Object getAuxiliaryInfo(String name) {
		return auxiliaryInfos.get(name);
	}
	public void setAuxiliaryInfo(String name, Object value) {
		if(value==null) {
			System.err.println("Attempt to set "+name+" to null");
			auxiliaryInfos.remove(name);
		} else {
			auxiliaryInfos.put(name, value);
		}
	}

	@Override
	public boolean equals(Object obj) {
		return obj==this;
	}

	public final void setInteractionAtomType(int atm, int type) {
		infos[atm][INFO_INTERACTION_ATOM_TYPE] = type;
	}
	
	public final int getInteractionAtomType(int atm) {
		return (int) infos[atm][INFO_INTERACTION_ATOM_TYPE];
	}

	@Override
	public void setMaxAtoms(int v) {
		super.setMaxAtoms(v);
		int u = atomFlags.length;
		atomFlags = Arrays.copyOf(atomFlags, v);
		partialCharges = Arrays.copyOf(partialCharges, v);
		infos = Arrays.copyOf(infos, v);
		for (int i=u; i<v; i++)
			infos[i] = new Object[MAX_INFOS];
	}

	@Override
	public int copyAtom(Molecule destMol, int sourceAtom, int esrGroupOffsetAND, int esrGroupOffsetOR) {
		int destAtom = super.copyAtom(destMol, sourceAtom, esrGroupOffsetAND, esrGroupOffsetOR);
		if (destMol instanceof Molecule3D) {

			//
			// If copy atom is called via the copy constructor the following arrays are not initialized.
			// They are filled later in the constructor.
			//
			if(((Molecule3D) destMol).atomFlags!=null) {
				((Molecule3D) destMol).atomFlags[destAtom] = atomFlags[sourceAtom];
			}

			if(((Molecule3D) destMol).partialCharges!=null){
				((Molecule3D) destMol).partialCharges[destAtom] = partialCharges[sourceAtom];
			}

			if(((Molecule3D) destMol).infos!=null) {
				((Molecule3D) destMol).infos[destAtom] = new Object[MAX_INFOS];//;= infos[sourceAtom];// = new Object[MAX_INFOS];
				for (int i = 0; i < infos[sourceAtom].length; i++) {
					((Molecule3D) destMol).infos[destAtom][i] = clone(infos[sourceAtom][i]);
				}
			}
		}
		return destAtom;
	}

	private Object clone(Object o) {	// since Object.clone() has protected access
		Object newObject = null;
		if(o!=null) {
		
			//System.out.println(o.getClass());
/*		if (o instanceof Cloneable) {	causes thousands of warning during obfuscation
			try {
				return o.getClass().getMethod("clone").invoke(o);
			} catch (Exception e) {}
		}*/

			if (o instanceof String)
				newObject = new String((String)o);
			else if(o instanceof Integer)
				newObject = new Integer((Integer)o);
			else if(o instanceof Double)
				newObject = new Double((Double)o);
			else {
				System.out.println("ERROR: unexpected Object type. Add support for new type: "+o);
			}
		}
		return newObject;
	}

	public void compact() {
		setMaxAtoms(getAllAtoms());
		setMaxBonds(getAllBonds());
	}

	@Override
	public void swapAtoms(int atom1, int atom2) {
		super.swapAtoms(atom1, atom2);
		int tempi = atomFlags[atom1];
		atomFlags[atom1] = atomFlags[atom2];
		atomFlags[atom2] = tempi;
		double tempd = partialCharges[atom1];
		partialCharges[atom1] = partialCharges[atom2];
		partialCharges[atom2] = tempd;
		Object[] tempo = infos[atom1];
		infos[atom1] = infos[atom2];
		infos[atom2] = tempo;
	}

	/**
	 * Add an atom with the given atomicNo
	 */
	public int addAtom(int atomicNo) {
		if(mAllAtoms > mMaxAtoms){
			setMaxAtoms(mMaxAtoms*2);
		}

		int atom = super.addAtom(atomicNo);

		atomFlags[atom] = 0;
		partialCharges[atom] = 0;
		infos[atom] = new Object[MAX_INFOS];

		/* TODO
		mm2AtomTypes[atom] = -1;
		mmffAtomTypes[atom] = -1;
		classIds[atom] = -1;
		*/

		return atom;
	}

	/**
	 * Add an atom by copying its properties from the given Molecule3D
	 * This has to be overriden by subclasses
	 * @param m
	 * @param i
	 * @return
	 */
	public int addAtom(Molecule3D m, int i) {
		int a = addAtom(m.getAtomicNo(i));
		atomFlags[a] = m.getAtomFlags(i);
		partialCharges[a] = m.getPartialCharge(i);
		infos[a] = m.infos[i].clone();
		setAtomX(i, getAtomX(i));
		setAtomY(i, getAtomY(i));
		setAtomZ(i, getAtomZ(i));
//		mm2AtomTypes[a] = m.getMM2AtomType(i);
//		mmffAtomTypes[a] = m.getMMFFAtomType(i);
//		classIds[a] = m.getAtomInteractionClass(i);
//		nMovables = -1;
		return a;
	}
	
	/**
	 * High level function for constructing a molecule.
	 * After the deletion the original order of atom and bond indexes is retained.
	 * @param atom
	 */
	@Override
	public void deleteAtom(int atom) {
		super.deleteAtom(atom);
		compressMolTable();
		mValidHelperArrays = cHelperNone;
		}
	
	@Override 
	
	protected int[] compressMolTable() {
		// neutralize charges, if after deletion of a bi-polar bond one atom remains
		for (int bnd=0; bnd<mAllBonds; bnd++) {
			if (mBondType[bnd] == cBondTypeDeleted) {
				int atom1 = mBondAtom[0][bnd];
				int atom2 = mBondAtom[1][bnd];
				if (mAtomicNo[atom1] == -1
				  ^ mAtomicNo[atom2] == -1) {
					if (mAtomCharge[atom1] != 0
					 && mAtomCharge[atom2] != 0) {
						if (mAtomCharge[atom1] < 0
						  ^ mAtomCharge[atom2] < 0) {
							if (mAtomCharge[atom1] < 0) {
								mAtomCharge[atom1]++;
								mAtomCharge[atom2]--;
								}
							else {
								mAtomCharge[atom1]--;
								mAtomCharge[atom2]++;
								}
							}
						}
					}
				}
			}
		
		int newAtmNo[] = new int[mAllAtoms];
		int atomDest = 0;
		for (int atom=0; atom<mAllAtoms; atom++) {
			if (mAtomicNo[atom] == -1) {
				newAtmNo[atom] = -1;
				continue;
				}
			if (atomDest < atom) {
				mAtomicNo[atomDest] = mAtomicNo[atom];
				mAtomCharge[atomDest] = mAtomCharge[atom];
				mAtomMass[atomDest] = mAtomMass[atom];
				mAtomFlags[atomDest] = mAtomFlags[atom];
				mAtomQueryFeatures[atomDest] = mAtomQueryFeatures[atom];
				mAtomMapNo[atomDest] = mAtomMapNo[atom];
				mCoordinates[atomDest].set(mCoordinates[atom]);
				atomFlags[atomDest] = atomFlags[atom];
				if(infos[atomDest]==null) {
					infos[atomDest] = new Object[MAX_INFOS];
				}
				for (int i = 0; i < infos[atom].length; i++) {
					infos[atomDest][i] = clone(infos[atom][i]);
				}
				partialCharges[atomDest] = partialCharges[atom];
				if (mAtomList != null)
					mAtomList[atomDest] = mAtomList[atom];
				if (mAtomCustomLabel != null)
					mAtomCustomLabel[atomDest] = mAtomCustomLabel[atom];
				}
			newAtmNo[atom] = atomDest;
			atomDest++;
			}
		mAllAtoms = atomDest;

		int bondDest = 0;
		for (int bnd=0; bnd<mAllBonds; bnd++) {
			if (mBondType[bnd] == cBondTypeDeleted) continue;
			mBondType[bondDest] = mBondType[bnd];
			mBondFlags[bondDest] = mBondFlags[bnd];
			mBondQueryFeatures[bondDest] = mBondQueryFeatures[bnd];
			mBondAtom[0][bondDest] = newAtmNo[mBondAtom[0][bnd]];
			mBondAtom[1][bondDest] = newAtmNo[mBondAtom[1][bnd]];
			bondDest++;
			}
		mAllBonds = bondDest;

		return newAtmNo;
		}

	///////////////////// UTILITY FUNCTIONS ////////////////////////////////////
	/**
	 * Reorganizes atom indexes, so that moveable atoms are first
	 * @return the number of moveable atoms
	 */
	public boolean reorderAtoms() {
		boolean changed = false;
		int N = getAllAtoms();
		int i = 0; 	 //index of the first moveable atom
		int j = N-1; //index of the last non-moveable atom

		if(N==0) {
			nMovables = 0;
			return false;
		}

		while(i<j) {
			//Move i to the first non-moveable atom
			while(i<j && !isAtomFlag(i, RIGID)) i++;

			//Move j to the last moveable atom
			while(i<j && isAtomFlag(j, RIGID)) j--;
			if(isAtomFlag(i, RIGID) && !isAtomFlag(j, RIGID)) {
				swapAtoms(i, j);
				changed = true;
			}
		}
		nMovables = isAtomFlag(i, RIGID)? i: i+1;

		return changed;
	}


	public void setCoordinates(int atom, Coordinates c){

		setAtomX(atom, c.x);
		setAtomY(atom, c.y);
		setAtomZ(atom, c.z);
	}

	/*
	@Override
	public void ensureHelperArrays(int level) {
		boolean doRings = (level & Molecule.cHelperBitRings) != 0 && (mValidHelperArrays & Molecule.cHelperBitRings) == 0;
		super.ensureHelperArrays(level);
		if (doRings) {
			RingCollection ringSet = getRingSet();
			mAtomToRings = new ArrayList[getAtoms()];
			for (int r=0; r<ringSet.getSize(); r++) {
				int[] ringAtom = ringSet.getRingAtoms(r);
				for (int atom:ringAtom) {
					if (mAtomToRings[atom] == null)
						mAtomToRings[atom] = new ArrayList();
					mAtomToRings[atom].add(r);
				}
			}
		}
	}

	public List<Integer>[] getAtomToRings() {
		return mAtomToRings;
	}
	*/

	/**
	 * @return the number of movable atoms (after reorderatoms has been called)
	 */
	public int getNMovables() {
		if(nMovables<0 || (getAllAtoms()>0 && !isAtomFlag(0, LIGAND))) reorderAtoms();
		return nMovables;
	}
	public boolean isMoleculeInOrder() {
		return reorderAtoms()==false;
	}

	public double getPartialCharge(int a) {
		return partialCharges[a];
	}
	public void setPartialCharge(int a, double v) {
		partialCharges[a] = v;
	}

	@Override
	public int compareTo(Molecule3D o) {
		return o==null? 1: getName().compareTo(o.getName());
	}

/*	@Override
	public int hashCode() {
		return name==null? 0 : name.hashCode();
	}*/

	/**
	 * Returns the Interaction Atom Type
	 * @param atm
	 * @param id
	 *
	public final void setAtomInteractionClass(int atm, int id) {
		classIds[atm] = id;
	}*/

	/**
	 * Set the MMFF atom type
	 * @param atm
	 * @param atomType
	 *
	public final void setMMFFAtomType(int atm, int atomType) {
		mmffAtomTypes[atm] = atomType;
	}*/

	/**
	 * Returns the MMFF Atom Type
	 * @param atm
	 * @return
	 *
	public final int getMMFFAtomType(int atm) {
		return mmffAtomTypes[atm];
	}*/

	/**
	 * Returns the MM2 Atom Type
	 * @param atm
	 * @return
	 *
	public final int getMM2AtomType(int atm) {
		return mm2AtomTypes[atm];
	}*/

	/**
	 * Set the MM2 atom type
	 * @param atm
	 * @param atomType
	 *
	public final void setMM2AtomType(int atm, int atomType) {
		mm2AtomTypes[atm] = atomType;
	}*/


/*	public final String getMM2AtomDescription(int atm) {
		return (String) infos[atm][INFO_MM2ATOMDESCRIPTION];
	}

	public final int getAtomInteractionClass(int atm) {
		return classIds[atm];
	}


	public void setMmffRingAtoms(Boolean[] ringAtoms) {
		this.mmffRingAtoms = ringAtoms;
	}
*/

	/**
	 * Determine if a ring is aromatic according to MMFF criteria. Only
	 * designed to work with rings of size 5 and 6. Returns the cached value.
	 *  @param r The ring index in the molecule.
	 *  @return True if the ring is aromatic, false otherwise.
	 *
	public boolean ringIsMMFFAromatic(int r) {
		return mmffRingAtoms[r] == Boolean.TRUE;
	}*/

	/**
	 * Returns true if the given ring has had its MMFF aromaticity flag set.
	 *  @param r The ring index in the molecule.
	 *  @return True if the ring has had its flag set, false otherwise.
	 *
	public boolean isSetRingMMFFAromaticity(int r) {
		return mmffRingAtoms[r] != null;
	}*/

}
