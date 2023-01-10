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
 * @author Thomas Sander
 */

package com.actelion.research.chem;

import com.actelion.research.util.Angle;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

/**
 * While the Molecule class covers all primary molecule information as atom and bond properties,
 * the atom connectivity and coordinates, its derived class ExtendedMolecule handles secondary,
 * i.e. calculated molecule information. Most important are the directly connected atoms and bonds
 * of every atom, information about rings and whether they are aromatic, and some atom properties
 * that depend on their near neighbors. This calculated information is cached in helper arrays and
 * stays valid as long as the molecule's primary information is not changed. High level methods,
 * e.g. getPath(), that require valid helper arrays take care of updating the cache themselves.
 * Low level methods, e.g. isAromaticAtom(), which typically are called very often, do not check
 * the validity of or update the helper arrays themselves for performance reasons. If you use low
 * level methods, then you need to make sure that the required helper array information is valid by
 * calling ensureHelperArrays().
 * Typically you will never instantiate an ExtendedMolecule, but rather use a StereoMolecule.
 */
public class ExtendedMolecule extends Molecule implements Serializable {
	static final long serialVersionUID = 0x2006CAFE;

	/**
	 * To interpret a stereo center as fisher projection, all non stereo
	 * bonds must be vertical and all stereo bonds must be horizontal.
	 * FISCHER_PROJECTION_LIMIT is the allowed tolerance (currently 5.0 degrees).
	 * In addition in large rings we don't assume Fisher projections,
	 * because coordinates, if just in a circle, may have subsequent almost vertical bonds.
	 */
	public static final float FISCHER_PROJECTION_LIMIT = (float)Math.PI / 36;
	public static final float FISCHER_PROJECTION_RING_LIMIT = 24;

	public static final float STEREO_ANGLE_LIMIT = (float)Math.PI / 36;   // 5 degrees

	public static final int cMaxConnAtoms = 16; // ExtendedMolecule is not restricted anymore
											   // However, this is a suggestion for editors and other classes
	transient private int mAtoms,mBonds;
	transient private RingCollection mRingSet;

	transient private int mPi[];
	transient private int mConnAtoms[];	// non-H neighbour counts
	transient private int mAllConnAtoms[];	// neighbour counts including explicit hydrogen
	transient private int mConnAtom[][];
	transient private int mConnBond[][];
	transient private int mConnBondOrder[][];

	public ExtendedMolecule() {
		}


	public ExtendedMolecule(int maxAtoms, int maxBonds) {
		super(maxAtoms, maxBonds);
		}

	public ExtendedMolecule(Molecule mol) {
		super(mol==null ? 256 : mol.getMaxAtoms(), mol==null ? 256 : mol.getMaxBonds());
		if (mol != null)
			mol.copyMolecule(this);
		}


	/**
	 * Clears destmol and then copies a part of this Molecule into destMol, being defined by a mask of atoms to be included.
	 * If not all atoms are copied, then destMol is set to be a substructure fragment.
	 * @param destMol receives the part of this Molecule
	 * @param includeAtom defines atoms to be copied; its size may be this.getAtoms() or this.getAllAtoms()
	 * @param recognizeDelocalizedBonds defines whether disconnected delocalized bonds will keep their
	 * single/double bond status or whether the query feature 'delocalized bond' will be set
	 * @param atomMap null or int[] not smaller than includeAtom.length; receives atom indices of dest molecule or -1 if not copied
	 */
	public void copyMoleculeByAtoms(ExtendedMolecule destMol, boolean[] includeAtom, boolean recognizeDelocalizedBonds, int[] atomMap) {
		if (recognizeDelocalizedBonds)
			ensureHelperArrays(cHelperRings);

		destMol.mAtomList = null;
		if (mIsFragment)
			destMol.setFragment(true);

		int atomCount = includeAtom.length;

		if (atomMap == null)
			atomMap = new int[atomCount];

		destMol.mAllAtoms = 0;
		for (int atom=0; atom<atomCount;atom++)
			atomMap[atom] = includeAtom[atom] ? copyAtom(destMol, atom, 0, 0) : -1;

		destMol.mAllBonds = 0;
		for (int bnd=0; bnd<mAllBonds;bnd++) {
			int atom1 = mBondAtom[0][bnd];
			int atom2 = mBondAtom[1][bnd];
			if (atom1<atomCount && atom2<atomCount) {
				if (includeAtom[atom1] && includeAtom[atom2])
					copyBond(destMol, bnd, 0, 0, atomMap, recognizeDelocalizedBonds);
				// if we keep only one atom and have a bi-polar bond, then we need to neutralize
				else if (mAtomCharge[atom1] != 0
					  && mAtomCharge[atom2] != 0
					  && (mAtomCharge[atom1] < 0
						^ mAtomCharge[atom2] < 0)) {
					if (includeAtom[atom1])
						destMol.mAtomCharge[atomMap[atom1]] += (mAtomCharge[atom1] < 0) ? 1 : -1;
					if (includeAtom[atom2])
						destMol.mAtomCharge[atomMap[atom2]] += (mAtomCharge[atom2] < 0) ? 1 : -1;
					}
				}
			}

		copyMoleculeProperties(destMol);
		destMol.mValidHelperArrays = cHelperNone;

		destMol.renumberESRGroups(cESRTypeAnd);
		destMol.renumberESRGroups(cESRTypeOr);

		if (destMol.mAllAtoms != atomCount)
			destMol.setFragment(true);

		if (recognizeDelocalizedBonds)
			new AromaticityResolver(destMol).locateDelocalizedDoubleBonds(null);
		}


	/**
	 * Clears destmol and then copies a part of this Molecule into destMol, being defined by a mask of bonds to be included.
	 * Bonds, whose atoms carry opposite charges are treated in the following manner: If only one of
	 * the two bond atoms is kept, then its absolute charge will be reduced by 1.
	 * @param destMol receives the part of this Molecule
	 * @param includeBond defines bonds to be copied
	 * @param recognizeDelocalizedBonds defines whether disconnected delocalized bonds will keep their
	 * single/double bond status or whether the query feature 'delocalized bond' will be set
	 * @param atomMap null or int[] not smaller than this.getAllAtoms()
	 * @return atom map from this to destMol with not copied atom's index being -1
	 */
	public int[] copyMoleculeByBonds(ExtendedMolecule destMol, boolean[] includeBond, boolean recognizeDelocalizedBonds, int[] atomMap) {
		if (recognizeDelocalizedBonds)
			ensureHelperArrays(cHelperRings);

		destMol.mAtomList = null;
		if (mIsFragment)
			destMol.setFragment(true);

		if (atomMap == null)
			atomMap = new int[mAllAtoms];

		destMol.mAllAtoms = 0;
		for (int atom=0; atom<mAllAtoms;atom++) {
			atomMap[atom] = -1;
			for (int i = 0; i< mConnAtoms[atom]; i++) {
				if (includeBond[mConnBond[atom][i]]) {
					atomMap[atom] = copyAtom(destMol, atom, 0, 0);
					break;
					}
				}
			}

		destMol.mAllBonds = 0;
		for (int bnd=0; bnd<mAllBonds;bnd++)
			if (includeBond[bnd]) {
				copyBond(destMol, bnd, 0, 0, atomMap, recognizeDelocalizedBonds);
				}
			// if we keep one atom of a removed bi-polar bonds, then we need to neutralize
			else {
				int atom1 = mBondAtom[0][bnd];
				int atom2 = mBondAtom[1][bnd];
				if (atomMap[atom1] == -1
				  ^ atomMap[atom2] == -1) {
					if (mAtomCharge[atom1] != 0
		  			 && mAtomCharge[atom2] != 0
					 && ((mAtomCharge[atom1] < 0)
					   ^ (mAtomCharge[atom2] < 0))) {
						if (atomMap[atom1] != -1)
							destMol.mAtomCharge[atomMap[atom1]] += (mAtomCharge[atom1] < 0) ? 1 : -1;
						if (atomMap[atom2] != -1)
							destMol.mAtomCharge[atomMap[atom2]] += (mAtomCharge[atom2] < 0) ? 1 : -1;
						}
					}
				}

		copyMoleculeProperties(destMol);
		destMol.mValidHelperArrays = cHelperNone;

		destMol.renumberESRGroups(cESRTypeAnd);
		destMol.renumberESRGroups(cESRTypeOr);

		if (destMol.mAllAtoms != mAllAtoms)
			destMol.setFragment(true);

		if (recognizeDelocalizedBonds)
			new AromaticityResolver(destMol).locateDelocalizedDoubleBonds(null);
 
		return atomMap;
		}


	/**
	 * The neighbours (connected atoms) of any atom are sorted by their relevance:<br>
	 * 1. non-hydrogen atoms (bond order 1 and above) and unusual hydrogen atoms (non-natural abundance isotopes, custom labelled hydrogen, etc.)<br>
	 * 2. plain-hydrogen atoms (natural abundance, bond order 1)<br>
	 * 3. loosely connected atoms (bond order 0, i.e. metal ligand bond)<br>
	 * Only valid after calling ensureHelperArrays(cHelperNeighbours or higher);
	 * @param atom
	 * @return count of category 1 & 2 neighbour atoms (excludes neighbours connected with zero bond order)
	 */
	public int getAllConnAtoms(int atom) {
		return mAllConnAtoms[atom];
		}


	/**
	 * @param atom
	 * @return the number of connected plain explicit and implicit hydrogen atoms (does not include 2H,3H, custom labelled H)
	 */
	public int getPlainHydrogens(int atom) {
		return getExplicitHydrogens(atom) + getImplicitHydrogens(atom);
		}


	/**
	 * @param atom
	 * @return the number of connected explicit and implicit hydrogen atoms (includes 2H,3H, custom labelled H)
	 */
	public int getAllHydrogens(int atom) {
		return mAllConnAtoms[atom] - getNonHydrogenNeighbourCount(atom) + getImplicitHydrogens(atom);
		}


	/**
	 * A validated molecule (after helper array creation) contains a sorted list of all atoms
	 * with the plain (negligible) hydrogen atoms at the end of the list. negligible hydrogen atoms
	 * a those that can be considered implicit, because they have no attached relevant information.
	 * Hydrogen atoms that cannot be neglected are special isotopes (mass != 0), if they carry a
	 * custom label, if they are connected to another atom with bond order different from 1, or
	 * if they are connected to another negligible hydrogen.<br>
	 * Only valid after calling ensureHelperArrays(cHelperNeighbours or higher);
	 * @return the number of relevant atoms not including negligible hydrogen atoms
	 */
	public int getAtoms() {
		return mAtoms;
		}


	/**
	 * @param atom
	 * @return count of neighbour atoms connected by a 0-order metal ligand bond
	 */
	public int getMetalBondedConnAtoms(int atom) {
		return mConnAtom[atom].length - mAllConnAtoms[atom];
		}


	/**
	 * This is different from the Hendrickson pi-value, which considers pi-bonds to carbons only.
	 * @param atom
	 * @return the number pi electrons at atom (the central atom of acetone would have 1)
	 */
	public int getAtomPi(int atom) {
		return mPi[atom];
	}


	/**
	 * Requires helper arrays state cHelperNeighbours.
	 * @param atom
	 * @return number of electronegative neighbour atoms
	 */
	public int getAtomElectronegativeNeighbours(int atom) {
		int e = 0;
		for (int i=0; i<mConnAtoms[atom]; i++)
			if (isElectronegative(mConnAtom[atom][i]) && !isBondBridge(mConnBond[atom][i]))
				e++;

		return e;
		}


	/**
	 * @param atom
	 * @return Hendrickson sigma-value, which is the number attached carbon atoms
	 *
	public int getAtomSigma(int atom) {
		int sigma = 0;
		for (int i=0; i<mConnAtoms[atom]; i++)
			if (mAtomicNo[mConnAtom[atom][i]] == 6)
				sigma++;
		return sigma;
		}*/


	/**
	 * Hendrickson Z-value, which is the sum of all bond orders to any attached hetero atoms.
	 * Delocalized bonds orders are considered 1.5. Since the z-value is an integer,
	 * one delocalized hetero atom results in z=1, two delocalized hetero atoms give z=3.
	 * Requires helper arrays state cHelperRings.
	 * @param atom
	 * @return Hendrickson Z-value
	 */
	public int getAtomZValue(int atom) {
		int z = 0;
		int arom = 0;
		for (int i=0; i<mConnAtoms[atom]; i++) {
			if (isElectronegative(mConnAtom[atom][i])) {
				if (isDelocalizedBond(mConnBond[atom][i])) {
					z++;
					arom++;
					}
				else {
					z += mConnBondOrder[atom][i];
					}
				}
			}
		return z + arom/2;
		}

	
	/**
	 * @param atom
	 * @return 0 or the size of the smallest ring that atom is a member of
	 */
	public int getAtomRingSize(int atom) {
		return (mRingSet != null && atom<mAtoms) ?
				mRingSet.getAtomRingSize(atom) : 0;
		}


	/**
	 * @param bond
	 * @return 0 or the size of the smallest ring that bond is a member of
	 */
	public int getBondRingSize(int bond) {
		return (mRingSet != null && bond<mBonds) ?
				mRingSet.getBondRingSize(bond) : 0;
		}


	/**
	 * The bond list is preprocessed such that all bonds leading to a plain hydrogen atom
	 * (natural abundance, no custom labels) are at the end of the list.
	 * Only valid after calling ensureHelperArrays(cHelperNeighbours or higher);
	 * @return count of bonds not including those connecting plain-H atoms
	 */
	public int getBonds() {
		return mBonds;
		}


	/**
	 * @return -1 or the bond that connects atom1 with atom2
	 */
	public int getBond(int atom1, int atom2) {
		for (int i=0; i<getAllConnAtomsPlusMetalBonds(atom1); i++)
			if (mConnAtom[atom1][i] == atom2)
				return mConnBond[atom1][i];
		return -1;
		}


	/**
	 * @return a copy of this with all arrays sized to just cover all existing atoms and bonds
	 */
	public ExtendedMolecule getCompactCopy() {
		ExtendedMolecule theCopy = new ExtendedMolecule(mAllAtoms, mAllBonds);
		copyMolecule(theCopy);
		return theCopy;
		}


	/**
	 * The neighbours (connected atoms) of any atom are sorted by their relevance:<br>
	 * 1. non-hydrogen atoms (bond order 1 and above) and unusual hydrogen atoms (non natural abundance isotops, custom labelled hydrogen, etc.)<br>
	 * 2. plain-hydrogen atoms (natural abundance, bond order 1)<br>
	 * 3. loosely connected atoms (bond order 0, i.e. metall ligand bond)<br>
	 * Only valid after calling ensureHelperArrays(cHelperNeighbours or higher);
	 * @param atom
	 * @param i index into sorted neighbour list
	 * @return the i-th neighbor atom of atom
	 */
	public int getConnAtom(int atom, int i) {
		return mConnAtom[atom][i];
		}


	/**
	 * The neighbours (connected atoms) of any atom are sorted by their relevance:<br>
	 * 1. non-hydrogen atoms (bond order 1 and above) and unusual hydrogen atoms (non natural abundance isotops, custom labelled hydrogen, etc.)<br>
	 * 2. plain-hydrogen atoms (natural abundance, bond order 1)<br>
	 * 3. loosely connected atoms (bond order 0, i.e. metall ligand bond)<br>
	 * Only valid after calling ensureHelperArrays(cHelperNeighbours or higher);
	 * @param atom
	 * @return count of category 1 neighbour atoms (excludes plain H and bond zero orders)
	 */
	public int getConnAtoms(int atom) {
		return mConnAtoms[atom];
		}


	/**
	 * The neighbours (connected atoms) of any atom are sorted by their relevance:<br>
	 * 1. non-hydrogen atoms (bond order 1 and above) and unusual hydrogen atoms (non natural abundance isotops, custom labelled hydrogen, etc.)<br>
	 * 2. plain-hydrogen atoms (natural abundance, bond order 1)<br>
	 * 3. loosely connected atoms (bond order 0, i.e. metall ligand bond)<br>
	 * Only valid after calling ensureHelperArrays(cHelperNeighbours or higher);
	 * @param atom
	 * @return count of category 1 & 2 & 3 neighbour atoms
	 */
	public int getAllConnAtomsPlusMetalBonds(int atom) {
		return mConnAtom[atom].length;
		}


	/**
	 * The neighbours (connected atoms) of any atom are sorted by their relevance:<br>
	 * 1. non-hydrogen atoms (bond order 1 and above) and unusual hydrogen atoms (non natural abundance isotops, custom labelled hydrogen, etc.)<br>
	 * 2. plain-hydrogen atoms (natural abundance, bond order 1)<br>
	 * 3. loosely connected atoms (bond order 0, i.e. metall ligand bond)<br>
	 * Only valid after calling ensureHelperArrays(cHelperNeighbours or higher);
	 * @param atom
	 * @param i index into sorted neighbour list
	 * @return index of bond connecting atom with its i-th neighbor
	 */
	public int getConnBond(int atom, int i) {
		return mConnBond[atom][i];
		}


	/**
	 * The neighbours (connected atoms) of any atom are sorted by their relevance:<br>
	 * 1. non-hydrogen atoms (bond order 1 and above) and unusual hydrogen atoms (non natural abundance isotops, custom labelled hydrogen, etc.)<br>
	 * 2. plain-hydrogen atoms (natural abundance, bond order 1)<br>
	 * 3. loosely connected atoms (bond order 0, i.e. metall ligand bond)<br>
	 * Only valid after calling ensureHelperArrays(cHelperNeighbours or higher);
	 * Orders of delocalized bonds, i.e. bonds in an aromatic 6-membered ring, are returned as 1.
	 * @param atom
	 * @param i index into sorted neighbour list
	 * @return order of bond connecting atom with its i-th neighbor
	 */
	public int getConnBondOrder(int atom, int i) {
		return mConnBondOrder[atom][i];
		}


	/**
	 * This method returns the non-hydrogen neighbour count of atom.
	 * It excludes any hydrogen atoms in contrast to getConnAtoms(), which only
	 * excludes plain hydrogen (not deuterium, tritium, custom labelled hydrogen, etc.).
	 * Don't use this method's return value for loops with getConnAtom(),
	 * getConnBond(), or getConnBondOrder().
	 * @param atom
	 * @return the number of non-hydrogen neighbor atoms
	 */
	public int getNonHydrogenNeighbourCount(int atom) {
		int count = mConnAtoms[atom];
		for (int i = 0; i< mConnAtoms[atom]; i++)
			if (mAtomicNo[mConnAtom[atom][i]] == 1)
				count--;
		return count;
		}


	/**
	 * This method returns the count of atom neighbours which are marked as being an exclude group.
	 * @param atom
	 * @return the number of non-hydrogen neighbor atoms
	 */
	public int getExcludedNeighbourCount(int atom) {
		int count = 0;
		for (int i=0; i<mConnAtoms[atom]; i++)
			if ((mAtomQueryFeatures[i] & Molecule.cAtomQFExcludeGroup) != 0)
				count++;
		return count;
		}


	/**
	 * Calculates and returns the mean bond length of all bonds including or not
	 * including hydrogen bonds.
	 * If there are no bonds, then the average distance between unconnected atoms is
	 * returned. If we have less than 2 atoms, cDefaultAverageBondLength is returned.
	 * @param nonHydrogenBondsOnly
	 * @return
	 */
	public double getAverageBondLength(boolean nonHydrogenBondsOnly) {
		if (nonHydrogenBondsOnly) {
			ensureHelperArrays(cHelperNeighbours);
			return getAverageBondLength(mAtoms, mBonds);
			}
		else {
			return getAverageBondLength(mAllAtoms, mAllBonds);
			}
		}


	/**
	 * Creates an array that maps connAtoms/connBonds sorted by atom indices.
	 * getConnAtom(atom, getSortedConnMap(atom, 0)) retrieves that neighbour
	 * of atom with the lowest atom index, i.e. that is the first in the atom table.
	 * @return neighbour index map
	 */
	private int[] getSortedConnMap(int atom) {
		int connAtoms = mAllConnAtoms[atom];
		int[] indexMap = new int[connAtoms];
		for (int i=0; i<connAtoms; i++)
			indexMap[i] = (mConnAtom[atom][i] << 16) + i;
		java.util.Arrays.sort(indexMap);
		for (int i=0; i<connAtoms; i++)
			indexMap[i] &= 0x0000FFFF;
		return indexMap;
		}


	/**
	 * The sum of bond orders of explicitly connected neighbour atoms including explicit hydrogen.
	 * In case of a fragment the occupied valence does not include bonds to atoms of which the cAtomQFExcludeGroup flag is set.
	 * Atom charge and radical states are not considered.
	 * @param atom
	 * @return explicitly used valence
	 */
	public int getOccupiedValence(int atom) {
		ensureHelperArrays(cHelperNeighbours);

		boolean piElectronsFound = false;
		boolean delocalizedBondFound = false;
		int valence = 0;
		for (int i=0; i<mAllConnAtoms[atom]; i++) {
			if (!mIsFragment || (mAtomQueryFeatures[mConnAtom[atom][i]] & cAtomQFExcludeGroup) == 0) {
				int order = mConnBondOrder[atom][i];
				valence += order;
				if (order > 1)
					piElectronsFound = true;

				int bond = mConnBond[atom][i];
				if (mBondType[bond] == cBondTypeDelocalized)
					delocalizedBondFound = true;
				}
			}

		if (delocalizedBondFound && !piElectronsFound)
			valence++;

		return valence;
		}


	/**
	 * The free valence is the number of potential additional single bonded
	 * neighbours to reach the atom's maximum valence. Atomic numbers that have
	 * multiple possible valences, the highest value is taken.
	 * Atom charges are considered. Implicit hydrogens are not considered.
	 * Thus, the oxygen in a R-O(-) has a free valence of 0, the nitrogen in R3N(+)
	 * has a free valence of 1. Chlorine in Cl(-) has a free valence of 6. If you need
	 * the free valence taking the lowest possible valence into account, use
	 * getLowestFreeValence(), which would return 0 for Cl(-).
	 * @param atom
	 * @return
	 */
	public int getFreeValence(int atom) {
		return getMaxValence(atom) - getOccupiedValence(atom);
		}


	/**
	 * The lowest free valence is the number of potential additional single bonded
	 * neighbours to reach the atom's lowest valence above or equal its current
	 * occupied valence. Atom charges are considered. Implicit hydrogens are not considered.
	 * Thus, the phosphor atoms in PF2 and PF4 both have a lowest free valence of 1.
	 * The oxygen in R-O(-) has a lowest free valence of 0, the nitrogen in R3N(+)
	 * has a free valence of 1. If you need the maximum possible free valence,
	 * use getFreeValence(), which would give 6 for Cl(-) and HCl.
	 * @param atom
	 * @return
	 */
	public int getLowestFreeValence(int atom) {
		int occupiedValence = getOccupiedValence(atom);
		int correction = getElectronValenceCorrection(atom, occupiedValence);

		int valence = getAtomAbnormalValence(atom);
		if (valence == -1) {
			byte[] valenceList = getAllowedValences(mAtomicNo[atom]);
			int i= 0;
			while ((occupiedValence > valenceList[i] + correction) && (i<valenceList.length-1))
				i++;
			valence = valenceList[i];
			}

		return valence + correction - occupiedValence;
		}


	/**
	 * If the explicitly attached neighbors cause an atom valence to exceed
	 * the lowest allowed valence for this atomic no, then this method returns
	 * the next higher allowed valence, e.g. O=P(-H)-OMe :<br>
	 * standard P valence is 3, used valence is 4, implicit abnormal valence is 5.
	 * The molecule is interpreted as O=PH2-OMe. Requires cHelperNeighbours!
	 * @param atom
	 * @param neglectExplicitHydrogen
	 * @return abnormal valence or -1 if valence doesn't exceed standard valence
	 */
	public int getImplicitHigherValence(int atom, boolean neglectExplicitHydrogen) {
		int occupiedValence = getOccupiedValence(atom);
		occupiedValence -= getElectronValenceCorrection(atom, occupiedValence);
		if (neglectExplicitHydrogen)
			occupiedValence -= mAllConnAtoms[atom] - mConnAtoms[atom];

		byte[] valences = getAllowedValences(mAtomicNo[atom]);
		if (occupiedValence <= valences[0])
			return -1;

		// stepwise try higher allowed valences
		for (int i=1; i<valences.length; i++)
			if (valences[i] >= occupiedValence)
				return valences[i];

		// if we have a compatible higher allowed valence, then use that, otherwise the occupied valence
		return occupiedValence;
		}

	/**
	 * Calculates for every non-H atom the mean value of all shortest routes (bonds in between)
	 * to any other atom of the same fragment.
	 * @return 
	 */
	public float[] getAverageTopologicalAtomDistance() {
		ensureHelperArrays(cHelperNeighbours);

		float[] meanDistance = new float[mAtoms];
		int[] graphAtom = new int[mAtoms];
		for (int startAtom=0; startAtom<mAtoms; startAtom++) {
			graphAtom[0] = startAtom;
			int[] graphLevel = new int[mAtoms];
			graphLevel[startAtom] = 1;

			int current = 0;
			int highest = 0;
			while (current <= highest) {
				for (int i = 0; i< mConnAtoms[graphAtom[current]]; i++) {
					int candidate = mConnAtom[graphAtom[current]][i];
					if (graphLevel[candidate] == 0) {
						graphLevel[candidate] = graphLevel[graphAtom[current]] + 1;
						graphAtom[++highest] = candidate;
						meanDistance[startAtom] += (graphLevel[candidate] - 1);
						}
					}
				current++;
				}
			meanDistance[startAtom] /= highest;
			}

		return meanDistance;
		}

	/**
	 * Calculates the length of the shortest path between atoms atom1 and atom2
	 * @param atom1
	 * @param atom2
	 * @return path length (no of bonds); -1 if there is no path
	 */
	public int getPathLength(int atom1, int atom2) {
		if (atom1 == atom2)
			return 0;

		ensureHelperArrays(cHelperNeighbours);

		int[] graphLevel = new int[mAllAtoms];
		int graphAtom[] = new int[mAllAtoms];

		graphAtom[0] = atom1;
		graphLevel[atom1] = 1;
		int current = 0;
		int highest = 0;
		while (current <= highest) {
			for (int i=0; i<mAllConnAtoms[graphAtom[current]]; i++) {
				int candidate = mConnAtom[graphAtom[current]][i];
				if (candidate == atom2)
					return graphLevel[graphAtom[current]];
				if (graphLevel[candidate] == 0) {
					graphAtom[++highest] = candidate;
					graphLevel[candidate] = graphLevel[graphAtom[current]]+1;
					}
				}
			current++;
			}
		return -1;
		}


	/**
	 * Calculates the length of the shortest path between atoms atom1 and atom2,
	 * which is not larger than maxLength and avoids atoms indicated by neglectAtom.
	 * @param atom1
	 * @param atom2
	 * @param maxLength paths larger than maxLength won't be detected
	 * @param neglectAtom null or atoms flagged which are forbidden path members
	 * @return path length (no of bonds); -1 if there is no path
	 */
	public int getPathLength(int atom1, int atom2, int maxLength, boolean[] neglectAtom) {
		if (atom1 == atom2)
			return 0;

		ensureHelperArrays(cHelperNeighbours);

		int[] graphLevel = new int[mAllAtoms];
		int graphAtom[] = new int[mAllAtoms];

		graphAtom[0] = atom1;
		graphLevel[atom1] = 1;
		int current = 0;
		int highest = 0;
		while (current <= highest && graphLevel[graphAtom[current]] <= maxLength) {
			for (int i=0; i<mAllConnAtoms[graphAtom[current]]; i++) {
				int candidate = mConnAtom[graphAtom[current]][i];
				if (candidate == atom2)
					return graphLevel[graphAtom[current]];

				if (graphLevel[candidate] == 0
				 && (neglectAtom == null || neglectAtom.length <= candidate || !neglectAtom[candidate])) {
					graphAtom[++highest] = candidate;
					graphLevel[candidate] = graphLevel[graphAtom[current]]+1;
					}
				}
			current++;
			}
		return -1;
		}

	/**
	 * Locates and returns the shortest path between atoms atom1 and atom2.
	 * If neglectBond is not null, then flagged neglectBonds may not be part of the path.
	 * @param pathAtom array large enough to hold all path atoms, i.e. maxLength+1
	 * @param atom1 first atom of path; ends up in pathAtom[0]
	 * @param atom2 last atom of path; ends up in pathAtom[pathLength]
	 * @param maxLength paths larger than maxLength won't be detected
	 * @param neglectBond null or bitmask of forbidden bonds
	 * @return number of bonds of path; -1 if there is no path
	 */
	public int getPath(int[] pathAtom, int atom1, int atom2, int maxLength, boolean[] neglectBond) {
		return getPath(pathAtom, atom1, atom2, maxLength, null, neglectBond);
		}

	/**
	 * Locates and returns the shortest path between atoms atom1 and atom2.
	 * If neglectAtom is not null, then flagged neglectAtoms may not be part of the path, except for atom1 and atom2.
	 * If neglectBond is not null, then flagged neglectBonds may not be part of the path.
	 * @param pathAtom array large enough to hold all path atoms, i.e. maxLength+1
	 * @param atom1 first atom of path; ends up in pathAtom[0]
	 * @param atom2 last atom of path; ends up in pathAtom[pathLength]
	 * @param maxLength paths larger than maxLength won't be detected
	 * @param neglectAtom null or bitmask of forbidden atoms
	 * @param neglectBond null or bitmask of forbidden bonds
	 * @return number of bonds of path; -1 if there is no path
	 */
	public int getPath(int[] pathAtom, int atom1, int atom2, int maxLength, boolean[] neglectAtom, boolean[] neglectBond) {
		if (atom1 == atom2) {
			pathAtom[0] = atom1;
			return 0;
			}

		ensureHelperArrays(cHelperNeighbours);

		int[] graphLevel = new int[mAllAtoms];
		int graphAtom[] = new int[mAllAtoms];
		int parentAtom[] = new int[mAllAtoms];

		graphAtom[0] = atom1;
		graphLevel[atom1] = 1;

		int current = 0;
		int highest = 0;
		while (current <= highest && graphLevel[graphAtom[current]] <= maxLength) {
			int parent = graphAtom[current];
			for (int i=0; i<mAllConnAtoms[parent]; i++) {
				if (neglectBond == null
				 || neglectBond.length <= mConnBond[parent][i]
				 || !neglectBond[mConnBond[parent][i]]) {
					int candidate = mConnAtom[parent][i];
					if (candidate == atom2) {
						int index = graphLevel[parent];
						pathAtom[index] = candidate;
						pathAtom[--index] = parent;
						while (index > 0) {
							pathAtom[index-1] = parentAtom[pathAtom[index]];
							index--;
							}
						return graphLevel[parent];
						}
	
					if (graphLevel[candidate] == 0
					 && (neglectAtom == null || neglectAtom.length <= candidate || !neglectAtom[candidate])) {
						graphAtom[++highest] = candidate;
						graphLevel[candidate] = graphLevel[parent]+1;
						parentAtom[candidate] = parent;
						}
					}
				}
			current++;
			}
		return -1;
		}


	/**
	 * Finds bonds of a path that is defined by an atom sequence.
	 * @param pathAtom pathAtom[0]...[pathLength] -> list of atoms on path 
	 * @param pathBond int array not smaller than pathLength
	 * @param pathLength no of path bonds == no of path atoms - 1
	 */
	public void getPathBonds(int[] pathAtom, int[] pathBond, int pathLength) {
		ensureHelperArrays(cHelperNeighbours);
		for (int i=0; i<pathLength; i++) {
			for (int j=0; j<mAllConnAtoms[pathAtom[i]]; j++) {
				if (mConnAtom[pathAtom[i]][j] == pathAtom[i+1]) {
					pathBond[i] = mConnBond[pathAtom[i]][j];
					break;
					}
				}
			}
		}


	/**
	 * @param atom1
	 * @param atom2
	 * @return whether there is a path of bonds leading from atom1 to atom2
	 */
	public boolean shareSameFragment(int atom1, int atom2) {
		return (getPathLength(atom1, atom2) != -1);
		}


	/**
	 * This adds a fragment from sourceMol to this molecule by first copying rootAtom and then
	 * all connected atoms and bonds by traversing the graph breadth first.
	 * @param sourceMol molecule from which the fragment is copied to this
	 * @param rootAtom
	 * @param atomMap null or int[] not smaller than sourceMol.mAllAtoms; receives atom indices of this molecule
	 */
	public void addFragment(ExtendedMolecule sourceMol, int rootAtom, int[] atomMap) {
		sourceMol.ensureHelperArrays(cHelperNeighbours);

		if (atomMap == null)
			atomMap = new int[sourceMol.mAllAtoms];

		int esrGroupCountAND = renumberESRGroups(Molecule.cESRTypeAnd);
		int esrGroupCountOR = renumberESRGroups(Molecule.cESRTypeOr);

		boolean[] isFragmentMember = new boolean[sourceMol.mAllAtoms];
		int graphAtom[] = new int[sourceMol.mAllAtoms];

		graphAtom[0] = rootAtom;
		isFragmentMember[rootAtom] = true;
		atomMap[rootAtom] = sourceMol.copyAtom(this, rootAtom, esrGroupCountAND, esrGroupCountOR);

		int current = 0;
		int highest = 0;
	 	while (current <= highest) {
			for (int i=0; i<sourceMol.getAllConnAtoms(graphAtom[current]); i++) {
				int candidate = sourceMol.mConnAtom[graphAtom[current]][i];
				if (!isFragmentMember[candidate]) {
					graphAtom[++highest] = candidate;
					isFragmentMember[candidate] = true;
					atomMap[candidate] = sourceMol.copyAtom(this, candidate, esrGroupCountAND, esrGroupCountOR);
					}
				}
			current++;
			}

		for (int bond=0; bond<sourceMol.mAllBonds; bond++)
		   	if (isFragmentMember[sourceMol.mBondAtom[0][bond]])
		   		sourceMol.copyBond(this, bond, esrGroupCountAND, esrGroupCountOR, atomMap, false);

		renumberESRGroups(cESRTypeAnd);
		renumberESRGroups(cESRTypeOr);

		mValidHelperArrays = cHelperNone;
		}


	/**
	 * Returns an array of all atoms for which a path of bonds leads to rootAtom
	 * not considering metal ligand bonds.
	 * @param rootAtom
	 * @return atoms being in the same fragment as rootAtom
	 */
	public int[] getFragmentAtoms(int rootAtom) {
		return getFragmentAtoms(rootAtom, false);
		}


	/**
	 * Returns an array of all atoms for which a path of bonds leads to rootAtom.
	 * Metal ligand bonds may or may not be considered a connection.
	 * @param rootAtom
	 * @param considerMetalBonds
	 * @return atoms being in the same fragment as rootAtom
	 */
	public int[] getFragmentAtoms(int rootAtom, boolean considerMetalBonds) {
		ensureHelperArrays(cHelperNeighbours);

		boolean[] isFragmentMember = new boolean[mAllAtoms];
		int graphAtom[] = new int[mAllAtoms];

		graphAtom[0] = rootAtom;
		isFragmentMember[rootAtom] = true;
		int current = 0;
		int highest = 0;
		int fragmentMembers = 1;
	 	while (current <= highest) {
			int connAtoms = considerMetalBonds ? getAllConnAtomsPlusMetalBonds(graphAtom[current])
											   : mAllConnAtoms[graphAtom[current]];
			for (int i=0; i<connAtoms; i++) {
				int candidate = mConnAtom[graphAtom[current]][i];
				if (!isFragmentMember[candidate]) {
					graphAtom[++highest] = candidate;
					isFragmentMember[candidate] = true;
					fragmentMembers++;
					}
				}
			current++;
			}

		int[] fragmentMember = new int[fragmentMembers];
		fragmentMembers = 0;
		for (int atom=0; atom<mAllAtoms; atom++)
			if (isFragmentMember[atom])
				fragmentMember[fragmentMembers++] = atom;

		return fragmentMember;
		}


	/**
	 * Locates all disconnected fragments in the molecule and assigns
	 * fragment numbers (starting from 0) to all atoms. Individual bonds can be marked to
	 * be skipped, i.e. to be treated as non-existing.
	 * Metal ligand bonds may or may not be considered a connection.
	 * @param fragmentNo array not smaller than getAllAtoms()
	 * @param neglectBond array not smaller than getAllBonds()
	 * @param considerMetalBonds
	 * @return array of atom's fragments indexes
	 */
	public int getFragmentNumbers(int[] fragmentNo, boolean[] neglectBond, boolean considerMetalBonds) {
		ensureHelperArrays(cHelperNeighbours);

		for (int atom=0; atom<mAllAtoms; atom++)
			fragmentNo[atom] = -1;

		int fragments = 0;
		for (int atom=0; atom<mAllAtoms; atom++) {
			if (fragmentNo[atom] == -1) {
				fragmentNo[atom] = fragments;
				int graphAtom[] = new int[mAllAtoms];
				graphAtom[0] = atom;
				int current = 0;
				int highest = 0;
				while (current <= highest) {
					int connAtoms = considerMetalBonds ? getAllConnAtomsPlusMetalBonds(graphAtom[current])
													   : mAllConnAtoms[graphAtom[current]];
					for (int i=0; i<connAtoms; i++) {
						int candidate = mConnAtom[graphAtom[current]][i];
						if (fragmentNo[candidate] == -1
						 && !neglectBond[mConnBond[graphAtom[current]][i]]) {
							graphAtom[++highest] = candidate;
							fragmentNo[candidate] = fragments;
							}
						}
					current++;
					}
				fragments++;
				}
			}
		return fragments;
		}


	/**
	 * Locates all unconnected fragments in the Molecule and assigns fragment indexes
	 * for every atom starting with 0. Optionally the fragment detection may be restricted to
	 * those atoms that have been previously marked with setAtomMarker(). In that case
	 * non-marked atoms receive the fragment number -1 and are not considered a connection between
	 * marked atoms potentially causing two marked atoms to end up in different fragments, despite
	 * sharing the same fragment.
	 * Metal ligand bonds may or may not be considered a connection.
	 * @param fragmentNo array at least mAllAtoms big to receive atom fragment indexes
	 * @param markedAtomsOnly if true, then only atoms marked with setAtomMarker() are considered
	 * @param considerMetalBonds
	 * @return number of disconnected fragments
	 */
	public int getFragmentNumbers(int[] fragmentNo, boolean markedAtomsOnly, boolean considerMetalBonds) {
		ensureHelperArrays(cHelperNeighbours);

		for (int atom=0; atom<mAllAtoms; atom++)
			fragmentNo[atom] = -1;

		int fragments = 0;
		for (int atom=0; atom<mAllAtoms; atom++) {
			if (fragmentNo[atom] == -1
			 && (!markedAtomsOnly || isMarkedAtom(atom))) {
				fragmentNo[atom] = fragments;
				int graphAtom[] = new int[mAllAtoms];
				graphAtom[0] = atom;
				int current = 0;
				int highest = 0;
				while (current <= highest) {
					int connAtoms = considerMetalBonds ? getAllConnAtomsPlusMetalBonds(graphAtom[current])
													   : mAllConnAtoms[graphAtom[current]];
					for (int i=0; i<connAtoms; i++) {
						int candidate = mConnAtom[graphAtom[current]][i];
						if (fragmentNo[candidate] == -1
						 && (!markedAtomsOnly || isMarkedAtom(candidate))) {
							graphAtom[++highest] = candidate;
							fragmentNo[candidate] = fragments;
							}
						}
					current++;
					}
				fragments++;
				}
			}
		return fragments;
		}


	/**
	 * Removes all unconnected fragments except for the largest one.
	 * If small fragments were removed, then canonizeCharge() is called to
	 * neutralize charges after potential removal of counter ions.
	 * Metal ligand bonds are not considered a connection.
	 * @return atom mapping from old to new index; null if no fragments were removed
	 */
	public int[] stripSmallFragments() {
		return stripSmallFragments(false);
		}


	/**
	 * Removes all unconnected fragments except for the largest one.
	 * If small fragments were removed, then canonizeCharge() is called to
	 * neutralize charges after potential removal of counter ions.
	 * Metal ligand bonds may or may not be considered a connection.
	 * @param considerMetalBonds
	 * @return atom mapping from old to new index; null if no fragments were removed
	 */
	public int[] stripSmallFragments(boolean considerMetalBonds) {
		int[] fragmentNo = new int[mAllAtoms];
		int fragmentCount = getFragmentNumbers(fragmentNo, false, considerMetalBonds);

		if (fragmentCount <= 1)
			return null;

		int[] fragmentSize = new int[fragmentCount];
		for (int atom=0; atom<mAtoms; atom++)
			fragmentSize[fragmentNo[atom]]++;

		int largestFragment = 0;
		int largestSize = fragmentSize[0];
		for (int i=1; i<fragmentCount; i++) {
			if (largestSize < fragmentSize[i]) {
				largestSize = fragmentSize[i];
				largestFragment = i;
				}
			}

		for (int atom=0; atom<mAllAtoms; atom++)
			if (fragmentNo[atom] != largestFragment)
				mAtomicNo[atom] = -1;				// mark for delete

		for (int bond=0; bond<mAllBonds; bond++)
			if ((!considerMetalBonds && mBondType[bond] == cBondTypeMetalLigand)
			 || fragmentNo[mBondAtom[0][bond]] != largestFragment)
				mBondType[bond] = cBondTypeDeleted;	// mark for delete

		int[] atomMap = compressMolTable();
		mValidHelperArrays = cHelperNone;

		try { canonizeCharge(true, true); } catch (Exception e) {}

		return atomMap;
		}


	/**
	 * Starting from startAtom this method locates a system of annelated or bridged ring systems
	 * with all members bonds being a ring bond. Detected member atoms and bonds are flagged
	 * accordingly.
	 * @param startAtom
	 * @param aromaticOnly if set then only aromatic atoms and bonds are considered
	 * @param isMemberAtom
	 * @param isMemberBond
	 */
	public void findRingSystem(int startAtom, boolean aromaticOnly, boolean[] isMemberAtom, boolean[] isMemberBond) {
		ensureHelperArrays(cHelperRings);
		
		if (!isRingAtom(startAtom) || (aromaticOnly && !isAromaticAtom(startAtom)))
			return;

		int[] graphAtom = new int[mAtoms];
			
		graphAtom[0] = startAtom;
		isMemberAtom[startAtom] = true;

		int current = 0;
		int highest = 0;
	 	while (current <= highest) {
			for (int i = 0; i< mConnAtoms[graphAtom[current]]; i++) {
				int candidateBond = mConnBond[graphAtom[current]][i];
				if (!isMemberBond[candidateBond]
				 && isRingBond(candidateBond)
				 && (!aromaticOnly || isAromaticBond(candidateBond))) {
					isMemberBond[candidateBond] = true;
					int candidateAtom = mConnAtom[graphAtom[current]][i];
					if (!isMemberAtom[candidateAtom]) {
						isMemberAtom[candidateAtom] = true;
						graphAtom[++highest] = candidateAtom;
						}
					}
				}
			current++;
			}
		}

	/**
	 * Determines all atoms of the substituent attached to coreAtom and starting
	 * with firstAtom. If isMemberAtom!=null, then all substituent member atoms
	 * will have the the respective index being flagged upon return. This includes
	 * firstAtom and excludes coreAtom.
	 * If substituent!=null, then it will contain the substituent as Molecule.
	 * At the position of the coreAtom substituent will contain a wildcard atom.
	 * If substituent!=null and atomMap!=null then atomMap receives atom index mapping from
	 * this to substituent with non-member atoms being -1.
	 * Returns -1 and an empty substituent if coreAtom and firstAtom share a ring
	 * @param coreAtom the atom to which the substituent is connected
	 * @param firstAtom the substituent's atom that is connected to coreAtom
	 * @param isMemberAtom may be null, otherwise set to contain atom membership mask
	 * @param substituent may be null, otherwise set to contain the substituent
	 * @param atomMap null or int[] not smaller than this.getAllAtoms()
	 * @return substituent atom count not counting coreAtom; -1 if coreAtom and firstAtom share a ring
	 */
	public int getSubstituent(int coreAtom, int firstAtom, boolean[] isMemberAtom, ExtendedMolecule substituent, int[] atomMap) {
		ensureHelperArrays(cHelperNeighbours);

		if (substituent != null) {
			substituent.clear();
			substituent.mIsFragment = false;
			}

		int[] graphAtom = new int[mAllAtoms];
		if (isMemberAtom == null)
			isMemberAtom = new boolean[mAllAtoms];
		else
			java.util.Arrays.fill(isMemberAtom, false);
			
		graphAtom[0] = coreAtom;
		graphAtom[1] = firstAtom;
		isMemberAtom[coreAtom] = true;
		isMemberAtom[firstAtom] = true;
		int current = 1;
		int highest = 1;
	 	while (current <= highest) {
//			int connAtoms = getAllConnAtomsPlusMetalBonds(graphAtom[current]);
//		    for (int i=0; i<connAtoms; i++) {     the default for graph methods is to neglect dative bonds; TLS 27 May 2020
// If dative bonds must be considered on a certain context, the we need to introduce a flag e.g. as in getFragmentAtoms()
			for (int i=0; i<mAllConnAtoms[graphAtom[current]]; i++) {
				int candidate = mConnAtom[graphAtom[current]][i];
				if (candidate == coreAtom) {
					if (current != 1)
						return -1;
					}
				if (!isMemberAtom[candidate]) {
					isMemberAtom[candidate] = true;
					graphAtom[++highest] = candidate;
					}
				}
			current++;
			}

	 	if (substituent != null) {
	 		if (atomMap == null)
	 			atomMap = new int[isMemberAtom.length];
	 		copyMoleculeByAtoms(substituent, isMemberAtom, false, atomMap);
			substituent.changeAtom(atomMap[coreAtom], 0, 0, -1, 0);
	 		}

		isMemberAtom[coreAtom] = false;
		return highest;
		}

	/**
	 * Counts the number of atoms of the substituent connected to coreAtom
	 * defined by firstAtom and not including the coreAtom.
	 * @param coreAtom
	 * @param firstAtom
	 * @return atom count of substituent or -1 if coreAtom and firstAtom are in the same ring
	 */
	public int getSubstituentSize(int coreAtom, int firstAtom) {
		ensureHelperArrays(cHelperNeighbours);
	
		int[] graphAtom = new int[mAtoms];
		boolean[] isMember = new boolean[mAtoms];
		graphAtom[0] = coreAtom;
		graphAtom[1] = firstAtom;
		isMember[coreAtom] = true;
		isMember[firstAtom] = true;
		int current = 1;
		int highest = 1;
		while (current <= highest) {
			for (int i = 0; i< mConnAtoms[graphAtom[current]]; i++) {
				int candidate = mConnAtom[graphAtom[current]][i];
				if (candidate == coreAtom) {
					if (current != 1)
						return -1;
					}
				if (!isMember[candidate]) {
					isMember[candidate] = true;
					graphAtom[++highest] = candidate;
					}
				}
			current++;
			}
		return highest;
		}

	/**
	 * Whether an atom may be considered to carry implicit hydrogen atoms depends
	 * on the atomicNo of that atom. Aluminum and all non/metal atoms except the
	 * nobel gases and except hydrogen itself are considered to carry implicit hydrogens
	 * to fill up their unoccupied valences. Atoms with an assigned unusual valence always
	 * support implicit hydrogens independent of their atomicNo.
	 * @param atom
	 * @return true if this atom's unoccupied valences are considered to be implicit hydrogens
	 */
	public boolean supportsImplicitHydrogen(int atom) {
		if ((mAtomFlags[atom] & cAtomFlagsValence) != 0)
			return true;
		if (mAtomicNo[atom] == 1)
			return false;
		return isOrganicAtom(atom)
				|| mAtomicNo[atom] == 13	// Al
				|| mAtomicNo[atom] >= 171;	// amino acids
		}

	/**
	 * Calculates and return the number of implicit hydrogens at atom.
	 * If atom is itself a hydrogen atom, a metal except Al, or a noble gas,
	 * then 0 is returned. For all other atom kinds the number of
	 * implicit hydrogens is basically the lowest typical valence that is compatible
	 * with the occupied valence, minus the occupied valence corrected by atom charge
	 * and radical state.
	 * @param atom
	 * @return
	 */
	public int getImplicitHydrogens(int atom) {
		if (mIsFragment
		 && (mAtomQueryFeatures[atom] & cAtomQFNoMoreNeighbours) == 0)
			return 0;

		// H, metals except Al, noble gases don't have implicit hydrogens
		if (!supportsImplicitHydrogen(atom))
			return 0;

		ensureHelperArrays(cHelperNeighbours);

		int occupiedValence = 0;
		for (int i=0; i<mAllConnAtoms[atom]; i++)
			occupiedValence += mConnBondOrder[atom][i];

		if (mIsFragment) {
			int delocalizedBonds = 1;
			for (int i = 0; i< mConnAtoms[atom]; i++)
				if (mBondType[mConnBond[atom][i]] == cBondTypeDelocalized)
					delocalizedBonds++;
			occupiedValence += delocalizedBonds >> 1;
			}

		occupiedValence -= getElectronValenceCorrection(atom, occupiedValence);
		int maxValence = getAtomAbnormalValence(atom);
		if (maxValence == -1) {
			byte[] valenceList = Molecule.getAllowedValences(mAtomicNo[atom]);
			maxValence = valenceList[0];
			for (int i=1; (maxValence<occupiedValence) && (i<valenceList.length); i++)
				maxValence = valenceList[i];
			}

		return Math.max(0, maxValence - occupiedValence);
		}

	/**
	 * @param atom
	 * @return number of explicit plain hydrogen atoms (does not include 2H, 3H, custom labelled H, etc)
	 */
	public int getExplicitHydrogens(int atom) {
		return mAllConnAtoms[atom] - mConnAtoms[atom];
		}

	/**
	 * Calculates a rounded mass of the molecule
	 * @return
	 */
	public int getMolweight() {
		ensureHelperArrays(cHelperNeighbours);
		int molweight = 0;
		for (int atom=0; atom<mAllAtoms; atom++) {
			int mass = mAtomMass[atom] != 0 ? mAtomMass[atom] : cRoundedMass[mAtomicNo[atom]];
			molweight += mass + getImplicitHydrogens(atom) * cRoundedMass[1];
			if (mAtomicNo[atom] >= 171 && mAtomicNo[atom] <= 190) {
				int connAtoms = mAllConnAtoms[atom];
				if (connAtoms > 2)
					molweight -= (connAtoms - 2) * cRoundedMass[1];
				}
			}

		return molweight;
		}

	/**
	 * Simple method to calculate rotatable bonds. This method counts all single
	 * bonds provided that they<br>
	 * - are not a terminal bond<br>
	 * - are not part of a ring<br>
	 * - are not an amide bond<br>
	 * - are not the second of two equivalent bonds next to the same triple bond<br>
	 * @return
	 */
	public int getRotatableBondCount() {
		int rCount = 0;
		ensureHelperArrays(Molecule.cHelperRings);
		for (int bond=0; bond<mBonds; bond++) {
			if (getBondOrder(bond) == 1 && !isRingBond(bond)) {
				boolean isRotatable = true;
				for (int i = 0; i < 2; i++) {
					int atom1 = mBondAtom[i][bond];
					if (mConnAtoms[atom1] == 1) {
						isRotatable = false;
						break;  // terminal bond
						}

					if (mAtomicNo[atom1] == 7 && !isAromaticAtom(atom1)) {
						int atom2 = mBondAtom[1 - i][bond];
						for (int j = 0; j < mConnAtoms[atom2]; j++) {
							int connAtom = mConnAtom[atom2][j];
							int connBond = mConnBond[atom2][j];
							if (connBond != bond
									&& getBondOrder(connBond) > 1
									&& !isAromaticAtom(connAtom)
									&& isElectronegative(connAtom)) {
								isRotatable = false;
								break;  // amid bond
								}
							}
						}
					}

				if (isRotatable && !isPseudoRotatableBond(bond))
					rCount++;
				}
			}
		return rCount;
		}

	/**
	 * In a consecutive sequence of sp-hybridized atoms multiple single bonds
	 * cause redundant torsions. Only that single bond with the smallest bond index
	 * is considered really rotatable; all other single bonds are pseudo rotatable.
	 * If one/both end(s) of the sp-atom sequence doesn't carry atoms
	 * outside of the straight line then no bond is considered rotatable.
	 * A simple terminal single bond
	 * @param bond
	 * @return true, if this bond is not considered rotatable because of a redundancy
	 */
	public boolean isPseudoRotatableBond(int bond) {
		if (getBondOrder(bond) != 1)
			return false;

		for (int i=0; i<2; i++) {
			int atom = mBondAtom[i][bond];
			int rearAtom = mBondAtom[1-i][bond];

			while (mPi[atom] == 2
				&& mConnAtoms[atom] == 2
				&& mAtomicNo[atom] < 10) {
				for (int j=0; j<2; j++) {
					int connAtom = mConnAtom[atom][j];
					if (connAtom != rearAtom) {
						if (mConnAtoms[connAtom] == 1)
							return true;

						int connBond = mConnBond[atom][j];
						if (getBondOrder(connBond) == 1
						 && connBond < bond)
							return true;

						rearAtom = atom;
						atom = connAtom;
						break;
						}
					}
				}

			if (mConnAtoms[atom] == 1)
				return true;
			}

		return false;
		}

	public int getAromaticRingCount() {
		ensureHelperArrays(cHelperRings);
		int count = 0;
		for (int i=0; i<mRingSet.getSize(); i++)
			if (mRingSet.isAromatic(i))
				count++;
		return count;
		}

	/**
	 * Calculates the number of independent rings of which 'atom' is a member.
	 * Any combination of two connected atoms to 'atom' is used for:
	 * - finding the shortest path connecting these two neighbors avoiding 'atom'
	 * - if such a path exists and at least one bonds of that path is not a member
	 *   of a path found earlier then count this path as an independent ring closure.
	 * @param atom
	 * @param maxRingSize
	 * @return number of independent rings
	 */
	public int getAtomRingCount(int atom, int maxRingSize) {
		ensureHelperArrays(cHelperRings);
		boolean[] bondTouched = new boolean[mBonds];
		boolean[] neglectBond = new boolean[mBonds];
		int[] ringAtom = new int[mAtoms];
		int count = 0;
		for (int i=1; i<mConnAtoms[atom]; i++) {
			int bond1 = mConnBond[atom][i];
			if (isRingBond(bond1)) {
				for (int j=0; j<i; j++) {
					int bond2 = mConnBond[atom][j];
					if (isRingBond(bond2)) {
						neglectBond[bond1] = true;
						neglectBond[bond2] = true;
						int pathLength = getPath(ringAtom, mConnAtom[atom][i], mConnAtom[atom][j], maxRingSize-2, neglectBond);
						neglectBond[bond1] = false;
						neglectBond[bond2] = false;
						if (pathLength != -1) {
							boolean isIndependentRing = false;
							int[] pathBond = new int[pathLength];
							getPathBonds(ringAtom, pathBond, pathLength);
							for (int k=0; k<pathLength; k++) {
								if (!bondTouched[pathBond[k]]) {
									bondTouched[pathBond[k]] = true;
									isIndependentRing = true;
									}
								}
							if (isIndependentRing)
								count++;
							}
						}
					}
				}
			}
		return count;
		}


	/**
	 * @return a RingCollection object, which contains a total set of small rings
	 */
	public RingCollection getRingSet() {
		ensureHelperArrays(cHelperRings);
		return mRingSet;
		}


	/**
	 * 
	 * @return a RingCollection object without aromaticity information
	 */
	
	public RingCollection getRingSetSimple() {
		ensureHelperArrays(cHelperRingsSimple);
		return mRingSet;
		}


	/**
	 * Locates that single bond which is the preferred one to be converted into up/down bond
	 * in order to define the atom chirality. This method does not consider, whether any
	 * bond is a stereo bond already.
	 * @param atom parity carrying atom, i.e. a tetrahedral stereocenter or central allene atom
	 * @return preferred bond or -1, if no single bond existing
	 */
	public int getAtomPreferredStereoBond(int atom) {
		ensureHelperArrays(cHelperRings);
		if (mPi[atom] == 2 && mConnAtoms[atom] == 2)
			return preferredAlleneStereoBond(atom, false);
		else
			return preferredTHStereoBond(atom, false);
		}


	/**
	 * Locates that single bond which is the preferred one to be converted into up/down bond
	 * in order to define the bond chirality. This method does not consider, whether any
	 * bond is a stereo bond already.
	 * @param bond BINAP type of chirality bond
	 * @return preferred bond or -1, if no single bond existing
	 */
	public int getBondPreferredStereoBond(int bond) {
		return preferredBinapStereoBond(bond);
		}


	private int getStereoBondScore(int bond, int atom) {
			// score used to select one bond to be the stereo bond
		if (getBondOrder(bond) != 1)
			return 0;

		return 16 - mAllConnAtoms[atom]
				  + (((mBondType[bond] & cBondTypeMaskStereo) == 0 || mBondAtom[0][bond] != atom) ? 32768 : 0)
				  + ((mAtomicNo[atom] == 1) ? 4096 : 0)
				  + ((mAllConnAtoms[atom] == 1) ? 2048 : 0)
				  + ((getAtomParity(atom) == 0) ? 1024 : 0)
				  + ((!isRingBond(bond)) ? 512 : 0)
				  + ((mAtomicNo[atom] != 6) ? 256 : 0);
		}


	/**
	 * @param atom
	 * @return whether the atom is in an allylic/benzylic position
	 */
	public boolean isAllylicAtom(int atom) {
		return (mAtomFlags[atom] & cAtomFlagAllylic) != 0;
		}


	public boolean isAromaticAtom(int atom) {
		return (atom < mAtoms) ? mRingSet.isAromaticAtom(atom) : false;
		}


	public boolean isHeteroAromaticAtom(int atom) {
		return (atom < mAtoms) ? mRingSet.isHeteroAromaticAtom(atom) : false;
	}


	/**
	 * @param atom
	 * @return whether the atom is a member of a delocalized ring (subset of aromatic rings)
	 */
	public boolean isDelocalizedAtom(int atom) {
		return (atom < mAtoms) ? mRingSet.isDelocalizedAtom(atom) : false;
	}


	public boolean isAromaticBond(int bond) {
		return (bond < mBonds) ? mRingSet.isAromaticBond(bond) : false;
	}


	public boolean isHeteroAromaticBond(int bond) {
		return (bond < mBonds) ? mRingSet.isHeteroAromaticBond(bond) : false;
	}


	/**
	 * A bond is considered delocalized, if it has different bond orders in
	 * different, but energetically equivalent mesomeric structures. Bonds in aromatic 6-membered
	 * rings typically are delocalized, while those in uncharged 5-membered aromatic rings are not.
	 * Indole has 6 delocalized bonds.
	 * This method also returns true, if the molecule is a fragment and if the bond is explicitly
	 * defined to be delocalized.
	 * @param bond
	 * @return
	 */
	@Override
	public boolean isDelocalizedBond(int bond) {
		return (bond < mBonds) ? mRingSet.isDelocalizedBond(bond) || mBondType[bond] == cBondTypeDelocalized : false;
		}


	public boolean isRingAtom(int atom) {
		return (mAtomFlags[atom] & cAtomFlagsRingBonds) != 0;
		}


	public boolean isRingBond(int bnd) {
		return (mBondFlags[bnd] & cBondFlagRing) != 0;
		}


	/**
	 * @param atom
	 * @return whether atom is a member of a ring not larger than 7 atoms
	 */
	public boolean isSmallRingAtom(int atom) {
		return (mAtomFlags[atom] & cAtomFlagSmallRing) != 0;
		}


	/**
	 * @param bond
	 * @return whether bond is a member of a ring not larger than 7 atoms
	 */
	public boolean isSmallRingBond(int bond) {
		return (mBondFlags[bond] & cBondFlagSmallRing) != 0;
		}


	/**
	 * @param atom
	 * @return whether atom has a neighbor that is connected through a double/triple bond to a hetero atom
	 */
	public boolean isStabilizedAtom(int atom) {
		return (mAtomFlags[atom] & cAtomFlagStabilized) != 0;
		}


	public int getAtomRingBondCount(int atom) {
		int flags = (mAtomFlags[atom] & cAtomFlagsRingBonds);
		return (flags == 0) ? 0
			 : (flags == cAtomFlags2RingBonds) ? 2
			 : (flags == cAtomFlags3RingBonds) ? 3 : 4;
		}


	public String getChiralText() {
		return null;
		}


	/**
	 * Checks whether at least one of the connected bonds is a stereo bond.
	 * If atom is the central atom of an allene, then its direct neighbours
	 * are checked, whether one of them has a stereo bond.
	 * @param atom
	 * @return the stereo bond or -1 if not found
	 */
	public int getStereoBond(int atom) {
		ensureHelperArrays(cHelperNeighbours);
		if (mConnAtoms[atom] == 2
		 && mConnBondOrder[atom][0] == 2
		 && mConnBondOrder[atom][1] == 2) {
			for (int i=0; i<2; i++)
				for (int j=0; j<mAllConnAtoms[mConnAtom[atom][i]]; j++)
					if (isStereoBond(mConnBond[mConnAtom[atom][i]][j], mConnAtom[atom][i]))
						return mConnBond[mConnAtom[atom][i]][j];
			}
		else {
			for (int i=0; i<mAllConnAtoms[atom]; i++)
				if (isStereoBond(mConnBond[atom][i], atom))
					return mConnBond[atom][i];
			}

		return -1;
		}


	/**
	 * Atom stereo parities and bond E/Z-parities are properties that are usually perceived
	 * from up/down-bonds and atom coordinates, respectively. This is done during the helper
	 * array calculation triggered by ensureHelperArrays(cHelperParities).<br>
	 * This method tells the molecule that current atom/bond parities are valid, even if the
	 * stereo perception not has been performed. In addition to the stereo parities one may
	 * declare CIP parities and/or symmetry ranks also to be valid (helperStereoBits != 0).
	 * setParitiesValid(0) should be called if no coordinates are available but the parities are valid
	 * nevertheless, e.g. after the IDCodeParser has parsed an idcode without coordinates.
	 * (Note: After idcode parsing unknown stereo centers have parities cAtomParityNone
	 * instead of cAtomParityUnknown. Thus, calling isStereoCenter(atom) returns false!!!)
	 * Declaring parities valid prevents the Canonizer to run the stereo recognition again when
	 * ensureHelperArrays(cHelperParities or higher) is called.<br>
	 * May also be called after filling valences with explicit hydrogen atoms, which have no
	 * coordinates, to tell the molecule that the earlier created stereo flags are still valid.
	 * @param helperStereoBits 0 or combinations of cHelperBitCIP,cHelperBitSymmetry...,cHelperBitIncludeNitrogenParities
	 */
	public void setParitiesValid(int helperStereoBits) {
		mValidHelperArrays |= (cHelperBitsStereo & (cHelperBitParities | helperStereoBits));
		}


	/**
	 * This converts one single bond per parity into a stereo up/down bond to
	 * correctly reflect the given parity. This works for tetrahedral and
	 * allene atom parities as well as for BINAP type of bond parities.
	 * Should only be called with valid TH and EZ parities and valid coordinates,
	 * e.g. after idcode parsing with coordinates or after coordinate generation.
	 */
	public void setStereoBondsFromParity() {
		assert((mValidHelperArrays & cHelperBitParities) != 0);
		ensureHelperArrays(cHelperRings); // in case we miss ring and neighbour information

		// Convert all stereo bonds into single bonds
		for (int bond=0; bond<mBonds; bond++)
			if (isStereoBond(bond))
				mBondType[bond] = cBondTypeSingle;

		for (int atom=0; atom<mAtoms; atom++)
			setStereoBondFromAtomParity(atom);

		for (int bond=0; bond<mBonds; bond++)
			setStereoBondFromBondParity(bond);

		for (int bond=0; bond<mBonds; bond++)
			if (mBondType[bond] == cBondTypeDouble
			 && getBondParity(bond) == Molecule.cBondParityUnknown)
				mBondType[bond] = cBondTypeCross;

		assert((mValidHelperArrays & cHelperBitParities) != 0);
		}


	/**
	 * Converts any stereo bond being attached with the pointed tip to the given potential stereocenter
	 * (TH or allene) into a single bond. For allenic stereo centers, stereo bonds at the allene end atoms
	 * are converted. Atoms with PI electrons are not touched, unless they are P,S,...
	 * @param atom
	 */
	public void convertStereoBondsToSingleBonds(int atom) {
		if (mPi[atom] == 2 && mConnAtoms[atom] == 2 && mConnBondOrder[atom][0] == 2) {
			for (int i=0; i<2; i++) {
				int alleneEnd = findAlleneEndAtom(atom, mConnAtom[atom][i]);
				if (alleneEnd != -1) {
					for (int j=0; j<mConnAtoms[alleneEnd];j++) {
						int connBond = mConnBond[alleneEnd][j];
						if (isStereoBond(connBond) && mBondAtom[0][connBond] == alleneEnd)
							mBondType[mConnBond[alleneEnd][j]] = cBondTypeSingle;
						}
					}
				}
			return;
			}

		if (mPi[atom] == 0 || mAtomicNo[atom] >= 15) {
			for (int i=0; i<mAllConnAtoms[atom]; i++) {
				int connBond = mConnBond[atom][i];
				if (isStereoBond(connBond, atom) && mBondAtom[0][connBond] == atom)
					mBondType[connBond] = cBondTypeSingle;
				}
			}
		}


	public void setStereoBondFromAtomParity(int atom) {
		convertStereoBondsToSingleBonds(atom);

			// set an optimal bond to up/down to reflect the atom parity
		if (getAtomParity(atom) == Molecule.cAtomParityNone
		 || getAtomParity(atom) == Molecule.cAtomParityUnknown)
			return;

		if (mPi[atom] == 2 && mConnAtoms[atom] == 2) {
			setAlleneStereoBondFromParity(atom);
			return;
			}

		if (mConnAtoms[atom] < 3 || mConnAtoms[atom] > 4) {
			setAtomParity(atom, cAtomParityNone, false);
			return;
			}

		int allConnAtoms = mAllConnAtoms[atom];

		// We may have a rare case without any single bond (e.g. O=S(=NH)=NMe) with parities assigned from 3D-coords
		boolean singleBondFound = false;
		for (int i=0; i<allConnAtoms; i++) {
			if (getBondOrder(mConnBond[atom][i]) == 1) {
				singleBondFound = true;
				break;
				}
			}
		if (!singleBondFound)
			return;

		int[] sortedConnMap = getSortedConnMap(atom);

		double angle[] = new double[allConnAtoms];
		for (int i=0; i<allConnAtoms; i++)
			angle[i] = getBondAngle(mConnAtom[atom][sortedConnMap[i]], atom);

		for (int i=0; i<allConnAtoms; i++)
			if (mBondAtom[0][mConnBond[atom][i]] == atom
			 && getBondOrder(mConnBond[atom][i]) == 1)
				mBondType[mConnBond[atom][i]] = cBondTypeSingle;

		// Don't allow Fisher projection for large rings, which may have multiple almost vertical bonds
		if (getAtomRingSize(atom) <= FISCHER_PROJECTION_RING_LIMIT
		 && setFisherProjectionStereoBondsFromParity(atom, sortedConnMap, angle))
			return;

		int preferredBond = preferredTHStereoBond(atom, true);

		if (mBondAtom[0][preferredBond] != atom) {
			mBondAtom[1][preferredBond] = mBondAtom[0][preferredBond];
			mBondAtom[0][preferredBond] = atom;
			}

		int preferredBondIndex = -1;
		for (int i=0; i<allConnAtoms; i++) {
			if (preferredBond == mConnBond[atom][sortedConnMap[i]]) {
				preferredBondIndex = i;
				break;
				}
			}

		final int[][] up_down = { { 2,1,2,1 },	// stereobond type to
								  { 1,2,2,1 },	// achieve parity = 1
								  { 1,1,2,2 },	// first dimension:
								  { 2,1,1,2 },	// one of the 6 angle orders
								  { 2,2,1,1 },	// second dimension:
								  { 1,2,1,2 } };// index of ConnAtom bearing stereobond

/*		int preferredBondIndex = -1;
		int sortedConn[] = new int[4];
		for (int i=0; i<mAllConnAtoms[atom]; i++) {
			int lowestConnAtom = Integer.MAX_VALUE;
			int lowestConnIndex = -1;
			for (int j=0; j<mAllConnAtoms[atom]; j++) {
				if ((i == 0 || mConnAtom[atom][j] > sortedConn[i-1]) && lowestConnAtom > mConnAtom[atom][j]) {
					lowestConnAtom = mConnAtom[atom][j];
					lowestConnIndex = j;
					}
				}

			sortedConn[i] = lowestConnAtom;
			if (mConnBond[atom][lowestConnIndex] == preferredBond)
				preferredBondIndex = i;
			}	*/

		for (int i=1; i<allConnAtoms; i++)
			if (angle[i] < angle[0])
				angle[i] += Math.PI*2;

		int bondType;
		if (allConnAtoms == 3) {
				// new handling!!! TLS 17.Oct.2005
				// Originally the parity depended solely on clockwise/anti-clockwise orientation
				// of three (priority sorted) angles, which included the angle of the stereo bond.
				// Now to handle strained rings with 'invalid' projections in an expected way
				// (all three bonds are within 180 degrees and the middle bond is the stereo bond
				// then treat this as if the stereo bond would be drawn 180 degrees rotated)
				// the angle of the stereobond is not considered anymore and the parity depends
				// solely on the question if the angle difference between higher priority bond
				// to lower priority bond is less or larger than 180 degrees.
			boolean inverted = false;
			switch (preferredBondIndex) {
			case 0:
				inverted = (((angle[1] < angle[2]) && (angle[2] - angle[1] < Math.PI))
						 || ((angle[1] > angle[2]) && (angle[1] - angle[2] > Math.PI)));
				break;
			case 1:
				inverted = (angle[2] - angle[0] > Math.PI);
				break;
			case 2:
				inverted = (angle[1] - angle[0] < Math.PI);
				break;
				}

			bondType = ((getAtomParity(atom) == cAtomParity1) ^ inverted) ?
						cBondTypeUp : cBondTypeDown;

/*				// original handling where parity depended solely on the clockwise/anti-clockwise of three angles
			bondType = ((getAtomParity(atom) == cAtomParity1)
					  ^ (angle[1] > angle[2])) ?
						cBondTypeDown : cBondTypeUp;	*/
			}
		else {
			int order = 0;
			if		(angle[1] <= angle[2] && angle[2] <= angle[3]) order = 0;
			else if (angle[1] <= angle[3] && angle[3] <= angle[2]) order = 1;
			else if (angle[2] <= angle[1] && angle[1] <= angle[3]) order = 2;
			else if (angle[2] <= angle[3] && angle[3] <= angle[1]) order = 3;
			else if (angle[3] <= angle[1] && angle[1] <= angle[2]) order = 4;
			else if (angle[3] <= angle[2] && angle[2] <= angle[1]) order = 5;
			bondType = ((getAtomParity(atom) == cAtomParity1)
					  ^ (up_down[order][preferredBondIndex] == 1)) ?
						cBondTypeDown : cBondTypeUp;
			}

		mBondType[preferredBond] = bondType;
		}


	private boolean setFisherProjectionStereoBondsFromParity(int atom, int[] sortedConnMap, double[] angle) {
		int allConnAtoms = mAllConnAtoms[atom];
		int[] direction = new int[allConnAtoms];
		int parity = getFisherProjectionParity(atom, sortedConnMap, angle, direction);
		if (parity == cAtomParityUnknown)
			return false;

		int bondType = (getAtomParity(atom) == parity) ? cBondTypeUp : cBondTypeDown;
		for (int i=0; i<allConnAtoms; i++) {
			if ((direction[i] & 1) == 1) {
				int bond = mConnBond[atom][sortedConnMap[i]];
				mBondType[bond] = bondType;
				if (mBondAtom[0][bond] != atom) {
					mBondAtom[1][bond] = mBondAtom[0][bond];
					mBondAtom[0][bond] = atom;
					}
				}
			}

		return true;
		}


	/**
	 * If the atom is a stereo center in fisher projection, then its
	 * tetrahedral parity is returned. If the horizontal bonds are plain
	 * single bonds, then they are interpreted as up-bonds.
	 * @param atom the stereo center
	 * @param sortedConnMap map of neighbours sorted by atom index
	 * @param angle bond angles sorted by neighbour atom index
	 * @param direction null or int[] large enough to receive bond directions
	 * @return cAtomParity1,cAtomParity2 or cAtomParityUnknown
	 */
	public int getFisherProjectionParity(int atom, int[] sortedConnMap, double[] angle, int[] direction) {
		if (getAtomRingSize(atom) > FISCHER_PROJECTION_RING_LIMIT)
			return cAtomParityUnknown;

		int allConnAtoms = mAllConnAtoms[atom];
		if (direction == null)
			direction = new int[allConnAtoms];
		if (!getFisherProjectionBondDirections(atom, sortedConnMap, angle, direction))
			return cAtomParityUnknown;

		int horizontalBondType = -1;
		for (int i=0; i<allConnAtoms; i++) {
			if ((direction[i] & 1) == 1) { // horizontal bond
				int bondType = mBondType[mConnBond[atom][sortedConnMap[i]]];
				if (horizontalBondType != -1
				 && horizontalBondType != bondType)
					return cAtomParityUnknown;
				horizontalBondType = bondType;
				}
			}

		// we determine the rotation for the first 3 substituents: cw/ccw
		int index = (Math.abs(direction[0] - direction[1]) == 2) ? 1 : 0;
		int dif = direction[index] - direction[index+1];
		boolean isClockwise = (Math.abs(dif) == 3) ^ (direction[index] < direction[index+1]);
		boolean is4thConnHorizontal = (allConnAtoms == 3 || ((direction[3] & 1) == 1));
		return (isClockwise ^ is4thConnHorizontal) ^ (horizontalBondType == cBondTypeDown) ?
				cAtomParity1 : cAtomParity2;
		}


	/**
	 * Checks whether we have two vertical non-stereo single bonds and
	 * two horizontal stereo single bonds. If these conditions
	 * are met, then an int array is filled defining the directions of all
	 * connected bonds (0:south; 1:west; 2:north; 3:east).
	 * @param atom
	 * @param sortedConnMap map of neighbours sorted by atom index
	 * @param angle bond angles sorted by neighbour atom index
	 * @param direction array large enough to receive bond directions
	 * @return false if fisher projection conditions are not met
	 */
	private boolean getFisherProjectionBondDirections(int atom, int[] sortedConnMap, double[] angle, int[] direction) {
		int allConnAtoms = mAllConnAtoms[atom];
		if (mPi[atom] != 0
		 || isAromaticAtom(atom)
		 || mConnAtoms[atom] < 3
		 || allConnAtoms > 4)
			return false;

		boolean[] isUsed = new boolean[4];
		for (int i=0; i<allConnAtoms; i++) {
			double a = Math.PI*5/4 - angle[i];
			if (Math.abs(Math.PI/4 - (a % (Math.PI/2))) > FISCHER_PROJECTION_LIMIT)
				return false;

			direction[i] = 3 & (int)(a / (Math.PI/2));
			if (isUsed[direction[i]])
				return false;
			isUsed[direction[i]] = true;

			if ((direction[i] & 1) == 0) {  // vertical bond
				if (mBondType[mConnBond[atom][sortedConnMap[i]]] != cBondTypeSingle)
					return false;
				}
			else {
				if (!isStereoBond(mConnBond[atom][sortedConnMap[i]], atom))
					return false;
				}
			}
		return isUsed[0] && isUsed[2];
		}

	
	private void setAlleneStereoBondFromParity(int atom) {
			// find preferred bond to serve as stereobond
		if (mConnAtoms[atom] != 2
		 || mConnBondOrder[atom][0] != 2
		 || mConnBondOrder[atom][1] != 2
		 || mConnAtoms[mConnAtom[atom][0]] < 2
		 || mConnAtoms[mConnAtom[atom][1]] < 2
		 || mPi[mConnAtom[atom][0]] != 1
		 || mPi[mConnAtom[atom][1]] != 1) {
			setAtomParity(atom, cAtomParityNone, false);
			return;
			}

		int preferredBond = -1;
		int preferredAtom = -1;
		int preferredAlleneAtom = -1;
		int oppositeAlleneAtom = -1;
		int bestScore = 0;
		for (int i=0; i<2; i++) {
			int alleneAtom = mConnAtom[atom][i];
			for (int j=0; j<mAllConnAtoms[alleneAtom]; j++) {
				int connAtom = mConnAtom[alleneAtom][j];
				if (connAtom != atom) {
					int connBond = mConnBond[alleneAtom][j];
					int score = getStereoBondScore(connBond,connAtom);
					if (bestScore < score) {
						bestScore = score;
						preferredAtom = connAtom;
						preferredBond = connBond;
						preferredAlleneAtom = alleneAtom;
						oppositeAlleneAtom = mConnAtom[atom][1-i];
						}
					}
				}
			}

		if (preferredAtom == -1)
			return;

		for (int i=0; i<2; i++) {
			int alleneAtom = mConnAtom[atom][i];
			for (int j = 0; j<mAllConnAtoms[alleneAtom]; j++) {
				int connAtom = mConnAtom[alleneAtom][j];
				int connBond = mConnBond[alleneAtom][j];
				if (connAtom != atom
				 && mBondAtom[0][connBond] == alleneAtom)
					mBondType[connBond] = cBondTypeSingle;
				}
			}

		if (mBondAtom[1][preferredBond] != preferredAtom) {
			mBondAtom[0][preferredBond] = mBondAtom[1][preferredBond];
			mBondAtom[1][preferredBond] = preferredAtom;
			}

		int highPriorityAtom = Integer.MAX_VALUE;
		for (int i = 0; i< mConnAtoms[preferredAlleneAtom]; i++) {
			int connAtom = mConnAtom[preferredAlleneAtom][i];
			if ((connAtom != atom) && (highPriorityAtom > connAtom))
				highPriorityAtom = connAtom;
			}

		int[] oppositeAtom = new int[2];
		int oppositeAtoms = 0;
		for (int i = 0; i< mConnAtoms[oppositeAlleneAtom]; i++) {
			int connAtom = mConnAtom[oppositeAlleneAtom][i];
			if (connAtom != atom)
				oppositeAtom[oppositeAtoms++] = connAtom;
			}

		double alleneAngle = getBondAngle(atom, oppositeAlleneAtom);
		double angleDif = 0.0;

		if (oppositeAtoms == 2) {
			if (oppositeAtom[0] > oppositeAtom[1]) {
				int temp = oppositeAtom[0];
				oppositeAtom[0] = oppositeAtom[1];
				oppositeAtom[1] = temp;
				}

			double hpAngleDif = getAngleDif(alleneAngle, getBondAngle(oppositeAlleneAtom, oppositeAtom[0]));
			double lpAngleDif = getAngleDif(alleneAngle, getBondAngle(oppositeAlleneAtom, oppositeAtom[1]));
			angleDif = hpAngleDif - lpAngleDif;
			}
		else {
			angleDif = getAngleDif(alleneAngle, getBondAngle(oppositeAlleneAtom, oppositeAtom[0]));
			}

		if ((angleDif < 0.0)
		  ^ (getAtomParity(atom) == cAtomParity1)
		  ^ (highPriorityAtom == preferredAtom))
			mBondType[preferredBond] = cBondTypeUp;
		else
			mBondType[preferredBond] = cBondTypeDown;
		}


	/**
	 * In case bond is a BINAP kind of chiral bond with defined parity,
	 * then the preferred neighbour single bond is converted into a
	 * stereo bond to correctly reflect its defined parity.
	 * This method assumes that none of the potential stereo bonds to indicate
	 * the axial configuration is a stereo bond already.
	 * @param bond
	 */
	public void setStereoBondFromBondParity(int bond) {
		// set an optimal bond to up/down to reflect the atom parity
		if (getBondParity(bond) == Molecule.cBondParityNone
		 || getBondParity(bond) == Molecule.cBondParityUnknown
		 || !isBINAPChiralityBond(bond))
			return;

		int preferredBond = -1;
		int preferredAtom = -1;
		int preferredBINAPAtom = -1;
		int oppositeBINAPAtom = -1;
		int bestScore = 0;
		for (int i=0; i<2; i++) {
			int atom = mBondAtom[i][bond];
			for (int j=0; j<mAllConnAtoms[atom]; j++) {
				int connBond = mConnBond[atom][j];
				if (connBond != bond && getBondOrder(connBond) == 1) {
					int connAtom = mConnAtom[atom][j];
					int score = getStereoBondScore(connBond, connAtom);
					if (bestScore < score) {
						bestScore = score;
						preferredAtom = connAtom;
						preferredBond = connBond;
						preferredBINAPAtom = atom;
						oppositeBINAPAtom = mBondAtom[1-i][bond];
   						}
   					}
   				}
   			}

		if (preferredAtom == -1)
			return;

		for (int i=0; i<2; i++) {
			for (int j=0; j<mAllConnAtoms[mBondAtom[i][bond]]; j++) {
				int connBond = mConnBond[mBondAtom[i][bond]][j];
				if (connBond != bond && getBondOrder(connBond) == 1)
					mBondType[connBond] = cBondTypeSingle;
				}
			}

		if (mBondAtom[1][preferredBond] != preferredAtom) {
			mBondAtom[0][preferredBond] = mBondAtom[1][preferredBond];
			mBondAtom[1][preferredBond] = preferredAtom;
			}

		int highPriorityAtom = Integer.MAX_VALUE;
		for (int i = 0; i< mConnAtoms[preferredBINAPAtom]; i++) {
			int connAtom = mConnAtom[preferredBINAPAtom][i];
			if ((mConnBond[preferredBINAPAtom][i] != bond) && (highPriorityAtom > connAtom))
				highPriorityAtom = connAtom;
	   		}

		int[] oppositeAtom = new int[2];
		int oppositeAtoms = 0;
		for (int i = 0; i< mConnAtoms[oppositeBINAPAtom]; i++)
			if (mConnBond[oppositeBINAPAtom][i] != bond)
				oppositeAtom[oppositeAtoms++] = mConnAtom[oppositeBINAPAtom][i];

		double binapAngle = getBondAngle(preferredBINAPAtom, oppositeBINAPAtom);
		double angleDif = 0.0;

		if (oppositeAtoms == 2) {
			if (oppositeAtom[0] > oppositeAtom[1]) {
				int temp = oppositeAtom[0];
				oppositeAtom[0] = oppositeAtom[1];
				oppositeAtom[1] = temp;
				}

			double hpAngleDif = getAngleDif(binapAngle, getBondAngle(oppositeBINAPAtom, oppositeAtom[0]));
			double lpAngleDif = getAngleDif(binapAngle, getBondAngle(oppositeBINAPAtom, oppositeAtom[1]));
			angleDif = hpAngleDif - lpAngleDif;
			}
		else {
			angleDif = getAngleDif(binapAngle, getBondAngle(oppositeBINAPAtom, oppositeAtom[0]));
			}

		if ((angleDif < 0.0)
		  ^ (getBondParity(bond) == cBondParityZor2)
	   	  ^ (highPriorityAtom == preferredAtom))
			mBondType[preferredBond] = cBondTypeUp;
	   	else
	   		mBondType[preferredBond] = cBondTypeDown;
		}


	protected boolean bondsAreParallel(double angle1, double angle2) {
		double angleDif = Math.abs(getAngleDif(angle1, angle2));
		return (angleDif < 0.08 || angleDif > Math.PI - 0.08);
		}


	private int preferredTHStereoBond(int atom, boolean excludeStereoBonds) {
		// If we have two (anti-)parallel bonds, then we need to select
		// that one of those that is closest to the other bonds.
		int allConnAtoms = mAllConnAtoms[atom];
		double[] angle = new double[allConnAtoms];
		for (int i=0; i<allConnAtoms; i++)
			angle[i] = getBondAngle(atom, mConnAtom[atom][i]);
		for (int i=1; i<allConnAtoms; i++) {
			for (int j=0; j<i; j++) {
				if (bondsAreParallel(angle[i], angle[j])) {
					float angleDistanceSum1 = 0;
					float angleDistanceSum2 = 0;
					for (int k=0; k<allConnAtoms; k++) {
						if (k != i && k != j) {
							angleDistanceSum1 += Math.abs(Angle.difference(angle[i], angle[k]));
							angleDistanceSum2 += Math.abs(Angle.difference(angle[j], angle[k]));
							}
						}
					int bond = (angleDistanceSum1 < angleDistanceSum2) ? mConnBond[atom][i] : mConnBond[atom][j];
					if (getBondOrder(bond) == 1
					 && (!excludeStereoBonds || !isStereoBond(bond)))
 						return bond;
					}
				}
			}

		// For every neighbour determine the closest left and the closest right neighbour.
		// If the angle between those two is smaller than 180 degrees, then give a strong preference
		// to the one in the middle to be used as stereo bond.
		boolean[] isPreferred = new boolean[allConnAtoms];
		for (int i=0; i<allConnAtoms; i++) {
			double closestLeftDif = -10.0;
			double closestRightDif = 10.0;
			for (int j=0; j<allConnAtoms; j++) {
				if (j != i) {
					double dif =Angle.difference(angle[i], angle[j]);
					if (dif < 0) {
						if (closestLeftDif < dif)
							closestLeftDif = dif;
						}
					else {
						if (closestRightDif > dif)
							closestRightDif = dif;
						}
					}
				isPreferred[i] = (closestRightDif - closestLeftDif < Math.PI);
				}
			}

		int preferredBond = -1;

		int bestScore = 0;
		for (int i=0; i<allConnAtoms; i++) {
			int connAtom = mConnAtom[atom][i];
			int connBond = mConnBond[atom][i];
			int score = getStereoBondScore(connBond, connAtom);
			if (isPreferred[i])
				score += 16384; // first priority
			if (bestScore < score
			 && (!excludeStereoBonds || !isStereoBond(connBond))) {
				bestScore = score;
				preferredBond = connBond;
				}
			}
		return preferredBond;
		}


	private int preferredAlleneStereoBond(int atom, boolean excludeStereoBonds) {
		int preferredBond = -1;
		int bestScore = 0;
		for (int i=0; i<2; i++) {
			int alleneAtom = mConnAtom[atom][i];
			for (int j=0; j<mAllConnAtoms[alleneAtom]; j++) {
				int connAtom = mConnAtom[alleneAtom][j];
				if (connAtom != atom) {
					int connBond = mConnBond[alleneAtom][j];
					int score = getStereoBondScore(connBond, connAtom);
					if (bestScore < score
					 && (!excludeStereoBonds || !isStereoBond(connBond))) {
						bestScore = score;
						preferredBond = connBond;
						}
					}
				}
			}
		return preferredBond;
		}


	private int preferredBinapStereoBond(int bond) {
		int preferredBond = -1;
		int bestScore = 0;
		for (int i=0; i<2; i++) {
			int atom = mBondAtom[i][bond];
			for (int j=0; j<mAllConnAtoms[atom]; j++) {
				int connAtom = mConnAtom[atom][j];
				if (connAtom != mBondAtom[1-i][bond]) {
					int connBond = mConnBond[atom][j];
					int score = getStereoBondScore(connBond, connAtom);
					if (bestScore < score) {
						bestScore = score;
						preferredBond = connBond;
						}
					}
				}
			}
		return preferredBond;
		}


	/**
	 * If atom is one of the two ends of an allene then returns allene center atom.
	 * @param atom
	 * @return allene center or -1
	 */
	public int findAlleneCenterAtom(int atom) {
		int center = -1;
		if (mPi[atom] == 1) {
			for (int i = 0; i< mConnAtoms[atom]; i++) {
				if (mConnBondOrder[atom][i] == 2) {
					int connAtom = mConnAtom[atom][i];
					if (mConnAtoms[connAtom] == 2 && mPi[connAtom] == 2) {
						for (int j=0; j<2; j++) {
							int endAtom = mConnAtom[connAtom][j];
							if (endAtom != atom && mPi[endAtom] == 1) {
								center = connAtom;
								break;
								}
							}
						}
					break;
					}
				}
			}
		return center;
		}

	/**
	 * Crawls along a chain of sp-hybridized atoms starting from atom2 (which may not be
	 * sp-hybridized) away from its sp-hybridized neighbour atom1. Returns the first atom
	 * that is either not sp-hybridized anymore or the last atom of the chain if that is still
	 * sp-hybridized. Returns -1 in case of an sp-hybridized cycle.
	 * @param atom1 sp-hybridized atom
	 * @param atom2 neighbour atom of atom1
	 * @return first non-sp-hybridized atom when crawling from atom2 away from atom1
	 */
	public int findAlleneEndAtom(int atom1, int atom2) {
		int startAtom = atom1;
		while (mConnAtoms[atom2] == 2 && mPi[atom2] == 2 && atom2 != startAtom) {
			int temp = atom2;
			atom2 = (mConnAtom[atom2][0] == atom1) ? mConnAtom[atom2][1] : mConnAtom[atom2][0];
			atom1 = temp;
			}
		return (atom2 == startAtom) ? -1 : atom2;
		}

		/**
		 * Checks whether atom is one of the two atoms of an axial chirality bond of BINAP type.
		 * Condition: non-aromatic single bond connecting two aromatic rings with 6 or more members
		 * that together bear at least three ortho substituents. A stereo bond indicating the
		 * chirality is not(!!!) a condition.
		 * @param atom to check, whether it is part of a bond, which has BINAP type of axial chirality
		 * @return opposite atom of axial chirality bond or -1 if axial chirality conditions are not met
		 */
	private int findBINAPOppositeAtom(int atom) {
		if (mConnAtoms[atom] == 3 && isAromaticAtom(atom) && getAtomRingSize(atom) >= 6)
			for (int i = 0; i< mConnAtoms[atom]; i++)
				if (isBINAPChiralityBond(mConnBond[atom][i]))
					return mConnAtom[atom][i];
		return -1;
		}


	/**
	 * Checks whether atom is one of the two atoms of an axial chirality bond of BINAP type.
	 * Condition: non-aromatic single bond connecting two aromatic rings with 6 or more members
	 * that together bear at least three ortho substituents. A stereo bond indicating the
	 * chirality is not(!!!) a condition.
	 * @param atom to check, whether it is part of a bond, which has BINAP type of axial chirality
	 * @return axial chirality bond or -1 if axial chirality conditions are not met 
	 */
	public int findBINAPChiralityBond(int atom) {
		if (mConnAtoms[atom] == 3 && isAromaticAtom(atom) && getAtomRingSize(atom) >= 5)
			for (int i = 0; i< mConnAtoms[atom]; i++)
				if (isBINAPChiralityBond(mConnBond[atom][i]))
					return mConnBond[atom][i];
		return -1;
		}


	/**
	 * Evaluates, whether bond is an amide bond, thio-amide, or amidine bond.
	 * @param bond
	 * @return
	 */
	public boolean isAmideTypeBond(int bond) {
		ensureHelperArrays(cHelperNeighbours);

		for (int i=0; i<2; i++) {
			int atom1 = mBondAtom[i][bond];
			if (mAtomicNo[atom1] == 7) {
				int atom2 = mBondAtom[1-i][bond];
				for (int j = 0; j< mConnAtoms[atom2]; j++) {
					int connAtom = mConnAtom[atom2][j];
					int connBond = mConnBond[atom2][j];
					if ((mAtomicNo[connAtom] == 7
					  || mAtomicNo[connAtom] == 8
					  || mAtomicNo[connAtom] == 16)
					 && getBondOrder(connBond) >= 2)
						return true;
					}
				}
			}

		return false;
		}


	/**
	 * @param atom
	 * @return whether this atom is the central atom of an allene
	 */
	public boolean isCentralAlleneAtom(int atom) {
		return mPi[atom] == 2
			&& mConnAtoms[atom] == 2
			&& mConnBondOrder[atom][0] == 2
			&& mConnBondOrder[atom][1] == 2
			&& mAtomicNo[atom] <= 7;
		}


	/**
	 * Checks whether this nitrogen atom is flat, because it has a double bond,
	 * is member of an aromatic ring or is part of amide, an enamine or
	 * in resonance with an aromatic ring. It is also checked that ortho
	 * substituents don't force the amine into a non-resonance torsion.
	 * State of helper arrays must be at least cHelperRings.
	 * @param atom
	 * @return
	 */
	public boolean isFlatNitrogen(int atom) {
		if (mAtomicNo[atom] != 7)
			return false;
		if (isAromaticAtom(atom) || mPi[atom] != 0 || (mAtomQueryFeatures[atom] & cAtomQFFlatNitrogen) != 0)
			return true;
		if (mAtomCharge[atom] == 1)
			return false;
		int heteroCount = 0;
		for (int i = 0; i< mConnAtoms[atom]; i++) {
			if (mConnBondOrder[atom][i] == 1) {
				int atomicNo = mAtomicNo[mConnAtom[atom][i]];
				if (atomicNo == 8 || atomicNo == 9 || atomicNo == 17)
					heteroCount++;
				}
			}
		if (heteroCount == 0) {
			for (int i = 0; i< mConnAtoms[atom]; i++) {
				int connAtom = mConnAtom[atom][i];
				if (mPi[connAtom] != 0) {
					if (isAromaticAtom(connAtom)) {
						if (getAtomRingSize(connAtom) >= 5) {
							int orthoSubstituentCount = 0;
							for (int j=0; j<mConnAtoms[connAtom]; j++) {
								// to be more exact, one might not count ortho-oxygen in case of 5-membered rings
								int ortho = mConnAtom[connAtom][j];
								if (ortho != atom && getNonHydrogenNeighbourCount(ortho) >= 3)
									orthoSubstituentCount++;
								}
							int nitrogenNeighbourCount = getNonHydrogenNeighbourCount(atom);
							if ((orthoSubstituentCount == 2 && nitrogenNeighbourCount >= 2)
							 || (orthoSubstituentCount == 1 && nitrogenNeighbourCount == 3))
								continue;  // the nitrogen is rotated out of PI-plane
							}
						return true;
						}

					// vinyloge amides, etc.
					for (int j = 0; j< mConnAtoms[connAtom]; j++) {
						if ((mConnBondOrder[connAtom][j] == 2 || isAromaticBond(mConnBond[connAtom][j])))
							return true;
						}
					}
				}
			}
		if (heteroCount < 2) {
			for (int i = 0; i< mConnAtoms[atom]; i++) {
				int connAtom = mConnAtom[atom][i];
				boolean isStabilized = false;
				boolean hasCompetitor = false;
				for (int j = 0; j< mConnAtoms[connAtom]; j++) {
					if (mConnAtom[connAtom][j] != atom) {
						if (mConnBondOrder[connAtom][j] != 1
						 && (mAtomicNo[mConnAtom[connAtom][j]] == 7
						  || mAtomicNo[mConnAtom[connAtom][j]] == 8
						  || mAtomicNo[mConnAtom[connAtom][j]] == 16))
							isStabilized = true;
	
						if (mConnBondOrder[connAtom][j] == 1
				   		 && mAtomicNo[mConnAtom[connAtom][j]] == 7)
							hasCompetitor = true;
						}
					}
				if (isStabilized && (!hasCompetitor || heteroCount == 0))
					return true;
				}
			}
		return false;
		}

	
	/**
	 * Checks whether bond is an axial chirality bond of the BINAP type.
	 * Condition: non-aromatic, non-small-ring (<= 7 members) single bond
	 * connecting two aromatic rings with 6 or more members each
	 * that together bear at least three ortho substituents. A stereo bond indicating the
	 * chirality is not(!!!) a condition.
	 * @param bond
	 * @return true if axial chirality conditions are met
	 */
	public boolean isBINAPChiralityBond(int bond) {
		if (mBondType[bond] != cBondTypeSingle
		 || isAromaticBond(bond)
		 || (isRingBond(bond) && getBondRingSize(bond) < 7))
			return false;

		int atom1 = mBondAtom[0][bond];
		if (!isAromaticAtom(atom1)
		 || getAtomRingSize(atom1) < 5)
			return false;

		int atom2 = mBondAtom[1][bond];
		if (!isAromaticAtom(atom2)
		 || getAtomRingSize(atom2) < 5)
			return false;

		int orthoSubstituentCount1 = getOrthoSubstituentCount(atom1, atom2);
		int orthoSubstituentCount2 = getOrthoSubstituentCount(atom2, atom1);

		// with aromatic 6-membered (or larger) rings only, we have a simple logic
		if (getAtomRingSize(atom1) > 5 && getAtomRingSize(atom2) > 5)
			return orthoSubstituentCount1 + orthoSubstituentCount2 > 2;

		int secondOrderOrthoSubstituentCount1 = getSecondOrderOrthoSubstituentCount(atom1, atom2);
		int secondOrderOrthoSubstituentCount2 = getSecondOrderOrthoSubstituentCount(atom2, atom1);

		// for 5-membered rings we need stronger constraints (currently we don't distignuid 5-5 and 5-6)

		if (orthoSubstituentCount1 == 2 && secondOrderOrthoSubstituentCount2 >= 1)
			return true;
		if (orthoSubstituentCount2 == 2 && secondOrderOrthoSubstituentCount1 >= 1)
			return true;

		if (secondOrderOrthoSubstituentCount1 == 2 && (orthoSubstituentCount2 >= 1 || secondOrderOrthoSubstituentCount2 >= 1))
			return true;
		if (secondOrderOrthoSubstituentCount2 == 2 && (orthoSubstituentCount1 >= 1 || secondOrderOrthoSubstituentCount1 >= 1))
			return true;

		return false;
		}

	private int getOrthoSubstituentCount(int atom, int otherBondAtom) {
		int count = 0;
		for (int i=0; i<mConnAtoms[atom]; i++) {
			int connAtom = mConnAtom[atom][i];
			if (connAtom != otherBondAtom && mConnAtoms[connAtom] > 2)
				count++;
			}
		return count;
		}

	private int getSecondOrderOrthoSubstituentCount(int atom, int otherBondAtom) {
		int count = 0;
		for (int i=0; i<mConnAtoms[atom]; i++) {
			int connAtom = mConnAtom[atom][i];
			if (connAtom != otherBondAtom) {
				int innerCount = 0;
				for (int j=0; j<mConnAtoms[connAtom]; j++) {
					int nextConnAtom = mConnAtom[connAtom][j];
					if (nextConnAtom != atom
					 && isAromaticBond(mConnBond[connAtom][j])
					 && mConnAtoms[nextConnAtom] > 2)
						innerCount++;
					}
				// if we have two next neighbours with 3 neighbours each, then one of them must be a second order ortho
				if (innerCount == 2)
					count++;
				}
			}
		return count;
		}

	protected boolean validateBondType(int bond, int type) {
		boolean ok = super.validateBondType(bond, type);

		if (ok && type == cBondTypeCross) {
			ensureHelperArrays(Molecule.cHelperRings);
			ok &= !isSmallRingBond(bond);
			}

		return ok;
		}

	
	public void validate() throws Exception {
		double avbl = getAverageBondLength();
		double minDistanceSquare = avbl * avbl / 16.0;
		for (int atom1=1; atom1<mAllAtoms; atom1++) {
			for (int atom2=0; atom2<atom1; atom2++) {
				double xdif = mCoordinates[atom2].x - mCoordinates[atom1].x;
				double ydif = mCoordinates[atom2].y - mCoordinates[atom1].y;
				double zdif = mCoordinates[atom2].z - mCoordinates[atom1].z;
				if ((xdif*xdif + ydif*ydif + zdif*zdif) < minDistanceSquare)
					throw new Exception("The distance between two atoms is too close.");
				}
			}

		ensureHelperArrays(cHelperNeighbours);
		int allCharge = 0;
		for (int atom=0; atom<mAtoms; atom++) {
			if (getOccupiedValence(atom) > getMaxValence(atom))
				throw new Exception("atom valence exceeded");
			allCharge += mAtomCharge[atom];
			}
		if (allCharge != 0)
			throw new Exception("unbalanced atom charge");
		}

	/**
	 * Normalizes different forms of functional groups (e.g. nitro)
	 * to a preferred one. This step should precede any canonicalization.
	 * @return true if the molecule was changed
	 */
	public boolean normalizeAmbiguousBonds() {
		ensureHelperArrays(cHelperNeighbours);

		// Molecules from Marvin Sketch, if generated from SMILES, may contain delocalized bonds
		// using the Daylight aromaticity model, e.g. with aromatic carbonyl carbons in case of pyridinone.
		// Before generating idcodes, we need to convert to cumulated double bonds and use our own
		// aromaticity model to locate delocalized bonds.
		normalizeExplicitlyDelocalizedBonds();

		boolean found = false;
		for (int atom=0; atom<mAtoms; atom++) {
			if (mAtomicNo[atom] == 7
			 && mAtomCharge[atom] == 0) {
				int valence = getOccupiedValence(atom);
				if (valence == 4) {
					// normalize -N(=O)-OH to -N(+)(=O)-O(-)
					for (int i=0; i<mConnAtoms[atom]; i++) {
						int connAtom = mConnAtom[atom][i];
						if (mConnBondOrder[atom][i] == 1
						 && mAtomicNo[connAtom] == 8
						 && mConnAtoms[connAtom] == 1
						 && mAtomCharge[connAtom] == 0) {
							found = true;
							mAtomCharge[atom]++;
							mAtomCharge[connAtom]--;
							break;
							}
						}
					}
				else if (valence == 5) {
					// normalize -N(=O)2 to -N(+)(=O)-O(-), -N#N to -N(+)=N(-),
					for (int i=0; i<mConnAtoms[atom]; i++) {
						int connAtom = mConnAtom[atom][i];
						int connBond = mConnBond[atom][i];
						if (mConnBondOrder[atom][i] == 2
						 && mAtomicNo[connAtom] == 8) {
							found = true;
							mAtomCharge[atom]++;
							mAtomCharge[connAtom]--;
							mBondType[connBond] = Molecule.cBondTypeSingle;
							break;
							}
						if (mConnBondOrder[atom][i] == 3
						 && mAtomicNo[connAtom] == 7) {
							found = true;
							mAtomCharge[atom]++;
							mAtomCharge[connAtom]--;
							mBondType[connBond] = Molecule.cBondTypeDouble;
							break;
							}
						}
					}
				}
			}

		// split covalent bonds between hetero atoms and one of (Li,Na,K,Mg,Ca,...)
		boolean bondDeleted = false;
		for (int bond=0; bond<mBonds; bond++) {
			for (int i=0; i<2; i++) {
				if (isElectronegative(mBondAtom[i][bond])) {
					int atom = mBondAtom[1-i][bond];
					if (isAlkaliMetal(atom) || isEarthAlkaliMetal(atom)) {
						if (getBondOrder(bond) == 1) {
							mAtomCharge[atom]++;
							mAtomCharge[mBondAtom[i][bond]]--;
							mBondType[bond] = cBondTypeDeleted;
							bondDeleted = true;
							}
						else if (mBondType[bond] == cBondTypeMetalLigand) {
							mBondType[bond] = cBondTypeDeleted;
							bondDeleted = true;
							}
						}
					break;
					}
				}
			}
		if (bondDeleted) {
			compressMolTable();
			found = true;
			}

		if (found)
			mValidHelperArrays = cHelperNone;

		return found;
		}

	/**
	 * Checks, converts any explicitly delocalized bonds with cBondTypeDelocalized into alternating
	 * single and double bonds.
	 * @return
	 */
	private boolean normalizeExplicitlyDelocalizedBonds() {
		for (int bond=0; bond<mBonds; bond++)
			if (mBondType[bond] == cBondTypeDelocalized)
				return new AromaticityResolver(this).locateDelocalizedDoubleBonds(null);

		return false;
		}

	/**
	 * @param atom
	 * @return whether atom is one of Li,Na,K,Rb,Cs
	 */
	public boolean isAlkaliMetal(int atom) {
		int atomicNo = mAtomicNo[atom];
		return atomicNo == 3	 // Li
			|| atomicNo == 11	// Na
			|| atomicNo == 19	// K
			|| atomicNo == 37	// Rb
			|| atomicNo == 55;	// Cs
		}

	/**
	 * @param atom
	 * @return whether atom is one of Mg,Ca,Sr,Ba
	 */
	public boolean isEarthAlkaliMetal(int atom) {
		int atomicNo = mAtomicNo[atom];
		return atomicNo == 12	// Mg
			|| atomicNo == 20	// Ca
			|| atomicNo == 38	// Sr
			|| atomicNo == 56;	// Ba
		}

	/**
	 * @param atom
	 * @return whether atom is one of N,P,As
	 */
	public boolean isNitrogenFamily(int atom) {
		int atomicNo = mAtomicNo[atom];
		return atomicNo == 7	// N
			|| atomicNo == 15	// P
			|| atomicNo == 33;	// As
		}

	/**
	 * @param atom
	 * @return whether atom is one of O,S,Se,Te
	 */
	public boolean isChalcogene(int atom) {
		int atomicNo = mAtomicNo[atom];
		return atomicNo == 8	// O
			|| atomicNo == 16	// S
			|| atomicNo == 34	// Se
			|| atomicNo == 52;	// Te
		}

	/**
	 * @param atom
	 * @return whether atom is one of F,Cl,Br,I
	 */
	public boolean isHalogene(int atom) {
		int atomicNo = mAtomicNo[atom];
		return atomicNo == 9	// F
			|| atomicNo == 17	// Cl
			|| atomicNo == 35	// Br
			|| atomicNo == 53;	// I
		}

	/**
	 * Normalizes charge distribution in single- and multifragment molecules.
	 * In a first step polar bonds (both atoms have opposite charge) are neutralized
	 * by removing both atom charges and increasing the bond order, provided that atom
	 * valences allow the change.
	 * Neutralizes positive and an equal amount of negative charges on electronegative atoms,
	 * provided these are not on 1,2-dipolar structures, in order to ideally achieve a neutral molecule.
	 * This method does not change the overall charge of the molecule. It does not change the number of
	 * explicit atoms or bonds or their connectivity except bond orders.
	 * This method does not deprotonate acidic groups to compensate for quarternary charged nitrogen.
	 * @param allowUnbalancedCharge throws an exception after polar bond neutralization, if overall charge is not zero
	 * @return remaining overall molecule charge
	 */
	public int canonizeCharge(boolean allowUnbalancedCharge) throws Exception {
		return canonizeCharge(allowUnbalancedCharge, false);
		}

	/**
	 * Normalizes charge distribution in single- and multifragment molecules.
	 * In a first step polar bonds (both atoms have opposite charge) are neutralized
	 * by removing both atom charges and increasing the bond order, provided that atom
	 * valences allow the change.
	 * Neutralizes positive and an equal amount of negative charges on electronegative atoms,
	 * provided these are not on 1,2-dipolar structures, in order to ideally achieve a neutral molecule.
	 * This method may change the overall charge of the molecule if doNeutralize==true.
	 * It does not change the number of explicit atoms or bonds or their connectivity except bond orders.
	 * This method does not deprotonate acidic groups to compensate for quarternary charged nitrogen.
	 * @param allowUnbalancedCharge throws an exception after polar bond neutralization, if overall charge is not zero
	 * @param doNeutralize if true, then tries to neutralize the molecule if overall charges are not zero
	 * @return remaining overall molecule charge
	 */
	public int canonizeCharge(boolean allowUnbalancedCharge, boolean doNeutralize) throws Exception {
		ensureHelperArrays(cHelperNeighbours);

		if (doNeutralize)
			allowUnbalancedCharge = true;

		for (int bond=0; bond<mAllBonds; bond++) {
			int bondOrder = getBondOrder(bond);
			if (bondOrder == 1 || bondOrder == 2) {
				int atom1,atom2;
				if (mAtomCharge[mBondAtom[0][bond]] > 0
				 && mAtomCharge[mBondAtom[1][bond]] < 0) {
					atom1 = mBondAtom[0][bond];
					atom2 = mBondAtom[1][bond];
					}
				else if (mAtomCharge[mBondAtom[0][bond]] < 0
					  && mAtomCharge[mBondAtom[1][bond]] > 0) {
					atom1 = mBondAtom[1][bond];
					atom2 = mBondAtom[0][bond];
					}
				else
					continue;

				if (isMetalAtom(atom1)
				 || isMetalAtom(atom2))
					continue;

				if ((mAtomicNo[atom1] < 9 && getOccupiedValence(atom1) > 3)
				 || (mAtomicNo[atom2] < 9 && getOccupiedValence(atom2) > 3))
						continue;

				boolean hasImplicitHydrogen = (getImplicitHydrogens(atom1) != 0);

				mAtomCharge[atom1] -= 1;
				mAtomCharge[atom2] += 1;

				// If the formerly positive atom had an implicit hydrogen, then removing
				// both charges mean shifting the implicit hydrogen to atom2.
				// If it had no implicit hydrogen, but an open valence, then we shift
				// the electron pair from the negative atom into the bond and increase the
				// bond order.
				if (!hasImplicitHydrogen) {
					int oldBondType = mBondType[bond];
					if (bondOrder == 1)
						mBondType[bond] = cBondTypeDouble;
					else
						mBondType[bond] = cBondTypeTriple;

					if (oldBondType == cBondTypeDown
					 || oldBondType == cBondTypeUp) {   // retain stereo information
						int stereoCenter = mBondAtom[0][bond];
						int newStereoBond = preferredTHStereoBond(stereoCenter, false);
						if (mBondAtom[0][newStereoBond] != stereoCenter) {
							mBondAtom[1][newStereoBond] = mBondAtom[0][newStereoBond];
							mBondAtom[1][newStereoBond] = stereoCenter;
							}
						}
					}

				mValidHelperArrays = cHelperNone;
				}
			}

		int overallCharge = 0;
		int negativeAtomCount = 0;
		int negativeAdjustableCharge = 0;
		for (int atom=0; atom<mAllAtoms; atom++) {
			overallCharge += mAtomCharge[atom];
			if (mAtomCharge[atom] < 0 && !hasPositiveNeighbour(atom)) {
				negativeAtomCount++;
				if (isElectronegative(atom))
					negativeAdjustableCharge -= mAtomCharge[atom];
				}
			}

		if (!allowUnbalancedCharge && overallCharge != 0)
			throw new Exception("molecule's overall charges are not balanced");

		ensureHelperArrays(cHelperNeighbours);
		int overallChargeChange = 0;
		int positiveChargeForRemoval = doNeutralize ?
					overallCharge + negativeAdjustableCharge : negativeAdjustableCharge;
		for (int atom=0; atom<mAllAtoms; atom++) {
			if (mAtomCharge[atom] > 0) {
				if (!hasNegativeNeighbour(atom) && isElectronegative(atom)) {
					int chargeReduction = Math.min(getImplicitHydrogens(atom), mAtomCharge[atom]);
					if (chargeReduction != 0 && positiveChargeForRemoval >= chargeReduction) {
						overallCharge -= chargeReduction;
						overallChargeChange -= chargeReduction;
						positiveChargeForRemoval -= chargeReduction;
						mAtomCharge[atom] -= chargeReduction;
						mValidHelperArrays &= cHelperNeighbours;
						}
					}
				}
			}

		int negativeChargeForRemoval = doNeutralize ? overallCharge : overallChargeChange;
		if (negativeChargeForRemoval < 0) {
			int[] negativeAtom = new int[negativeAtomCount];
			negativeAtomCount = 0;
			for (int atom=0; atom<mAllAtoms; atom++) {
				if (mAtomCharge[atom] < 0 && !hasPositiveNeighbour(atom)) {
					// Ideally priorities for negative atom protonation should
					// be done based on atom priorities from Canonizer.
					negativeAtom[negativeAtomCount++] = (mAtomicNo[atom] << 22)
													  + atom;
					}
				}
			java.util.Arrays.sort(negativeAtom);
			for (int i=negativeAtom.length-1; (negativeChargeForRemoval < 0) && (i>=negativeAtom.length-negativeAtomCount); i--) {
				int atom = negativeAtom[i] & 0x003FFFFF;
				if (isElectronegative(atom)) {
					int chargeReduction = Math.min(-negativeChargeForRemoval, -mAtomCharge[atom]);
					overallCharge += chargeReduction;
					negativeChargeForRemoval += chargeReduction;
					mAtomCharge[atom] += chargeReduction;
					mValidHelperArrays &= cHelperNeighbours;
					}
				}
			}

		return overallCharge;
		}


	/**
	 * @param atom
	 * @return true if atom has a neighbour with a negative charge
	 */
	private boolean hasNegativeNeighbour(int atom) {
		for (int i = 0; i< mConnAtoms[atom]; i++)
			if (mAtomCharge[mConnAtom[atom][i]] < 0)
				return true;
		return false;
		}


	/**
	 * @param atom
	 * @return true if atom has a neighbour with a positive charge
	 */
	private boolean hasPositiveNeighbour(int atom) {
		for (int i = 0; i< mConnAtoms[atom]; i++)
			if (mAtomCharge[mConnAtom[atom][i]] > 0)
				return true;
		return false;
		}


	/**
	 * Provided that the bond parity of a double bond is available,
	 * this method determines, whether connAtom has a counterpart with
	 * Z- (cis) configuration at the other end of the double bond.
	 * If there is no Z-counterpart, then -1 is returned.
	 * Requires cHelperParities.
	 * @param connAtom directly connected to one of the double bond atoms
	 * @param bond double bond with available bond parity
	 * @return -1 or counterpart to connAtom in Z-configuration
	 */
	public int getZNeighbour(int connAtom, int bond) {
		if (getBondOrder(bond) != 2 && !isAromaticBond(bond))
			return -1;
		int parity = getBondParity(bond);
		if (parity != cBondParityEor1 && parity != cBondCIPParityZorM)
			return -1;

		for (int i=0; i<2; i++) {
			int atom1 = mBondAtom[i][bond];
			int atom2 = mBondAtom[1-i][bond];
			int other1 = -1;
			boolean found = false;
			for (int j = 0; j< mConnAtoms[atom1]; j++) {
				int conn = mConnAtom[atom1][j];
				if (conn != atom2) {
					if (conn == connAtom)
						found = true;
					else
						other1 = conn;
					}
				}
			if (found) {
				int lowConn = -1;
				int highConn = -1;
				for (int j = 0; j< mConnAtoms[atom2]; j++) {
					int conn = mConnAtom[atom2][j];
					if (conn != atom1) {
						if (lowConn == -1)
							lowConn = conn;
						else if (conn > lowConn)
							highConn = conn;
						else {
							highConn = lowConn;
							lowConn = conn;
							}
						}
					}

				if (mConnAtoms[atom1] == 2) {
					if (mConnAtoms[atom2] == 2)
						return parity == cBondCIPParityZorM ? lowConn : -1;
					return (parity == cBondCIPParityZorM) ? lowConn : highConn;
					}
				else {
					if (mConnAtoms[atom2] == 2)
						return (parity == cBondCIPParityZorM) ^ (connAtom < other1) ? -1 : lowConn;
					return (parity == cBondCIPParityZorM) ^ (connAtom < other1) ? highConn : lowConn;
					}
				}
			}
		return -1;
		}

	public int getHelperArrayStatus() {
		return mValidHelperArrays;
		}

	/**
	 * While the Molecule class covers all primary molecule information, its derived class
	 * ExtendedMolecule handles secondary, i.e. calculated molecule information, which is cached
	 * in helper arrays and stays valid as long as the molecule's primary information is not changed.
	 * Most methods of ExtendedMolecule require some of the helper array's information. High level
	 * methods, e.g. getPath(), take care of updating an outdated cache themselves. Low level methods,
	 * e.g. isAromaticAtom(), which typically are called very often, do not check for validity
	 * nor update the helper arrays themselves. If you use low level methods, then you need to make
	 * sure that the needed helper array information is valid by this method.<br>
	 * For performance reasons there are <b>distinct levels of helper information</b>. (A higher
	 * level always includes all properties of the previous level):<br>
	 * <i>cHelperNeighbours:</i> explicit hydrogen atoms are moved to the end of the atom table and
	 * bonds leading to them are moved to the end of the bond table. This way algorithms can skip
	 * hydrogen atoms easily. For every atom directly connected atoms and bonds (with and without
	 * hydrogens) are determined. The number of pi electrons is counted.<br>
	 * <i>cHelperRings</i>: Aromatic and non-aromatic rings are detected. Atom and bond ring
	 * properties are set and a ring collection provides a total set of small rings (7 or less atoms).
	 * Atoms being in allylic/benzylic or stabilized (neighbor of a carbonyl or similar group) position
	 * are flagged as such.<br>
	 * <i>cHelperParities</i>: Atom (tetrahedral or axial) and bond (E/Z or atrop) parities are calculated
	 * from the stereo configurations.<br>
	 * <i>cHelperCIP</i>: Cahn-Ingold-Prelog stereo information for atoms and bonds.<br>
	 * <br>cHelperParities and cHelperCIP require a StereoMolecule!!!<br>
	 * @param required one of cHelperNeighbours,cHelperRings,cHelperParities,cHelperCIP
	 * @return true if the molecule was changed
	 */
	public void ensureHelperArrays(int required) {
		// cHelperNeighbours: mConnAtoms,mConnBonds,mPi for all atoms
		// cHelperRings: rings,aromaticity/allylic/stabilized for non-H-atoms only
		// cHelperParities: stereo parities for non-H-atoms/bonds only
		// cHelperCIP: mCanonizer, stereo parities for non-H-atoms/bonds only

		if ((required & ~mValidHelperArrays) == 0)
			return;

		if ((mValidHelperArrays & cHelperBitNeighbours) == 0) {
			handleHydrogens();
			calculateNeighbours();

			mValidHelperArrays |= cHelperBitNeighbours;

			if (convertHydrogenToQueryFeatures()) {
				handleHydrogens();
				calculateNeighbours();
				}
			}

		if ((required & ~mValidHelperArrays) == 0)
			return;

		if ((mValidHelperArrays & ~(cHelperBitRingsSimple | cHelperBitRings)) != 0) {
			for (int atom = 0; atom < mAtoms; atom++)
				mAtomFlags[atom] &= ~cAtomFlagsHelper2;
			for (int bond = 0; bond < mBonds; bond++)
				mBondFlags[bond] &= ~cBondFlagsHelper2;

			// if we are asked to only detect small rings and skip aromaticity, allylic and stabilized detection
			if ((required & cHelperBitRings) == 0) {
				findRings(RingCollection.MODE_SMALL_RINGS_ONLY);
				mValidHelperArrays |= cHelperBitRingsSimple;
				return;
				}

			findRings(RingCollection.MODE_SMALL_AND_LARGE_RINGS_AND_AROMATICITY);

			for (int atom=0; atom<mAtoms; atom++) {	// allylic & stabilized flags
				for (int i = 0; i< mConnAtoms[atom]; i++) {
					int connBond = mConnBond[atom][i];
					if (isAromaticBond(connBond))
						continue;

					int connAtom = mConnAtom[atom][i];
					for (int j = 0; j< mConnAtoms[connAtom]; j++) {
						if (mConnBond[connAtom][j] == connBond)
							continue;

						if (mConnBondOrder[connAtom][j] > 1) {
							if (mAtomicNo[mConnAtom[connAtom][j]] == 6)
								mAtomFlags[atom] |= cAtomFlagAllylic;
							else {
								if (!isAromaticBond(mConnBond[connAtom][j])
								 && isElectronegative(mConnAtom[connAtom][j]))
									mAtomFlags[atom] |= cAtomFlagStabilized;
								}
							}
						}
					}
				}

				// propagate stabilized flags to vinylic positions
			while (true) {
				boolean found = false;
				for (int atom=0; atom<mAtoms; atom++) {
						 // for non-aromatic stabilized atoms with pi-electrons
					if (mPi[atom] > 0
					 && (mAtomFlags[atom] & cAtomFlagStabilized) != 0
					 && !mRingSet.isAromaticAtom(atom)) {
						for (int i = 0; i< mConnAtoms[atom]; i++) {
							if (mConnBondOrder[atom][i] > 1) {
								int connAtom = mConnAtom[atom][i];
								int connBond = mConnBond[atom][i];
								for (int j = 0; j< mConnAtoms[connAtom]; j++) {
									if (mConnBond[connAtom][j] != connBond) {
										int candidate = mConnAtom[connAtom][j];
										if ((mAtomFlags[candidate] & cAtomFlagStabilized) == 0) {
											mAtomFlags[candidate] |= cAtomFlagStabilized;
											found = true;
											}
										}
									}
								}
							}
						}
					}
				if (!found)
					break;
				}
			}

		mValidHelperArrays |= (cHelperBitRingsSimple | cHelperBitRings);
		}


	/**
	 * Usually explicit hydrogen atoms can be removed without changing a molecule,
	 * because the removal just converts explicit into implicit hydrogen atoms.
	 * Exceptions are hydrogen with isotop information, hydrogens not connected to a non-H atom,
	 * hydrogens carrying a custom label, hydrogens whose existence implicitly defines a neighbour
	 * atom to have an abnormal valence, hydrogens with a different bonding environment than exactly
	 * one single bond, and hydrogen atoms connected to metal atoms.<br>
	 * This method moves all simple hydrogen atoms and associated bonds to the end of the atom/bond tables.
	 * It sets mAtoms to exclude simple hydrogen atoms and mBonds to exclude bonds leading to them.
	 * Simple hydrogens are not deleted, though. They are always displayed and the stereo perception
	 * still considers up/down bonds leading to hydrogen atoms. Most other functions, however, can
	 * happily neglect them.<br>
	 * If non-simple hydrogen atoms (D,T and custom labelled hydrogen) exist, then these are positioned
	 * after the non-hydrogen atoms and just before the simple hydrogen atoms. mAtoms and mBonds include
	 * these atoms.<br>
	 * mConnAtoms/mConnBonds/mConnBondOrder are neither used nor updated.<br>
	 * <b>Note:</b> This method changes the order among the non-hydrogen atoms. To translate to the original
	 * order use getHandleHydrogenMap() before calling ensureHelperArrays() if the original atom order is relevant.
	 */
	private void handleHydrogens() {
		// find all hydrogens that are connected to a non-H atom and therefore can be implicit		
		boolean[] isHydrogen = findSimpleHydrogens();

					// move all hydrogen atoms to end of atom table
		int lastNonHAtom = mAllAtoms;
		do lastNonHAtom--;
			while ((lastNonHAtom >= 0) && isHydrogen[lastNonHAtom]);

		for (int atom=0; atom<lastNonHAtom; atom++) {
			if (isHydrogen[atom]) {
				swapAtoms(atom, lastNonHAtom);

				// swap simple H flags also
				boolean temp = isHydrogen[atom];
				isHydrogen[atom] = isHydrogen[lastNonHAtom];
				isHydrogen[lastNonHAtom] = temp;

				do lastNonHAtom--;
					while (isHydrogen[lastNonHAtom]);
				}
			}
		mAtoms = lastNonHAtom + 1;

					// move all bonds to hydrogen to end of bond table
		if (mAllAtoms == mAtoms) {
			mBonds = mAllBonds;
			return;
			}

		boolean isHydrogenBond[] = new boolean[mAllBonds];
		for (int bond=0; bond<mAllBonds; bond++) {	// mark all bonds to hydrogen
			int atom1 = mBondAtom[0][bond];
			int atom2 = mBondAtom[1][bond];
			if (isHydrogen[atom1]
			 || isHydrogen[atom2])
				isHydrogenBond[bond] = true;
			}

		int lastNonHBond = mAllBonds;
		do lastNonHBond--; while ((lastNonHBond >= 0) && isHydrogenBond[lastNonHBond]);

		for (int bond=0; bond<lastNonHBond; bond++) {
			if (isHydrogenBond[bond]) {
				swapBonds(bond, lastNonHBond);
/*				int tempInt = mBondAtom[0][bond];
				mBondAtom[0][bond] = mBondAtom[0][lastNonHBond];
				mBondAtom[0][lastNonHBond] = tempInt;
				tempInt = mBondAtom[1][bond];
				mBondAtom[1][bond] = mBondAtom[1][lastNonHBond];
				mBondAtom[1][lastNonHBond] = tempInt;
				tempInt = mBondType[bond];
				mBondType[bond] = mBondType[lastNonHBond];
				mBondType[lastNonHBond] = tempInt;  */
				isHydrogenBond[bond] = false;
				do lastNonHBond--;
					while (isHydrogenBond[lastNonHBond]);
				}
			}
		mBonds = lastNonHBond + 1;
		}

	/**
	 * If ensureHelperArrays() (and with it handleHydrogens()) was not called yet
	 * on a fresh molecule and if the molecule contains simple hydrogen atoms within
	 * non-hydrogens atoms, then this function returns a map from current atom indexes
	 * to those new atom indexes that would result from a call to handleHydrogens.
	 * @return
	 */
	public int[] getHandleHydrogenMap() {
		return getHandleHydrogenAtomMap(findSimpleHydrogens());
		}

	/**
	 * If ensureHelperArrays() (and with it handleHydrogens()) was not called yet
	 * on a fresh molecule and if the molecule contains simple hydrogen atoms within
	 * non-hydrogens atoms, then this function returns a map from current atom indexes
	 * to those new atom indexes that would result from a call to handleHydrogens.
	 * @return
	 */
	public int[] getHandleHydrogenAtomMap(boolean[] isSimpleHydrogen) {
		int[] map = new int[mAllAtoms];
		for (int i=0; i<mAllAtoms; i++)
			map[i] = i;

		int lastNonHAtom = mAllAtoms;
		do lastNonHAtom--;
		while ((lastNonHAtom >= 0) && isSimpleHydrogen[lastNonHAtom]);

		for (int i=0; i<lastNonHAtom; i++) {
			if (isSimpleHydrogen[map[i]]) {
				int tempIndex = map[i];
				map[i] = map[lastNonHAtom];
				map[lastNonHAtom] = tempIndex;

				// swap simple H flags also
				boolean temp = isSimpleHydrogen[i];
				isSimpleHydrogen[i] = isSimpleHydrogen[lastNonHAtom];
				isSimpleHydrogen[lastNonHAtom] = temp;

				do lastNonHAtom--;
				while (isSimpleHydrogen[lastNonHAtom]);
				}
			}

		return map;
		}

	/**
	 * If ensureHelperArrays() (and with it handleHydrogens()) was not called yet
	 * on a fresh molecule and if the molecule contains simple hydrogen atoms within
	 * non-hydrogens atoms, then this function returns a map from current bond indexes
	 * to those new bond indexes that would result from a call to handleHydrogens.
	 * @return
	 */
	public int[] getHandleHydrogenBondMap() {
		boolean[] isSimpleHydrogen = findSimpleHydrogens();
		int[] map = new int[mAllBonds];
		for (int i=0; i<mAllBonds; i++)
			map[i] = i;
		
		boolean isHydrogenBond[] = new boolean[mAllBonds];
		for (int bond=0; bond<mAllBonds; bond++)	// mark all bonds to hydrogen
			isHydrogenBond[bond] = isSimpleHydrogen[mBondAtom[0][bond]]
								|| isSimpleHydrogen[mBondAtom[1][bond]];

		int lastNonHBond = mAllBonds;
		do lastNonHBond--; while ((lastNonHBond >= 0) && isHydrogenBond[lastNonHBond]);

		for (int i=0; i<lastNonHBond; i++) {
			if (isHydrogenBond[map[i]]) {
				int tempIndex = map[i];
				map[i] = map[lastNonHBond];
				map[lastNonHBond] = tempIndex;

				do lastNonHBond--;
					while (isHydrogenBond[map[lastNonHBond]]);
				}
			}

		return map;
		}

	private boolean[] findSimpleHydrogens() {
		boolean[] isSimpleHydrogen = new boolean[mAllAtoms];
		for (int atom=0; atom<mAllAtoms; atom++)
			isSimpleHydrogen[atom] = isSimpleHydrogen(atom);

		// unflag simple hydrogens that have a bond with order != 1
		// or that have more than one bond or that are connected to metal atoms
		boolean[] oneBondFound = new boolean[mAllAtoms];
		for (int bond=0; bond<mAllBonds; bond++) {
			int atom1 = mBondAtom[0][bond];
			int atom2 = mBondAtom[1][bond];

			if (getBondOrder(bond) != 1) {
				isSimpleHydrogen[atom1] = false;
				isSimpleHydrogen[atom2] = false;
				continue;
				}

			if (oneBondFound[atom1])
				isSimpleHydrogen[atom1] = false;
			if (oneBondFound[atom2])
				isSimpleHydrogen[atom2] = false;

			if (isSimpleHydrogen[atom1] && isMetalAtom(atom2) && mAtomicNo[atom2] != 13)    // Al
				isSimpleHydrogen[atom1] = false;
			if (isSimpleHydrogen[atom2] && isMetalAtom(atom1) && mAtomicNo[atom1] != 13)
				isSimpleHydrogen[atom2] = false;

			oneBondFound[atom1] = true;
			oneBondFound[atom2] = true;
			}

		// unflag simple hydrogens within an H2 molecule
		for (int bond=0; bond<mAllBonds; bond++) {
			if (isSimpleHydrogen[mBondAtom[0][bond]] && isSimpleHydrogen[mBondAtom[1][bond]]) {
				isSimpleHydrogen[mBondAtom[0][bond]] = false;
				isSimpleHydrogen[mBondAtom[1][bond]] = false;
				}
			}

			// unflag simple hydrogens that have no connection to another atom
		for (int atom=0; atom<mAllAtoms; atom++)
			if (!oneBondFound[atom])
				isSimpleHydrogen[atom] = false;

		return isSimpleHydrogen;
		}

	/**
	 * Uncharged hydrogen atoms with no isotop information nor with an attached custom label
	 * are considered simple and can usually be suppressed, effectively converting them from an
	 * explicit to an implicit hydrogen atom.<br>
	 * <b>Note:</b> This method returns true for uncharged, natural abundance hydrogens without
	 * custom labels even if they have a non-standard bonding situation (everything being different
	 * from having one single bonded non-simple-hydrogen neighbour, e.g. unbonded hydrogen, H2,
	 * a metal ligand bond to another atom, two single bonds, etc.)
	 * If unusual bonding needs to be considered, check for that independently from this method.
	 * @param atom
	 * @return
	 */
	public boolean isSimpleHydrogen(int atom) {
		return mAtomicNo[atom] == 1 && mAtomMass[atom] == 0 && mAtomCharge[atom] == 0 // && mAtomMapNo[atom] == 0
			&& (mAtomCustomLabel == null || mAtomCustomLabel[atom] == null);
		// Since a mapNo is not part of an idcode, a mapped but otherwise simple H must be considered simple; TLS 29Aug2020
		}

	/**
	 * Removes all plain explicit hydrogens atoms from the molecule, converting them effectively to implicit ones.
	 * Assuming that this molecule has 2D-coordinates, this method perceives stereo configurations from up/down-bonds
	 * to explicit hydrogens before deleting them and turns another bond into a stereo bond to indicate the proper configuration.
	 * If the removal of a hydrogen atom would change an atom's implicit valence, the atom's abnormal valence is set accordingly.
	 */
	public void removeExplicitHydrogens() {
		removeExplicitHydrogens(true, false);
		}

	/**
	 * Removes all plain explicit hydrogens atoms from the molecule, converting them effectively to implicit ones.
	 * If the molecules has 2D-coordinates (hasValid2DCoords==true), then this method initially perceives stereo
	 * configurations E/Z/R/S from coordinates and up/down-bonds to explicit hydrogens before deleting them
	 * and turns another bond into a stereo bond to indicate the proper configuration.
	 * If the removal of a hydrogen atom would change an atom's implicit valence, the atom's abnormal valence is set accordingly.
	 * @param hasValid2DCoords pass true, if molecule has 2D-atom-coordinates
	 * @param hasValid3DCoords pass true, if molecule has 3D-atom-coordinates
	 */
	public void removeExplicitHydrogens(boolean hasValid2DCoords, boolean hasValid3DCoords) {
		// calculate EZ and TH stereo parities if 2D-atom coordinates are available
		ensureHelperArrays(hasValid2DCoords ? cHelperParities : cHelperNeighbours);
		mAllAtoms = mAtoms;
		mAllBonds = mBonds;
		for (int atom=0; atom<mAtoms; atom++) {
			if (mAllConnAtoms[atom] != mConnAtoms[atom]) {
				// If we have an abnormal valence implicitly defined by explicit
				// hydrogens, we need to explicitly define that abnormal valence!
				int abnormalValence = getImplicitHigherValence(atom, false);

				mAllConnAtoms[atom] = mConnAtoms[atom];

				if (abnormalValence != -1) {
					int newAbnormalValence = getImplicitHigherValence(atom, true);
					if (abnormalValence != newAbnormalValence) {
						int explicitAbnormalValence = getAtomAbnormalValence(atom);
						if (explicitAbnormalValence == -1 || explicitAbnormalValence < abnormalValence)
							setAtomAbnormalValence(atom, abnormalValence);
						}
					}
				}
			}

		if (hasValid2DCoords)
			setStereoBondsFromParity();

		mValidHelperArrays = cHelperNone;
		}


	private void calculateNeighbours() {
		mConnAtoms = new int[mAllAtoms];
		mAllConnAtoms = new int[mAllAtoms];
		mConnAtom = new int[mAllAtoms][];
		mConnBond = new int[mAllAtoms][];
		mConnBondOrder = new int[mAllAtoms][];
		mPi = new int[mAtoms];

		int[] connCount = new int[mAllAtoms];
		for(int bnd=0; bnd<mAllBonds; bnd++) {
			connCount[mBondAtom[0][bnd]]++;
			connCount[mBondAtom[1][bnd]]++;
			}

		for(int atom=0; atom<mAllAtoms; atom++) {
			mConnAtom[atom] = new int[connCount[atom]];
			mConnBond[atom] = new int[connCount[atom]];
			mConnBondOrder[atom] = new int[connCount[atom]];
			}

		boolean metalBondFound = false;
		for(int bnd=0; bnd<mBonds; bnd++) {
			int order = getBondOrder(bnd);
			if (order == 0) {
				metalBondFound = true;
				continue;
				}

			for (int i = 0; i < 2; i++) {
				int atom = mBondAtom[i][bnd];
				int allConnAtoms = mAllConnAtoms[atom];
				mConnBondOrder[atom][allConnAtoms] = order;
				mConnAtom[atom][allConnAtoms] = mBondAtom[1 - i][bnd];
				mConnBond[atom][allConnAtoms] = bnd;
				mAllConnAtoms[atom]++;	// all non metal-bonded neighbours (non-H, H)
				mConnAtoms[atom]++;	// non-H and non-metal-bonded neighbours
				if (atom < mAtoms) {
					if (order > 1)
						mPi[atom] += order - 1;
					else if (mBondType[bnd] == cBondTypeDelocalized)
						mPi[atom] = 1;
					}
				}
			}

		for(int bnd=mBonds; bnd<mAllBonds; bnd++) {
			int order = getBondOrder(bnd);
			if (order == 0) {
				metalBondFound = true;
				continue;
				}

			for (int i = 0; i < 2; i++) {
				int atom = mBondAtom[i][bnd];
				int allConnAtoms = mAllConnAtoms[atom];
				mConnBondOrder[atom][allConnAtoms] = order;
				mConnAtom[atom][allConnAtoms] = mBondAtom[1 - i][bnd];
				mConnBond[atom][allConnAtoms] = bnd;
				mAllConnAtoms[atom]++;	// all non metal-bonded neighbours (non-H, H)
				if (mBondAtom[1-i][bnd] < mAtoms)
					mConnAtoms[atom]++;
				}
			}

		if (metalBondFound) {
			int[] allConnAtoms = new int[mAllAtoms];
			for(int atom=0; atom<mAllAtoms; atom++)
				allConnAtoms[atom] = mAllConnAtoms[atom];

			for(int bnd=0; bnd<mAllBonds; bnd++) {
				int order = getBondOrder(bnd);
				if (order == 0) {
					for (int i = 0; i < 2; i++) {
						int atom = mBondAtom[i][bnd];
						mConnBondOrder[atom][allConnAtoms[atom]] = order;
						mConnAtom[atom][allConnAtoms[atom]] = mBondAtom[1 - i][bnd];
						mConnBond[atom][allConnAtoms[atom]] = bnd;
						allConnAtoms[atom]++;		// all neighbours (non-H, H, metal-bonded)
						}
					}
				}
			}
		}


	private void findRings(int mode) {
		mRingSet = new RingCollection(this, mode);

		int[] atomRingBondCount = new int[mAtoms];
		for (int bond=0; bond<mBonds; bond++) {
			if (mRingSet.getBondRingSize(bond) != 0) {
				mBondFlags[bond] |= cBondFlagRing;
				atomRingBondCount[mBondAtom[0][bond]]++;
				atomRingBondCount[mBondAtom[1][bond]]++;
				}
			}
		for (int atom=0; atom<mAtoms; atom++) {
			if (atomRingBondCount[atom] == 2)
				mAtomFlags[atom] |= cAtomFlags2RingBonds;
			else if (atomRingBondCount[atom] == 3)
				mAtomFlags[atom] |= cAtomFlags3RingBonds;
			else if (atomRingBondCount[atom] > 3)
				mAtomFlags[atom] |= cAtomFlags4RingBonds;
			}

		for (int ringNo=0; ringNo<mRingSet.getSize(); ringNo++) {
			int ringAtom[] = mRingSet.getRingAtoms(ringNo);
			int ringBond[] = mRingSet.getRingBonds(ringNo);
			int ringAtoms = ringAtom.length;
			for (int i=0; i<ringAtoms; i++) {
				mAtomFlags[ringAtom[i]] |= cAtomFlagSmallRing;
				mBondFlags[ringBond[i]] |= cBondFlagSmallRing;

				if (mBondType[ringBond[i]] == cBondTypeCross)
					mBondType[ringBond[i]] = cBondTypeDouble;
				}
			}
		}


	/**
	 * @return true if hydrogens were deleted and, thus, mConnAtoms are invalid
	 */
	private boolean convertHydrogenToQueryFeatures() {
		if (!mIsFragment)
			return false;

		// if an atom has no free valence then cAtomQFNoMoreNeighbours is not necessary
		// and cAtomQFMoreNeighbours is not possible
		// unless it is an uncharged N- or O-family atom that could be e.g. methylated
		for (int atom=0; atom<mAllAtoms; atom++) {
			if (getFreeValence(atom) <= 0
			 && !(mAtomCharge[atom]==0 && (mAtomicNo[atom]==5 || isNitrogenFamily(atom) || isChalcogene(atom))))
				mAtomQueryFeatures[atom] &= ~(cAtomQFNoMoreNeighbours | cAtomQFMoreNeighbours);
			}

			// approximate explicit hydrogens by query features
			// and remove explicit hydrogens except those with stereo bonds
		boolean deleteHydrogens = false;
		for (int atom=0; atom<mAtoms; atom++) {
			int explicitHydrogens = getExplicitHydrogens(atom);
			if (!mProtectHydrogen && explicitHydrogens > 0) {
				if ((mAtomQueryFeatures[atom] & cAtomQFNoMoreNeighbours) == 0) {
						// add query feature hydrogen to explicit hydrogens
					int queryFeatureHydrogens =
						(mAtomQueryFeatures[atom] & cAtomQFHydrogen) == (cAtomQFNot0Hydrogen | cAtomQFNot1Hydrogen | cAtomQFNot2Hydrogen) ? 3
					  : (mAtomQueryFeatures[atom] & cAtomQFHydrogen) == (cAtomQFNot0Hydrogen | cAtomQFNot1Hydrogen) ? 2
					  : (mAtomQueryFeatures[atom] & cAtomQFNot0Hydrogen) == cAtomQFNot0Hydrogen ? 1 : 0;

					// For atoms with no specific charge definition a potential charge is not considered by getFreeValence().
					// Non-carbon atoms with a charge may have an increased valence. We need to consider that.
					int freeValence = getFreeValence(atom);
					if (mAtomCharge[atom] == 0
					 && (mAtomQueryFeatures[atom] & cAtomQFCharge) == 0
					 && mAtomicNo[atom] != 6)
						freeValence++;

					int queryFeatureShift = explicitHydrogens;
					if (queryFeatureShift > 3 - queryFeatureHydrogens)
						queryFeatureShift = 3 - queryFeatureHydrogens;
					if (queryFeatureShift > freeValence + explicitHydrogens - queryFeatureHydrogens)
						queryFeatureShift = freeValence + explicitHydrogens - queryFeatureHydrogens;

					if (queryFeatureShift > 0) {
						long queryFeatures = (queryFeatureHydrogens == 0) ?  // purge 'less than' options
								0 : (mAtomQueryFeatures[atom] & cAtomQFHydrogen) << queryFeatureShift;
						queryFeatures |= (queryFeatureShift == 3 ? 7 : explicitHydrogens == 2 ? 3 : 1) << cAtomQFHydrogenShift;

						mAtomQueryFeatures[atom] &= ~cAtomQFHydrogen;
						mAtomQueryFeatures[atom] |= (cAtomQFHydrogen & queryFeatures);
						}
					}

				for (int i=mConnAtoms[atom]; i<mAllConnAtoms[atom]; i++) {
					int connBond = mConnBond[atom][i];
					if (mBondType[connBond] == cBondTypeSingle) {	// no stereo bond
						mAtomicNo[mConnAtom[atom][i]] = -1;
						mBondType[connBond] = cBondTypeDeleted;
						deleteHydrogens = true;
						}
					}
				}
			}
		if (deleteHydrogens)
			compressMolTable();

		return deleteHydrogens;
		}


	/**
	 * Normalizes atom query features by removing redundant ones that are implicitly defined
	 * by the structure itself and those that are contradictory.
	 * This method may be used before creating a canonical representation of a substructure.
	 */
	public void validateAtomQueryFeatures() {
		if (!mIsFragment)
			return;

		ensureHelperArrays(cHelperRings);

		for (int atom=0; atom<mAtoms; atom++) {
			if (isRingAtom(atom))
				mAtomQueryFeatures[atom] &= ~(cAtomQFRingSize0 | cAtomQFNotChain);  // forbidden

			if (isAromaticAtom(atom))
				mAtomQueryFeatures[atom] &= ~cAtomQFAromState;   // redundant or impossible
			else if ((mAtomQueryFeatures[atom] & cAtomQFAromatic) != 0)
				mAtomQueryFeatures[atom] &= ~cAtomQFNotAromatic;

			if ((mAtomQueryFeatures[atom] & cAtomQFSmallRingSize) != 0)
				mAtomQueryFeatures[atom] &= ~cAtomQFNotChain;   // redundant

			if (mAtomCharge[atom] != 0)	// explicit charge supersedes query features
				mAtomFlags[atom] &= ~cAtomQFCharge;
			}
		}

	/**
	 * Normalizes bond query features by removing redundant ones that are implicitly defined
	 * by the structure itself and those that are contradictory.
	 * This method may be used before creating a canonical representation of a substructure.
	 */
	public void validateBondQueryFeatures() {
		if (!mIsFragment)
			return;

		ensureHelperArrays(cHelperRings);

		for (int bond=0; bond<mBonds; bond++) {
			if (isDelocalizedBond(bond))
				mBondQueryFeatures[bond] &= ~cBondQFDelocalized;

			int bondType = mBondType[bond] & cBondTypeMaskSimple;
			if (bondType == cBondTypeSingle)
				mBondQueryFeatures[bond] &= ~cBondQFSingle;
			else if (bondType == cBondTypeDouble)
				mBondQueryFeatures[bond] &= ~cBondQFDouble;
			else if (bondType == cBondTypeTriple)
				mBondQueryFeatures[bond] &= ~cBondQFTriple;
			else if (bondType == cBondTypeQuadruple)
				mBondQueryFeatures[bond] &= ~cBondQFQuadruple;
			else if (bondType == cBondTypeQuintuple)
				mBondQueryFeatures[bond] &= ~cBondQFQuintuple;
			else if (bondType == cBondTypeMetalLigand)
				mBondQueryFeatures[bond] &= ~cBondQFMetalLigand;
			else if (bondType == cBondTypeDelocalized)
				mBondQueryFeatures[bond] &= ~cBondQFDelocalized;
			}
		}


	private void writeObject(ObjectOutputStream stream) throws IOException {}
	private void readObject(ObjectInputStream stream) throws IOException {}


	public final static Coordinates getCenterGravity(ExtendedMolecule mol) {

		int n = mol.getAllAtoms();

		int [] indices = new int [n];

		for (int i = 0; i < indices.length; i++) {
			indices[i]=i;
		}

		return getCenterGravity(mol, indices);
	}

	public final static Coordinates getCenterGravity(ExtendedMolecule mol, int[] indices) {

		Coordinates c = new Coordinates();
		for (int i = 0; i < indices.length; i++) {
			c.x += mol.getAtomX(indices[i]);
			c.y += mol.getAtomY(indices[i]);
			c.z += mol.getAtomZ(indices[i]);
		}
		c.x /= indices.length;
		c.y /= indices.length;
		c.z /= indices.length;

		return c;
	}

}
