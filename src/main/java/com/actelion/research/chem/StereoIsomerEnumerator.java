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

import com.actelion.research.chem.coords.CoordinateInventor;

import java.util.TreeMap;

public class StereoIsomerEnumerator {
	private boolean mSkipEnantiomers;
	private StereoMolecule mMol;
	private int[][] mAtomGroupList,mBondGroupList,mUnknownDoubleBondList;
	private boolean[][] mAtomIsParity1,mBondIsParity1;

	/**
	 * If the passed molecule has stereo-chemically undefined configurations
	 * (double bonds, stereo centers) or/and one or more AND/OR groups of
	 * defined relative stereo configurations, then it represents multiple
	 * stereo isomers. The StereoIsomerEnumerator generates all individual
	 * stereo isomers of the passed molecule. If the passed molecule does
	 * not include absolute stereo centers (or atrop isomeric configuration),
	 * but unknown stereo centers or groups with defined relative configuration,
	 * then we have pairs of enantiomers. In this case the StereoIsomerEnumerator
	 * may either generate one or both stereo isomers of each enantiomeric pair.
	 * @param mol WARNING
	 * @param skipEnantiomers
	 */
	public StereoIsomerEnumerator(StereoMolecule mol, boolean skipEnantiomers) {
		mMol = mol;
		mMol.ensureHelperArrays(Molecule.cHelperParities);
		TreeMap<String,int[]> groupMap = new TreeMap<String,int[]>();
		mSkipEnantiomers = skipEnantiomers;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.getAtomParity(atom) != Molecule.cAtomParityNone && !mMol.isAtomParityPseudo(atom)) {
				if (mMol.getAtomParity(atom) == Molecule.cAtomParityUnknown) {
					int[] atomList = new int[1];
					atomList[0] = atom;
					groupMap.put("U"+atom, atomList);
					}
				else {
					if (mMol.getAtomESRType(atom) == Molecule.cESRTypeAbs) {
						mSkipEnantiomers = false;
						}
					else {
						String key = (mMol.getAtomESRType(atom) == Molecule.cESRTypeOr ? "O" : "A") + mMol.getAtomESRGroup(atom);
						int[] atomList = groupMap.get(key);
						if (atomList == null) {
							atomList = new int[1];
							atomList[0] = atom;
							groupMap.put(key, atomList);
							}
						else {
							int[] newAtomList = new int[atomList.length+1];
							for (int i=0; i<atomList.length; i++)
								newAtomList[i] = atomList[i];
							newAtomList[atomList.length] = atom;
							groupMap.put(key, newAtomList);
							}
						}
					}
				}
			}
		mAtomGroupList = groupMap.values().toArray(new int[0][]);
		mAtomIsParity1 = new boolean[mAtomGroupList.length][];
		for (int i=0; i<mAtomGroupList.length; i++) {
			mAtomIsParity1[i] = new boolean[mAtomGroupList[i].length];
			for (int j=0; j<mAtomGroupList[i].length; j++)
				mAtomIsParity1[i][j] = (mMol.getAtomParity(mAtomGroupList[i][j]) == Molecule.cAtomParity1);
			}

		// handle now atrop isomeric bonds (BINAP type)
		groupMap.clear();
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondParity(bond) != Molecule.cBondParityNone && !mMol.isBondParityPseudo(bond) && mMol.getBondOrder(bond) == 1) {
				if (mMol.getBondParity(bond) == Molecule.cBondParityUnknown) {	// atrop and double bonds
					int[] bondList = new int[1];
					bondList[0] = bond;
					groupMap.put("U"+bond, bondList);
					}
				else {	// handle ESR for atrop bonds
					if (mMol.getBondESRType(bond) == Molecule.cESRTypeAbs) {
						mSkipEnantiomers = false;
						}
					else {
						String key = (mMol.getBondESRType(bond) == Molecule.cESRTypeOr ? "O" : "A") + mMol.getBondESRGroup(bond);
						int[] bondList = groupMap.get(key);
						if (bondList == null) {
							bondList = new int[1];
							bondList[0] = bond;
							groupMap.put(key, bondList);
							}
						else {
							int[] newBondList = new int[bondList.length+1];
							for (int i=0; i<bondList.length; i++)
								newBondList[i] = bondList[i];
							newBondList[bondList.length] = bond;
							groupMap.put(key, newBondList);
							}
						}
					}
				}
			}
		mBondGroupList = groupMap.values().toArray(new int[0][]);
		mBondIsParity1 = new boolean[mBondGroupList.length][];
		for (int i=0; i<mBondGroupList.length; i++) {
			mBondIsParity1[i] = new boolean[mBondGroupList[i].length];
			for (int j=0; j<mBondGroupList[i].length; j++)
				mBondIsParity1[i][j] = (mMol.getBondParity(mBondGroupList[i][j]) == Molecule.cBondParityEor1);
			}

		if (mAtomGroupList.length == 0 && mBondGroupList.length == 0)
			mSkipEnantiomers = false;

		// handle now double bonds with unknown configuration
		groupMap.clear();
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondParity(bond) == Molecule.cBondParityUnknown
			 && !mMol.isBondParityPseudo(bond)
			 && mMol.getBondOrder(bond) == 2) {
				int[] bondList = new int[1];
				bondList[0] = bond;
				groupMap.put("U"+bond, bondList);
				}
			}
		mUnknownDoubleBondList = groupMap.values().toArray(new int[0][]);
		}

	/**
	 * Checks and returns, whether skipEnantiomers was chosen in the constructor
	 * and whether the complete set of isomers would contain enantiomeric sets,
	 * because of the absence of absolute stereo centers and the presence of
	 * unknown or relative tetrahedral or atrop stereo configurations.
	 * @return whether we are skipping enantiomers
	 */
	public boolean isSkippingEnantiomers() {
		return mSkipEnantiomers;
		}

	/**
	 * This calculates the stereo isomer count that may be requested by the getStereoIsomer() method.
	 * If isSkippingEnantiomers() returns true, then the enumerated stereo isomers don't
	 * contain enatiomers, i.e. in order to complete the set of stereo isomers, one needs
	 * to duplicate and invert externally in this case.
	 * @return stereo isomer count that may be requested by the getStereoIsomer() method.
	 */
	public int getStereoIsomerCount() {
		return 1 << (mAtomGroupList.length+mBondGroupList.length+mUnknownDoubleBondList.length-(mSkipEnantiomers?1:0));
		}

	/**
	 * Creates and returns the molecule in a specific atom and bond parity permutation state,
	 * i.e. a specific stereo isomer. If stereo bonds are permuted (as opposed to tetrahedral centers),
	 * then new atom coordinates are generated to reflect the specific bond configuration.
	 * Otherwise, only atom parities are set and up/down stereo bonds are set to reflect the
	 * expected correct atom stereo configuration.<br>
	 * This method does not sort out stereo isomers that are impossible because of geometric strain.
	 * @param index parity permutation
	 * @return
	 */
	public StereoMolecule getStereoIsomer(int index) {
		StereoMolecule mol = new StereoMolecule(mMol);

		// if we produce one enantiomer only, we don't permute the first atom group
		// (or bond group if there are no atom groups)
		boolean skipFirst = mSkipEnantiomers;
		for (int i=0; i<mAtomGroupList.length; i++) {
			boolean invertGroup = !skipFirst && ((index & 1) != 0);
			for (int j=0; j<mAtomGroupList[i].length; j++) {
				int parity = (invertGroup ^ mAtomIsParity1[i][j]) ? Molecule.cAtomParity2 : Molecule.cAtomParity1;
				int atom = mAtomGroupList[i][j];
				mol.setAtomParity(atom, parity, false);
				mol.setAtomESR(atom, Molecule.cESRTypeAbs, 0);
				}
			if (skipFirst) {
				skipFirst = false;
				continue;
				}
			index >>= 1;
			}

		for (int i=0; i<mBondGroupList.length; i++) {
			boolean invertGroup = !skipFirst && ((index & 1) != 0);
			for (int j=0; j<mBondGroupList[i].length; j++) {
				int parity = (invertGroup ^ mBondIsParity1[i][j]) ? Molecule.cBondParityZor2 : Molecule.cBondParityEor1;
				int bond = mBondGroupList[i][j];
				mol.setBondParity(bond, parity, false);
				mol.setBondESR(bond, Molecule.cESRTypeAbs, 0);
				}
			if (skipFirst) {
				skipFirst = false;
				continue;
				}
			index >>= 1;
			}

		for (int i=0; i<mUnknownDoubleBondList.length; i++) {
			int parity = ((index & 1) != 0) ? Molecule.cBondParityZor2 : Molecule.cBondParityEor1;
			int bond = mUnknownDoubleBondList[i][0];
			mol.setBondType(bond, Molecule.cBondTypeDouble);
			mol.setBondParity(bond, parity, false);
			index >>= 1;
			}

		mol.setParitiesValid(0);
		mol.setStereoBondsFromParity();

		// if we have unknown double bonds, we need to generate new 2D-coordinates for E and Z configuration
		if (mUnknownDoubleBondList.length != 0) {
			mol.setParitiesValid(0);
			new CoordinateInventor(0).invent(mol);
			}

		return mol;
		}

	/**
	 * Checks, whether the given molecule's atom and bond parities are the expected one's
	 * for the given stereo isomer index. This may be useful when creating conformers
	 * from the stereo isomer returned by getStereoIsomer() to check, whether conformers
	 * still have the correct stereo configurations. If the correct stereo configuration is
	 * an impossible one due to ring strains, then conformer generators or force field may
	 * change the configuration in order to create a conformer at all.
	 * @param mol
	 * @param index
	 * @return true, if the given molecule's atom and bond parities are correct the given index
	 */
	public boolean isCorrectStereoIsomer(StereoMolecule mol, int index) {
		mol. ensureHelperArrays(Molecule.cHelperParities);

		// if we produce one enantiomer only, we don't permute the first atom group
		// (or bond group if there are no atom groups)
		boolean skipFirst = mSkipEnantiomers;
		for (int i=0; i<mAtomGroupList.length; i++) {
			boolean invertGroup = !skipFirst && ((index & 1) != 0);
			for (int j=0; j<mAtomGroupList[i].length; j++) {
				int parity = (invertGroup ^ mAtomIsParity1[i][j]) ? Molecule.cAtomParity2 : Molecule.cAtomParity1;
				int atom = mAtomGroupList[i][j];
				if (mol.getAtomParity(atom) != parity)
					return false;
				}
			if (skipFirst) {
				skipFirst = false;
				continue;
				}
			index >>= 1;
			}

		for (int i=0; i<mBondGroupList.length; i++) {
			boolean invertGroup = !skipFirst && ((index & 1) != 0);
			for (int j=0; j<mBondGroupList[i].length; j++) {
				int parity = (invertGroup ^ mBondIsParity1[i][j]) ? Molecule.cBondParityZor2 : Molecule.cBondParityEor1;
				int bond = mBondGroupList[i][j];
				if (mol.getBondParity(bond) != parity)
					return false;
				}
			if (skipFirst) {
				skipFirst = false;
				continue;
				}
			index >>= 1;
			}

		for (int i=0; i<mUnknownDoubleBondList.length; i++) {
			int parity = ((index & 1) != 0) ? Molecule.cBondParityZor2 : Molecule.cBondParityEor1;
			int bond = mUnknownDoubleBondList[i][0];
			if (mol.getBondParity(bond) !=  parity)
				return false;
			index >>= 1;
			}

		return true;
		}
	}

