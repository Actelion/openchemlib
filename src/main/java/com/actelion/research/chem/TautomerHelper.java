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
*/

package com.actelion.research.chem;

import com.actelion.research.util.IntArrayComparator;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeSet;

public class TautomerHelper {
	private static final int MAX_TAUTOMERS = 100000;

	private StereoMolecule mOriginalMol;
	private boolean[] mIsTautomerBond;
	private boolean[] mHasFreeValence;
	private int[] mRegionPiCount;
	private int[] mRegionDCount;
	private int[] mRegionTCount;

	private Iterator<BondOrders> mBondOrderIterator;
	private TreeSet<BondOrders> mBondOrderSet;
	private ArrayDeque<BondOrders> mBondOrderDeque;

	/**
	 * @param mol
	 */
	public TautomerHelper(StereoMolecule mol) {
		mOriginalMol = mol;
		mOriginalMol.ensureHelperArrays(Molecule.cHelperParities);

		mIsTautomerBond = new boolean[mOriginalMol.getBonds()];
		mHasFreeValence = new boolean[mOriginalMol.getAtoms()];
		for (int i=0; i<mOriginalMol.getAtoms(); i++) {
			int atomicNo = mOriginalMol.getAtomicNo(i);
			int valence = (atomicNo < Molecule.cAtomValence.length && Molecule.cAtomValence[atomicNo] != null) ?
					Molecule.cAtomValence[atomicNo][0] : Molecule.cDefaultAtomValence;
			mHasFreeValence[i] = (mOriginalMol.getNonHydrogenNeighbourCount(i) < valence);
			}

		createAllTautomers();
		}

	/**
	 * @param tautomer null or another tautomer of same molecule, which receives the new bond orders
	 * @return
	 */
	public StereoMolecule getNextTautomer(StereoMolecule tautomer) {
		if (mBondOrderIterator == null)
			mBondOrderIterator = mBondOrderSet.iterator();

		if (!mBondOrderIterator.hasNext())
			return null;

		if (tautomer == null)
			tautomer = mOriginalMol.getCompactCopy();

		mBondOrderIterator.next().copyToTautomer(tautomer);

		return tautomer;
		}

	public int getTautomerCount() {
		return mBondOrderSet.size();
		}

	/**
	 * Identifies connected tautomeric regions and assign region numbers to all atoms.
	 * Atoms sharing the same region share the same number.<br>
	 * 0: not member of a tautomer region; 1 and above: region number
	 * @param atomRegionNo int[mol.getAtoms()] filled with 0
	 * @return region count
	 */
	public int getAtomRegionNumbers(int[] atomRegionNo) {
		if (mBondOrderSet.size() == 1)
			return 0;

		int regionCount = assignRegionNumbers(atomRegionNo);
		return regionCount;
		}

	/**
	 * Considers connected tautomer bonds to belong to a tautomer region.
	 * All independent tautomer regions are located and member atoms assigned to them.
	 * mAtomRegionNo[] is set accordingly.
	 * 0: not member of a tautomer region; 1 and above: region number
	 * @param atomRegionNo int[mol.getAtoms()] filled with 0
	 * @return number of found tautomer regions
	 */
	private int assignRegionNumbers(int[] atomRegionNo) {
		int[] graphAtom = new int[mOriginalMol.getAtoms()];
		boolean[] bondWasSeen = new boolean[mOriginalMol.getBonds()];
		int region = 0;

		for (int bond=0; bond<mOriginalMol.getBonds(); bond++) {
			if (!bondWasSeen[bond] && mIsTautomerBond[bond]) {
				region++;
				atomRegionNo[mOriginalMol.getBondAtom(0, bond)] = region;
				atomRegionNo[mOriginalMol.getBondAtom(1, bond)] = region;
				bondWasSeen[bond] = true;
				for (int i=0; i<2; i++) {
					int atom = mOriginalMol.getBondAtom(i, bond);
					atomRegionNo[atom] = region;
					int current = 0;
					int highest = 0;
					graphAtom[0] = atom;
					while (current <= highest) {
						for (int j=0; j<mOriginalMol.getConnAtoms(graphAtom[current]); j++) {
							int connBond = mOriginalMol.getConnBond(graphAtom[current], j);
							if (!bondWasSeen[connBond] && mIsTautomerBond[connBond]) {
								bondWasSeen[connBond] = true;
								int connAtom = mOriginalMol.getConnAtom(graphAtom[current], j);
								if (atomRegionNo[connAtom] == 0) {
									atomRegionNo[connAtom] = region;
									graphAtom[++highest] = connAtom;
									}
								}
							}
						current++;
						}
					}
				}
			}

		return region;
		}

	/**
	 * Counts for every region: pi-electrons, deuterium atoms, tritium atoms.
	 * Must be called after assignRegionNumbers().
	 * @param atomRegionNo array with valid region numbers
	 */
	private void compileRegionCounts(int[] atomRegionNo, int regionCount) {
		mRegionPiCount = new int[regionCount];
		mRegionDCount = new int[regionCount];
		mRegionTCount = new int[regionCount];

		for (int atom=0; atom<mOriginalMol.getAtoms(); atom++) {
			if (atomRegionNo[atom] != 0) {
				int regionIndex = atomRegionNo[atom]-1;
				for (int i=0; i<mOriginalMol.getConnAtoms(atom); i++) {
					int connAtom = mOriginalMol.getConnAtom(atom, i);
					if (mOriginalMol.getAtomicNo(connAtom) == 1) {
						if (mOriginalMol.getAtomMass(connAtom) == 2)
							mRegionDCount[regionIndex]++;
						if (mOriginalMol.getAtomMass(connAtom) == 3)
							mRegionTCount[regionIndex]++;
						}
					}
				}
			}
		for (int bond=0; bond<mOriginalMol.getBonds(); bond++) {
			if (mIsTautomerBond[bond]
			 && mOriginalMol.getBondOrder(bond) == 2) {
				mRegionPiCount[atomRegionNo[mOriginalMol.getBondAtom(0, bond)]-1] += 2;
				}
			}
		}

	private boolean hasAcidicHydrogen(StereoMolecule mol, int atom) {
		if (mol.getAllHydrogens(atom) <= 0)
			return false;

		if (mol.isElectronegative(atom))
			return true;

		if (mol.getAtomPi(atom) != 0)
			return false;

		return true;
		}

	/**
	 * If no tautomers can be formed then a copy of the original molecule is returned.
	 * Otherwise a copy of the original molecule is normalized to create a
	 * generic tautomer structure. The original molecule is not touched.
	 * Different tautomers of the same molecule should always result in the same
	 * generic tautomer structure. A generic tautomer contains one or more regions
	 * indicated by bond query features cBondQFSingle & cBondQFDouble. Bond types
	 * of all bonds of any tautomer region are cBondTypeSingle.
	 * The highest ranking atom in every region carries a label defining the number
	 * of double bonds and D and T atoms. The returned molecule has the fragment bit set.
	 * Canonicalizing the returned molecule with Canonizer mode ENCODE_ATOM_CUSTOM_LABELS
	 * produces the same idcode from any tautomer.
	 * If keepStereoCenters is true, then stereo centers with parity 1 or 2, if they are
	 * absolute or if they are part of an AND/OR group with more than one member,
	 * are considered stable (non racemising) and, thus, their proton is not considered
	 * being able to take part in a tautomeric transition.
	 * @return generic tautomer with normalized tautomer regions and custom label to encode pi,D,T counts
	 */
	public StereoMolecule createGenericTautomer() {
		if (mBondOrderSet.size() == 1)
			return mOriginalMol;

		StereoMolecule mol = mOriginalMol.getCompactCopy();
		mol.ensureHelperArrays(Molecule.cHelperRings);

		mol.setFragment(true);
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (mIsTautomerBond[bond]) {
				mol.setBondType(bond, Molecule.cBondTypeSingle);
				mol.setBondQueryFeature(bond, Molecule.cBondQFSingle | Molecule.cBondQFDouble, true);
				}
			}

		int[] atomRegionNo = new int[mol.getAtoms()];
		int regionCount = assignRegionNumbers(atomRegionNo);

		// remove stereo information from bonds that indicate a stereo center at one of the tautomer region atoms
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (atomRegionNo[atom] != 0
			 && mOriginalMol.getNonHydrogenNeighbourCount(atom) < 4) {   // keep terternary S or P
				mol.convertStereoBondsToSingleBonds(atom);
				mol.setAtomConfigurationUnknown(atom,false);
				mol.setAtomESR(atom, Molecule.cESRTypeAbs, -1);
				}
			}

		// find highest ranking atom in every region
		int[] maxAtom = new int[regionCount];
		int[] maxRank = new int[regionCount];
		int[] atomRank = new Canonizer(mol).getFinalRank();
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (atomRegionNo[atom] != 0) {
				int regionIndex = atomRegionNo[atom]-1;
				if (maxRank[regionIndex] < atomRank[atom]) {
					maxRank[regionIndex] = atomRank[atom];
					maxAtom[regionIndex] = atom;
					}
				}
			}

		// attach label with region counts to highest ranking atoms
		compileRegionCounts(atomRegionNo, regionCount);
		for (int i=0; i<regionCount; i++) {
			String label = ""+mRegionPiCount[i]+"|"+mRegionDCount[i]+"|"+mRegionTCount[i];
			mol.setAtomCustomLabel(maxAtom[i], label);
			}

		return mol;
		}

	private void createAllTautomers() {
		mBondOrderSet = new TreeSet<>();
		mBondOrderDeque = new ArrayDeque<>();

		addTautomerIfNew(new BondOrders(mOriginalMol));

		// This serves as recycled molecule container
		StereoMolecule tautomer = mOriginalMol.getCompactCopy();

		while (!mBondOrderDeque.isEmpty()) {
			mBondOrderDeque.poll().copyToTautomer(tautomer);
			addAllTautomers(tautomer);

			if (mBondOrderSet.size() >= MAX_TAUTOMERS) {
				System.out.println("Tautomer count exceeds maximum: "+new Canonizer(mOriginalMol).getIDCode());
				break;
				}
			}

		// TODO racemize stereo centers
		}

	/**
	 * Find all HX-Y=Z / X=Y-ZH type 3-atom sequences and recursively vinylogous sequences.
	 */
	private void addAllTautomers(StereoMolecule mol) {
		ArrayList<Integer> bondList = new ArrayList<>();
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		boolean [] isUsedAtom = new boolean[mol.getAtoms()];   // atom use buffer for recursive methods

		for (int atom1=0; atom1<mol.getAtoms(); atom1++) {
			if (isValidHeteroAtom(atom1)) {
				isUsedAtom[atom1] = true;
				for (int i=0; i<mol.getConnAtoms(atom1); i++) {
					int atom2 = mol.getConnAtom(atom1, i);
					int bond12 = mol.getConnBond(atom1, i);
					int order12 = mol.getConnBondOrder(atom1, i);
					if (mol.getAtomPi(atom2) != 0
					 && mol.getAtomPi(atom1) < order12) {   // make sure atom1 has no other double bond
						isUsedAtom[atom2] = true;
						bondList.add(bond12);
						for (int j=0; j<mol.getConnAtoms(atom2); j++) {
							int atom3 = mol.getConnAtom(atom2, j);
							if (!isUsedAtom[atom3]) {
								isUsedAtom[atom3] = true;
								int bond23 = mol.getConnBond(atom2, j);
								int order23 = mol.getConnBondOrder(atom2, j);
								if (mol.getAtomPi(atom2) + 2 == order12 + order23) {
									bondList.add(bond23);

									// if atom3 has other double bond we check for vinylogous donor atoms
									if (order12 >= order23) {
										if (mol.getAtomPi(atom3) < order23) {   // if atom3 has no other double bond
											if (hasAcidicHydrogen(mol, atom3))
												addVinylogousTautomers(mol, atom3, true, false, isUsedAtom, bondList);
											}
										else {
											addVinylogousTautomers(mol, atom3, true, true, isUsedAtom, bondList);
											}
										}

									if (order23 >= order12 && hasAcidicHydrogen(mol, atom1)) {
										addVinylogousTautomers(mol, atom3, false, false, isUsedAtom, bondList);
										}

									if (isValidDonorAtom(atom3)
									 && mol.getAtomPi(atom3) < order23) {   // make sure atom3 has no other double bond
										if (order12<=2 && order23>=2 && hasAcidicHydrogen(mol, atom1)) {
											addDirectTautomer(mol, bond12, bond23);
											}

										if (order12>=2 && order23<=2 && hasAcidicHydrogen(mol, atom3)) {
											addDirectTautomer(mol, bond23, bond12);
											}
										}

									bondList.remove(bondList.size()-1);
									}
								isUsedAtom[atom3] = false;
								}
							}
						bondList.remove(bondList.size()-1);
						isUsedAtom[atom2] = false;
						}
					}
				isUsedAtom[atom1] = false;
				}
			}
		}

	private boolean isValidDonorAtom(int atom) {
		return mHasFreeValence[atom]
			&& (mOriginalMol.getAtomicNo(atom) == 5
			 || mOriginalMol.getAtomicNo(atom) == 6
			 || mOriginalMol.getAtomicNo(atom) == 7
			 || mOriginalMol.getAtomicNo(atom) == 8
			 || mOriginalMol.getAtomicNo(atom) == 16);
		}

	private boolean isValidHeteroAtom(int atom) {
		return mHasFreeValence[atom]
			&& (mOriginalMol.getAtomicNo(atom) == 7
			 || mOriginalMol.getAtomicNo(atom) == 8
			 || mOriginalMol.getAtomicNo(atom) == 16);
		}

	private void addDirectTautomer(StereoMolecule mol, int bondSToD, int bondDToS) {
		BondOrders bondOrders = new BondOrders(mol);
		bondOrders.setBond(bondSToD, mol.getBondOrder(bondSToD) == 1 ? 2 : 3);
		bondOrders.setBond(bondDToS, mol.getBondOrder(bondDToS) == 2 ? 1 : 2);
		mIsTautomerBond[bondSToD] = true;
		mIsTautomerBond[bondDToS] = true;
		addTautomerIfNew(bondOrders);
		}

	private boolean addVinylogousTautomers(StereoMolecule mol, int atom1, boolean firstBondIsDouble, boolean thirdBondIsDouble, boolean[] isUsedAtom, ArrayList<Integer> bondList) {
		for (int i=0; i<mol.getConnAtoms(atom1); i++) {
			int atom2 = mol.getConnAtom(atom1, i);
			if (!isUsedAtom[atom2]) {
				int bond12 = mol.getConnBond(atom1, i);
				int order12 = mol.getBondOrder(bond12);
				if ((firstBondIsDouble && order12 >= 2)
				 || (!firstBondIsDouble && order12 <= 2)) {
					isUsedAtom[atom2] = true;
					bondList.add(bond12);
					for (int j=0; j<mol.getConnAtoms(atom2); j++) {
						int atom3 = mol.getConnAtom(atom2, j);
						if (!isUsedAtom[atom3]) {
							int bond23 = mol.getConnBond(atom2, j);
							int order23 = mol.getBondOrder(bond23);
							if (mol.getAtomPi(atom2) + 2 == order12 + order23
							 && ((firstBondIsDouble && order23 <= 2)
							  || (!firstBondIsDouble && order23 >= 2))) {
								isUsedAtom[atom3] = true;
								bondList.add(bond23);
								if (isValidDonorAtom(atom3) && (!firstBondIsDouble || hasAcidicHydrogen(mol, atom3))) {
									BondOrders bondOrders = new BondOrders(mol);
									for (int k=0; k<bondList.size(); k++) {
										int bond = bondList.get(k);
										boolean makeDouble = (k < 2) ? (firstBondIsDouble ^ (k & 1) == 0)
												: (thirdBondIsDouble ^ (k & 1) == 0);
										if (makeDouble)
											bondOrders.setBond(bond, mol.getBondOrder(bond) == 1 ? 2 : 3);
										else
											bondOrders.setBond(bond, mol.getBondOrder(bond) == 2 ? 1 : 2);
										mIsTautomerBond[bond] = true;
										}
									addTautomerIfNew(bondOrders);
									}
								else {
									addVinylogousTautomers(mol, atom3, firstBondIsDouble, thirdBondIsDouble, isUsedAtom, bondList);
									}
								bondList.remove(bondList.size()-1);
								isUsedAtom[atom3] = false;
								}
							}
						}
					bondList.remove(bondList.size()-1);
					isUsedAtom[atom2] = false;
					}
				}
			}
		return false;
		}

	private void addTautomerIfNew(BondOrders bondOrders) {
		if (mBondOrderSet.add(bondOrders))
			mBondOrderDeque.add(bondOrders);
		}

	class BondOrders implements Comparable<BondOrders> {
		private int[] encoding;

		public BondOrders(StereoMolecule mol) {
			encoding = new int[(mOriginalMol.getBonds()+15) / 16];
			for (int i=0; i<mOriginalMol.getBonds(); i++)
				encoding[i >> 4] |= (mol.getBondOrder(i) << (2*(i & 15)));
			}

		@Override
		public int compareTo(BondOrders o) {
			return new IntArrayComparator().compare(encoding, o.encoding);
			}

		public void setBond(int bond, int order) {
			int high = bond >> 4;
			int shift = 2 * (bond & 15);
			encoding[high] &= ~(3 << shift);
			encoding[high] |= (order << shift);
			}

		public void copyToTautomer(StereoMolecule tautomer) {
			for (int i=0; i<mOriginalMol.getBonds(); i++) {
				int bo = 3 & (encoding[i >> 4] >> (2*(i & 15)));
				tautomer.setBondType(i, bo == 1 ? Molecule.cBondTypeSingle
						: bo == 2 ? Molecule.cBondTypeDouble
						: bo == 3 ? Molecule.cBondTypeTriple
						: Molecule.cBondTypeMetalLigand);
				}
			}
		}
	}
