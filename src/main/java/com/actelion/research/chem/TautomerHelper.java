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

import com.actelion.research.util.IntArrayComparator;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeSet;

public class TautomerHelper {
	private static final int MAX_TAUTOMERS = 100000;
	private static int sMaxTautomers = MAX_TAUTOMERS;
	private static boolean sSuppressWarning = false;

	private StereoMolecule mOriginalMol;
	private boolean[] mIsTautomerBond;
	private boolean[] mHasFreeValence;
	private int[] mRegionPiCount;
	private int[] mRegionDCount;
	private int[] mRegionTCount;
	private int[] mAtomRegionNo;
	private int[] mAtomDCount;
	private int[] mAtomTCount;
	private int mRegionCount;

	private Iterator<BondOrders> mBondOrderIterator;
	private TreeSet<BondOrders> mBondOrderSet;
	private ArrayDeque<BondOrders> mBondOrderDeque;

	public static void setMaxTautomers(int maxTautomers) {
		sMaxTautomers = maxTautomers;
		}

	public static void setSuppressWarning(boolean suppressWarning) {
		sSuppressWarning = suppressWarning;
		}

	/**
	 * @param mol
	 */
	public TautomerHelper(StereoMolecule mol) {
		mOriginalMol = mol.getCompactCopy();
		moveDeuteriumAndTritiumToTableEnd();
		mOriginalMol.ensureHelperArrays(Molecule.cHelperRings);

		mIsTautomerBond = new boolean[mOriginalMol.getBonds()];
		mHasFreeValence = new boolean[mOriginalMol.getAtoms()];
		for (int i=0; i<mOriginalMol.getAtoms(); i++) {
			int valence = Molecule.getAllowedValences(mOriginalMol.getAtomicNo(i))[0];
			mHasFreeValence[i] = (mOriginalMol.getNonHydrogenNeighbourCount(i) < valence);
			}

		createAllTautomers();

		assignRegionNumbers();
		countAndRemoveDAndT();
		compileRegionCounts();
		}

	/**
	 *
	 * @param tautomer null or molecule container that received the new tautomer
	 * @return
	 */
	public StereoMolecule getNextTautomer(StereoMolecule tautomer) {
		if (mBondOrderIterator == null)
			mBondOrderIterator = mBondOrderSet.iterator();

		if (!mBondOrderIterator.hasNext())
			return null;

		if (tautomer != null) {
			tautomer.clear();
			mOriginalMol.copyMolecule(tautomer);
			}
		else {
			tautomer = mOriginalMol.getCompactCopy();
			}

		mBondOrderIterator.next().copyToTautomer(tautomer);

		// Now we add originally removed D and T at arbitrary positions
		if (mAtomDCount != null || mAtomTCount != null) {
			tautomer.ensureHelperArrays(Molecule.cHelperNeighbours);
			int[] freeValence = new int[tautomer.getAtoms()];
			for (int atom=0; atom<tautomer.getAtoms(); atom++)
				freeValence[atom] = tautomer.getFreeValence(atom);

			// Add D and T in two steps. First try to attach then to the atoms, where they have been originally.
			// If limited free valence prevents that (because now there is a double bond, where it was single),
			// then add remaining D and to other atoms in same region, which gained a free valence.
			if (mAtomDCount != null) {
				int[] regionDCount = mRegionDCount.clone();
				for (int atom=0; atom<tautomer.getAtoms(); atom++) {
					int dCount = Math.min(freeValence[atom], mAtomDCount[atom]);
					for (int i=0; i<dCount; i++) {
						addHydrogen(tautomer, atom, 2);
						freeValence[atom]--;
						regionDCount[mAtomRegionNo[atom]-1]--;
						}
					}
				for (int atom=0; atom<tautomer.getAtoms(); atom++) {
					if (mAtomRegionNo[atom] != 0 && regionDCount[mAtomRegionNo[atom]-1] > 0
					 && tautomer.getAtomPi(atom) < mOriginalMol.getAtomPi(atom)
					 && freeValence[atom] != 0) {
						addHydrogen(tautomer, atom, 2);
						freeValence[atom]--;
						regionDCount[mAtomRegionNo[atom]-1]--;
						}
					}
				}

			if (mAtomTCount != null) {
				int[] regionTCount = mRegionTCount.clone();
				for (int atom=0; atom<tautomer.getAtoms(); atom++) {
					int tCount = Math.min(freeValence[atom], mAtomTCount[atom]);
					for (int i=0; i<tCount; i++) {
						addHydrogen(tautomer, atom, 3);
						freeValence[atom]--;
						regionTCount[mAtomRegionNo[atom]-1]--;
						}
					}
				for (int atom=0; atom<tautomer.getAtoms(); atom++) {
					if (mAtomRegionNo[atom] != 0 && regionTCount[mAtomRegionNo[atom]-1] > 0
					 && tautomer.getAtomPi(atom) < mOriginalMol.getAtomPi(atom)
					 && freeValence[atom] != 0) {
						addHydrogen(tautomer, atom, 3);
						freeValence[atom]--;
						regionTCount[mAtomRegionNo[atom]-1]--;
						}
					}
				}
			}

		// Double bonds, which are single bonds in other tautomers, and which carry two different
		// substituents on both double bond atoms, may have two configurations (E and Z),
		// unless the double bond is in a small ring.
		// In these cases we use a cross-bond instead of creating both configurations,
		// to avoid another dimension of tautomer count explosion.
		// Here we need to convert erroneous cross-bond assignments back to double, if the bond isn't a stereo bond.
		tautomer.ensureHelperArrays(Molecule.cHelperParities);
		for (int bond=0; bond<tautomer.getBonds(); bond++)
			if (tautomer.getBondType(bond) == Molecule.cBondTypeCross
			 && tautomer.getBondParity(bond) == Molecule.cBondParityNone)
				tautomer.setBondType(bond, Molecule.cBondTypeDouble);

		return tautomer;
		}

	private void moveDeuteriumAndTritiumToTableEnd() {
		mOriginalMol.ensureHelperArrays(Molecule.cHelperNeighbours);

		int lastNonHAtom = mOriginalMol.getAtoms();
		do lastNonHAtom--;
		while ((lastNonHAtom >= 0) && mOriginalMol.getAtomicNo(lastNonHAtom) == 1);

		for (int atom=0; atom<lastNonHAtom; atom++) {
			if (mOriginalMol.getAtomicNo(atom) == 1) {
				mOriginalMol.swapAtoms(atom, lastNonHAtom);

				do lastNonHAtom--;
				while (mOriginalMol.getAtomicNo(lastNonHAtom) == 1);
				}
			}

		if (lastNonHAtom == mOriginalMol.getAtoms()-1)
			return;

		boolean isHydrogenBond[] = new boolean[mOriginalMol.getBonds()];
		for (int bond=0; bond<mOriginalMol.getBonds(); bond++) {	// mark all bonds to hydrogen
			int atom1 = mOriginalMol.getBondAtom(0, bond);
			int atom2 = mOriginalMol.getBondAtom(1, bond);
			if (mOriginalMol.getAtomicNo(atom1) == 1
			 || mOriginalMol.getAtomicNo(atom2) == 1)
				isHydrogenBond[bond] = true;
			}

		int lastNonHBond = mOriginalMol.getBonds();
		do lastNonHBond--; while ((lastNonHBond >= 0) && isHydrogenBond[lastNonHBond]);

		for (int bond=0; bond<lastNonHBond; bond++) {
			if (isHydrogenBond[bond]) {
				mOriginalMol.swapBonds(bond, lastNonHBond);
				isHydrogenBond[bond] = false;
				do lastNonHBond--;
				while (isHydrogenBond[lastNonHBond]);
				}
			}
		}

	private void addHydrogen(StereoMolecule mol, int atom, int mass) {
		int hydrogen = mol.addAtom(1);
		mol.setAtomMass(hydrogen, mass);
		mol.addBond(atom, hydrogen, Molecule.cBondTypeSingle);
		}

	public int getTautomerCount() {
		return mBondOrderSet.size();
		}

	private void countAndRemoveDAndT() {
		for (int bond=0; bond<mOriginalMol.getAllBonds(); bond++) {
			for (int i=0; i<2; i++) {
				int atom1 = mOriginalMol.getBondAtom(i, bond);
				int atom2 = mOriginalMol.getBondAtom(1-i, bond);
				if (mOriginalMol.getAtomicNo(atom1) == 1
				 && mOriginalMol.getAtomMass(atom1) > 1
				 && mOriginalMol.getAtomicNo(atom2) > 1
				 && mAtomRegionNo[atom2] != 0) {
					if (mOriginalMol.getAtomMass(atom1) == 2) {
						if (mAtomDCount == null)
							mAtomDCount = new int[mOriginalMol.getAllAtoms()];
						mAtomDCount[atom2]++;
						}
					else {
						if (mAtomTCount == null)
							mAtomTCount = new int[mOriginalMol.getAllAtoms()];
						mAtomTCount[atom2]++;
						}
					mOriginalMol.markAtomForDeletion(atom1);
					}
				}
			}
		if (mAtomDCount != null || mAtomTCount != null)
			mOriginalMol.deleteMarkedAtomsAndBonds();
		}

	/**
	 * @return number of tautomeric regions, i.e. atoms being connected by tautomeric bonds
	 */
	public int getAtomRegionCount() {
		return mRegionCount;
		}

	/**
	 * Returns region numbers for all atoms, i.e. numbers identifying the respective connected
	 * tautomeric regions the atom belongs to. Atoms sharing the same region share the same number.<br>
	 * 0: not member of a tautomer region; 1 and above: region number
	 * @return region count
	 */
	public int[] getAtomRegionNumbers() {
		return mAtomRegionNo;
		}

	/**
	 * Considers connected tautomer bonds to belong to a tautomer region.
	 * All independent tautomer regions are located and member atoms assigned to them.
	 * mAtomRegionNo[] is set accordingly.
	 * 0: not member of a tautomer region; 1 and above: region number
	 */
	private void assignRegionNumbers() {
		mAtomRegionNo = new int[mOriginalMol.getAtoms()];
		int[] graphAtom = new int[mOriginalMol.getAtoms()];
		boolean[] bondWasSeen = new boolean[mOriginalMol.getBonds()];
		int region = 0;

		for (int bond=0; bond<mOriginalMol.getBonds(); bond++) {
			if (!bondWasSeen[bond] && mIsTautomerBond[bond]) {
				region++;
				mAtomRegionNo[mOriginalMol.getBondAtom(0, bond)] = region;
				mAtomRegionNo[mOriginalMol.getBondAtom(1, bond)] = region;
				bondWasSeen[bond] = true;
				for (int i=0; i<2; i++) {
					int atom = mOriginalMol.getBondAtom(i, bond);
					mAtomRegionNo[atom] = region;
					int current = 0;
					int highest = 0;
					graphAtom[0] = atom;
					while (current <= highest) {
						for (int j=0; j<mOriginalMol.getConnAtoms(graphAtom[current]); j++) {
							int connBond = mOriginalMol.getConnBond(graphAtom[current], j);
							if (!bondWasSeen[connBond] && mIsTautomerBond[connBond]) {
								bondWasSeen[connBond] = true;
								int connAtom = mOriginalMol.getConnAtom(graphAtom[current], j);
								if (mAtomRegionNo[connAtom] == 0) {
									mAtomRegionNo[connAtom] = region;
									graphAtom[++highest] = connAtom;
									}
								}
							}
						current++;
						}
					}
				}
			}

		mRegionCount = region;
		}

	/**
	 * Counts for every region: pi-electrons, deuterium atoms, tritium atoms.
	 * Must be called after assignRegionNumbers().
	 */
	private void compileRegionCounts() {
		mRegionPiCount = new int[mRegionCount];
		mRegionDCount = new int[mRegionCount];
		mRegionTCount = new int[mRegionCount];

		for (int atom=0; atom<mOriginalMol.getAtoms(); atom++) {
			if (mAtomRegionNo[atom] != 0) {
				int regionIndex = mAtomRegionNo[atom]-1;
				if (mAtomDCount != null)
					mRegionDCount[regionIndex] += mAtomDCount[atom];
				if (mAtomTCount != null)
					mRegionTCount[regionIndex] += mAtomTCount[atom];
				}
			}
		for (int bond=0; bond<mOriginalMol.getBonds(); bond++) {
			if (mIsTautomerBond[bond]
			 && mOriginalMol.getBondOrder(bond) == 2) {
				mRegionPiCount[mAtomRegionNo[mOriginalMol.getBondAtom(0, bond)]-1] += 2;
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
	 * of pi-electrons and D and T atoms. The returned molecule has the fragment bit set.
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
		mol.setFragment(true);

		mol.ensureHelperArrays(Molecule.cHelperRings);

		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (mIsTautomerBond[bond]) {
				mol.setBondType(bond, Molecule.cBondTypeSingle);
				mol.setBondQueryFeature(bond, Molecule.cBondTypeSingle | Molecule.cBondTypeDouble, true);
				}
			}

		// remove stereo information from bonds that indicate a stereo center at one of the tautomer region atoms
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mAtomRegionNo[atom] != 0
			 && mOriginalMol.getNonHydrogenNeighbourCount(atom) < 4) {   // keep terternary S or P
				mol.convertStereoBondsToSingleBonds(atom);
				mol.setAtomConfigurationUnknown(atom,false);
				mol.setAtomESR(atom, Molecule.cESRTypeAbs, -1);
				}
			}

		// find highest ranking atom in every region
		int[] maxAtom = new int[mRegionCount];
		int[] maxRank = new int[mRegionCount];
		int[] atomRank = new Canonizer(mol).getFinalRank();
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mAtomRegionNo[atom] != 0) {
				int regionIndex = mAtomRegionNo[atom]-1;
				if (maxRank[regionIndex] < atomRank[atom]) {
					maxRank[regionIndex] = atomRank[atom];
					maxAtom[regionIndex] = atom;
					}
				}
			}

		// attach label with region counts to highest ranking atoms
		for (int i=0; i<mRegionCount; i++) {
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

			if (mBondOrderSet.size() >= sMaxTautomers) {
				if (!sSuppressWarning)
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
			 || mOriginalMol.getAtomicNo(atom) == 16
			 || mOriginalMol.getAtomicNo(atom) == 34
			 || mOriginalMol.getAtomicNo(atom) == 52);
		}

	private boolean isValidHeteroAtom(int atom) {
		return mHasFreeValence[atom]
			&& (mOriginalMol.getAtomicNo(atom) == 7
			 || mOriginalMol.getAtomicNo(atom) == 8
			 || mOriginalMol.getAtomicNo(atom) == 16
			 || mOriginalMol.getAtomicNo(atom) == 34
			 || mOriginalMol.getAtomicNo(atom) == 52);
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
				encoding[i >> 4] |= (Math.min(3, mol.getBondOrder(i)) << (2*(i & 15)));
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
				if (mIsTautomerBond[i]) {
					int bo = 3 & (encoding[i >> 4] >> (2*(i & 15)));
					tautomer.setBondType(i, bo == 1 ? Molecule.cBondTypeSingle
							: bo == 2 ? (mIsTautomerBond[i] && !mOriginalMol.isSmallRingBond(i) ? Molecule.cBondTypeCross : Molecule.cBondTypeDouble)
							: bo == 3 ? Molecule.cBondTypeTriple
							: Molecule.cBondTypeMetalLigand);
					}
				}
			}
		}
	}
