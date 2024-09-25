/*
 * @(#)MoleculeShredder.java
 *
 * Copyright 1997-2001 Actelion Ltd., Inc. All Rights Reserved.
 *
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.
 *
 * @author Thomas Sander
 */

package com.actelion.research.chem.shredder;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.*;

public class MoleculeShredder {
	public static final int cModeSetAtomQFNoMoreNeighbours = 1;
	public static final int cModeExcludeCoreFragments = 64;

		// mode may contain one or none of those:
	public static final int cModeRetainSubstituentTypes = 2;
	public static final int cModeDiscardSubstitutionInfo = 128;

		// mode may contain one or none of these:
	public static final int cModeCoreFragmentsOnly = 4;
	public static final int cModeOneSpereAroundCoreOnly = 32;

		// mode may contain one or none of the 'Cut' options:
	public static final int cModeCutTerminalBondsAlso = 8;
	public static final int cModeCutSpacersOffRings = 16;

	public static final int cColorFragmentsInMolecule = 256;

	private static final int cMaxCoreFragments = 8;

	private final StereoMolecule mMol;
	private boolean[]			mIsCuttableBond;
	private boolean[][]			mAreNeighbours;
	private final int			mMode;
	private int					mCuttableBonds;
	private int					mCoreFragments;
	private int[]				mCoreFragmentNo;
	private ArrayList<Long>		mFragmentList;

	public MoleculeShredder(StereoMolecule mol, int mode) {
		mMol = mol;
		mMol.ensureHelperArrays(Molecule.cHelperRings);
		mMode = mode;
		locateCuttableBonds();
		locateCoreFragments();
		if ((mode & cModeCoreFragmentsOnly) == 0)
			createCoreNeighbourMatrix();
		createFragmentList();
		}


	public int getFragments() {
		return mFragmentList.size();
		}


	public int getCoreFragments() {
		return mCoreFragments;
		}


	public StereoMolecule getFragment(StereoMolecule theFragment, int fragmentNo) {
		boolean[] isMember = new boolean[mMol.getAtoms()];
		long fragmentBits = mFragmentList.get(fragmentNo);
		for (int coreFragment=0; coreFragment<mCoreFragments; coreFragment++)
			if ((fragmentBits & ((long)1 << coreFragment)) != 0)
				for (int atom=0; atom<mMol.getAtoms(); atom++)
					if (mCoreFragmentNo[atom] == coreFragment)
						isMember[atom] = true;

		if (theFragment == null) {
			int atomCount = 0;
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (isMember[atom])
					atomCount++;

			int bondCount = 0;
			for (int sourceBond=0; sourceBond<mMol.getBonds(); sourceBond++) {
				int atom1 = mMol.getBondAtom(0, sourceBond);
				int atom2 = mMol.getBondAtom(1, sourceBond);
				if (isMember[atom1] && isMember[atom2])
					bondCount++;
				else if ((isMember[atom1] || isMember[atom2])
					  && (mMode & cModeRetainSubstituentTypes) != 0)
					bondCount++;
				}

			theFragment = new StereoMolecule(atomCount, bondCount);
			theFragment.setFragment(true);
			}
		else {
			theFragment.clear();
			theFragment.setFragment(true);
			}

		int[] newAtom = new int[mMol.getAtoms()];
		for (int sourceAtom=0; sourceAtom<mMol.getAtoms(); sourceAtom++) {
			if (isMember[sourceAtom])
				newAtom[sourceAtom] = mMol.copyAtom(theFragment, sourceAtom, 0, 0);

			if ((mMode & cColorFragmentsInMolecule) != 0) {
				if (isMember[sourceAtom])
					mMol.setAtomColor(sourceAtom, Molecule.cAtomColorDarkRed);
				else
					mMol.setAtomColor(sourceAtom, Molecule.cAtomColorNone);
				}
			}

		if ((mMode & cModeSetAtomQFNoMoreNeighbours) != 0)
			for (int atom=0; atom<theFragment.getAllAtoms(); atom++)
				theFragment.setAtomQueryFeature(atom, Molecule.cAtomQFNoMoreNeighbours, true);

		for (int sourceBond=0; sourceBond<mMol.getBonds(); sourceBond++) {
			int atom1 = mMol.getBondAtom(0, sourceBond);
			int atom2 = mMol.getBondAtom(1, sourceBond);
			if (isMember[atom1] && isMember[atom2])
			    mMol.copyBond(theFragment, sourceBond, 0, 0, newAtom, false);
			else if (mIsCuttableBond[sourceBond] && (mMode & cModeDiscardSubstitutionInfo) == 0) {
				if (isMember[atom1]) {
					if ((mMode & cModeRetainSubstituentTypes) == 0) {
						theFragment.setAtomQueryFeature(newAtom[atom1], Molecule.cAtomQFNoMoreNeighbours, false);
						theFragment.setAtomQueryFeature(newAtom[atom1], Molecule.cAtomQFMoreNeighbours, true);
						}
					else {	// distinguish between carbon and non-carbon substituents
						int substituent = theFragment.addAtom(6);
						theFragment.addBond(newAtom[atom1], substituent, mMol.getBondType(sourceBond));
						if (mMol.getAtomicNo(atom2) != 6) {	// set atom list to exclude carbon
                            int[] list = new int[1];
							list[0] = 6;
							theFragment.setAtomList(substituent, list, true);
							}
						}
					}
				if (isMember[atom2]) {
					if ((mMode & cModeRetainSubstituentTypes) == 0) {
						theFragment.setAtomQueryFeature(newAtom[atom2], Molecule.cAtomQFNoMoreNeighbours, false);
						theFragment.setAtomQueryFeature(newAtom[atom2], Molecule.cAtomQFMoreNeighbours, true);
						}
					else {
						int substituent = theFragment.addAtom(6);
						theFragment.addBond(newAtom[atom2], substituent, mMol.getBondType(sourceBond));
						if (mMol.getAtomicNo(atom1) != 6) {
                            int[] list = new int[1];
                            list[0] = 6;
                            theFragment.setAtomList(substituent, list, true);
							}
						}
					}
				}
			}

		return theFragment;
		}


	private void locateCuttableBonds() {
		mCuttableBonds = 0;
		mIsCuttableBond = new boolean[mMol.getBonds()];
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if ((mMode & cModeCutTerminalBondsAlso) != 0) {
				// cut all non-ring-single-bonds
				if (!mMol.isRingBond(bond)
				 && mMol.getBondOrder(bond) == 1) {
					mIsCuttableBond[bond] = true;
					mCuttableBonds++;
					}
				}
			else if ((mMode & cModeCutSpacersOffRings) != 0) {
				// cut all exocyclic bonds
				if (!mMol.isRingBond(bond)
				 && (mMol.isRingAtom(mMol.getBondAtom(0, bond))
				  || mMol.isRingAtom(mMol.getBondAtom(1, bond)))) {
					mIsCuttableBond[bond] = true;
					mCuttableBonds++;
					}
				}
			else {
				// cut all non-terminal-non-ring-single-bonds
				if (!mMol.isRingBond(bond)
				 && mMol.getBondOrder(bond) == 1
				 && ((mMol.getConnAtoms(mMol.getBondAtom(0, bond)) > 1
				   && mMol.getConnAtoms(mMol.getBondAtom(1, bond)) > 1))) {
					mIsCuttableBond[bond] = true;
					mCuttableBonds++;
					}
				}
			}
		}


	private void locateCoreFragments() {
		mCoreFragmentNo = new int[mMol.getAtoms()];
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			mCoreFragmentNo[atom] = -1;

		mCoreFragments = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mCoreFragmentNo[atom] == -1) {
				mCoreFragmentNo[atom] = mCoreFragments;
				int[] graphAtom = new int[mMol.getAtoms()];
				graphAtom[0] = atom;
				int current = 0;
				int highest = 0;
				while (current <= highest) {
					for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
						if (!mIsCuttableBond[mMol.getConnBond(graphAtom[current], i)]) {
							int candidate = mMol.getConnAtom(graphAtom[current], i);
							if (mCoreFragmentNo[candidate] == -1) {
								graphAtom[++highest] = candidate;
								mCoreFragmentNo[candidate] = mCoreFragments;
								}
							}
						}
					current++;
					}
				mCoreFragments++;
				}
			}
		}


	private void createCoreNeighbourMatrix() {
		mAreNeighbours = new boolean[mCoreFragments][mCoreFragments];
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mIsCuttableBond[bond]) {
				int f1 = mCoreFragmentNo[mMol.getBondAtom(0, bond)];
				int f2 = mCoreFragmentNo[mMol.getBondAtom(1, bond)];
				mAreNeighbours[f1][f2] = true;
				mAreNeighbours[f2][f1] = true;
				}
			}
		}


	private void createFragmentList() {
		mFragmentList = new ArrayList<Long>();

		if (mCoreFragments == 1)	// if the only fragment is the entire molecule then return no fragments
			return;

		if (mCoreFragments > 64)	// fragment combinations cannot be coded in long -> don't shredder
			return;

		if ((mMode & cModeCoreFragmentsOnly) == 0
		 && mCoreFragments > cMaxCoreFragments)	// do not create fragments of to flexible molecules
			return;

		for (int fragment=0; fragment<mCoreFragments; fragment++)
			mFragmentList.add((long)1 << fragment);

		if ((mMode & cModeCoreFragmentsOnly) != 0)
			return;

		if ((mMode & cModeOneSpereAroundCoreOnly) != 0)
			createOneSphereFragmentList();
		else
			createFullPermutationFragmentList();

		if ((mMode & cModeExcludeCoreFragments) != 0)
			for (int fragment=0; fragment<mCoreFragments; fragment++)
				mFragmentList.remove(0);
		}

	private void createFullPermutationFragmentList() {
		int current = 0;
		for (int fragmentSize=2; fragmentSize<mCoreFragments; fragmentSize++) {
			int firstOfThisSize = mFragmentList.size();
			while (current < firstOfThisSize) {
				long currentFragment = mFragmentList.get(current);
				for (int coreFragment=0; coreFragment<mCoreFragments; coreFragment++) {
					if ((currentFragment & ((long)1 << coreFragment)) != 0) {
						for (int coreNeighbour=0; coreNeighbour<mCoreFragments; coreNeighbour++) {
							if (mAreNeighbours[coreFragment][coreNeighbour]
							 && (currentFragment & ((long)1 << coreNeighbour)) == 0) {
								long proposedFragment = currentFragment | ((long)1 << coreNeighbour);
								boolean found = false;
								for (int i=firstOfThisSize; i<mFragmentList.size(); i++) {
									if (proposedFragment == mFragmentList.get(i)) {
										found = true;
										break;
										}
									}

								if (!found)
									mFragmentList.add(proposedFragment);
								}
							}
						}
					}
				current++;
				}
			}
		}

	private void createOneSphereFragmentList() {
		int current = 0;

		ArrayList<Integer> centralCoreList = new ArrayList<Integer>();
		for (int fragment=0; fragment<mCoreFragments; fragment++)
			centralCoreList.add(fragment);

		for (int fragmentSize=2; fragmentSize<mCoreFragments; fragmentSize++) {
			int firstOfThisSize = mFragmentList.size();
			while (current < firstOfThisSize) {
				long currentFragment = mFragmentList.get(current);
				int coreFragment = centralCoreList.get(current);
				for (int coreNeighbour=0; coreNeighbour<mCoreFragments; coreNeighbour++) {
					if (mAreNeighbours[coreFragment][coreNeighbour]
					 && (currentFragment & ((long)1 << coreNeighbour)) == 0) {
						long proposedFragment = currentFragment | ((long)1 << coreNeighbour);
						boolean found = false;
						for (int i=firstOfThisSize; i<mFragmentList.size(); i++) {
							if (proposedFragment == mFragmentList.get(i)) {
								found = true;
								break;
								}
							}

						if (!found) {
							mFragmentList.add(proposedFragment);
							centralCoreList.add(centralCoreList.get(current));
							}
						}
					}
				current++;
				}
			}
		}
	}
