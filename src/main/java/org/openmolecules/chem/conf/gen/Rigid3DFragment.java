/*
 * @(#)Rigid3DFragment.java
 *
 * Copyright 2013 openmolecules.org, Inc. All Rights Reserved.
 *
 * NOTICE: All information contained herein is, and remains the property
 * of openmolecules.org.  The intellectual and technical concepts contained
 * herein are proprietary to openmolecules.org.
 * Actelion Pharmaceuticals Ltd. is granted a non-exclusive, non-transferable
 * and timely unlimited usage license.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.gen;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import org.openmolecules.chem.conf.so.ConformationSelfOrganizer;
import org.openmolecules.chem.conf.so.SelfOrganizedConformer;

import java.util.ArrayList;

public class Rigid3DFragment {
	private static final int MAX_CONFORMERS = 16;
	private static long sRandomSeed = 0;	 // no specific seed

	private int mCoreAtomCount;
	private int[] mFragmentToOriginalAtom,mCoreToFragmentAtom,
				  mExtendedToFragmentAtom,mOriginalToExtendedAtom;
	private double[] mConformerLikelyhood;
	private ConformationSelfOrganizer mSelfOrganizer;
	private SelfOrganizedConformer[] mConformerList;

private StereoMolecule[] mFragment;	// TODO remove
public StereoMolecule[] getFragment() { return mFragment; };

	public static void setRandomSeed(long seed) {
		sRandomSeed = seed;
		}

	/**
	 * Creates a conformer of a part of the passed molecule that is supposed to
	 * be a rigid fragment. Atoms are considered part of the fragment, if they have
	 * the specified fragmentIndex. Member atoms that have non-member neighbors
	 * will have incomplete valences in the rigid fragment.
	 * Helper arrays must provide valid parities, which are considered when
	 * generating 3D coordinates for all member atoms.
	 * @param mol
	 * @param fragmentNo
	 * @param fragmentIndex
	 */
	public Rigid3DFragment(StereoMolecule mol, int[] fragmentNo, int fragmentIndex) {
		mCoreAtomCount = 0;
		int atomCount = 0;
		int extendedAtomCount = 0;

		// mark all atoms with specified fragmentNo and two layers around it
		boolean[] includeAtom = new boolean[mol.getAllAtoms()];
		boolean[] isOuterShellAtom = new boolean[mol.getAllAtoms()];
		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			if (fragmentNo[atom] == fragmentIndex) {
				includeAtom[atom] = true;
				atomCount++;
				mCoreAtomCount++;
				for (int i=0; i<mol.getConnAtoms(atom); i++) {
					int connAtom = mol.getConnAtom(atom, i);
					if (fragmentNo[connAtom] != fragmentIndex) {
						includeAtom[connAtom] = true;
						atomCount++;
						extendedAtomCount++;
						for (int j=0; j<mol.getAllConnAtoms(connAtom); j++) {
							int connconn = mol.getConnAtom(connAtom, j);
							if (fragmentNo[connconn] != fragmentIndex) {
								isOuterShellAtom[connconn] = true;
								includeAtom[connconn] = true;
								atomCount++;
								}
							}
						}
					}
				}
			}

		int bondCount = 0;
		for (int bond=0; bond<mol.getAllBonds(); bond++)
			if (includeAtom[mol.getBondAtom(0, bond)]
			 && includeAtom[mol.getBondAtom(1, bond)])
				bondCount++;

		StereoMolecule fragment = new StereoMolecule(atomCount, bondCount);
		mol.copyMoleculeByAtoms(fragment, includeAtom, false, null);
		fragment.setFragment(false); // if mFragment is a fragment, then H-atoms are converted to query features!!!
		fragment.setParitiesValid(0);

		mCoreToFragmentAtom = new int[mCoreAtomCount];
		mFragmentToOriginalAtom = new int[atomCount];
		mExtendedToFragmentAtom = new int[mCoreAtomCount+extendedAtomCount];
		mOriginalToExtendedAtom = new int[mol.getAllAtoms()];
		int coreAtom = 0;
		int fragmentAtom = 0;
		int extendedAtom = 0;
		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			if (includeAtom[atom]) {
// this would convert the molecule to a fragment
//				if (mol.isFlatNitrogen(atom))
//				    fragment.setAtomQueryFeature(fragmentAtom, Molecule.cAtomQFFlatNitrogen, true);

				if (isOuterShellAtom[atom])
					fragment.setAtomMarker(fragmentAtom, true);

				if (fragmentNo[atom] == fragmentIndex || !isOuterShellAtom[atom]) {
					mExtendedToFragmentAtom[extendedAtom] = fragmentAtom;
					mOriginalToExtendedAtom[atom] = extendedAtom;
					extendedAtom++;
					}
				if (fragmentNo[atom] == fragmentIndex) {
					mCoreToFragmentAtom[coreAtom] = fragmentAtom;
					coreAtom++;
					}
				mFragmentToOriginalAtom[fragmentAtom] = atom;
				fragmentAtom++;
				}
			}

		mSelfOrganizer = new ConformationSelfOrganizer(fragment, true);
		mSelfOrganizer.initializeConformers(sRandomSeed, MAX_CONFORMERS);

		// Generate multiple low constrain conformers
		ArrayList<SelfOrganizedConformer> conformerList = new ArrayList<SelfOrganizedConformer>();
		SelfOrganizedConformer bestConformer = mSelfOrganizer.getNextConformer();
		conformerList.add(bestConformer);
		SelfOrganizedConformer conformer = mSelfOrganizer.getNextConformer();
		while (conformer != null) {
			conformerList.add(conformer);
			conformer = mSelfOrganizer.getNextConformer();
			}
		mConformerList = conformerList.toArray(new SelfOrganizedConformer[0]);

		mConformerLikelyhood = new double[mConformerList.length];
		double likelyhoodSum = 0.0;
		for (int i=0; i<mConformerList.length; i++) {
			mConformerLikelyhood[i] = mConformerList[i].getLikelyhood();
			likelyhoodSum += mConformerLikelyhood[i];
			}
		if (likelyhoodSum != 0.0) {
			for (int i=0; i<mConformerList.length; i++)
				mConformerLikelyhood[i] /= likelyhoodSum;
			}

if (ConformerGenerator.WRITE_DW_FRAGMENT_FILE) {
 mFragment = new StereoMolecule[mConformerList.length];
 for (int i=0; i<mConformerList.length; i++) {
  mFragment[i] = fragment.getCompactCopy();
  mConformerList[i].toMolecule(mFragment[i]);
 }}
	    }

	public int getConformerCount() {
		return mConformerList.length;
		}

	public double getConformerLikelyhood(int i) {
		return mConformerLikelyhood[i];
		}

	/**
	 * Calculates a random conformer index giving conformers with lower strains
	 * a higher chance to be selected. With a progress value of 0.0 selection
	 * likelyhoods are proportional to conformer likelyhoods due to lower strains.
	 * With increasing progress value the higher frequent conformers get less
	 * and less preferred until 1.0 without any preference.
	 * @param random
	 * @param progress 0...1 
	 */
	public int getLikelyRandomConformerIndex(double random, double progress) {
		double sum = 0;
		for (int t=0; t<mConformerLikelyhood.length; t++) {
			double contribution = (1f-progress)*mConformerLikelyhood[t] + progress/mConformerLikelyhood.length;
			sum += contribution;
			if (random <= sum)
				return t;
			}
		return mConformerLikelyhood.length-1;  // should never reach this
		}

	public int originalToExtendedAtom(int originalAtom) {
		return mOriginalToExtendedAtom[originalAtom];
		}

	public int coreToOriginalAtom(int atom) {
		return mFragmentToOriginalAtom[mCoreToFragmentAtom[atom]];
		}

	public int extendedToOriginalAtom(int atom) {
		return mFragmentToOriginalAtom[mExtendedToFragmentAtom[atom]];
		}

	public int getConnectionPointCount() {
		return mExtendedToFragmentAtom.length - mCoreToFragmentAtom.length;
		}

	/**
	 * @return count of core atoms, i.e. atoms inside of rotatable bonds
	 */
	public int getCoreSize() {
		return mCoreToFragmentAtom.length;
		}

	/**
	 * @return count of core and 1st shell atoms, i.e. core and rotatable bond atoms
	 */
	public int getExtendedSize() {
		return mExtendedToFragmentAtom.length;
		}

	public Coordinates getCoreCoordinates(int conformer, int atom) {
		return mConformerList[conformer].getCoordinates(mCoreToFragmentAtom[atom]);
		}

	public Coordinates getExtendedCoordinates(int conformer, int atom) {
		return mConformerList[conformer].getCoordinates(mExtendedToFragmentAtom[atom]);
		}
	}
