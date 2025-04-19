package com.actelion.research.chem.shredder;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.*;
import com.actelion.research.util.IntArrayComparator;

import java.util.ArrayList;
import java.util.TreeSet;

public class Fragmenter3D {
	private final int mMinAtoms,mMaxAtoms,mMaxBondFlexibilitySum,mMinExits,mMaxExits;
	private final ArrayList<Fragment3D> mFragmentList;

	/**
	 * 
	 * @param minAtoms
	 * @param maxAtoms
	 * @param maxBondFlexibilitySum
	 * @param minExits
	 * @param maxExits
	 */
	public Fragmenter3D(int minAtoms, int maxAtoms, int maxBondFlexibilitySum, int minExits, int maxExits) {
		mMinAtoms = minAtoms;
		mMaxAtoms = maxAtoms;
		mMaxBondFlexibilitySum = maxBondFlexibilitySum;
		mMinExits = minExits;
		mMaxExits = maxExits;
		mFragmentList = new ArrayList<>();
		}

	/**
	 * Applying the constraints passed to the Fragmenter3D constructor, this method shredders
	 * the given 3D-molecule and returns all generated 3D-fragments as an ArrayList.
	 * The list is re-used by subsequent calls to this nethod. Thus, process/consume the
	 * fragment list before calling this method again.
	 * @param mol
	 * @param withHydrogen whether built fragments shall include hydrogen atoms
	 * @return Fragment3D list of given molecule
	 */
	public ArrayList<Fragment3D> buildFragments(StereoMolecule mol, boolean withHydrogen) {
		mFragmentList.clear();

		mol.stripSmallFragments();
		if (withHydrogen)
			new AtomAssembler(mol).addImplicitHydrogens();
		else
			mol.removeExplicitHydrogens(false);

		boolean[] isRotatableBond = new boolean[mol.getAllBonds()];
		TorsionDB.findRotatableBonds(mol, true, isRotatableBond);
		float[] bondFlexibility = new MolecularFlexibilityCalculator().calculateBondFlexibilities(mol, isRotatableBond);

		int[] fragmentNo = new int[mol.getAllAtoms()];
		int fragmentCount = mol.getFragmentNumbers(fragmentNo, isRotatableBond, true);

		int[] atomCount = new int[fragmentCount];
		for (int atom=0; atom<mol.getAtoms(); atom++)
			atomCount[fragmentNo[atom]]++;

		BaseFragmentInfo[] fragmentData = new BaseFragmentInfo[fragmentCount];
		for (int i=0; i<fragmentCount; i++)
			fragmentData[i] = new BaseFragmentInfo(mol, fragmentNo, i, atomCount[i]);

		TreeSet<int[]> baseFragmentCombinationSet = new TreeSet<>(new IntArrayComparator());
		boolean[] isMemberFragment = new boolean[fragmentCount];
		for (int i=0; i<fragmentCount; i++) {
			isMemberFragment[i] = true;

			addNewFragments(mol, fragmentNo, bondFlexibility, 0f, isMemberFragment, i, 1,
					fragmentData[i].atomCount, baseFragmentCombinationSet, fragmentData);

			isMemberFragment[i] = false;
			}

		return mFragmentList;
		}

	/**
	 * Recursively adds one baseFragment forming a larger joined Fragment3D, which is added to the set,
	 * provided that the size criteria are fullfilled.
	 * @param mol
	 */
	private void addNewFragments(StereoMolecule mol, int[] fragmentNo, float[] bondFlexibility, float bondFlexibilitySum,
	                             boolean[] isMemberFragment, int previousBaseFragment, int usedBaseFragmentCount, int atomCount,
	                             TreeSet<int[]> baseFragmentCombinationSet, BaseFragmentInfo[] baseFragmentInfo) {
		if (atomCount > mMaxAtoms)
			return;

		int[] baseFragmentList = new int[usedBaseFragmentCount];
		int index = 0;
		for (int i=0; i<isMemberFragment.length; i++)
			if (isMemberFragment[i])
				baseFragmentList[index++] = i;

		if (baseFragmentCombinationSet.contains(baseFragmentList))
			return;

		baseFragmentCombinationSet.add(baseFragmentList);

		// build the Fragment3D from the base fragment members and add it to the queue
		if (atomCount >= mMinAtoms)
			addFragment(mol, fragmentNo, isMemberFragment);

		BaseFragmentInfo previousBaseFragmentInfo = baseFragmentInfo[previousBaseFragment];
		for (int i=0; i<previousBaseFragmentInfo.neighbourCount; i++) {
			int neighbour = previousBaseFragmentInfo.neighbourFragment[i];
			int neighbourBond = previousBaseFragmentInfo.neighbourBond[i];
			if (!isMemberFragment[neighbour]
			 && bondFlexibilitySum + bondFlexibility[neighbourBond] <= mMaxBondFlexibilitySum) {
				isMemberFragment[neighbour] = true;

				addNewFragments(mol, fragmentNo, bondFlexibility, bondFlexibilitySum + bondFlexibility[neighbourBond],
						isMemberFragment, neighbour, usedBaseFragmentCount+1,
						atomCount+baseFragmentInfo[neighbour].atomCount, baseFragmentCombinationSet, baseFragmentInfo);

				isMemberFragment[neighbour] = false;
				}
			}
		}

	private void addFragment(StereoMolecule mol, int[] fragmentNo, boolean[] isMemberFragment) {
		int atomCount = 0;
		int extendedAtomCount = 0;

		// mark all atoms with specified fragmentNo and two layers around it
		boolean[] includeAtom = new boolean[mol.getAllAtoms()];
		boolean[] isExitAtom = new boolean[mol.getAllAtoms()];
		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			if (isMemberFragment[fragmentNo[atom]]) {
				includeAtom[atom] = true;
				atomCount++;
				for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
					int connAtom = mol.getConnAtom(atom, i);
					if (!isMemberFragment[fragmentNo[connAtom]]) {
						includeAtom[connAtom] = true;
						isExitAtom[connAtom] = true;
						atomCount++;
						extendedAtomCount++;
						}
					}
				}
			}

		if (extendedAtomCount < mMinExits || extendedAtomCount > mMaxExits)
			return;

		int bondCount = 0;
		for (int bond=0; bond<mol.getAllBonds(); bond++)
			if (includeAtom[mol.getBondAtom(0, bond)]
			 && includeAtom[mol.getBondAtom(1, bond)])
				bondCount++;

		int[] atomMap = new int[mol.getAllAtoms()];

		StereoMolecule fragment = new StereoMolecule(atomCount, bondCount);
		mol.copyMoleculeByAtoms(fragment, includeAtom, false, atomMap);

		fragment.setFragment(false);
		fragment.center();

		int[] exitAtom = new int[extendedAtomCount];
		int exitAtomIndex = 0;
		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			if (isExitAtom[atom]) {
				fragment.setAtomCustomLabel(atomMap[atom],"*");
				exitAtom[exitAtomIndex++] = atomMap[atom];
				}
			}

		Canonizer canonizer = new Canonizer(fragment, Canonizer.ENCODE_ATOM_CUSTOM_LABELS);
		String idcode = canonizer.getIDCode();
		String coords = canonizer.getEncodedCoordinates();
		TorsionDescriptor td = new TorsionDescriptorHelper(canonizer.getCanMolecule()).getTorsionDescriptor();

		for (int i=0; i<exitAtom.length; i++)
			exitAtom[i] = canonizer.getGraphIndexes()[exitAtom[i]];

		mFragmentList.add(new Fragment3D(idcode, coords, td, exitAtom));
		}
	}

class BaseFragmentInfo {
	int[] neighbourFragment,neighbourAtom,neighbourBond;
	int neighbourCount,atomCount;

	public BaseFragmentInfo(StereoMolecule mol, int[] fragmentNo, int fragmentIndex, int atomCount) {
		neighbourCount = 0;
		for (int bond=0; bond<mol.getBonds(); bond++)
			if (fragmentNo[mol.getBondAtom(0, bond)] == fragmentIndex
					^ fragmentNo[mol.getBondAtom(1, bond)] == fragmentIndex)
				neighbourCount++;

		neighbourFragment = new int[neighbourCount];
		neighbourAtom = new int[neighbourCount];
		neighbourBond = new int[neighbourCount];

		int neighbourIndex = 0;
		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (fragmentNo[mol.getBondAtom(0, bond)] == fragmentIndex
			  ^ fragmentNo[mol.getBondAtom(1, bond)] == fragmentIndex) {
				neighbourAtom[neighbourIndex] = mol.getBondAtom(fragmentNo[mol.getBondAtom(0, bond)] == fragmentIndex ? 1 : 0, bond);
				neighbourBond[neighbourIndex] = bond;
				neighbourFragment[neighbourIndex] = fragmentNo[neighbourAtom[neighbourIndex]];
				neighbourIndex++;
				}
			}

		this.atomCount = atomCount;
		}
	}