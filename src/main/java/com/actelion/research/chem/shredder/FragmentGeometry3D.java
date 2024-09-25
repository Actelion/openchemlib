package com.actelion.research.chem.shredder;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;
import java.util.Arrays;

public class FragmentGeometry3D {
	public static final int MODE_SELECTED_ATOMS = 1;	// molecule with fragment atoms selected, exit vector atoms not selected
	public static final int MODE_FRAGMENT_WITH_EXIT_VECTORS = 2;	// defined as atom custom label "*"

	private StereoMolecule mMol;
	private ExitVector[] mExitVector;
	private String mFootPrint;	// canonical String describing atomic numbers or exit vectors
	private int[][] mPermutations;
	private Coordinates[] mCoordinates;

	/**
	 * Creates a FragmentGeometry3D from a StereoMolecule with the mode defining the situation.
	 * This creates a canonical exit vector footprint. For geometries with the same footprint
	 * it provides geometry comparisons regarding alignment and exit vector similarity.
	 * A transformation mask is provided for the best alignment one object to another.
	 * Helper functions to attach properly aligned substituents are given.
	 * @param mol
	 */
	public FragmentGeometry3D(StereoMolecule mol, int mode) {
		mMol = mol;
		mMol.ensureHelperArrays(Molecule.cHelperNeighbours);

		switch (mode) {
			case MODE_SELECTED_ATOMS: initMoleculeWithSelection();
			case MODE_FRAGMENT_WITH_EXIT_VECTORS: initFragmentWithExitVectors();
		}

		Arrays.sort(mExitVector);

		StringBuilder footprint = new StringBuilder();
		for (ExitVector ev : mExitVector)
			footprint.append(ev.atomicNo);
		mFootPrint = footprint.toString();
	}

	public void initMoleculeWithSelection() {
		ArrayList<ExitVector> exitVectorList = new ArrayList<>();
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			if (mMol.isSelectedAtom(atom)) {
				for (int i=0; i<mMol.getAllConnAtoms(atom); i++) {
					int connAtom = mMol.getConnAtom(atom, i);
					if (!mMol.isSelectedAtom(connAtom))
						exitVectorList.add(new ExitVector(atom, connAtom, mMol.getAtomicNo(connAtom)));
				}
			}
		}

		mExitVector = exitVectorList.toArray(new ExitVector[0]);
	}

	private void initFragmentWithExitVectors() {
		ArrayList<ExitVector> exitVectorList = new ArrayList<>();
		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			if ("*".equals(mMol.getAtomCustomLabel(atom)))
				exitVectorList.add(new ExitVector(mMol.getConnAtom(atom, 0), atom, mMol.getAtomicNo(atom)));

		mExitVector = exitVectorList.toArray(new ExitVector[0]);
	}

	/**
	 * @param geometry
	 * @return whether it has the same number and kind of exit vectors
	 */
	public boolean equals(FragmentGeometry3D geometry) {
		return mFootPrint.equals(geometry.mFootPrint);
	}

	public float align(FragmentGeometry3D geometry, int[] permutation) {
		if (mCoordinates == null) {
			mCoordinates = new Coordinates[mExitVector.length];
			for (int i=0; i<mExitVector.length; i++)
				mCoordinates[i] = mMol.getCoordinates(mExitVector[i].rootAtom);
		}

		Coordinates[] coords = new Coordinates[mExitVector.length];
		for (int i=0; i<mExitVector.length; i++)
			coords[i] = geometry.mMol.getCoordinates(mExitVector[permutation[i]].rootAtom);

		// TODO align...
		return 0;
	}

	public int[][] getPermutations() {
		if (mPermutations == null) {
			int[] sameCount = new int[mExitVector.length];
			int[] permCount = new int[mExitVector.length];
			int index = 0;
			int totalPermCount = 1;
			while (index<mExitVector.length) {
				sameCount[index] = 1;
				permCount[index] = 1;
				while (index+sameCount[index]<mExitVector.length
					&& mExitVector[index].compareTo(mExitVector[index+sameCount[index]]) == 0) {
					sameCount[index]++;
					permCount[index] *= sameCount[index];
				}
				totalPermCount *= sameCount[index];
				index += sameCount[index];
			}

			mPermutations = new int[totalPermCount][mExitVector.length];
			int multiplier = 1;
			for (int i=0; i<mExitVector.length; i++) {
				if (sameCount[i] <= 1) {
					for (int p=0; p<totalPermCount; p++)
						mPermutations[p][i] = i;
				}
				else {
					int[][] basePerm = createPermutations(sameCount[i]);
					int tp = 0;
					while (tp < totalPermCount) {
						for (int p=0; p<basePerm.length; p++) {
							for (int j=0; j<sameCount[i]; j++) {
								for (int m=0; m<multiplier; m++) {
									mPermutations[tp+m][i+j] = i+basePerm[p][j];
								}
							}
						}
						tp += multiplier;
					}
					multiplier *= permCount[i];
				}
			}
		}

		return mPermutations;
	}

	private int[][] createPermutations(int objectCount) {
		ArrayList<int[]> list = new ArrayList<>();
		addToPermutation(new int[1], objectCount, list);
		return list.toArray(new int[0][]);
	}

	private void addToPermutation(int[] input, int max, ArrayList<int[]> list) {
		int size = input.length + 1;
		for (int pos=0; pos<size; pos++) {
			int[] output = new int[size];
			int source = 0;
			int target = 0;
			while (target < size) {
				if (target == pos)
					output[target] = input.length;
				else
					output[target] = input[source++];
				target++;
			}

			if (size == max)
				list.add(output);
			else
				addToPermutation(output, max, list);
		}
	}

	private static class ExitVector implements Comparable<ExitVector> {
		int rootAtom,exitAtom,atomicNo;

		public ExitVector(int rootAtom, int exitAtom, int atomicNo) {
			this.rootAtom = rootAtom;
			this.exitAtom = exitAtom;
			this.atomicNo = atomicNo;
		}

		@Override
		public int compareTo(ExitVector o) {
			return Integer.compare(o.atomicNo, atomicNo);
		}
	}
}
