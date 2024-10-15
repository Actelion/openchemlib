package com.actelion.research.chem.shredder;

import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.SingularValueDecomposition;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;
import java.util.Arrays;

public class FragmentGeometry3D {
	public static final int MODE_SELECTED_ATOMS = 1;	// molecule with fragment atoms selected, exit vector atoms not selected
	public static final int MODE_FRAGMENT_WITH_EXIT_VECTORS = 2;	// defined as atom custom label "*"

	private final StereoMolecule mMol;
	private ExitVector[] mExitVector;
	private final String mFootPrint;	// canonical String describing atomic numbers or exit vectors
	private int[][] mPermutation;
	private final Coordinates[] mAlignmentCoords;
	private final Coordinates mAlignmentCOG;

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
		case MODE_SELECTED_ATOMS:
			initMoleculeWithSelection();
			break;
		case MODE_FRAGMENT_WITH_EXIT_VECTORS:
			initFragmentWithExitVectors();
			break;
		}

		Arrays.sort(mExitVector);

		StringBuilder footprint = new StringBuilder();
		for (ExitVector ev : mExitVector)
			footprint.append(Molecule.cAtomLabel[ev.atomicNo]);
		mFootPrint = footprint.toString();

		// compile coordinates of all root atoms and calculate their center of gravity
		mAlignmentCoords = determineAlignmentCoords();
		mAlignmentCOG = centerOfGravity(mAlignmentCoords);
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

	private Coordinates[] determineAlignmentCoords() {
		Coordinates[] coords = new Coordinates[2*mExitVector.length];
		for (int i=0; i<mExitVector.length; i++) {
			coords[i] = mMol.getCoordinates(mExitVector[i].rootAtom);
			coords[mExitVector.length + i] = mMol.getCoordinates(mExitVector[i].exitAtom);

			// for lonely hydrogens (single selected H connecting to non-selected exit atom)
			// we need to place the root coord (hydrogen) further away from the exit atom,
			// to reflect the longer bond length of a C-C compared to H-C
			if (mMol.getAtomicNo(mExitVector[i].rootAtom) == 1)
				coords[i] = (coords[mExitVector.length + i]).subC(coords[i]).scale(0.3).add(coords[i]);
		}

		return coords;
	}

	private void initFragmentWithExitVectors() {
		ArrayList<ExitVector> exitVectorList = new ArrayList<>();
		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			if ("*".equals(mMol.getAtomCustomLabel(atom)))
				exitVectorList.add(new ExitVector(mMol.getConnAtom(atom, 0), atom, mMol.getAtomicNo(atom)));

		mExitVector = exitVectorList.toArray(new ExitVector[0]);
	}

	public int getExitVectorCount() {
		return mExitVector.length;
	}

	public int getExitAtom(int i) {
		return mExitVector[i].exitAtom;
	}

	public int getRootAtom(int i) {
		return mExitVector[i].rootAtom;
	}

	public int getPermutedExitVectorIndex(int permutation, int i) {
		return mPermutation[permutation][i];
	}

	/**
	 * @param geometry
	 * @return whether it has the same number and kind of exit vectors
	 */
	public boolean equals(FragmentGeometry3D geometry) {
		return mFootPrint.equals(geometry.mFootPrint);
	}

	/**
	 * Aligns all root atoms of the passed geometry to this geometry. If the RMSD of all aligned
	 * root atoms is lower than maxRMSD, then the rotations matrix is returned. Otherwise, null
	 * is returned.
	 * @param geometry the geometry to be aligned with this
	 * @param permutation permutation index for equivalent root atoms of the passed geometry
	 * @param rmsdHolder null or double[1] to receive RMSD value
	 * @return rotation matrix or null, depending on whether alignment is acceptable
	 */
	public double[][] alignRootAndExitAtoms(FragmentGeometry3D geometry, int permutation, double[] rmsdHolder, double maxRMSD) {
		ExitVector[] geomEV = geometry.mExitVector;
		Coordinates[] coords = new Coordinates[2*geomEV.length];
		for (int i=0; i<geomEV.length; i++)
			coords[i] = new Coordinates(geometry.mAlignmentCoords[mPermutation[permutation][i]]);
		for (int i=0; i<geomEV.length; i++)
			coords[geomEV.length+i] = new Coordinates(geometry.mAlignmentCoords[geomEV.length+mPermutation[permutation][i]]);

		double[][] matrix = kabschAlign(mAlignmentCoords, coords, mAlignmentCOG, geometry.mAlignmentCOG);

		for (Coordinates c : coords) {
			c.sub(geometry.mAlignmentCOG);
			c.rotate(matrix);
			c.add(mAlignmentCOG);
		}

		double rmsd = Coordinates.getRmsd(mAlignmentCoords, coords);
		if (rmsdHolder != null)
			rmsdHolder[0] = rmsd;

		return rmsd > maxRMSD ? null : matrix;
	}

	public boolean hasMatchingExitVectors(FragmentGeometry3D geometry, Coordinates[] coords, int permutation, double[][] angleHolder, double maxAngleDivergence) {
		double[] angleDif = new double[mExitVector.length];
		boolean qualifies = true;
		for (int i=0; i<mExitVector.length; i++) {
			Coordinates v1 = mAlignmentCoords[mExitVector.length+i].subC(mAlignmentCoords[i]);
			ExitVector ev2 = geometry.mExitVector[mPermutation[permutation][i]];
			Coordinates v2 = coords[ev2.exitAtom].subC(coords[ev2.rootAtom]);
			angleDif[i] = v1.getAngle(v2) * 180 / Math.PI;
			if (angleDif[i] > maxAngleDivergence)
				qualifies = false;
		}
		if (angleHolder != null)
			angleHolder[0] = angleDif;
		return qualifies;
	}

	public int getPermutationCount() {
		if (mPermutation == null) {
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

			mPermutation = new int[totalPermCount][mExitVector.length];
			int multiplier = 1;
			for (int i=0; i<mExitVector.length; i++) {
				if (sameCount[i] <= 1) {
					for (int p=0; p<totalPermCount; p++)
						mPermutation[p][i] = i;
				}
				else {
					int[][] basePerm = createPermutations(sameCount[i]);
					int tp = 0;
					while (tp < totalPermCount) {
						for (int p=0; p<basePerm.length; p++) {
							for (int j=0; j<sameCount[i]; j++) {
								for (int m=0; m<multiplier; m++) {
									mPermutation[tp+m][i+j] = i+basePerm[p][j];
								}
							}
						}
						tp += multiplier;
					}
					multiplier *= permCount[i];
				}
			}
		}

		return mPermutation.length;
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

	public Coordinates getAlignmentCOG() {
		return mAlignmentCOG;
	}

	public static Coordinates centerOfGravity(Coordinates[] coords) {
		int counter = 0;
		Coordinates cog = new Coordinates();
		for(Coordinates c:coords) {
			cog.add(c);
			counter++;
		}
		cog.scale(1.0/counter);
		return cog;
	}

	/**
	 * Calculates the rotation matrix to rigidly align coords2 on coords1.
	 * coords1 and coords2 are not touched.<br>
	 * To actually perform the alignment with any set of coordinates do:<br>
	 * for (Coordinates c : anyCoords) { c.sub(cog2); c.rotate(matrix); c.add(cog1); }
	 * @param coords1
	 * @param coords2
	 * @param cog1
	 * @param cog2
	 * @return
	 */
	public static double[][] kabschAlign(Coordinates[] coords1, Coordinates[] coords2,
										 Coordinates cog1, Coordinates cog2) {
		double[][] m = new double[3][3];
		double[][] c1 = Arrays.stream(coords1).map(e -> new double[] {e.x-cog1.x,e.y-cog1.y,e.z-cog1.z}).toArray(double[][]::new);
		double[][] c2 = Arrays.stream(coords2).map(e -> new double[] {e.x-cog2.x,e.y-cog2.y,e.z-cog2.z}).toArray(double[][]::new);
		for(int i=0;i<3;i++) {
			for(int j=0;j<3;j++) {
				double rij = 0.0;
				for(int a=0; a<c1.length; a++)
					rij+= c2[a][i]* c1[a][j];
				m[i][j] = rij;
			}
		}

		SingularValueDecomposition svd = new SingularValueDecomposition(m,null,null);

		Matrix u = new Matrix(svd.getU());
		Matrix v = new Matrix(svd.getV());
		Matrix ut = u.getTranspose();
		Matrix vut = v.multiply(ut);
		double det = vut.det();

		Matrix ma = new Matrix(3,3);
		ma.set(0,0,1.0);
		ma.set(1,1,1.0);
		ma.set(2,2,det);

		Matrix rot = ma.multiply(ut);
		rot = v.multiply(rot);
		assert(rot.det()>0.0);
		rot = rot.getTranspose();
		return rot.getArray();
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
