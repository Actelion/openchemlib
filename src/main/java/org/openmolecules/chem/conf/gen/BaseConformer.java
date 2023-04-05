package org.openmolecules.chem.conf.gen;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.TorsionDB;

import java.util.ArrayList;
import java.util.Random;

/**
 * A BaseConformer is an extended Conformer that uses a specific set of rigid fragment conformations.
 * A BaseConformer knows all rotatable bonds, which are those bonds connecting the rigid fragments.
 * For every rotatable bond, the BaseConformer calculates a set of possible torsion angles and their
 * likelyhoods considering the specific conformations of the pair of connected rigid fragments.
 * For every rotatable bond the number and order of torsion angles matches the one of the TorsionDB for
 * matching bond environments, however, torsion values and their likelyhoods may be adjusted to reflect
 * the particular sterical situation of the fragments being connected.
 * Typically, multiple Conformer instances are derived from the same BaseConformer by creating random or
 * likely torsion sets using the deriveConformer() method.
 * A BaseConformer also keeps track of EliminationRules, which is a dictionary of bond torsion combinations
 * that cause severe atom clashes. The dictionary grows when trying to assemble new Conformers based on
 * torsion sets, which cause clashes.
 */
public class BaseConformer extends Conformer {
	private static final double ANGLE_TOLERANCE = 0.001f;	// limit for considering bonds as parallel

	// For every optimum torsion we check for collisions and in case we try, whether left or right
	// range ends are substantially better.
	private static final double ACCEPTABLE_CENTER_STRAIN = 0.05f;	// if we have less strain than this we don't check whether range edges improve strain
	private static final double NECESSARY_EDGE_STRAIN_IMPROVEMENT_FACTOR = 0.8f; // necessary strain improvement to use edge torsion rather than center torsion
	private static final double MAXIMUM_CENTER_STRAIN = 0.2f;		// if center torsion is above this limit we refuse that torsion

	// If no acceptable torsion remains after checking and attempting to use range edge values,
	// we try to modify the least bad torsion stepwise until it is acceptable.
	private static final int ESCAPE_ANGLE = 8;  // degrees to rotate the rotatable bond to escape collisions
	private static final int ESCAPE_STEPS = 4;	// how often we apply this rotation trying to solve the collision
	private static final double MIN_ESCAPE_GAIN_PER_STEP = 0.05;

	private short[][] mTorsion;
	private short[][] mFrequency;
	private double[][] mLikelyhood; // considering directly connected rigid fragments (frequency and collision strain)
	private int[] mBestTorsionIndex;
	private RotatableBond[] mRotatableBond;
	private ArrayList<TorsionSetEliminationRule> mEliminationRuleList;
	private Random mRandom;
	private int mConformerCount;

	public BaseConformer(StereoMolecule mol, RigidFragment[] rigidFragment, RotatableBond[] rotatableBond, int[] fragmentPermutation, Random random) {
		super(mol);

		mRandom = random;
		mRotatableBond = rotatableBond;

		mTorsion = new short[rotatableBond.length][];
		mFrequency = new short[rotatableBond.length][];
		mLikelyhood = new double[rotatableBond.length][];
		short[][][] defaultTorsionRange = new short[rotatableBond.length][][];
		for (int bondIndex=0; bondIndex<rotatableBond.length; bondIndex++) {
			mTorsion[bondIndex] = rotatableBond[bondIndex].getDefaultTorsions().clone();
			mFrequency[bondIndex] = rotatableBond[bondIndex].getDefaultFrequencies().clone();
			mLikelyhood[bondIndex] = new double[mTorsion[bondIndex].length];
			defaultTorsionRange[bondIndex] = rotatableBond[bondIndex].getDefaultTorsionRanges();
		}

		mBestTorsionIndex = new int[rotatableBond.length];
		boolean[] isAttached = new boolean[rigidFragment.length];
		for (int bondIndex=0; bondIndex<rotatableBond.length; bondIndex++)
			connectFragments(rotatableBond[bondIndex], bondIndex, isAttached, fragmentPermutation, defaultTorsionRange[bondIndex]);

		// for separated fragments without connection points we need to get coordinates
		for (int i=0; i<rigidFragment.length; i++)
			if (!isAttached[i])
				for (int j=0; j<rigidFragment[i].getCoreSize(); j++)
					setCoordinates(rigidFragment[i].coreToOriginalAtom(j),
							rigidFragment[i].getCoreCoordinates(fragmentPermutation[i], j));

		mEliminationRuleList = new ArrayList<>();
	}

	public RotatableBond[] getRotatableBonds() {
		return mRotatableBond;
	}

	public int getTorsion(int bondIndex, int torsionIndex) {
		return mTorsion[bondIndex][torsionIndex];
	}

	public double getTorsionLikelyhood(int bondIndex, int torsionIndex) {
		return mLikelyhood[bondIndex][torsionIndex];
	}

	public int[] getMostLikelyTorsionIndexes() {
		return mBestTorsionIndex;
	}

	public Conformer deriveConformer(int[] torsionIndex, String name) {
		mConformerCount++;

		Conformer conformer = new Conformer(this);
		conformer.setName(name);
		for (int rb=mRotatableBond.length-1; rb>=0; rb--)
			rotateToIndex(conformer, mRotatableBond[rb], rb, torsionIndex[rb]);
		return conformer;
	}

	/**
	 * @return number of derived (not necessarily successful) conformers
	 */
	public int getDerivedConformerCount() {
		return mConformerCount;
	}

	public ArrayList<TorsionSetEliminationRule> getEliminationRules() {
		return mEliminationRuleList;
	}

	/**
	 * Calculates a random torsion index giving torsions with higher likelyhoods
	 * (i.e. frequencies and collision strains) a higher chance to be selected.
	 * With a progress value of 0.0 selection likelyhoods are proportional to
	 * the torsion frequencies. With increasing progress value the higher frequent
	 * torsions are less and less preferred until 1.0 without any preference.
	 * @param random
	 * @param progress 0...1
	 */
	public int getLikelyRandomTorsionIndex(int bondIndex, double random, double progress) {
		double sum = 0;
		for (int t=0; t<mTorsion[bondIndex].length; t++) {
			double contribution = (1f-progress)*mLikelyhood[bondIndex][t] + progress/mTorsion[bondIndex].length;
			sum += contribution;
			if (random <= sum)
				return t;
		}
		return mTorsion[bondIndex].length-1;  // should never reach this
	}

	/**
	 * Checks both rigid fragments that are connected by this bond, whether they have
	 * been attached to the conformer yet, i.e. whether their local coordinates have been
	 * copied to conformer and transformed according to the connecting torsion.
	 * If one was already attached, then the other's coordinates are transformed according
	 * to the torsion and copied to the conformer. A likelihood is calculated for every torsion
	 * value based on its frequency and the collision strain of the two fragments' atoms.
	 * If both fragments were not attached yet, then the larger one's coordinates are
	 * copied and the smaller one's coordinates are translated and then copied.
	 * Unlikely torsions, i.e. where collisions strain outweighs frequency, are removed from torsion list.
	 */
	private void connectFragments(RotatableBond rotatableBond, int bondIndex, boolean[] isAttached, int[] fragmentPermutation, short[][] defaultTorsionRange) {
		int fragmentNo1 = rotatableBond.getFragmentNo(0);
		int fragmentNo2 = rotatableBond.getFragmentNo(1);
		RigidFragment fragment1 = rotatableBond.getFragment(0);
		RigidFragment fragment2 = rotatableBond.getFragment(1);
		if (!isAttached[fragmentNo1] && !isAttached[fragmentNo2]) {
			RigidFragment largerFragment = (fragment1.getCoreSize() > fragment2.getCoreSize()) ? fragment1 : fragment2;
			int largerFragmentNo = (fragment1.getCoreSize() > fragment2.getCoreSize()) ? fragmentNo1 : fragmentNo2;
			isAttached[largerFragmentNo] = true;
			int fragmentConformer = (fragmentPermutation == null) ? 0 : fragmentPermutation[largerFragmentNo];
			for (int i=0; i<largerFragment.getExtendedSize(); i++) {
				int atom = largerFragment.extendedToOriginalAtom(i);
				setCoordinates(atom, largerFragment.getExtendedCoordinates(fragmentConformer, i));
			}
		}

		assert(isAttached[fragmentNo1] ^ isAttached[fragmentNo2]);

		int rootAtom,rearAtom,fragmentNo,bondAtomIndex;
		RigidFragment fragment;
		if (isAttached[fragmentNo1]) {
			fragmentNo = fragmentNo2;
			fragment = fragment2;
			bondAtomIndex = rotatableBond.areBondAtomsInFragmentOrder() ? 1 : 0;
		}
		else {
			fragmentNo = fragmentNo1;
			fragment = fragment1;
			bondAtomIndex = rotatableBond.areBondAtomsInFragmentOrder() ? 0 : 1;
		}

		int bond = rotatableBond.getBond();
		rootAtom = getMolecule().getBondAtom(bondAtomIndex, bond);
		rearAtom = getMolecule().getBondAtom(1-bondAtomIndex, bond);

		int fragmentConformer = (fragmentPermutation == null) ? 0 : fragmentPermutation[fragmentNo];

		int fRootAtom = fragment.originalToExtendedAtom(rootAtom);
		int fRearAtom = fragment.originalToExtendedAtom(rearAtom);

		Coordinates froot = fragment.getExtendedCoordinates(fragmentConformer, fRootAtom);
		Coordinates root = getCoordinates(rootAtom);
		Coordinates fuv = froot.subC(fragment.getExtendedCoordinates(fragmentConformer, fRearAtom)).unit();
		Coordinates uv = root.subC(getCoordinates(rearAtom)).unit();
		double alpha = fuv.getAngle(uv);

		Coordinates[] fcoords = new Coordinates[fragment.getExtendedSize()];

		if (alpha < ANGLE_TOLERANCE) { // special case both axes parallel: we only need to translate
			for (int i=0; i<fragment.getExtendedSize(); i++) {
				int atom = fragment.extendedToOriginalAtom(i);
				if (atom != rootAtom && atom != rearAtom)
					fcoords[i] = (alpha > Math.PI / 2) ?
							froot.subC(fragment.getExtendedCoordinates(fragmentConformer, i))
							: fragment.getExtendedCoordinates(fragmentConformer, i).subC(froot);
			}
		}
		else {
			Coordinates rotationAxis;
			if (alpha < Math.PI - ANGLE_TOLERANCE) {    // normal case, just rotate around the plane orthogonal
				rotationAxis = uv.cross(fuv);
			}
			else {    // special case both axes anti-parallel: for cross-product we need any vector being different from uv
				if (Math.abs(uv.x) >= Math.abs(uv.y) && Math.abs(uv.x) >= Math.abs(uv.z))
					rotationAxis = new Coordinates(-(uv.y + uv.z) / uv.x, 1.0, 1.0);
				else if (Math.abs(uv.y) >= Math.abs(uv.x) && Math.abs(uv.y) >= Math.abs(uv.z))
					rotationAxis = new Coordinates(1.0, -(uv.x + uv.z) / uv.y, 1.0);
				else
					rotationAxis = new Coordinates(1.0, 1.0, -(uv.x + uv.y) / uv.z);
			}

			double[][] m = getRotationMatrix(rotationAxis.unit(), alpha);

			for (int i=0; i<fragment.getExtendedSize(); i++) {
				int atom = fragment.extendedToOriginalAtom(i);
				if (atom != rootAtom && atom != rearAtom)
					fcoords[i] = fragment.getExtendedCoordinates(fragmentConformer, i).subC(froot).rotate(m);
			}
		}
		isAttached[fragmentNo] = true;

		// we need to restore valid coordinates for rootAtom and neighbors
		// to correctly calculate torsion from atom sequence
		for (int i=0; i<getMolecule().getConnAtoms(rootAtom); i++) {
			int connAtom = getMolecule().getConnAtom(rootAtom, i);
			if (connAtom != rearAtom) {
				int fAtom = fragment.originalToExtendedAtom(connAtom);
				setCoordinatesReplace(connAtom, fcoords[fAtom].addC(root));
			}
		}

		double startTorsion = TorsionDB.calculateTorsionExtended(this, rotatableBond.getTorsionAtoms());

		short currentTorsion = -1;
		int bestTorsionEdgeUsed = 0;	// 0:peakCenter, 1:leftEdge, 2:rightEdge
		double bestTorsionStrain = 0;
		for (int t=0; t<mTorsion[bondIndex].length; t++) {
			currentTorsion = mTorsion[bondIndex][t];
			double[][] m = getRotationMatrix(uv, Math.PI * currentTorsion / 180 - startTorsion);
			for (int i=0; i<fragment.getExtendedSize(); i++) {
				int atom = fragment.extendedToOriginalAtom(i);
				if (atom != rootAtom && atom != rearAtom) {
					setCoordinates(atom, fcoords[i]);
					getCoordinates(atom).rotate(m).add(root);
				}
			}

			double centerStrain = calculateCollisionStrain(rotatableBond);
			double usedStrain = centerStrain;	// default
			int torsionEdgeUsed = 0;	// default
			// if the strain is above a certain limit, we investigate whether we should use one
			// of the torsion range limits rather than the central torsion value, which has the
			// highest frequency.
			// If one of the torsion range limits causes a lower strain, then we use that instead
			// of the center value and apply the original center frequency value for it.
			if (centerStrain < ACCEPTABLE_CENTER_STRAIN) {
				double relativeStrain = centerStrain / MAXIMUM_CENTER_STRAIN;
				mLikelyhood[bondIndex][t] = mFrequency[bondIndex][t] * (1.0 - relativeStrain * relativeStrain);
			}
			else {
				int centerTorsion = mTorsion[bondIndex][t];
				double usedRangeStrain = Double.MAX_VALUE;
				for (int r=0; r<2; r++) {
					currentTorsion = defaultTorsionRange[t][r];
					if (currentTorsion != centerTorsion) {
						m = getRotationMatrix(uv, Math.PI * currentTorsion / 180 - startTorsion);
						for (int i=0; i<fragment.getExtendedSize(); i++) {
							int atom = fragment.extendedToOriginalAtom(i);
							if (atom != rootAtom && atom != rearAtom) {
								setCoordinates(atom, fcoords[i]);
								getCoordinates(atom).rotate(m).add(root);
							}
						}

						double rangeStrain = calculateCollisionStrain(rotatableBond);
						if (rangeStrain < centerStrain * NECESSARY_EDGE_STRAIN_IMPROVEMENT_FACTOR
						 && rangeStrain < usedRangeStrain) {
							mTorsion[bondIndex][t] = defaultTorsionRange[t][r];
							double relativeStrain = rangeStrain / (NECESSARY_EDGE_STRAIN_IMPROVEMENT_FACTOR * MAXIMUM_CENTER_STRAIN);
							mLikelyhood[bondIndex][t] = mFrequency[bondIndex][t] * (1f - relativeStrain * relativeStrain);
							usedStrain = rangeStrain;
							torsionEdgeUsed = r+1;
							usedRangeStrain = rangeStrain;
						}
					}
				}
				if (usedRangeStrain == Double.MAX_VALUE
				 && centerStrain < MAXIMUM_CENTER_STRAIN) {
					// Range strains aren't better and center strain is below allowed max.
					// Calculate likelyhood from center strain.
					double relativeStrain = centerStrain / MAXIMUM_CENTER_STRAIN;
					mLikelyhood[bondIndex][t] = mFrequency[bondIndex][t] * (1f - relativeStrain * relativeStrain);
				}
			}
			if (mLikelyhood[bondIndex][mBestTorsionIndex[bondIndex]] < mLikelyhood[bondIndex][t]) {
				mBestTorsionIndex[bondIndex] = t;
				bestTorsionEdgeUsed = torsionEdgeUsed;
				bestTorsionStrain = usedStrain;	// this is the strain with the highest likelyhood (not necessarily the lowest strain)
			}
		}

		double totalLikelyhood = 0f;
		for (int t=0; t<mTorsion[bondIndex].length; t++)
			if (mLikelyhood[bondIndex][t] > 0f)
				totalLikelyhood += mLikelyhood[bondIndex][t];

		// make sure, we have at least one torsion with positive likelyhood, because only those are considered later
		if (mLikelyhood[bondIndex][mBestTorsionIndex[bondIndex]] <= 0f) {
			mLikelyhood[bondIndex][mBestTorsionIndex[bondIndex]] = 1.0f;
			int angle = bestTorsionEdgeUsed == 1 ? -ESCAPE_ANGLE
					: bestTorsionEdgeUsed == 2 ? ESCAPE_ANGLE
					: (mRandom.nextDouble() < 0.5) ? -ESCAPE_ANGLE : ESCAPE_ANGLE;
			for (int step=1; step<=ESCAPE_STEPS; step++) {
				currentTorsion = (short)(mTorsion[bondIndex][mBestTorsionIndex[bondIndex]]+angle*step);
				double[][] m = getRotationMatrix(uv, Math.PI * currentTorsion / 180 - startTorsion);
				for (int i=0; i<fragment.getExtendedSize(); i++) {
					int atom = fragment.extendedToOriginalAtom(i);
					if (atom != rootAtom && atom != rearAtom) {
						setCoordinates(atom, fcoords[i]);
						getCoordinates(atom).rotate(m).add(root);
					}
				}

				double escapeStrain = calculateCollisionStrain(rotatableBond);
				if (bestTorsionStrain - escapeStrain < MIN_ESCAPE_GAIN_PER_STEP)
					break;

				mTorsion[bondIndex][mBestTorsionIndex[bondIndex]] = currentTorsion;
			}
		}
		else {
			for (int t=0; t<mTorsion[bondIndex].length; t++)
				mLikelyhood[bondIndex][t] /= totalLikelyhood;
		}

//		bestTorsionIndex = buildSortedTorsionMap(bondIndex, bestTorsionIndex);

//System.out.print("connectFragments() applied torsions:"); for (int t=0; t<mTorsion.length; t++) System.out.print(mTorsion[t]+" "); System.out.println();
//System.out.print("connectFragments() applied likelyhoods:"); for (int t=0; t<mTorsion.length; t++) System.out.print(mLikelyhood[t]+" "); System.out.println();
//System.out.println();

		if (currentTorsion != mTorsion[bondIndex][mBestTorsionIndex[bondIndex]]) {
			double[][] m = getRotationMatrix(uv, Math.PI * mTorsion[bondIndex][mBestTorsionIndex[bondIndex]] / 180 - startTorsion);
			for (int i=0; i<fragment.getExtendedSize(); i++) {
				int atom = fragment.extendedToOriginalAtom(i);
				if (atom != rootAtom && atom != rearAtom) {
					setCoordinates(atom, fcoords[i]);
					getCoordinates(atom).rotate(m).add(root);
				}
			}
		}

		setBondTorsion(bond, mTorsion[bondIndex][mBestTorsionIndex[bondIndex]]);
	}

	private double calculateCollisionStrain(RotatableBond rotatableBond) {
		RigidFragment fragment1 = rotatableBond.getFragment(0);
		RigidFragment fragment2 = rotatableBond.getFragment(1);
		int bond = rotatableBond.getBond();
		double panalty = 0f;
		int bondAtom1 = getMolecule().getBondAtom(0, bond);
		int bondAtom2 = getMolecule().getBondAtom(1, bond);
		for (int i=0; i<fragment1.getCoreSize(); i++) {
			int atom1 = fragment1.coreToOriginalAtom(i);
			if (atom1 != bondAtom1 && atom1 != bondAtom2) {
				double vdwr1 = ConformerGenerator.getToleratedVDWRadius(getMolecule().getAtomicNo(atom1));
				for (int j=0; j<fragment2.getCoreSize(); j++) {
					int atom2 = fragment2.coreToOriginalAtom(j);
					if (atom2 != bondAtom1 && atom2 != bondAtom2) {
						double minDistance = vdwr1 + ConformerGenerator.getToleratedVDWRadius(
								getMolecule().getAtomicNo(atom2));
						double dx = Math.abs(getX(atom1) - getX(atom2));
						if (dx < minDistance) {
							double dy = Math.abs(getY(atom1) - getY(atom2));
							if (dy < minDistance) {
								double dz = Math.abs(getZ(atom1) - getZ(atom2));
								if (dz < minDistance) {
									double distance = Math.sqrt(dx*dx+dy*dy+dz*dz);
									if (distance < minDistance) {
										double p = (minDistance - distance) / minDistance;
										panalty += p*p;
									}
								}
							}
						}
					}
				}
			}
		}
		return panalty;
	}

	private static double[][] getRotationMatrix(Coordinates n, double alpha) {
		double sinAlpha = Math.sin(alpha);
		double cosAlpha = Math.cos(alpha);
		double invcosAlpha = 1.0-cosAlpha;

		// rotation matrix is:  m11 m12 m13
		//					    m21 m22 m23
		//					    m31 m32 m33
		double[][] m = new double[3][3];
		m[0][0] = n.x*n.x*invcosAlpha+cosAlpha;
		m[1][1] = n.y*n.y*invcosAlpha+cosAlpha;
		m[2][2] = n.z*n.z*invcosAlpha+cosAlpha;
		m[0][1] = n.x*n.y*invcosAlpha-n.z*sinAlpha;
		m[1][2] = n.y*n.z*invcosAlpha-n.x*sinAlpha;
		m[2][0] = n.z*n.x*invcosAlpha-n.y*sinAlpha;
		m[0][2] = n.x*n.z*invcosAlpha+n.y*sinAlpha;
		m[1][0] = n.y*n.x*invcosAlpha+n.z*sinAlpha;
		m[2][1] = n.z*n.y*invcosAlpha+n.x*sinAlpha;
		return m;
	}

	/**
	 * Rotate the smaller side of the molecule around this bond
	 * to reach the torsion angle defined by torsionIndex.
	 * @param torsionIndex
	 */
	public void rotateToIndex(Conformer conformer, RotatableBond rotatableBond, int bondIndex, int torsionIndex) {
		rotateTo(conformer, rotatableBond, mTorsion[bondIndex][torsionIndex]);
	}

	/**
	 * Rotate the smaller side of the molecule around this bond
	 * to reach the defined torsion angle.
	 * @param torsion final torsion in degrees
	 */
	public void rotateTo(Conformer conformer, RotatableBond rotatableBond, short torsion) {
		while (torsion < 0)
			torsion += 360;
		while (torsion >= 360)
			torsion -= 360;

		int bond = rotatableBond.getBond();
		if (torsion != conformer.getBondTorsion(bond)) {
			int deltaTorsion = torsion - conformer.getBondTorsion(bond);
			rotateSmallerSide(conformer, rotatableBond, Math.PI * deltaTorsion / 180.0);
			conformer.setBondTorsion(bond, torsion);
		}
	}

	/**
	 * Rotates all atoms in atomList around the axis leading from atom1 through atom2
	 * by angle. The coordinates of atom1 and atom2 are not touched.
	 * @param alpha
	 */
	private void rotateSmallerSide(Conformer conformer, RotatableBond rotatableBond, double alpha) {
		int[] torsionAtom = rotatableBond.getTorsionAtoms();
		int rotationCenter = rotatableBond.getRotationCenter();
		Coordinates t2 = conformer.getCoordinates(torsionAtom[2]);
		Coordinates unit = t2.subC(conformer.getCoordinates(torsionAtom[1])).unit();

		double[][] m = getRotationMatrix(unit, (rotationCenter == torsionAtom[1]) ? alpha : -alpha);

		for (int atom:rotatableBond.getSmallerSideAtoms()) {
			if (atom != rotationCenter) {
				double x = conformer.getX(atom) - t2.x;
				double y = conformer.getY(atom) - t2.y;
				double z = conformer.getZ(atom) - t2.z;
				conformer.setX(atom, x*m[0][0]+y*m[0][1]+z*m[0][2] + t2.x);
				conformer.setY(atom, x*m[1][0]+y*m[1][1]+z*m[1][2] + t2.y);
				conformer.setZ(atom, x*m[2][0]+y*m[2][1]+z*m[2][2] + t2.z);
			}
		}
	}
}
