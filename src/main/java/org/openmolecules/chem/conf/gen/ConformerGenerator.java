/*
 * @(#)ConformerGenerator.java
 *
 * Copyright 2013 openmolecules.org, Inc. All Rights Reserved.
 * 
 * NOTICE: All information contained herein is, and remains the property
 * of openmolecules.org. The intellectual and technical concepts contained
 * herein are proprietary to openmolecules.org.
 * Actelion Pharmaceuticals Ltd. is granted a non-exclusive, not transferable
 * and timely unlimited usage license.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.gen;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.chem.conf.VDWRadii;
import com.actelion.research.util.DoubleFormat;
import com.actelion.research.util.IntArrayComparator;
import org.openmolecules.chem.conf.so.ConformationSelfOrganizer;
import org.openmolecules.chem.conf.so.SelfOrganizedConformer;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;
import java.util.TreeMap;

/**
 * This class generates 3D-conformers of a given molecule using the following strategy:
 * <li>All rotatable, non-ring bonds are determined.
 * <li>The fragments between rotatable bonds are considered rigid.
 * <li>For every fragment the relative atom coordinates are determined using a self organization based algorithm.
 * <li>For every rotatable bond a list of preferred torsion angles is determined based on from a CSD statistics of similar bond environments.
 * <li>For individual torsion values likelihoods are estimated based on frequency and atom collisions of vicinal fragments.
 * <li>A dedicated (systematic, biased or random) torsion set strategy delivers collision-free torsion sets, i.e. conformers.
 */
public class ConformerGenerator {
	public static final int STRATEGY_LIKELY_SYSTEMATIC = 1;
	public static final int STRATEGY_PURE_RANDOM = 2;
	public static final int STRATEGY_LIKELY_RANDOM = 3;
	public static final int STRATEGY_ADAPTIVE_RANDOM = 4;

	protected static final double VDW_TOLERANCE_HYDROGEN = 0.80;  // factor on VDW radii for minimum tolerated non bound atom distances
	protected static final double VDW_TOLERANCE_OTHER = 0.80;     // factor on VDW radii for minimum tolerated non bound atom distances

	private static final int ESCAPE_ANGLE = 8;  // degrees to rotate two adjacent rotatable bonds to escape collisions
	private static final int ESCAPE_STEPS = 4;	// how often we apply this rotation trying to solve the collision
	private static final double MIN_ESCAPE_GAIN_PER_STEP = 0.05;

	private StereoMolecule		mMolecule;
	private TreeMap<int[],Conformer> mBaseConformerMap;
	private RotatableBond[]		mRotatableBond;
	private Rigid3DFragment[]	mRigidFragment;
	private ConformationSelfOrganizer mSelfOrganizer;
	private TorsionSetStrategy	mTorsionSetStrategy;
	private TorsionSet			mTorsionSet;
	private long				mRandomSeed;
	private int					mDisconnectedFragmentCount,mConformerCount;
	private boolean				mUseSelfOrganizerIfAllFails;
	private double				mContribution;
	private int[]				mFragmentNo,mDisconnectedFragmentNo,mDisconnectedFragmentSize;
	private boolean[][]			mSkipCollisionCheck;
	private Random				mRandom;

	public static final boolean PRINT_TORSION_AND_FRAGMENT_LIKELYHOODS = false;
	public static final boolean PRINT_DEBUG_INDEXES = false;
	public static final boolean PRINT_EXIT_REASON = false;
	public static final boolean PRINT_ELIMINATION_RULES_WITH_STRUCTURES = false;

public static final String DW_FRAGMENTS_FILE = "/home/thomas/data/debug/conformationGeneratorFragments.dwar";
public String mDiagnosticCollisionString,mDiagnosticTorsionString;	// TODO get rid of this
public int[] mDiagnosticCollisionAtoms;	// TODO get rid of this
public static boolean WRITE_DW_FRAGMENT_FILE = false;

	/**
	 * Adds explicit hydrogen atoms where they are implicit by filling valences
	 * and adapting for atom charges. New hydrogen atoms receive new 2D-coordinates
	 * by equally locating them between those two neighbors with the widest angle between
	 * their bonds. Any stereo configurations deducible from 2D-coordinates are retained.
	 * @param mol
	 */
	public static void addHydrogenAtoms(StereoMolecule mol) {
		// We may have parities but empty coordinates. In this case we need to protect parities.
		boolean paritiesValid = (mol.getHelperArrayStatus() & Molecule.cHelperBitParities) != 0;

		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		int[] implicitHydrogen = new int[mol.getAtoms()];
		for (int atom=0; atom<mol.getAtoms(); atom++)
			implicitHydrogen[atom] = mol.getImplicitHydrogens(atom);

		double avbl = mol.getAverageBondLength();
		for (int atom=0; atom<implicitHydrogen.length; atom++)
			if (implicitHydrogen[atom] != 0)
				for (int i=0; i<implicitHydrogen[atom]; i++)
					mol.addBond(atom, mol.addAtom(1), Molecule.cBondTypeSingle);

		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		for (int atom=0; atom<implicitHydrogen.length; atom++)
			if (implicitHydrogen[atom] != 0)
				setHydrogenLocations(mol, atom, implicitHydrogen[atom], avbl);

		// addAtom() and addBond() clear the helper status, i.e. flag all helpers as invalid.
		// Adding hydrogens does not destroy parities. Though, we may flag them to be valid again.
		if (paritiesValid)
			mol.setParitiesValid(0);
		}

	/**
	 * Finds the widest open angle between all connected non-stereo bonds of atom, divides this angle
	 * into hydrogenCount+1 equals parts and sets atom coordinates of hydrogenCount new hydrogen atoms
	 * such, that they equally occupy the space. Helper arrays are assumed to have the state before
	 * adding any new hydrogen atoms.
	 * @param mol
	 * @param atom
	 * @param newHydrogenCount new hydrogen atoms added to atom
	 * @param avbl
	 */
	private static void setHydrogenLocations(StereoMolecule mol, int atom, int newHydrogenCount, double avbl) {

		int firstNewHydrogenNeighbour = mol.getAllConnAtoms(atom) - newHydrogenCount;

		int stereoBondCount = 0;
		for (int i=0; i<firstNewHydrogenNeighbour; i++)
			if (mol.isStereoBond(mol.getConnBond(atom,i)))
				stereoBondCount++;

		double[] angle = null;
		if (stereoBondCount < firstNewHydrogenNeighbour) {
			angle = new double[firstNewHydrogenNeighbour-stereoBondCount];
			int bond = 0;
			for (int i=0; i<firstNewHydrogenNeighbour; i++)
				if (!mol.isStereoBond(mol.getConnBond(atom, i)))
					angle[bond++] = mol.getBondAngle(atom, mol.getConnAtom(atom, i));
	
			Arrays.sort(angle);
			}

		double angleIncrement = 2.0*Math.PI/newHydrogenCount;
		double startAngle = 0.0;

		if (angle != null) {
			double biggestAngleDif = 0.0;
			for (int i=0; i<angle.length; i++) {
				double a1 = (i == 0) ? angle[angle.length-1] - Math.PI*2.0 : angle[i-1];
				double a2 = angle[i];
				if (biggestAngleDif < a2 - a1) {
					biggestAngleDif = a2 - a1;
					startAngle = a1;
					}
				}
			angleIncrement = biggestAngleDif / (newHydrogenCount + 1);
			}

		for (int i=0; i<newHydrogenCount; i++) {
			startAngle += angleIncrement;
			int newHydrogen = mol.getConnAtom(atom, firstNewHydrogenNeighbour+i);
			mol.setAtomX(newHydrogen, mol.getAtomX(atom) + avbl * Math.sin(startAngle));
			mol.setAtomY(newHydrogen, mol.getAtomY(atom) + avbl * Math.cos(startAngle));
			}
		}

	public static double getToleratedVDWRadius(int atomicNo) {
		float vdwr = VDWRadii.VDW_RADIUS[atomicNo < VDWRadii.VDW_RADIUS.length ? atomicNo : 6];
		return vdwr * (atomicNo == 1 ? VDW_TOLERANCE_HYDROGEN : VDW_TOLERANCE_OTHER);
		}

	/**
	 * Instantiates a ConformerGenerator for creating not reproducible conformers.
	 */
	public ConformerGenerator() {
		this(0L);
		}

	/**
	 * @param seed != 0L if conformers shall be created in a reproducible way
	 */
	public ConformerGenerator(long seed) {
		TorsionDB.initialize(TorsionDB.MODE_ANGLES);
		mRandomSeed = seed;
		mRandom = (seed == 0) ? new Random() : new Random(seed);
		Rigid3DFragment.setRandomSeed(seed);
		}

	/**
	 * Fills all free valences of mol with explicit hydrogens and tries to
	 * create a reasonable conformer by starting with the most likely torsion set.
	 * If there are collisions, then less likely torsions are tried to find
	 * a collision free conformer. If it succeeds, mol receives the modified
	 * atom coordinates and mol is returned. If the conformer generation fails,
	 * then null is returned. The torsion strategy used is STRATEGY_ADAPTIVE_RANDOM.
	 * New 3D-coordinates correctly reflect E/Z and R/S bond/atom parities.
	 * This is a convenience method that does not require any initialization.
	 * @param mol the molecule that receive new 3D coordinates in place
	 * @return original molecule with new 3D-coordinates or null
	 */
	public StereoMolecule getOneConformerAsMolecule(StereoMolecule mol) {
		Conformer conformer = getOneConformer(mol);
		return (conformer == null) ? null : conformer.toMolecule(mol);
		}

	/**
	 * Fills all free valences of mol with explicit hydrogens and tries to
	 * create a reasonable conformer by starting with the most likely torsion set.
	 * If there are collisions, then less likely torsions are tried to find
	 * a collision free conformer. If this fails, then null is returned.
	 * The torsion strategy used is STRATEGY_ADAPTIVE_RANDOM.
	 * New 3D-coordinates correctly reflect E/Z and R/S bond/atom parities.
	 * This is a convenience method that does not require any initialization.
	 * @param mol the molecule from which to create the conformer
	 */
	public Conformer getOneConformer(StereoMolecule mol) {
try {	// TODO remove try catch
		if (!initialize(mol, false))
			return null;

		if (mRotatableBond != null) {
			mTorsionSetStrategy = new TorsionSetStrategyAdaptiveRandom(mRotatableBond, mRigidFragment, true, true, mRandomSeed);
			mTorsionSetStrategy.setMaxTotalCount(400);
			mBaseConformerMap = new TreeMap<int[],Conformer>(new IntArrayComparator());
			return getNextConformer();
			}
		else {
			ConformationSelfOrganizer sampler = new ConformationSelfOrganizer(mol, true);
			Conformer conformer = sampler.generateOneConformer(mRandomSeed);
			separateDisconnectedFragments(conformer);
			return conformer;
			}
} catch (Exception e) { e.printStackTrace(); return null; }
		}

		/**
		 * Adds implicit hydrogens to the molecule and determines all rotatable bonds,
		 * which are not part of a ring. Generates rigid fragments between rotatable bonds.
		 * The base conformer is constructed by connecting all rigid fragments using the
		 * most frequent torsions. Atoms of different fragments in the base conformer may
		 * collide. In order to obtain collision free conformers choose a TorsionStrategy
		 * and call getNextConformer() at least once.
		 * @param mol
		 * @param use60degreeSteps use 60 degree steps for every rotatable bond instead of torsion DB
		 */
	private boolean initialize(StereoMolecule mol, boolean use60degreeSteps) {
		mSelfOrganizer = null;

		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		for (int atom=0; atom<mol.getAtoms(); atom++)
			if (mol.getOccupiedValence(atom) > mol.getMaxValence(atom))
				return false;

		addHydrogenAtoms(mol);

		// we need to protect explicit hydrogens in the fragments that we create from mol
		mol.setHydrogenProtection(true);

		int[] atomParity = null;
		int[] bondParity = null;
		boolean[] atomParityIsPseudo = null;
		boolean[] bondParityIsPseudo = null;
		// if we have parities and no atom coordinates, parities cannot be correctly recreated
		// by the Canonizer and, thus, need to be cached and copied back after the symmetry detection.
		if ((mol.getHelperArrayStatus() & Molecule.cHelperBitParities) != 0) {
			atomParity = new int[mol.getAtoms()];
			atomParityIsPseudo = new boolean[mol.getAtoms()];
			for (int atom=0; atom<mol.getAtoms(); atom++) {
				atomParity[atom] = mol.getAtomParity(atom);
				atomParityIsPseudo[atom] = mol.isAtomParityPseudo(atom);
				}
			bondParity = new int[mol.getBonds()];
			bondParityIsPseudo = new boolean[mol.getBonds()];
			for (int bond=0; bond<mol.getBonds(); bond++) {
				bondParity[bond] = mol.getBondParity(bond);
				bondParityIsPseudo[bond] = mol.isBondParityPseudo(bond);
				}
			}

		// we need symmetry ranks for detecting equivalent torsions
		mol.ensureHelperArrays(Molecule.cHelperSymmetrySimple);

		if (atomParity != null) {
			for (int atom=0; atom<mol.getAtoms(); atom++)
				mol.setAtomParity(atom, atomParity[atom], atomParityIsPseudo[atom]);
			for (int bond=0; bond<mol.getBonds(); bond++)
				mol.setBondParity(bond, bondParity[bond], bondParityIsPseudo[bond]);
		}

		mMolecule = mol;
		mConformerCount = 0;
		mTorsionSet = null;
		mRotatableBond = null;

		// check, whether we have disconnected fragments
		mDisconnectedFragmentNo = new int[mol.getAllAtoms()];
		mDisconnectedFragmentCount = mol.getFragmentNumbers(mDisconnectedFragmentNo, false, true);
		mDisconnectedFragmentSize = new int[mDisconnectedFragmentCount];
		for (int atom=0; atom<mol.getAllAtoms(); atom++)
			mDisconnectedFragmentSize[mDisconnectedFragmentNo[atom]]++;

		boolean[] isRotatableBond = new boolean[mol.getAllBonds()];
		int count = TorsionDB.findRotatableBonds(mol, true, isRotatableBond);
		if (count == 0)
			return true;

		locateInitialFragments(isRotatableBond);

if (WRITE_DW_FRAGMENT_FILE) {
 try {
  BufferedWriter writer = new BufferedWriter(new FileWriter(DW_FRAGMENTS_FILE));
  writer.write("<column properties>");
  writer.newLine();
  writer.write("<columnName=\"Structure\">");
  writer.newLine();
  writer.write("<columnProperty=\"specialType\tidcode\">");
  writer.newLine();
  writer.write("<columnName=\"coords\">");
  writer.newLine();
  writer.write("<columnProperty=\"specialType\tidcoordinates3D\">");
  writer.newLine();
  writer.write("<columnProperty=\"parent\tStructure\">");
  writer.newLine();
  writer.write("</column properties>");
  writer.newLine();
  writer.write("Structure\tcoords");
  writer.newLine();
  for (Rigid3DFragment f:mRigidFragment) {
   StereoMolecule[] fragment = f.getFragment();
   for (int i=0; i<fragment.length; i++) {
	Canonizer canonizer = new Canonizer(fragment[i]);
	String idcode = canonizer.getIDCode();
	String coords = canonizer.getEncodedCoordinates();
	writer.write(idcode + "\t" + coords + "\t");
	writer.newLine();
    }
   }
  writer.close();
  }
 catch (IOException ioe) {}
 }

		mRotatableBond = new RotatableBond[count];
		int rotatableBond = 0;
		for (int bond=0; bond<mol.getBonds(); bond++)
			if (isRotatableBond[bond])
				mRotatableBond[rotatableBond++] = new RotatableBond(mol, bond, mFragmentNo,
						mDisconnectedFragmentNo, mDisconnectedFragmentSize[mDisconnectedFragmentNo[mol.getBondAtom(0, bond)]],
						mRigidFragment, mRandom, use60degreeSteps);

		// sort by descending atom count of smaller side, i.e. we want those bond dividing into equal parts first!
		Arrays.sort(mRotatableBond, new Comparator<RotatableBond>() {
			@Override
			public int compare(RotatableBond b1, RotatableBond b2) {
				int c1 = b1.getSmallerSideAtomCount();
				int c2 = b2.getSmallerSideAtomCount();
				return (c1 == c2) ? 0 : (c1 < c2) ? 1 : -1;
				}
			});

		// TODO this is actually only used with random conformers. Using a conformer mode
		// may save some time, if systematic conformers are created.
		initializeCollisionCheck();

		return true;
		}

	private Conformer getBaseConformer(int[] fragmentPermutation) {
		Conformer baseConformer = mBaseConformerMap.get(fragmentPermutation);
		if (baseConformer != null)
			return baseConformer;

		baseConformer = new Conformer(mMolecule);
		boolean[] isAttached = new boolean[mRigidFragment.length];
		for (RotatableBond rb:mRotatableBond)
			rb.connectFragments(baseConformer, isAttached, fragmentPermutation);

		// for separated fragments without connection points we need to get coordinates
		for (int i=0; i<mRigidFragment.length; i++)
			if (!isAttached[i])
				for (int j=0; j<mRigidFragment[i].getCoreSize(); j++)
					baseConformer.setCoordinates(mRigidFragment[i].coreToOriginalAtom(j),
							mRigidFragment[i].getCoreCoordinates(fragmentPermutation[i], j));

		mBaseConformerMap.put(fragmentPermutation, baseConformer);

		if (PRINT_TORSION_AND_FRAGMENT_LIKELYHOODS) {
			for (int rb=0; rb<mRotatableBond.length; rb++) {
				System.out.print("RotBond["+rb+"]("+mRotatableBond[rb].getBond()+"):");
				for (int i=0; i<mRotatableBond[rb].getTorsionCount(); i++) {
					System.out.print(" "+mRotatableBond[rb].getTorsion(i)+"("+(int)(100*mRotatableBond[rb].getTorsionLikelyhood(i))+"%)");
					}
				System.out.println();
				}
			for (int rb=0; rb<mRigidFragment.length; rb++)
				System.out.println("Fragment["+rb+"] conformers:"+mRigidFragment[rb].getConformerCount());
			}

		return baseConformer;
		}

	/**
	 * Creates the next random, likely or systematic new(!) conformer of the molecule
	 * that was passed when calling initializeConformers(). A new conformer is one,
	 * whose combination of torsion angles was not used in a previous conformer
	 * created by this function since the last call of initializeConformers().
	 * Parameter mol may be null or recycle the original molecule to receive new 3D coordinates.
	 * If it is null, then a fresh copy of the original molecule with new atom coordinates is returned.
	 * Every call of this method creates a new collision-free conformer until the employed torsion set
	 * strategy decides that it cannot generate any more suitable torsion sets.
	 * @param mol null or molecule used during initialization or a copy of it
	 * @return conformer or null, if all/maximum torsion permutations have been tried
	 */
	public StereoMolecule getNextConformerAsMolecule(StereoMolecule mol) {
		Conformer conformer = getNextConformer();
		return (conformer == null) ? null : conformer.toMolecule(mol);
		}

	/**
	 * Creates the next random, likely or systematic new(!) conformer of the molecule
	 * that was passed when calling initializeConformers(). A new conformer is one,
	 * whose combination of torsion angles was not used in a previous conformer
	 * created by this function since the last call of initializeConformers().
	 * Every call of this method creates a new collision-free conformer until the employed torsion set
	 * strategy decides that it cannot generate any more suitable torsion sets.
	 * @return conformer or null, if all/maximum torsion permutations have been tried
	 */
	public Conformer getNextConformer() {
		if (mRotatableBond == null && mSelfOrganizer == null)
			return null;

		if (mSelfOrganizer != null) {
			SelfOrganizedConformer conformer = mSelfOrganizer.getNextConformer();
			if (conformer != null) {
				separateDisconnectedFragments(conformer);
				mConformerCount++;
				return conformer;
				}

			return null;
			}

		if (mBaseConformerMap == null)
			return null;

		// create a base conformer from first set of fragments and calculate torsion likelyhoods
		if (mBaseConformerMap.size() == 0)
			getBaseConformer(new int[mRigidFragment.length]);

		mTorsionSet = mTorsionSetStrategy.getNextTorsionSet(mTorsionSet);
		while (mTorsionSet != null) {
/*
System.out.println("---- new torsion and conformer index set: -----");
for (int i=0; i<mRotatableBond.length; i++) System.out.println("rb:"+i+" index:"+torsionSet.getTorsionIndexes()[i]
+" torsion:"+mRotatableBond[i].getTorsion(torsionSet.getTorsionIndexes()[i])
+" likelyhood:"+mRotatableBond[i].getTorsionLikelyhood(torsionSet.getTorsionIndexes()[i]));
for (int i=0; i<mRigidFragment.length; i++) System.out.println("rf:"+i+" index:"+torsionSet.getConformerIndexes()[i]
+" likelyhood:"+mRigidFragment[i].getConformerLikelyhood(torsionSet.getConformerIndexes()[i]));
*/

			Conformer conformer = new Conformer(getBaseConformer(mTorsionSet.getConformerIndexes()));

			for (int j=mRotatableBond.length-1; j>=0; j--)
				mRotatableBond[j].rotateToIndex(conformer, mTorsionSet.getTorsionIndexes()[j]);

if (WRITE_DW_FRAGMENT_FILE) {
 mDiagnosticTorsionString = ""+mTorsionSet.getTorsionIndexes()[0];
 for (int i=1; i<mTorsionSet.getTorsionIndexes().length; i++)
  mDiagnosticTorsionString = mDiagnosticTorsionString + ":" + mTorsionSet.getTorsionIndexes()[i];
 mDiagnosticTorsionString = mDiagnosticTorsionString + "<->" + mTorsionSet.getConformerIndexes()[0];
 for (int i=1; i<mTorsionSet.getConformerIndexes().length; i++)
  mDiagnosticTorsionString = mDiagnosticTorsionString + ":" + mTorsionSet.getConformerIndexes()[i];
 }

			// If the torsionSet has already a collision value, then it is a second choice torsion set
			// with a collision value below the tolerance. No need to check again.
			if (mTorsionSet.getCollisionIntensitySum() == 0.0)
				checkCollision(conformer, mTorsionSet);

			// Even if this is a second choice torsion set, we need to fix collisions again,
			// because the torsion set stores indexes, not real torsions.
			if (mTorsionSet.getCollisionIntensitySum() != 0.0)
				tryFixCollisions(conformer, mTorsionSet);

			if (PRINT_DEBUG_INDEXES) {
				String collisionString = (mTorsionSet.getCollisionIntensitySum() != 0.0) ?
						" collides:"+ DoubleFormat.toString(mTorsionSet.getCollisionIntensitySum()) : "";
				System.out.println(mTorsionSet.toString() + collisionString);
				}

			if (mTorsionSet.getCollisionIntensitySum() > mTorsionSetStrategy.calculateCollisionTolerance()) {
//System.out.println("COLLIDES!");

// TODO remove
String idcode,coords;
int elimRules;
if (PRINT_ELIMINATION_RULES_WITH_STRUCTURES) {
Canonizer can = new Canonizer(conformer.toMolecule(null));
idcode = can.getIDCode();
coords = can.getEncodedCoordinates();
elimRules = mTorsionSetStrategy.getEliminationRuleList().size();
}

				mTorsionSet = mTorsionSetStrategy.getNextTorsionSet(mTorsionSet);

// TODO remove
if (PRINT_ELIMINATION_RULES_WITH_STRUCTURES) {
if (elimRules != mTorsionSetStrategy.getEliminationRuleList().size()) {
	if (elimRules == 0)
		System.out.println("idcode\tidcoords\telimRules\tstrategy");
	StringBuilder sb = new StringBuilder();
	for (int i = elimRules; i < mTorsionSetStrategy.getEliminationRuleList().size(); i++) {
		TorsionSetEliminationRule rule = mTorsionSetStrategy.getEliminationRuleList().get(i);
		sb.append(mTorsionSetStrategy.eliminationRuleString(rule));
		sb.append("  m:" + Long.toHexString(rule.getMask()[0]));
		sb.append(" d:" + Long.toHexString(rule.getData()[0]));
		sb.append("<NL>");
	}
	System.out.println(idcode + "\t" + coords + "\t" + sb.toString() + "\tlikelyRandom");
}}


				if (mTorsionSet != null || mConformerCount != 0)
					continue;

				if (mUseSelfOrganizerIfAllFails) {
					mSelfOrganizer = new ConformationSelfOrganizer(mMolecule, true);
					mSelfOrganizer.initializeConformers(mRandomSeed, -1);
					conformer = mSelfOrganizer.getNextConformer();
					if (conformer != null) {
						separateDisconnectedFragments(conformer);
						mConformerCount++;
						return conformer;
						}
					}

				// we didn't get any torsion set that didn't collide; take the best we had
				mTorsionSet = mTorsionSetStrategy.getBestCollidingTorsionIndexes();
				conformer = new Conformer(getBaseConformer(mTorsionSet.getConformerIndexes()));

				for (int j=mRotatableBond.length-1; j>=0; j--)
					mRotatableBond[j].rotateToIndex(conformer, mTorsionSet.getTorsionIndexes()[j]);

				mBaseConformerMap = null;	// we are finished with conformers
				}

//System.out.println("passed collision check! "+mTorsionSet.toString());
			separateDisconnectedFragments(conformer);
			mContribution = mTorsionSetStrategy.getContribution(mTorsionSet);
			mTorsionSet.setUsed();
			mConformerCount++;
			return conformer;
			}

		return null;
		}

	/**
	 * @return count of valid delivered conformers
	 */
	public int getConformerCount() {
		return mConformerCount;
		}

	/**
	 * Calculates the potential count of conformers by multiplying degrees of freedom
	 * (torsions per rotatable bond & rigid fragment multiplicities).
	 * Cannot be called before calling initializeConformers().
	 * @return
	 */
	public int getPotentialConformerCount() {
		return (mTorsionSetStrategy == null) ? 1 : mTorsionSetStrategy.getPermutationCount();
		}

	/**
	 * With best current knowledge about colliding torsion combinations
	 * and based on the individual frequencies of currently active torsions
	 * this method returns the conformers's overall contribution to the
	 * total set of non colliding conformers.
	 * @return this conformer's contribution to all conformers
	 */
	public double getPreviousConformerContribution() {
		return mRotatableBond == null ? 1f : mContribution;
		}

	/**
	 * If a molecule has at least one rotatable bond if all permutations
	 * of torsions collide beyond a tolerated strain, then the standard
	 * behaviour of this class is to return that clashing conformer with
	 * the lowest strain.<br>
	 * If passing true to this method, the ConformerGenerator will use
	 * the ConformerSelfOrganizer in these cases to generate conformers.
	 * getNextConformer will then deliver conformers until the self organization
	 * fails to create new conformers.
	 * @param b
	 */
	public void setUseSelfOrganizerIfAllFails(boolean b) {
		mUseSelfOrganizerIfAllFails = b;
		}

	/**
	 * One of the initializeConformers() methods needs to be called, before getting individual
	 * conformers of the same molecule by getNextConformer().
	 * Open valences of the passed molecule are filled with hydrogen atoms.
	 * The passed molecule may repeatedly be used as container for a new conformer's atom
	 * coordinates, if it is passed as parameter to getNextConformer().
	 * This method uses the STRATEGY_LIKELY_RANDOM strategy with a maximum of 100.000 distinct
	 * torsion sets and uses torsions from crystallographic data.
	 * @param mol will be saturated with hydrogen atoms
	 * @return false if there is a structure problem
	 */
	public boolean initializeConformers(StereoMolecule mol) {
		return initializeConformers(mol, STRATEGY_LIKELY_RANDOM, 100000, false);
		}

	/**
	 * One of the initializeConformers() methods needs to be called, before getting individual
	 * conformers of the same molecule by getNextConformer().
	 * Open valences of the passed molecule are filled with hydrogen atoms.
	 * The passed molecule may repeatedly be used as container for a new conformer's atom
	 * coordinates, if it is passed as parameter to getNextConformer().
	 * @param mol will be saturated with hydrogen atoms
	 * @param strategy one of the STRATEGY_ constants
	 * @param maxTorsionSets maximum number of distinct torsion sets the strategy will try (default 100000)
	 * @param use60degreeSteps use 60 degree steps for every rotatable bond instead of torsion DB
	 * @return false if there is a structure problem
	 */
	public boolean initializeConformers(StereoMolecule mol, int strategy, int maxTorsionSets, boolean use60degreeSteps) {
		if (!initialize(mol, use60degreeSteps))
			return false;

		if (mRotatableBond == null) {
			mSelfOrganizer = new ConformationSelfOrganizer(mol, true);
			mSelfOrganizer.initializeConformers(mRandomSeed, -1);
			}
		else {
			switch(strategy) {
			case STRATEGY_PURE_RANDOM:
				mTorsionSetStrategy = new TorsionSetStrategyRandom(mRotatableBond, mRigidFragment, false, mRandomSeed);
				break;
			case STRATEGY_LIKELY_RANDOM:
				mTorsionSetStrategy = new TorsionSetStrategyRandom(mRotatableBond, mRigidFragment, true, mRandomSeed);
				break;
			case STRATEGY_ADAPTIVE_RANDOM:
				mTorsionSetStrategy = new TorsionSetStrategyAdaptiveRandom(mRotatableBond, mRigidFragment, true, true, mRandomSeed);
				break;
			case STRATEGY_LIKELY_SYSTEMATIC:
				mTorsionSetStrategy = new TorsionSetStrategyLikelySystematic(mRotatableBond, mRigidFragment);
				break;
				}
			mTorsionSetStrategy.setMaxTotalCount(maxTorsionSets);
			mBaseConformerMap = new TreeMap<int[],Conformer>(new IntArrayComparator());
			}

		return true;
		}

	/**
	 * Moves disconnected fragments along Z-axis such that there is an
	 * empty z-space SEPARATION_DISTANCE thick between the fragments.
	 * @param conformer
	 */
	private void separateDisconnectedFragments(Conformer conformer) {
		final double SEPARATION_DISTANCE = 3.0;
		if (mDisconnectedFragmentCount > 1) {
			double[] meanX = new double[mDisconnectedFragmentCount];
			double[] meanY = new double[mDisconnectedFragmentCount];
			double[] minZ = new double[mDisconnectedFragmentCount];
			double[] maxZ = new double[mDisconnectedFragmentCount];
			for (int i=0; i<mDisconnectedFragmentCount; i++) {
				minZ[i] =  1000000000f;
				maxZ[i] = -1000000000f;
				}
			for (int atom=0; atom<conformer.getSize(); atom++) {
				meanX[mDisconnectedFragmentNo[atom]] += conformer.getX(atom);
				meanY[mDisconnectedFragmentNo[atom]] += conformer.getY(atom);
				if (minZ[mDisconnectedFragmentNo[atom]] > conformer.getZ(atom))
					minZ[mDisconnectedFragmentNo[atom]] = conformer.getZ(atom);
				if (maxZ[mDisconnectedFragmentNo[atom]] < conformer.getZ(atom))
					maxZ[mDisconnectedFragmentNo[atom]] = conformer.getZ(atom);
				}
			for (int i=0; i<mDisconnectedFragmentCount; i++) {
				meanX[i] /= mDisconnectedFragmentSize[i];
				meanY[i] /= mDisconnectedFragmentSize[i];
				}
			double[] shiftX = new double[mDisconnectedFragmentCount];
			double[] shiftY = new double[mDisconnectedFragmentCount];
			double[] shiftZ = new double[mDisconnectedFragmentCount];
			for (int i=1; i<mDisconnectedFragmentCount; i++) {
				shiftX[i] = meanX[0] - meanX[i];
				shiftY[i] = meanY[0] - meanY[i];
				shiftZ[i] = shiftZ[i-1] + maxZ[i-1] - minZ[i] + SEPARATION_DISTANCE;
				}
			for (int atom=0; atom<conformer.getSize(); atom++) {
				if (mDisconnectedFragmentNo[atom] != 0) {
					Coordinates c = conformer.getCoordinates(atom);
					c.x += shiftX[mDisconnectedFragmentNo[atom]];
					c.y += shiftY[mDisconnectedFragmentNo[atom]];
					c.z += shiftZ[mDisconnectedFragmentNo[atom]];
					}
				}
			}
		}

	private void locateInitialFragments(boolean[] isRotatableBond) {
		mFragmentNo = new int[mMolecule.getAllAtoms()];
		int fragmentCount = mMolecule.getFragmentNumbers(mFragmentNo, isRotatableBond, true);
		mRigidFragment = new Rigid3DFragment[fragmentCount];
		for (int i=0; i<fragmentCount; i++)
			mRigidFragment[i] = new Rigid3DFragment(mMolecule, mFragmentNo, i);
		}

	/**
	 * Creates a atom matrix flagging all atom pairs, which should be skipped
	 * during collision check, because they are members of the same fragment,
	 * they are members of two adjacent fragments or because the number of
	 * bonds between them is smaller than 3.
	 */
	private void initializeCollisionCheck() {
		mSkipCollisionCheck = new boolean[mMolecule.getAllAtoms()][];
		for (int atom=1; atom<mMolecule.getAllAtoms(); atom++)
			mSkipCollisionCheck[atom] = new boolean[atom];

		// skip collision check for two atoms in adjacent fragments
/* This we cannot do, if we don't sort out colliding torsions in the RotatableBond class,
   which is not straight forward in case we have multiple fragment conformers.
		for (RotatableBond rb:mRotatableBond) {
			Rigid3DFragment f1 = rb.getFragment(0);
			Rigid3DFragment f2 = rb.getFragment(1);
			for (int i=0; i<f1.getCoreSize(); i++) {
				int atom1 = f1.coreToOriginalAtom(i);
				for (int j=0; j<f2.getCoreSize(); j++)
					skipCollisionCheck(atom1, f2.coreToOriginalAtom(j));
				}
			}*/

		// skip collision check for two atoms of the same fragment
		for (Rigid3DFragment rf:mRigidFragment)
			for (int i=1; i<rf.getExtendedSize(); i++)
				for (int j=0; j<i; j++)
					skipCollisionCheck(rf.extendedToOriginalAtom(i), rf.extendedToOriginalAtom(j));

		// skip collision check for atom pairs with 2 bonds in between
		for (int atom=0; atom<mMolecule.getAtoms(); atom++)
			for (int i=1; i<mMolecule.getAllConnAtoms(atom); i++)
				for (int j=0; j<i; j++)
					skipCollisionCheck(mMolecule.getConnAtom(atom, i), mMolecule.getConnAtom(atom, j));

		// skip collision check for any two atoms that belong to different disconnected fragments
		if (mDisconnectedFragmentNo != null)
			for (int atom1=1; atom1<mMolecule.getAllAtoms(); atom1++)
				for (int atom2=0; atom2<atom1; atom2++)
					if (mDisconnectedFragmentNo[atom1] != mDisconnectedFragmentNo[atom2])
						mSkipCollisionCheck[atom1][atom2] = true;
		}

	private void skipCollisionCheck(int atom1, int atom2) {
		if (atom1 < atom2)
			mSkipCollisionCheck[atom2][atom1] = true;
		else
			mSkipCollisionCheck[atom1][atom2] = true;
		}

	private boolean checkCollision(Conformer conformer, TorsionSet torsionSet) {
mDiagnosticCollisionString = "";
mDiagnosticCollisionAtoms = null;
		double collisionIntensitySum = 0;
		double[][] collisionIntensityMatrix = null;
		StereoMolecule mol = conformer.getMolecule();
		for (int atom1=1; atom1<mol.getAllAtoms(); atom1++) {
			double vdwr1 = getToleratedVDWRadius(mol.getAtomicNo(atom1));
			for (int atom2=0; atom2<atom1; atom2++) {
				if (!mSkipCollisionCheck[atom1][atom2]) {
					double minDistance = vdwr1+getToleratedVDWRadius(mol.getAtomicNo(atom2));
					double dx = Math.abs(conformer.getX(atom1) - conformer.getX(atom2));
					if (dx < minDistance) {
						double dy = Math.abs(conformer.getY(atom1) - conformer.getY(atom2));
						if (dy < minDistance) {
							double dz = Math.abs(conformer.getZ(atom1) - conformer.getZ(atom2));
							if (dz < minDistance) {
								double distance = Math.sqrt(dx*dx+dy*dy+dz*dz);
								if (distance < minDistance) {
									double relativeCollision = (minDistance - distance) / minDistance;
									double collisionIntensity = relativeCollision * relativeCollision;
									collisionIntensitySum += collisionIntensity;
//System.out.println("a1:"+atom1+" f1:"+mFragmentNo[atom1]+" a2:"+atom2+" f2:"+mFragmentNo[atom2]+" distance:"+distance+" min:"+minDistance);
if (WRITE_DW_FRAGMENT_FILE) {
 if (mDiagnosticCollisionString.length() != 0) mDiagnosticCollisionString = mDiagnosticCollisionString + "<NL>";
  mDiagnosticCollisionString = mDiagnosticCollisionString+"a1:"+atom1+" f1:"+mFragmentNo[atom1]+" a2:"+atom2+" f2:"+mFragmentNo[atom2]+" distance:"+distance+" min:"+minDistance;
 if (mDiagnosticCollisionAtoms == null) {
  mDiagnosticCollisionAtoms = new int[2];
  mDiagnosticCollisionAtoms[0] = atom1;
  mDiagnosticCollisionAtoms[1] = atom2;
 }
}
									if (collisionIntensityMatrix == null)
										collisionIntensityMatrix = new double[mRigidFragment.length][];
									int f1 = mFragmentNo[atom1];
									int f2 = mFragmentNo[atom2];
									if (f1 < f2) {
										if (collisionIntensityMatrix[f2] == null)
											collisionIntensityMatrix[f2] = new double[f2];
										collisionIntensityMatrix[f2][f1] += collisionIntensity;
										}
									else {
										if (collisionIntensityMatrix[f1] == null)
											collisionIntensityMatrix[f1] = new double[f1];
										collisionIntensityMatrix[f1][f2] += collisionIntensity;
										}
									continue;
									}
								}
							}
						}
					}
				}
			}
		torsionSet.setCollisionIntensity(collisionIntensitySum, collisionIntensityMatrix);
		return (collisionIntensitySum != 0);
		}

	/**
	 * If we have collisions between two fragments that have two rotatable bonds in between,
	 * then we try to rotate both rotatable bonds in any of these situations stepwise
	 * by a small angle to reduce the collision intensity sum. If we achieve a sufficient
	 * improvement without hampering the rest of the molecule, the passed conformer is modified.
	 * @param origConformer
	 * @param torsionSet
	 * @return true if successful; false if no fix possible that reduces total collision intensity below tolerance
	 */
	private boolean tryFixCollisions(Conformer origConformer, TorsionSet torsionSet) {
		double[][] origCollisionIntensityMatrix = torsionSet.getCollisionIntensityMatrix();
		double remainingCollisionIntensity = 0;
		for (int f1=1; f1<origCollisionIntensityMatrix.length; f1++)
			if (origCollisionIntensityMatrix[f1] != null)
				for (int f2=0; f2<f1; f2++)
					if (origCollisionIntensityMatrix[f1][f2] != 0f
					 && mTorsionSetStrategy.getBondsBetweenFragments(f1, f2).length != 2)
						remainingCollisionIntensity += origCollisionIntensityMatrix[f1][f2];

		// if other collision strength is already inacceptable
		if (remainingCollisionIntensity > TorsionSetStrategy.MAX_ALLOWED_COLLISION_INTENSITY)
			return false;

		boolean changeDone = false;
		Conformer conformer = null;
		Conformer backup = null;
		double origCollisionIntensitySum = torsionSet.getCollisionIntensitySum();
		for (int f1=1; f1<origCollisionIntensityMatrix.length; f1++) {
			if (origCollisionIntensityMatrix[f1] != null) {
				for (int f2=0; f2<f1; f2++) {
					if (origCollisionIntensityMatrix[f1][f2] > MIN_ESCAPE_GAIN_PER_STEP
					 && mTorsionSetStrategy.getBondsBetweenFragments(f1, f2).length == 2) {
						int[] rotatableBond = mTorsionSetStrategy.getBondsBetweenFragments(f1, f2);
						int angle = (mRandom.nextDouble() < 0.5) ? -ESCAPE_ANGLE : ESCAPE_ANGLE;
						double origCollisionIntensity = origCollisionIntensityMatrix[f1][f2];
						for (int i=1; i<=ESCAPE_STEPS; i++) {
							// make sure we keep the original coordinates
							if (conformer == null)
								conformer = new Conformer(origConformer);
							else if (backup == null)
								backup = new Conformer(conformer);
							else
								backup.copyFrom(conformer);

							for (int bond : rotatableBond)
								mRotatableBond[bond].rotateTo(conformer, (short) (conformer.getBondTorsion(bond) + i*angle));

							double localCollisionIntensity = calculateCollisionIntensity(conformer, mRigidFragment[f1], mRigidFragment[f2]);
							if (localCollisionIntensity < origCollisionIntensity - MIN_ESCAPE_GAIN_PER_STEP) {
								origCollisionIntensity = localCollisionIntensity;
								changeDone = true;

								// not enough collision intensity left for correction
								if (localCollisionIntensity < MIN_ESCAPE_GAIN_PER_STEP)
									break;
								}
							else {
								if (backup != null)
									conformer.copyFrom(backup);
								else
									conformer.copyFrom(origConformer);
								break;
								}
							}
						}
					}
				}
			}

		if (!changeDone)
			return false;

		checkCollision(conformer, torsionSet);
		if (torsionSet.getCollisionIntensitySum() >= origCollisionIntensitySum) {
			// local improvement causes worsening somewhere else
			torsionSet.setCollisionIntensity(origCollisionIntensitySum, origCollisionIntensityMatrix);
			return false;
			}

		origConformer.copyFrom(conformer);
		return true;
		}

	private double calculateCollisionIntensity(Conformer conformer, Rigid3DFragment f1, Rigid3DFragment f2) {
		double collisionIntensitySum = 0;
		for (int i=0; i<f1.getCoreSize(); i++) {
			int atom1 = f1.coreToOriginalAtom(i);
			double vdwr1 = getToleratedVDWRadius(mMolecule.getAtomicNo(atom1));
			for (int j=0; j<f2.getCoreSize(); j++) {
				int atom2 = f2.coreToOriginalAtom(j);
				double minDistance = vdwr1+getToleratedVDWRadius(mMolecule.getAtomicNo(atom2));
				double dx = Math.abs(conformer.getX(atom1) - conformer.getX(atom2));
				if (dx < minDistance) {
					double dy = Math.abs(conformer.getY(atom1) - conformer.getY(atom2));
					if (dy < minDistance) {
						double dz = Math.abs(conformer.getZ(atom1) - conformer.getZ(atom2));
						if (dz < minDistance) {
							double distance = Math.sqrt(dx*dx+dy*dy+dz*dz);
							if (distance < minDistance) {
								double relativeCollision = (minDistance - distance) / minDistance;
								double collisionIntensity = relativeCollision * relativeCollision;
								collisionIntensitySum += collisionIntensity;
								}
							}
						}
					}
				}
			}
		return collisionIntensitySum;
		}
	}