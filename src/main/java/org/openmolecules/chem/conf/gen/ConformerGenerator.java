/*
 * Copyright 2013-2020 Thomas Sander, openmolecules.org
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.gen;

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.chem.conf.VDWRadii;
import com.actelion.research.util.DoubleFormat;
import com.actelion.research.util.IntArrayComparator;
import org.openmolecules.chem.conf.so.ConformationSelfOrganizer;
import org.openmolecules.chem.conf.so.SelfOrganizedConformer;

import java.util.Arrays;
import java.util.Random;
import java.util.TreeMap;

/**
 * This class generates 3D-conformers of a given molecule using the following strategy:
 * <li>All rotatable, non-ring bonds are determined.
 * <li>Fragments separated by rotatable bonds are considered rigid, but there may be more than one possible
 * fragment conformer, e.g. chair- and boat conformers of a saturated 6-membered ring.
 * <li>These atom coordinate sets of rigid fragments are handed out by a dedicated RigidFragmentProvider instance,
 * which either generates them using a self organization algorithm, or which takes it from a cache.
 * <li>For every rotatable bond a list of preferred torsion angles is determined based on from a COD statistics
 * of similar bond environments.
 * <li>For individual torsion values likelihoods are estimated based on frequency and atom collisions of vicinal fragments.
 * <li>A dedicated (systematic, biased or random) torsion set strategy delivers collision-free torsion sets, i.e. conformers.
 * <br><br>
 * For generating conformers in multiple threads, every thread needs its own ConformerGenerator instance.
 * If they use a RigidFragmentCache, then the cache is shared among all ConformerGenerators.<br>
 * Important: Input molecules should contain absolute stereo centers. If they contain undefined or ESR type '&' or 'or'
 * stereo centers, then a ConformerGenerator randomly takes one of the possible stereo isomers and generates conformers
 * for that. If you want conformers for all possible stereo isomers of a molecules with non-absolute stereo centers,
 * you should use a StereoIsomerEnumerator to produce all possible stereo isomers and then produce conformers for every
 * one of them. If half of a set of stereo isomers consists of the enantiomers of the other half, then it is advisable
 * to generate conformes for one half only and to generate the second half by just mirroring the first halfs coordinates.
 * To do that use option skipEnantiomers==true create a mirrored set of conformers, if isSkippingEnantiomers() of the
 * StereoIsomerEnumerator returns true.
 */
public class ConformerGenerator {
	public static final int STRATEGY_LIKELY_SYSTEMATIC = 1;
	public static final int STRATEGY_PURE_RANDOM = 2;
	public static final int STRATEGY_LIKELY_RANDOM = 3;
	public static final int STRATEGY_ADAPTIVE_RANDOM = 4;

	protected static final double VDW_TOLERANCE_HYDROGEN = 0.90;  // factor on VDW radii for minimum tolerated non bound atom distances
	protected static final double VDW_TOLERANCE_OTHER = 0.90;     // factor on VDW radii for minimum tolerated non bound atom distances

	private static final int ESCAPE_ANGLE = 8;  // degrees to rotate two adjacent rotatable bonds to escape collisions
	private static final int ESCAPE_STEPS = 3;	// how often we apply this rotation trying to solve the collision
	private static final double MIN_ESCAPE_GAIN_PER_STEP = 1.0;

	// We try to translate arbitrary collision values into a kcal/mol energy scale
	public static final double COLLISION_STRAIN_TO_ENERGY_FACTOR = 20;

	private StereoMolecule mMolecule;
	private TreeMap<int[],BaseConformer> mBaseConformerMap;
	private RotatableBond[] mRotatableBond;
	private RigidFragment[] mRigidFragment;
	private ConformationSelfOrganizer mSelfOrganizer;
	private final RigidFragmentProvider mRigidFragmentProvider;
	private TorsionSetStrategy mTorsionSetStrategy;
	private TorsionSet mTorsionSet;
	private final long mRandomSeed;
	private long mTimeOut,mStopMillis;
	private int mDisconnectedFragmentCount,mAllConformerCount,mReturnedConformerCount;
	private boolean mUseSelfOrganizerIfAllFails,mIsDiagnosticsMode,mIsFinished;
	private int[] mFragmentNo,mDisconnectedFragmentNo,mDisconnectedFragmentSize;
	private boolean[][] mSkipCollisionCheck;
	private final Random mRandom;
	private ThreadMaster mThreadMaster;
	private ConformerSetDiagnostics mDiagnostics;

	/**
	 * Assuming that the given molecule has 2D-coordinates, this method
	 * converts all implicit hydrogen atoms into explicit ones by filling valences
	 * and adapting for atom charges. New hydrogen atoms receive new 2D-coordinates
	 * by equally locating them between those two neighbors with the widest angle between
	 * their bonds. Any stereo configurations deducible from 2D-coordinates are retained.
	 * @param mol
	 */
	public static void addHydrogenAtoms(StereoMolecule mol) {
		// We may have parities but empty coordinates. In this case we need to protect parities.
		int oldStereoHelperBits = mol.getHelperArrayStatus() & Molecule.cHelperBitsStereo;

		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		int[] implicitHydrogen = new int[mol.getAtoms()];
		for (int atom=0; atom<mol.getAtoms(); atom++)
			implicitHydrogen[atom] = mol.getImplicitHydrogens(atom);

		double hydrogenBondLength = 0.8 * mol.getAverageBondLength();

		for (int atom=0; atom<implicitHydrogen.length; atom++)
			if (implicitHydrogen[atom] != 0)
				for (int i=0; i<implicitHydrogen[atom]; i++)
					mol.addBond(atom, mol.addAtom(1), Molecule.cBondTypeSingle);

		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		for (int atom=0; atom<implicitHydrogen.length; atom++)
			if (implicitHydrogen[atom] != 0)
				setHydrogenLocations(mol, atom, implicitHydrogen[atom], hydrogenBondLength);

		// addAtom() and addBond() clear the helper status, i.e. flag all helpers as invalid.
		// Adding hydrogens does not destroy parities. Though, we may flag them to be valid again.
		if ((oldStereoHelperBits & Molecule.cHelperBitParities) != 0)
			mol.setParitiesValid(oldStereoHelperBits);
		}

	/**
	 * Finds the widest open angle between all connected non-stereo bonds of atom, divides this angle
	 * into hydrogenCount+1 equal parts and sets atom coordinates of hydrogenCount new hydrogen atoms
	 * such, that they equally occupy the space and not collide with a potential stereo-bonded neighbour.
	 * @param mol
	 * @param atom
	 * @param newHydrogenCount new hydrogen atoms added to atom
	 * @param avbl
	 */
	private static void setHydrogenLocations(StereoMolecule mol, int atom, int newHydrogenCount, double avbl) {
		int firstNewHydrogenConnIndex = mol.getAllConnAtoms(atom) - newHydrogenCount;

		int stereoBondIndex = -1;
		double stereoBondAngle = Double.NaN;
		for (int i=0; i<firstNewHydrogenConnIndex; i++) {
			if (mol.isStereoBond(mol.getConnBond(atom, i), atom)) {
				stereoBondIndex = i;
				stereoBondAngle = mol.getBondAngle(atom, mol.getConnAtom(atom, i));
				break;
				}
			}
		int stereoBondCount = (stereoBondIndex == -1) ? 0 : 1;

		double[] angle = null;
		if (stereoBondCount < firstNewHydrogenConnIndex) {
			angle = new double[firstNewHydrogenConnIndex-stereoBondCount];
			int bond = 0;
			for (int i=0; i<firstNewHydrogenConnIndex; i++)
				if (i != stereoBondIndex)
					angle[bond++] = mol.getBondAngle(atom, mol.getConnAtom(atom, i));
	
			Arrays.sort(angle);
			}

		double angleIncrement = 2.0*Math.PI/Math.max(newHydrogenCount, 3);
		double startAngle = 0.0;

		if (angle != null) {
			if (angle.length == 1 && newHydrogenCount == 1) {
				startAngle = angle[0];
				}
			else {
				double biggestAngleDif = 0.0;
				for (int i=0; i<angle.length; i++) {
					double a1 = (i == 0) ? angle[angle.length - 1] - Math.PI * 2.0 : angle[i - 1];
					double a2 = angle[i];
					if (biggestAngleDif < a2 - a1) {
						biggestAngleDif = a2 - a1;
						startAngle = a1;
						}
					}
				angleIncrement = biggestAngleDif / (newHydrogenCount + 1);
				}
			}

		for (int i=0; i<newHydrogenCount; i++) {
			startAngle += angleIncrement;
			double hydrogenAngle = startAngle;
			if (stereoBondCount != 0) {
				double dif = Molecule.getAngleDif(hydrogenAngle, stereoBondAngle);
				if (Math.abs(dif) < angleIncrement/2) {
					if (dif < 0)
						hydrogenAngle = stereoBondAngle - angleIncrement/2;
					else
						hydrogenAngle = stereoBondAngle + angleIncrement/2;
					}
				}
			int newHydrogen = mol.getConnAtom(atom, firstNewHydrogenConnIndex+i);
			mol.setAtomX(newHydrogen, mol.getAtomX(atom) + avbl * Math.sin(hydrogenAngle));
			mol.setAtomY(newHydrogen, mol.getAtomY(atom) + avbl * Math.cos(hydrogenAngle));
			}
		}

	public static double getToleratedVDWRadius(int atomicNo) {
		float vdwr = VDWRadii.VDW_RADIUS[atomicNo < VDWRadii.VDW_RADIUS.length ? atomicNo : 6];
		return vdwr * (atomicNo == 1 ? VDW_TOLERANCE_HYDROGEN : VDW_TOLERANCE_OTHER);
		}

	/**
	 * Instantiates a ConformerGenerator for creating non-reproducible conformers.
	 * Uses the default rigid fragment cache, which is initially empty unless
	 * RigidFragmentCache.getDefaultInstance().loadDefaultCache() was called before.
	 * While rigid fragments, when found in the cache are always energy  minimized,
	 * the ones not found in the cache are generated by self organization and not minimized afterawrds.
	 */
	public ConformerGenerator() {
		this(0L, RigidFragmentCache.getDefaultInstance(), false);
		}

	/**
	 * Instantiates a ConformerGenerator for creating reproducible conformers unless seed is 0L.
	 * Uses the default rigid fragment cache, which is initially empty unless
	 * RigidFragmentCache.getDefaultInstance().loadDefaultCache() was called before.
	 * While rigid fragments, when found in the cache, are always energy  minimized,
	 * the ones not found in the cache are generated by self organization and energy
	 * minimized before cached and being used, if optimizeRigidFragments is true.
	 * @param seed != 0L if conformers shall be created in a reproducible way
	 * @param optimizeRigidFragments if true, then all rigid fragments will be energy minimized by the MMFF94s+ forcefield
	 */
	public ConformerGenerator(long seed, boolean optimizeRigidFragments) {
		this(seed, RigidFragmentCache.getDefaultInstance(), optimizeRigidFragments);
	}

	/**
	 * Instantiates a ConformerGenerator for creating non-reproducible conformers.
	 * Uses the default rigid fragment cache, which is initially empty unless
	 * RigidFragmentCache.getDefaultInstance().loadDefaultCache() was called before.
	 * While rigid fragments, when found in the cache, are always energy minimized,
	 * the ones not found in the cache are generated by self organization and energy
	 * minimized before cached and being used, if optimizeRigidFragments is true.
	 * @param optimizeRigidFragments if true, then all rigid fragments will be energy minimized by the MMFF94s+ forcefield
	 */
	public ConformerGenerator(boolean optimizeRigidFragments) {
		this(0L, RigidFragmentCache.getDefaultInstance(), optimizeRigidFragments);
	}

	/**
	 * Instantiates a ConformerGenerator for creating reproducible conformers unless seed is 0L.
	 * Uses a custom rigid fragment cache or no cache at all.
	 * Rigid fragments, when found in the cache, are always energy  minimized.
	 * The ones not found in the cache are generated by self organization. They are energy
	 * minimized before being cached and being used, if optimizeRigidFragments is true.
	 * @param seed != 0L if conformers shall be created in a reproducible way
	 * @param cache may be null for generating all rigid fragment conformers on the fly
	 * @param optimizeRigidFragments if true, then all rigid fragments will be energy minimized by the MMFF94s+ forcefield
	 */
	public ConformerGenerator(long seed, RigidFragmentCache cache, boolean optimizeRigidFragments) {
		TorsionDB.initialize(TorsionDB.MODE_ANGLES);
		mRandomSeed = seed;
		mRandom = (seed == 0) ? new Random() : new Random(seed);
		mRigidFragmentProvider = new RigidFragmentProvider(seed, cache, optimizeRigidFragments);
	}

	/**
	 * Instantiates a ConformerGenerator for creating reproducible conformers unless seed is 0L.
	 * Uses a custom RigidFragmentProvider, which typically should be done if your fragments shall be
	 * energy minimized with something else than the MMFF94s+ forcefield, e.g. a QM method.
	 * @param seed != 0L if conformers shall be created in a reproducible way
	 * @param rfp
	 */
	public ConformerGenerator(long seed, RigidFragmentProvider rfp) {
		TorsionDB.initialize(TorsionDB.MODE_ANGLES);
		mRandomSeed = seed;
		mRandom = (seed == 0) ? new Random() : new Random(seed);
		mRigidFragmentProvider = rfp;
	}

	/**
	 * In diagnostic mode the ConformerGenerator records statistics about rigid conformer indexes and
	 * torsion indexes used when building conformers. It also stores rigid conformer and torsion likelyhoods
	 * applied when building conformers. Diagnostics information for the most recently used molecule can
	 * be retrieved after generating conformer with the getDiagnostics() methods.
	 * @param b
	 */
	public void setDiagnosticMode(boolean b) {
		mIsDiagnosticsMode = b;
	}

	public ConformerSetDiagnostics getDiagnostics() {
		return mDiagnostics;
	}

	/**
	 * If the conformer generation shall be stopped after a certain time
	 * of unsuccessfully trying generating conformers, then set a timeout.
	 * @param millis milli-seconds after which to stop the conformer generation or 0L to remove the timeout
	 */
	public void setTimeOut(long millis) {
		mTimeOut = millis;
		applyTimeOut();
	}

	public void applyTimeOut() {
		mStopMillis = (mTimeOut == 0L) ? 0L : System.currentTimeMillis() + mTimeOut;
		mRigidFragmentProvider.setStopTime(mStopMillis);
	}

	/**
	 * If the conformer generation must be stopped from outside, for instance because of user
	 * intervention or because of a defined timeout, then provide a ThreadMaster with this method.
	 * @param tm
	 */
	public void setThreadMaster(ThreadMaster tm) {
		mThreadMaster = tm;
		mRigidFragmentProvider.setThreadMaster(tm);
	}

	private boolean mustStop() {
		return (mThreadMaster != null && mThreadMaster.threadMustDie())
			|| (mStopMillis != 0 && System.currentTimeMillis() > mStopMillis);
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
		if (!initialize(mol, false))
			return null;

		if (mRotatableBond != null) {
			mTorsionSetStrategy = new TorsionSetStrategyAdaptiveRandom(this, true, true, mRandomSeed);
			mTorsionSetStrategy.setMaxTotalCount(400);
			mBaseConformerMap = new TreeMap<>(new IntArrayComparator());
			return getNextConformer();
			}
		else {
			ConformationSelfOrganizer sampler = new ConformationSelfOrganizer(mol, true);
			sampler.setThreadMaster(mThreadMaster);
			sampler.setStopTime(mStopMillis);
			Conformer conformer = sampler.generateOneConformer(mRandomSeed);
			separateDisconnectedFragments(conformer);
			conformer.setName("SO#1");
			return conformer;
			}
		}

		/**
		 * Don't call this method directly. Rather call initializeConformers() !!!
		 * Adds implicit hydrogens to the molecule and determines all rotatable bonds,
		 * which are not part of a ring. Generates rigid fragments between rotatable bonds.
		 * The base conformer is constructed by connecting all rigid fragments using the
		 * most frequent torsions. Atoms of different fragments in the base conformer may
		 * collide. In order to obtain collision free conformers choose a TorsionStrategy
		 * and call getNextConformer() at least once.
		 * @param mol
		 * @param use60degreeSteps use 60 degree steps for every rotatable bond instead of torsion DB
		 */
	protected boolean initialize(StereoMolecule mol, boolean use60degreeSteps) {
		applyTimeOut();
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
		mAllConformerCount = 0;
		mReturnedConformerCount = 0;
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

		if (!locateInitialFragments(isRotatableBond))
			return false;

		mRotatableBond = new RotatableBond[count];
		int rotatableBond = 0;
		for (int bond=0; bond<mol.getBonds(); bond++)
			if (isRotatableBond[bond])
				mRotatableBond[rotatableBond++] = new RotatableBond(mol, bond, mFragmentNo,
						mDisconnectedFragmentNo, mDisconnectedFragmentSize[mDisconnectedFragmentNo[mol.getBondAtom(0, bond)]],
						mRigidFragment, use60degreeSteps);

		// sort by descending atom count of smaller side, i.e. we want those bond dividing into equal parts first!
		Arrays.sort(mRotatableBond, (b1, b2) ->
			Integer.compare(b2.getSmallerSideAtoms().length, b1.getSmallerSideAtoms().length) );

		if (mIsDiagnosticsMode)
			mDiagnostics = new ConformerSetDiagnostics(this);

		// TODO this is actually only used with random conformers. Using a conformer mode
		// may save some time, if systematic conformers are created.
		initializeCollisionCheck();

		return true;
		}

	/**
	 * After calling initializeConformers() this method returns the number of rotatable bonds,
	 * which are used to separate the molecule into rigid fragments.
	 * @return
	 */
	public int getRotatableBondCount() {
		return mRotatableBond == null ? 0 : mRotatableBond.length;
		}

	public TreeMap<int[],BaseConformer> getBaseConformerMap() {
		return mBaseConformerMap;
		}

	public TorsionSetStrategy getTorsionSetStrategy() {
		return mTorsionSetStrategy;
		}

	public BaseConformer getBaseConformer(int[] fragmentPermutation) {
		BaseConformer baseConformer = mBaseConformerMap.get(fragmentPermutation);
		if (baseConformer != null)
			return baseConformer;

		baseConformer = new BaseConformer(mMolecule, mRigidFragment, mRotatableBond, fragmentPermutation, mRandom);
//printDebugConformers(baseConformer);
		mBaseConformerMap.put(fragmentPermutation.clone(), baseConformer);
		return baseConformer;
		}

	/**
	 * Creates the next random, likely or systematic new(!) conformer of the molecule
	 * that was passed when calling initializeConformers(). A new conformer is one,
	 * whose combination of torsion angles was not used in a previous conformer
	 * created by this function since the last call of initializeConformers().
	 * If parameter mol is null, then a compact copy of the original molecule with the new
	 * coordinates is returned. You may pass the original molecule or a copy of it to
	 * recycle the original molecule to receive new 3D coordinates.
	 * Every call of this method creates a new collision-free conformer until the employed torsion set
	 * strategy decides that it cannot generate any more suitable torsion sets.
	 * @param mol null or molecule used during initialization or a copy of it
	 * @return conformer or null, if all/maximum torsion permutations have been tried
	 */
	public StereoMolecule getNextConformerAsMolecule(StereoMolecule mol) {
		Conformer conformer = getNextConformer();
		return (conformer == null) ? null : conformer.toMolecule(mol);
		}

	public Conformer getNextConformer() {
		return getNextConformer(null);
	}

	/**
	 * Creates the next random, likely or systematic new(!) conformer of the molecule
	 * that was passed when calling initializeConformers(). A new conformer is one,
	 * whose combination of torsion angles was not used in a previous conformer
	 * created by this function since the last call of initializeConformers().
	 * Every call of this method creates a new collision-free conformer until the employed torsion set
	 * strategy decides that it cannot generate any more suitable torsion sets.
	 * @param torsionSetHolder : will contain the TorsionSet that gave rise to the conformer.
	 * @return conformer or null, if all/maximum torsion permutations have been tried
	 */
	public Conformer getNextConformer(TorsionSet[] torsionSetHolder) {
		if (mRotatableBond == null && mSelfOrganizer == null)
			return null;

		if (mSelfOrganizer != null) {
			SelfOrganizedConformer conformer = mSelfOrganizer.getNextConformer();
			if (conformer != null) {
				separateDisconnectedFragments(conformer);
				mReturnedConformerCount++;
				conformer.setName("SO#"+(++mAllConformerCount));
				return conformer;
				}

			return null;
			}

		if (mIsFinished)
			return null;

		// create a base conformer from first set of fragments and calculate torsion likelihoods
		if (mBaseConformerMap.isEmpty())
			getBaseConformer(new int[mRigidFragment.length]);

		mTorsionSet = mTorsionSetStrategy.getNextTorsionSet(mTorsionSet, mDiagnostics);
		while (mTorsionSet != null && !mustStop()) {
			BaseConformer baseConformer = getBaseConformer(mTorsionSet.getConformerIndexes());
			if (mTorsionSet.getConformer() == null) {
				mTorsionSet.setConformer(baseConformer.deriveConformer(mTorsionSet.getTorsionIndexes(), "#" + (++mAllConformerCount)));
				if (mIsDiagnosticsMode)
					mDiagnostics.addNew(mTorsionSet);
				}

			// If the torsionSet has already a collision value, then it is a second choice torsion set
			// with a collision value below the tolerance. No need to check or fix again.
			if (mTorsionSet.getCollisionStrainSum() == 0.0) {
				calculateCollision(mTorsionSet, mTorsionSet.getConformer());

				if (mTorsionSet.getCollisionStrainSum() != 0.0)
					tryFixCollisions(baseConformer, mTorsionSet);

				if (mIsDiagnosticsMode) {
					mDiagnostics.get(mTorsionSet).setConformer(baseConformer, mTorsionSet.getConformer());
					mDiagnostics.get(mTorsionSet).setCollisionStrain(mTorsionSet.getCollisionStrainSum());
					}
				}

			if (mTorsionSet.getCollisionStrainSum() > mTorsionSetStrategy.calculateCollisionTolerance()) {
				mTorsionSet = mTorsionSetStrategy.getNextTorsionSet(mTorsionSet, mDiagnostics);
				if (mTorsionSet != null || mReturnedConformerCount != 0)
					continue;

				if (mUseSelfOrganizerIfAllFails) {
					// We couldn't create torsion strategy based conformers: switch to self organizer!
					mSelfOrganizer = new ConformationSelfOrganizer(mMolecule, true);
					mSelfOrganizer.setThreadMaster(mThreadMaster);
					mSelfOrganizer.setStopTime(mStopMillis);
					mSelfOrganizer.initializeConformers(mRandomSeed, -1);
					SelfOrganizedConformer conformer = mSelfOrganizer.getNextConformer();
					if (conformer != null) {
						separateDisconnectedFragments(conformer);
						mReturnedConformerCount++;
						conformer.setName("SO#"+(++mAllConformerCount));
						return conformer;
						}
					}

				// we didn't get any torsion set that didn't collide; take the best we had
				mTorsionSet = mTorsionSetStrategy.getBestCollidingTorsionIndexes();

				mIsFinished = true;	// we are finished with conformers
				}

			if (mTorsionSet != null) {
				separateDisconnectedFragments(mTorsionSet.getConformer());
				mTorsionSet.setUsed();
				mReturnedConformerCount++;

				if(torsionSetHolder != null)
					torsionSetHolder[0] = new TorsionSet(mTorsionSet);

				if (mIsDiagnosticsMode)
					mDiagnostics.get(mTorsionSet).setSuccess(true);

				return mTorsionSet.getConformer();
				}
// This poses the danger to produce highly strained unrealistic conformers, e.g. for impossible compounds...
//			else if (mReturnedConformerCount == 0) {
//				// We couldn't create torsion strategy based conformers: try creating one conformer with self organizer!
//				ConformationSelfOrganizer sampler = new ConformationSelfOrganizer(mMolecule, true);
//	    		sampler.setThreadMaster(mThreadMaster);
//  			sampler.setStopTime(mStopMillis);
//				Conformer conformer = sampler.generateOneConformer(mRandomSeed);
//				separateDisconnectedFragments(conformer);
//				mReturnedConformerCount++;
//				conformer.setName("SO#1");
//				return conformer;
//				}
			}

		return null;
		}

/*private void printDebugConformers(BaseConformer bc) {
 int degreesA = bc.getBondTorsion(mRotatableBond[0].getBond());
 int degreesB = bc.getBondTorsion(mRotatableBond[1].getBond());
 System.out.println("Startangles A:"+degreesA+" B:"+degreesB);
 System.out.println("idcode\tidcoords\tdegrees A\tdegrees B");
 for (int degrees=0; degrees<100; degrees+=20) {
  Conformer co = new Conformer(bc);
  bc.rotateTo(co, mRotatableBond[0], (short)(degreesA + degrees));
  Canonizer ca = new Canonizer(co.toMolecule(null));
  System.out.println(ca.getIDCode() + "\t" + ca.getEncodedCoordinates() + "\t" + degrees +"\t0");
 }
 for (int degrees=0; degrees<100; degrees+=20) {
  Conformer co = new Conformer(bc);
  bc.rotateTo(co, mRotatableBond[1], (short)(degreesB + degrees));
  Canonizer ca = new Canonizer(co.toMolecule(null));
  System.out.println(ca.getIDCode() + "\t" + ca.getEncodedCoordinates() + "\t0\t" + degrees);
 }
}*/

	/**
	 * @return count of valid delivered conformers
	 */
	public int getConformerCount() {
		return mReturnedConformerCount;
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
			mSelfOrganizer.setThreadMaster(mThreadMaster);
			mSelfOrganizer.setStopTime(mStopMillis);
			mSelfOrganizer.initializeConformers(mRandomSeed, -1);
			}
		else {
			mBaseConformerMap = new TreeMap<>(new IntArrayComparator());
			switch(strategy) {
			case STRATEGY_PURE_RANDOM:
				mTorsionSetStrategy = new TorsionSetStrategyRandom(this, false, mRandomSeed);
				break;
			case STRATEGY_LIKELY_RANDOM:
				mTorsionSetStrategy = new TorsionSetStrategyRandom(this, true, mRandomSeed);
				break;
			case STRATEGY_ADAPTIVE_RANDOM:
				mTorsionSetStrategy = new TorsionSetStrategyAdaptiveRandom(this, true, true, mRandomSeed);
				break;
			case STRATEGY_LIKELY_SYSTEMATIC:
				mTorsionSetStrategy = new TorsionSetStrategyLikelySystematic(this);
				break;
				}
			mTorsionSetStrategy.setMaxTotalCount(maxTorsionSets);
			}

		return true;
		}

	public RotatableBond[] getRotatableBonds() {
		return mRotatableBond;
		}

	public RigidFragment[] getRigidFragments() {
		return mRigidFragment;
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

	private boolean locateInitialFragments(boolean[] isRotatableBond) {
		mFragmentNo = new int[mMolecule.getAllAtoms()];
		int fragmentCount = mMolecule.getFragmentNumbers(mFragmentNo, isRotatableBond, true);
		mRigidFragment = new RigidFragment[fragmentCount];
		boolean ok = true;
		for (int i=0; i<fragmentCount; i++)
			ok &= ((mRigidFragment[i] = mRigidFragmentProvider.createFragment(mMolecule, mFragmentNo, i)) != null);
		return ok;
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
		for (RigidFragment rf:mRigidFragment)
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

	private boolean calculateCollision(TorsionSet torsionSet, Conformer conformer) {
		if (mIsDiagnosticsMode)
			mDiagnostics.get(torsionSet).setCollisionAtoms(null);

		double collisionStrainSum = 0;
		double[][] collisionStrainMatrix = null;
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
									double collisionStrain = COLLISION_STRAIN_TO_ENERGY_FACTOR * relativeCollision * relativeCollision;
									collisionStrainSum += collisionStrain;

									if (mIsDiagnosticsMode) {
										mDiagnostics.get(torsionSet).writeCollisionLog("a1:" + atom1 + " f1:" + mFragmentNo[atom1] + " a2:" + atom2 + " f2:" + mFragmentNo[atom2] + " distance:" + DoubleFormat.toString(distance) + " min:" + DoubleFormat.toString(minDistance));
										if (mDiagnostics.get(torsionSet).getCollisionAtoms() == null) {
											int[] atoms = new int[2];
											atoms[0] = atom1;
											atoms[1] = atom2;
											mDiagnostics.get(torsionSet).setCollisionAtoms(atoms);
											}
										}

									if (collisionStrainMatrix == null)
										collisionStrainMatrix = new double[mRigidFragment.length][];
									int f1 = mFragmentNo[atom1];
									int f2 = mFragmentNo[atom2];
									if (f1 < f2) {
										if (collisionStrainMatrix[f2] == null)
											collisionStrainMatrix[f2] = new double[f2];
										collisionStrainMatrix[f2][f1] += collisionStrain;
										}
									else {
										if (collisionStrainMatrix[f1] == null)
											collisionStrainMatrix[f1] = new double[f1];
										collisionStrainMatrix[f1][f2] += collisionStrain;
										}
									}
								}
							}
						}
					}
				}
			}
		torsionSet.setCollisionStrain(collisionStrainSum, collisionStrainMatrix);
		return (collisionStrainSum != 0);
		}

	/**
	 * If we have collisions between two fragments that have two rotatable bonds in between,
	 * then we try to rotate both rotatable bonds in any of these situations stepwise
	 * by a small angle to reduce the collision strain sum. If we achieve a sufficient
	 * improvement without hampering the rest of the molecule, the passed conformer is modified.
	 * @param torsionSet
	 * @return true if successful; false if no fix possible that reduces total collision strain below tolerance
	 */
	private boolean tryFixCollisions(BaseConformer baseConformer, TorsionSet torsionSet) {
		if (mustStop())
			return false;

		double[][] origCollisionStrainMatrix = torsionSet.getCollisionStrainMatrix();
		double remainingCollisionStrain = 0;
		for (int f1=1; f1<origCollisionStrainMatrix.length; f1++)
			if (origCollisionStrainMatrix[f1] != null)
				for (int f2=0; f2<f1; f2++)
					if (origCollisionStrainMatrix[f1][f2] != 0f
					 && mTorsionSetStrategy.getBondsBetweenFragments(f1, f2).length != 2)
						remainingCollisionStrain += origCollisionStrainMatrix[f1][f2];

		// if other collision strength is already unacceptable
		if (remainingCollisionStrain > TorsionSetStrategy.MAX_COLLISION_STRAIN)
			return false;

		Conformer conformer = null;
		Conformer localBest = null;
		Conformer globalBest = null;
		double bestCollisionStrain = Double.MAX_VALUE;
		double origCollisionStrainSum = torsionSet.getCollisionStrainSum();
		for (int f1=1; f1<origCollisionStrainMatrix.length; f1++) {
			if (origCollisionStrainMatrix[f1] != null) {
				for (int f2=0; f2<f1; f2++) {
					if (origCollisionStrainMatrix[f1][f2] > MIN_ESCAPE_GAIN_PER_STEP
					 && mTorsionSetStrategy.getBondsBetweenFragments(f1, f2).length == 2) {
						int[] rbIndex = mTorsionSetStrategy.getBondsBetweenFragments(f1, f2);
						int[] startTorsion = new int[2];
						for (int i=0; i<2; i++) {
							RotatableBond rotatableBond = mRotatableBond[rbIndex[i]];
							startTorsion[i] = torsionSet.getConformer().getBondTorsion(rotatableBond.getBond());
//							short[] torsionRange = rotatableBond.getDefaultTorsionRanges()[torsionSet.getTorsionIndexes()[rbIndex[i]]];
							}

						int[] torsionDif = new int[2];
						for (int r1=-1; r1<=1; r1+=2) {
							torsionDif[0] = r1 * ESCAPE_ANGLE;
							for (int r2=-1; r2<=1; r2+=2) {
								torsionDif[1] = r2 * ESCAPE_ANGLE;

								double collisionStrain = origCollisionStrainMatrix[f1][f2];
								for (int step=1; step<=ESCAPE_STEPS; step++) {
									if (conformer == null)
										conformer = new Conformer(torsionSet.getConformer());
									else
										conformer.copyFrom(torsionSet.getConformer());

									for (int i=0; i<2; i++)
										baseConformer.rotateTo(conformer, mRotatableBond[rbIndex[i]], (short)(startTorsion[i] + step * torsionDif[i]));

									double strain = calculateCollisionStrain(conformer, mRigidFragment[f1], mRigidFragment[f2]);
									if (strain < collisionStrain - MIN_ESCAPE_GAIN_PER_STEP) {
										collisionStrain = strain;
										if (localBest == null)
											localBest = new Conformer(conformer);
										else
											localBest.copyFrom(conformer);

										// not enough collision strain left for correction
										if (collisionStrain < MIN_ESCAPE_GAIN_PER_STEP)
											break;
										}
									else {
										break;
										}
									}
								if (collisionStrain < origCollisionStrainMatrix[f1][f2]
								 && collisionStrain < bestCollisionStrain) {
									bestCollisionStrain = collisionStrain;
									if (globalBest == null)
										globalBest = new Conformer(localBest);
									else
										globalBest.copyFrom(localBest);
									}
								}
							}
						}
					}
				}
			}

		if (globalBest == null)
			return false;

		calculateCollision(torsionSet, globalBest);
		if (torsionSet.getCollisionStrainSum() >= origCollisionStrainSum) {
			// local improvement causes worsening somewhere else
			torsionSet.setCollisionStrain(origCollisionStrainSum, origCollisionStrainMatrix);
			return false;
			}

		torsionSet.getConformer().copyFrom(globalBest);

		return true;
		}

	private double calculateCollisionStrain(Conformer conformer, RigidFragment f1, RigidFragment f2) {
		double collisionStrainSum = 0;
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
								double collisionStrain = relativeCollision * relativeCollision;
								collisionStrainSum += collisionStrain;
								}
							}
						}
					}
				}
			}
		return COLLISION_STRAIN_TO_ENERGY_FACTOR * collisionStrainSum;
		}


	/**
	 * Computes a conformer from a TorsionSet object.
 	 * @param torsionset
	 * @return the generated conformer.
	 *
	public Conformer generateConformerFromTorsionSet(TorsionSet torsionset) {
		return generateConformerFromTorsionSet(torsionset.getConformerIndexes(), torsionset.getTorsionIndexes());
	}

	private Conformer generateConformerFromTorsionSet(int[] conformerIndexes, int[] torsionIndexes) {
		Conformer conformer = getBaseConformer(conformerIndexes).deriveConformer(torsionIndexes, "#"+(++mAllConformerCount));
		separateDisconnectedFragments(conformer);
		return conformer;
	}*/
}
