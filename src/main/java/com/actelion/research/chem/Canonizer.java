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
* 3. Neither the name of the copyright holder nor the
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

// Restriction: - Although bond query features are encoded into the idcode, idcodes of
//				fragments with bond query features are not necessarily unique.
//			  - Atom query features are(!) considered during atom priority assignment
//				and therefore fragments with atom query features get perfectly unique
//				idcode. The only exception are atoms with assigned atom lists, which
//				are encoded into idcodes but may result in non-unique idcode.

package com.actelion.research.chem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

public class Canonizer {
	public static final int CREATE_SYMMETRY_RANK = 1;

	// The following options CONSIDER_DIASTEREOTOPICITY, CONSIDER_ENANTIOTOPICITY
	// and CONSIDER_STEREOHETEROTOPICITY have no influence on the idcode,
	// i.e. the idcode is the same whether or not one of these options is
	// used. However, if you require e.g. a pro-E atom always to appear
	// before the pro-Z, e.g. because you can distinguish them from the
	// encoded coordinates, then you need to use one of these options.
	// Of course, pro-R or pro-S can only be assigned, if one of the bonds
	// of the pro-chiral center is a stereo bond, which, of course,
	// will be recognized as an over-specified stereo feature.

	// Consider both diastereotopic and enantiotopic atoms as being different when creating the symmetry rank
	public static final int CONSIDER_STEREOHETEROTOPICITY = 2;

	// Consider custom atom labels for atom ranking and encode them into idcodes
	public static final int ENCODE_ATOM_CUSTOM_LABELS = 8;

	// Consider the atom selection for atom ranking and encode it into idcodes
	public static final int ENCODE_ATOM_SELECTION = 16;

	// Assign parities to tetrahedral nitrogen (useful in crystals or at low temp, when N inversion is frozen out)
	public static final int ASSIGN_PARITIES_TO_TETRAHEDRAL_N = 32;

	// Enforces creation of 3D-coordinate encoding even if all z-coords are 0.0
	public static final int COORDS_ARE_3D = 64;

	// The Canonizer normalizes pseudo stereo parities within any rigid fragments, where
	// two representations of relative stereo features encode equal stereo configurations
	// (e.g. cis-1,4-dimethyl-cyclohexane). If the CREATE_PSEUDO_STEREO_GROUPS bit is set,
	// then one can use getPseudoTHGroups() and getPseudoEZGroups(), which return stereo
	// feature group numbers that are shared among all related relative stereo features
	// (tetrahedral and E/Z-bonds).
	public static final int CREATE_PSEUDO_STEREO_GROUPS = 128;

	// If two molecules are identical except for an inverted configuration of stereo centers
	// within one OR group, then they receive the same idcode, because 'this or the other'
	// and 'the other or this' are effectively the same. If we know, however, that we have
	// both enantiomers, but cannot assign which is which, we may want to create different
	// idcodes for either enantiomer. If the mode includes DISTINGUISH_RACEMIC_OR_GROUPS,
	// then the normalization of tetrahedral stereo centers is skipped for OR groups retaining
	// the given configuration within all OR groups.
	public static final int DISTINGUISH_RACEMIC_OR_GROUPS = 256;

	// If we have fragments instead of molecules, then unused valences are not meant to be
	// filled by implicit hydrogen atoms. That means that atoms may have free valences.
	// If two otherwise equivalent (symmetrical) atoms have
	// free valences, then these may differ in the context of a super-structure match.
	// A stereo center in the super-structure may not be a stereo center in the fragment alone.
	// Same is true for stereo bonds. To discover all potential stereo features within a
	// substructure fragment use more CONSIDER_FREE_VALENCES, which breaks the ties
	// between equivalent atoms that have free valences.
	public static final int TIE_BREAK_FREE_VALENCE_ATOMS = 512;

	// ENCODE_ATOM_CUSTOM_LABELS (above) causes customs labels to encode into the idcode in
	// a canonical way, i.e. the label is considered when ranking. Two otherwise symmetrical
	// atoms are considered different, if one has a custom label, or both have different custom
	// labels. The ENCODE_ATOM_CUSTOM_LABELS_WITHOUT_RANKING option does not consider such
	// atoms being different. This option can be used to mark an atom witghout influencing
	// ranking and graph generation. This is typically used for diagnostics.
	public static final int ENCODE_ATOM_CUSTOM_LABELS_WITHOUT_RANKING = 1024;

	public static final int NEGLECT_ANY_STEREO_INFORMATION = 2048;

	protected static final int cIDCodeVersion2 = 8;
		// productive version till May 2006 based on the molfile version 2

	protected static final int cIDCodeVersion3 = 9;
		// productive version since May 2006 based on the molfile version 3
		// being compatible with MDL's "Enhanced Stereo Representation"

	public static final int cIDCodeCurrentVersion = cIDCodeVersion3;
		// New version numbers start with 9 because they are stored in the
		// idcode's first 4 bits which were originally used to store the number
		// of bits for atom indices which did never exceed the value 8.

		// In addition to the four TH/EZ parities from the Molecule
		// idcodes may contain four more types encoding the enhanced
		// stereo representation modes AND and OR.
	protected static final int cParity1And = 4;
	protected static final int cParity2And = 5;
	protected static final int cParity1Or = 6;
	protected static final int cParity2Or = 7;

	private static final int CALC_PARITY_MODE_PARITY = 1;
	private static final int CALC_PARITY_MODE_PRO_PARITY = 2;
	private static final int CALC_PARITY_MODE_IN_SAME_FRAGMENT = 3;

//	public static final int mAtomBits = 16;
//	public static final int MAX_ATOMS = 0xFFFF;
//	public static final int MAX_BONDS = 0xFFFF;

	private StereoMolecule mMol;
	private int[] mCanRank;
	private int[] mCanRankBeforeTieBreaking;
	private int[] mPseudoTHGroup;
	private int[] mPseudoEZGroup;
	private byte[] mTHParity;	        // is tetrahedral parity based on atom symmetry ranks
	private byte[] mEZParity;           // is double bond parity based on atom symmetry ranks
	private byte[] mTHConfiguration;	// is tetrahedral parity based on atom indexes in canonical graph
	private byte[] mEZConfiguration;	// is double bond parity based on atom indexes in canonical graph
	private byte[] mTHCIPParity;
	private byte[] mEZCIPParity;
	private byte[] mTHESRType;
	private byte[] mTHESRGroup;
	private byte[] mEZESRType;
	private byte[] mEZESRGroup;
	private byte[] mAbnormalValence;
	private CanonizerBaseValue[] mCanBase;
	private CanonizerMesoHelper mMesoHelper;
	private boolean mIsMeso,mStereoCentersFound;
	private boolean[] mIsStereoCenter;  // based on extended stereo ranking, i.e. considering ESR type and group
	private boolean[] mTHParityIsMesoInverted;  // whether the atom's parity must be inverted in the idcode because of meso fragment parity normalization
	private boolean[] mTHParityNeedsNormalization;
	private boolean[] mTHESRTypeNeedsNormalization;
	private boolean[] mTHParityRoundIsOdd;
	private boolean[] mEZParityRoundIsOdd;
	private boolean[] mTHParityIsPseudo;
	private boolean[] mEZParityIsPseudo;
	private boolean[] mProTHAtomsInSameFragment;
	private boolean[] mProEZAtomsInSameFragment;
	private boolean[] mNitrogenQualifiesForParity;
	private ArrayList<CanonizerFragment> mFragmentList;
	private ArrayList<int[]> mTHParityNormalizationGroupList;
	private int mMode,mNoOfRanks,mNoOfPseudoGroups;
	private boolean mIsOddParityRound;
	private boolean mZCoordinatesAvailable;
	private boolean mCIPParityNoDistinctionProblem;
	private boolean mEncodeAvoid127;

	private boolean mGraphGenerated;
	private int mGraphRings,mFeatureBlock;
	private int[] mGraphAtom;
	private int[] mGraphIndex;
	private int[] mGraphBond;
	private int[] mGraphFrom;
	private int[] mGraphClosure;

	private String		    mIDCode, mEncodedCoords,mMapping;
	private StringBuilder	mEncodingBuffer;
	private	int				mEncodingBitsAvail,mEncodingTempData,mAtomBits,mMaxConnAtoms;

	/**
	 * Runs a canonicalization procedure for the given molecule that creates unique atom ranks,
	 * which takes stereo features, ESR settings and query features into account.
	 * @param mol
	 */
	public Canonizer(StereoMolecule mol) {
		this(mol, 0);
		}


	/**
	 * Runs a canonicalization procedure for the given molecule that creates unique atom ranks,
	 * which takes stereo features, ESR settings and query features into account.
	 * If mode includes ENCODE_ATOM_CUSTOM_LABELS, than custom atom labels are
	 * considered for the atom ranking and are encoded into the idcode.<br>
	 * If mode includes COORDS_ARE_3D, then getEncodedCoordinates() always returns
	 * a 3D-encoding even if all z-coordinates are 0.0. Otherwise coordinates are
	 * encoded in 3D only, if at least one of the z-coords is not 0.0.
	 * @param mol
	 * @param mode 0 or one or more of CONSIDER...TOPICITY, CREATE..., ENCODE_ATOM_CUSTOM_LABELS, ASSIGN_PARITIES_TO_TETRAHEDRAL_N, COORDS_ARE_3D
	 */
	public Canonizer(StereoMolecule mol, int mode) {
//		if (mol.getAllAtoms()>MAX_ATOMS)
//			throw new IllegalArgumentException("Cannot canonize a molecule having more than "+MAX_ATOMS+" atoms");
//		if (mol.getAllBonds()>MAX_BONDS)
//			throw new IllegalArgumentException("Cannot canonize a molecule having more than "+MAX_BONDS+" bonds");

		mMol = mol;
		mMode = mode;

		mMol.ensureHelperArrays(Molecule.cHelperRings);
		mAtomBits = getNeededBits(mMol.getAtoms());

		if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0)
			canFindNitrogenQualifyingForParity();

		mZCoordinatesAvailable = ((mode & COORDS_ARE_3D) != 0) || mMol.is3D();

		if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0) {
			mTHParity = new byte[mMol.getAtoms()];
			mTHParityIsPseudo = new boolean[mMol.getAtoms()];
			mTHParityRoundIsOdd = new boolean[mMol.getAtoms()];
			mEZParity = new byte[mMol.getBonds()];
			mEZParityRoundIsOdd = new boolean[mMol.getBonds()];
			mEZParityIsPseudo = new boolean[mMol.getBonds()];
			}

		mCIPParityNoDistinctionProblem = false;

		canInitializeRanking();
		if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0)
			canRankStereo();
		canRankFinal();
//		if (mCIPParityNoDistinctionProblem)
//			System.out.println("No distinction applying CIP rules: "+getIDCode()+" "+getEncodedCoordinates());
		}

	public boolean hasCIPParityDistinctionProblem() {
		return mCIPParityNoDistinctionProblem;
		}

	/**
	 * Locate those tetrahedral nitrogen atoms with at least 3 neighbors that
	 * qualify for tetrahedral parity calculation because:<br>
	 * - being a quarternary nitrogen atom<br>
	 * - being an aziridin nitrogen atom<br>
	 * - the configuration inversion is hindered in a polycyclic structure<br>
	 * - flag ASSIGN_PARITIES_TO_TETRAHEDRAL_N is set
	 */
	private void canFindNitrogenQualifyingForParity() {
		mNitrogenQualifiesForParity = new boolean[mMol.getAtoms()];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.getAtomicNo(atom) == 7) {
				if (mMol.getConnAtoms(atom) == 4) {
					mNitrogenQualifiesForParity[atom] = true;
					continue;
					}
				if (mMol.getConnAtoms(atom) == 3) {
					if (mMol.getAtomRingSize(atom) == 3) {
						mNitrogenQualifiesForParity[atom] = true;
						continue;
						}

					if (mMol.getAtomCharge(atom) == 1) {
						mNitrogenQualifiesForParity[atom] = true;
						continue;
						}

					if  (mMol.isFlatNitrogen(atom))
						continue;

					if ((mMode & ASSIGN_PARITIES_TO_TETRAHEDRAL_N) != 0) {
						mNitrogenQualifiesForParity[atom] = true;
						continue;
						}

					if (mMol.getAtomRingBondCount(atom) != 3)
						continue;

					// For any neutral nitrogen with three ring bonds
					// we find the smallest ring. If this is not larger than 7 members
					// then we find that connBond, which does not belong to that ring.
					// The we find the bridge size (atom count) to where it touches the smallest ring.
					// We also find the path length from the touch point on the smallest ring back
					// to the nitrogen atom.
					int smallRingSize = mMol.getAtomRingSize(atom);
					if (smallRingSize > 7)
						continue;

					RingCollection ringSet = mMol.getRingSet();
					int smallRingNo = 0;
					while (smallRingNo<ringSet.getSize()) {
						if (ringSet.getRingSize(smallRingNo) == smallRingSize
						 && ringSet.isAtomMember(smallRingNo, atom))
							break;

						smallRingNo++;
						}

					if (smallRingNo >= RingCollection.MAX_SMALL_RING_COUNT
					 && smallRingNo == ringSet.getSize())
						continue;   // very rare case, but found with wrongly highly bridged CSD entry JORFAZ

					int firstBridgeAtom = -1;
					int firstBridgeBond = -1;
					for (int i=0; i<3; i++) {
						int connBond = mMol.getConnBond(atom, i);
						if (!ringSet.isBondMember(smallRingNo, connBond)) {
							firstBridgeAtom = mMol.getConnAtom(atom, i);
							firstBridgeBond = connBond;
							break;
							}
						}
					
					boolean[] neglectBond = new boolean[mMol.getBonds()];
					neglectBond[firstBridgeBond] = true;
					int[] pathAtom = new int[11];
					int pathLength = mMol.getPath(pathAtom, firstBridgeAtom, atom, 10, neglectBond);
					if (pathLength == -1)
						continue;

					int bridgeAtomCount = 1;
					while (!ringSet.isAtomMember(smallRingNo, pathAtom[bridgeAtomCount]))
						bridgeAtomCount++;

					int bondCountToBridgeHead = pathLength - bridgeAtomCount;

					int bridgeHead = pathAtom[bridgeAtomCount];

					if (smallRingSize == 6 && bondCountToBridgeHead == 2 && bridgeAtomCount == 3) {
						if (mMol.getAtomRingBondCount(pathAtom[1]) >= 3) {
							boolean isAdamantane = false;
							int[] ringAtom = ringSet.getRingAtoms(smallRingNo);
							for (int i=0; i<6; i++) {
								if (atom == ringAtom[i]) {
									int potentialOtherBridgeHeadIndex = ringSet.validateMemberIndex(smallRingNo,
											(bridgeHead == ringAtom[ringSet.validateMemberIndex(smallRingNo, i+2)]) ? i-2 : i+2);
									int potentialOtherBridgeHead = ringAtom[potentialOtherBridgeHeadIndex];
									if (mMol.getAtomRingBondCount(potentialOtherBridgeHead) >= 3
									 && mMol.getPathLength(pathAtom[1], potentialOtherBridgeHead, 2, null) == 2)
										isAdamantane = true;
									break;
									}
								}
							if (isAdamantane) {
								mNitrogenQualifiesForParity[atom] = true;
								continue;
								}
							}
						}

					boolean bridgeHeadIsFlat = (mMol.getAtomPi(bridgeHead) == 1
											 || mMol.isAromaticAtom(bridgeHead)
											 || mMol.isFlatNitrogen(bridgeHead));
					boolean bridgeHeadMayInvert = !bridgeHeadIsFlat
												&& mMol.getAtomicNo(bridgeHead) == 7
												&& mMol.getAtomCharge(bridgeHead) != 1;

					if (bondCountToBridgeHead == 1) {
						if (!bridgeHeadIsFlat && !bridgeHeadMayInvert && smallRingSize <= 4 && bridgeAtomCount <= 3)
							mNitrogenQualifiesForParity[atom] = true;
						continue;
						}

					switch (smallRingSize) {
					// case 3 is fully handled
					case 4:	// must be bondCountToBridgeHead == 2
						if (!bridgeHeadIsFlat && !bridgeHeadMayInvert) {
							if (bridgeAtomCount <= 4)
								mNitrogenQualifiesForParity[atom] = true;
							}
						break;
					case 5:	// must be bondCountToBridgeHead == 2
						if (bridgeHeadMayInvert) {
							if (bridgeAtomCount <= 3)
								mNitrogenQualifiesForParity[atom] = true;
							}
						else if (!bridgeHeadIsFlat) {
							if (bridgeAtomCount <= 4)
								mNitrogenQualifiesForParity[atom] = true;
							}
						break;
					case 6:
						if (bondCountToBridgeHead == 2) {
							if (bridgeHeadIsFlat) {
								if (bridgeAtomCount <= 4)
									mNitrogenQualifiesForParity[atom] = true;
								}
							else if (!bridgeHeadMayInvert) {
								if (bridgeAtomCount <= 3)
									mNitrogenQualifiesForParity[atom] = true;
								}
							}
						else if (bondCountToBridgeHead == 3) {
							if (bridgeHeadIsFlat) {
								if (bridgeAtomCount <= 6)
									mNitrogenQualifiesForParity[atom] = true;
								}
							else {
								if (bridgeAtomCount <= 4)
									mNitrogenQualifiesForParity[atom] = true;
								}
							}
						break;
					case 7:
						if (bondCountToBridgeHead == 3) {
							if (bridgeAtomCount <= 3)
								mNitrogenQualifiesForParity[atom] = true;
							}
						break;
						}
					}
				}
			}
		}

	/**
	 * Calculates and stores implicit abnormal valences due to 
	 * explicit hydrogen attachments.
	 * @param atom
	 * @return
	 */
	private int canCalcImplicitAbnormalValence(int atom) {
		int explicitAbnormalValence = mMol.getAtomAbnormalValence(atom);
		int implicitHigherValence = mMol.getImplicitHigherValence(atom, false);
		int newImplicitHigherValence = mMol.getImplicitHigherValence(atom, true);

		int valence = -1;
		if (implicitHigherValence != newImplicitHigherValence) {
			if (explicitAbnormalValence != -1
			 && explicitAbnormalValence > implicitHigherValence)
				valence = (byte)explicitAbnormalValence;
			else
				valence = (byte)implicitHigherValence;
			}
		else if (explicitAbnormalValence != -1) {
			if (explicitAbnormalValence > newImplicitHigherValence
			 || (explicitAbnormalValence < newImplicitHigherValence
			  && explicitAbnormalValence >= mMol.getOccupiedValence(atom)))
				valence = (byte)explicitAbnormalValence;
			}
		else if (!mMol.supportsImplicitHydrogen(atom)
			  && mMol.getExplicitHydrogens(atom) != 0) {
			valence = mMol.getOccupiedValence(atom);
			valence -= mMol.getElectronValenceCorrection(atom, valence);
			}

		canSetAbnormalValence(atom, valence);
		return valence;
		}

	private void canSetAbnormalValence(int atom, int valence) {
    	if (mAbnormalValence == null) {
			mAbnormalValence = new byte[mMol.getAtoms()];
			Arrays.fill(mAbnormalValence, (byte)-1);
    		}
    	mAbnormalValence[atom] = (byte)valence;
		}

	private void canRankStereo() {
		// Store ranking state before considering stereo information
		int noOfRanksWithoutStereo = mNoOfRanks;
		int[] canRankWithoutStereo = Arrays.copyOf(mCanRank, mMol.getAtoms());

		// Calculate the Cahn-Ingold-Prelog stereo assignments based
		// on drawn stereo bonds neglecting any ESR group assignments
		if (!mMol.isFragment()) {
			canRecursivelyFindCIPParities();
//System.out.println("after CIP parity ranking");

			// rollback stereo information
			initializeParities(noOfRanksWithoutStereo, canRankWithoutStereo);
			}

		mTHESRType = new byte[mMol.getAtoms()];
		mTHESRGroup = new byte[mMol.getAtoms()];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			mTHESRType[atom] = (byte)mMol.getAtomESRType(atom);
			mTHESRGroup[atom] = (byte)mMol.getAtomESRGroup(atom);
			}
		mEZESRType = new byte[mMol.getBonds()];
		mEZESRGroup = new byte[mMol.getBonds()];
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			mEZESRType[bond] = (byte)mMol.getBondESRType(bond);
			mEZESRGroup[bond] = (byte)mMol.getBondESRGroup(bond);
			}

   		// In meso fragments we may have arbitrary ESR group assignments,
		// i.e. there are different possible ESR groupings defining the same
		// meso fragment. Thus, we need to normalize ESR grouping before
		// considering this information at all.
		// Stereo ranking is based on parities, which in turn can only be
		// relied on, if they belong to ABS stereo centers (rather than OR or AND).
		// Thus, before we do the real stereo ranking we need to run a limited
		// ranking first, that enables us to detect meso fragments. Then we
		// roll back, normalize the ESR grouping of the meso fragments and do
		// discover all stereo centers a second time, now based on normalized
		// ESR groups.

		canRecursivelyFindAllParities();
//System.out.println("after meso parity ranking");

		// indicate all stereo centers
		mStereoCentersFound = false;
		mIsStereoCenter = new boolean[mMol.getAtoms()];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mTHParity[atom] != Molecule.cAtomParityNone) {
				mIsStereoCenter[atom] = true;
				mStereoCentersFound = true;
				}
			}
/*
System.out.print("all stereo centers:");
for (int atom=0; atom<mMol.getAtoms(); atom++)
 if (mIsStereoCenter[atom])
  System.out.print(" "+atom+",p:"+mTHParity[atom]+",r:"+mCanRank[atom]);
System.out.println();
*/
		// remove all ESR group assignments from non-stereo-centers
		canRemoveOverspecifiedESRGroups();

		mMesoHelper = null;
		mTHESRTypeNeedsNormalization = new boolean[mMol.getAtoms()];
		if (mStereoCentersFound) {
			mMesoHelper = new CanonizerMesoHelper(mMol,
												  canRankWithoutStereo,
												  mIsStereoCenter,
												  mTHParity,
												  mEZParity,
												  mTHESRType,
												  mTHESRGroup,
												  mEZESRType,
												  mEZESRGroup,
												  mTHParityRoundIsOdd,
												  mEZParityRoundIsOdd,
												  mTHESRTypeNeedsNormalization);
			// For every meso fragment
			// put all ABS stereo centers into new AND/OR group.
			mMesoHelper.normalizeESRGroups();
			}

		// schedule all atoms of any ESR group (AND and OR) to be normalized
		// concerning their parities, in order to be able to consider the
		// parities for ranking and ,thus, for recursive parity determination.
		mTHParityIsMesoInverted = new boolean[mMol.getAtoms()];
		mTHParityNeedsNormalization = new boolean[mMol.getAtoms()];
		mTHParityNormalizationGroupList = new ArrayList<>();

		if (mMesoHelper != null)
			mMesoHelper.normalizeESRGroupSwappingAndRemoval(mCanRank);

		canMarkESRGroupsForParityNormalization();

		// rollback stereo information
		initializeParities(noOfRanksWithoutStereo, canRankWithoutStereo);

		// Find real stereo features based on initial ranking and then
		// in a loop check for stereo features that depend on the configuration
		// of other stereo features already found and rank atoms again

		canRecursivelyFindCanonizedParities();
//System.out.println("after parity ranking");

		if (mMesoHelper != null)
			mIsMeso = mMesoHelper.isMeso();

		// at this point all real stereo features have been found

		determineChirality(canRankWithoutStereo);
		}


	private void canRankFinal() {
		// Atoms, which still share an equal rank, can now be considered symmetrical
		// regarding connectivity and stereo features excluding stereo hetero-topicity.
		if ((mMode & CREATE_SYMMETRY_RANK) != 0
		 && (mMode & CONSIDER_STEREOHETEROTOPICITY) == 0) {
			mCanRankBeforeTieBreaking = Arrays.copyOf(mCanRank, mMol.getAtoms());
			}

		if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0) {
			// locate atom differences due to pro-chiral or pro-E/Z location and
			// detect for every proTH- or proEZ-parity whether pro-atoms are
			// in same fragment as the pro-chiral-center or double-bond, respectively
			mProTHAtomsInSameFragment = new boolean[mMol.getAtoms()];
			mProEZAtomsInSameFragment = new boolean[mMol.getBonds()];

			// If we consider stereo information then we try in a reproducible way to distinguish symmetrical atoms
			// considering enantio- or diastereo-topicity, i.e. we check for pro-R or -S atoms with hitherto equal ranks.
			if (mNoOfRanks < mMol.getAtoms()) {
				canBreakTiesByHeteroTopicity();

				if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0) {
					canNormalizeGroupParities();
					if (mMesoHelper != null)
						mMesoHelper.normalizeESRGroupSwappingAndRemoval(mCanRank);
					}
				}
			}

		// Atoms, which still share an equal rank now, can now be considered symmetrical
		// regarding connectivity and stereo features including stereo hetero-topicity.
		if (mCanRankBeforeTieBreaking == null
		 && (mMode & CREATE_SYMMETRY_RANK) != 0
		 && (mMode & CONSIDER_STEREOHETEROTOPICITY) != 0)
			mCanRankBeforeTieBreaking = Arrays.copyOf(mCanRank, mMol.getAtoms());

		// ############### begin tie breaking ##############
		// i.e. if not all atoms have a different rank yet, then
		// select randomly one atom of those atoms that share the
		// lowest rank and consider it higher ranked. Propagate the
		// new ranking asymmetry through the molecule (performFullRanking()).
		// Repeat these steps until all atoms have a different rank.

//System.out.println("start of tie breaking");
/*
for (int atom=0; atom<mMol.getAtoms(); atom++)
System.out.println("mTHParity["+atom+"] = "+mTHParity[atom]+" roundIsOdd = "+mTHParityRoundIsOdd+" pseudo = "+mTHParityIsPseudo[atom]);
for (int bond=0; bond<mMol.getBonds(); bond++)
System.out.println("mEZParity["+bond+"] = "+mEZParity[bond]);
*/
		// If no further resolution can be done using different hetero-topicity, then we promote one atom of lowest shared rank randomly.
		while (mNoOfRanks < mMol.getAtoms()) {
			canBreakTiesRandomly();

			if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0) {
				canNormalizeGroupParities();
				if (mMesoHelper != null)
					mMesoHelper.normalizeESRGroupSwappingAndRemoval(mCanRank);
				}
			}

			// normalize those group's parities that could not be normalized before tie-breaking
		if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0) {
			canNormalizeGroupParities();

			// detect all not yet discovered pseudo-parities
			canFindPseudoParities();
			flagStereoProblems();
			}
/*
for (int atom=0; atom<mMol.getAtoms(); atom++)
System.out.println("mTHParity["+atom+"] = "+mTHParity[atom]+" roundIsOdd = "+mTHParityRoundIsOdd+" pseudo = "+mTHParityIsPseudo[atom]);
for (int bond=0; bond<mMol.getBonds(); bond++)
System.out.println("mEZParity["+bond+"] = "+mEZParity[bond]);
*/
		}


/*	private boolean canBreakTiesByHeteroTopicity() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			mCanBase[atom].init(atom);
			mCanBase[atom].add((2*mAtomBits+4), (long)mCanRank[atom] << (mAtomBits+4));
			}

		boolean found = false;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			found |= canCalcTHParity(atom, true);
		for (int bond=0; bond<mMol.getBonds(); bond++)
			found |= canCalcEZParity(bond, true);

		if (found)
			mNoOfRanks = canPerformRanking();

		return found;
		}*/

	private boolean canBreakTiesByHeteroTopicity() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			mCanBase[atom].init(atom);
			mCanBase[atom].add((2*mAtomBits+4), (long)mCanRank[atom] << (mAtomBits+4));
			}

		// create in-same-fragment information for all heterotopic atoms
		boolean found = false;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			found |= canCalcTHParity(atom, CALC_PARITY_MODE_IN_SAME_FRAGMENT);
		for (int bond=0; bond<mMol.getBonds(); bond++)
			found |= canCalcEZParity(bond, CALC_PARITY_MODE_IN_SAME_FRAGMENT);

		if (!found)
			return false;

		while (mNoOfRanks < mMol.getAtoms()) {
			found = canInnerBreakTiesByHeteroTopicity();
			if (!found)
				break;

			canNormalizeGroupParities();
			if (mMesoHelper != null)
				mMesoHelper.normalizeESRGroupSwappingAndRemoval(mCanRank);
			}

		return true;
		}


	private boolean canInnerBreakTiesByHeteroTopicity() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			mCanBase[atom].init(atom);
			mCanBase[atom].add((2*mAtomBits+4), (long)mCanRank[atom] << (mAtomBits+4));
			}

		for (int rank=1; rank<=mNoOfRanks; rank++) {
			boolean found = false;
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (mCanRank[atom] == rank)
					found |= canCalcTHParity(atom, CALC_PARITY_MODE_PRO_PARITY);

			if (found) {
				int oldRanks = mNoOfRanks;
				mNoOfRanks = canPerformRanking();
				if (mNoOfRanks != oldRanks)
					return true;

				for (int atom=0; atom<mMol.getAtoms(); atom++) {
					mCanBase[atom].init(atom);
					mCanBase[atom].add((2*mAtomBits+4), (long)mCanRank[atom] << (mAtomBits+4));
					}
				}
			}

		CanonizerBond[] rankedBond = new CanonizerBond[mMol.getBonds()];
		for (int i = 0; i<rankedBond.length; i++)
			rankedBond[i] = new CanonizerBond(mCanRank[mMol.getBondAtom(0, i)], mCanRank[mMol.getBondAtom(1, i)], i);
		Arrays.sort(rankedBond);

		for (int i=0; i<rankedBond.length; i++) {
			if (canCalcEZParity(rankedBond[i].bond, CALC_PARITY_MODE_PRO_PARITY)) {
				while (i+1 < rankedBond.length && rankedBond[i].compareTo(rankedBond[i+1]) == 0)
					canCalcEZParity(rankedBond[++i].bond, CALC_PARITY_MODE_PRO_PARITY);

				int oldRanks = mNoOfRanks;
				mNoOfRanks = canPerformRanking();
				if (mNoOfRanks != oldRanks)
					return true;

				for (int atom=0; atom<mMol.getAtoms(); atom++) {
					mCanBase[atom].init(atom);
					mCanBase[atom].add((2*mAtomBits+4), (long)mCanRank[atom] << (mAtomBits+4));
					}
				}
			}

		return false;
		}

	private void canBreakTiesRandomly() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			mCanBase[atom].init(atom);
			mCanBase[atom].add(mAtomBits+1, 2*mCanRank[atom]);
			}

		// promote randomly one atom of lowest shared rank.
		int[] rankCount = new int[mNoOfRanks+1];
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			rankCount[mCanRank[atom]]++;
		int rank = 1;
		while (rankCount[rank] == 1)
			rank++;

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mCanRank[atom] == rank) {
				mCanBase[atom].add(1);
				break;
				}
			}

		mNoOfRanks = canPerformRanking();
		}

	private void initializeParities(int noOfRanksWithoutStereo, int[] canRankWithoutStereo) {
		mNoOfRanks = noOfRanksWithoutStereo;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			mCanRank[atom] = canRankWithoutStereo[atom];
			mTHParity[atom] = 0;
			mTHParityRoundIsOdd[atom] = false;
			}
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			mEZParity[bond] = 0;
			mEZParityRoundIsOdd[bond] = false;
			}
		}

	/**
	 * This ranks all atoms based on inherent atom properties, their neighbourhood
	 * and connecting bond types until all atom and connectivity differences have
	 * caused different atom ranks. Only stereo chemistry is not considered in this step.
	 */
	private void canInitializeRanking() {
		boolean bondQueryFeaturesPresent = false;
		if (mMol.isFragment()) {
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (mMol.getBondQueryFeatures(bond) != 0) {
					bondQueryFeaturesPresent = true;
					break;
					}
				}
			}

		mMaxConnAtoms = 2;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			mMaxConnAtoms = Math.max(mMaxConnAtoms, mMol.getConnAtoms(atom)+mMol.getMetalBondedConnAtoms(atom));
		int baseValueSize = Math.max(2, bondQueryFeaturesPresent ?
				(62 + mAtomBits + mMaxConnAtoms * (mAtomBits+Molecule.cBondQFNoOfBits)) / 63
			  : (62 + mAtomBits + mMaxConnAtoms * (mAtomBits+5)) / 63);

		mCanRank = new int[mMol.getAllAtoms()];
		mCanBase = new CanonizerBaseValue[mMol.getAtoms()];
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			mCanBase[atom] = new CanonizerBaseValue(baseValueSize);

		boolean atomListFound = false;

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			mCanBase[atom].init(atom);
			if (((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0)
			 || mMol.getAtomList(atom) != null)
				mCanBase[atom].add(8, 6);
			else
				mCanBase[atom].add(8, mMol.getAtomicNo(atom));
			mCanBase[atom].add(8, mMol.getAtomMass(atom));
			mCanBase[atom].add(2, mMol.getAtomPi(atom));
			mCanBase[atom].add(4, mMol.getConnAtoms(atom)+mMol.getMetalBondedConnAtoms(atom));
			if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0)
				mCanBase[atom].add(4, 8);
			else
				mCanBase[atom].add(4, 8 + mMol.getAtomCharge(atom));
			mCanBase[atom].add(5, Math.min(31, mMol.getAtomRingSize(atom)));
			mCanBase[atom].add(4, canCalcImplicitAbnormalValence(atom)+1);
			mCanBase[atom].add(2, mMol.getAtomRadical(atom) >> Molecule.cAtomRadicalStateShift);

			if (mMol.isFragment()) {
				mCanBase[atom].add(Molecule.cAtomQFNoOfBits, mMol.getAtomQueryFeatures(atom));
				if (mMol.getAtomList(atom) != null)
					atomListFound = true;
				}
			}

		mNoOfRanks = canPerformRanking();

		// In very rare cases we need to consider the bond ring size (we neglect rings caused by metal ligand bonds)
		if (mNoOfRanks < mMol.getAtoms()) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(mAtomBits, mCanRank[atom]);
				int[] bondRingSize = new int[mMol.getConnAtoms(atom)];
				for (int i=0; i<mMol.getConnAtoms(atom); i++) {
					bondRingSize[i] = mCanRank[mMol.getConnAtom(atom, i)] << 5;
					bondRingSize[i] |= Math.min(31, mMol.getBondRingSize(mMol.getConnBond(atom, i)));
					}
				Arrays.sort(bondRingSize);
				for (int i=mMaxConnAtoms; i>bondRingSize.length; i--)
					mCanBase[atom].add(mAtomBits+5, 0);
				for (int i=bondRingSize.length-1; i>=0; i--)
					mCanBase[atom].add(mAtomBits+5, bondRingSize[i]);
				}
			mNoOfRanks = canPerformRanking();
			}

		if (atomListFound && mNoOfRanks < mMol.getAtoms()) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(mAtomBits, mCanRank[atom]);
				int[] atomList = mMol.getAtomList(atom);
				int listLength = (atomList == null) ? 0 : Math.min(12, atomList.length);
				for (int i=12; i>listLength; i--)
					mCanBase[atom].add(8, 0);
				for (int i=listLength-1; i>=0; i--)
					mCanBase[atom].add(8, atomList[i]);
				}

			mNoOfRanks = canPerformRanking();
			}

		if (bondQueryFeaturesPresent && mNoOfRanks < mMol.getAtoms()) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(mAtomBits, mCanRank[atom]);
				long[] bondQFList = new long[mMol.getConnAtoms(atom)+mMol.getMetalBondedConnAtoms(atom)];
				int index = 0;
				for (int i=0; i<mMol.getAllConnAtomsPlusMetalBonds(atom); i++) {
					if (i< mMol.getConnAtoms(atom) || i>=mMol.getAllConnAtoms(atom)) {
						bondQFList[index] = mCanRank[mMol.getConnAtom(atom, i)];
						bondQFList[index] <<= Molecule.cBondQFNoOfBits;
						bondQFList[index] |= mMol.getBondQueryFeatures(mMol.getConnBond(atom, i));
						index++;
						}
					}
				Arrays.sort(bondQFList);
				for (int i=mMaxConnAtoms; i>bondQFList.length; i--)
					mCanBase[atom].add(mAtomBits+Molecule.cBondQFNoOfBits, 0);
				for (int i=bondQFList.length-1; i>=0; i--)
					mCanBase[atom].add(mAtomBits+Molecule.cBondQFNoOfBits, bondQFList[i]);
				}
			mNoOfRanks = canPerformRanking();
			}

		if ((mMode & ENCODE_ATOM_CUSTOM_LABELS) != 0 && mNoOfRanks < mMol.getAtoms()) {
			SortedStringList list = new SortedStringList();
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (mMol.getAtomCustomLabel(atom) != null)
					list.addString(mMol.getAtomCustomLabel(atom));

			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				int rank = (mMol.getAtomCustomLabel(atom) == null) ?
						   0 : 1+list.getListIndex(mMol.getAtomCustomLabel(atom));
				mCanBase[atom].init(atom);
				mCanBase[atom].add(mAtomBits, mCanRank[atom]);
				mCanBase[atom].add(mAtomBits, rank);
				}

			mNoOfRanks = canPerformRanking();
			}

		if ((mMode & ENCODE_ATOM_SELECTION) != 0 && mNoOfRanks < mMol.getAtoms()) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(mAtomBits, mCanRank[atom]);
				mCanBase[atom].add(1, mMol.isSelectedAtom(atom) ? 1 : 0);
				}

			mNoOfRanks = canPerformRanking();
			}

		if ((mMode & TIE_BREAK_FREE_VALENCE_ATOMS) != 0 && mMol.isFragment())
			canBreakFreeValenceAtomTies();

		//System.out.println("after initial ranking");
		}


	private void canBreakFreeValenceAtomTies() {
		while (true) {
			boolean[] isFreeValenceRank = new boolean[mNoOfRanks+1];
			int highestSharedFreeValenceRank = -1;
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mMol.getLowestFreeValence(atom) != 0) {
					if (isFreeValenceRank[mCanRank[atom]] && highestSharedFreeValenceRank < mCanRank[atom])
						highestSharedFreeValenceRank = mCanRank[atom];
					isFreeValenceRank[mCanRank[atom]] = true;
					}
				}

			if (highestSharedFreeValenceRank == -1)
				break;

			int increment = 0;
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				int value = 0;
				if (mCanRank[atom] == highestSharedFreeValenceRank)
					value = ++increment;
				mCanBase[atom].init(atom);
				mCanBase[atom].add(mAtomBits, mCanRank[atom]);
				mCanBase[atom].add(8, value);
				}

			mNoOfRanks = canPerformRanking();
			}
		}


	private void canRemoveOverspecifiedESRGroups() {
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (!mIsStereoCenter[atom]
			 || mTHParity[atom] == Molecule.cAtomParityUnknown)
				mTHESRType[atom] = Molecule.cESRTypeAbs;
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mMol.getBondType(bond) != Molecule.cBondTypeSingle
			 || mEZParity[bond] == Molecule.cBondParityNone
			 || mEZParity[bond] == Molecule.cBondParityUnknown)
				mEZESRType[bond] = Molecule.cESRTypeAbs;
		}


	/**
	 * Detects all stereo centers and their parities while considering
	 * absolute ESR group numbers and the drawn configuration of stereo
	 * centers within OR or AND groups. The calculated parities are not
	 * canonized and only intended for meso fragment detection and complete
	 * detection of all stereo centers.
	 */
	private void canRecursivelyFindAllParities() {
		mIsOddParityRound = true;

		// Find absolute stereo features based on current ranking
		boolean paritiesFound = canFindParities(false);

			// in a loop check for stereo features that depend on the configuration
			// of other stereo features already found and rank atoms again
		final int parityInfoBits = 2 + 2 + Molecule.cESRGroupBits;
		while ((mNoOfRanks < mMol.getAtoms()) && paritiesFound) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(mAtomBits, mCanRank[atom]);

					// In case of stereo centers with ESR type!=ABS then consider
					// group number and type in addition to the parity.
				int thParityInfo = mTHParity[atom] << (2 + Molecule.cESRGroupBits);
				if ((mTHParity[atom] == Molecule.cAtomParity1
				  || mTHParity[atom] == Molecule.cAtomParity2)
				 && mTHESRType[atom] != Molecule.cESRTypeAbs) {
					thParityInfo |= (mTHESRType[atom] << Molecule.cESRGroupBits);
					thParityInfo |= mTHESRGroup[atom];
					}

				mCanBase[atom].add(2 * parityInfoBits, thParityInfo << parityInfoBits); // generate space for bond parity
				}

			for (int bond=0; bond<mMol.getBonds(); bond++) {
				int ezParityInfo = mEZParity[bond] << (2 + Molecule.cESRGroupBits);
				if ((mEZParity[bond] == Molecule.cBondParityEor1
				  || mEZParity[bond] == Molecule.cBondParityZor2)
				  && mMol.getBondType(bond) == Molecule.cBondTypeSingle
				  && mEZESRType[bond] != Molecule.cESRTypeAbs) {
					ezParityInfo |= (mEZESRType[bond] << Molecule.cESRGroupBits);
					ezParityInfo |= mEZESRGroup[bond];
					}

				mCanBase[mMol.getBondAtom(0, bond)].add(ezParityInfo);
				mCanBase[mMol.getBondAtom(1, bond)].add(ezParityInfo);
				}
	
			int newNoOfRanks = canPerformRanking();
			if (mNoOfRanks == newNoOfRanks)
				break;

			mNoOfRanks = newNoOfRanks;
			paritiesFound = canFindParities(false);
			}
		}


	/**
	 * Detects stereo center's parities considering the drawn configuration.
	 * Any ESR group assignments are neglected. From these parities the CIP
	 * assignments are calculated. These may be for display purposes only.
	 */
	private void canRecursivelyFindCIPParities() {
		mIsOddParityRound = true;

		// Find absolute stereo features based on current ranking
		mTHCIPParity = new byte[mMol.getAtoms()];
		mEZCIPParity = new byte[mMol.getBonds()];
		boolean paritiesFound = canFindParities(true);

			// in a loop check for stereo features that depend on the configuration
			// of other stereo features already found and rank atoms again
		while ((mNoOfRanks < mMol.getAtoms()) && paritiesFound) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(mAtomBits+4, (mCanRank[atom] << 4)
											  | (mTHParity[atom] << 2));
				}

			for (int bond=0; bond<mMol.getBonds(); bond++) {
				mCanBase[mMol.getBondAtom(0,bond)].add(mEZParity[bond]);
				mCanBase[mMol.getBondAtom(1,bond)].add(mEZParity[bond]);
				}
	
			int newNoOfRanks = canPerformRanking();
			if (mNoOfRanks == newNoOfRanks)
				break;

			mNoOfRanks = newNoOfRanks;
			paritiesFound = canFindParities(true);
			}
		}


	/**
	 * Detects parities for full normalization and idcode generation
	 * (does not consider absolute group numbers nor drawn configuration of
	 *  stereo centers within OR or AND groups, before they are normalized)
	 */
	private void canRecursivelyFindCanonizedParities() {
		mIsOddParityRound = true;

		int[][][] esrGroupMember = compileESRGroupMembers();

		if (mMesoHelper != null
		 && mMesoHelper.normalizeESRGroupSwappingAndRemoval(mCanRank))
			esrGroupMember = compileESRGroupMembers();

			// Find absolute stereo features based on current ranking.
			// If stereo features could be found try to normalize ambiguously
			// definable parities within ESR groups.
		if (canFindParities(false))
			canNormalizeGroupParities();

			// in a loop check for stereo features that depend on the configuration
			// of other stereo features already found and rank atoms again
		boolean newStereoInfoAvailable = true;
		while ((mNoOfRanks < mMol.getAtoms()) && newStereoInfoAvailable) {
			int[][] groupRank = canGetESRGroupRank(esrGroupMember);

			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(mAtomBits, mCanRank[atom]);
				mCanBase[atom].add(20, 0);

				// Certain groups of parities require normalization before
				// the definite parity values can be considered here. These
				// groups of stereo centers consist of all ESR groups of type
				// OR and AND, but also of independent ABS atoms in meso
				// fragments.
				// If a stereo center is marked to require normalization
				// then only consider the ESR type and the rank of the group.

				if (!mTHESRTypeNeedsNormalization[atom]
						&& mTHESRType[atom] != Molecule.cESRTypeAbs)
					mCanBase[atom].add((mTHESRType[atom] << 18)
							+ (groupRank[(mTHESRType[atom] == Molecule.cESRTypeAnd) ? 0 : 1][mTHESRGroup[atom]] << 8));

//				if (!mTHParityNeedsNormalization[atom]) {
					int parity = mTHParity[atom];
					if (mTHParityIsMesoInverted[atom]) {
						if (parity == Molecule.cAtomParity1)
							parity = Molecule.cAtomParity2;
						else if (parity == Molecule.cAtomParity2)
							parity = Molecule.cAtomParity1;
						}
					mCanBase[atom].add(parity << 4);
					}
//				}

// TODO consider groupRank for bonds
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				mCanBase[mMol.getBondAtom(0,bond)].add(mEZParity[bond]);
				mCanBase[mMol.getBondAtom(1,bond)].add(mEZParity[bond]);
				}
/*
for (int atom=0; atom<mMol.getAtoms(); atom++)
System.out.println("mCanBaseValue["+atom+"] = "+Long.toHexString(mCanBase[atom].mValue[0])+" parity:"+mTHParity[atom]+" needsNorm:"+mTHParityNeedsNormalization[atom]+" ESRType:"+mTHESRType[atom]);
*/

			int newNoOfRanks = canPerformRanking();
			if (mNoOfRanks == newNoOfRanks)
				break;

			mNoOfRanks = newNoOfRanks;

			newStereoInfoAvailable = false;

			if (mMesoHelper != null
			 && mMesoHelper.normalizeESRGroupSwappingAndRemoval(mCanRank)) {
				newStereoInfoAvailable = true;
				esrGroupMember = compileESRGroupMembers();
				}

			if (canFindParities(false)) {
				newStereoInfoAvailable = true;
				canNormalizeGroupParities();
				}
			}
		}


	private int[][][] compileESRGroupMembers() {
		int[][][] esrGroupMember = new int[2][Molecule.cESRMaxGroups][];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mIsStereoCenter[atom]) {
				if (mTHESRType[atom] == Molecule.cESRTypeAnd)
					esrGroupMember[0][mTHESRGroup[atom]] = CanonizerMesoHelper.addToIntArray(
							esrGroupMember[0][mTHESRGroup[atom]], atom);
				else if (mTHESRType[atom] == Molecule.cESRTypeOr)
					esrGroupMember[1][mTHESRGroup[atom]] = CanonizerMesoHelper.addToIntArray(
							esrGroupMember[0][mTHESRGroup[atom]], atom);
				}
			}
		return esrGroupMember;
		}

	/**
	 * This determines for relative parities within any ESR group
	 * and for any meso fragment group, whether the
	 * parity should be inverted for normalization in the idcode.
	 * The rule is that the highest ranking atom or bond gets parity2
	 * and the others are adapted accordingly.
	 * @return
	 */
	private boolean canNormalizeGroupParities() {
// TODO take also BINAP parities into account
		boolean groupNormalized = false;
		for (int i=0; i<mTHParityNormalizationGroupList.size(); i++) {
			int[] groupAtom = mTHParityNormalizationGroupList.get(i);
			boolean allParitiesDetermined = true;
			int maxRank = -1;
			boolean invertParities = false;
			for (int j=0; j<groupAtom.length; j++) {
				int atom = groupAtom[j];
				if (mTHParity[atom] == Molecule.cAtomParityNone) {
					allParitiesDetermined = false;
					break;
					}
				if (mTHParity[atom] != Molecule.cAtomParityUnknown) {
					boolean isUniqueRank = true;
					for (int k=0; k<groupAtom.length; k++) {
						if (k != j && mCanRank[atom] == mCanRank[groupAtom[k]]) {
							isUniqueRank = false;
							break;
							}
						}
					if (isUniqueRank
					 && maxRank < mCanRank[atom]) {
						maxRank = mCanRank[atom];
						invertParities = (mTHParity[atom] == Molecule.cAtomParity1);
						}
					}
				}
			if (allParitiesDetermined
			 && maxRank != -1) {
				for (int atom:groupAtom) {
					if ((mTHParity[atom] == Molecule.cAtomParity1
					  || mTHParity[atom] == Molecule.cAtomParity2))
						mTHParityIsMesoInverted[atom] = invertParities;
					mTHParityNeedsNormalization[atom] = false;
					}
				mTHParityNormalizationGroupList.remove(groupAtom);
				groupNormalized = true;
				i--;
				}
//System.out.println(" allParitiesDetermined:"+allParitiesDetermined+" maxRank:"+maxRank+" invertParities:"+invertParities);
			}
		return groupNormalized;
		}


	private void canMarkESRGroupsForParityNormalization() {
		int count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mTHESRType[atom] != Molecule.cESRTypeAbs
			 && (mTHESRType[atom] != Molecule.cESRTypeOr || (mMode & DISTINGUISH_RACEMIC_OR_GROUPS) == 0))
				count++;

		if (count == 0)
			return;

		int[] parity = new int[count];
		count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mTHESRType[atom] != Molecule.cESRTypeAbs
			 && (mTHESRType[atom] != Molecule.cESRTypeOr || (mMode & DISTINGUISH_RACEMIC_OR_GROUPS) == 0)) {
				parity[count] = (mTHESRType[atom] << 29)
							  | (mTHESRGroup[atom] << 24)
							  | (mCanRank[atom] << 12)
							  | atom;
				count++;
				}
			}

		Arrays.sort(parity);
		int groupBase = 0;
		int nextGroupBase = 0;
		int groupID = parity[0] & 0xff000000;
		while (true) {
			nextGroupBase++;
			if (nextGroupBase == parity.length
			 || groupID != (parity[nextGroupBase] & 0xff000000)) {
				int[] atomList = new int[nextGroupBase-groupBase];
				for (int i=groupBase; i<nextGroupBase; i++) {
					int atom = parity[i] & 0x00000fff;
					atomList[i-groupBase] = atom;
					mTHParityNeedsNormalization[atom] = true;
					}
				mTHParityNormalizationGroupList.add(atomList);

				if (nextGroupBase == parity.length)
					break;

				groupID = (parity[nextGroupBase] & 0xff000000);
				groupBase = nextGroupBase;
				}
			}
		}

	
	private int[][] canGetESRGroupRank(int[][][] groupMember) {
		int[][] groupRank = new int[2][Molecule.cESRMaxGroups];
		for (int groupTypeIndex=0; groupTypeIndex<2; groupTypeIndex++) {
			int[][] atomRank = new int[Molecule.cESRMaxGroups][];
			int rankCount = 0;
			for (int group=0; group<Molecule.cESRMaxGroups; group++) {
				if (groupMember[groupTypeIndex][group] != null) {
					int memberCount = groupMember[groupTypeIndex][group].length;
					atomRank[group] = new int[memberCount];
					for (int i=0; i<memberCount; i++)
						atomRank[group][i] = mCanRank[groupMember[groupTypeIndex][group][i]];
	
					Arrays.sort(atomRank[group]);
					rankCount++;
					}
				}
			for (int rank=rankCount; rank>0; rank--) {
				int maxGroup = 0;
				int[] maxAtomRank = null;
				for (int group=0; group<Molecule.cESRMaxGroups; group++) {
					if (atomRank[group] != null) {
						if (maxAtomRank == null
						 || maxAtomRank.length < atomRank[group].length) {
							maxAtomRank = atomRank[group];
							maxGroup = group;
							}
						else if (maxAtomRank.length == atomRank[group].length) {
							for (int i=maxAtomRank.length-1; i>=0; i--) {
								if (maxAtomRank[i] < atomRank[group][i]) {
									maxAtomRank = atomRank[group];
									maxGroup = group;
									break;
									}
								}
							}
						}
					}
				groupRank[groupTypeIndex][maxGroup] = rank;
				atomRank[maxGroup] = null;
				}
			}
		return groupRank;
		}


	private boolean canFindParities(boolean doCIP) {
		boolean ezFound = false;
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (canCalcEZParity(bond, CALC_PARITY_MODE_PARITY)) {
//System.out.println("found EZ parity:"+mEZParity[bond]+" bond:"+bond+" roundIsOdd:"+mIsOddParityRound);
				mEZParityRoundIsOdd[bond] = mIsOddParityRound;
				if (doCIP)
					cipCalcEZParity(bond);
				ezFound = true;
				}
		boolean thFound = false;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (canCalcTHParity(atom, CALC_PARITY_MODE_PARITY)) {
				mTHParityRoundIsOdd[atom] = mIsOddParityRound;
				if (doCIP)
					cipCalcTHParity(atom);
				thFound = true;
				}
		if (thFound)
			mIsOddParityRound = !mIsOddParityRound;

		return ezFound || thFound;
		}


	private void determineChirality(int[] canRankWithoutStereo) {
		int stereoCenters = 0;
		int stereoCentersUnknown = 0;
		int stereoCentersTypeAbs = 0;
		int stereoCentersTypeAbsInMesoFragment = 0;
		int stereoCentersTypeAndGroup0 = 0;
		int stereoCentersTypeOrGroup0 = 0;
		int typeAndGroups = 0;
		boolean typeAndInMesoFragmentFound = false;
		boolean[] andGroupUsed = new boolean[Molecule.cESRMaxGroups];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mTHParity[atom] != Molecule.cAtomParityNone) {
				stereoCenters++;
				if (mTHParity[atom] == Molecule.cAtomParityUnknown) {
					stereoCentersUnknown++;
					}
				else {
					if (mTHESRType[atom] == Molecule.cESRTypeAbs) {
						stereoCentersTypeAbs++;
						if (mMesoHelper != null && mMesoHelper.isInMesoFragment(atom))
							stereoCentersTypeAbsInMesoFragment++;
						}
					else if (mTHESRType[atom] == Molecule.cESRTypeOr) {
						if (mTHESRGroup[atom] == 0)
							stereoCentersTypeOrGroup0++;
						}
					else if (mTHESRType[atom] == Molecule.cESRTypeAnd) {
						int group = mTHESRGroup[atom];
						if (!andGroupUsed[group]) {
							typeAndGroups++;
							andGroupUsed[group] = true;
							}
							
						if (mTHESRGroup[atom] == 0)
							stereoCentersTypeAndGroup0++;
						if (mMesoHelper != null && mMesoHelper.isInMesoFragment(atom))
							typeAndInMesoFragmentFound = true;
						}
					}
				}
			}
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mEZParity[bond] != Molecule.cBondParityNone && mMol.getBondType(bond) == Molecule.cBondTypeSingle) {
				stereoCenters++;
				if (mEZParity[bond] == Molecule.cBondParityUnknown) {
					stereoCentersUnknown++;
					}
				else {
					if (mEZESRType[bond] == Molecule.cESRTypeAbs) {
						stereoCentersTypeAbs++;
						if (mMesoHelper != null
						 && mMesoHelper.isInMesoFragment(mMol.getBondAtom(0, bond))
						 && mMesoHelper.isInMesoFragment(mMol.getBondAtom(1, bond)))
							stereoCentersTypeAbsInMesoFragment++;
						}
					else if (mEZESRType[bond] == Molecule.cESRTypeOr) {
						if (mEZESRGroup[bond] == 0)
							stereoCentersTypeOrGroup0++;
						}
					else if (mEZESRType[bond] == Molecule.cESRTypeAnd) {
						int group = mEZESRGroup[bond];
						if (!andGroupUsed[group]) {
							typeAndGroups++;
							andGroupUsed[group] = true;
							}
							
						if (mEZESRGroup[bond] == 0)
							stereoCentersTypeAndGroup0++;
						if (mMesoHelper != null
						 && mMesoHelper.isInMesoFragment(mMol.getBondAtom(0, bond))
						 && mMesoHelper.isInMesoFragment(mMol.getBondAtom(1, bond)))
							typeAndInMesoFragmentFound = true;
						}
					}
				}
			}

		if (stereoCenters == 0) {
			mMol.setChirality(Molecule.cChiralityNotChiral);
			return;
			}

		if (stereoCentersUnknown != 0) {
			mMol.setChirality(Molecule.cChiralityUnknown);
			return;
			}

		if (mIsMeso) {
			mMol.setChirality(Molecule.cChiralityMeso+(1<<typeAndGroups));
			return;
			}

		if (stereoCentersTypeAndGroup0 + stereoCentersTypeAbsInMesoFragment == stereoCenters
		 && !typeAndInMesoFragmentFound) {
			// If we have found type AND stereo centers in a meso fragment, then
			// the same fragment must contain also ABS atoms; otherwise AND would
			// have been normalized to ABS. ABS and AND in the same meso fragment
			// cannot be racemic.
			mMol.setChirality(Molecule.cChiralityRacemic);
			}
		else if (stereoCentersTypeAbs == stereoCenters) {
			mMol.setChirality(Molecule.cChiralityKnownEnantiomer);
			}
		else if (stereoCentersTypeOrGroup0 == stereoCenters) {
			mMol.setChirality(Molecule.cChiralityUnknownEnantiomer);
			}
		else if (stereoCentersTypeAbs == stereoCenters-1
			  && stereoCentersTypeAndGroup0 == 1) {
			mMol.setChirality(Molecule.cChiralityEpimers);
			}
		else {
			mMol.setChirality(Molecule.cChiralityDiastereomers+(1<<typeAndGroups));
			}
		}
	
	private boolean canFindPseudoParities() {
		boolean[] isFreshPseudoParityAtom = new boolean[mMol.getAtoms()];
		boolean[] isFreshPseudoParityBond = new boolean[mMol.getBonds()];
		int anyPseudoParityCount = 0;
		boolean pseudoParity1Or2Found = false;

		if ((mMode & CREATE_PSEUDO_STEREO_GROUPS) != 0) {
			mPseudoTHGroup = new int[mMol.getAtoms()];
			mPseudoEZGroup = new int[mMol.getBonds()];
			}

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mProTHAtomsInSameFragment[atom]) {
				if (!mTHParityIsPseudo[atom]) {
					if (canCalcTHParity(atom, CALC_PARITY_MODE_PARITY)) {
						mTHParityIsPseudo[atom] = true;
						isFreshPseudoParityAtom[atom] = true;
						anyPseudoParityCount++;
						}
					}
				}
			}
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mProEZAtomsInSameFragment[bond]) {
				if (!mEZParityIsPseudo[bond]) {
					if (canCalcEZParity(bond, CALC_PARITY_MODE_PARITY)) {
						mEZParityIsPseudo[bond] = true;
						isFreshPseudoParityBond[bond] = true;
						anyPseudoParityCount++;
						}
					}
				}
			}

		if (anyPseudoParityCount == 1) {	// delete meaningless pseudo stereo feature single in molecule
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (isFreshPseudoParityAtom[atom]) {
					mTHParity[atom] = Molecule.cAtomParityNone;
					break;
					}
				}
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (isFreshPseudoParityBond[bond]) {
					mEZParity[bond] = Molecule.cBondParityNone;
					break;
					}
				}
			}
		else if (anyPseudoParityCount > 1) {
			canEnsureFragments();
			mNoOfPseudoGroups = 0;
			for (CanonizerFragment f:mFragmentList) {
				int pseudoParitiesInGroup = 0;
				int pseudoParity1Or2InGroup = 0;
				int highRankingTHAtom = 0;
				int highRankingEZBond = 0;
				int highTHAtomRank = -1;
				int highEZBondRank = -1;
				for (int i=0; i<f.atom.length; i++) {
					if (isFreshPseudoParityAtom[f.atom[i]]) {
						pseudoParitiesInGroup++;

						if (mTHParity[f.atom[i]] == Molecule.cAtomParity1
						 || mTHParity[f.atom[i]] == Molecule.cAtomParity2) {

							pseudoParity1Or2InGroup++;
							pseudoParity1Or2Found = true;

							if (highTHAtomRank < mCanRank[f.atom[i]]) {
								highTHAtomRank = mCanRank[f.atom[i]];
								highRankingTHAtom = f.atom[i];
								}
							}
						}
					}

				for (int i=0; i<f.bond.length; i++) {
					if (isFreshPseudoParityBond[f.bond[i]]) {
						pseudoParitiesInGroup++;

						int rank1 = mCanRank[mMol.getBondAtom(0,f.bond[i])];
						int rank2 = mCanRank[mMol.getBondAtom(1,f.bond[i])];
						int higherRank = (rank1 > rank2) ? (rank1 << 16) + rank2 : (rank2 << 16) + rank1;
						if (mEZParity[f.bond[i]] == Molecule.cBondParityEor1
						 || mEZParity[f.bond[i]] == Molecule.cBondParityZor2) {

							pseudoParity1Or2InGroup++;
							pseudoParity1Or2Found = true;

							if (highEZBondRank < higherRank) {
								highEZBondRank = higherRank;
								highRankingEZBond = f.bond[i];
								}
							}
						}
					}

				if (pseudoParitiesInGroup == 0)
					continue;

				if (pseudoParitiesInGroup == 1) {	// delete meaningless pseudo stereo feature single in fragment
					for (int i=0; i<f.atom.length; i++)
						if (isFreshPseudoParityAtom[f.atom[i]])
							mTHParity[f.atom[i]] = Molecule.cAtomParityNone;
					for (int i=0; i<f.bond.length; i++)
						if (isFreshPseudoParityBond[f.bond[i]])
							mEZParity[f.bond[i]] = Molecule.cBondParityNone;
					}
				else {
					if (pseudoParity1Or2InGroup == 1) {
						for (int i=0; i<f.atom.length; i++)
							if (isFreshPseudoParityAtom[f.atom[i]])
								mTHParity[f.atom[i]] = Molecule.cAtomParityUnknown;
						for (int i=0; i<f.bond.length; i++)
							if (isFreshPseudoParityBond[f.bond[i]])
								mEZParity[f.bond[i]] = Molecule.cBondParityUnknown;
						}
					else {	// canonize relative stereo features of fragment
						if ((mMode & CREATE_PSEUDO_STEREO_GROUPS) != 0) {
							mNoOfPseudoGroups++;
							for (int i=0; i<f.atom.length; i++)
								if (isFreshPseudoParityAtom[f.atom[i]])
									mPseudoTHGroup[f.atom[i]] = mNoOfPseudoGroups;
							for (int i=0; i<f.bond.length; i++)
								if (isFreshPseudoParityBond[f.bond[i]])
									mPseudoEZGroup[f.bond[i]] = mNoOfPseudoGroups;
							}

						boolean invertFragmentsStereoFeatures = false;
						if (highTHAtomRank != -1) {
							if (mTHParity[highRankingTHAtom] == Molecule.cAtomParity2)
								invertFragmentsStereoFeatures = true;
							}
						else {
							if (mEZParity[highRankingEZBond] == Molecule.cBondParityZor2)
								invertFragmentsStereoFeatures = true;
							}

						if (invertFragmentsStereoFeatures) {
							for (int i=0; i<f.atom.length; i++) {
								if (isFreshPseudoParityAtom[f.atom[i]]) {
									switch (mTHParity[f.atom[i]]) {
									case Molecule.cAtomParity1:
										mTHParity[f.atom[i]] = Molecule.cAtomParity2;
										break;
									case Molecule.cAtomParity2:
										mTHParity[f.atom[i]] = Molecule.cAtomParity1;
										break;
										}
									}
								}
							for (int i=0; i<f.bond.length; i++) {
								if (isFreshPseudoParityBond[f.bond[i]]) {
									switch (mEZParity[f.bond[i]]) {
									case Molecule.cBondParityEor1:
										mEZParity[f.bond[i]] = Molecule.cBondParityZor2;
										break;
									case Molecule.cBondParityZor2:
										mEZParity[f.bond[i]] = Molecule.cBondParityEor1;
										break;
										}
									}
								}
							}
						}
					}
				}
			}
		return pseudoParity1Or2Found;
		}


	private void canEnsureFragments() {
		if (mFragmentList != null)
			return;

		mFragmentList = new ArrayList<>();

		int fragmentCount = 0;
		int[] fragmentNo = new int[mMol.getAtoms()];
		int[] fragmentAtom = new int[mMol.getAtoms()];
		int[] fragmentBond = new int[mMol.getBonds()];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (fragmentNo[atom] == 0
			 && (mMol.isRingAtom(atom) || mMol.getAtomPi(atom) == 1)) {
				fragmentAtom[0] = atom;
				int fragmentAtoms = 1;
				int fragmentBonds = 0;
				fragmentNo[atom] = ++fragmentCount;
				boolean[] bondHandled = new boolean[mMol.getBonds()];
				for (int current=0; current<fragmentAtoms; current++) {
					for (int i=0; i<mMol.getConnAtoms(fragmentAtom[current]); i++) {
						int connBond = mMol.getConnBond(fragmentAtom[current],i);
						if (mMol.isRingBond(connBond) || mMol.getBondOrder(connBond) == 2 || mMol.isBINAPChiralityBond(connBond)) {
							int connAtom = mMol.getConnAtom(fragmentAtom[current],i);
	
							if (!bondHandled[connBond]) {
								fragmentBond[fragmentBonds++] = connBond;
								bondHandled[connBond] = true;
								}
							if (fragmentNo[connAtom] == 0) {
								fragmentAtom[fragmentAtoms++] = connAtom;
								fragmentNo[connAtom] = fragmentCount;
								}
							}
						}
					}
				mFragmentList.add(new CanonizerFragment(fragmentAtom, fragmentAtoms,
														fragmentBond, fragmentBonds));
				}
			}
		}


	private int canPerformRanking() {
		// Does ranking based on sorted rank list of neighbour atoms and an attached
		// flag if neighbour is connected via a non-aromatic double bond.
		// This needs to be called at least once before tie-breaking in order to
		// distinguish atoms connected in a symmetrical way to to anti-aromatic rings
		// with different distribution of pi-bonds.
		int oldNoOfRanks,newNoOfRanks;

		newNoOfRanks = canConsolidate();
		do {
			oldNoOfRanks = newNoOfRanks;
			canCalcNextBaseValues();
			newNoOfRanks = canConsolidate();
			} while (oldNoOfRanks != newNoOfRanks);

		return newNoOfRanks;
		}


	private void canCalcNextBaseValues() {
		int[] connRank = new int[mMaxConnAtoms];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
								// generate sorted list of ranks of neighbours
			int neighbours = mMol.getConnAtoms(atom)+mMol.getMetalBondedConnAtoms(atom);
			int neighbour = 0;
			for (int i=0; i<mMol.getAllConnAtomsPlusMetalBonds(atom); i++) {
				if (i<mMol.getConnAtoms(atom) || i>=mMol.getAllConnAtoms(atom)) {
					int rank = 2 * mCanRank[mMol.getConnAtom(atom, i)];
					int connBond = mMol.getConnBond(atom, i);
					if (mMol.getBondOrder(connBond) == 2)
						if (!mMol.isAromaticBond(connBond))
							rank++;        // set a flag for non-aromatic double bond
					int j;
					for (j = 0; j < neighbour; j++)
						if (rank < connRank[j])
							break;
					for (int k = neighbour; k > j; k--)
						connRank[k] = connRank[k - 1];
					connRank[j] = rank;
					neighbour++;
					}
				}

			mCanBase[atom].init(atom);
			mCanBase[atom].add(mAtomBits, mCanRank[atom]);
			for (int i=neighbours; i<mMaxConnAtoms; i++)
				mCanBase[atom].add(mAtomBits + 1, 0);
			for (int i=0; i<neighbours; i++)
				mCanBase[atom].add(mAtomBits + 1, connRank[i]);
			}
		}


	private int canConsolidate() {
		int canRank = 0;
		Arrays.sort(mCanBase);
		for (int i=0; i<mCanBase.length; i++) {
			if (i == 0 || mCanBase[i].compareTo(mCanBase[i-1]) != 0)
				canRank++;
			mCanRank[mCanBase[i].getAtom()] = canRank;
			}
/*
System.out.print("Ranks:");
for (int i=0; i<mMol.getAtoms(); i++)
System.out.print(mCanRank[i]+",");
System.out.println("noOfRanks:"+canRank);
*/
		return canRank;
		}


	/**
	 * Calculate tetrahedral parity of a stereo center based on neighbours ranks and
	 * 3D-orientation. Axial (allene type) are also calculated.
	 * If calcProParity is true then this function calculates proparities of those atoms
	 * adjacent to the stereo center rather than a parity of the stereo center itself.
	 * @param atom
	 * @param mode whether to calculate parities, pro-parities, or determine whether heterotopic atoms are in same fragment
	 * @return false if atom's mTHParity != 0; otherwise true if it finds and sets a parity other than 0
	 */
	private boolean canCalcTHParity(int atom, int mode) {
		if (mTHParity[atom] != 0)
			return false;

		if (mMol.getAtomicNo(atom) != 5
		 && mMol.getAtomicNo(atom) != 6
		 && mMol.getAtomicNo(atom) != 7
		 && mMol.getAtomicNo(atom) != 14
		 && mMol.getAtomicNo(atom) != 15
		 && mMol.getAtomicNo(atom) != 16)
			return false;

		if (mMol.getAtomPi(atom) != 0) {
			if (mMol.isCentralAlleneAtom(atom))
				return canCalcAlleneParity(atom, mode);

			if (mMol.getAtomicNo(atom) != 15
			 && mMol.getAtomicNo(atom) != 16)
				return false;
			}

		if (mMol.getConnAtoms(atom) < 3 || mMol.getAllConnAtoms(atom) > 4)
			return false;

		// no carbenium
		if (mMol.getAtomCharge(atom) > 0 && mMol.getAtomicNo(atom) == 6)
			return false;

		// no trivalent boron
		if (mMol.getAtomicNo(atom) == 5 && mMol.getAllConnAtoms(atom) != 4)
			return false;

		// don't consider tetrahedral nitrogen, unless found to qualify for parity calculation
		if (mMol.getAtomicNo(atom) == 7
		 && !mNitrogenQualifiesForParity[atom])
			return false;

				// create array to remap connAtoms according to canRank order
		int[] remappedConn = new int[4];
		int[] remappedRank = new int[4];
		boolean[] neighbourUsed = new boolean[4];
		for (int i=0; i<mMol.getAllConnAtoms(atom); i++) {
			int highestRank = -1;
			int highestConn = 0;
			for (int j=0; j<mMol.getAllConnAtoms(atom); j++) {
				if (!neighbourUsed[j]) {
					if (highestRank < mCanRank[mMol.getConnAtom(atom,j)]) {
						highestRank = mCanRank[mMol.getConnAtom(atom,j)];
						highestConn = j;
						}
					}
				}
			remappedConn[i] = highestConn;
			remappedRank[i] = highestRank;
			neighbourUsed[highestConn] = true;
			}

		if (mMol.getAllConnAtoms(atom) == 4
		 && remappedRank[0] == remappedRank[1]
		 && remappedRank[2] == remappedRank[3])
			return false;			// two pairs of equal ranking atoms

		if ((mMol.getAllConnAtoms(atom) == 4)
		 && (remappedRank[0] == remappedRank[2]
		  || remappedRank[1] == remappedRank[3]))
			return false;			//three equal ranking atoms out of 4

		if (mMol.getAllConnAtoms(atom) == 3
		 && remappedRank[0] == remappedRank[2])
			return false;			//three equal ranking atoms out of 3

		int proTHAtom1 = 0;
		int proTHAtom2 = 0;
		boolean proTHAtomsFound = false;
		for (int i=1; i<mMol.getAllConnAtoms(atom); i++) {
			if (remappedRank[i-1] == remappedRank[i]) {
				if (mode == CALC_PARITY_MODE_PARITY || remappedRank[i] == 0)	// don't calculate pro-parities or hydrogens
					return false;	// two equal ranking atoms -> no stereo center

				proTHAtom1 = mMol.getConnAtom(atom, remappedConn[i-1]);
				proTHAtom2 = mMol.getConnAtom(atom, remappedConn[i]);
				if (mode == CALC_PARITY_MODE_IN_SAME_FRAGMENT
				 && mMol.isRingBond(mMol.getConnBond(atom,remappedConn[i])))
					mProTHAtomsInSameFragment[atom] = true;
				proTHAtomsFound = true;
				}
			}

		if (mode != CALC_PARITY_MODE_PARITY && !proTHAtomsFound)
			return false;

		byte atomTHParity = (mZCoordinatesAvailable) ?
							canCalcTHParity3D(atom, remappedConn)
						  : canCalcTHParity2D(atom, remappedConn);

		if (mode == CALC_PARITY_MODE_PARITY) {
			mTHParity[atom] = atomTHParity;
			}
		else if (mode == CALC_PARITY_MODE_PRO_PARITY) {
			if (atomTHParity == Molecule.cAtomParity1) {
				mCanBase[proTHAtom1].add(mCanRank[atom]);
				}
			else if (atomTHParity == Molecule.cAtomParity2) {
				mCanBase[proTHAtom2].add(mCanRank[atom]);
				}
			}

		return true;
		}


	private byte canCalcTHParity2D(int atom, int[] remappedConn) {
		final int[][] up_down = { { 2,1,2,1 },	// direction of stereobond
								  { 1,2,2,1 },	// for parity = 1
								  { 1,1,2,2 },	// first dimension: order of
								  { 2,1,1,2 },	// angles to connected atoms
								  { 2,2,1,1 },	// second dimension: number of
								  { 1,2,1,2 } };// mMol.getConnAtom that has stereobond

		double[] angle = new double[mMol.getAllConnAtoms(atom)];
		for (int i=0; i<mMol.getAllConnAtoms(atom); i++)
			angle[i] = mMol.getBondAngle(mMol.getConnAtom(atom, remappedConn[i]),atom);

		byte parity = (byte)mMol.getFisherProjectionParity(atom, remappedConn, angle, null);
		if (parity != Molecule.cAtomParityUnknown)
			return parity;

		int stereoBond = 0;	// number of last via stereobond connected Atom
		int stereoType = 0;	// bond geometry (up = 2 / down = 1)
		for (int i=0; i<mMol.getAllConnAtoms(atom); i++) {
			int bnd = mMol.getConnBond(atom,remappedConn[i]);
			if (mMol.getBondAtom(0,bnd) == atom) {
				if (mMol.getBondType(bnd) == Molecule.cBondTypeDown) {
					if (stereoType != 0)
						mMol.setStereoProblem(atom);
					stereoBond = i;
					stereoType = 1;
					}
				if (mMol.getBondType(bnd) == Molecule.cBondTypeUp) {
					if (stereoType != 0)
						mMol.setStereoProblem(atom);
					stereoBond = i;
					stereoType = 2;
					}
				}
			}
		if (stereoType == 0)		// no stereobond at center
			return Molecule.cAtomParityUnknown;

		for (int i=1; i<mMol.getAllConnAtoms(atom); i++)
			if (angle[i] < angle[0])
				angle[i] += Math.PI*2;

		if (mMol.getAllConnAtoms(atom) == 3) {
				// new handling!!! TLS 17.Oct.2005
				// Originally the parity depended solely on clockwise/anti-clockwise orientation
				// of three (priority sorted) angles, which included the angle of the stereo bond.
				// Now to handle strained rings with 'invalid' projections in an expected way
				// (all three bonds are within 180 degrees and the middle bond is the stereo bond
				// then treat this as if the stereo bond would be drawn 180 degrees rotated)
				// the angle of the stereobond is not considered anymore and the parity depends
				// solely on the question if the angle difference between higher priority bond
				// to lower priority bond is less or larger than 180 degrees.
			switch (stereoBond) {
			case 0:
				if (((angle[1] < angle[2]) && (angle[2] - angle[1] < Math.PI))
				 || ((angle[1] > angle[2]) && (angle[1] - angle[2] > Math.PI)))
					stereoType = 3 - stereoType;
				break;
			case 1:
				if (angle[2] - angle[0] > Math.PI)
					stereoType = 3 - stereoType;
				break;
			case 2:
				if (angle[1] - angle[0] < Math.PI)
					stereoType = 3 - stereoType;
				break;
				}

			return (stereoType == 1) ? (byte)Molecule.cAtomParity2 : Molecule.cAtomParity1;
			}

		int order = 0;
		if		(angle[1] <= angle[2] && angle[2] <= angle[3]) order = 0;
		else if (angle[1] <= angle[3] && angle[3] <= angle[2]) order = 1;
		else if (angle[2] <= angle[1] && angle[1] <= angle[3]) order = 2;
		else if (angle[2] <= angle[3] && angle[3] <= angle[1]) order = 3;
		else if (angle[3] <= angle[1] && angle[1] <= angle[2]) order = 4;
		else if (angle[3] <= angle[2] && angle[2] <= angle[1]) order = 5;
		return (up_down[order][stereoBond] == stereoType) ?
				(byte)Molecule.cAtomParity2 : Molecule.cAtomParity1;
		}


	private byte canCalcTHParity3D(int atom, int[] remappedConn) {
		int[] atomList = new int[4];
		for (int i=0; i<mMol.getAllConnAtoms(atom); i++)
			atomList[i] = mMol.getConnAtom(atom, remappedConn[i]);
		if (mMol.getAllConnAtoms(atom) == 3)
			atomList[3] = atom;

		double[][] coords = new double[3][3];
		for (int i=0; i<3; i++) {
			coords[i][0] = mMol.getAtomX(atomList[i+1]) - mMol.getAtomX(atomList[0]);
			coords[i][1] = mMol.getAtomY(atomList[i+1]) - mMol.getAtomY(atomList[0]);
			coords[i][2] = mMol.getAtomZ(atomList[i+1]) - mMol.getAtomZ(atomList[0]);
			}

		// calculate the normal vector (vector product of coords[0] and coords[1])
		double[] n = new double[3];
		n[0] = coords[0][1]*coords[1][2]-coords[0][2]*coords[1][1];
		n[1] = coords[0][2]*coords[1][0]-coords[0][0]*coords[1][2];
		n[2] = coords[0][0]*coords[1][1]-coords[0][1]*coords[1][0];

		// calculate cos(angle) of coords[2] to normal vector
		double cosa = (coords[2][0]*n[0]+coords[2][1]*n[1]+coords[2][2]*n[2])
					/ (Math.sqrt(coords[2][0]*coords[2][0]+coords[2][1]*coords[2][1]+coords[2][2]*coords[2][2])
					 * Math.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]));

		return (cosa > 0.0) ? (byte)Molecule.cAtomParity1 : Molecule.cAtomParity2;
		}


	private boolean canCalcAlleneParity(int atom, int mode) {
		if (mMol.getAtomicNo(atom) != 6
		 && mMol.getAtomicNo(atom) != 7)
				   return false;

		int atom1 = mMol.getConnAtom(atom,0);
		int atom2 = mMol.getConnAtom(atom,1);

		if (mMol.getAtomPi(atom1) != 1 || mMol.getAtomPi(atom2) != 1)
			return false;

		if (mMol.getConnAtoms(atom1) == 1 || mMol.getConnAtoms(atom2) == 1)
			return false;

		if ((mMol.getAllConnAtoms(atom1) > 3) || (mMol.getAllConnAtoms(atom2) > 3))
			return false;

		EZHalfParity halfParity1 = new EZHalfParity(mMol, mCanRank, atom, atom1);
		if (halfParity1.mRanksEqual && mode == CALC_PARITY_MODE_PARITY)
			return false;
		EZHalfParity halfParity2 = new EZHalfParity(mMol, mCanRank,atom, atom2);
		if (halfParity2.mRanksEqual && mode == CALC_PARITY_MODE_PARITY)
			return false;

		if (halfParity1.mRanksEqual && halfParity2.mRanksEqual)
			return false;	// both ends of DB bear equal substituents

		if (mode == CALC_PARITY_MODE_IN_SAME_FRAGMENT) {
			if (halfParity1.mRanksEqual && halfParity1.mInSameFragment)
				mProTHAtomsInSameFragment[atom] = true;
			if (halfParity2.mRanksEqual && halfParity2.mInSameFragment)
				mProTHAtomsInSameFragment[atom] = true;
			}

		byte alleneParity = mZCoordinatesAvailable ? canCalcAlleneParity3D(halfParity1, halfParity2)
												  : canCalcAlleneParity2D(halfParity1, halfParity2);

		if (mode == CALC_PARITY_MODE_PARITY) {	// increment mProParity[] for atoms that are Pro-Parity1
			mTHParity[atom] = alleneParity;
			}
		else if (mode == CALC_PARITY_MODE_PRO_PARITY) {
			if (halfParity1.mRanksEqual) {
				if (alleneParity == Molecule.cAtomParity1) {
					mCanBase[halfParity1.mHighConn].add(mCanRank[atom1]);
					}
				else {
					mCanBase[halfParity1.mLowConn].add(mCanRank[atom1]);
					}
				}
			if (halfParity2.mRanksEqual) {
				if (alleneParity == Molecule.cAtomParity2) {
					mCanBase[halfParity2.mHighConn].add(mCanRank[atom2]);
					}
				else {
					mCanBase[halfParity2.mLowConn].add(mCanRank[atom2]);
					}
				}
			}

		return true;
		}


	private byte canCalcAlleneParity2D(EZHalfParity halfParity1, EZHalfParity halfParity2) {
		int hp1 = halfParity1.getValue();
		int hp2 = halfParity2.getValue();
		if (hp1 == -1 || hp2 == -1 || ((hp1 + hp2) & 1) == 0)
			return Molecule.cAtomParityUnknown;

		byte alleneParity = 0;
		switch (hp1 + hp2) {
			case 3:
			case 7:
				alleneParity = Molecule.cAtomParity2;
				break;
			case 5:
				alleneParity = Molecule.cAtomParity1;
				break;
			}
		return alleneParity;
		}


	private byte canCalcAlleneParity3D(EZHalfParity halfParity1, EZHalfParity halfParity2) {
		int[] atom = new int[4];
		atom[0] = halfParity1.mHighConn;
		atom[1] = halfParity1.mCentralAxialAtom;
		atom[2] = halfParity2.mCentralAxialAtom;
		atom[3] = halfParity2.mHighConn;
		double torsion = mMol.calculateTorsion(atom);
		// if the torsion is not significant (less than ~10 degrees) then return cAtomParityUnknown
		if (Math.abs(torsion) < 0.3 || Math.abs(torsion) > Math.PI-0.3)
			return Molecule.cAtomParityUnknown;
		if (torsion < 0)
			return Molecule.cAtomParity2;
		else
			return Molecule.cAtomParity1;
		}


	private boolean canCalcBINAPParity(int bond, int mode) {
		if (!mMol.isBINAPChiralityBond(bond))
			return false;

		int atom1 = mMol.getBondAtom(0, bond);
		int atom2 = mMol.getBondAtom(1, bond);

		EZHalfParity halfParity1 = new EZHalfParity(mMol, mCanRank, atom1, atom2);
		if (halfParity1.mRanksEqual && mode == CALC_PARITY_MODE_PARITY)
			return false;
		EZHalfParity halfParity2 = new EZHalfParity(mMol, mCanRank, atom2, atom1);
		if (halfParity2.mRanksEqual && mode == CALC_PARITY_MODE_PARITY)
			return false;

		if (halfParity1.mRanksEqual && halfParity2.mRanksEqual)
			return false;	// both ends of DB bear equal substituents

		if (mode == CALC_PARITY_MODE_IN_SAME_FRAGMENT) {
			if (halfParity1.mRanksEqual)	// this is a hack, we should find a better solution considering other proparities as well
				mProEZAtomsInSameFragment[bond] = hasSecondBINAPBond(atom2);
			if (halfParity2.mRanksEqual)
				mProEZAtomsInSameFragment[bond] = hasSecondBINAPBond(atom1);
			}

		byte axialParity = mZCoordinatesAvailable ? canCalcBINAPParity3D(halfParity1, halfParity2)
												  : canCalcBINAPParity2D(halfParity1, halfParity2);

		if (mode == CALC_PARITY_MODE_PARITY) {	// increment mProParity[] for atoms that are Pro-E
			mEZParity[bond] = axialParity;
			}
		else if (mode == CALC_PARITY_MODE_PRO_PARITY) {
			if (halfParity1.mRanksEqual) {
				if (axialParity == Molecule.cBondParityZor2) {
					mCanBase[halfParity1.mHighConn].add(mCanRank[atom2]);
					}
				else {
					mCanBase[halfParity1.mLowConn].add(mCanRank[atom2]);
					}
				}
			if (halfParity2.mRanksEqual) {
				if (axialParity == Molecule.cBondParityZor2) {
					mCanBase[halfParity2.mHighConn].add(mCanRank[atom1]);
					}
				else {
					mCanBase[halfParity2.mLowConn].add(mCanRank[atom1]);
					}
				}
			}

		return true;
		}


	private byte canCalcBINAPParity2D(EZHalfParity halfParity1, EZHalfParity halfParity2) {
		int hp1 = halfParity1.getValue();
		int hp2 = halfParity2.getValue();
		if (hp1 == -1 || hp2 == -1 || ((hp1 + hp2) & 1) == 0)
			return Molecule.cBondParityUnknown;

		byte axialParity = 0;
		switch (hp1 + hp2) {
		case 3:
		case 7:
			axialParity = Molecule.cBondParityEor1;
			break;
		case 5:
			axialParity = Molecule.cBondParityZor2;
			break;
			}
		return axialParity;
		}


	private byte canCalcBINAPParity3D(EZHalfParity halfParity1, EZHalfParity halfParity2) {
		int[] atom = new int[4];
		atom[0] = halfParity1.mHighConn;
		atom[1] = halfParity1.mCentralAxialAtom;
		atom[2] = halfParity2.mCentralAxialAtom;
		atom[3] = halfParity2.mHighConn;
		double torsion = mMol.calculateTorsion(atom);
		// if the torsion is not significant (less than ~10 degrees) then return cBondParityUnknown
		if (Math.abs(torsion) < 0.3 || Math.abs(torsion) > Math.PI-0.3)
			return Molecule.cBondParityUnknown;
		if (torsion < 0)
			return Molecule.cBondParityEor1;
		else
			return Molecule.cBondParityZor2;
		}


	private boolean hasSecondBINAPBond(int atom) {
		RingCollection ringSet = mMol.getRingSet();
		for (int i=0; i<ringSet.getSize(); i++) {
			if (ringSet.isAromatic(i) && ringSet.isAtomMember(i, atom)) {
				for (int j:ringSet.getRingAtoms(i))
					if (j != atom)
						for (int k=0; k<mMol.getConnAtoms(j); k++)
							if (mMol.isBINAPChiralityBond(mMol.getConnBond(j, k)))
								return true;
				return false;
				}
			}
		return false;
		}


	/**
	 * Calculate E/Z parity of a double bond based on their neighbours ranks.
	 * If calcProParity is true then the pro-parities of the allylic atoms are calculated
	 * rather than EZ-parities of double bonds.
	 * @param bond
	 * @param mode whether to calculate parities, pro-parities, or determine whether heterotopic atoms are in same fragment
	 * @return false if bond's mEZParity != 0; otherwise true if it finds and sets a parity other than 0
	 */
	private boolean canCalcEZParity(int bond, int mode) {
			// depending on the flag calcProParity this routine either
			//	- updates EZ-parities of double bonds
			// or - calculates EZ-pro-parities of allylic atoms
		if (mEZParity[bond] != 0)
			return false;

		if (mMol.getBondOrder(bond) == 1)
			return canCalcBINAPParity(bond, mode);

		if (mMol.getBondOrder(bond) != 2)
			return false;

//		if (mMol.isSmallRingBond(bond))		// don't generate EZParity within small rings
//			return false;

		if (mMol.isAromaticBond(bond))	// different mesomeric structures would have EZParities
			return false;				// assigned to different bonds and therefore break the symmetry

		int dbAtom1 = mMol.getBondAtom(0,bond);
		int dbAtom2 = mMol.getBondAtom(1,bond);

		if (mMol.getConnAtoms(dbAtom1) == 1 || mMol.getConnAtoms(dbAtom2) == 1)
			return false;

		if ((mMol.getConnAtoms(dbAtom1) > 3) || (mMol.getConnAtoms(dbAtom2) > 3))
			return false;

		if (mMol.getAtomPi(dbAtom1) == 2 || mMol.getAtomPi(dbAtom2) == 2)	// allene
			return false;

		EZHalfParity halfParity1 = new EZHalfParity(mMol, mCanRank, dbAtom2, dbAtom1);
		if (halfParity1.mRanksEqual && mode == CALC_PARITY_MODE_PARITY)
			return false;
		EZHalfParity halfParity2 = new EZHalfParity(mMol, mCanRank, dbAtom1, dbAtom2);
		if (halfParity2.mRanksEqual && mode == CALC_PARITY_MODE_PARITY)
			return false;

		if (halfParity1.mRanksEqual && halfParity2.mRanksEqual)
			return false;	// both ends of DB bear equal substituents

		if (mode == CALC_PARITY_MODE_IN_SAME_FRAGMENT) {
			if (halfParity1.mRanksEqual && halfParity1.mInSameFragment)
				mProEZAtomsInSameFragment[bond] = true;
			if (halfParity2.mRanksEqual && halfParity2.mInSameFragment)
				mProEZAtomsInSameFragment[bond] = true;
			}

		byte bondDBParity = mMol.isBondParityUnknownOrNone(bond) ? Molecule.cBondParityUnknown
							  		  : (mZCoordinatesAvailable) ? canCalcEZParity3D(halfParity1, halfParity2)
							  				  					 : canCalcEZParity2D(halfParity1, halfParity2);

		if (mode == CALC_PARITY_MODE_PARITY) {
			mEZParity[bond] = bondDBParity;
			}
		else if (mode == CALC_PARITY_MODE_PRO_PARITY) {
			if (halfParity1.mRanksEqual) {
				if (bondDBParity == Molecule.cBondParityEor1) {
					mCanBase[halfParity1.mHighConn].add(mCanRank[dbAtom1]);
					}
				else if (bondDBParity == Molecule.cBondParityZor2) {
					mCanBase[halfParity1.mLowConn].add(mCanRank[dbAtom1]);
					}
				}
			if (halfParity2.mRanksEqual) {
				if (bondDBParity == Molecule.cBondParityEor1) {
					mCanBase[halfParity2.mHighConn].add(mCanRank[dbAtom2]);
					}
				else if (bondDBParity == Molecule.cBondParityZor2) {
					mCanBase[halfParity2.mLowConn].add(mCanRank[dbAtom2]);
					}
				}
			}

		return true;
		}


	private byte canCalcEZParity2D(EZHalfParity halfParity1, EZHalfParity halfParity2) {
		if (halfParity1.getValue() == -1 || halfParity2.getValue() == -1)
			return Molecule.cBondParityUnknown;

		if (((halfParity1.getValue() | halfParity2.getValue()) & 1) != 0)
			return Molecule.cBondParityUnknown;

		return (halfParity1.getValue() == halfParity2.getValue()) ?
				(byte)Molecule.cBondParityEor1 : Molecule.cBondParityZor2;
		}


	private byte canCalcEZParity3D(EZHalfParity halfParity1, EZHalfParity halfParity2) {
		double[] db = new double[3];
		db[0] = mMol.getAtomX(halfParity2.mCentralAxialAtom) - mMol.getAtomX(halfParity1.mCentralAxialAtom);
		db[1] = mMol.getAtomY(halfParity2.mCentralAxialAtom) - mMol.getAtomY(halfParity1.mCentralAxialAtom);
		db[2] = mMol.getAtomZ(halfParity2.mCentralAxialAtom) - mMol.getAtomZ(halfParity1.mCentralAxialAtom);

		double[] s1 = new double[3];
		s1[0] = mMol.getAtomX(halfParity1.mHighConn) - mMol.getAtomX(halfParity1.mCentralAxialAtom);
		s1[1] = mMol.getAtomY(halfParity1.mHighConn) - mMol.getAtomY(halfParity1.mCentralAxialAtom);
		s1[2] = mMol.getAtomZ(halfParity1.mHighConn) - mMol.getAtomZ(halfParity1.mCentralAxialAtom);

		double[] s2 = new double[3];
		s2[0] = mMol.getAtomX(halfParity2.mHighConn) - mMol.getAtomX(halfParity2.mCentralAxialAtom);
		s2[1] = mMol.getAtomY(halfParity2.mHighConn) - mMol.getAtomY(halfParity2.mCentralAxialAtom);
		s2[2] = mMol.getAtomZ(halfParity2.mHighConn) - mMol.getAtomZ(halfParity2.mCentralAxialAtom);

		// calculate the normal vector n1 of plane from db and s1 (vector product)
		double[] n1 = new double[3];
		n1[0] = db[1]*s1[2]-db[2]*s1[1];
		n1[1] = db[2]*s1[0]-db[0]*s1[2];
		n1[2] = db[0]*s1[1]-db[1]*s1[0];

		// calculate the normal vector n2 of plane from db and n1 (vector product)
		double[] n2 = new double[3];
		n2[0] = db[1]*n1[2]-db[2]*n1[1];
		n2[1] = db[2]*n1[0]-db[0]*n1[2];
		n2[2] = db[0]*n1[1]-db[1]*n1[0];

		// calculate cos(angle) of s1 and normal vector n2
		double cosa = (s1[0]*n2[0]+s1[1]*n2[1]+s1[2]*n2[2])
					/ (Math.sqrt(s1[0]*s1[0]+s1[1]*s1[1]+s1[2]*s1[2])
					 * Math.sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]));

		// calculate cos(angle) of s2 and normal vector n2
		double cosb = (s2[0]*n2[0]+s2[1]*n2[1]+s2[2]*n2[2])
					/ (Math.sqrt(s2[0]*s2[0]+s2[1]*s2[1]+s2[2]*s2[2])
					 * Math.sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]));

		return ((cosa < 0.0) ^ (cosb < 0.0)) ? (byte)Molecule.cBondParityEor1 : Molecule.cBondParityZor2;
		}


	private void flagStereoProblems() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			// if stereo center is declared unknown and no recognized as such; or vice versa
			if (mTHParity[atom] == Molecule.cAtomParityUnknown
			 && !mMol.isAtomConfigurationUnknown(atom))
				mMol.setStereoProblem(atom);

			// if no parity found, but atom was assigned to AND or OR group
			if ((mMol.getAtomESRType(atom) == Molecule.cESRTypeAnd
			  || mMol.getAtomESRType(atom) == Molecule.cESRTypeOr)
			 && (mTHParity[atom] == Molecule.cAtomParityUnknown))
				mMol.setStereoProblem(atom);

			if (mMol.isAtomConfigurationUnknown(atom)
			 && mTHParity[atom] != Molecule.cAtomParityUnknown
			 && !isUnknownBINAPBondAtom(atom))
				mMol.setStereoProblem(atom);
			}

		// if a stereo bond does not cause the recognition of a parity (tetrahedral, BINAP or allene type)
		for (int bond=0; bond<mMol.getAllBonds(); bond++)
			if (mMol.isStereoBond(bond) && !isJustifiedStereoBond(bond))
				mMol.setStereoProblem(mMol.getBondAtom(0,bond));

		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondOrder(bond) == 2) {
				// colorize double bonds with undefined EZParity
				// except crossbonds which are intentionally undefined
				if (mMol.isBondParityUnknownOrNone(bond)
				 && (mEZParity[bond] == Molecule.cBondParityEor1
				  || mEZParity[bond] == Molecule.cBondParityZor2)) {
					mEZParity[bond] = Molecule.cBondParityUnknown;
					mMol.setBondType(bond, ExtendedMolecule.cBondTypeCross);
					}
	
				if (mEZParity[bond] == Molecule.cBondParityUnknown
				 && !mEZParityIsPseudo[bond]) {
					if (mMol.getBondType(bond) != ExtendedMolecule.cBondTypeCross) {
						mMol.setStereoProblem(mMol.getBondAtom(0,bond));
						mMol.setStereoProblem(mMol.getBondAtom(1,bond));
						}
					}
				}

			if (mMol.getBondType(bond) == Molecule.cBondTypeSingle
			 && mEZParity[bond] == Molecule.cBondParityUnknown
			 && !mMol.isAtomConfigurationUnknown(mMol.getBondAtom(0, bond))
			 && !mMol.isAtomConfigurationUnknown(mMol.getBondAtom(1, bond))) {
				mMol.setStereoProblem(mMol.getBondAtom(0, bond));
				mMol.setStereoProblem(mMol.getBondAtom(1, bond));
				}

			if ((mMol.getBondESRType(bond) == Molecule.cESRTypeAnd
			  || mMol.getBondESRType(bond) == Molecule.cESRTypeOr)
			 && (mMol.getBondType(bond) != Molecule.cBondTypeSingle
			  || (mEZParity[bond] != Molecule.cBondParityEor1
			   && mEZParity[bond] != Molecule.cBondParityZor2))) {
				mMol.setStereoProblem(mMol.getBondAtom(0, bond));
				mMol.setStereoProblem(mMol.getBondAtom(1, bond));
				}
			}
		}


	private boolean isUnknownBINAPBondAtom(int atom) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (mEZParity[mMol.getConnBond(atom, i)] == Molecule.cBondParityUnknown
			 && mMol.getConnBondOrder(atom, i) == 1)
				return true;

		return false;
		}


	private boolean isJustifiedStereoBond(int bond) {
		int atom = mMol.getBondAtom(0, bond);
		if (atom >= mMol.getAtoms())
			return false;

		if (mTHParity[atom] == Molecule.cAtomParity1
		 || mTHParity[atom] == Molecule.cAtomParity2)
			return true;
		if (mTHParity[atom] == Molecule.cAtomParityUnknown)
			return false;

		int binapBond = mMol.findBINAPChiralityBond(atom);
		if (binapBond != -1)
			return mEZParity[binapBond] == Molecule.cBondParityEor1
				|| mEZParity[binapBond] == Molecule.cBondParityZor2;

		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			if (mMol.getConnBondOrder(atom, i) == 2) {
				if (mTHParity[mMol.getConnAtom(atom, i)] == Molecule.cAtomParity1
				 || mTHParity[mMol.getConnAtom(atom, i)] == Molecule.cAtomParity2)
					return true;
				}
			}

		return false;
		}


	private void generateGraph() {
		if (mMol.getAtoms() == 0)
			return;
		if (mGraphGenerated)
			return;

		mGraphRings = 0;

		int startAtom = 0;
		for (int atom=1; atom<mMol.getAtoms(); atom++)
			if (mCanRank[atom] > mCanRank[startAtom])
				startAtom = atom;

		boolean[] atomHandled = new boolean[mMol.getAtoms()];
		boolean[] bondHandled = new boolean[mMol.getBonds()];
		mGraphIndex = new int[mMol.getAtoms()];
		mGraphAtom = new int[mMol.getAtoms()];
		mGraphFrom = new int[mMol.getAtoms()];
		mGraphBond = new int[mMol.getBonds()];
		mGraphAtom[0] = startAtom;
		mGraphIndex[startAtom] = 0;
		atomHandled[startAtom] = true;

		int atomsWithoutParents = 1;	// the startatom has no parent
		int firstUnhandled = 0;
		int firstUnused = 1;
		int graphBonds = 0;
		while (firstUnhandled < mMol.getAtoms()) {
			if (firstUnhandled < firstUnused) {	// attach neighbours in rank order to unhandled
				while (true) {
					int highestRankingConnAtom = 0;
					int highestRankingConnBond = 0;
					int highestRank = -1;
					int atom = mGraphAtom[firstUnhandled];
					for (int i=0; i<mMol.getAllConnAtomsPlusMetalBonds(atom); i++) {
						if (i<mMol.getConnAtoms(atom) || i>=mMol.getAllConnAtoms(atom)) {
							int connAtom = mMol.getConnAtom(atom, i);
							if (!atomHandled[connAtom] && mCanRank[connAtom] > highestRank) {
								highestRankingConnAtom = connAtom;
								highestRankingConnBond = mMol.getConnBond(atom, i);
								highestRank = mCanRank[connAtom];
								}
							}
						}

					if (highestRank == -1)
						break;

					mGraphIndex[highestRankingConnAtom] = firstUnused;
					mGraphFrom[firstUnused] = firstUnhandled;
					mGraphAtom[firstUnused++] = highestRankingConnAtom;
					mGraphBond[graphBonds++] = highestRankingConnBond;
					atomHandled[highestRankingConnAtom] = true;
					bondHandled[highestRankingConnBond] = true;
					}
				firstUnhandled++;
				}
			else {
				int highestRankingAtom = 0;
				int highestRank = -1;
				for (int atom=0; atom<mMol.getAtoms(); atom++) {
					if (!atomHandled[atom] && mCanRank[atom] > highestRank) {
						highestRankingAtom = atom;
						highestRank = mCanRank[atom];
						}
					}
				atomsWithoutParents++;
				mGraphIndex[highestRankingAtom] = firstUnused;
				mGraphFrom[firstUnused] = -1;	// no parent atom in graph tree
				mGraphAtom[firstUnused++] = highestRankingAtom;
				atomHandled[highestRankingAtom] = true;
				}
			}

		mGraphClosure = new int[2 * (mMol.getBonds() - graphBonds)];

		while (true) {	// add ring closure bonds (those with lowest new atom numbers first)
			int lowAtomNo1 = mMol.getMaxAtoms();
			int lowAtomNo2 = mMol.getMaxAtoms();
			int lowBond = -1;
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				int loAtom,hiAtom;
				if (!bondHandled[bond]) {
					if (mGraphIndex[mMol.getBondAtom(0,bond)]
					  < mGraphIndex[mMol.getBondAtom(1,bond)]) {
						loAtom = mGraphIndex[mMol.getBondAtom(0,bond)];
						hiAtom = mGraphIndex[mMol.getBondAtom(1,bond)];
						}
					else {
						loAtom = mGraphIndex[mMol.getBondAtom(1,bond)];
						hiAtom = mGraphIndex[mMol.getBondAtom(0,bond)];
						}
					if (loAtom < lowAtomNo1
					 || (loAtom == lowAtomNo1 && hiAtom < lowAtomNo2)) {
						lowAtomNo1 = loAtom;
						lowAtomNo2 = hiAtom;
						lowBond = bond;
						}
					}
				}

			if (lowBond == -1)
				break;

			bondHandled[lowBond] = true;
			mGraphBond[graphBonds++] = lowBond;
			mGraphClosure[2*mGraphRings] = lowAtomNo1;
			mGraphClosure[2*mGraphRings+1] = lowAtomNo2;
			mGraphRings++;
			}

		mGraphGenerated = true;
		}


	public StereoMolecule getCanMolecule() {
		return getCanMolecule(false);
		}


	public StereoMolecule getCanMolecule(boolean includeExplicitHydrogen) {
		generateGraph();

		StereoMolecule mol = new StereoMolecule(mMol.getAtoms(), mMol.getBonds());

		mol.setFragment(mMol.isFragment());	// to allow copying of atom/bond query features

		for(int i=0; i<mMol.getAtoms(); i++) {
			mMol.copyAtom(mol, mGraphAtom[i], 0, 0);
			mol.setAtomESR(i, mTHESRType[mGraphAtom[i]], mTHESRGroup[mGraphAtom[i]]);
			}

		for(int i=0; i<mMol.getBonds(); i++) {
			mMol.copyBond(mol, mGraphBond[i], 0, 0, mGraphIndex, false);
			if (!mol.isStereoBond(i) && mol.getBondAtom(0, i) > mol.getBondAtom(1, i)) {
				int temp = mol.getBondAtom(0, i);
				mol.setBondAtom(0, i, mol.getBondAtom(1, i));
				mol.setBondAtom(1, i, temp);
				}
			mol.setBondESR(i, mEZESRType[mGraphBond[i]], mEZESRGroup[mGraphBond[i]]);
			}

		if (includeExplicitHydrogen) {
			for (int i=0; i<mMol.getAtoms(); i++) {
				int atom = mGraphAtom[i];
				for (int j=mMol.getConnAtoms(atom); j<mMol.getAllConnAtoms(atom); j++) {
					int hydrogen = mMol.copyAtom(mol, mMol.getConnAtom(atom, j), 0, 0);
					mMol.copyBond(mol, mMol.getConnBond(atom, j), 0, 0, mGraphIndex[atom], hydrogen, false);
					}
				}
			}

		for (int bond=0; bond<mol.getAllBonds(); bond++) {
			int atom = mol.getBondAtom(0, bond);
			if (mTHParityIsMesoInverted[mGraphAtom[atom]]) {
				if (mol.getBondType(bond) == Molecule.cBondTypeUp)
					mol.setBondType(bond, Molecule.cBondTypeDown);
				else if (mol.getBondType(bond) == Molecule.cBondTypeDown)
					mol.setBondType(bond, Molecule.cBondTypeUp);
				}
			}

		mMol.copyMoleculeProperties(mol);
		mMol.invalidateHelperArrays(Molecule.cHelperBitParities);

		return mol;
		}


	/**
	 * Sets all atoms with TH-parity 'unknown' to explicitly defined 'unknown'.
	 * Sets all double bonds with EZ-parity 'unknown' to cross bonds.
	 * Sets the first bond atom of all BINAP type bonds with parity 'unknown' to explicitly defined 'unknown' parity.
	 */
	public void setUnknownParitiesToExplicitlyUnknown() {
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (!mMol.isAtomConfigurationUnknown(atom)
			 && mTHParity[atom] == Molecule.cAtomParityUnknown)
				mMol.setAtomConfigurationUnknown(atom, true);
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mEZParity[bond] == Molecule.cBondParityUnknown) {
				int order = mMol.getBondOrder(bond);
				if (order == 1) {
					mMol.setAtomConfigurationUnknown(mMol.getBondAtom(0, bond), true);
					}
				else if (order == 2) {
			   		mMol.setBondType(bond, Molecule.cBondTypeCross);
					}
				}
			}
		}


	/**
	 * If the molecule contains exactly one stereo center and if that has unknown configuration,
	 * than assume that the configuration is meant to be racemic and update molecule accordingly.
	 * If stereo configuration is ill defined with a stereo bond whose pointed tip is not at the
	 * stereo center, then the molecule is not touched and the stereo center kept as undefined.
	 * @return whether a stereo center was converted to be racemic
	 */
	public boolean setSingleUnknownAsRacemicParity() {
		int unknownTHParities = 0;
		int knownTHParities = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mTHParity[atom] != Molecule.cAtomParityNone
			 && !mTHParityIsPseudo[atom]) {
				if (mTHParity[atom] == Molecule.cAtomParityUnknown)
					unknownTHParities++;
				else
					knownTHParities++;
				}
			}
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondType(bond) == Molecule.cBondTypeSingle
			 && mEZParity[bond] != Molecule.cAtomParityNone
			 && !mEZParityIsPseudo[bond]) {
				if (mEZParity[bond] == Molecule.cAtomParityUnknown)
					unknownTHParities++;
				else
					knownTHParities++;
				}
			}

		if (knownTHParities == 0 && unknownTHParities == 1) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mTHParity[atom] == Molecule.cAtomParityUnknown
				 && !mTHParityIsPseudo[atom]) {
					// The parity may be unknown, because we have inconclusive
					// stereobonds, e.g. drawn in a wrong direction, or more
					// than one stereo bond in an inconsistent way.
					// In this case we keep the parity as unknown.
					if (mMol.getAtomPi(atom) == 2 && mMol.getConnAtoms(atom) == 2) {
						// chiral allene case
						for (int i=0; i<2; i++) {
							int connAtom = mMol.getConnAtom(atom, i);
							for (int j=0; j<mMol.getConnAtoms(connAtom); j++)
								if (mMol.isStereoBond(mMol.getConnBond(connAtom, j)))
									return false;
							}
						}
					else {
						// tetrahedral chirality
						for (int i=0; i<mMol.getConnAtoms(atom); i++)
							if (mMol.isStereoBond(mMol.getConnBond(atom, i)))
								return false;
						}

						// This must be parity2 (rather than 1) because
						// group parities are normalized at this stage
						// and parity2 is the normalized parity of the
						// highest ranking atom in a ESR group.
					mTHParity[atom] = Molecule.cAtomParity2;
					mTHESRType[atom] = Molecule.cESRTypeAnd;
					mTHESRGroup[atom] = 0;

					mMol.setAtomParity(atom, Molecule.cAtomParity2, false);
					mMol.setAtomESR(atom, Molecule.cESRTypeAnd, 0);

					int stereoBond = mMol.getAtomPreferredStereoBond(atom);
					mMol.setBondType(stereoBond, Molecule.cBondTypeUp);
					if (mMol.getBondAtom(1, stereoBond) == atom) {
						int connAtom = mMol.getBondAtom(0, stereoBond);
						mMol.setBondAtom(0, stereoBond, atom);
						mMol.setBondAtom(1, stereoBond, connAtom);
						}
					return true;
					}
				}
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (mMol.getBondType(bond) == Molecule.cBondTypeSingle
				 && mEZParity[bond] != Molecule.cAtomParityNone
				 && !mEZParityIsPseudo[bond]) {
					// The parity may be unknown, because we have inconclusive
					// stereobonds, e.g. drawn in a wrong direction, or more
					// than one stereo bond in an inconsistent way.
					// In this case we keep the parity as unknown.

					// BINAP chirality bond case
					for (int i=0; i<2; i++) {
						int atom = mMol.getBondAtom(i, bond);
						for (int j=0; j<mMol.getConnAtoms(atom); j++)
							if (mMol.isStereoBond(mMol.getConnBond(atom, j)))
								return false;
						}
	
						// This must be parity2 (rather than 1) because
						// group parities are normalized at this stage
						// and parity2 is the normalized parity of the
						// highest ranking atom in a ESR group.
					mEZParity[bond] = Molecule.cBondParityZor2;
					mEZESRType[bond] = Molecule.cESRTypeAnd;
					mEZESRGroup[bond] = 0;
	
					mMol.setBondParity(bond, Molecule.cBondParityZor2, false);
					mMol.setAtomESR(bond, Molecule.cESRTypeAnd, 0);
	
					int stereoBond = mMol.getBondPreferredStereoBond(bond);
					mMol.setBondType(stereoBond, Molecule.cBondTypeUp);
					if (mMol.getBondAtom(1, stereoBond) == mMol.getBondAtom(0, bond)
					 || mMol.getBondAtom(1, stereoBond) == mMol.getBondAtom(1, bond)) {
						int connAtom = mMol.getBondAtom(0, stereoBond);
						mMol.setBondAtom(0, stereoBond, mMol.getBondAtom(1, stereoBond));
						mMol.setBondAtom(1, stereoBond, connAtom);
						}
					return true;
					}
				}
			}
		return false;
		}


	public String getIDCode() {
		if (mIDCode == null) {
			generateGraph();
			if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0) {
				idGenerateConfigurations();
				idNormalizeESRGroupNumbers();
				}
			idCodeCreate();
			}

		return mIDCode;
		}


	public int[] getFinalRank() {
		// this is the final mCanRank after all tie breaking steps
		return mCanRank;
		}


	/**
	 * Returns the symmetry rank before tie breaking. For this the Canonizer
	 * mode must contain the CREATE_SYMMETRY_RANK option. If ranking
	 * shall reflect atom diastereotopicity or even enantiotopicity, use
	 * mode CONSIDER_DIASTEREOTOPICITY or CONSIDER_STEREOHETEROTOPICITY,
	 * respectively.
	 * @param atom
	 * @return rank
	 */
	public int getSymmetryRank(int atom) {
		return (mCanRankBeforeTieBreaking == null) ? -1 : mCanRankBeforeTieBreaking[atom];
		}


	/**
	 * Returns the symmetry ranks before tie breaking. For this the Canonizer
	 * mode must contain the CREATE_SYMMETRY_RANK option. If ranking
	 * shall reflect atom diastereotopicity or even enantiotopicity, use
	 * mode CONSIDER_DIASTEREOTOPICITY or CONSIDER_STEREOHETEROTOPICITY,
	 * respectively.
	 * @return ranks
	 */
	public int[] getSymmetryRanks() {
		return mCanRankBeforeTieBreaking;
		}


	private void idCodeCreate() {
		encodeBitsStart(false);
		encodeBits(cIDCodeVersion3, 4);
		int nbits = Math.max(getNeededBits(mMol.getAtoms()),
							 getNeededBits(mMol.getBonds()));
		encodeBits(nbits, 4);

		if (nbits == 0) {
			encodeBits(mMol.isFragment() ? 1 : 0, 1);	// query fragment ?
			encodeBits(0, 1);
			mIDCode = encodeBitsEnd();
			return;
			}

		int nitrogens,oxygens,otherAtoms,chargedAtoms;
		nitrogens = oxygens = otherAtoms = chargedAtoms = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) == 0) {
				switch (mMol.getAtomicNo(atom)) {
				case 6:
					break;
				case 7:
					nitrogens++;
					break;
				case 8:
					oxygens++;
					break;
				default:
					otherAtoms++;
					break;
					}
				if (mMol.getAtomCharge(atom) != 0)
					chargedAtoms++;
				}
			}

		encodeBits(mMol.getAtoms(), nbits);
		encodeBits(mMol.getBonds(), nbits);
		encodeBits(nitrogens, nbits);
		encodeBits(oxygens, nbits);
		encodeBits(otherAtoms, nbits);
		encodeBits(chargedAtoms, nbits);

/*
System.out.println("Atoms:" + mMol.getAtoms());
System.out.println("Bonds:" + mMol.getBonds());
System.out.print("Nitrogens:");
for (int atom=0; atom<mMol.getAtoms(); atom++)
if (mMol.getAtomicNo(mGraphAtom[atom]) == 7)
System.out.print(atom + ";");
System.out.println();
System.out.print("Oxygens:");
for (int atom=0; atom<mMol.getAtoms(); atom++)
if (mMol.getAtomicNo(mGraphAtom[atom]) == 8)
System.out.print(atom + ";");
System.out.println();
System.out.print("Others:");
for (int atom=0; atom<mMol.getAtoms(); atom++)
if (mMol.getAtomicNo(mGraphAtom[atom]) != 6
 && mMol.getAtomicNo(mGraphAtom[atom]) != 7
 && mMol.getAtomicNo(mGraphAtom[atom]) != 8)
System.out.print("(" + atom + "," + mMol.getAtomicNo(mGraphAtom[atom]) + ");");
System.out.println();
System.out.print("Charged:");
for (int atom=0; atom<mMol.getAtoms(); atom++)
if (mMol.getAtomCharge(mGraphAtom[atom]) != 0)
System.out.print("(" + atom + "," + mMol.getAtomCharge(mGraphAtom[atom]) + ");");
System.out.println();
*/

		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomicNo(mGraphAtom[atom]) == 7
			 && (mMol.getAtomQueryFeatures(mGraphAtom[atom]) & Molecule.cAtomQFAny) == 0)
				encodeBits(atom, nbits);
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomicNo(mGraphAtom[atom]) == 8
			 && (mMol.getAtomQueryFeatures(mGraphAtom[atom]) & Molecule.cAtomQFAny) == 0)
				encodeBits(atom, nbits);
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomicNo(mGraphAtom[atom]) != 6
			 && mMol.getAtomicNo(mGraphAtom[atom]) != 7
			 && mMol.getAtomicNo(mGraphAtom[atom]) != 8
			 && (mMol.getAtomQueryFeatures(mGraphAtom[atom]) & Molecule.cAtomQFAny) == 0) {
				encodeBits(atom, nbits);
				encodeBits(mMol.getAtomicNo(mGraphAtom[atom]), 8);
				}
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomCharge(mGraphAtom[atom]) != 0
			 && (mMol.getAtomQueryFeatures(mGraphAtom[atom]) & Molecule.cAtomQFAny) == 0) {
				encodeBits(atom, nbits);
				encodeBits(8 + mMol.getAtomCharge(mGraphAtom[atom]), 4);
				}

		int maxdif = 0;
		int base = 0;
		for (int atom=1; atom<mMol.getAtoms(); atom++) {
			int dif;
			if (mGraphFrom[atom] == -1) {
				dif = 0;
				}
			else {
				dif = 1 + mGraphFrom[atom] - base;
				base = mGraphFrom[atom];
				}
			if (maxdif < dif)
				maxdif = dif;
			}
/*
System.out.println("D-Bits:" + idGetNeededBits(maxdif));
System.out.print("GraphAtomDifs:");
*/
		int dbits = getNeededBits(maxdif);
		encodeBits(dbits, 4);
		base = 0;
		for (int atom=1; atom<mMol.getAtoms(); atom++) {
			int dif;
			if (mGraphFrom[atom] == -1) {
				dif = 0;
				}
			else {
				dif = 1 + mGraphFrom[atom] - base;
				base = mGraphFrom[atom];
				}
			encodeBits(dif, dbits);

//System.out.print(dif + ";");
			}
//System.out.println();

		for (int i=0; i<2*mGraphRings; i++)
			encodeBits(mGraphClosure[i], nbits);
/*
System.out.print("ClosureAtoms:");
for (int i=0; i<2*mGraphRings; i++)
System.out.print(mGraphClosure[i] + ";");
System.out.println();
System.out.print("BondOrders:");
*/
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			int bondOrder = ((mMol.getBondQueryFeatures(mGraphBond[bond]) & Molecule.cBondQFBridge) != 0
						  || mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeMetalLigand) ?
							1 : (mMol.isDelocalizedBond(mGraphBond[bond])) ?
							0 : Math.min(3, mMol.getBondOrder(mGraphBond[bond]));
			encodeBits(bondOrder, 2);
//System.out.print(bondOrder + ";");
			}
//System.out.println();

		int THCount = 0;
		if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0) {
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (mTHConfiguration[mGraphAtom[atom]] != Molecule.cAtomParityNone
				 && mTHConfiguration[mGraphAtom[atom]] != Molecule.cAtomParityUnknown)
					THCount++;
			}
		encodeBits(THCount, nbits);
		if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mTHConfiguration[mGraphAtom[atom]] != Molecule.cAtomParityNone
				 && mTHConfiguration[mGraphAtom[atom]] != Molecule.cAtomParityUnknown) {
					encodeBits(atom, nbits);
					if (mTHESRType[mGraphAtom[atom]] == Molecule.cESRTypeAbs) {
						encodeBits(mTHConfiguration[mGraphAtom[atom]], 3);
						}
					else {
						int parity = (mTHConfiguration[mGraphAtom[atom]] == Molecule.cAtomParity1) ?
								((mTHESRType[mGraphAtom[atom]] == Molecule.cESRTypeAnd) ?
										cParity1And : cParity1Or)
							  : ((mTHESRType[mGraphAtom[atom]] == Molecule.cESRTypeAnd) ?
										cParity2And : cParity2Or);

						encodeBits(parity, 3);
						encodeBits(mTHESRGroup[mGraphAtom[atom]], 3);
						}
					}
				}
			}

		int EZCount = 0;
		if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0) {
			for (int bond=0; bond<mMol.getBonds(); bond++)	// parity of all double bonds
				if (mEZConfiguration[mGraphBond[bond]] != 0
				 && mEZConfiguration[mGraphBond[bond]] != Molecule.cBondParityUnknown
				 && (!mMol.isSmallRingBond(mGraphBond[bond]) || mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeSingle))
					EZCount++;
			}
		encodeBits(EZCount, nbits);
		if ((mMode & NEGLECT_ANY_STEREO_INFORMATION) == 0) {
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (mEZConfiguration[mGraphBond[bond]] != 0
				 && mEZConfiguration[mGraphBond[bond]] != Molecule.cBondParityUnknown
				 && (!mMol.isSmallRingBond(mGraphBond[bond]) || mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeSingle)) {
					encodeBits(bond, nbits);
					if (mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeSingle) {	// BINAP type of axial chirality
						if (mEZESRType[mGraphBond[bond]] == Molecule.cESRTypeAbs) {
							encodeBits(mEZConfiguration[mGraphBond[bond]], 3);
							}
						else {
							int parity = (mEZConfiguration[mGraphBond[bond]] == Molecule.cBondParityEor1) ?
									((mEZESRType[mGraphBond[bond]] == Molecule.cESRTypeAnd) ?
											cParity1And : cParity1Or)
								  : ((mEZESRType[mGraphBond[bond]] == Molecule.cESRTypeAnd) ?
											cParity2And : cParity2Or);

							encodeBits(parity, 3);
							encodeBits(mEZESRGroup[mGraphBond[bond]], 3);
							}
						}
					else {
						encodeBits(mEZConfiguration[mGraphBond[bond]], 2);	// for compatibility reasons we use only two bits for double bonds
						}
					}
				}
			}
/*
THCount = 0;
for (int atom=0; atom<mMol.getAtoms(); atom++)
if (mTHConfiguration[mGraphAtom[atom]] != Molecule.cAtomParityNone
&& mTHConfiguration[mGraphAtom[atom]] != Molecule.cAtomParityUnknown)
THCount++;
System.out.print("THCount:"+THCount+"  ");
for (int atom=0; atom<mMol.getAtoms(); atom++)
if (mTHConfiguration[mGraphAtom[atom]] != Molecule.cAtomParityNone
&& mTHConfiguration[mGraphAtom[atom]] != Molecule.cAtomParityUnknown) {
if (mMol.getAtomESRType(mGraphAtom[atom]) == Molecule.cAtomESRTypeAbs) {
System.out.print("("+atom+":"+mTHConfiguration[mGraphAtom[atom]]+")");
}
else {
int parity = (mTHConfiguration[mGraphAtom[atom]] == Molecule.cAtomParity1) ?
((mMol.getAtomESRType(mGraphAtom[atom]) == Molecule.cAtomESRTypeAnd) ?
cTHParity1And : cTHParity1Or)
: ((mMol.getAtomESRType(mGraphAtom[atom]) == Molecule.cAtomESRTypeAnd) ?
cTHParity2And : cTHParity2Or);
System.out.print("("+atom+":"+parity+"-"+mTHESRGroup[mGraphAtom[atom]]+")");
}}
System.out.println();

EZCount = 0;
for (int bond=0; bond<mMol.getBonds(); bond++)	// parity of all double bonds
if (mEZConfiguration[mGraphBond[bond]] != 0
&& mEZConfiguration[mGraphBond[bond]] != Molecule.cBondParityUnknown
&& !mMol.isSmallRingBond(mGraphBond[bond]))
EZCount++;
System.out.print("EZCount:"+EZCount+"  ");
for (int bond=0; bond<mMol.getBonds(); bond++)
if (mEZConfiguration[mGraphBond[bond]] != 0
&& mEZConfiguration[mGraphBond[bond]] != Molecule.cBondParityUnknown
&& !mMol.isSmallRingBond(mGraphBond[bond])) {
System.out.print("("+bond+":"+mEZConfiguration[mGraphBond[bond]]+")");
}
System.out.println();
*/

		encodeBits(mMol.isFragment() ? 1 : 0, 1);	// query fragment ?

		int count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomMass(mGraphAtom[atom]) != 0)
				count++;
		if (count != 0) {
			encodeFeatureNo(1);	//	1 = datatype 'isotope'
			encodeBits(count, nbits);
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mMol.getAtomMass(mGraphAtom[atom]) != 0) {
					encodeBits(atom, nbits);
					encodeBits(mMol.getAtomMass(mGraphAtom[atom]), 8);
					}
				}
			}

		mFeatureBlock = 0;

		if (mMol.isFragment()) {	// QueryFeatures and fragment specific properties
			addAtomQueryFeatures(0, nbits, Molecule.cAtomQFNoMoreNeighbours, 1, -1);

			addAtomQueryFeatures(3, nbits, Molecule.cAtomQFMoreNeighbours, 1, -1);

			addAtomQueryFeatures(4, nbits,
								 Molecule.cAtomQFRingState,
								 Molecule.cAtomQFRingStateBits,
								 Molecule.cAtomQFRingStateShift);

			addAtomQueryFeatures(5, nbits,
								 Molecule.cAtomQFAromState,
								 Molecule.cAtomQFAromStateBits,
								 Molecule.cAtomQFAromStateShift);

			addAtomQueryFeatures(6, nbits, Molecule.cAtomQFAny, 1, -1);

			addAtomQueryFeatures(7, nbits,
								 Molecule.cAtomQFHydrogen,
								 Molecule.cAtomQFHydrogenBits,
								 Molecule.cAtomQFHydrogenShift);

			count = 0;
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (mMol.getAtomList(mGraphAtom[atom]) != null)
					count++;
			if (count > 0) {
				encodeFeatureNo(8);	//	8 = datatype 'AtomList'
				encodeBits(count, nbits);
				for (int atom=0; atom<mMol.getAtoms(); atom++) {
					int[] atomList = mMol.getAtomList(mGraphAtom[atom]);
					if (atomList != null) {
						encodeBits(atom, nbits);
						encodeBits(atomList.length, 4);
						for (int a:atomList)
							encodeBits(a, 8);
						}
					}
				}

			addBondQueryFeatures(9, nbits,
								 Molecule.cBondQFRingState,
								 Molecule.cBondQFRingStateBits,
								 Molecule.cBondQFRingStateShift);

			addBondQueryFeatures(10, nbits,
								 Molecule.cBondQFBondTypes,
								 Molecule.cBondQFBondTypesBits,
								 Molecule.cBondQFBondTypesShift);

			addAtomQueryFeatures(11, nbits, Molecule.cAtomQFMatchStereo, 1, -1);

			addBondQueryFeatures(12, nbits,
								 Molecule.cBondQFBridge,
								 Molecule.cBondQFBridgeBits,
								 Molecule.cBondQFBridgeShift);

			addAtomQueryFeatures(13, nbits,
								 Molecule.cAtomQFPiElectrons,
								 Molecule.cAtomQFPiElectronBits,
								 Molecule.cAtomQFPiElectronShift);

			addAtomQueryFeatures(14, nbits,
								 Molecule.cAtomQFNeighbours,
								 Molecule.cAtomQFNeighbourBits,
								 Molecule.cAtomQFNeighbourShift);

			addAtomQueryFeatures(16, nbits,
								 Molecule.cAtomQFSmallRingSize,
								 Molecule.cAtomQFSmallRingSizeBits,
								 Molecule.cAtomQFSmallRingSizeShift);
			}

		count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mAbnormalValence != null && mAbnormalValence[mGraphAtom[atom]] != -1)
				count++;
		if (count != 0) {
			encodeFeatureNo(17);   //  17 = datatype 'AtomAbnormalValence'
			encodeBits(count, nbits);
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mAbnormalValence != null && mAbnormalValence[mGraphAtom[atom]] != -1) {
					encodeBits(atom, nbits);
					encodeBits(mAbnormalValence[mGraphAtom[atom]], 4);
					}
				}
			}

		if ((mMode & ENCODE_ATOM_CUSTOM_LABELS) != 0
		 ||	(mMode & ENCODE_ATOM_CUSTOM_LABELS_WITHOUT_RANKING) != 0) {
			count = 0;
			int maxLength = 0;
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				String label = mMol.getAtomCustomLabel(mGraphAtom[atom]);
				if (label != null) {
					count++;
					maxLength = Math.max(maxLength, label.length());
					}
				}
			if (count != 0) {
				int lbits = getNeededBits(maxLength);
				encodeFeatureNo(18);   //  18 = datatype 'AtomCustomLabel'
				encodeBits(count, nbits);
				encodeBits(lbits, 4);
				for (int atom=0; atom<mMol.getAtoms(); atom++) {
					String customLabel = mMol.getAtomCustomLabel(mGraphAtom[atom]);
					if (customLabel != null) {
						encodeBits(atom, nbits);
						encodeBits(customLabel.length(), lbits);
						for (int i=0; i<customLabel.length(); i++)
							encodeBits(customLabel.charAt(i), 7);
						}
					}
				}
			}

		if (mMol.isFragment()) {	// more QueryFeatures and fragment specific properties
			addAtomQueryFeatures(19, nbits,
					 			 Molecule.cAtomQFCharge,
					 			 Molecule.cAtomQFChargeBits,
					 			 Molecule.cAtomQFChargeShift);

			addBondQueryFeatures(20, nbits,
								 Molecule.cBondQFRingSize,
								 Molecule.cBondQFRingSizeBits,
								 Molecule.cBondQFRingSizeShift);
			}

		count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomRadical(mGraphAtom[atom]) != 0)
				count++;
		if (count != 0) {
			encodeFeatureNo(21);   //  21 = datatype 'AtomRadicalState'
			encodeBits(count, nbits);
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mMol.getAtomRadical(mGraphAtom[atom]) != 0) {
					encodeBits(atom, nbits);
					encodeBits(mMol.getAtomRadical(mGraphAtom[atom]) >> Molecule.cAtomRadicalStateShift, 2);
					}
				}
			}

		if (mMol.isFragment()) {	// more QueryFeatures and fragment specific properties
			addAtomQueryFeatures(22, nbits, Molecule.cAtomQFFlatNitrogen, 1, -1);
			addBondQueryFeatures(23, nbits, Molecule.cBondQFMatchStereo, 1, -1);
			addBondQueryFeatures(24, nbits,
					 							Molecule.cBondQFAromState,
					 							Molecule.cBondQFAromStateBits,
					 							Molecule.cBondQFAromStateShift);
			}

		if ((mMode & ENCODE_ATOM_SELECTION) != 0) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mMol.isSelectedAtom(mGraphAtom[atom])) {
					encodeFeatureNo(25);   //  25 = datatype 'AtomSelection'
					for (int a=0; a<mMol.getAtoms(); a++)
						encodeBits(mMol.isSelectedAtom(mGraphAtom[a]) ? 1:0, 1);
					break;
					}
				}
			}

		boolean[] isAromaticSPBond = getAromaticSPBonds();
		if (isAromaticSPBond != null) {
			count = 0;
			for (int bond=0; bond<mMol.getBonds(); bond++)
				if (isAromaticSPBond[mGraphBond[bond]])
					count++;

			encodeFeatureNo(26);   //  26 = datatype 'delocalized high order bond'
			encodeBits(count, nbits);
			for (int bond=0; bond<mMol.getBonds(); bond++)
				if (isAromaticSPBond[mGraphBond[bond]])
					encodeBits(bond, nbits);
			}

		if (mMol.isFragment())	// 27 = datatype 'part of an exclude-group'
			addAtomQueryFeatures(27, nbits, Molecule.cAtomQFExcludeGroup, 1, -1);

		count = 0;
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeMetalLigand)
				count++;
		if (count != 0) {
			encodeFeatureNo(28);    // 28 = datatype 'dative (0-order) bond'
			encodeBits(count, nbits);
			for (int bond=0; bond<mMol.getBonds(); bond++)
				if (mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeMetalLigand)
					encodeBits(bond, nbits);
			}

		if (mMol.isFragment()) {    // 29 = datatype 'reaction parity hint'
			addAtomQueryFeatures(29, nbits, Molecule.cAtomQFRxnParityHint, Molecule.cAtomQFRxnParityBits, Molecule.cAtomQFRxnParityShift);
			addAtomQueryFeatures(30, nbits, Molecule.cAtomQFNewRingSize, Molecule.cAtomQFNewRingSizeBits, Molecule.cAtomQFNewRingSizeShift);
			addAtomQueryFeatures(32, nbits, Molecule.cAtomQFStereoState, Molecule.cAtomQFStereoStateBits, Molecule.cAtomQFStereoStateShift);
			addAtomQueryFeatures(33, nbits, Molecule.cAtomQFENeighbours, Molecule.cAtomQFENeighbourBits, Molecule.cAtomQFENeighbourShift);
			addAtomQueryFeatures(34, nbits, Molecule.cAtomQFHeteroAromatic, 1, -1);
			addBondQueryFeatures(35, nbits, Molecule.cBondQFMatchFormalOrder, 1, -1);
			addBondQueryFeatures(36, nbits, Molecule.cBondQFRareBondTypes, Molecule.cBondQFRareBondTypesBits, Molecule.cBondQFRareBondTypesShift);
			}

		count = 0;
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if (mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeQuadruple
			 || mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeQuintuple)
				count++;
		if (count != 0) {
			encodeFeatureNo(37);    // 37 = datatype 'rare order bond'
			encodeBits(count, nbits);
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeQuadruple
				 || mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeQuintuple) {
					encodeBits(bond, nbits);
					encodeBits(mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeQuadruple ? 0 : 1, 1);
					}
				}
			}

		encodeBits(0, 1);
		mIDCode = encodeBitsEnd();
		}


	private void addAtomQueryFeatures(int codeNo, int nbits, long qfMask, int qfBits, int qfShift) {
		int count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if ((mMol.getAtomQueryFeatures(mGraphAtom[atom]) & qfMask) != 0)
				count++;

		if (count == 0)
			return;

		encodeFeatureNo(codeNo);
		encodeBits(count, nbits);
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			long feature = mMol.getAtomQueryFeatures(mGraphAtom[atom]) & qfMask;
			if (feature != 0) {
				encodeBits(atom, nbits);
				if (qfBits != 1)
					encodeBits(feature >> qfShift, qfBits);
				}
			}
		}


	private void addBondQueryFeatures(int codeNo, int nbits, int qfMask, int qfBits, int qfShift) {
		int count = 0;
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if ((mMol.getBondQueryFeatures(mGraphBond[bond]) & qfMask) != 0)
				count++;

		if (count == 0)
			return;

		encodeFeatureNo(codeNo);
		encodeBits(count, nbits);
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			int feature = mMol.getBondQueryFeatures(mGraphBond[bond]) & qfMask;
			if (feature != 0) {
				encodeBits(bond, nbits);
				if (qfBits != 1)
					encodeBits(feature >> qfShift, qfBits);
				}
			}
		}


	private boolean[] getAromaticSPBonds() {
		boolean[] isAromaticSPBond = null;
		RingCollection ringSet = mMol.getRingSet();
		for (int r=0; r<ringSet.getSize(); r++) {
			if (ringSet.isDelocalized(r)) {
				int count = 0;
				int[] ringAtom = ringSet.getRingAtoms(r);
				for (int atom:ringAtom)
					if (hasTwoAromaticPiElectrons(atom))
						count++;
				if (count != 0) {
					int[] ringBond = ringSet.getRingBonds(r);
					if (isAromaticSPBond == null)
						isAromaticSPBond = new boolean[mMol.getBonds()];
					if (count == ringAtom.length) {
						// We have two alternatives for marking every second bond in the ring
						// We choose the one with the bond with the lowest graphIndex.
						int minIndex = -1;
						int minValue = Integer.MAX_VALUE;
						for (int i=0; i<ringAtom.length; i++) {
							if (minValue > mGraphAtom[ringBond[i]]) {
								minValue = mGraphAtom[ringBond[i]];
								minIndex = i;
								}
							}
						while (count > 0) {
							isAromaticSPBond[ringBond[minIndex]] = true;
							minIndex = validateCyclicIndex(minIndex+2, ringAtom.length);
							count -= 2;
							}
						}
					else {
						int index = 0;
						while (hasTwoAromaticPiElectrons(ringAtom[index]))
							index++;
						while (!hasTwoAromaticPiElectrons(ringAtom[index]))
							index = validateCyclicIndex(index+1, ringAtom.length);
						while (count > 0) {
							isAromaticSPBond[ringBond[index]] = true;
							index = validateCyclicIndex(index+2, ringAtom.length);
							count -= 2;
							while (!hasTwoAromaticPiElectrons(ringAtom[index]))
								index = validateCyclicIndex(index+1, ringAtom.length);
							}
						}
					}
				}
			}
		return isAromaticSPBond;
		}


	private boolean hasTwoAromaticPiElectrons(int atom) {
		if (mMol.getAtomPi(atom) < 2)
			return false;
		if (mMol.getConnAtoms(atom) == 2)
			return true;
		int aromaticPi = 0;
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int connBond = mMol.getConnBond(atom, i);
			if (mMol.isAromaticBond(connBond))
				aromaticPi += mMol.getBondOrder(connBond) - 1;
			}
		return aromaticPi > 1;
		}


	private int validateCyclicIndex(int index, int limit) {
		return (index < limit) ? index : index - limit;
		}

	public void invalidateCoordinates() {
		mEncodedCoords = null;
		}

	/**
	 * Encodes the molecule's atom coordinates into a compact String. Together with the
	 * idcode the coordinate string can be passed to the IDCodeParser to recreate the
	 * original molecule including coordinates.<br>
	 * If the molecule's coordinates are 2D, then coordinate encoding will be relative,
	 * i.e. scale and absolute positions get lost during the encoding.
	 * 3D-coordinates, however, are encoded retaining scale and absolute positions.<br>
	 * If the molecule has 3D-coordinates and if there are no implicit hydrogen atoms,
	 * i.e. all hydrogen atoms are explicitly available with their coordinates, then
	 * hydrogen 3D-coordinates are also encoded despite the fact that the idcode itself does
	 * not contain hydrogen atoms, because it must be canonical.
	 * @return
	 */
	public String getEncodedCoordinates() {
		return getEncodedCoordinates(mZCoordinatesAvailable);
		}

	/**
	 * Encodes the molecule's atom coordinates into a compact String. Together with the
	 * idcode the coordinate string can be passed to the IDCodeParser to recreate the
	 * original molecule including coordinates.<br>
	 * If keepPositionAndScale==false, then coordinate encoding will be relative,
	 * i.e. scale and absolute positions get lost during the encoding.
	 * Otherwise the encoding retains scale and absolute positions.<br>
	 * If the molecule has 3D-coordinates and if there are no implicit hydrogen atoms,
	 * i.e. all hydrogen atoms are explicitly available with their coordinates, then
	 * hydrogen 3D-coordinates are also encoded despite the fact that the idcode itself does
	 * not contain hydrogen atoms, because it must be canonical.
	 * @param keepPositionAndScale if false, then coordinates are scaled to an average bond length of 1.5 units
	 * @return
	 */
	public String getEncodedCoordinates(boolean keepPositionAndScale) {
		if (mEncodedCoords == null) {
			generateGraph();
			encodeCoordinates(keepPositionAndScale, mMol.getAtomCoordinates());
			}

		return mEncodedCoords;
		}

	/**
	 * Encodes the molecule's atom coordinates into a compact String. Together with the
	 * idcode the coordinate string can be passed to the IDCodeParser to recreate the
	 * original molecule including coordinates.<br>
	 * If keepPositionAndScale==false, then coordinate encoding will be relative,
	 * i.e. scale and absolute positions get lost during the encoding.
	 * Otherwise the encoding retains scale and absolute positions.<br>
	 * If the molecule has 3D-coordinates and if there are no implicit hydrogen atoms,
	 * i.e. all hydrogen atoms are explicitly available with their coordinates, then
	 * hydrogen 3D-coordinates are also encoded despite the fact that the idcode itself does
	 * not contain hydrogen atoms, because it must be canonical.
	 * @param keepPositionAndScale if false, then coordinates are scaled to an average bond length of 1.5 units
	 * @param atomCoordinates external atom coordinate set for the same molecule, e.g. from a Conformer
	 * @return
	 */
	public String getEncodedCoordinates(boolean keepPositionAndScale, Coordinates[] atomCoordinates) {
		if (mEncodedCoords == null) {
			generateGraph();
			encodeCoordinates(keepPositionAndScale, atomCoordinates);
			}

		return mEncodedCoords;
		}

	private void encodeCoordinates(boolean keepPositionAndScale, Coordinates[] coords) {
		if (mMol.getAtoms() == 0) {
			mEncodedCoords = "";
			return;
			}

		// if we have 3D-coords and explicit hydrogens and if all hydrogens are explicit then encode hydrogen coordinates
		boolean includeHydrogenCoordinates = false;
		if (mZCoordinatesAvailable
		 && mMol.getAllAtoms() > mMol.getAtoms()
		 && !mMol.isFragment()) {
			includeHydrogenCoordinates = true;
			for (int i=0; i<mMol.getAtoms(); i++) {
				if (mMol.getImplicitHydrogens(i) != 0) {
					includeHydrogenCoordinates = false;
					break;
					}
				}
			}

		int resolutionBits = mZCoordinatesAvailable ? 16 : 8;	// must be an even number
		encodeBitsStart(true);
		mEncodingBuffer.append(includeHydrogenCoordinates ? '#' : '!');
		encodeBits(mZCoordinatesAvailable ? 1 : 0, 1);
		encodeBits(keepPositionAndScale ? 1 : 0, 1);
		encodeBits(resolutionBits/2, 4);	// resolution bits devided by 2

		double maxDelta = 0.0;
		for (int i=1; i<mMol.getAtoms(); i++)
			maxDelta = getMaxDelta(mGraphAtom[i], (mGraphFrom[i] == -1) ? -1 : mGraphAtom[mGraphFrom[i]], maxDelta, coords);
		if (includeHydrogenCoordinates) {
			for (int i=0; i<mMol.getAtoms(); i++) {
				int atom = mGraphAtom[i];
				for (int j=mMol.getConnAtoms(atom); j<mMol.getAllConnAtoms(atom); j++)
					maxDelta = getMaxDelta(mMol.getConnAtom(atom, j), atom, maxDelta, coords);
				}
			}

		if (mMol.getAtoms() > 1 && maxDelta == 0.0) {
			mEncodedCoords = "";
			return;
			}

		int binCount = (1 << resolutionBits);
		double increment = maxDelta / (binCount / 2.0 - 1);
		double maxDeltaPlusHalfIncrement = maxDelta + increment / 2.0;

		for (int i=1; i<mMol.getAtoms(); i++)
			encodeCoords(mGraphAtom[i], (mGraphFrom[i] == -1) ? -1 : mGraphAtom[mGraphFrom[i]], maxDeltaPlusHalfIncrement, increment, resolutionBits, coords);
		if (includeHydrogenCoordinates) {
			for (int i=0; i<mMol.getAtoms(); i++) {
				int atom = mGraphAtom[i];
				for (int j=mMol.getConnAtoms(atom); j<mMol.getAllConnAtoms(atom); j++)
					encodeCoords(mMol.getConnAtom(atom, j), atom, maxDeltaPlusHalfIncrement, increment, resolutionBits, coords);
				}
			}

		if (keepPositionAndScale) {
			double avblDefault = mZCoordinatesAvailable ? 1.5 : Molecule.getDefaultAverageBondLength();
			double avbl = mMol.getAverageBondLength(
					includeHydrogenCoordinates ? mMol.getAllAtoms() : mMol.getAtoms(),
					includeHydrogenCoordinates ? mMol.getAllBonds() : mMol.getBonds(),
					avblDefault, coords);
			encodeBits(encodeABVL(avbl, binCount), resolutionBits);

			encodeBits(encodeShift(coords[mGraphAtom[0]].x / avbl, binCount), resolutionBits);
			encodeBits(encodeShift(coords[mGraphAtom[0]].y / avbl, binCount), resolutionBits);

			if (mZCoordinatesAvailable)
				encodeBits(encodeShift(coords[mGraphAtom[0]].z / avbl, binCount), resolutionBits);
			}

		mEncodedCoords = encodeBitsEnd();
		}

	private double getMaxDelta(int atom, int from, double maxDelta, Coordinates[] coords) {
		double deltaX = (from == -1) ?
						Math.abs(coords[atom].x - coords[mGraphAtom[0]].x) / 8.0
					  : Math.abs(coords[atom].x - coords[from].x);
		if (maxDelta < deltaX)
			maxDelta = deltaX;

		double deltaY = (from == -1) ?
						Math.abs(coords[atom].y - coords[mGraphAtom[0]].y) / 8.0
					  : Math.abs(coords[atom].y - coords[from].y);
		if (maxDelta < deltaY)
			maxDelta = deltaY;

		if (mZCoordinatesAvailable) {
			double deltaZ = (from == -1) ?
							Math.abs(coords[atom].z - coords[mGraphAtom[0]].z) / 8.0
						  : Math.abs(coords[atom].z - coords[from].z);
			if (maxDelta < deltaZ)
				maxDelta = deltaZ;
			}

		return maxDelta;
		}

	private void encodeCoords(int atom, int from, double maxDeltaPlusHalfIncrement, double increment, int resolutionBits, Coordinates[] coords) {
		double deltaX = (from == -1) ?
						(coords[atom].x - coords[mGraphAtom[0]].x) / 8.0
					   : coords[atom].x - coords[from].x;

		double deltaY = (from == -1) ?
						(coords[atom].y - coords[mGraphAtom[0]].y) / 8.0
					   : coords[atom].y - coords[from].y;

		encodeBits((int)((maxDeltaPlusHalfIncrement + deltaX) / increment), resolutionBits);
		encodeBits((int)((maxDeltaPlusHalfIncrement + deltaY) / increment), resolutionBits);

		if (mZCoordinatesAvailable) {
			double deltaZ = (from == -1) ?
							(coords[atom].z - coords[mGraphAtom[0]].z) / 8.0
						   : coords[atom].z - coords[from].z;

			encodeBits((int)((maxDeltaPlusHalfIncrement + deltaZ) / increment), resolutionBits);
			}
		}

	/**
	 * Encode a floating point value into an integer with precision proportional to the value itself.
	 * @param value
	 * @return
	 */
	private int encodeABVL(double value, int binCount) {
		return Math.min(binCount-1, Math.max(0, (int)(0.5 + Math.log10(value/0.1) / Math.log10(200/0.1) * (binCount-1))));
		}

	private int encodeShift(double value, int binCount) {
		int halfBinCount = binCount / 2;
		boolean isNegative =  (value < 0);
		value = Math.abs(value);
		double steepness = binCount/32;
		int intValue = Math.min(halfBinCount-1, (int)Math.round(value * halfBinCount / (value + steepness)));
		return isNegative ? halfBinCount + intValue : intValue;
		}

	public String getEncodedMapping() {
		if (mMapping == null) {
			generateGraph();
			encodeMapping();
			}

		return mMapping;
		}


	private void encodeMapping() {
		if (mMol.getAtoms() == 0) {
			mMapping = "";
			return;
			}

		int maxMapNo = 0;
		boolean autoMappingFound = false;
		boolean manualMappingFound = false;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (maxMapNo < mMol.getAtomMapNo(atom))
				maxMapNo = mMol.getAtomMapNo(atom);
			if (mMol.isAutoMappedAtom(atom))
				autoMappingFound = true;
			else
				manualMappingFound = true;
			}

		if (maxMapNo == 0) {
			mMapping = "";
			return;
			}

		int nbits = getNeededBits(maxMapNo);
		encodeBitsStart(true);
		encodeBits(nbits, 4);
		encodeBits(autoMappingFound ? 1 : 0, 1);
		encodeBits(manualMappingFound ? 1 : 0, 1);
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			encodeBits(mMol.getAtomMapNo(mGraphAtom[atom]), nbits);
			if (autoMappingFound && manualMappingFound)
				encodeBits(mMol.isAutoMappedAtom(mGraphAtom[atom]) ? 1 : 0, 1);
			}

		mMapping = encodeBitsEnd();
		}


	/**
	 * Creates parities based on atom indices in graph rather than on priority values.
	 * These values are more meaningful to be written into idcodes, because they allow
	 * to create coordinates or running Configuration aware substructure searches on
	 * molecules creates from idcode without the necessity to recreate the priority values.
	 */
	private void idGenerateConfigurations() {
		mTHConfiguration = new byte[mMol.getAtoms()];

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mTHParity[atom] == Molecule.cAtomParity1
			 || mTHParity[atom] == Molecule.cAtomParity2) {
				boolean inversion = mTHParityIsMesoInverted[atom];
				if (mMol.isCentralAlleneAtom(atom)) {
					for (int i=0; i<mMol.getConnAtoms(atom); i++) {
						int connAtom = mMol.getConnAtom(atom,i);
						int neighbours = 0;
						int[] neighbour = new int[3];
						for (int j=0; j<mMol.getConnAtoms(connAtom); j++) {
							neighbour[neighbours] = mMol.getConnAtom(connAtom,j);
							if (neighbour[neighbours] != atom)
								neighbours++;
							}
						if (neighbours == 2
						 && ((mCanRank[neighbour[0]] > mCanRank[neighbour[1]])
							^(mGraphIndex[neighbour[0]] < mGraphIndex[neighbour[1]])))
							inversion = !inversion;
						}
					}
				else {
					for (int i=1; i<mMol.getConnAtoms(atom); i++) {
						for (int j=0; j<i; j++) {
							int connAtom1 = mMol.getConnAtom(atom,i);
							int connAtom2 = mMol.getConnAtom(atom,j);
							if (mCanRank[connAtom1] > mCanRank[connAtom2])
								inversion = !inversion;
							if (mGraphIndex[connAtom1] < mGraphIndex[connAtom2])
								inversion = !inversion;
							}
						}
					}

				mTHConfiguration[atom] = ((mTHParity[atom] == Molecule.cAtomParity1) ^ inversion) ?
											(byte)Molecule.cAtomParity1 : Molecule.cAtomParity2;
				}
			else {
				mTHConfiguration[atom] = mTHParity[atom];
				}
			}

		mEZConfiguration = new byte[mMol.getBonds()];

		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mEZParity[bond] == Molecule.cBondParityEor1
			 || mEZParity[bond] == Molecule.cBondParityZor2) {
				boolean inversion = false;
				for (int i=0; i<2; i++) {
					int bondAtom = mMol.getBondAtom(i,bond);
					if (mMol.getConnAtoms(bondAtom) == 3) {
						int[] neighbour = new int[2];
						int neighbours = 0;
						for (int j=0; j<3; j++)
							if (mMol.getConnAtom(bondAtom,j) != mMol.getBondAtom(1-i,bond))
								neighbour[neighbours++] = mMol.getConnAtom(bondAtom,j);
						if (mCanRank[neighbour[0]] > mCanRank[neighbour[1]])
							inversion = !inversion;
						if (mGraphIndex[neighbour[0]] < mGraphIndex[neighbour[1]])
							inversion = !inversion;
						}
					}

				mEZConfiguration[bond] = ((mEZParity[bond] == Molecule.cBondParityEor1) ^ inversion) ?
										   (byte)Molecule.cBondParityEor1 : Molecule.cBondParityZor2;
				}
			else {
				mEZConfiguration[bond] = mEZParity[bond];
				}
			}
		}


	private void idNormalizeESRGroupNumbers() {
		idNormalizeESRGroupNumbers(Molecule.cESRTypeAnd);
		idNormalizeESRGroupNumbers(Molecule.cESRTypeOr);
		}

	
	private void idNormalizeESRGroupNumbers(int type) {
		int[] groupRank = new int[Molecule.cESRMaxGroups];
		int groups = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if ((mTHConfiguration[atom] == Molecule.cAtomParity1
			  || mTHConfiguration[atom] == Molecule.cAtomParity2)
			 && mTHESRType[atom] == type) {
				int group = mTHESRGroup[atom];
				if (groupRank[group] < mCanRank[atom]) {
					if (groupRank[group] == 0)
						groups++;
					groupRank[group] = mCanRank[atom];
					}
				}
			}
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if ((mEZConfiguration[bond] == Molecule.cBondParityEor1
			  || mEZConfiguration[bond] == Molecule.cBondParityZor2)
			 && mEZESRType[bond] == type
			 && mMol.getBondType(bond) == Molecule.cBondTypeSingle) {
				int group = mEZESRGroup[bond];
				int rank = Math.max(mCanRank[mMol.getBondAtom(0, bond)], mCanRank[mMol.getBondAtom(1, bond)]);
				if (groupRank[group] < rank) {
					if (groupRank[group] == 0)
						groups++;
					groupRank[group] = rank;
					}
				}
			}
		byte[] canGroup = new byte[Molecule.cESRMaxGroups];
		for (int i=0; i<groups; i++) {
			int maxGroup = -1;
			int maxRank = 0;
			for (int j=0; j<Molecule.cESRMaxGroups; j++) {
				if (maxRank < groupRank[j]) {
					maxRank = groupRank[j];
					maxGroup = j;
					}
				}
			groupRank[maxGroup] = 0;
			canGroup[maxGroup] = (byte)i;
			}
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if ((mTHConfiguration[atom] == Molecule.cAtomParity1
			  || mTHConfiguration[atom] == Molecule.cAtomParity2)
			 && mTHESRType[atom] == type)
				mTHESRGroup[atom] = canGroup[mTHESRGroup[atom]];
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if ((mEZConfiguration[bond] == Molecule.cBondParityEor1
			  || mEZConfiguration[bond] == Molecule.cBondParityZor2)
			 && mEZESRType[bond] == type
			 && mMol.getBondType(bond) == Molecule.cBondTypeSingle)
				mEZESRGroup[bond] = canGroup[mEZESRGroup[bond]];
		}


	private void encodeBitsStart(boolean avoid127) {
		mEncodingBuffer = new StringBuilder();
		mEncodingBitsAvail = 6;
		mEncodingTempData = 0;
		mEncodeAvoid127 = avoid127;
		}


	private void encodeFeatureNo(int codeNo) {
		for (int i=0; i<mFeatureBlock; i++)
			codeNo -= 16;

		if (codeNo < 0)
			System.out.println("ERROR in Canonizer: Code unexpectedly low.");

		while (codeNo > 15) {
			encodeBits(1, 1);   //  more data to come
			encodeBits(15, 4);  //  15 = datatype 'start next feature set'
			codeNo -= 16;
			mFeatureBlock++;
			}

		encodeBits(1, 1);   // more features to come
		encodeBits(codeNo, 4);
		}

	private void encodeBits(long data, int bits) {
//System.out.println(bits+" bits:"+data+"  mode="+mode);
		while (bits != 0) {
			if (mEncodingBitsAvail == 0) {
				if (!mEncodeAvoid127 || mEncodingTempData != 63)
					mEncodingTempData += 64;
				mEncodingBuffer.append((char)mEncodingTempData);
				mEncodingBitsAvail = 6;
				mEncodingTempData = 0;
				}
			mEncodingTempData <<= 1;
			mEncodingTempData |= (data & 1);
			data >>= 1;
			bits--;
			mEncodingBitsAvail--;
			}
		}


	private String encodeBitsEnd() {
		mEncodingTempData <<= mEncodingBitsAvail;
		if (!mEncodeAvoid127 || mEncodingTempData != 63)
			mEncodingTempData += 64;
		mEncodingBuffer.append((char)mEncodingTempData);
		return mEncodingBuffer.toString();
		}


	/**
	 * @param maxNo highest possible index of some kind
	 * @return number of bits needed to represent numbers up to maxNo
	 */
	public static int getNeededBits(int maxNo) {
		int bits = 0;
		while (maxNo > 0) {
			maxNo >>= 1;
			bits++;
			}
		return bits;
		}


	/**
	 * Returns the absolute tetrahedral parity, which is based on priority ranks.
	 * @param atom
	 * @return one of the Molecule.cAtomParityXXX constants
	 */
	public int getTHParity(int atom) {
		return mTHParity[atom];
		}


	/**
	 * Returns the absolute bond parity, which is based on priority ranks.
	 * @param bond
	 * @return one of the Molecule.cBondParityXXX constants
	 */
	public int getEZParity(int bond) {
		return mEZParity[bond];
		}


	/**
	 * If mMode includes CREATE_PSEUDO_STEREO_GROUPS, then this method returns
	 * the number of independent relative stereo feature groups. A relative stereo
	 * feature group always contains more than one pseudo stereo features (TH or EZ),
	 * which only in combination define a certain stereo configuration.
	 * @return
	 */
	public int getPseudoStereoGroupCount() {
		return mNoOfPseudoGroups;
		}

	/**
	 * If mMode includes CREATE_PSEUDO_STEREO_GROUPS, then this method returns
	 * this bond's relative stereo feature group number provided this bond is a
	 * pseudo stereo bond, i.e. its stereo configuration only is relevant in
	 * combination with other pseudo stereo features.
	 * If this bond is not a pseudo stereo bond, then this method returns 0.
	 * @param bond
	 * @return
	 */
	public int getPseudoEZGroup(int bond) {
		return mPseudoEZGroup[bond];
		}


	/**
	 * If mMode includes CREATE_PSEUDO_STEREO_GROUPS, then this method returns
	 * this atom's relative stereo feature group number provided this atom is a
	 * pseudo stereo center, i.e. its stereo configuration only is relevant in
	 * combination with other pseudo stereo features.
	 * If this atom is not a pseudo stereo center, then this method returns 0.
	 * @param atom
	 * @return
	 */
	public int getPseudoTHGroup(int atom) {
		return mPseudoTHGroup[atom];
		}


	/**
	 * This normalizes all absolute tetrahedral-, allene- and atrop-parities within the molecule.
	 * This is done by finding the lowest atom rank that is shared by an odd number of
	 * atoms with determines parities, not counting unknown and none.
	 * If there number of parity2 atoms is higher than parity1 atoms of that rank, then
	 * all parities are inverted.<br>
	 * You may call this method before creating the idcode from this Canonizer to convert
	 * internal parity information to the noermalized enantiomer. When calling getIDCode()
	 * afterwards, the idcode represents the normalized enantiomer. Stereo information of
	 * the underlying molecule is not touched.
	 * @return true, if all internal parities were inverted
	 */
	public boolean normalizeEnantiomer() {
		int[] parityCount = new int[mNoOfRanks + 1];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.getAtomESRType(atom) == Molecule.cESRTypeAbs) {
				if (mTHParity[atom] == Molecule.cAtomParity1)
					parityCount[mCanRank[atom]]++;
				else if (mTHParity[atom] == Molecule.cAtomParity2)
					parityCount[mCanRank[atom]]--;
				}
			}
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondOrder(bond) == 1 && mMol.getBondESRType(bond) == Molecule.cESRTypeAbs) {
				if (mEZParity[bond] == Molecule.cBondParityEor1) {
					parityCount[mCanRank[mMol.getBondAtom(0, bond)]]++;
					parityCount[mCanRank[mMol.getBondAtom(1, bond)]]++;
					}
				else if (mEZParity[bond] == Molecule.cBondParityZor2) {
					parityCount[mCanRank[mMol.getBondAtom(0, bond)]]--;
					parityCount[mCanRank[mMol.getBondAtom(1, bond)]]--;
					}
				}
			}
		for (int rank=1; rank<=mNoOfRanks; rank++) {
			if (parityCount[rank] != 0) {
				boolean invert = (parityCount[rank] < 0);
				if (invert) {
					for (int atom=0; atom<mMol.getAtoms(); atom++) {
						if (mMol.getAtomESRType(atom) == Molecule.cESRTypeAbs) {
							if (mTHParity[atom] == Molecule.cAtomParity1)
								mTHParity[atom] = Molecule.cAtomParity2;
							else if (mTHParity[atom] == Molecule.cAtomParity2)
								mTHParity[atom] = Molecule.cAtomParity1;
							}
						}
					for (int bond=0; bond<mMol.getBonds(); bond++) {
						if (mMol.getBondOrder(bond) == 1 && mMol.getBondESRType(bond) == Molecule.cESRTypeAbs) {
							if (mEZParity[bond] == Molecule.cBondParityEor1)
								mEZParity[bond] = Molecule.cBondParityZor2;
							else if (mEZParity[bond] == Molecule.cBondParityZor2)
								mEZParity[bond] = Molecule.cBondParityEor1;
							}
						}
					}
				return invert;
				}
			}

		return false;
		}

	/**
	 * Creates parities based on atom indices of original molecule and
	 * copies them back into that molecule. It also sets the stereo center flag.
	 * These atom parities are not normalized (e.g. in racemic or meso groups)
	 * and reflect the original molecule's up/down bonds.
	 */
	public void setParities() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mTHParity[atom] == Molecule.cAtomParity1
			 || mTHParity[atom] == Molecule.cAtomParity2) {
				boolean inversion = false;
				if (mMol.isCentralAlleneAtom(atom)) {   // allene parities
					for (int i=0; i<mMol.getConnAtoms(atom); i++) {
						int connAtom = mMol.getConnAtom(atom,i);
						int neighbours = 0;
						int[] neighbour = new int[3];
						for (int j=0; j<mMol.getConnAtoms(connAtom); j++) {
							neighbour[neighbours] = mMol.getConnAtom(connAtom,j);
							if (neighbour[neighbours] != atom)
								neighbours++;
							}
						if (neighbours == 2
						 && ((mCanRank[neighbour[0]] > mCanRank[neighbour[1]])
							^(neighbour[0] < neighbour[1])))
							inversion = !inversion;
						}
					}
				else {
					for (int i=1; i<mMol.getConnAtoms(atom); i++) {
						for (int j=0; j<i; j++) {
							int connAtom1 = mMol.getConnAtom(atom,i);
							int connAtom2 = mMol.getConnAtom(atom,j);
							if (mCanRank[connAtom1] > mCanRank[connAtom2])
								inversion = !inversion;
							if (connAtom1 < connAtom2)
								inversion = !inversion;
							}
						}
					}

				mMol.setAtomParity(atom, ((mTHParity[atom] == Molecule.cAtomParity1) ^ inversion) ?
										   Molecule.cAtomParity1 : Molecule.cAtomParity2,
										   mTHParityIsPseudo[atom]);
				}
			else {
				mMol.setAtomParity(atom, mTHParity[atom], mTHParityIsPseudo[atom]);
				}
			}

		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mEZParity[bond] == Molecule.cBondParityEor1
			 || mEZParity[bond] == Molecule.cBondParityZor2) {
				boolean inversion = false;
				for (int i=0; i<2; i++) {
					int bondAtom = mMol.getBondAtom(i,bond);
					if (mMol.getConnAtoms(bondAtom) == 3) {
						int[] neighbour = new int[2];
						int neighbours = 0;
						for (int j=0; j<3; j++)
							if (mMol.getConnAtom(bondAtom,j) != mMol.getBondAtom(1-i,bond))
								neighbour[neighbours++] = mMol.getConnAtom(bondAtom,j);
						if (mCanRank[neighbour[0]] > mCanRank[neighbour[1]])
							inversion = !inversion;
						if (neighbour[0] < neighbour[1])
							inversion = !inversion;
						}
					}

				mMol.setBondParity(bond, ((mEZParity[bond] == Molecule.cBondParityEor1) ^ inversion) ?
										   Molecule.cBondParityEor1 : Molecule.cBondParityZor2,
										   mEZParityIsPseudo[bond]);
				}
			else {
				mMol.setBondParity(bond, mEZParity[bond], mEZParityIsPseudo[bond]);
				}
			}
		}


	protected void setStereoCenters() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			mMol.setAtomStereoCenter(atom, mIsStereoCenter[atom]);
			}
		}


	protected void setCIPParities() {
		if (mTHCIPParity != null)
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				mMol.setAtomCIPParity(atom, mTHCIPParity[atom]);
		if (mEZCIPParity != null)
			for (int bond=0; bond<mMol.getBonds(); bond++)
				mMol.setBondCIPParity(bond, mEZCIPParity[bond]);
		}


	private void cipCalcTHParity(int atom) {
		if ((mTHParity[atom] == Molecule.cAtomParity1
		  || mTHParity[atom] == Molecule.cAtomParity2)) {
			boolean invertedOrder = false;

			if (mMol.getAtomPi(atom) == 2) {	// allene
				try {
					for (int i=0; i<2; i++) {
						int alleneAtom = mMol.getConnAtom(atom,i);
						if (mMol.getConnAtoms(alleneAtom) == 3) {
							int[] connAtom = new int[2];
							int count = 0;
							for (int j=0; j<mMol.getConnAtoms(alleneAtom); j++)
								if (mMol.getConnBondOrder(alleneAtom,j) == 1)
									connAtom[count++] = mMol.getConnAtom(alleneAtom,j);
							if ((mCanRank[connAtom[0]] > mCanRank[connAtom[1]])
							  ^ cipComparePriority(alleneAtom,connAtom[0],connAtom[1]))
								invertedOrder = !invertedOrder;
							}
						}
					}
				catch (Exception e) {
					mTHCIPParity[atom] = Molecule.cAtomCIPParityProblem; // to indicate assignment problem
					return;
					}
				}
			else {
				int[] cipConnAtom;
				try {
					cipConnAtom = cipGetOrderedConns(atom);
					}
				catch (Exception e) {
					mTHCIPParity[atom] = Molecule.cAtomCIPParityProblem; // to indicate assignment problem
					return;
					}
				for (int i=1; i<cipConnAtom.length; i++)
					for (int j=0; j<i; j++)
						if (mCanRank[cipConnAtom[i]] < mCanRank[cipConnAtom[j]])
							invertedOrder = !invertedOrder;
				}

			if ((mTHParity[atom] == Molecule.cAtomParity1) ^ invertedOrder)
				mTHCIPParity[atom] = Molecule.cAtomCIPParityRorM;
			else
				mTHCIPParity[atom] = Molecule.cAtomCIPParitySorP;
			}
		}


	private void cipCalcEZParity(int bond) {
		if ((mEZParity[bond] == Molecule.cBondParityEor1
		  || mEZParity[bond] == Molecule.cBondParityZor2)
		 && !mMol.isSmallRingBond(bond)) {
			boolean invertedOrder = false;

		   	try {
				for (int i=0; i<2; i++) {
					int bondAtom = mMol.getBondAtom(i,bond);
					if (mMol.getConnAtoms(bondAtom) == 3) {
						int[] connAtom = new int[2];
						int count = 0;
						for (int j=0; j<mMol.getConnAtoms(bondAtom); j++)
							if (mMol.getConnBond(bondAtom,j) != bond)
								connAtom[count++] = mMol.getConnAtom(bondAtom,j);
						if ((mCanRank[connAtom[0]] > mCanRank[connAtom[1]])
						  ^ cipComparePriority(bondAtom,connAtom[0],connAtom[1]))
							invertedOrder = !invertedOrder;
						}
					}
				}
			catch (Exception e) {
				mEZCIPParity[bond] = Molecule.cBondCIPParityProblem; // to indicate assignment problem
				return;
				}

			if ((mEZParity[bond] == Molecule.cBondParityEor1) ^ invertedOrder)
				mEZCIPParity[bond] = Molecule.cBondCIPParityEorP;
			else
				mEZCIPParity[bond] = Molecule.cBondCIPParityZorM;
			}
		}


	private int[] cipGetOrderedConns(int atom) throws Exception {
		int noOfConns = mMol.getAllConnAtoms(atom);
		int[] orderedConn = new int[noOfConns];
		for (int i=0; i<noOfConns; i++)
			orderedConn[i] = mMol.getConnAtom(atom,i);
		for (int i=noOfConns; i>1; i--) {
			boolean found = false;
			for (int j=1; j<i; j++) {
				if (cipComparePriority(atom, orderedConn[j-1], orderedConn[j])) {
					found = true;
					int temp = orderedConn[j-1];
					orderedConn[j-1] = orderedConn[j];
					orderedConn[j] = temp;
					}
				}
			if (!found)
				break;
			}
		return orderedConn;
		}

	/**
	 * Returns true if atom1's CIP-priority is higher than atom2's CIP-priority.
	 * 10.03.2009 MvK added emergency break in while loop.
	 * @param rootAtom
	 * @param atom1
	 * @param atom2
	 * @return
	 * @throws Exception
	 */
	private boolean cipComparePriority(int rootAtom, int atom1, int atom2) throws Exception {
		if (mMol.getAtomicNo(atom1) != mMol.getAtomicNo(atom2))
			return (mMol.getAtomicNo(atom1) > mMol.getAtomicNo(atom2));
		if (mMol.getAtomMass(atom1) != mMol.getAtomMass(atom2)) {
			int mass1 = mMol.isNaturalAbundance(atom1) ? Molecule.cRoundedMass[mMol.getAtomicNo(atom1)] : mMol.getAtomMass(atom1);
			int mass2 = mMol.isNaturalAbundance(atom2) ? Molecule.cRoundedMass[mMol.getAtomicNo(atom2)] : mMol.getAtomMass(atom2);
			return (mass1 > mass2);
			}

		int graphSize = mMol.getAtoms();

		int[] graphAtom = new int[graphSize];
		int[] graphParent = new int[graphSize];
		int[] graphRank = new int[graphSize];
		boolean[] graphIsPseudo = new boolean[graphSize];

		boolean[] atomUsed = new boolean[mMol.getAllAtoms()];

		graphAtom[0] = rootAtom;
		graphAtom[1] = atom1;
		graphAtom[2] = atom2;

		graphParent[0] = -1;
		graphParent[1] = 0;
		graphParent[2] = 0;

		atomUsed[rootAtom] = true;
		atomUsed[atom1] = true;
		atomUsed[atom2] = true;

		int current = 1;
		int highest = 2;

		int[] levelStart = new int[64];
		levelStart[1] = 1;
		levelStart[2] = 3;
		int currentLevel = 2;
		while (current <= highest) {
			while (current < levelStart[currentLevel]) {
				int currentAtom = graphAtom[current];

								// do not consider neighbours of pseudo atoms
				if (!graphIsPseudo[current]) {
					int delocalizedBondCount = 0;
					int delocalizedMeanAtomicNo = 0;

					for (int i=0; i<mMol.getConnAtoms(currentAtom); i++) {
						int candidate = mMol.getConnAtom(currentAtom, i);

						if (highest+mMol.getConnBondOrder(currentAtom, i)+1 >= graphSize) {
							graphSize += mMol.getAtoms();
							graphAtom = resize(graphAtom, graphSize);
							graphParent = resize(graphParent, graphSize);
							graphRank = resize(graphRank, graphSize);
							graphIsPseudo = resize(graphIsPseudo, graphSize);
							}

						if (mMol.isDelocalizedBond(mMol.getConnBond(currentAtom, i))) {
							// if candidate is part of a delocalized ring, than we need to add
							// one pseudo atom with mean atomic no between all delocalized neighbors
							delocalizedBondCount++;
							delocalizedMeanAtomicNo += mMol.getAtomicNo(candidate);
							}
						else {
							// add pseudo atoms for double and triple bonds
							for (int j=1; j<mMol.getConnBondOrder(currentAtom, i); j++) {
								highest++;
								graphAtom[highest] = candidate;
								graphParent[highest] = current;
								graphIsPseudo[highest] = true;
								}
							}

						int parentGraphIndex = graphParent[current];
						if (candidate == graphAtom[parentGraphIndex])
							continue;

						boolean atomInParentChain = false;
						if (atomUsed[candidate]) {
							int parent = graphParent[parentGraphIndex];
							while (parent != -1) {
								if (candidate == graphAtom[parent]) {
									atomInParentChain = true;
									break;
									}
								parent = graphParent[parent];
								}
							}

						// add all atoms moving away from rootAtom
						if (atomInParentChain) {
							// if we have a cycle closure then add closure atom
							// a second time as pseudo atom
							highest++;
							graphAtom[highest] = candidate;
							graphParent[highest] = current;
							graphIsPseudo[highest] = true;
							}
						else {
							highest++;
							graphAtom[highest] = candidate;
							graphParent[highest] = current;
							atomUsed[candidate] = true;
							}
						}

					if (delocalizedBondCount != 0) {
						highest++;
						graphRank[highest] = (delocalizedMeanAtomicNo << 2) / delocalizedBondCount;
						graphParent[highest] = current;
						graphIsPseudo[highest] = true;
						}
					}
				current++;

				if(current == 10000) {
/* TODO fix this
					try {
						java.io.BufferedWriter w = new java.io.BufferedWriter(new java.io.FileWriter("/Users/sandert/canonizerEmergencyRoot"+rootAtom+".txt"));
						w.write("index\tgraphAtom\tgraphParent\tgraphRank\tgraphIsPseudo\tatomUsed");
						w.newLine();
						for (int i=0; i<10000; i++) {
							w.write(""+i+"\t"+graphAtom[i]+"\t"+graphParent[i]+"\t"+graphRank[i]+"\t"+graphIsPseudo[i]+"\t"+atomUsed[graphAtom[i]]);
							w.newLine();
							}
						w.close();
						}
					catch (Exception e) {}
*/
					throw new Exception("Emergency break in while loop.");
					}
				}

			if (levelStart.length == currentLevel+1)
				levelStart = resize(levelStart, levelStart.length+64);

			levelStart[currentLevel+1] = highest + 1;

			// compile initial ranks of current level
			for (int i=levelStart[currentLevel]; i<levelStart[currentLevel+1]; i++) {
				if (graphRank[i] == 0)
					graphRank[i] = (mMol.getAtomicNo(graphAtom[i]) == 151 ? 1   // D
								  : mMol.getAtomicNo(graphAtom[i]) == 152 ? 1   // T
								  : mMol.getAtomicNo(graphAtom[i])) << 2;
				graphRank[i] += (graphRank[graphParent[i]] << 16);
				}

			// Locate parent atoms with equal ranks that may be distinguished
			// by attachments on current level. Adjust parent ranks accordingly
			// and propagate rank distinctions upwards the hierarchy.
			cipUpdateParentRanking(graphIsPseudo, graphRank, graphParent, graphAtom, levelStart, currentLevel);

/*
System.out.println("Ranking atomic numbers on currentLevel:"+currentLevel+" levelStart:"+levelStart[currentLevel]+" nextLevelStart:"+levelStart[currentLevel+1]);
System.out.print(" graphAtoms:");
for (int i=0; i<levelStart[currentLevel]; i++) System.out.print(graphAtom[i]+" ");
System.out.print(" >> ");
for (int i=levelStart[currentLevel]; i<levelStart[currentLevel+1]; i++) System.out.print(graphAtom[i]+" ");
System.out.println("");
System.out.print(" graphRank1:");
for (int i=0; i<levelStart[currentLevel]; i++) System.out.print(Integer.toHexString(graphRank[i])+" ");
System.out.print(" >> ");
for (int i=levelStart[currentLevel]; i<levelStart[currentLevel+1]; i++) System.out.print(Integer.toHexString(graphRank[i])+" ");
System.out.println("");
if (graphRank[1] != graphRank[2]) System.out.println("@@@ atom"+atom1+((graphRank[1]>graphRank[2])?" > ":" < ")+"atom"+atom2);
*/
			// check whether both substituents are now recognized to be different
			if (graphRank[1] != graphRank[2])
				return (graphRank[1] > graphRank[2]);

			// generate relative ranking for current level by adding all relative
			// ranks of parent levels with higher significance the higher a parent
			// is in the tree.
			if (currentLevel > 1)
//{
				cipCompileRelativeRanks(graphRank, graphParent, levelStart, currentLevel);
/*
System.out.print(" graphRank2:");
for (int i=0; i<levelStart[currentLevel]; i++) System.out.print(Integer.toHexString(graphRank[i])+" ");
System.out.print(" >> ");
for (int i=levelStart[currentLevel]; i<levelStart[currentLevel+1]; i++) System.out.print(Integer.toHexString(graphRank[i])+" ");
System.out.println("");
}
*/
			currentLevel++;
			}

		// consider isotop differences
		int[] cipRank = new int[mMol.getAtoms()];
		boolean isotopDataFound = false;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (atomUsed[atom] && !mMol.isNaturalAbundance(atom)) {
				isotopDataFound = true;
				break;
				}
			}
		if (isotopDataFound) {
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				cipRank[atom] = mMol.isNaturalAbundance(atom) ? Molecule.cRoundedMass[mMol.getAtomicNo(atom)] : mMol.getAtomMass(atom);

			if (cipTryDistinguishBranches(graphIsPseudo, graphRank, graphParent, graphAtom, cipRank, levelStart, currentLevel))
				return (graphRank[1] > graphRank[2]);
			}

		// consider E/Z configuration differences
		Arrays.fill(cipRank, 0);
		boolean ezDataFound = false;
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (atomUsed[mMol.getBondAtom(0, bond)]
			 || atomUsed[mMol.getBondAtom(1, bond)]) {
				if (mEZCIPParity[bond] == Molecule.cBondCIPParityEorP) {
					cipRank[mMol.getBondAtom(0,bond)] = 1;
					cipRank[mMol.getBondAtom(1,bond)] = 1;
					ezDataFound = true;
					}
				else if (mEZCIPParity[bond] == Molecule.cBondCIPParityZorM) {
					cipRank[mMol.getBondAtom(0,bond)] = 2;
					cipRank[mMol.getBondAtom(1,bond)] = 2;
					ezDataFound = true;
					}
				}
			}
		if (ezDataFound
		 && cipTryDistinguishBranches(graphIsPseudo, graphRank, graphParent, graphAtom, cipRank, levelStart, currentLevel))
			return (graphRank[1] > graphRank[2]);

		// TODO consider relative configuration differences RR/SS higher than RS/SR etc.

		// consider R/S configuration differences
		Arrays.fill(cipRank, 0);
		boolean rsDataFound = false;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (atomUsed[atom]) {
				if (mTHCIPParity[atom] == Molecule.cAtomCIPParitySorP) {
					cipRank[atom] = 1;
					rsDataFound = true;
					}
				else if (mTHCIPParity[atom] == Molecule.cAtomCIPParityRorM) {
					cipRank[atom] = 2;
					rsDataFound = true;
					}
				}
			}
		if (rsDataFound
		 && cipTryDistinguishBranches(graphIsPseudo, graphRank, graphParent,
									  graphAtom, cipRank, levelStart, currentLevel))
			return (graphRank[1] > graphRank[2]);

		mCIPParityNoDistinctionProblem = true;
		throw new Exception("no distinction applying CIP rules");
		}

	private boolean cipTryDistinguishBranches(boolean[] graphIsPseudo,
											  int[] graphRank,
											  int[] graphParent,
											  int[] graphAtom,
											  int[] cipRank,
											  int[] levelStart,
											  int currentLevel) {
		for (int level=1; level<currentLevel; level++) {
			for (int i=levelStart[level]; i<levelStart[level+1]; i++)
				graphRank[i] = cipRank[graphAtom[i]]
							 + (graphRank[graphParent[i]] << 8);

			cipUpdateParentRanking(graphIsPseudo, graphRank, graphParent,
								   graphAtom, levelStart, level);
/*
System.out.println("2nd ranking atomic numbers on level:"+level+" levelStart:"+levelStart[level]+" nextLevelStart:"+levelStart[level+1]);
System.out.print(" graphAtoms:");
for (int i=0; i<levelStart[level]; i++) System.out.print(graphAtom[i]+" ");
System.out.print(" >> ");
for (int i=levelStart[level]; i<levelStart[level+1]; i++) System.out.print(graphAtom[i]+" ");
System.out.println("");
System.out.print(" graphRank1:");
for (int i=0; i<levelStart[level]; i++) System.out.print(Integer.toHexString(graphRank[i])+" ");
System.out.print(" >> ");
for (int i=levelStart[level]; i<levelStart[level+1]; i++) System.out.print(Integer.toHexString(graphRank[i])+" ");
System.out.println("");
if (graphRank[1] != graphRank[2]) System.out.println("@@@ atom"+graphAtom[1]+((graphRank[1]>graphRank[2])?" > ":" < ")+"atom"+graphAtom[2]);
*/
			if (graphRank[1] != graphRank[2])
				return true;

			if (level > 1)
//{
				cipCompileRelativeRanks(graphRank, graphParent, levelStart, level);
/*
System.out.print(" graphRank2:");
for (int i=0; i<levelStart[level]; i++) System.out.print(Integer.toHexString(graphRank[i])+" ");
System.out.print(" >> ");
for (int i=levelStart[level]; i<levelStart[level+1]; i++) System.out.print(Integer.toHexString(graphRank[i])+" ");
System.out.println("");
}*/
			}

		return false;
		}

	private int[] resize(int[] array, int newSize) {
		int [] copy = new int[newSize];
		System.arraycopy( array, 0, copy, 0, array.length );
		return copy;
		}

	private boolean[] resize(boolean[] array, int newSize) {
		boolean [] copy = new boolean[newSize];
		System.arraycopy( array, 0, copy, 0, array.length );
		return copy;
		}

	private void cipUpdateParentRanking(boolean[] graphIsPseudo,
										int[] graphRank,
										int[] graphParent,
										int[] graphAtom,
										int[] levelStart,
										int currentLevel) {
		class RankObject {
			int parentIndex;
			int parentRank;
			int parentHCount;
			int[] childRank;
			}

		for (int level=currentLevel; level>1; level--) {
			int parentCount = levelStart[level] - levelStart[level-1];
			RankObject[] rankObject = new RankObject[parentCount];
			int baseIndex = levelStart[level];
			for (int parent=0; parent<parentCount; parent++) {
				int parentIndex = levelStart[level-1] + parent;
				int nextBaseIndex = baseIndex;
				while (nextBaseIndex<levelStart[level+1]
					&& graphParent[nextBaseIndex] == parentIndex)
					nextBaseIndex++;
				rankObject[parent] = new RankObject();
				rankObject[parent].parentIndex = parentIndex;
				rankObject[parent].parentRank = graphRank[parentIndex];
				rankObject[parent].parentHCount = graphIsPseudo[parentIndex] ? 0
												: mMol.getPlainHydrogens(graphAtom[parentIndex]);
				rankObject[parent].childRank = new int[nextBaseIndex-baseIndex];
				for (int i=baseIndex; i<nextBaseIndex; i++)
					rankObject[parent].childRank[i-baseIndex] = graphRank[i];
				Arrays.sort(rankObject[parent].childRank);
				baseIndex = nextBaseIndex;
				}

			Comparator<RankObject> comparator = new Comparator<RankObject>() {
				public int compare(RankObject r1, RankObject r2) {
					if (r1.parentRank != r2.parentRank)
						return (r1.parentRank > r2.parentRank) ? 1 : -1;
					int i1 = r1.childRank.length;
					int i2 = r2.childRank.length;
					int count = Math.min(i1, i2);
					for (int i=0; i<count; i++) {
						i1--;
						i2--;
						if (r1.childRank[i1] != r2.childRank[i2])
							return (r1.childRank[i1] > r2.childRank[i2]) ? 1 : -1;
						}
					if (i1 != i2)
						return (i1 > i2) ? 1 : -1;
					if (r1.parentHCount != r2.parentHCount)
						return (r1.parentHCount > r2.parentHCount) ? 1 : -1;
					return 0;
					}
				};
			Arrays.sort(rankObject, comparator);

			int consolidatedRank = 1;
			for (int parent=0; parent<parentCount; parent++) {
				graphRank[rankObject[parent].parentIndex] = consolidatedRank;
				if (parent != parentCount-1
				 && comparator.compare(rankObject[parent], rankObject[parent+1]) != 0)
					consolidatedRank++;
				}
			}
		}

	private void cipCompileRelativeRanks(int[] graphRank, int[] graphParent, int[] levelStart, int currentLevel) {
		class RankObject {
			int rank;
			int parent;
			int index;
			}

		int levelOffset = levelStart[currentLevel];
		int count = levelStart[currentLevel+1] - levelOffset;
		RankObject[] rankObject = new RankObject[count];
		for (int i=0; i<count; i++) {
			rankObject[i] = new RankObject();
			rankObject[i].rank = graphRank[i+levelOffset];
			rankObject[i].parent = graphParent[i+levelOffset];
			rankObject[i].index = i+levelOffset;
			}

		Comparator<RankObject> comparator = new Comparator<RankObject>() {
			public int compare(RankObject r1, RankObject r2) {
				if (r1.rank != r2.rank)
					return (r1.rank > r2.rank) ? 1 : -1;
				return 0;
				}
			};

		for (int level=currentLevel; level>1; level--) {
			for (int i=0; i<count; i++) {
				rankObject[i].rank += (graphRank[rankObject[i].parent] << 16);
				rankObject[i].parent = graphParent[rankObject[i].parent];
				}
			Arrays.sort(rankObject, comparator);

			int consolidatedRank = 1;
			for (int i=0; i<count; i++) {
				graphRank[rankObject[i].index] = consolidatedRank;
				if (i != count-1
				 && comparator.compare(rankObject[i], rankObject[i+1]) != 0)
					consolidatedRank++;
				}
			}
		}
	
	/**
	 * @return an int[] giving all atom indexes in the order as they appear in the graph
	 */
	public int[] getGraphAtoms() {
		generateGraph();
		return mGraphAtom;
		}

	/**
	 * @return an int[] giving the relationship between new atom numbers and old atom numbers
	 */
	public int[] getGraphIndexes() {
		generateGraph();
		return mGraphIndex;
		}

	class ESRGroup implements Comparable<ESRGroup> {
		int[] atomList;
		int[] rankList;

		protected ESRGroup(int type, int group) {
			int count = 0;
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (mTHESRType[atom] == type
				 && mTHESRGroup[atom] == group)
					count++;
			atomList = new int[count];
			rankList = new int[count];
			count = 0;
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mTHESRType[atom] == type
				 && mTHESRGroup[atom] == group) {
					atomList[count] = atom;
					rankList[count++] = mCanRank[atom];
					}
				}
			Arrays.sort(rankList);
			}

		public int compareTo(ESRGroup g) {
			if (rankList.length != g.rankList.length)
				return rankList.length < g.rankList.length ? -1 : 1;
			for (int i=0; i<rankList.length; i++)
				if (rankList[i] != g.rankList[i])
					return rankList[i] < g.rankList[i] ? -1 : 1;
			return 0;
			}
		}
	}


class EZHalfParity {
	ExtendedMolecule mMol;
	int 	mCentralAxialAtom;
	int 	mRemoteAxialAtom;
	int		mHighConn;
	int		mLowConn;
	int		mValue;
	boolean	mStereoBondFound;
	boolean	mRanksEqual;
	boolean mInSameFragment;

	protected EZHalfParity(ExtendedMolecule mol, int[] rank, int atom1, int atom2) {
		mMol = mol;
		mRemoteAxialAtom = atom1;
		mCentralAxialAtom = atom2;
		int highRank = -1;
		for (int i=0; i<mMol.getAllConnAtoms(mCentralAxialAtom); i++) {
			int connAtom = mMol.getConnAtom(mCentralAxialAtom,i);
			int connBond = mMol.getConnBond(mCentralAxialAtom,i);

			if (connAtom == mRemoteAxialAtom) {
				if (mMol.getBondType(connBond) == Molecule.cBondTypeCross)
					mValue = -1;
				continue;
				}

			if (mMol.isStereoBond(connBond, mCentralAxialAtom)) {
				if (mStereoBondFound)
					mol.setStereoProblem(atom2);
				mStereoBondFound = true;
				}

			if (highRank == rank[connAtom]) {
				mLowConn = connAtom;	// two symmetrical atoms connected to one double bond atom
				mRanksEqual = true;
				mInSameFragment = mMol.isRingBond(connBond);
				}
			else if (highRank < rank[connAtom]) {
				highRank = rank[connAtom];
				mLowConn = mHighConn;
				mHighConn = connAtom;
				}
			else {
				mLowConn = connAtom;
				}
			}
		}


	protected int getValue() {
		// return values: -1:no assignment possible; 1:down; 2:left; 3:up; 4:right
		if (mValue != 0)
			return mValue;

		if (mStereoBondFound
		 && mMol.getAtomicNo(mCentralAxialAtom) != 15
		 && mMol.getAtomicNo(mCentralAxialAtom) != 16) {
			for (int i=0; i<mMol.getAllConnAtoms(mCentralAxialAtom); i++) {
				int connBond = mMol.getConnBond(mCentralAxialAtom,i);
				if (mMol.isStereoBond(connBond, mCentralAxialAtom)) {
					if (mMol.getConnAtom(mCentralAxialAtom,i) == mHighConn)
						mValue = (mMol.getBondType(connBond) == Molecule.cBondTypeUp) ? 3 : 1;
					else
						mValue = (mMol.getBondType(connBond) == Molecule.cBondTypeUp) ? 1 : 3;
					return mValue;
					}
				}
			}

		double angleDB = mMol.getBondAngle(mCentralAxialAtom,mRemoteAxialAtom);
		double angleHigh = mMol.getBondAngle(mCentralAxialAtom,mHighConn);
		if (angleHigh < angleDB)
			angleHigh += Math.PI*2;

		if (mMol.getAllConnAtoms(mCentralAxialAtom) == 2) {
			double angleDif = angleHigh - angleDB;
			if ((angleDif > Math.PI - 0.05) && (angleDif < Math.PI + 0.05)) {
				mValue = -1;	// less than 3 degrees different from double bond
				return mValue;	// is counted as non-stereo-specified double bond
				}
			mValue = (angleDif < Math.PI) ? 4 : 2;
			return mValue;
			}
		else {
			double angleOther = mMol.getBondAngle(mCentralAxialAtom,mLowConn);
			if (angleOther < angleDB)
				angleOther += Math.PI*2;
			mValue = (angleOther < angleHigh) ? 2 : 4;
			return mValue;
			}
		}
	}

class CanonizerBond implements Comparable<CanonizerBond> {
	int maxAtomRank,minAtomRank,bond;

	protected CanonizerBond(int atomRank1, int atomRank2, int bond) {
		maxAtomRank = Math.max(atomRank1, atomRank2);
		minAtomRank = Math.min(atomRank1, atomRank2);
		this.bond = bond;
		}

	public int compareTo(CanonizerBond cb) {
		if (maxAtomRank != cb.maxAtomRank)
			return maxAtomRank > cb.maxAtomRank ? -1 : 1;
		if (minAtomRank != cb.minAtomRank)
			return minAtomRank > cb.minAtomRank ? -1 : 1;	 // we want high ranks first
		return 0;
		}
	}

class CanonizerFragment {
	int[] atom;
	int[] bond;

	protected CanonizerFragment(int[] atom, int atoms, int[] bond, int bonds) {
		this.atom = Arrays.copyOf(atom, atoms);
		this.bond = Arrays.copyOf(bond, bonds);
		}
	}

class CanonizerParity implements Comparable<CanonizerParity> {
	int atom;
	int group;
	int rank;
	
	protected CanonizerParity(int atom, int group, int type, int rank) {
		this.atom = atom;
		this.group = group + (type << 8);
		this.rank = rank;
		}

	public int compareTo(CanonizerParity p) {
		if (group != p.group)
			return group < p.group ? -1 : 1;
		if (rank != p.rank)
			return rank > p.rank ? -1 : 1;	 // we want high ranks first
		return 0;
		}
	}
