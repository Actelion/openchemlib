/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
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

import java.util.*;

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

	// Consider diastereotopic/enantiotopic atoms uniquely for atom ranking
	public static final int CONSIDER_DIASTEREOTOPICITY = 2;
	public static final int CONSIDER_ENANTIOTOPICITY = 4;

	// Consider both diastereotopic and enantiotopic atoms uniquely for atom ranking
	public static final int CONSIDER_STEREOHETEROTOPICITY = CONSIDER_DIASTEREOTOPICITY | CONSIDER_ENANTIOTOPICITY;
	
	// Consider custom atom labels for atom ranking and encode them into idcodes
	public static final int ENCODE_ATOM_CUSTOM_LABELS = 8;

	// Consider the atom selection for atom ranking and encode it into idcodes
	public static final int ENCODE_ATOM_SELECTION = 16;

	// Assign parities to tetrahedral nitrogen (useful in crystals or at low temp, when N inversion is frozen out)
	public static final int ASSIGN_PARITIES_TO_TETRAHEDRAL_N = 32;

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

	public static final int ATOM_BITS = 16;
	public static final int MAX_ATOMS = 0xFFFF;
	public static final int MAX_BONDS = 0xFFFF;

	private ExtendedMolecule mMol;
	private int[] mCanRank;
	private int[] mCanRankBeforeTieBreaking;
	private byte[] mTHParity;
	private byte[] mEZParity;
	private byte[] mTHConfiguration;	// is tetrahedral parity based on atom numbers in graph
	private byte[] mEZConfiguration;	// is double bond parity based on atom numbers in graph
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
	private int mMode,mNoOfRanks;
	private boolean mIsOddParityRound;
	private boolean mZCoordinatesAvailable;
	private boolean mCIPParityNoDistinctionProblem;

	private boolean mGraphGenerated;
	private int mGraphRings;
	private int[] mGraphAtom;
	private int[] mGraphIndex;
	private int[] mGraphBond;
	private int[] mGraphFrom;
	private int[] mGraphClosure;

	private String		  mIDCode,mCoordinates,mMapping;
	private StringBuilder	mEncodingBuffer;
	private	int				mEncodingBitsAvail,mEncodingTempData;

	/**
	 * Runs a canonicalization process molecule creating a unique atom ranking
	 * taking stereo features, ESR settings and query features into account.
	 * @param mol
	 */
	public Canonizer(ExtendedMolecule mol) {
		this(mol, 0);
		}


	/**
	 * Runs a canonicalization process molecule creating a unique atom ranking
	 * taking stereo features, ESR settings and query features into account.
	 * If mode includes ENCODE_ATOM_CUSTOM_LABELS, than custom atom labels are
	 * used for atom ranking and are encoded into the idcode.
	 * 
	 * @param mol
	 * @param mode 0 or one or more of CONSIDER...TOPICITY, CREATE_SYMMETRY_RANK, ENCODE_ATOM_CUSTOM_LABELS, ASSIGN_PARITIES_TO_TETRAHEDRAL_N
	 */
	public Canonizer(ExtendedMolecule mol, int mode) {
		if (mol.getAllAtoms()>MAX_ATOMS)
			throw new IllegalArgumentException("Cannot canonize a molecule having more than "+MAX_ATOMS+" atoms");
		if (mol.getAllBonds()>MAX_BONDS)
			throw new IllegalArgumentException("Cannot canonize a molecule having more than "+MAX_BONDS+" bonds");

		mMol = mol;
		mMode = mode;

		mMol.ensureHelperArrays(Molecule.cHelperRings);
		canFindNitrogenQualifyingForParity();

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.getAtomZ(atom) != 0.0) {
				mZCoordinatesAvailable = true;
				break;
				}
			}

		mCanRank = new int[mMol.getAllAtoms()];
		mCanBase = new CanonizerBaseValue[mMol.getAtoms()];
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			mCanBase[atom] = new CanonizerBaseValue();
			
		mTHParity = new byte[mMol.getAtoms()];
		mTHParityIsPseudo = new boolean[mMol.getAtoms()];
		mTHParityRoundIsOdd = new boolean[mMol.getAtoms()];
		mEZParity = new byte[mMol.getBonds()];
		mEZParityRoundIsOdd = new boolean[mMol.getBonds()];
		mEZParityIsPseudo = new boolean[mMol.getBonds()];

		mCIPParityNoDistinctionProblem = false;

		canInitializeRanking();
		canRankStereo();
		canRankFinal();

		if (mCIPParityNoDistinctionProblem)
			System.out.println("No distinction applying CIP rules: "+getIDCode()+" "+getEncodedCoordinates());
		}


	/**
	 * Locate those tetrahedral nitrogen atoms with at least 3 neighbors that
	 * qualify for tetrahedral parity calculation because:<br>
	 * - they are quarternary nitrogen atoms<br>
	 * or - their configuration inversion is hindered in a polycyclic structure<br>
	 * or - flag ASSIGN_PARITIES_TO_TETRAHEDRAL_N is set
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
		int[] canRankWithoutStereo = new int[mMol.getAtoms()];
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			canRankWithoutStereo[atom] = mCanRank[atom];

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
		mTHParityNeedsNormalization = new boolean[mMol.getAtoms()];
		mTHParityNormalizationGroupList = new ArrayList<int[]>();
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
			// locate atom differences due to pro-chiral or pro-E/Z location and
			// detect for every proTH- or proEZ-parity whether pro-atoms are
			// in same fragment as the pro-chiral-center or double-bond, respectively
		mProTHAtomsInSameFragment = new boolean[mMol.getAtoms()];
		mProEZAtomsInSameFragment = new boolean[mMol.getBonds()];

		if ((mMode & CONSIDER_STEREOHETEROTOPICITY) != 0) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(ATOM_BITS+4, mCanRank[atom] << 12);
				}
			}
		if (mNoOfRanks < mMol.getAtoms()) {
			int proParities = 0;
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (canCalcTHParity(atom, true))
					proParities++;
			for (int bond=0; bond<mMol.getBonds(); bond++)
				if (canCalcEZParity(bond, true))
					proParities++;
			}
		if ((mMode & CONSIDER_STEREOHETEROTOPICITY) != 0) {
			mNoOfRanks = canPerformRanking();
			}

			// Atoms, which still share an equal rank, can now be
			// considered symmetrical considering their connectivity
			// and stereo features.
		if ((mMode & CREATE_SYMMETRY_RANK) != 0) {
			mCanRankBeforeTieBreaking = new int[mMol.getAtoms()];
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				mCanRankBeforeTieBreaking[atom] = mCanRank[atom];
			}


			// ############### begin tie breaking ##############
			// i.e. if not all atoms have a different rank yet, then
			// select randomly one atom of those atoms that share the
			// lowest rank and consider it higher ranked. Propagate the
			// new ranking assymetry through the molecule (performFullRanking()).
			// Repeat these steps until all atoms have a different rank.

//System.out.println("start of tie breaking");
/*
for (int atom=0; atom<mMol.getAtoms(); atom++)
System.out.println("mTHParity["+atom+"] = "+mTHParity[atom]+" round = "+mTHParityRound[atom]+" pseudo = "+mTHParityIsPseudo[atom]);
for (int bond=0; bond<mMol.getBonds(); bond++)
System.out.println("mEZParity["+bond+"] = "+mEZParity[bond]);
*/
		while (mNoOfRanks < mMol.getAtoms()) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(ATOM_BITS+1, 2*mCanRank[atom]);
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
			canNormalizeGroupParities();
			if (mMesoHelper != null)
				mMesoHelper.normalizeESRGroupSwappingAndRemoval(mCanRank);
			}

			// normalize those group's parities that could not be normalized before tie-breaking
		canNormalizeGroupParities();

			// detect all not yet discovered pseudo-parities
		canFindPseudoParities();
		flagStereoProblems();
/*
for (int atom=0; atom<mMol.getAtoms(); atom++)
System.out.println("mTHParity["+atom+"] = "+mTHParity[atom]+" round = "+mTHParityRound[atom]+" pseudo = "+mTHParityIsPseudo[atom]);
for (int bond=0; bond<mMol.getBonds(); bond++)
System.out.println("mEZParity["+bond+"] = "+mEZParity[bond]);
*/
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


	private void canInitializeRanking() {
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
			mCanBase[atom].add(3, mMol.getConnAtoms(atom));
			if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0)
				mCanBase[atom].add(4, 8);
			else
				mCanBase[atom].add(4, 8 + mMol.getAtomCharge(atom));
			mCanBase[atom].add(5, mMol.getAtomRingSize(atom));
			mCanBase[atom].add(4, canCalcImplicitAbnormalValence(atom)+1);
			mCanBase[atom].add(2, mMol.getAtomRadical(atom) >> Molecule.cAtomRadicalStateShift);
			mCanBase[atom].add(Molecule.cAtomQFNoOfBits, mMol.getAtomQueryFeatures(atom));
			if (mMol.getAtomList(atom) != null)
				atomListFound = true;
			}

		mNoOfRanks = canPerformRanking();

		if (atomListFound) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(ATOM_BITS, mCanRank[atom]);
				int[] atomList = mMol.getAtomList(atom);
				int listLength = (atomList == null) ? 0 : atomList.length;
				mCanBase[atom].add((6 - listLength)*8, 0);
				for (int i=listLength-1; i>=0; i--)
					mCanBase[atom].add(8, atomList[i]);
				}

			mNoOfRanks = canPerformRanking();
			}

		if (mMol.isFragment()) {
			boolean bondQueryFeaturesPresent = false;
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				if (mMol.getBondQueryFeatures(bond) != 0) {
					bondQueryFeaturesPresent = true;
					break;
					}
				}
			if (bondQueryFeaturesPresent) {
				for (int atom=0; atom<mMol.getAtoms(); atom++) {
					mCanBase[atom].init(atom);
					mCanBase[atom].add(ATOM_BITS, mCanRank[atom]);
					int[] bondQFList = new int[mMol.getConnAtoms(atom)];
					for (int i=0; i<mMol.getConnAtoms(atom); i++)
						bondQFList[i] = mMol.getBondQueryFeatures(mMol.getConnBond(atom, i));
					Arrays.sort(bondQFList);
					mCanBase[atom].add((6 - bondQFList.length)*Molecule.cBondQFNoOfBits, 0);
					for (int i=bondQFList.length-1; i>=0; i--)
						mCanBase[atom].add(Molecule.cBondQFNoOfBits, bondQFList[i]);
					}
				mNoOfRanks = canPerformRanking();
				}
			}

		if ((mMode & ENCODE_ATOM_CUSTOM_LABELS) != 0) {
			SortedStringList list = new SortedStringList();
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (mMol.getAtomCustomLabel(atom) != null)
					list.addString(mMol.getAtomCustomLabel(atom));

			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				int rank = (mMol.getAtomCustomLabel(atom) == null) ?
						   0 : 1+list.getListIndex(mMol.getAtomCustomLabel(atom));
				mCanBase[atom].init(atom);
				mCanBase[atom].add(ATOM_BITS, mCanRank[atom]);
				mCanBase[atom].add(ATOM_BITS, rank);
				}

			mNoOfRanks = canPerformRanking();
			}

		if ((mMode & ENCODE_ATOM_SELECTION) != 0) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(ATOM_BITS, mCanRank[atom]);
				mCanBase[atom].add(1, mMol.isSelectedAtom(atom) ? 1 : 0);
				}

			mNoOfRanks = canPerformRanking();
			}
//System.out.println("after initial ranking");
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


	private void canRecursivelyFindAllParities() {
		// detects all stereo centers and their parities while considering 
		// absolute ESR group numbers and the drawn configuration of stereo
		// centers within OR or AND groups. The calculated parities are not
		// canonized and only intended for meso fragment detection and complete
		// detection of all stereo centers.
		mIsOddParityRound = true;

		// Find absolute stereo features based on current ranking
		boolean paritiesFound = canFindParities(false);

			// in a loop check for stereo features that depend on the configuration
			// of other stereo features already found and rank atoms again
		final int parityInfoBits = 2 + 2 + Molecule.cESRGroupBits;
		while ((mNoOfRanks < mMol.getAtoms()) && paritiesFound) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				mCanBase[atom].init(atom);
				mCanBase[atom].add(ATOM_BITS, mCanRank[atom]);

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


	private void canRecursivelyFindCIPParities() {
		// Detects stereo center's parities considering the drawn configuration.
		// Any ESR group assignments are neglected. From these parities the CIP
		// assignments are calculated. These may be for display purposes only.

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
				mCanBase[atom].add(ATOM_BITS+4, (mCanRank[atom] << 4)
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


	private void canRecursivelyFindCanonizedParities() {
			// detects parities for full normalization and idcode generation
			// (does not consider absolute group numbers nor drawn configuration of
			//  stereo centers within OR or AND groups, before they are normalized)
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
				mCanBase[atom].add(ATOM_BITS, mCanRank[atom]);
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
									+  (groupRank[(mTHESRType[atom] == Molecule.cESRTypeAnd) ? 0 : 1][mTHESRGroup[atom]] << 8));

//				if (!mTHParityNeedsNormalization[atom])
					mCanBase[atom].add(mTHParity[atom] << 4);
				}

// TODO consider groupRank for bonds
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				mCanBase[mMol.getBondAtom(0,bond)].add(mEZParity[bond]);
				mCanBase[mMol.getBondAtom(1,bond)].add(mEZParity[bond]);
				}
/*
for (int atom=0; atom<mMol.getAtoms(); atom++)
System.out.println("mCanBaseValue["+atom+"] = "+Long.toHexString(mCanBaseValue[atom])+" parity:"+mTHParity[atom]+" needsNorm:"+mTHParityNeedsNormalization[atom]+" ESRType:"+mTHESRType[atom]);
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
	 * This normalizes relative parities within any ESR group
	 * such that the highest ranking atom or bond gets parity2
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
				for (int j=0; j<groupAtom.length; j++) {
					int atom = groupAtom[j];
					if (invertParities) {
						if (mTHParity[atom] == Molecule.cAtomParity1)
							mTHParity[atom] = Molecule.cAtomParity2;
						else if (mTHParity[atom] == Molecule.cAtomParity2)
							mTHParity[atom] = Molecule.cAtomParity1;
						}
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
			if (mTHESRType[atom] != Molecule.cESRTypeAbs)
				count++;

		if (count == 0)
			return;

		int[] parity = new int[count];
		count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mTHESRType[atom] != Molecule.cESRTypeAbs) {
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
			if (canCalcEZParity(bond, false)) {
//System.out.println("found EZ parity:"+mEZParity[bond]+" bond:"+bond+" round:"+(mRealParityRounds+1));
				mEZParityRoundIsOdd[bond] = mIsOddParityRound;
				if (doCIP)
					cipCalcEZParity(bond);
				ezFound = true;
				}
		boolean thFound = false;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (canCalcTHParity(atom, false)) {
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
		boolean isFreshPseudoParityAtom[] = new boolean[mMol.getAtoms()];
		boolean isFreshPseudoParityBond[] = new boolean[mMol.getBonds()];
		int anyPseudoParityCount = 0;
		boolean pseudoParity1Or2Found = false;

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mProTHAtomsInSameFragment[atom]) {
				if (!mTHParityIsPseudo[atom]) {
					if (canCalcTHParity(atom, false)) {
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
					if (canCalcEZParity(bond, false)) {
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

		mFragmentList = new ArrayList<CanonizerFragment>();

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
				boolean bondHandled[] = new boolean[mMol.getBonds()];
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
		int	connRank[] = new int[ExtendedMolecule.cMaxConnAtoms];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
								// generate sorted list of ranks of neighbours
			for (int i=0; i<mMol.getConnAtoms(atom); i++) {
				int rank = 2 * mCanRank[mMol.getConnAtom(atom,i)];
				int connBond = mMol.getConnBond(atom,i);
				if (mMol.getBondOrder(connBond) == 2)
					if (!mMol.isAromaticBond(connBond))
						rank++;		// set a flag for non-aromatic double bond
				int j;
				for (j=0; j<i; j++)
					if (rank < connRank[j])
						break;
				for (int k=i; k>j; k--)
					connRank[k] = connRank[k-1];
				connRank[j] = rank;
				}

			int neighbours = Math.min(6, mMol.getConnAtoms(atom));
			mCanBase[atom].init(atom);
			mCanBase[atom].add(ATOM_BITS, mCanRank[atom]);
			mCanBase[atom].add((6 - neighbours)*(ATOM_BITS + 1), 0);
			for (int i=0; i<neighbours; i++)
				mCanBase[atom].add(ATOM_BITS + 1, connRank[i]);
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
	 * @param calcProParity if true then calculate pro-parities
	 * @return false if atom's mTHParity != 0; otherwise true if it finds and sets a parity other than 0
	 */
	private boolean canCalcTHParity(int atom, boolean calcProParity) {

		//#		if (!calcProParity && mTHParity[atom] != 0)
		if (mTHParity[atom] != 0)
			return false;

		if (mMol.getAtomicNo(atom) != 6
		 && mMol.getAtomicNo(atom) != 7
		 && mMol.getAtomicNo(atom) != 14
		 && mMol.getAtomicNo(atom) != 15
		 && mMol.getAtomicNo(atom) != 16)
			return false;

		if (mMol.getAtomPi(atom) != 0) {
			if (isCentralAlleneAtom(atom))
				return canCalcAlleneParity(atom, calcProParity);

			if (mMol.getAtomicNo(atom) != 15
			 && mMol.getAtomicNo(atom) != 16)
				return false;
			}

		if (mMol.getConnAtoms(atom) < 3 || mMol.getAllConnAtoms(atom) > 4)
			return false;

		// don't tetrahedral nitrogen, unless they were found to qualify for parity calculation
		if (mMol.getAtomicNo(atom) == 7
		 && !mNitrogenQualifiesForParity[atom])
			return false;

				// create array to remap connAtoms according to canRank order
		int remappedConn[] = new int[4];
		int remappedRank[] = new int[4];
		boolean neighbourUsed[] = new boolean[4];
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
				if (!calcProParity || remappedRank[i] == 0)	// don't calculate pro-parities or hydrogens
					return false;	// two equal ranking atoms -> no stereo center

				proTHAtom1 = mMol.getConnAtom(atom, remappedConn[i-1]);
				proTHAtom2 = mMol.getConnAtom(atom, remappedConn[i]);
				if (mMol.isRingBond(mMol.getConnBond(atom,remappedConn[i])))
					mProTHAtomsInSameFragment[atom] = true;
				proTHAtomsFound = true;
				}
			}

		if (calcProParity && !proTHAtomsFound)
			return false;

		byte atomTHParity = (mZCoordinatesAvailable) ?
							canCalcTHParity3D(atom, remappedConn)
						  : canCalcTHParity2D(atom, remappedConn);

		if (!calcProParity) {
			mTHParity[atom] = atomTHParity;
			}
		else if ((mStereoCentersFound && (mMode & CONSIDER_DIASTEREOTOPICITY) != 0)
			  || (!mStereoCentersFound && (mMode & CONSIDER_ENANTIOTOPICITY) != 0)) {
			// increment mCanBase[] for atoms that are Pro-Parity1
			if (atomTHParity == Molecule.cAtomParity1) {
				mCanBase[proTHAtom1].add(0x0400);
				mCanBase[proTHAtom2].add(0x0100);
				}
			else if (atomTHParity == Molecule.cAtomParity2) {
				mCanBase[proTHAtom1].add(0x0100);
				mCanBase[proTHAtom2].add(0x0400);
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

		float angle[] = new float[mMol.getAllConnAtoms(atom)];
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

		float[][] coords = new float[3][3];
		for (int i=0; i<3; i++) {
			coords[i][0] = mMol.getAtomX(atomList[i+1]) - mMol.getAtomX(atomList[0]);
			coords[i][1] = mMol.getAtomY(atomList[i+1]) - mMol.getAtomY(atomList[0]);
			coords[i][2] = mMol.getAtomZ(atomList[i+1]) - mMol.getAtomZ(atomList[0]);
			}

		// calculate the normal vector (vector product of coords[0] and coords[1])
		float[] n = new float[3];
		n[0] = coords[0][1]*coords[1][2]-coords[0][2]*coords[1][1];
		n[1] = coords[0][2]*coords[1][0]-coords[0][0]*coords[1][2];
		n[2] = coords[0][0]*coords[1][1]-coords[0][1]*coords[1][0];

		// calculate cos(angle) of coords[2] to normal vector
		float cosa = (coords[2][0]*n[0]+coords[2][1]*n[1]+coords[2][2]*n[2])
					/ ((float)Math.sqrt(coords[2][0]*coords[2][0]+coords[2][1]*coords[2][1]+coords[2][2]*coords[2][2])
					 * (float)Math.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]));

		return (cosa > 0.0) ? (byte)Molecule.cAtomParity1 : Molecule.cAtomParity2;
		}


	private boolean canCalcAlleneParity(int atom, boolean calcProParity) {
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
		if (halfParity1.mRanksEqual && !calcProParity)
			return false;
		EZHalfParity halfParity2 = new EZHalfParity(mMol, mCanRank,atom, atom2);
		if (halfParity2.mRanksEqual && !calcProParity)
			return false;

		if (halfParity1.mRanksEqual && halfParity2.mRanksEqual)
			return false;	// both ends of DB bear equal substituents

		if (calcProParity) {
			if (halfParity1.mRanksEqual && halfParity1.mInSameFragment)
				mProTHAtomsInSameFragment[atom] = true;
			if (halfParity2.mRanksEqual && halfParity2.mInSameFragment)
				mProTHAtomsInSameFragment[atom] = true;
			}

		int hp1 = halfParity1.getValue();
		int hp2 = halfParity2.getValue();
		if (hp1 == -1 || hp2 == -1 || ((hp1 + hp2) & 1) == 0) {
			if (!calcProParity) {
				mTHParity[atom] = Molecule.cAtomParityUnknown;
				}
			return true;
			}

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

		if (!calcProParity) {	// increment mProParity[] for atoms that are Pro-Parity1
			mTHParity[atom] = alleneParity;
			}
		else if ((mStereoCentersFound && (mMode & CONSIDER_DIASTEREOTOPICITY) != 0)
			  || (!mStereoCentersFound && (mMode & CONSIDER_ENANTIOTOPICITY) != 0)) {
			if (halfParity1.mRanksEqual) {
				if (alleneParity == Molecule.cAtomParity1) {
					mCanBase[halfParity1.mHighConn].add(0x0040);
					mCanBase[halfParity1.mLowConn].add(0x0010);
					}
				else {
					mCanBase[halfParity1.mHighConn].add(0x0010);
					mCanBase[halfParity1.mLowConn].add(0x0040);
					}
				}
			if (halfParity2.mRanksEqual) {
				if (alleneParity == Molecule.cAtomParity2) {
					mCanBase[halfParity2.mHighConn].add(0x0040);
					mCanBase[halfParity2.mLowConn].add(0x0010);
					}
				else {
					mCanBase[halfParity2.mHighConn].add(0x0010);
					mCanBase[halfParity2.mLowConn].add(0x0040);
					}
				}
			}

		return true;
		}


	private boolean canCalcBINAPParity(int bond, boolean calcProParity) {
		if (!mMol.isBINAPChiralityBond(bond))
			return false;

		EZHalfParity halfParity1 = new EZHalfParity(mMol, mCanRank, mMol.getBondAtom(0, bond), mMol.getBondAtom(1, bond));
		if (halfParity1.mRanksEqual && !calcProParity)
			return false;
		EZHalfParity halfParity2 = new EZHalfParity(mMol, mCanRank, mMol.getBondAtom(1, bond), mMol.getBondAtom(0, bond));
		if (halfParity2.mRanksEqual && !calcProParity)
			return false;

		if (halfParity1.mRanksEqual && halfParity2.mRanksEqual)
			return false;	// both ends of DB bear equal substituents

		if (calcProParity) {
			if (halfParity1.mRanksEqual || halfParity2.mRanksEqual)
				mProEZAtomsInSameFragment[bond] = true;
			}

		int hp1 = halfParity1.getValue();
		int hp2 = halfParity2.getValue();
		if (hp1 == -1 || hp2 == -1 || ((hp1 + hp2) & 1) == 0) {
			if (!calcProParity) {
				mEZParity[bond] = Molecule.cBondParityUnknown;
				}
			return true;
			}

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

		if (!calcProParity) {	// increment mProParity[] for atoms that are Pro-E
			mEZParity[bond] = axialParity;
			}
		else if ((mStereoCentersFound && (mMode & CONSIDER_DIASTEREOTOPICITY) != 0)
			  || (!mStereoCentersFound && (mMode & CONSIDER_ENANTIOTOPICITY) != 0)) {
			if (halfParity1.mRanksEqual) {
				if (axialParity == Molecule.cBondParityZor2) {
					mCanBase[halfParity1.mHighConn].add(0x0004);
					mCanBase[halfParity1.mLowConn].add(0x0001);
					}
				else {
					mCanBase[halfParity1.mHighConn].add(0x0001);
					mCanBase[halfParity1.mLowConn].add(0x0004);
					}
				}
			if (halfParity2.mRanksEqual) {
				if (axialParity == Molecule.cBondParityZor2) {
					mCanBase[halfParity2.mHighConn].add(0x0004);
					mCanBase[halfParity2.mLowConn].add(0x0001);
					}
				else {
					mCanBase[halfParity2.mHighConn].add(0x0001);
					mCanBase[halfParity2.mLowConn].add(0x0004);
					}
				}
			}

		return true;
		}


	/**
	 * Calculate E/Z parity of a double bond based on their neighbours ranks.
	 * If calcProParity is true then the pro-parities of the allylic atoms are calculated
	 * rather than EZ-parities of double bonds.
	 * @param bond
	 * @param calcProParity if true then calculate pro-parities
	 * @return false if bond's mEZParity != 0; otherwise true if it finds and sets a parity other than 0
	 */
	private boolean canCalcEZParity(int bond, boolean calcProParity) {
			// depending on the flag calcProParity this routine either
			//	- updates EZ-parities of double bonds
			// or - calculates EZ-pro-parities of allylic atoms
//#		if (!calcProParity && mEZParity[bond] != 0)
		if (mEZParity[bond] != 0)
			return false;

		if (mMol.getBondOrder(bond) == 1)
			return canCalcBINAPParity(bond, calcProParity);

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
		if (halfParity1.mRanksEqual && !calcProParity)
			return false;
		EZHalfParity halfParity2 = new EZHalfParity(mMol, mCanRank, dbAtom1, dbAtom2);
		if (halfParity2.mRanksEqual && !calcProParity)
			return false;

		if (halfParity1.mRanksEqual && halfParity2.mRanksEqual)
			return false;	// both ends of DB bear equal substituents

		if (calcProParity) {
			if (halfParity1.mRanksEqual && halfParity1.mInSameFragment)
				mProEZAtomsInSameFragment[bond] = true;
			if (halfParity2.mRanksEqual && halfParity2.mInSameFragment)
				mProEZAtomsInSameFragment[bond] = true;
			}

		byte bondDBParity = mMol.isBondParityUnknownOrNone(bond) ? Molecule.cBondParityUnknown
							  		  : (mZCoordinatesAvailable) ? canCalcEZParity3D(halfParity1, halfParity2)
							  				  					 : canCalcEZParity2D(halfParity1, halfParity2);

		if (!calcProParity) {
			mEZParity[bond] = bondDBParity;
			}
		else if ((mMode & CONSIDER_DIASTEREOTOPICITY) != 0) {
			// increment mProParity[] for atoms that are Pro-E
			if (halfParity1.mRanksEqual) {
				if (bondDBParity == Molecule.cBondParityEor1) {
					mCanBase[halfParity1.mHighConn].add(0x0004);
					mCanBase[halfParity1.mLowConn].add(0x0001);
					}
				else if (bondDBParity == Molecule.cBondParityZor2) {
					mCanBase[halfParity1.mHighConn].add(0x0001);
					mCanBase[halfParity1.mLowConn].add(0x0004);
					}
				}
			if (halfParity2.mRanksEqual) {
				if (bondDBParity == Molecule.cBondParityEor1) {
					mCanBase[halfParity2.mHighConn].add(0x0004);
					mCanBase[halfParity2.mLowConn].add(0x0001);
					}
				else if (bondDBParity == Molecule.cBondParityZor2) {
					mCanBase[halfParity2.mHighConn].add(0x0001);
					mCanBase[halfParity2.mLowConn].add(0x0004);
					}
				}
			}

		return true;
		}


	private boolean isCentralAlleneAtom(int atom) {
		return mMol.getConnAtoms(atom) == 2
			&& mMol.getConnBondOrder(atom,0) == 2
			&& mMol.getConnBondOrder(atom,1) == 2;
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
		float[] db = new float[3];
		db[0] = mMol.getAtomX(halfParity2.mCentralAxialAtom) - mMol.getAtomX(halfParity1.mCentralAxialAtom);
		db[1] = mMol.getAtomY(halfParity2.mCentralAxialAtom) - mMol.getAtomY(halfParity1.mCentralAxialAtom);
		db[2] = mMol.getAtomZ(halfParity2.mCentralAxialAtom) - mMol.getAtomZ(halfParity1.mCentralAxialAtom);

		float[] s1 = new float[3];
		s1[0] = mMol.getAtomX(halfParity1.mHighConn) - mMol.getAtomX(halfParity1.mCentralAxialAtom);
		s1[1] = mMol.getAtomY(halfParity1.mHighConn) - mMol.getAtomY(halfParity1.mCentralAxialAtom);
		s1[2] = mMol.getAtomZ(halfParity1.mHighConn) - mMol.getAtomZ(halfParity1.mCentralAxialAtom);

		float[] s2 = new float[3];
		s2[0] = mMol.getAtomX(halfParity2.mHighConn) - mMol.getAtomX(halfParity2.mCentralAxialAtom);
		s2[1] = mMol.getAtomY(halfParity2.mHighConn) - mMol.getAtomY(halfParity2.mCentralAxialAtom);
		s2[2] = mMol.getAtomZ(halfParity2.mHighConn) - mMol.getAtomZ(halfParity2.mCentralAxialAtom);

		// calculate the normal vector n1 of plane from db and s1 (vector product)
		float[] n1 = new float[3];
		n1[0] = db[1]*s1[2]-db[2]*s1[1];
		n1[1] = db[2]*s1[0]-db[0]*s1[2];
		n1[2] = db[0]*s1[1]-db[1]*s1[0];

		// calculate the normal vector n2 of plane from db and n1 (vector product)
		float[] n2 = new float[3];
		n2[0] = db[1]*n1[2]-db[2]*n1[1];
		n2[1] = db[2]*n1[0]-db[0]*n1[2];
		n2[2] = db[0]*n1[1]-db[1]*n1[0];

		// calculate cos(angle) of s1 and normal vector n2
		float cosa = (s1[0]*n2[0]+s1[1]*n2[1]+s1[2]*n2[2])
					/ ((float)Math.sqrt(s1[0]*s1[0]+s1[1]*s1[1]+s1[2]*s1[2])
					 * (float)Math.sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]));

		// calculate cos(angle) of s2 and normal vector n2
		float cosb = (s2[0]*n2[0]+s2[1]*n2[1]+s2[2]*n2[2])
					/ ((float)Math.sqrt(s2[0]*s2[0]+s2[1]*s2[1]+s2[2]*s2[2])
					 * (float)Math.sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]));

		return ((cosa < 0.0) ^ (cosb < 0.0)) ? (byte)Molecule.cBondParityEor1 : Molecule.cBondParityZor2;
		}


	private void flagStereoProblems() {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			// if stereo center is declared unknown and no recognized as such; or vice versa
			if (mMol.isAtomConfigurationUnknown(atom)
			  ^ mTHParity[atom] == Molecule.cAtomParityUnknown)
				mMol.setStereoProblem(atom);

			// if no parity found, but atom was assigned to AND or OR group
			if ((mMol.getAtomESRType(atom) == Molecule.cESRTypeAnd
			  || mMol.getAtomESRType(atom) == Molecule.cESRTypeOr)
			 && (!mIsStereoCenter[atom]
			  || mTHParity[atom] == Molecule.cAtomParityUnknown))
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
			 && mEZParity[bond] == Molecule.cBondParityUnknown) {
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

		boolean atomHandled[] = new boolean[mMol.getAtoms()];
		boolean bondHandled[] = new boolean[mMol.getBonds()];
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
					for (int i=0; i<mMol.getConnAtoms(mGraphAtom[firstUnhandled]); i++) {
						int connAtom = mMol.getConnAtom(mGraphAtom[firstUnhandled],i);
						if (!atomHandled[connAtom] && mCanRank[connAtom] > highestRank) {
							highestRankingConnAtom = connAtom;
							highestRankingConnBond = mMol.getConnBond(mGraphAtom[firstUnhandled],i);
							highestRank = mCanRank[connAtom];
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
		generateGraph();

		StereoMolecule mol = new StereoMolecule(mMol.getAtoms(), mMol.getBonds());

		for(int i=0; i<mMol.getAtoms(); i++) {
			mMol.copyAtom(mol, mGraphAtom[i], 0, 0);
			mol.setAtomESR(i, mTHESRType[mGraphAtom[i]], mTHESRGroup[mGraphAtom[i]]);
			}

		for(int i=0; i<mMol.getBonds(); i++) {
			mMol.copyBond(mol, mGraphBond[i], 0, 0, mGraphIndex, false);
			mol.setBondESR(i, mEZESRType[mGraphBond[i]], mEZESRGroup[mGraphBond[i]]);
			}

		mMol.copyMoleculeProperties(mol);
		mMol.invalidateHelperArrays(Molecule.cHelperBitParities);

		return mol;
		}


	/**
	 * Sets all atoms with TH-parity 'unknown' to explicitly defined 'unknown'.
	 * Sets all double bonds with EZ-parity 'unknown' to cross bonds.
	 */
	public void setUnknownParitiesToExplicitlyUnknown() {
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (!mMol.isAtomConfigurationUnknown(atom)
			 && mTHParity[atom] == Molecule.cAtomParityUnknown)
				mMol.setAtomConfigurationUnknown(atom, true);
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mEZParity[bond] == Molecule.cBondParityUnknown) {
				if (mMol.getBondOrder(bond) == 2) {
			   		mMol.setBondType(bond, Molecule.cBondTypeCross);
					}
				// TODO once there is a BINAP config unknown, adapt here
				}
			}
		}


	public void setSingleUnknownAsRacemicParity() {
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
									return;
							}
						}
					else {
						// tetrahedral chirality
						for (int i=0; i<mMol.getConnAtoms(atom); i++)
							if (mMol.isStereoBond(mMol.getConnBond(atom, i)))
								return;
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
					return;
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
								return;
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
					return;
					}
				}
			}
		}


	public String getIDCode() {
		if (mIDCode == null) {
			generateGraph();
			idGenerateConfigurations();
//			idNormalizeConfigurations();
			idNormalizeESRGroupNumbers();
			idCodeCreate();
			}

		return mIDCode;
		}


	public int[] getFinalRank() {
		// this is the final mCanRank after all tie breaking steps
		return mCanRank;
		}


	/**
	 * Returned symmetry rank before tie breaking. For this the Canonizer
	 * mode must contain the CREATE_SYMMETRY_RANK option. If ranking
	 * shall reflect atom diastereotopicity or even enantiotopicity, use
	 * mode CONSIDER_DIASTEREOTOPICITY or CONSIDER_STEREOHETEROTOPICITY,
	 * respectively.
	 * @param atom
	 * @return rank
	 */
	public int getSymmetryRank(int atom) {
		// requires to use the 'createSymmetryRank' option
		return (mCanRankBeforeTieBreaking == null) ? -1 : mCanRankBeforeTieBreaking[atom];
		}


	private void idCodeCreate() {
		encodeBitsStart();
		encodeBits(cIDCodeVersion3, 4);
		int nbits = Math.max(idGetNeededBits(mMol.getAtoms()),
							 idGetNeededBits(mMol.getBonds()));
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
		int dbits = idGetNeededBits(maxdif);
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
			int bondOrder = ((mMol.getBondQueryFeatures(bond) & Molecule.cBondQFBridge) != 0) ?
							1 : (mMol.isDelocalizedBond(mGraphBond[bond])) ?
							0 : mMol.getBondOrder(mGraphBond[bond]);
			encodeBits(bondOrder, 2);
//System.out.print(bondOrder + ";");
			}
//System.out.println();

		int THCount = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mTHConfiguration[mGraphAtom[atom]] != Molecule.cAtomParityNone
			 && mTHConfiguration[mGraphAtom[atom]] != Molecule.cAtomParityUnknown)
				THCount++;
		encodeBits(THCount, nbits);
		for (int atom=0; atom<mMol.getAtoms(); atom++)
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

		int EZCount = 0;
		for (int bond=0; bond<mMol.getBonds(); bond++)	// parity of all double bonds
			if (mEZConfiguration[mGraphBond[bond]] != 0
			 && mEZConfiguration[mGraphBond[bond]] != Molecule.cBondParityUnknown
			 && (!mMol.isSmallRingBond(mGraphBond[bond]) || mMol.getBondType(mGraphBond[bond]) == Molecule.cBondTypeSingle))
				EZCount++;
		encodeBits(EZCount, nbits);
		for (int bond=0; bond<mMol.getBonds(); bond++)
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
			encodeBits(1, 1);	//	more data to come
			encodeBits(1, 4);	//	1 = datatype 'isotope'
			encodeBits(count, nbits);
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mMol.getAtomMass(mGraphAtom[atom]) != 0) {
					encodeBits(atom, nbits);
					encodeBits(mMol.getAtomMass(mGraphAtom[atom]), 8);
					}
				}
			}

		boolean isSecondFeatureBlock = false;

		if (mMol.isFragment()) {	// QueryFeatures and fragment specific properties
			addAtomQueryFeatures(0, false, nbits, Molecule.cAtomQFNoMoreNeighbours, 1, -1);

			addBondQueryFeatures(2, false, nbits, Molecule.cBondTypeDelocalized, 1, -1);

			addAtomQueryFeatures(3, false, nbits, Molecule.cAtomQFMoreNeighbours, 1, -1);

			addAtomQueryFeatures(4, false, nbits,
								 Molecule.cAtomQFRingState,
								 Molecule.cAtomQFRingStateBits,
								 Molecule.cAtomQFRingStateShift);

			addAtomQueryFeatures(5, false, nbits,
								 Molecule.cAtomQFAromState,
								 Molecule.cAtomQFAromStateBits,
								 Molecule.cAtomQFAromStateShift);

			addAtomQueryFeatures(6, false, nbits, Molecule.cAtomQFAny, 1, -1);

			addAtomQueryFeatures(7, false, nbits,
								 Molecule.cAtomQFHydrogen,
								 Molecule.cAtomQFHydrogenBits,
								 Molecule.cAtomQFHydrogenShift);

			count = 0;
			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (mMol.getAtomList(mGraphAtom[atom]) != null)
					count++;
			if (count > 0) {
				encodeBits(1, 1);	//	more data to come
				encodeBits(8, 4);	//	8 = datatype 'AtomList'
				encodeBits(count, nbits);
				for (int atom=0; atom<mMol.getAtoms(); atom++) {
					int[] atomList = mMol.getAtomList(mGraphAtom[atom]);
					if (atomList != null) {
						encodeBits(atom, nbits);
						encodeBits(atomList.length, 4);
						for (int i=0; i<atomList.length; i++)
							encodeBits(atomList[i], 8);
						}
					}
				}

			addBondQueryFeatures(9, false, nbits,
								 Molecule.cBondQFRingState,
								 Molecule.cBondQFRingStateBits,
								 Molecule.cBondQFRingStateShift);

			addBondQueryFeatures(10, false, nbits,
								 Molecule.cBondQFBondTypes,
								 Molecule.cBondQFBondTypesBits,
								 Molecule.cBondQFBondTypesShift);

			addAtomQueryFeatures(11, false, nbits, Molecule.cAtomQFMatchStereo, 1, -1);

			addBondQueryFeatures(12, false, nbits,
								 Molecule.cBondQFBridge,
								 Molecule.cBondQFBridgeBits,
								 Molecule.cBondQFBridgeShift);

			addAtomQueryFeatures(13, false, nbits,
								 Molecule.cAtomQFPiElectrons,
								 Molecule.cAtomQFPiElectronBits,
								 Molecule.cAtomQFPiElectronShift);

			addAtomQueryFeatures(14, false, nbits,
								 Molecule.cAtomQFNeighbours,
								 Molecule.cAtomQFNeighbourBits,
								 Molecule.cAtomQFNeighbourShift);

			isSecondFeatureBlock |= addAtomQueryFeatures(16, isSecondFeatureBlock, nbits,
														 Molecule.cAtomQFRingSize,
														 Molecule.cAtomQFRingSizeBits,
														 Molecule.cAtomQFRingSizeShift);
			}

		count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mAbnormalValence != null && mAbnormalValence[mGraphAtom[atom]] != -1)
				count++;
		if (count != 0) {
			isSecondFeatureBlock = ensureSecondFeatureBlock(isSecondFeatureBlock);
			encodeBits(1, 1);   //  more data to come
			encodeBits(1, 4);   //  (17-offset) 17 = datatype 'AtomAbnormalValence'
			encodeBits(count, nbits);
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mAbnormalValence != null && mAbnormalValence[mGraphAtom[atom]] != -1) {
					encodeBits(atom, nbits);
					encodeBits(mAbnormalValence[mGraphAtom[atom]], 4);
					}
				}
			}

		if ((mMode & ENCODE_ATOM_CUSTOM_LABELS) != 0) {
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
				isSecondFeatureBlock = ensureSecondFeatureBlock(isSecondFeatureBlock);
				int lbits = idGetNeededBits(maxLength);
				encodeBits(1, 1);   //  more data to come
				encodeBits(2, 4);   //  (18-offset) 18 = datatype 'AtomCustomLabel'
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
			isSecondFeatureBlock |= addAtomQueryFeatures(19, isSecondFeatureBlock, nbits,
					 									 Molecule.cAtomQFCharge,
					 									 Molecule.cAtomQFChargeBits,
					 									 Molecule.cAtomQFChargeShift);

			isSecondFeatureBlock |= addBondQueryFeatures(20, isSecondFeatureBlock, nbits,
														 Molecule.cBondQFRingSize,
														 Molecule.cBondQFRingSizeBits,
														 Molecule.cBondQFRingSizeShift);
			}

		count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mMol.getAtomRadical(mGraphAtom[atom]) != 0)
				count++;
		if (count != 0) {
			isSecondFeatureBlock = ensureSecondFeatureBlock(isSecondFeatureBlock);
			encodeBits(1, 1);   //  more data to come
			encodeBits(5, 4);   //  (21-offset) 21 = datatype 'AtomRadicalState'
			encodeBits(count, nbits);
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mMol.getAtomRadical(mGraphAtom[atom]) != 0) {
					encodeBits(atom, nbits);
					encodeBits(mMol.getAtomRadical(mGraphAtom[atom]) >> Molecule.cAtomRadicalStateShift, 2);
					}
				}
			}

		if (mMol.isFragment()) {	// more QueryFeatures and fragment specific properties
			isSecondFeatureBlock |= addAtomQueryFeatures(22, isSecondFeatureBlock, nbits, Molecule.cAtomQFFlatNitrogen, 1, -1);
			isSecondFeatureBlock |= addBondQueryFeatures(23, isSecondFeatureBlock, nbits, Molecule.cBondQFMatchStereo, 1, -1);
			isSecondFeatureBlock |= addBondQueryFeatures(24, isSecondFeatureBlock, nbits,
					 									 Molecule.cBondQFAromState,
					 									 Molecule.cBondQFAromStateBits,
					 									 Molecule.cBondQFAromStateShift);
			}

		if ((mMode & ENCODE_ATOM_SELECTION) != 0) {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (mMol.isSelectedAtom(mGraphAtom[atom])) {
					isSecondFeatureBlock = ensureSecondFeatureBlock(isSecondFeatureBlock);
					encodeBits(1, 1);   //  more data to come
					encodeBits(9, 4);   //  (25-offset) 25 = datatype 'AtomSelection'
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

			isSecondFeatureBlock = ensureSecondFeatureBlock(isSecondFeatureBlock);
			encodeBits(1, 1);   //  more data to come
			encodeBits(10, 4);   //  (26-offset) 26 = datatype 'delocalized high order bond'
			encodeBits(count, nbits);
			for (int bond=0; bond<mMol.getBonds(); bond++)
				if (isAromaticSPBond[mGraphBond[bond]])
					encodeBits(bond, nbits);
			}

		encodeBits(0, 1);
		mIDCode = encodeBitsEnd();
		}


	private boolean ensureSecondFeatureBlock(boolean isSecondFeatureBlock) {
		if (!isSecondFeatureBlock) {
			encodeBits(1, 1);   //  more data to come
			encodeBits(15, 4);   //  15 = datatype 'start second query feature set'
			}
		return true;
		}


	private boolean addAtomQueryFeatures(int codeNo, boolean isSecondFeatureBlock, int nbits,
										 int qfMask, int qfBits, int qfShift) {
		int count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if ((mMol.getAtomQueryFeatures(mGraphAtom[atom]) & qfMask) != 0)
				count++;

		if (count == 0)
			return false;

		if (codeNo > 15) {
			ensureSecondFeatureBlock(isSecondFeatureBlock);
			codeNo -= 16;
			}

		encodeBits(1, 1);		   //  more data to come
		encodeBits(codeNo, 4);	  //  datatype
		encodeBits(count, nbits);
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			int feature = mMol.getAtomQueryFeatures(mGraphAtom[atom]) & qfMask;
			if (feature != 0) {
				encodeBits(atom, nbits);
				if (qfBits != 1)
					encodeBits(feature >> qfShift, qfBits);
				}
			}

		return true;
		}


	private boolean addBondQueryFeatures(int codeNo, boolean isSecondFeatureBlock, int nbits, int qfMask, int qfBits, int qfShift) {
		int count = 0;
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if ((mMol.getBondQueryFeatures(mGraphBond[bond]) & qfMask) != 0)
				count++;

		if (count == 0)
			return false;

		if (codeNo > 15) {
			ensureSecondFeatureBlock(isSecondFeatureBlock);
			codeNo -= 16;
			}

		encodeBits(1, 1);		   //  more data to come
		encodeBits(codeNo, 4);	  //  datatype
		encodeBits(count, nbits);
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			int feature = mMol.getBondQueryFeatures(mGraphBond[bond]) & qfMask;
			if (feature != 0) {
				encodeBits(bond, nbits);
				if (qfBits != 1)
					encodeBits(feature >> qfShift, qfBits);
				}
			}

		return true;
		}


	private boolean[] getAromaticSPBonds() {
		boolean[] isAromaticSPBond = null;
		RingCollection ringSet = mMol.getRingSet();
		for (int r=0; r<ringSet.getSize(); r++) {
			if (ringSet.isDelocalized(r)) {
				int[] ringAtom = ringSet.getRingAtoms(r);
				int count = 0;
				for (int atom:ringAtom)
					if (mMol.getAtomPi(atom) == 2)
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
						while (mMol.getAtomPi(ringAtom[index]) == 2)
							index++;
						while (mMol.getAtomPi(ringAtom[index]) != 2)
							index = validateCyclicIndex(index+1, ringAtom.length);
						while (count > 0) {
							isAromaticSPBond[ringBond[index]] = true;
							index = validateCyclicIndex(index+2, ringAtom.length);
							count -= 2;
							while (mMol.getAtomPi(ringAtom[index]) != 2)
								index = validateCyclicIndex(index+1, ringAtom.length);
							}
						}
					}
				}
			}
		return isAromaticSPBond;
		}


	private int validateCyclicIndex(int index, int limit) {
		return (index < limit) ? index : index - limit;
		}


	public String getEncodedCoordinates() {
		return getEncodedCoordinates(mZCoordinatesAvailable);
		}


	public String getEncodedCoordinates(boolean keepAbsoluteValues) {
		if (mCoordinates == null) {
			generateGraph();
			encodeCoordinates(keepAbsoluteValues);
			}

		return mCoordinates;
		}


/*	private void encodeCoordinates(boolean keepAbsoluteValues) {
		if (mMol.getAtoms() == 0) {
			mCoordinates = "";
			return;
			}

		float maxDelta = 0.0f;
		for (int i=1; i<mMol.getAtoms(); i++) {
			int atom = mGraphAtom[i];
			int from = (mGraphFrom[i] == -1) ? -1 : mGraphAtom[mGraphFrom[i]];

			float deltaX = (from == -1) ?
							Math.abs(mMol.getAtomX(atom) - mMol.getAtomX(mGraphAtom[0])) / 8.0f
						  : Math.abs(mMol.getAtomX(atom) - mMol.getAtomX(from));
			if (maxDelta < deltaX)
				maxDelta = deltaX;

			float deltaY = (from == -1) ?
							Math.abs(mMol.getAtomY(atom) - mMol.getAtomY(mGraphAtom[0])) / 8.0f
						  : Math.abs(mMol.getAtomY(atom) - mMol.getAtomY(from));
			if (maxDelta < deltaY)
				maxDelta = deltaY;

			if (mZCoordinatesAvailable) {
				float deltaZ = (from == -1) ?
								Math.abs(mMol.getAtomZ(atom) - mMol.getAtomZ(mGraphAtom[0])) / 8.0f
							  : Math.abs(mMol.getAtomZ(atom) - mMol.getAtomZ(from));
				if (maxDelta < deltaZ)
					maxDelta = deltaZ;
				}
			}

		if (maxDelta == 0.0) {
			mCoordinates = "";
			return;
			}

		float increment = maxDelta / 43.0f;
		float halfIncrement = increment / 2.0f;

		StringBuilder coordinateBuffer = new StringBuilder();

		for (int i=1; i<mMol.getAtoms(); i++) {
			int atom = mGraphAtom[i];
			int from = (mGraphFrom[i] == -1) ? -1 : mGraphAtom[mGraphFrom[i]];

			float deltaX = (from == -1) ?
							(mMol.getAtomX(atom) - mMol.getAtomX(mGraphAtom[0])) / 8.0f
						   : mMol.getAtomX(atom) - mMol.getAtomX(from);

							float deltaY = (from == -1) ?
							(mMol.getAtomY(atom) - mMol.getAtomY(mGraphAtom[0])) / 8.0f
						   : mMol.getAtomY(atom) - mMol.getAtomY(from);

			coordinateBuffer.append((char)(40 + (int)((maxDelta + deltaX + halfIncrement) / increment)));
			coordinateBuffer.append((char)(40 + (int)((maxDelta + deltaY + halfIncrement) / increment)));
			}

		if (mZCoordinatesAvailable) {
			for (int i=1; i<mMol.getAtoms(); i++) {
				int atom = mGraphAtom[i];
				int from = (mGraphFrom[i] == -1) ? -1 : mGraphAtom[mGraphFrom[i]];

				float deltaZ = (from == -1) ?
								(mMol.getAtomZ(atom) - mMol.getAtomZ(mGraphAtom[0])) / 8.0f
							   : mMol.getAtomZ(atom) - mMol.getAtomZ(from);

				coordinateBuffer.append((char)(40 + (int)((maxDelta + deltaZ + halfIncrement) / increment)));
				}
			}

		if (keepAbsoluteValues) {
			coordinateBuffer.append('&');	// old faulty encoding started with "'"

			int avblInt = encodeABVL(mMol.getAverageBondLength(), 7396);
			coordinateBuffer.append((char)(40 + avblInt/86));
			coordinateBuffer.append((char)(40 + avblInt%86));

			int xInt = encodeShift(mMol.getAtomX(mGraphAtom[0]), 7396);
			coordinateBuffer.append((char)(40 + xInt/86));
			coordinateBuffer.append((char)(40 + xInt%86));

			int yInt = encodeShift(mMol.getAtomY(mGraphAtom[0]), 7396);
			coordinateBuffer.append((char)(40 + yInt/86));
			coordinateBuffer.append((char)(40 + yInt%86));
			if (mZCoordinatesAvailable) {
				int zInt = encodeShift(mMol.getAtomZ(mGraphAtom[0]), 7396);
				coordinateBuffer.append((char)(40 + zInt/86));
				coordinateBuffer.append((char)(40 + zInt%86));
				}
			}
		
		mCoordinates = coordinateBuffer.toString();
		}	*/

	private void encodeCoordinates(boolean keepAbsoluteValues) {
		if (mMol.getAtoms() == 0) {
			mCoordinates = "";
			return;
			}

		int resolutionBits = mZCoordinatesAvailable ? 16 : 8;	// must be an even number
		encodeBitsStart();
		mEncodingBuffer.append('!');
		encodeBits(mZCoordinatesAvailable ? 1 : 0, 1);
		encodeBits(keepAbsoluteValues ? 1 : 0, 1);
		encodeBits(resolutionBits/2, 4);	// resolution bits devided by 2

		float maxDelta = 0.0f;
		for (int i=1; i<mMol.getAtoms(); i++) {
			int atom = mGraphAtom[i];
			int from = (mGraphFrom[i] == -1) ? -1 : mGraphAtom[mGraphFrom[i]];

			float deltaX = (from == -1) ?
							Math.abs(mMol.getAtomX(atom) - mMol.getAtomX(mGraphAtom[0])) / 8.0f
						  : Math.abs(mMol.getAtomX(atom) - mMol.getAtomX(from));
			if (maxDelta < deltaX)
				maxDelta = deltaX;

			float deltaY = (from == -1) ?
							Math.abs(mMol.getAtomY(atom) - mMol.getAtomY(mGraphAtom[0])) / 8.0f
						  : Math.abs(mMol.getAtomY(atom) - mMol.getAtomY(from));
			if (maxDelta < deltaY)
				maxDelta = deltaY;

			if (mZCoordinatesAvailable) {
				float deltaZ = (from == -1) ?
								Math.abs(mMol.getAtomZ(atom) - mMol.getAtomZ(mGraphAtom[0])) / 8.0f
							  : Math.abs(mMol.getAtomZ(atom) - mMol.getAtomZ(from));
				if (maxDelta < deltaZ)
					maxDelta = deltaZ;
				}
			}

		if (maxDelta == 0.0) {
			mCoordinates = "";
			return;
			}

		int binCount = (1 << resolutionBits);
		float increment = maxDelta / (binCount / 2.0f - 1);
		float halfIncrement = increment / 2.0f;

		for (int i=1; i<mMol.getAtoms(); i++) {
			int atom = mGraphAtom[i];
			int from = (mGraphFrom[i] == -1) ? -1 : mGraphAtom[mGraphFrom[i]];

			float deltaX = (from == -1) ?
							(mMol.getAtomX(atom) - mMol.getAtomX(mGraphAtom[0])) / 8.0f
						   : mMol.getAtomX(atom) - mMol.getAtomX(from);

			float deltaY = (from == -1) ?
							(mMol.getAtomY(atom) - mMol.getAtomY(mGraphAtom[0])) / 8.0f
						   : mMol.getAtomY(atom) - mMol.getAtomY(from);

			encodeBits((int)((maxDelta + deltaX + halfIncrement) / increment), resolutionBits);
			encodeBits((int)((maxDelta + deltaY + halfIncrement) / increment), resolutionBits);

			if (mZCoordinatesAvailable) {
				float deltaZ = (from == -1) ?
								(mMol.getAtomZ(atom) - mMol.getAtomZ(mGraphAtom[0])) / 8.0f
							   : mMol.getAtomZ(atom) - mMol.getAtomZ(from);

				encodeBits((int)((maxDelta + deltaZ + halfIncrement) / increment), resolutionBits);
				}
			}

		if (keepAbsoluteValues) {
			encodeBits(encodeABVL(mMol.getAverageBondLength(true), binCount), resolutionBits);

			encodeBits(encodeShift(mMol.getAtomX(mGraphAtom[0]), binCount), resolutionBits);
			encodeBits(encodeShift(mMol.getAtomY(mGraphAtom[0]), binCount), resolutionBits);

			if (mZCoordinatesAvailable)
				encodeBits(encodeShift(mMol.getAtomZ(mGraphAtom[0]), binCount), resolutionBits);
			}

		mCoordinates = encodeBitsEnd();
		}

	/**
	 * Encode a floating point value into an integer with precision proportional to the value itself.
	 * @param value
	 * @return
	 */
	private int encodeABVL(float value, int binCount) {
		return Math.min(binCount-1, Math.max(0, (int)(0.5 + Math.log10(value/0.1) / Math.log10(200/0.1) * (binCount-1))));
		}

	private int encodeShift(float value, int binCount) {
		int halfBinCount = binCount / 2;
		boolean isNegative =  (value < 0);
		value = Math.abs(value);
		float steepness = (float)binCount/100f;
		int intValue = (int)(0.5 + value * (halfBinCount-1) / (value + steepness));
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

		int nbits = idGetNeededBits(maxMapNo);
		encodeBitsStart();
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


	private void idGenerateConfigurations() {
		// Creates parities based on atom indices in graph rather than on priority values.
		// These values are more meaningful to be written into idcodes, because they allow
		// to create coordinates or running Configuration aware substructure searches on
		// molecules creates from idcode without the necessity to recreate the priority values.

		mTHConfiguration = new byte[mMol.getAtoms()];

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mTHParity[atom] == Molecule.cAtomParity1
			 || mTHParity[atom] == Molecule.cAtomParity2) {
				boolean inversion = false;
				if (isCentralAlleneAtom(atom)) {
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


/*	private void idNormalizeConfigurations() {
		// Atom TH-parities of type ABS in meso fragments can be defined in
		// two degenerate ways. Normalize them based on the current mCanRank
		if (mMesoHelper != null) {
			mMesoHelper.normalizeFragmentsAbsAtoms(mCanRank, mTHConfiguration);
			}

		// Atom TH-parities of ESR-AND or ESR-OR groups are up to now
		// arbitrary values, i.e. one of two possible ways of encoding
		// the group's relative parity information.
		// First we create for every ESR group a list of its atoms.
		// All parities of a list's atoms are then normalized by inverting
		// all parities of a group if the parity of the highest ranking
		// group member is parity2.

		int count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (canIsMemberOfESRGroup(atom))
				count++;

		if (count == 0)
			return;

		int[] parity = new int[count];
		count = 0;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (canIsMemberOfESRGroup(atom)) {
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
				for (int i=groupBase; i<nextGroupBase; i++)
					atomList[i-groupBase] = parity[i] & 0x00000fff;
				idNormalizeConfigurations(atomList);

				if (nextGroupBase == parity.length)
					break;

				groupID = (parity[nextGroupBase] & 0xff000000);
				groupBase = nextGroupBase;
				}
			}
		}


	private void idNormalizeConfigurations(int[] atomList) {
			// atomList is sorted by mCanRank with the lowest rank first
		if (mTHParity[atomList[atomList.length-1]] == Molecule.cAtomParity2) {
			for (int i=0; i<atomList.length; i++) {
				int atom = atomList[i];
				if (mTHConfiguration[atom] == Molecule.cAtomParity1)
					mTHConfiguration[atom] = Molecule.cAtomParity2;
				else if (mTHConfiguration[atom] == Molecule.cAtomParity2)
					mTHConfiguration[atom] = Molecule.cAtomParity1;
				}
			}
		}*/

	
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


	private void encodeBitsStart() {
		mEncodingBuffer = new StringBuilder();
		mEncodingBitsAvail = 6;
		mEncodingTempData = 0;
		}


	private void encodeBits(int data, int bits) {
//System.out.println(bits+" bits:"+data+"  mode="+mode);
		while (bits != 0) {
			if (mEncodingBitsAvail == 0) {
				mEncodingBuffer.append((char)(mEncodingTempData + 64));
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
		mEncodingBuffer.append((char)(mEncodingTempData + 64));
		return mEncodingBuffer.toString();
		}


	private int idGetNeededBits(int no) {
		int bits = 0;
		while (no > 0) {
			no >>= 1;
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


	protected void setParities() {
		// Creates parities based on atom indices of original molecule and
		// stores them in molecule. It also set the stereo center flag.
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mTHParity[atom] == Molecule.cAtomParity1
			 || mTHParity[atom] == Molecule.cAtomParity2) {
				boolean inversion = false;
				if (mMol.getAtomPi(atom) != 0
				 && mMol.getConnAtoms(atom) == 2
				 && mMol.getConnBondOrder(atom,0) == 2
				 && mMol.getConnBondOrder(atom,1) == 2) {   // allene parities
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
							int connAtom[] = new int[2];
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
						int connAtom[] = new int[2];
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
		int orderedConn[] = new int[noOfConns];
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

		int graphAtom[] = new int[graphSize];
		int graphParent[] = new int[graphSize];
		int graphRank[] = new int[graphSize];
		boolean graphIsPseudo[] = new boolean[graphSize];

		boolean atomUsed[] = new boolean[mMol.getAllAtoms()];

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
			cipUpdateParentRanking(graphIsPseudo, graphRank, graphParent,
								   graphAtom, levelStart, currentLevel);
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

	private boolean cipTryDistinguishBranches(boolean graphIsPseudo[],
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

	private void cipUpdateParentRanking(boolean graphIsPseudo[],
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
												: mMol.getAllHydrogens(graphAtom[parentIndex]);
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

	protected EZHalfParity(ExtendedMolecule mol, int rank[], int atom1, int atom2) {
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
				continue;
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

		float angleDB = mMol.getBondAngle(mCentralAxialAtom,mRemoteAxialAtom);
		float angleHigh = mMol.getBondAngle(mCentralAxialAtom,mHighConn);
		if (angleHigh < angleDB)
			angleHigh += Math.PI*2;

		if (mMol.getAllConnAtoms(mCentralAxialAtom) == 2) {
			float angleDif = angleHigh - angleDB;
			if ((angleDif > Math.PI - 0.05) && (angleDif < Math.PI + 0.05)) {
				mValue = -1;	// less than 3 degrees different from double bond
				return mValue;	// is counted as non-stereo-specified double bond
				}
			mValue = (angleDif < Math.PI) ? 4 : 2;
			return mValue;
			}
		else {
			float angleOther = mMol.getBondAngle(mCentralAxialAtom,mLowConn);
			if (angleOther < angleDB)
				angleOther += Math.PI*2;
			mValue = (angleOther < angleHigh) ? 2 : 4;
			return mValue;
			}
		}
	}


class CanonizerFragment {
	int		atom[];
	int		bond[];

	protected CanonizerFragment(int[] atom, int atoms, int[] bond, int bonds) {
		this.atom = new int[atoms];
		this.bond = new int[bonds];
		for (int a=0; a<atoms; a++)
			this.atom[a] = atom[a];
		for (int b=0; b<bonds; b++)
			this.bond[b] = bond[b];
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
