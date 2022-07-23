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
 * 3. Neither the name of the the copyright holder nor the
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


// This class handles the meso detection for the Canonizer class.

package com.actelion.research.chem;

import java.util.*;

public class CanonizerMesoHelper {
	private static final int REMOVE_ESR_GROUP = 1;
	private static final int SWAP_ESR_GROUPS = 2;

	private ExtendedMolecule mMol;
	private int[] mCanRankWithoutStereo;
	private byte[] mTHParity;
	private byte[] mEZParity;
	private byte[] mTHESRType;
	private byte[] mTHESRGroup;
//	private byte[] mEZESRType;
//	private byte[] mEZESRGroup;
	private int[][] mMesoFragmentAtom;
	private boolean[] mIsStereoCenter;  // based on meso ranking
	private boolean[] mIsMesoFragmentMember;
	private boolean[] mTHParityRoundIsOdd;
	private boolean[] mEZParityRoundIsOdd;
	private boolean[] mTHESRTypeNeedsNormalization;
	private ArrayList<ESRGroupNormalizationInfo> mESRGroupNormalizationInfoList;

	protected CanonizerMesoHelper(ExtendedMolecule mol,
								  int[] canRankWithoutStereo,
								  boolean[] isStereoCenter,
								  byte[] thParity,
								  byte[] ezParity,
								  byte[] thESRType,
								  byte[] thESRGroup,
								  byte[] ezESRType,
								  byte[] ezESRGroup,
								  boolean[] thParityRoundIsOdd,
								  boolean[] ezParityRoundIsOdd,
								  boolean[] esrTypeNeedsNormalization) {
		mMol = mol;
		mCanRankWithoutStereo = canRankWithoutStereo;
		mIsStereoCenter = isStereoCenter;
		mTHParity = thParity;
		mEZParity = ezParity;
		mTHESRType = thESRType;
		mTHESRGroup = thESRGroup;
//		mEZESRType = ezESRType;
//		mEZESRGroup = ezESRGroup;
		mTHParityRoundIsOdd = thParityRoundIsOdd;
		mEZParityRoundIsOdd = ezParityRoundIsOdd;
		mTHESRTypeNeedsNormalization = esrTypeNeedsNormalization;

		findMesoFragments();
		}

	protected boolean isMeso() {
		// Must be called with mTHParities and mEZParities
		// freshly calculated after ESR pre-normalization.

		boolean meso = true;
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
//System.out.println("isMeso() atom:"+atom+" parity:"+mTHParity[atom]+" mesoFragment:"+mMesoFragmentNo[atom]+" atomicNo:"+mMol.getAtomicNo(atom));
			if (mTHParity[atom] != Molecule.cAtomParityNone
			 && !mIsMesoFragmentMember[atom]) {
				meso = false;
				break;
				}
			}
		return meso;
		}


	protected boolean isInMesoFragment(int atom) {
		return mIsMesoFragmentMember[atom];
		}


	private boolean mayBeMirrorAtoms(int atom1, int atom2) {
		if (atom1 == atom2)
			return false;

			// a necessary condition is equal symmetry rank
		if (mCanRankWithoutStereo[atom1] != mCanRankWithoutStereo[atom2])
			return false;

			// if our atoms are stereo centers ...
		if (mTHParity[atom1] != Molecule.cAtomParityNone) {
				// ... they must have known configuration
			if (mTHParity[atom1] == Molecule.cAtomParityUnknown
			 || mTHParity[atom2] == Molecule.cAtomParityUnknown)
				return false;

				// ... they must have inverse configuration
				// (or matching configuration if the parity detection round was even)
			if (mTHParityRoundIsOdd[atom1]
			 ^ (mTHParity[atom1] != mTHParity[atom2]))
				return false;

				// ... they must have matching ESR type and group
			if (mTHESRType[atom1] != mTHESRType[atom2]
			 || mTHESRGroup[atom1] != mTHESRGroup[atom2])
				return false;
			}

			// if we have a direct BINAP bond or E-double bond we refuse
		int bond = mMol.getBond(atom1, atom2);
		if (bond != -1) {
			if (mMol.getBondOrder(bond) == 1
			 && mEZParity[bond] != Molecule.cBondParityNone)
				return false;
			if (mMol.getBondOrder(bond) == 2
			 && mEZParity[bond] == Molecule.cBondParityEor1)
				return false;
			}

			// if our atoms are part of a non-aromatic double bond ...
		if (mMol.getAtomPi(atom1) == 1 && !mMol.isAromaticAtom(atom1)) {
			int bond1 = -1;
			for (int i=0; i<mMol.getConnAtoms(atom1); i++) {
				if (mMol.getConnAtom(atom1, i) != atom2
				 && mMol.getConnBondOrder(atom1, i) == 2) {
					bond1 = mMol.getConnBond(atom1, i);
					break;
					}
				}
			int bond2 = -1;
			for (int i=0; i<mMol.getConnAtoms(atom2); i++) {
				if (mMol.getConnAtom(atom2, i) != atom1
				 && mMol.getConnBondOrder(atom2, i) == 2) {
					bond2 = mMol.getConnBond(atom2, i);
					break;
					}
				}

			if (bond1 != -1 // for symmetry reasons bond2 must not be -1 then
			 && mEZParity[bond1] != Molecule.cBondParityNone
			 && (mEZParityRoundIsOdd[bond1] ^ (mEZParity[bond1] == mEZParity[bond2])))
				return false;
			}

		return true;
		}


	private void findMesoFragments() {
		TreeSet<MesoFragmentMembers> mesoFragmentList = new TreeSet<>();

		// Detect mirror planes by finding a seed atom (an atom with
		// at least 2 neighbours sharing the same canRankWithoutStereo)
		// or seed bond (a bond connecting to atoms sharing the same
		// canRankWithoutStereo).
		for (int seedAtom=0; seedAtom<mMol.getAtoms(); seedAtom++) {
			if (mMol.getAtomPi(seedAtom) < 2  // exclude allenes but include sulfones
			 || mMol.getConnAtoms(seedAtom) > 2) {
				for (int i=1; i<mMol.getConnAtoms(seedAtom); i++) {
					int atom1 = mMol.getConnAtom(seedAtom, i);
					for (int j=0; j<i; j++) {
						int atom2 = mMol.getConnAtom(seedAtom, j);
						if (mayBeMirrorAtoms(atom1, atom2))
							tryAddNewMesoFragment(atom1, atom2, mesoFragmentList);
						}
					}
				}
			}
		for (int seedBond=0; seedBond<mMol.getBonds(); seedBond++) {
			if (mEZParity[seedBond] != Molecule.cBondParityNone) {
				// don't consider single bonds with axial chirality or double bonds with unknown or E-configuration
				if (mMol.getBondOrder(seedBond) != 2
				 || mEZParity[seedBond] != Molecule.cBondParityZor2)
					continue;
				}
			
			int atom1 = mMol.getBondAtom(0, seedBond);
			int atom2 = mMol.getBondAtom(1, seedBond);
			if (mayBeMirrorAtoms(atom1, atom2))
				tryAddNewMesoFragment(atom1, atom2, mesoFragmentList);
			}

		mMesoFragmentAtom = new int[mesoFragmentList.size()][];
		mIsMesoFragmentMember = new boolean[mMol.getAtoms()];
		int fragmentNo = 0;
		for (MesoFragmentMembers members:mesoFragmentList) {
			mMesoFragmentAtom[fragmentNo++] = members.memberAtom;
			for (int i=0; i<members.memberAtom.length; i++)
				mIsMesoFragmentMember[members.memberAtom[i]] = true;
			}
/*
System.out.println("--mesofragments found:--------------------------------------------");
for (int i=0; i<mMesoFragmentAtom.length; i++) {
 System.out.print("fragment["+i+"]:");
 for (int j=0; j<mMesoFragmentAtom[i].length; j++) {
  int atom = mMesoFragmentAtom[i][j];
  System.out.print(" "+atom);
  if (mTHParity[atom] != 0)
   System.out.print(",p:"+mTHParity[atom]+",o:"+mTHParityRoundIsOdd[atom]);
  }
 System.out.println();
 }
*/
		}


	private void tryAddNewMesoFragment(int atom1, int atom2, TreeSet<MesoFragmentMembers> mesoFragmentList) {
		MesoFragmentMembers members = tryFindMesoFragment(atom1, atom2);
		if (members != null && members.hasStereoCenters(mIsStereoCenter))
			mesoFragmentList.add(members);
		}


	/**
	 * Tries to find a symmetry plane starting from symmetrical atom1 and atom2.
	 * Stops atom matching at non-ring single bonds, thus locates and assigns one meso fragment.		
	 * @param atom1
	 * @param atom2
	 * @return a MesoFragmentMembers object, if a meso fragment was found
	 */
	private MesoFragmentMembers tryFindMesoFragment(int atom1, int atom2) {
		int[] graphAtom = new int[mMol.getAtoms()];
		int[] matchAtom = new int[mMol.getAtoms()];
		boolean[] isOrthogonal = new boolean[mMol.getAtoms()];
		boolean[] hasOrthogonality = new boolean[mMol.getAtoms()];
		MesoFragmentBranch[] branch = new MesoFragmentBranch[mMol.getAtoms()];
		MesoFragmentMembers members = new MesoFragmentMembers(mMol.getAtoms());

		graphAtom[0] = atom1;
		matchAtom[atom1] = atom2;
		matchAtom[atom2] = -2;  // -2 := on mirror side
		members.add(atom1);
		members.add(atom2);

		int current = 0;
		int highest = 0;
		while (current <= highest) {
			int currentAtom = graphAtom[current];

				// if currentAtom is in mirror plane
			if (matchAtom[currentAtom] == currentAtom) {
				for (int i=0; i<mMol.getConnAtoms(currentAtom); i++) {
					int candidate = mMol.getConnAtom(currentAtom, i);
					if (!members.isMember[candidate]) {

						// if candidate is a period 2 element (C,N,O rather than Si,P,S)
						// and if candidateBond is double-bond then these belong
						// to the fragment. The candidate's orthogonality is inverted
						// with respect to the currentAtom.
						if (mMol.getConnBondOrder(currentAtom, i) == 2
						 && mMol.getAtomicNo(candidate) < 10) {
							graphAtom[++highest] = candidate;
							matchAtom[candidate] = candidate;
							hasOrthogonality[candidate] = hasOrthogonality[currentAtom] || (mMol.getAtomPi(candidate) == 2);
							isOrthogonal[candidate] = hasOrthogonality[currentAtom] ?
													  !isOrthogonal[currentAtom] : false;
							members.add(candidate);
							}

						else if (hasOrthogonality[currentAtom]
							  && isOrthogonal[currentAtom]) {
							int opponent = findMirrorAtom(candidate, matchAtom[currentAtom], members.isMember);
							if (opponent == -1)
								return null;
							
							graphAtom[++highest] = candidate;
							matchAtom[candidate] = opponent;
							matchAtom[opponent] = -2;
							hasOrthogonality[candidate] = false;
							members.add(candidate);
							members.add(opponent);
							}

						else if (mMol.isRingBond(mMol.getConnBond(currentAtom, i))) {
							graphAtom[++highest] = candidate;
							matchAtom[candidate] = candidate;
							hasOrthogonality[candidate] = false;
							members.add(candidate);

							// Tetrahedral atoms with more than two neighbours
							// must have two symmetrical atoms sticking out of
							// the mirror plane in opposite directions.
							if (isTetrahedral(candidate)
							 && mMol.getConnAtoms(candidate) > 2) {
								boolean found = false;
								for (int j=1; j<mMol.getConnAtoms(candidate); j++) {
									int symAtom1 = mMol.getConnAtom(candidate, j);
									if (!members.isMember[symAtom1]) {
										for (int k=0; k<j; k++) {
											int symAtom2 = mMol.getConnAtom(candidate, k);
											if (!members.isMember[symAtom2]) {
												if (mayBeMirrorAtoms(symAtom1, symAtom2)) {
													graphAtom[++highest] = symAtom1;
													matchAtom[symAtom1] = symAtom2;
													matchAtom[symAtom2] = -2;
													hasOrthogonality[symAtom1] = false;
													members.add(symAtom1);
													members.add(symAtom2);
													found = true;
													}
												}
											}
										}
									}
								if (!found)
									return null;
								}
							}
						}
					}
				}
			else {  // for current atoms that are not in the mirror plane

				// mark those neighbours of the currentAtom which are in the symmetry plane
				boolean[] connIsOnMirrorPlane = new boolean[mMol.getConnAtoms(currentAtom)];
				for (int i=0; i<mMol.getConnAtoms(currentAtom); i++) {
					int candidate = mMol.getConnAtom(currentAtom, i);
					if (members.isMember[candidate]) {
						connIsOnMirrorPlane[i] = (matchAtom[candidate] == candidate);
						}
					else {
						for (int j=0; j<mMol.getConnAtoms(candidate); j++) {
							if (mMol.getConnAtom(candidate, j) == matchAtom[currentAtom]) {
								connIsOnMirrorPlane[i] = true;
								break;
								}
							}
						}
					}

					// add first those neighbours which are in the mirror plane
				for (int i=0; i<mMol.getConnAtoms(currentAtom); i++) {
					if (connIsOnMirrorPlane[i]) {
						int candidate = mMol.getConnAtom(currentAtom, i);
						if (members.isMember[candidate]) {
							// If we have a ring closure to an atom on the mirror plane, then
							// make sure the current atom's match atom also connects to the candidate
							if (mMol.getBond(candidate, matchAtom[currentAtom]) == -1)
								return null;
							}
						else {
							graphAtom[++highest] = candidate;
							matchAtom[candidate] = candidate;
							isOrthogonal[candidate] = false;
							hasOrthogonality[candidate] = true;
							members.add(candidate);
							}
						}
					}

					// add now those neighbours and opponents not being in the mirror plane
				MesoFragmentBranch b = branch[currentAtom];
				for (int i=(b == null) ? 0 : b.neighbourIndex; i<mMol.getConnAtoms(currentAtom); i++) {
					if (!connIsOnMirrorPlane[i]) {
						int candidate = mMol.getConnAtom(currentAtom, i);
						if (!members.isMember[candidate]) {
							int opponent = findMirrorAtom(candidate, matchAtom[currentAtom], members.isMember);
							if (opponent == -1)
								return null;

							graphAtom[++highest] = candidate;
							matchAtom[candidate] = opponent;
							matchAtom[opponent] = -2;
							hasOrthogonality[candidate] = false;
							members.add(candidate);
							members.add(opponent);
							}
						}
					}

				}
			current++;
			}

//printFragmentMsg("exit3", atom1, atom2, current, highest, graphAtom, matchAtom, isOrthogonal, hasOrthogonality);

/*
System.out.println("--mesofragment found----atom1:"+atom1+"--atom2:"+atom2+"------------------------------------");
for (int atom=0; atom<mMol.getAtoms(); atom++)
if (isFragmentMember[atom]) {
if (matchAtom[atom] == atom)
	System.out.println(" "+atom+",p:"+mTHParity[atom]+",o:"+mTHParityRoundIsOdd[atom]);
else if (matchAtom[atom] != -2)
	System.out.println(" "+atom+",p:"+mTHParity[atom]+",o:"+mTHParityRoundIsOdd[atom]
					  +"<->"+matchAtom[atom]+",p:"+mTHParity[matchAtom[atom]]+",o:"+mTHParityRoundIsOdd[matchAtom[atom]]);
}
System.out.print("graphAtom:");
for (int j=0; j<=highest; j++)
	System.out.print(" "+graphAtom[j]);
System.out.println();
*/
		return members;
		}
/*
private void printFragmentMsg(String msg, int atom1, int atom2, int current, int highest, int[] graphAtom, int[] matchAtom, boolean[] isOrthogonal, boolean[] hasOrthogonality) {
System.out.println("##fragmentMessage:"+msg+" atom1:"+atom1+" atom2:"+atom2+" current:"+current+" highest:"+highest);
System.out.print("graphAtom:");
for (int j=0; j<graphAtom.length; j++)
 System.out.print(" "+graphAtom[j]);
System.out.println();
System.out.print("matchAtom:");
for (int j=0; j<matchAtom.length; j++)
 System.out.print(" "+matchAtom[j]);
System.out.println();
System.out.print("isOrthogonal:");
for (int j=0; j<isOrthogonal.length; j++)
 System.out.print(" "+(isOrthogonal[j] ? 't' : 'f'));
System.out.println();
System.out.print("hasOrthogonality:");
for (int j=0; j<hasOrthogonality.length; j++)
 System.out.print(" "+(hasOrthogonality[j] ? 't' : 'f'));
System.out.println();
}
*/
	private boolean isTetrahedral(int atom) {
		return (mMol.getAtomicNo(atom) == 6 && mMol.getAtomPi(atom) == 0)
			|| (mMol.getAtomicNo(atom) == 7 && mMol.getAtomCharge(atom) == 1)
			|| (mMol.getAtomicNo(atom) == 14)
			|| (mMol.getAtomicNo(atom) == 15 && mMol.getConnAtoms(atom) > 2)
			|| (mMol.getAtomicNo(atom) == 16 && mMol.getConnAtoms(atom) > 2);
		}

	/**
	 * temporarilly
	 * @param atom
	 * @param parentOfMirrorAtom
	 * @param isFragmentMember
	 * @return
	 */
	private int findMirrorAtom(int atom, int parentOfMirrorAtom, boolean[] isFragmentMember) {
		int[] candidate = new int[mMol.getConnAtoms(parentOfMirrorAtom)];
		int index = 0;
		for (int i=0; i<mMol.getConnAtoms(parentOfMirrorAtom); i++) {
			candidate[index] = mMol.getConnAtom(parentOfMirrorAtom, i);
			if (!isFragmentMember[candidate[index]]
			 && mayBeMirrorAtoms(atom, candidate[index]))
				index++;
			}
		if (index == 0)
			return -1;
		if (index == 1)
			return candidate[0];

		// if we have multiple candidates, then take that one that is topologically closer to atom
		int lowCandidate = -1;
		int lowPathLength = Integer.MAX_VALUE;
		for (int i=0; i<index; i++) {
			int pathLength = mMol.getPathLength(atom, candidate[i], Integer.MAX_VALUE, isFragmentMember);
			if (pathLength < lowPathLength) {
				lowPathLength = pathLength;
				lowCandidate = candidate[i];
				}
			}
		return lowCandidate;
		}


	/**
	 * If a meso fragment contains stereo atoms belonging to ESR groups
	 * then there may be an alternative way of specifying the same
	 * meso fragment because of the symmetry of the fragment.
	 *
	 * The procedure to normalize a meso fragment's ESR definition
	 * depends on whether it contains ESR groups that have members
	 * outside of the fragment or not.
	 *
	 * If we have group dependency cycles, i.e. some meso fragments contain
	 * at least two groups each, as f1:g1,g2 f2:g2,g3 f3:g3,g1, then we need
	 * to convert one group of the cycle into ABS atoms.
	 *
	 * To be precise we can determine here which situation to normalize,
	 * however the actual normalization should be postponed until we have
	 * canonization ranks that don't depend on the original grouping.
	 */
	protected void normalizeESRGroups() {
//System.out.println("normalizeESRGroups() mMesoFragmentCount:"+mMesoFragmentAtom.length);
		if (mMesoFragmentAtom != null) {
			ESRGroupFragmentMatrix matrix = new ESRGroupFragmentMatrix();
			mESRGroupNormalizationInfoList = new ArrayList<>();

			for (int fragment=0; fragment<mMesoFragmentAtom.length; fragment++) {
				int dependentGroupCount = matrix.getDependentGroupCount(fragment);
//System.out.println("normalizeESRGroups() dependentGroupCount:"+dependentGroupCount);
				if (dependentGroupCount == 0) {
					// IF we have 1 or more OR groups in the fragment
					//	 IF we have ABS atoms
					//		 1.) create new OR group and put all ABS atoms into it
					//	 ENDIF
					//	 2.) convert all atoms of one of the OR groups into ABS atoms
					//		 (select the OR group in a nomalized way)
					// ELSE IF we have one or more AND groups in the fragment
					//	 IF we have ABS atoms
					//		 1.) create new AND group and put all ABS atoms into it
					//	 ENDIF
					//	 2.) convert all atoms of one of the AND groups into ABS atoms
					//		 (select the AND group in a normalized way)
					// ENDIF
					matrix.cutTiesOfIndependentGroups(fragment);
					int orCount = countESRGroups(fragment, Molecule.cESRTypeOr);
					int andCount = countESRGroups(fragment, Molecule.cESRTypeAnd);
					boolean containsABS = containsTypeABSParity1Or2(fragment);
					if (orCount == 1 && andCount == 1 && !containsABS) {
						putORAtomsIntoANDGroup(fragment, matrix.newESRGroup(Molecule.cESRTypeAnd));

						// after stereo ranking convert lowest ranking group of esrType into ABS atoms
						mESRGroupNormalizationInfoList.add(new ESRGroupNormalizationInfo(fragment,
																REMOVE_ESR_GROUP, -1, -1));
						}
					if (orCount > 0) {
						// put temporarily all fragment's ABS atoms int new group of esrType OR
						if (containsABS) {
							putABSAtomsIntoESRGroup(fragment, matrix.newESRGroup(Molecule.cESRTypeOr), Molecule.cESRTypeOr);
							orCount++;
							}

						// after stereo ranking convert lowest ranking group of esrType into ABS atoms
						mESRGroupNormalizationInfoList.add(new ESRGroupNormalizationInfo(fragment,
																REMOVE_ESR_GROUP, -1, -1));
						}
					else if (andCount > 0) {
						// put temporarily all fragment's ABS atoms int new group of esrType AND
						if (containsABS)
							putABSAtomsIntoESRGroup(fragment, matrix.newESRGroup(Molecule.cESRTypeAnd), Molecule.cESRTypeAnd);

						// after stereo ranking convert lowest ranking group of esrType into ABS atoms
						mESRGroupNormalizationInfoList.add(new ESRGroupNormalizationInfo(fragment,
																REMOVE_ESR_GROUP, -1, -1));
						}
					else if (containsABS) {
						putABSAtomsIntoESRGroup(fragment, matrix.newESRGroup(Molecule.cESRTypeAnd), Molecule.cESRTypeAnd);
						
						// after stereo ranking convert lowest ranking group of esrType into ABS atoms
						mESRGroupNormalizationInfoList.add(new ESRGroupNormalizationInfo(fragment,
																REMOVE_ESR_GROUP, -1, -1));
						}
					}
				else if (dependentGroupCount == 1) {
					// IF we have ABS atoms
					//	 1.) Swapping ESR roles of ABS atoms and dependent group atoms
					//		 doesn't change the molecule. Therefore select one of these
					//		 two representations in a normalized way.
					// ELSE
					//	 1.) Groups are independent from the rest of the molecule, so
					//		 give them independent group numbers.
					//	 2.) We can convert one of the groups into ABS atoms without
					//		 changing the molecule. 
					// ENDIF
					if (containsTypeABSParity1Or2(fragment)) {
						int group = matrix.getDependentGroup(fragment);
						int type = matrix.getDependentType(fragment);
						mESRGroupNormalizationInfoList.add(new ESRGroupNormalizationInfo(fragment,
																SWAP_ESR_GROUPS, group, type));
						}
					else {
						matrix.cutTiesOfIndependentGroups(fragment);
						mESRGroupNormalizationInfoList.add(new ESRGroupNormalizationInfo(fragment,
																REMOVE_ESR_GROUP, -1, -1));
						}
					}
				}
			}
		}


	private boolean containsTypeABSParity1Or2(int fragment) {
		for (int i=0; i<mMesoFragmentAtom[fragment].length; i++) {
			int atom = mMesoFragmentAtom[fragment][i];
			if (hasParity1or2(atom)
			 && mTHESRType[atom] == Molecule.cESRTypeAbs)
				return true;
			}
		
		return false;
		}


	private int countESRGroups(int fragment, int esrType) {
		int count = 0;
		int groupBits = 0;
		for (int i=0; i<mMesoFragmentAtom[fragment].length; i++) {
			int atom = mMesoFragmentAtom[fragment][i];
			if (mTHESRType[atom] == esrType) {
				int groupBit = (1 << mTHESRGroup[atom]);
				if ((groupBits & groupBit) == 0) {
					groupBits |= groupBit;
					count++;
					}
				}
			}
		
		return count;
		}


	private boolean hasParity1or2(int atom) {
		return mIsStereoCenter[atom]
			 && (mTHParity[atom] == Molecule.cAtomParity1
			  || mTHParity[atom] == Molecule.cAtomParity2);
		}


	private void putABSAtomsIntoESRGroup(int fragment, int esrGroup, int esrType) {
		for (int j=0; j<mMesoFragmentAtom[fragment].length; j++) {
			int atom = mMesoFragmentAtom[fragment][j];
			if (hasParity1or2(atom)
			 && mTHESRType[atom] == Molecule.cESRTypeAbs) {
				mTHESRType[atom] = (byte)esrType;
				mTHESRGroup[atom] = (byte)esrGroup;
				}
			}
		}

	
	private void putORAtomsIntoANDGroup(int fragment, int esrGroup) {
		for (int j=0; j<mMesoFragmentAtom[fragment].length; j++) {
			int atom = mMesoFragmentAtom[fragment][j];
			if (mTHESRType[atom] == Molecule.cESRTypeOr) {
				mTHESRType[atom] = Molecule.cESRTypeAnd;
				mTHESRGroup[atom] = (byte)esrGroup;
				}
			}
		}

	
	protected boolean normalizeESRGroupSwappingAndRemoval(int[] canRank) {
		if (mESRGroupNormalizationInfoList == null)
			return false;

		boolean doneAny = false;
		for (int i=mESRGroupNormalizationInfoList.size()-1; i>=0; i--) {
			boolean done = false;
			ESRGroupNormalizationInfo info = mESRGroupNormalizationInfoList.get(i);
			if (info.action == SWAP_ESR_GROUPS) {
				done = normalizeESRGroupSwapping(info.fragment, info.group, info.type, canRank);
				}
			else if (info.action == REMOVE_ESR_GROUP) {
				done = removeESRGroupFromFragment(info.fragment, canRank);
				}
			if (done) {
				mESRGroupNormalizationInfoList.remove(info);
				for (int j=0; j<mMesoFragmentAtom[info.fragment].length; j++) {
					int atom = mMesoFragmentAtom[info.fragment][j];
					mTHESRTypeNeedsNormalization[atom] = false;
					}
				doneAny = true;
				}
			}
		return doneAny;
		}


	private boolean normalizeESRGroupSwapping(int fragment, int group, int type, int[] canRank) {
		int[] groupAtom = null;
		int[] absAtom = null;
		for (int i=0; i<mMesoFragmentAtom[fragment].length; i++) {
			int atom = mMesoFragmentAtom[fragment][i];
			if (hasParity1or2(atom)) {
				if (mTHESRType[atom] == Molecule.cESRTypeAbs)
					absAtom = addToIntArray(absAtom, (canRank[atom] << 16) + atom);
				else if (mTHESRType[atom] == type
					  && mTHESRGroup[atom] == group)
					groupAtom = addToIntArray(groupAtom, (canRank[atom] << 16) + atom);
				}
			}

		int comparison = new CanonizerRankListComparator().compare(groupAtom, absAtom);
		if (comparison == 0)
			return false;

		if (comparison < 0) {
			for (int i=0; i<mMesoFragmentAtom[fragment].length; i++) {
				int atom = mMesoFragmentAtom[fragment][i];
				if (hasParity1or2(atom)) {
					if (mTHESRType[atom] == Molecule.cESRTypeAbs) {
						mTHESRType[atom] = (byte)type;
						mTHESRGroup[atom] = (byte)group;
						}
					else if (mTHESRType[atom] == type
						  && mTHESRGroup[atom] == group) {
						mTHESRType[atom] = Molecule.cESRTypeAbs;
						mTHESRGroup[atom] = -1;
						}
					}
				}
			}

		return true;
		}


	/**
	 * Changes the lowest ranking ESR group of the fragment to ABS atoms.
	 * Checks first if we have OR groups. If this is the case we must convert one of those.
	 * @param fragment
	 * @param canRank
	 * @return
	 */
	private boolean removeESRGroupFromFragment(int fragment, int[] canRank) {
//System.out.println("removeESRGroupFromFragment() entry");
		int[] fragmentAtom = mMesoFragmentAtom[fragment];

		int esrType = Molecule.cESRTypeAnd; // default
		for (int i=0; i<fragmentAtom.length; i++) {
			int atom = fragmentAtom[i];
			if (mIsStereoCenter[atom]
			 && mTHESRType[atom] == Molecule.cESRTypeOr) {
				esrType = Molecule.cESRTypeOr;
				break;
				}
			}

		int[][] groupMember = new int[Molecule.cESRMaxGroups][];
		for (int i=0; i<fragmentAtom.length; i++) {
			int atom = fragmentAtom[i];
			if (mIsStereoCenter[atom]
			 && mTHESRType[atom] == esrType)
				groupMember[mTHESRGroup[atom]] = addToIntArray(groupMember[mTHESRGroup[atom]],
						(canRank[atom] << 16) + atom);
			}
		for (int i=0; i<Molecule.cESRMaxGroups; i++)
			if (groupMember[i] != null)
				Arrays.sort(groupMember[i]);
		Arrays.sort(groupMember, new CanonizerRankListComparator());

		if (new CanonizerRankListComparator().compare(groupMember[0], groupMember[1]) == 0)
			return false;

//System.out.print("removeESRGroupFromFragment() removing");
		for (int i=0; i<groupMember[0].length; i++) {
			int atom = groupMember[0][i] & 0x0000ffff;
//System.out.print(" "+atom);
			mTHESRType[atom] = Molecule.cESRTypeAbs;
			mTHESRGroup[atom] = -1;
			}
//System.out.println();

		return true;
		}


	static protected int[] addToIntArray(int[] intArray, int intValue) {
		int[] newArray = new int[(intArray == null) ? 1 : intArray.length+1];
		for (int i=0; i<newArray.length-1; i++)
			newArray[i] = intArray[i];
		newArray[newArray.length-1] = intValue;
		return newArray;
		}


	private class ESRGroupFragmentMatrix {
			// Compiles a membership matrix of ESR groups and ABS atoms
			// versus meso fragment numbers and non-fragment-space.
			// From this it derives conditions for meso group normalization.
		private int mAndGroupCount,mOrGroupCount,mGroupCount,
					mNewAndGroupCount,mNewOrGroupCount;
		private boolean[][] mMatrix;
		private int[] mGroupDependence;
		private int[][] mGroupNeighbour;

		public ESRGroupFragmentMatrix() {
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (hasParity1or2(atom)) {
					if (mTHESRType[atom] == Molecule.cESRTypeAnd) {
						if (mAndGroupCount <= mTHESRGroup[atom])
							mAndGroupCount = 1 + mTHESRGroup[atom];
						}
					else if (mTHESRType[atom] == Molecule.cESRTypeOr) {
						if (mOrGroupCount <= mTHESRGroup[atom])
							mOrGroupCount = 1 + mTHESRGroup[atom];
						}
					}
				}

			mGroupCount = mAndGroupCount + mOrGroupCount;
			mMatrix = new boolean[mGroupCount+1][mMesoFragmentAtom.length+1];

			for (int atom=0; atom<mMol.getAtoms(); atom++)
				if (hasParity1or2(atom)
				 && !mIsMesoFragmentMember[atom])
					mMatrix[groupIndex(atom)][mMesoFragmentAtom.length] = true;

			for (int fragment=0; fragment<mMesoFragmentAtom.length; fragment++) {
				for (int j=0; j<mMesoFragmentAtom[fragment].length; j++) {
					int atom = mMesoFragmentAtom[fragment][j];
					if (hasParity1or2(atom))
						mMatrix[groupIndex(atom)][fragment] = true;
					}
				}
			
			// compile group neighbours, i.e. groups that share any fragment
			// of a given group
			mGroupNeighbour = new int[mGroupCount][];
			for (int fragment=0; fragment<mMesoFragmentAtom.length; fragment++) {
				for (int group1=1; group1<mGroupCount; group1++) {
					if (mMatrix[group1][fragment]) {
						for (int group2=0; group2<group1; group2++) {
							if (mMatrix[group2][fragment]) {
								mGroupNeighbour[group1] = addToIntArray(mGroupNeighbour[group1], group2);
								mGroupNeighbour[group2] = addToIntArray(mGroupNeighbour[group2], group1);
								}
							}
						}
					}
				}

			// Find for every ESR group/type whether it is a dependent group/type
			// combination. 
			// i.e. whether its atoms are at the same time members of
			//		 - at least one meso fragment that contains ABS atoms
			//	 and   the remaining atoms outside of any meso fragment
			// or	  - at least two meso fragments that contain ABS atoms
			// or 

			mGroupDependence = new int[mGroupCount+1];
				// mGroupDependence[] stores for every AND/OR group number
				//  -3 -> is a dependent group
				//  -2 -> is free group, i.e. is neither member of an ABS atoms
				//		containing fragment nor of the outside area
				//  -1 -> is member of outside area
				// >=0 -> fragment# with ABS atoms to which all group atoms belong
			for (int group=0; group<mGroupCount; group++) {
				if (mMatrix[group][mMesoFragmentAtom.length])
					mGroupDependence[group] = -1;
				else
					mGroupDependence[group] = -2;
				}

			for (int fragment=0; fragment<mMesoFragmentAtom.length; fragment++) {
				if (mMatrix[mGroupCount][fragment]) { // if fragment has ABS stereo centers
					for (int group=0; group<mGroupCount; group++) {
						if (mMatrix[group][fragment]
						 && mGroupDependence[group] != fragment) {
							if (mGroupDependence[group] == -2)
								mGroupDependence[group] = fragment;
							else
								mGroupDependence[group] = -3;
							}
						}
					}
				}

			for (int anchorGroup=0; anchorGroup<mGroupCount; anchorGroup++) {
				if (mGroupDependence[anchorGroup] >= -1) {
					int[] chainMemberLevel = new int[mGroupCount];
					if (extendAnchorChain(chainMemberLevel, anchorGroup)) {
						for (int group=0; group<mGroupCount; group++) {
							if (chainMemberLevel[group] != 0)
								mGroupDependence[group] = -3;
							}
						}
					}
				}

			// Find dependency cycles:
			// For every pair of non-dependent groups within every fragment
			// try to find a cyclic connection
			for (int fragment=0; fragment<mMesoFragmentAtom.length-1; fragment++) {
				for (int group1=1; group1<mGroupCount; group1++) {
					if (mMatrix[group1][fragment]
					 && mGroupDependence[group1] != -3) {
						for (int group2=0; group2<group1; group2++) {
							if (mMatrix[group2][fragment]
							 && mGroupDependence[group2] != -3) {
								int[] cycle = getDependencyCycle(group1, group2, fragment);
								if (cycle != null) {
									for (int i=0; i<cycle.length; i++)
										mGroupDependence[cycle[i]] = -3;

									removeOneGroupFromCycle(cycle);

									break;
									}
								}
							}
						}
					}
				}
			}

		private boolean extendAnchorChain(int[] chainMemberLevel, int anchorGroup) {
			boolean secondAnchorFound = false;
			int level = 1;
			chainMemberLevel[anchorGroup] = level;
			boolean chainExtentionFound = true;
			while (chainExtentionFound) {
				chainExtentionFound = false;
				for (int group1=0; group1<mGroupCount; group1++) {
					if (chainMemberLevel[group1] == level) {
						for (int group2=0; group2<mGroupCount; group2++) {
							if (chainMemberLevel[group2] == 0
							 && groupsShareFragment(group1, group2)) {
								if (mGroupDependence[group2] == -2) {
									chainMemberLevel[group2] = level+1;
									chainExtentionFound = true;
									}
								else if (mGroupDependence[group2] != mGroupDependence[anchorGroup]) {
									chainMemberLevel[group2] = level+1;
									secondAnchorFound = true;
									}
								}
							}
						}
					}
				level++;
				}
			return secondAnchorFound;
			}

		private int[] getDependencyCycle(int group1, int group2, int startFragment) {
			// check first, if there is a direct back connection in another fragment
			for (int fragment=startFragment+1; fragment<mMesoFragmentAtom.length; fragment++) {
				if (fragment != startFragment
				 && mMatrix[group1][fragment]
				 && mMatrix[group2][fragment]) {
					int[] cycle = new int[2];
					cycle[0] = group2;
					cycle[1] = group1;
					return cycle;
					}
				}

			int[] parentGroup = new int[mGroupCount];
			int[] graphLevel = new int[mGroupCount];
			int[] graphGroup = new int[mGroupCount];

			int current = 0;
			int highest = 0;

			graphGroup[0] = group1;
			graphLevel[group1] = 1;

			while (current <= highest) {
				for (int i=0; i<mGroupNeighbour[graphGroup[current]].length; i++) {
					int candidate = mGroupNeighbour[graphGroup[current]][i];

					if (candidate == group2) {
						if (current == 0)
							continue;

						int cycleLength = graphLevel[graphGroup[current]] + 1;
						int[] cycle = new int[cycleLength];
						cycle[0] = candidate;
						cycle[1] = graphGroup[current];
						for (int j=2; j<cycleLength; j++)
							cycle[j] = parentGroup[cycle[j-1]];

						return cycle;
						}

					if (graphLevel[candidate] == 0
					 && mGroupDependence[candidate] != -3) {
						graphLevel[candidate] = graphLevel[graphGroup[current]] + 1;
						graphGroup[++highest] = candidate;
						parentGroup[candidate] = graphGroup[current];
						}
					}

				current++;
				}

			return null;
			}

		private void removeOneGroupFromCycle(int[] cycle) {
			int minRank = Integer.MAX_VALUE;
			int minGroup = -1;
			int minType = -1;
			int minGroupIndex = -1;
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (hasParity1or2(atom)
				 && mTHESRType[atom] != Molecule.cESRTypeAbs) {
					for (int i=0; i<cycle.length; i++) {
						int esrGroup = this.getESRGroup(cycle[i]);
						int esrType = this.getESRType(cycle[i]);
						if (mTHESRType[atom] == esrType
						 && mTHESRGroup[atom] == esrGroup) {
							if (minRank > mCanRankWithoutStereo[atom] + ((esrType == Molecule.cESRTypeAnd) ? 0x10000 : 0)) {
								minRank = mCanRankWithoutStereo[atom] + ((esrType == Molecule.cESRTypeAnd) ? 0x10000 : 0);
								minGroup = esrGroup;
								minType = esrType;
								minGroupIndex = cycle[i];
								}
							}
						}
					}
				}
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (hasParity1or2(atom)
				 && mTHESRType[atom] == minType
				 && mTHESRGroup[atom] == minGroup) {
					mTHESRType[atom] = Molecule.cESRTypeAbs;
					mTHESRGroup[atom] = -1;
					}
				}

			for (int fragment=0; fragment<mMesoFragmentAtom.length; fragment++)
				 mMatrix[minGroupIndex][fragment] = false;
			}

		private boolean groupsShareFragment(int group1, int group2) {
			for (int fragment=0; fragment<mMesoFragmentAtom.length; fragment++)
				if (mMatrix[group1][fragment]
				 && mMatrix[group2][fragment])
					return true;

			return false;
			}

		private int getDependentGroupCount(int fragment) {
			int count = 0;
			for (int group=0; group<mGroupCount; group++)
				if (mMatrix[group][fragment]
				 && mGroupDependence[group] == -3)
					count++;

			return count;
			}

		private int getDependentType(int fragment) {
			for (int group=0; group<mGroupCount; group++)
				if (mMatrix[group][fragment]
				 && mGroupDependence[group] == -3)
					return getESRType(group);
			return -1;
			}

		private int getDependentGroup(int fragment) {
			for (int group=0; group<mGroupCount; group++)
				if (mMatrix[group][fragment]
				 && mGroupDependence[group] == -3)
					return getESRGroup(group);
			return -1;
			}

		private int getESRType(int group) {
				// extract ESR type from internal group/type index
			return (group < mAndGroupCount) ? Molecule.cESRTypeAnd
				 : (group < mGroupCount)	? Molecule.cESRTypeOr
											: Molecule.cESRTypeAbs;
				}

		private int getESRGroup(int group) {
				// extract ESR group from internal group/type index
			return (group < mAndGroupCount) ? group
				 : (group < mGroupCount)	? group - mAndGroupCount
											: -1;
			}

		private int groupIndex(int atom) {
			int type = mTHESRType[atom];
			int group = mTHESRGroup[atom];
			return (type == Molecule.cESRTypeAbs) ? mGroupCount
				 : (type == Molecule.cESRTypeAnd) ? group
				 : mAndGroupCount + group;
			}

		private void cutTiesOfIndependentGroups(int fragment) {
			for (int group=0; group<mGroupCount; group++) {
				if (mMatrix[group][fragment]
				 && mGroupDependence[group] != -3) {
					for (int f=0; f<=mMesoFragmentAtom.length; f++) {
						if (f != fragment && mMatrix[group][f]) {
							mMatrix[group][fragment] = false;
							int oldESRGroup = getESRGroup(group);
							int newESRGroup = newESRGroup(getESRType(group));
							for (int i=0; i<mMesoFragmentAtom[fragment].length; i++) {
								int atom = mMesoFragmentAtom[fragment][i];
								if (hasParity1or2(atom)
								 && mTHESRGroup[atom] == oldESRGroup)
									mTHESRGroup[atom] = (byte)newESRGroup;
								}
							}
						}
					}
				}
			}

		private int newESRGroup(int esrType) {
			return (esrType == Molecule.cESRTypeAnd) ?
					mAndGroupCount + mNewAndGroupCount++
					: mOrGroupCount + mNewOrGroupCount++;
			}
		}
	}


class MesoFragmentBranch {
	private int[] mirrorAtom;
	private int mirrorAtomIndex;
	public int neighbourIndex;
	public int current;

	public MesoFragmentBranch(int[] mirrorAtom, int neighbourIndex, int current) {
		this.mirrorAtom = mirrorAtom;
		this.neighbourIndex = neighbourIndex;
		this.current = current;
		this.mirrorAtomIndex = 1;
		}

	public int getNextMirrorAtom() {
		return mirrorAtomIndex < mirrorAtom.length ? mirrorAtom[mirrorAtomIndex++] : -1;
		}

	public boolean hasNextMirrorAtom() {
		return mirrorAtomIndex < mirrorAtom.length;
		}
	}

class MesoFragmentMembers implements Comparable<MesoFragmentMembers> {
	public boolean[] isMember;
	public int[] memberAtom;

	public MesoFragmentMembers(int atoms) {
		isMember = new boolean[atoms];
		}

	public void add(int atom) {
		isMember[atom] = true;
		}

	private void consolidate() {
		int count = 0;
		for (boolean is:isMember)
			if (is)
				count++;
		memberAtom = new int[count];
		count = 0;
		for (int atom=0; atom<isMember.length; atom++)
			if (isMember[atom])
				memberAtom[count++] = atom;
		}

	public boolean hasStereoCenters(boolean[] isStereoCenter) {
		consolidate();

		for (int j=0; j<memberAtom.length; j++)
			if (isStereoCenter[memberAtom[j]])
				return true;

		return false;
		}

	@Override
	public int compareTo(MesoFragmentMembers members) {
		if (memberAtom.length != members.memberAtom.length)
			return (memberAtom.length < members.memberAtom.length) ? -1 : 1;

		for (int i=0; i<memberAtom.length; i++)
			if (memberAtom[i] != members.memberAtom[i])
				return (memberAtom[i] < members.memberAtom[i]) ? -1 : 1;

		return 0;
		}
	}

class ESRGroupNormalizationInfo {
	public int fragment;
	public int action;
	public int group;
	public int type;

	public ESRGroupNormalizationInfo(int fragment, int action, int group, int type) {
		this.fragment = fragment;
		this.action = action;
		this.group = group;
		this.type = type;
		}
	}

class CanonizerRankListComparator implements Comparator<int[]> {
	public int compare(int[] o1, int[] o2) {
		if (o1 == null) // put null arrays at the end of the list values
			return (o2 == null) ? 0 : 1;
		if (o2 == null)
			return -1;
		int count = Math.min(o1.length, o2.length);
		for (int i=0; i<count; i++)
			if ((o1[i] & 0xffff0000) != (o2[i] & 0xffff0000))
				return ((o1[i] & 0xffff0000) < (o2[i] & 0xffff0000)) ? -1 : 1;
		return (o1.length == o2.length) ? 0 : (o1.length < o2.length) ? -1 : 1;
		}
	}
