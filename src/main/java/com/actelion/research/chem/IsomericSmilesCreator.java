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
*/

package com.actelion.research.chem;

import com.actelion.research.chem.reaction.Reaction;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IsomericSmilesCreator {
	public static final int MODE_CREATE_SMARTS = 1;
	public static final int MODE_INCLUDE_MAPPING = 2;

	private StereoMolecule mMol;
	private Canonizer mCanonizer;
	private String mSmiles;
	private int mMode;
	private int[] mAtomRank;
	private int[] mClosureNumber;
	private int[] mSmilesIndex;
	private int[] mClosureBuffer;
	private int[][] mKnownTHCountInESRGroup;
	private List<SmilesAtom> mGraphAtomList;
	private boolean[] mAtomUsed;
	private boolean[] mBondUsed;
	private boolean[] mClosureOpened;
	private boolean[] mPseudoStereoGroupInversion;
	private boolean[] mPseudoStereoGroupInitialized;

	/**
	 * Convenience method to generate a canonical and isomeric SMILES from a given molecule.
	 * @param mol
	 * @return
	 */
	public static String createSmiles(StereoMolecule mol) {
		return new IsomericSmilesCreator(mol, 0).getSmiles();
		}

	/**
	 * Convenience method to generate a canonical and isomeric SMILES from a given molecule.
	 * @param mol
	 * @return
	 */
	public static String createSmarts(StereoMolecule mol) {
		return new IsomericSmilesCreator(mol, MODE_CREATE_SMARTS).getSmiles();
	}

	public static String createReactionSmiles(Reaction rxn) {
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<rxn.getReactants(); i++) {
			if (i != 0)
				sb.append('.');
			sb.append(new IsomericSmilesCreator(rxn.getReactant(i), MODE_INCLUDE_MAPPING).getSmiles());
		}
		sb.append('>');
		for (int i=0; i<rxn.getCatalysts(); i++) {
			if (i != 0)
				sb.append('.');
			sb.append(new IsomericSmilesCreator(rxn.getCatalyst(i)).getSmiles());
		}
		sb.append('>');
		for (int i=0; i<rxn.getProducts(); i++) {
			if (i != 0)
				sb.append('.');
			sb.append(new IsomericSmilesCreator(rxn.getProduct(i), MODE_INCLUDE_MAPPING).getSmiles());
		}
		return sb.toString();
	}

	/**
	 * Creates an IsomericSmilesCreator, which doesn't include atom mapping into generated smiles.
	 * @param mol
	 */
	public IsomericSmilesCreator(StereoMolecule mol) {
		this(mol, false);
	}

	/**
	 * Creates an IsomericSmilesCreator with the given mode.
	 * @param mol
	 * @param mode combination of MODE_... flags
	 */
	public IsomericSmilesCreator(StereoMolecule mol, int mode) {
		mMol = mol;
		mMode = mode;
	}

	/**
	 * Creates an IsomericSmilesCreator, which may include atom mapping numbers into generated smiles.
	 * @param mol
	 * @param includeAtomMapping
	 */
	@Deprecated
	public IsomericSmilesCreator(StereoMolecule mol, boolean includeAtomMapping) {
		this(mol, includeAtomMapping ? MODE_INCLUDE_MAPPING : 0);
	}

	public String getSmiles() {
		if (mSmiles == null)
			mSmiles = createSmiles();

		return mSmiles;
	}

	private String createSmiles() {
		if (mMol == null || mMol.getAllAtoms() == 0)
			return "";

		mMol.ensureHelperArrays(Molecule.cHelperParities);

		mCanonizer = new Canonizer(mMol, Canonizer.CREATE_SYMMETRY_RANK | Canonizer.CREATE_PSEUDO_STEREO_GROUPS);

		int groupCount = mCanonizer.getPseudoStereoGroupCount();
		mPseudoStereoGroupInversion = new boolean[groupCount+1];
		mPseudoStereoGroupInitialized = new boolean[groupCount+1];

		mKnownTHCountInESRGroup = new int[2][Molecule.cESRMaxGroups];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			int type = mMol.getAtomESRType(atom) - 1;
			if (type != -1)
				mKnownTHCountInESRGroup[type][mMol.getAtomESRGroup(atom)]++;
		}

		generateCanonicalTree();
		findRingClosures();
		calculateEZBonds();

		StringBuilder builder = new StringBuilder();
		StringBuilder buffer = new StringBuilder();
		boolean isFirst = true;
		for (SmilesAtom smilesAtom:mGraphAtomList) {
			if (smilesAtom.parent == -1) {
				if (isFirst)
					isFirst = false;
				else
					builder.append('.');
			}
			addAtomString(smilesAtom, builder, buffer);
		}

		return builder.toString();
	}

	private void generateCanonicalTree() {
		mAtomRank = mCanonizer.getFinalRank();
		mAtomUsed = new boolean[mMol.getAtoms()];
		mBondUsed = new boolean[mMol.getBonds()];
		mGraphAtomList = new ArrayList<>();

		int atom = findUnusedStartAtom();
		while (atom != -1) {
			int graphIndex = mGraphAtomList.size();
			addToGraph(new SmilesAtom(atom, -1, false, false), graphIndex);
			if (mMol.getConnAtoms(atom) != 0) {
				addHighestRankingChain(graphIndex, false);

				while (graphIndex < mGraphAtomList.size() - 1) {
					while (hasUnusedNeighborAtom(mGraphAtomList.get(graphIndex).atom))
						addHighestRankingChain(graphIndex, true);
					graphIndex++;
				}
			}

			// we may have multiple unconnected fragments
			atom = findUnusedStartAtom();
		}

		// assign a smiles position index to every atom (needed for stereo assignment)
		mSmilesIndex = new int[mMol.getAtoms()];
		int index = 0;
		for (SmilesAtom smilesAtom:mGraphAtomList)
			mSmilesIndex[smilesAtom.atom] = index++;
	}

	private void findRingClosures() {
		boolean[] closureNumberUsed = new boolean[mMol.getBonds()];
		mClosureNumber = new int[mMol.getBonds()];

		for (SmilesAtom smilesAtom:mGraphAtomList) {
			for (int i=0; i<mMol.getConnAtoms(smilesAtom.atom); i++) {
				int bond = mMol.getConnBond(smilesAtom.atom, i);
				closureNumberUsed[mClosureNumber[bond]] = false;
			}

			int index = getUnusedConnBondIndex(smilesAtom.atom);
			while (index != -1) {
				int closureBond = mMol.getConnBond(smilesAtom.atom, index);
				mBondUsed[closureBond] = true;

				int closureNumber = 1;
				while (closureNumberUsed[closureNumber])
					closureNumber++;

				mClosureNumber[closureBond] = closureNumber;
				closureNumberUsed[closureNumber] = true;

				index = getUnusedConnBondIndex(smilesAtom.atom);
			}
		}

		mClosureOpened = new boolean[mMol.getBonds()];
		mClosureBuffer = new int[8];
	}

	private void calculateEZBonds() {
		for (SmilesAtom currentSA:mGraphAtomList) {
			if (currentSA.parent != -1) {
				int bond = mMol.getBond(currentSA.atom, currentSA.parent);
				if (!mMol.isBINAPChiralityBond(bond)
				 && !mMol.isSmallRingBond(bond)
				 && (mMol.getBondParity(bond) == Molecule.cBondParityEor1
				  || mMol.getBondParity(bond) == Molecule.cBondParityZor2)) {
					SmilesAtom parentSA = mGraphAtomList.get(mSmilesIndex[currentSA.parent]);

					int halfParity1 = parentSA.ezHalfParity;
					if (halfParity1 == 0)
						halfParity1 = parentSA.ezHalfParity = 1;

					int halfParity = halfParity1;	// we assume an E-double bond

					if (mMol.getConnAtoms(parentSA.atom) == 3) {
						for (int i=0; i<mMol.getConnAtoms(parentSA.atom); i++) {
							int connAtom = mMol.getConnAtom(parentSA.atom, i);
							if (connAtom != parentSA.parent && connAtom != currentSA.atom) {
								// If we have a branch at the left double bond atom, we need to give it the same
								// half-parity given to the parent atom (other atom at the left double bond end),
								// because this is translated into the same sign ('/' or '\'), which
								// actually means on the different side of the double bond, because
								// "The 'visual interpretation' of the 'up-ness' or 'down-ness' of each single bond
								//  is relative to the carbon atom, not the double bond" (opensmiles.org).
								SmilesAtom branchSA = mGraphAtomList.get(mSmilesIndex[connAtom]);
								if (branchSA.parent == parentSA.atom)
									branchSA.ezHalfParity = halfParity1;	// same half-parity

								if (connAtom < parentSA.parent)
									halfParity = 3 - halfParity;	// invert
								break;
							}
						}
					}

					if (mMol.getBondParity(bond) == Molecule.cBondParityZor2)
						halfParity = 3 - halfParity;

					for (int i=0; i<mMol.getConnAtoms(currentSA.atom); i++) {
						int childAtom = mMol.getConnAtom(currentSA.atom, i);
						if (childAtom != currentSA.parent) {
							int halfParity2 = halfParity;
							if (mMol.getConnAtoms(currentSA.atom) == 3) {
								for (int j=0; j<mMol.getConnAtoms(currentSA.atom); j++) {
									int connAtom = mMol.getConnAtom(currentSA.atom, j);
									if (connAtom != currentSA.parent && connAtom != childAtom) {
										if (connAtom < childAtom)
											halfParity2 = 3 - halfParity2;	// invert
										break;
									}
								}
							}

							if (mMol.isBondParityPseudo(bond)) {
								int group = mCanonizer.getPseudoEZGroup(bond);
								if (!mPseudoStereoGroupInitialized[group]) {
									mPseudoStereoGroupInitialized[group] = true;
									mPseudoStereoGroupInversion[group] = (halfParity2 == 2);
								}
								if (mPseudoStereoGroupInversion[group])
									halfParity2 = 3 - halfParity2;	// invert
							}

							SmilesAtom childSA = mGraphAtomList.get(mSmilesIndex[childAtom]);
							if (childSA.parent == currentSA.atom)
								childSA.ezHalfParity = halfParity2;
//							else we have a ring closure
						}
					}
				}
			}
		}
	}

	/**
	 * Of the not yet used atoms find that atom with the lowest number of neighbour atoms.
	 * Exception: Atom with zero neighbors are only considered if no other atoms are left.
	 * If more than one atom with the same neighbor count exist, then the one of them with
	 * the lowest symmetry rank is returned.
	 * @return -1 if all atoms were already used
	 */
	private int findUnusedStartAtom() {
		int startAtom = -1;
		int startRank = Integer.MAX_VALUE;

		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (!mAtomUsed[atom]) {
				int rank = mAtomRank[atom];

				// prefer non-exclude group atoms
				if ((mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFExcludeGroup) != 0)
					rank += 0x40000000;

				// prefer lower neighbour count except zero neighbours
				if (mMol.getConnAtoms(atom) == 0)
					rank += 0x3F000000;
				else
					rank += mMol.getConnAtoms(atom) << 24;

				if (startRank > rank) {
					startRank = rank;
					startAtom = atom;
				}
			}
		}

		return startAtom;
	}

	private boolean hasUnusedNeighborAtom(int atom) {
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (!mAtomUsed[mMol.getConnAtom(atom, i)])
				return true;
		return false;
	}

	private void addHighestRankingChain(int graphIndex, boolean isSideChain) {
		boolean isFirst = true;
		int parent = mGraphAtomList.get(graphIndex).atom;
		int neighborIndex = getUnusedConnAtomIndex(parent);
		while (neighborIndex != -1) {
			int atom = mMol.getConnAtom(parent, neighborIndex);
			int bond =  mMol.getConnBond(parent, neighborIndex);
			neighborIndex = getUnusedConnAtomIndex(atom);
			addToGraph(new SmilesAtom(atom, parent,
					isSideChain && isFirst,
					isSideChain && neighborIndex == -1), ++graphIndex);
			parent = atom;
			isFirst = false;
		}
	}

	/**
	 * Find and return the preferred unused neighbor index.
	 * Higher bond orders are preferred to ensure that finally all remaining
	 * ring closures are single bonds. Within equal bond order neighbors
	 * those with a higher atom rank are preferred.
	 * @param atom
	 * @return highest ranking unused neighbor index
	 */
	private int getUnusedConnAtomIndex(int atom) {
		int bestIndex = -1;
		int bestRank = -1;
		for (int i=0; i<mMol.getConnAtoms(atom); i++) {
			int connAtom = mMol.getConnAtom(atom, i);
			int rank = (mMol.getConnBondOrder(atom, i) << 24) + mAtomRank[connAtom];
			if (!mAtomUsed[connAtom]
			 && (bestIndex == -1 || bestRank < rank)) {
				bestIndex = i;
				bestRank = rank;
			}
		}
		return bestIndex;
	}

	private int getUnusedConnBondIndex(int atom) {
		int bestIndex = -1;
		for (int i=0; i<mMol.getConnAtoms(atom); i++)
			if (!mBondUsed[mMol.getConnBond(atom, i)]
			 && (bestIndex == -1 || mAtomRank[mMol.getConnAtom(atom, bestIndex)] < mAtomRank[mMol.getConnAtom(atom, i)]))
				bestIndex = i;

		return bestIndex;
	}

	private void addToGraph(SmilesAtom smilesAtom, int listIndex) {
		mGraphAtomList.add(listIndex, smilesAtom);
		mAtomUsed[smilesAtom.atom] = true;
		if (smilesAtom.parent != -1)
			mBondUsed[mMol.getBond(smilesAtom.atom, smilesAtom.parent)] = true;
	}

	private void addAtomString(SmilesAtom smilesAtom, StringBuilder builder, StringBuilder buffer) {
		int atom = smilesAtom.atom;
		int parent = smilesAtom.parent;

		boolean isAnyAtom = (mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFAny) != 0;
		int[] atomList = mMol.getAtomList(atom);

		String label = (atomList != null) ? createMultiAtomLabel(atom, atomList, buffer) : isAnyAtom ? "*" : mMol.getAtomLabel(atom);
		if (!isAnyAtom
		 && atomList == null
		 && mMol.isAromaticAtom(atom)
		 && (mMol.getAtomPi(atom) != 0
		  || (mMol.getAtomAbnormalValence(atom) == -1
		   && mMol.getAtomRadical(atom) == Molecule.cAtomRadicalStateNone)))
			label = label.toLowerCase();

		if (smilesAtom.isSideChainStart)
			builder.append('(');

		if (parent != -1)
			appendBondOrderSymbol(smilesAtom, builder);

		int charge = mMol.getAtomCharge(atom);
		if (charge == 0 && (mMode & MODE_CREATE_SMARTS) != 0) {
			// Because SMARTS don't know a charge query feature, we set the atom charge in case of neg/pos charge required.
			// There is no way to require an uncharged atom (Molecule.cAtomQFNotChargeNeg | Molecule.cAtomQFNotChargePos).
			int chargeFeatures = mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFCharge;
			if (chargeFeatures == (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargePos))
				charge = -1;    // 'require negative charge' is translated into charge := -1
			else if (chargeFeatures == (Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargeNeg))
				charge = 1;    // 'require positive charge' is translated into charge := +1
			}

		int isotop = mMol.getAtomMass(atom);
		int mapNo = (mMode & MODE_INCLUDE_MAPPING) != 0 ? mMol.getAtomMapNo(atom) : 0;

		String smartsFeatures = (mMode & MODE_CREATE_SMARTS) != 0 ? getAtomSMARTSFeatures(atom, buffer) : null;

		boolean useBrackets =
				(!isAnyAtom && !isOrganic(mMol.getAtomicNo(atom)))
				|| atomList != null
				|| qualifiesForAtomParity(atom)
				|| (mMol.isAromaticAtom(atom) && mMol.getAtomPi(atom) == 0)
				|| charge != 0
				|| isotop != 0
				|| mapNo != 0
				|| mMol.getAtomAbnormalValence(atom) != -1
				|| mMol.getAtomRadical(atom) != Molecule.cAtomRadicalStateNone
				|| smartsFeatures != null;

		if (useBrackets)
			builder.append('[');

		if (isotop != 0)
			builder.append(isotop);

		builder.append(label);

		if (qualifiesForAtomParity(atom))
			builder.append(getAtomParitySymbol(atom, parent));

		if ((mMode & MODE_CREATE_SMARTS) == 0 && useBrackets) {
			int hCount = mMol.getPlainHydrogens(atom);
			if (hCount == 1)
				builder.append("H");
			else if (hCount > 1)
				builder.append("H"+hCount);
			}

		if (charge != 0) {
			builder.append(charge > 0 ? '+' : '-');
			if (Math.abs(charge) > 1)
				builder.append(Math.abs(charge));
			}

		if (smartsFeatures != null)
			builder.append(smartsFeatures);

		if (mapNo != 0) {
			builder.append(':');
			builder.append(mapNo);
			}

		if (useBrackets)
			builder.append(']');

		appendClosureBonds(smilesAtom, builder);

		if (smilesAtom.isSideChainEnd)
			builder.append(')');
	}

	private String createMultiAtomLabel(int atom, int[] atomicNo, StringBuilder buffer) {
		buffer.setLength(0);
		boolean isAromatic = mMol.isAromaticAtom(atom);

		for (int a:atomicNo) {
			if (buffer.length() != 0)
				buffer.append(',');
			String label = Molecule.cAtomLabel[a];
				buffer.append(isAromatic ? label.toLowerCase() : label);
			}

		return buffer.toString();
		}

	private String getAtomSMARTSFeatures(int atom, StringBuilder buffer) {
		buffer.setLength(0);

		int queryFeatures = mMol.getAtomQueryFeatures(atom);

		// SMARTS don't distinguish between a charged atom and atom charge query features
		int chargeFeatures = queryFeatures & Molecule.cAtomQFCharge;
		switch (chargeFeatures) {
			case Molecule.cAtomQFNotChargeNeg | Molecule.cAtomQFNotChargePos:
				buffer.append("+0");   // 'require negative charge' is translated into charge := -1
				break;
			case Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargePos:
				if (mMol.getAtomCharge(atom) == 0)
					buffer.append("-");   // 'require negative charge' is translated into charge := -1
				break;
			case Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargeNeg:
				if (mMol.getAtomCharge(atom) == 0)
					buffer.append("+");   // 'require positive charge' is translated into charge := +1
				break;
			}

		int aromState = queryFeatures & Molecule.cAtomQFAromState;
		if (aromState == Molecule.cAtomQFAromatic)
			buffer.append(";a");
		else if (aromState == Molecule.cAtomQFNotAromatic)
			buffer.append(";A");

		int hydrogenQueryFeatures = queryFeatures & Molecule.cAtomQFHydrogen;
		if (hydrogenQueryFeatures != 0) {
			if (hydrogenQueryFeatures == (Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen))
				buffer.append(";H0");
			else if (hydrogenQueryFeatures == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen))
				buffer.append(";H1");
			else if (hydrogenQueryFeatures == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot3Hydrogen))
				buffer.append(";H2");
			else if (hydrogenQueryFeatures == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen | Molecule.cAtomQFNot2Hydrogen))
				buffer.append(";H3");
			else if (hydrogenQueryFeatures == Molecule.cAtomQFNot0Hydrogen)
				buffer.append(";!H0");
			else if (hydrogenQueryFeatures == (Molecule.cAtomQFNot0Hydrogen | Molecule.cAtomQFNot1Hydrogen))
				buffer.append(";!H0;!H1");
			else if (hydrogenQueryFeatures == (Molecule.cAtomQFNot2Hydrogen | Molecule.cAtomQFNot3Hydrogen))
				buffer.append(";!H2;!H3");
			else if (hydrogenQueryFeatures == Molecule.cAtomQFNot3Hydrogen)
				buffer.append(";!H3");
			}

		// Atom membership of SSSR rings cannot be directly translated into number of ring bonds.
		// We try to get close...
		int ringState = queryFeatures & Molecule.cAtomQFRingState;
		switch (ringState) {
			case Molecule.cAtomQFNotChain:
				buffer.append(";!R0");
				break;
			case Molecule.cAtomQFNot2RingBonds:
				buffer.append(";!R1");
				break;
			case Molecule.cAtomQFNot3RingBonds:
				buffer.append(";!R2");
				break;
			case Molecule.cAtomQFNot4RingBonds:
				buffer.append(";!R3");
				break;
			case Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds:
				buffer.append(";R0");
				break;
			case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds:
				buffer.append(";R1");
				break;
			case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot4RingBonds:
				buffer.append(";R2");
				break;
			case Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds:
				buffer.append(";R3");
				break;
		}

		int ringSize = (queryFeatures & Molecule.cAtomQFRingSize) >> Molecule.cAtomQFRingSizeShift;
		if (ringSize != 0)
			buffer.append(";r"+ringSize);

		int neighbourFeatures = queryFeatures & Molecule.cAtomQFNeighbours;
		switch (neighbourFeatures) {
			case Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot1Neighbour:
				buffer.append(";D1");   // exactly 1
				break;
			case Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot2Neighbours:
				buffer.append(";D2");   // exactly 2
				break;
			case Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot3Neighbours:
				buffer.append(";D3");   // exactly 3
				break;
			case Molecule.cAtomQFNot3Neighbours | Molecule.cAtomQFNot4Neighbours:
				buffer.append(";!D3;!D4");   // less than 3
				break;
			case Molecule.cAtomQFNot4Neighbours:
				buffer.append(";!D4");   // less than 4
				break;
			case Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour:
				buffer.append(";!D0;!D1");   // more than 1
				break;
			case Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour | Molecule.cAtomQFNot2Neighbours:
				buffer.append(";!D0;!D1;!D2");   // more than 2
				break;
			case Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot4Neighbours:
				buffer.append(";!D0;!D1;!D2;!D3");   // more than 3
				break;
		}

		if ((queryFeatures & Molecule.cAtomQFNoMoreNeighbours) != 0)
			buffer.append(";D"+mMol.getConnAtoms(atom));   // Convert into exact explicit neighbour count 'D'

		if ((queryFeatures & Molecule.cAtomQFMoreNeighbours) != 0)
			buffer.append(";!D"+mMol.getConnAtoms(atom));  // Convert into exact explicit neighbour count 'D'

		return buffer.length() == 0 ? null : buffer.toString();
		}

	/**
	 * Don't store nitrogen parities in SMILES unless we have a quarternary nitrogen.
	 * Some software packages seem to have promblems when decoding parities on nitrogen atoms.
	 * @param atom
	 * @return
	 */
	private boolean qualifiesForAtomParity(int atom) {
		return (mMol.getAtomParity(atom) == Molecule.cAtomParity1
			 || mMol.getAtomParity(atom) == Molecule.cAtomParity2)
			&& !isSingleKnownStereoCenterInESRGroup(atom)
			&& (mMol.getAtomicNo(atom) != 7 || mMol.getAtomCharge(atom) > 0);
	}

	private boolean isSingleKnownStereoCenterInESRGroup(int atom) {
		int type = mMol.getAtomESRType(atom) - 1;
		return type == -1 ? false : mKnownTHCountInESRGroup[type][mMol.getAtomESRGroup(atom)] <= 1;
	}

	private void appendClosureBonds(SmilesAtom smilesAtom, StringBuilder builder) {
		int closureCount = 0;
		for (int i=0; i<mMol.getConnAtoms(smilesAtom.atom); i++) {
			int bond = mMol.getConnBond(smilesAtom.atom, i);
			if (mClosureNumber[bond] != 0) {
				int isOpenFlag = mClosureOpened[bond] ? 0 : 0x40000000;
				mClosureBuffer[closureCount++] = isOpenFlag | (mClosureNumber[bond] << 20) | bond;
			}
		}
		if (closureCount != 0) {
			// when sorting, then put and handle open closures first
			Arrays.sort(mClosureBuffer, 0, closureCount); // we must sort to be canonical
			for (int i=0; i<closureCount; i++) {
				int bond = mClosureBuffer[i] & 0x0003FFFF;
				int closureNumber = ((mClosureBuffer[i] & 0x3FFC0000) >> 20);
				if (!mClosureOpened[bond]) {
					mClosureOpened[bond] = true;
					appendBondOrderSymbol(bond, smilesAtom.atom, builder);
				}
				if (closureNumber > 9)
					builder.append('%');
				builder.append(closureNumber);
			}
		}
	}

	private void appendBondOrderSymbol(SmilesAtom smilesAtom, StringBuilder builder) {
		if (smilesAtom.ezHalfParity != 0) {
			builder.append(smilesAtom.ezHalfParity == 1 ? '/' : '\\');
			return;
		}
		appendBondOrderSymbol(mMol.getBond(smilesAtom.atom, smilesAtom.parent), smilesAtom.parent, builder);
	}

	private void appendBondOrderSymbol(int bond, int parentAtom, StringBuilder builder) {
		int startLength = builder.length();
		if (mMode == MODE_CREATE_SMARTS) {
			int bondTypes = mMol.getBondQueryFeatures(Molecule.cBondQFBondTypes);
			if (bondTypes != 0) {
				if ((bondTypes & Molecule.cBondTypeSingle) != 0)
					builder.append('-');
				if ((bondTypes & Molecule.cBondTypeDouble) != 0) {
					if (builder.length() != startLength)
						builder.append(',');
					builder.append('=');
					}
				if ((bondTypes & Molecule.cBondTypeTriple) != 0) {
					if (builder.length() != startLength)
						builder.append(',');
					builder.append('#');
					}
				if ((bondTypes & Molecule.cBondTypeDelocalized) != 0) {
					if (builder.length() != startLength)
						builder.append(',');
					builder.append(':');
					}
				if ((bondTypes & Molecule.cBondTypeMetalLigand) != 0) {
					if (builder.length() != startLength)
						builder.append(',');
					builder.append(mMol.isMetalAtom(parentAtom) ? "<-" : "->");
					}
				}
			}
		if (startLength == builder.length() && !mMol.isAromaticBond(bond)) {
			int bondType = mMol.getBondType(bond) & Molecule.cBondTypeMaskSimple;
			if (bondType == Molecule.cBondTypeSingle) {
				if (mMol.isAromaticAtom(mMol.getBondAtom(0, bond))
				 && mMol.isAromaticAtom(mMol.getBondAtom(1, bond)))
					builder.append('-');
			} else if (bondType == Molecule.cBondTypeDouble)
				builder.append('=');
			else if (bondType == Molecule.cBondTypeTriple)
				builder.append('#');
//				else if (order == Molecule.cBondTypeQuadruple)	// not supported by Molecule
//					builder.append('$');
			else if (bondType == Molecule.cBondTypeDelocalized)
				builder.append(':');
			else if (bondType == Molecule.cBondTypeMetalLigand)
				builder.append(mMol.isMetalAtom(parentAtom) ? "<-" : "->");
		}
		if (mMode == MODE_CREATE_SMARTS) {
			String gap = (startLength == builder.length()) ? "" : ";";
			int ringState = mMol.getBondQueryFeatures(bond) & Molecule.cBondQFRingState;
			if (ringState == Molecule.cBondQFRing)
				builder.append(gap+"@");
			else if (ringState == Molecule.cBondQFNotRing)
				builder.append(gap+"!@");
		}
	}

	private String getAtomParitySymbol(int atom, int parent) {
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
						&& ((mSmilesIndex[neighbour[0]] < mSmilesIndex[neighbour[1]])
						^(neighbour[0] < neighbour[1])))
					inversion = !inversion;
			}
		}
		else {
			// we map 3 neighbors of stereo center to 3 neighbors in smiles
			// (excluding from-atom and including implicit H, if exists)
			int[] neighborAtom = new int[3];
			int[] neighborRank = new int[3];
			int index = 0;
			for (int i=0; i<mMol.getConnAtoms(atom); i++) {
				int connAtom = mMol.getConnAtom(atom, i);
				if (connAtom != parent) {
					neighborAtom[index] = connAtom;
					neighborRank[index++] = getSmilesRank(atom, i);
				}
			}
			if (index == 2) {	// we have an implicit hydrogen (or lone pair in case of chiral sulfur)
				if (mMol.getImplicitHydrogens(atom) == 0)    // no hydrogen
					neighborAtom[index] = atom; // we use the central atom's index for the lone pair
				else
					neighborAtom[index] = Integer.MAX_VALUE; // index must be higher than any non-hydrogen neighbor
				neighborRank[index++] = 8 * mSmilesIndex[atom]; // smaller than any potential closure neighbor
			}

			if (neighborRank[0] > neighborRank[1])
				inversion = !inversion;
			if (neighborRank[0] > neighborRank[2])
				inversion = !inversion;
			if (neighborRank[1] > neighborRank[2])
				inversion = !inversion;

			if (neighborAtom[0] > neighborAtom[1])
				inversion = !inversion;
			if (neighborAtom[0] > neighborAtom[2])
				inversion = !inversion;
			if (neighborAtom[1] > neighborAtom[2])
				inversion = !inversion;

			for (int i=0; i<3; i++)
				if (parent > neighborAtom[i])
					inversion = !inversion;
		}

		boolean isClockwise = ((mMol.getAtomParity(atom) == Molecule.cAtomParity1) ^ inversion);

		if (mMol.isAtomParityPseudo(atom)) {
			int group = mCanonizer.getPseudoTHGroup(atom);
			if (!mPseudoStereoGroupInitialized[group]) {
				mPseudoStereoGroupInitialized[group] = true;
				mPseudoStereoGroupInversion[group] = isClockwise;
			}
			if (mPseudoStereoGroupInversion[group])
				isClockwise = !isClockwise;
		}

		return isClockwise ? "@@" : "@";
	}

	/**
	 * @param atom for which to return a neighbor's smiles rank
	 * @param neighborIndex index for getConnAtoms() to get neighbor atom and bond
	 * @return neighbor's position rank in smiles from a perspective of atom using closure digit positions for closure neighbors
	 */
	private int getSmilesRank(int atom, int neighborIndex) {
		int bond = mMol.getConnBond(atom, neighborIndex);
		if (mClosureNumber[bond] != 0) {
			// if neighbor is attached via a closure digit, then the rank is based primarily on atom's position
			// in the smiles and secondary on the count of other closures at atom that precede this closure
			int rank = 8 * mSmilesIndex[atom] + 1;
			for (int i=0; i<neighborIndex; i++)
				if (mClosureNumber[mMol.getConnAtom(atom, i)] != 0)
					rank++;
			return rank;
			}

		// if the neighbor is not a closure return a rank based on its atom index in the smiles
		return 8 * mSmilesIndex[mMol.getConnAtom(atom, neighborIndex)];
	}

	private boolean isOrganic (int atomicNo) {
		return (atomicNo >= 5 && atomicNo <= 9)     // B,C,N,O,F
			|| (atomicNo >= 15 && atomicNo <= 17)   // P,S,Cl
			|| atomicNo == 35                       // Br
			|| atomicNo == 53;                      // I
	}
}

class SmilesAtom {
	public int atom,parent,ezHalfParity;
	public boolean isSideChainStart,isSideChainEnd;

	public SmilesAtom(int atom, int parent, boolean isSideChainStart, boolean isSideChainEnd) {
		this.atom = atom;
		this.parent = parent;
		this.isSideChainStart = isSideChainStart;
		this.isSideChainEnd = isSideChainEnd;
	}
}