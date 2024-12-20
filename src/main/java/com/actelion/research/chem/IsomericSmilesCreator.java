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

package com.actelion.research.chem;

import com.actelion.research.chem.reaction.Reaction;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IsomericSmilesCreator {
	public static final int MODE_CREATE_SMARTS = 1;
	public static final int MODE_INCLUDE_MAPPING = 2;
	public static final int MODE_KEKULIZED_OUTPUT = 4;  // no lower case atom labels and single/double bonds to represent aromaticity

	private StereoMolecule mMol;
	private Canonizer mCanonizer;
	private String mSmiles;
	private int mMode;
	private int[] mAtomRank;
	private int[] mClosureNumber;
	private int[] mSmilesIndex;
	private int[][] mKnownTHCountInESRGroup;
	private List<SmilesAtom> mGraphAtomList;
	private boolean[] mAtomUsed;
	private boolean[] mBondUsed;
	private boolean[] mPseudoStereoGroupInversion;
	private boolean[] mPseudoStereoGroupInitialized;
	private int[] mEZHalfParity;

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

	public static String createReactionSmarts(Reaction rxn) {
		return createReactionSmiles(rxn, MODE_INCLUDE_MAPPING | MODE_CREATE_SMARTS);
	}

	public static String createReactionSmiles(Reaction rxn) {
		return createReactionSmiles(rxn, MODE_INCLUDE_MAPPING);
	}

	public static String createReactionSmiles(Reaction rxn, int mode) {
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<rxn.getReactants(); i++) {
			if (i != 0)
				sb.append('.');
			sb.append(new IsomericSmilesCreator(rxn.getReactant(i), mode).getSmiles());
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
			sb.append(new IsomericSmilesCreator(rxn.getProduct(i), mode).getSmiles());
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
			addToGraph(new SmilesAtom(atom, -1, -1, false, false), graphIndex);
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
		// find closure neighbours of every atom and put them in canonical order. i.e. order of appearance in SMILES
		for (SmilesAtom smilesAtom:mGraphAtomList) {
			int closureCount = 0;
			for (int i=0; i<mMol.getConnAtoms(smilesAtom.atom); i++)
				if (!mBondUsed[mMol.getConnBond(smilesAtom.atom, i)])
					closureCount++;

			if (closureCount != 0) {
				smilesAtom.closureNeighbour = new int[closureCount];

				closureCount = 0;
				for (int i=0; i<mMol.getConnAtoms(smilesAtom.atom); i++) {
					if (!mBondUsed[mMol.getConnBond(smilesAtom.atom, i)]) {
						int neighbour = mMol.getConnAtom(smilesAtom.atom, i);
						smilesAtom.closureNeighbour[closureCount++] = (mSmilesIndex[neighbour] << 16) | neighbour;
					}
				}

				Arrays.sort(smilesAtom.closureNeighbour);
				for (int i=0; i<smilesAtom.closureNeighbour.length; i++)
					smilesAtom.closureNeighbour[i] = 0x0000FFFF & smilesAtom.closureNeighbour[i];
			}
		}

		// assign closure digits to closure bonds
		boolean[] closureNumberUsed = new boolean[mMol.getBonds()];
		mClosureNumber = new int[mMol.getBonds()];
		for (SmilesAtom smilesAtom:mGraphAtomList) {
			if (smilesAtom.closureNeighbour != null) {
				smilesAtom.closureOpens = new boolean[smilesAtom.closureNeighbour.length];
				for (int i=0; i<smilesAtom.closureNeighbour.length; i++) {
					for (int j=0; j<mMol.getConnAtoms(smilesAtom.atom); j++) {
						if (smilesAtom.closureNeighbour[i] == mMol.getConnAtom(smilesAtom.atom, j)) {
							int bond = mMol.getConnBond(smilesAtom.atom, j);
							if (!mBondUsed[bond]) { // opening closure bond
								mBondUsed[bond] = true;
								smilesAtom.closureOpens[i] = true;

								// assign and allocation closure digit
								mClosureNumber[bond] = 1;
								while (closureNumberUsed[mClosureNumber[bond]])
									mClosureNumber[bond]++;
								closureNumberUsed[mClosureNumber[bond]] = true;
							}
							else {  // closing closure bond: release closure digit
								closureNumberUsed[mClosureNumber[bond]] = false;
							}
						}
					}
				}
			}
		}
	}

	private void calculateEZBonds() {
		ArrayList<int[]> relativeBondParityList = new ArrayList<>();
		for (SmilesAtom currentSA:mGraphAtomList) {
			if (currentSA.parent != -1) {
				int ezBond = mMol.getBond(currentSA.atom, currentSA.parent);
				if (!mMol.isBINAPChiralityBond(ezBond)
				 && !mMol.isSmallRingBond(ezBond)
				 && (mMol.getBondParity(ezBond) == Molecule.cBondParityEor1
				  || mMol.getBondParity(ezBond) == Molecule.cBondParityZor2)) {
					SmilesAtom parentSA = mGraphAtomList.get(mSmilesIndex[currentSA.parent]);

					// Here we collect for every stereo double bond a list of relative halfParities,
					// which collects the relative bond symbol dependency of all connected single bonds.
					int[] bondWithHalfParity = new int[mMol.getConnAtoms(currentSA.atom)+mMol.getConnAtoms(parentSA.atom)-2];
					int halfParityIndex = 0;

					// halfParities translate 1-to-1 into the bond symbol ('/' or '\').
					// The value of halfParity1 is arbitrary, but for canonical SMILES its value must be reproducible.
					boolean parity = false;
					if (parentSA.parent != -1) {
						// If there is a bond leading to the first double bond atom, we take it as the reference.
						bondWithHalfParity[halfParityIndex++] = parentSA.bond;  // reference bond with halfParity flag unset
						}
					else {
						// If the graph starts with one of the double bond atoms (rare, but possible), then this atom has no parent atom.
						// Instead, it may be added in forward direction as a branch or it may be connected as a ring closure.
						// If that single bonded neighbour, which has the lower SMILES index, is effectively the second atom of the '/' or '\'
						// bond in the SMILES, we need to invert its halfParity.
						int firstNeighbourIndex = -1;
						int secondNeighbourIndex = -1;;
						int firstSmilesIndex = Integer.MAX_VALUE;
						for (int i=0; i<mMol.getConnAtoms(parentSA.atom); i++) {
							int connAtom = mMol.getConnAtom(parentSA.atom, i);
							if (connAtom != currentSA.atom) {
								if (firstNeighbourIndex == -1) {
									firstNeighbourIndex = i;
									firstSmilesIndex = mSmilesIndex[connAtom];
									}
								else {
									if (firstSmilesIndex < mSmilesIndex[connAtom]) {
										secondNeighbourIndex = i;
										}
									else {
										secondNeighbourIndex = firstNeighbourIndex;
										firstNeighbourIndex = i;
										}
									}
								}
							}
						if (secondNeighbourIndex == -1) {     // one neighbour atom
							int neighbourAtom = mMol.getConnAtom(parentSA.atom, firstNeighbourIndex);
							int neighbourBond = mMol.getConnBond(parentSA.atom, firstNeighbourIndex);
							bondWithHalfParity[halfParityIndex++] = neighbourBond | (isBondFromTo(parentSA.atom, neighbourAtom) ? 0x40000000 : 0);
							}
						else {      // two neighbour atoms
							int connAtom1 = mMol.getConnAtom(parentSA.atom, firstNeighbourIndex);
							int connBond1 = mMol.getConnBond(parentSA.atom, firstNeighbourIndex);
							int connAtom2 = mMol.getConnAtom(parentSA.atom, secondNeighbourIndex);
							int connBond2 = mMol.getConnBond(parentSA.atom, secondNeighbourIndex);
							bondWithHalfParity[halfParityIndex++] = connBond1 | (isBondFromTo(parentSA.atom, connAtom1) ? 0x40000000 : 0);
							bondWithHalfParity[halfParityIndex++] = connBond2 | (isBondFromTo(parentSA.atom, connAtom2) ? 0 : 0x40000000);
							}
						}

					if (mMol.getConnAtoms(parentSA.atom) == 3 && parentSA.parent != -1) {
						for (int i=0; i<mMol.getConnAtoms(parentSA.atom); i++) {
							int connAtom = mMol.getConnAtom(parentSA.atom, i);
							if (connAtom != parentSA.parent && connAtom != currentSA.atom) {
								// If we have a branch at the left double bond atom, we need to give it the same
								// half-parity given to the parent atom (other atom at the left double bond end),
								// because this is translated into the same sign ('/' or '\'), which
								// actually means on the different side of the double bond, because
								// "The 'visual interpretation' of the 'up-ness' or 'down-ness' of each single bond
								//  is relative to the carbon atom, not the double bond" (opensmiles.org).
								int branchBond = mMol.getConnBond(parentSA.atom, i);
								bondWithHalfParity[halfParityIndex++] = branchBond | (isBondFromTo(parentSA.atom, connAtom) ? 0x40000000 : 0);

								if (connAtom < parentSA.parent) // the other neighbour is the reference in OpenChemLib
									parity = !parity;	// invert
								break;
							}
						}
					}

					if (mMol.getBondParity(ezBond) == Molecule.cBondParityZor2)
						parity = !parity;

					for (int i=0; i<mMol.getConnAtoms(currentSA.atom); i++) {
						int childAtom = mMol.getConnAtom(currentSA.atom, i);
						if (childAtom != currentSA.parent) {
							boolean halfParity2 = parity;
							if (mMol.getConnAtoms(currentSA.atom) == 3) {
								for (int j=0; j<mMol.getConnAtoms(currentSA.atom); j++) {
									int connAtom = mMol.getConnAtom(currentSA.atom, j);
									if (connAtom != currentSA.parent && connAtom != childAtom) {
										if (connAtom < childAtom)
											halfParity2 = !halfParity2;
										break;
									}
								}
							}

							if (mMol.isBondParityPseudo(ezBond)) {
								int group = mCanonizer.getPseudoEZGroup(ezBond);
								if (!mPseudoStereoGroupInitialized[group]) {
									mPseudoStereoGroupInitialized[group] = true;
									mPseudoStereoGroupInversion[group] = halfParity2;
								}
								if (mPseudoStereoGroupInversion[group])
									halfParity2 = !halfParity2;
							}

						// If the graph continues normally and includes the childAtom down the line
						// or if we have a closure bond, which connects from current atom to child atom,
						// then we use the halfParity2 as generated.
						// If, however, the closure leads from the childAtom to the current atom, we
						// need to invert the halfParity2, since we put the '/' or '\' on the closing closure bond.
						int childBond = mMol.getBond(currentSA.atom, childAtom);
						bondWithHalfParity[halfParityIndex++] = childBond | (halfParity2 ^ isBondFromTo(currentSA.atom, childAtom) ? 0 : 0x40000000);
						}
					}

				relativeBondParityList.add(bondWithHalfParity);
				}
			}
		}

		mEZHalfParity = new int[mMol.getBonds()];
		if (relativeBondParityList.size() != 0)
			addRelativeBondHalfParities(relativeBondParityList.remove(0), false);

		while (relativeBondParityList.size() != 0) {
			int startSize = relativeBondParityList.size();
			for (int i=relativeBondParityList.size()-1; i>=0; i--) {
				int[] bondWithHalfParity = relativeBondParityList.get(i);
				int overlapCount = 0;
				boolean inverted = false;
				boolean collides = false;
				for (int bwhp:bondWithHalfParity) {
					int bond = bwhp & 0x3FFFFFFF;
					if (mEZHalfParity[bond] != 0) {
						boolean differs = ((bwhp & 0x40000000) != 0) ^ (mEZHalfParity[bond] == 2);
						if (overlapCount == 0)
							inverted = differs;
						else if (inverted != differs)
							collides = true;
						overlapCount++;
						}
					}
				if (overlapCount != 0) {
					bondWithHalfParity = relativeBondParityList.remove(i);
					// collisions, i.e. incompatible halfParity constraints should be very rare, but not impossible
					if (!collides)
						addRelativeBondHalfParities(bondWithHalfParity, inverted);
					}
				}

			// if we haven't found a constrained set, we just add the first one
			if (startSize == relativeBondParityList.size())
				addRelativeBondHalfParities(relativeBondParityList.remove(0), false);
			}
		}

	private void addRelativeBondHalfParities(int[] bondWithHalfParity, boolean inverted) {
		for (int bwhp:bondWithHalfParity) {
			mEZHalfParity[bwhp & 0x3FFFFFFF] = ((bwhp & 0x40000000) != 0) ^ inverted ? 2 : 1;
		}
	}

	private boolean isBondFromTo(int atom1, int atom2) {
		SmilesAtom sa1 = mGraphAtomList.get(mSmilesIndex[atom1]);
		if (sa1.parent == atom2)
			return false;
		SmilesAtom sa2 = mGraphAtomList.get(mSmilesIndex[atom2]);
		if (sa2.parent == atom1)
			return true;
		return sa2.isOpeningClosureTo(atom1);
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
			addToGraph(new SmilesAtom(atom, bond, parent,
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
		 && (mMode & MODE_KEKULIZED_OUTPUT) == 0
		 && (mMol.getAtomPi(atom) != 0
		  || (mMol.getAtomAbnormalValence(atom) == -1
		   && mMol.getAtomRadical(atom) == Molecule.cAtomRadicalStateNone)))
			label = label.toLowerCase();

		if (smilesAtom.isSideChainStart)
			builder.append('(');

		if (parent != -1)
			appendBondOrderSymbol(mMol.getBond(smilesAtom.atom, smilesAtom.parent), smilesAtom.parent, builder);

		int charge = mMol.getAtomCharge(atom);
		if (charge == 0 && (mMode & MODE_CREATE_SMARTS) != 0) {
			// Because SMARTS don't know a charge query feature, we set the atom charge in case of neg/pos charge required.
			// There is no way to require an uncharged atom (Molecule.cAtomQFNotChargeNeg | Molecule.cAtomQFNotChargePos).
			long chargeFeatures = mMol.getAtomQueryFeatures(atom) & Molecule.cAtomQFCharge;
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
				|| (mMol.isAromaticAtom(atom) && mMol.getAtomPi(atom) == 0 && (mMode & MODE_KEKULIZED_OUTPUT) == 0)
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
		boolean isLowerCase = mMol.isAromaticAtom(atom) && (mMode & MODE_KEKULIZED_OUTPUT) == 0;

		for (int a:atomicNo) {
			if (buffer.length() != 0)
				buffer.append(',');
			String label = Molecule.cAtomLabel[a];
				buffer.append(isLowerCase ? label.toLowerCase() : label);
			}

		return buffer.toString();
		}

	private String getAtomSMARTSFeatures(int atom, StringBuilder buffer) {
		buffer.setLength(0);

		long queryFeatures = mMol.getAtomQueryFeatures(atom);

		// SMARTS don't distinguish between a charged atom and atom charge query features
		int chargeFeatures = (int)((queryFeatures & Molecule.cAtomQFCharge) >> Molecule.cAtomQFChargeBits);
		switch (chargeFeatures) {
			case (int)((Molecule.cAtomQFNotChargeNeg | Molecule.cAtomQFNotChargePos) >> Molecule.cAtomQFChargeBits):
				buffer.append("+0");   // 'require negative charge' is translated into charge := -1
				break;
			case (int)((Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargePos) >> Molecule.cAtomQFChargeBits):
				if (mMol.getAtomCharge(atom) == 0)
					buffer.append("-");   // 'require negative charge' is translated into charge := -1
				break;
			case (int)((Molecule.cAtomQFNotCharge0 | Molecule.cAtomQFNotChargeNeg) >> Molecule.cAtomQFChargeBits):
				if (mMol.getAtomCharge(atom) == 0)
					buffer.append("+");   // 'require positive charge' is translated into charge := +1
				break;
			}

		long aromState = queryFeatures & Molecule.cAtomQFAromState;
		if (aromState == Molecule.cAtomQFAromatic)
			buffer.append(";a");
		else if (aromState == Molecule.cAtomQFNotAromatic)
			buffer.append(";A");

		long hydrogenQueryFeatures = queryFeatures & Molecule.cAtomQFHydrogen;
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
		long ringState = queryFeatures & Molecule.cAtomQFRingState;
		if (ringState == Molecule.cAtomQFNotChain)
			buffer.append(";!R0");
		else if (ringState == Molecule.cAtomQFNot2RingBonds)
			buffer.append(";!R1");
		else if (ringState == Molecule.cAtomQFNot3RingBonds)
			buffer.append(";!R2");
		else if (ringState == Molecule.cAtomQFNot4RingBonds)
			buffer.append(";!R3");
		else if (ringState == (Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds))
			buffer.append(";R0");
		else if (ringState == (Molecule.cAtomQFNotChain | Molecule.cAtomQFNot3RingBonds | Molecule.cAtomQFNot4RingBonds))
			buffer.append(";R1");
		else if (ringState == (Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot4RingBonds))
			buffer.append(";R2");
		else if (ringState == (Molecule.cAtomQFNotChain | Molecule.cAtomQFNot2RingBonds | Molecule.cAtomQFNot3RingBonds))
			buffer.append(";R3");

		long ringSize = queryFeatures & Molecule.cAtomQFNewRingSize;
		if (ringSize == Molecule.cAtomQFRingSize0)
			buffer.append(";!r" + ringSize);
		else if (ringSize == (Molecule.cAtomQFNewRingSize & ~Molecule.cAtomQFRingSize0))
			buffer.append(";r" + ringSize);
		else if (ringSize != 0) {
			if ((ringSize & Molecule.cAtomQFRingSizeLarge) != 0) {  // negative logic
				if ((ringSize & Molecule.cAtomQFRingSize0) == 0)
					buffer.append(";!r0" + ringSize);
				if ((ringSize & Molecule.cAtomQFRingSize3) == 0)
					buffer.append(";!r3" + ringSize);
				if ((ringSize & Molecule.cAtomQFRingSize4) == 0)
					buffer.append(";!r4" + ringSize);
				if ((ringSize & Molecule.cAtomQFRingSize5) == 0)
					buffer.append(";!r5" + ringSize);
				if ((ringSize & Molecule.cAtomQFRingSize6) == 0)
					buffer.append(";!r6" + ringSize);
				if ((ringSize & Molecule.cAtomQFRingSize7) == 0)
					buffer.append(";!r7" + ringSize);
			}
			else {
				buffer.append(";");
				if ((ringSize & Molecule.cAtomQFRingSize0) != 0)
					buffer.append("r0," + ringSize);
				if ((ringSize & Molecule.cAtomQFRingSize3) != 0)
					buffer.append("r3," + ringSize);
				if ((ringSize & Molecule.cAtomQFRingSize4) != 0)
					buffer.append("r4," + ringSize);
				if ((ringSize & Molecule.cAtomQFRingSize5) != 0)
					buffer.append("r5," + ringSize);
				if ((ringSize & Molecule.cAtomQFRingSize6) != 0)
					buffer.append("r6," + ringSize);
				if ((ringSize & Molecule.cAtomQFRingSize7) != 0)
					buffer.append("r7," + ringSize);
				buffer.setLength(buffer.length()-1);
			}
		}

		if (ringSize == 0) {
			ringSize = (queryFeatures & Molecule.cAtomQFSmallRingSize) >> Molecule.cAtomQFSmallRingSizeShift;
			if (ringSize != 0)
				buffer.append(";r" + ringSize);
		}

		long neighbourFeatures = queryFeatures & Molecule.cAtomQFNeighbours;
		if (neighbourFeatures == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot1Neighbour))
			buffer.append(";D1");   // exactly 1
		if (neighbourFeatures == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot2Neighbours))
			buffer.append(";D2");   // exactly 2
		if (neighbourFeatures == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot3Neighbours))
			buffer.append(";D3");   // exactly 3
		if (neighbourFeatures == (Molecule.cAtomQFNot3Neighbours | Molecule.cAtomQFNot4Neighbours))
			buffer.append(";!D3;!D4");   // less than 3
		if (neighbourFeatures == Molecule.cAtomQFNot4Neighbours)
			buffer.append(";!D4");   // less than 4
		if (neighbourFeatures == (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour))
			buffer.append(";!D0;!D1");   // more than 1
		if (neighbourFeatures == (Molecule.cAtomQFNot0Neighbours | Molecule.cAtomQFNot1Neighbour | Molecule.cAtomQFNot2Neighbours))
			buffer.append(";!D0;!D1;!D2");   // more than 2
		if (neighbourFeatures == (Molecule.cAtomQFNeighbours & ~Molecule.cAtomQFNot4Neighbours))
			buffer.append(";!D0;!D1;!D2;!D3");   // more than 3

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
		if (smilesAtom.closureNeighbour != null) {
			for (int i=0; i<smilesAtom.closureNeighbour.length; i++) {
				for (int j=0; j<mMol.getConnAtoms(smilesAtom.atom); j++) {
					if (smilesAtom.closureNeighbour[i] == mMol.getConnAtom(smilesAtom.atom, j)) {
						int bond = mMol.getConnBond(smilesAtom.atom, j);
						if (!smilesAtom.closureOpens[i])
							appendBondOrderSymbol(bond, smilesAtom.atom, builder);
						if (mClosureNumber[bond] > 9)
							builder.append('%');
						builder.append(mClosureNumber[bond]);
					}
				}
			}
		}
	}

	private void appendBondOrderSymbol(int bond, int parentAtom, StringBuilder builder) {
		int startLength = builder.length();
		if (mEZHalfParity[bond] != 0)
			builder.append(mEZHalfParity[bond] == 1 ? '/' : '\\');
		if (mMode == MODE_CREATE_SMARTS) {
			int bondQFTypes = mMol.getBondQueryFeatures(bond) & (Molecule.cBondQFBondTypes | Molecule.cBondQFRareBondTypes);
			if (bondQFTypes != 0) {
				if ((bondQFTypes & Molecule.cBondTypeSingle) != 0 && mEZHalfParity[bond] == 0) {
					builder.append('-');
					}
				if ((bondQFTypes & Molecule.cBondTypeDouble) != 0) {
					if (builder.length() != startLength)
						builder.append(',');
					builder.append('=');
					}
				if ((bondQFTypes & Molecule.cBondTypeTriple) != 0) {
					if (builder.length() != startLength)
						builder.append(',');
					builder.append('#');
					}
				if ((bondQFTypes & Molecule.cBondTypeQuadruple) != 0) {
					if (builder.length() != startLength)
						builder.append(',');
					builder.append('$');
					}
				if ((bondQFTypes & Molecule.cBondTypeQuintuple) != 0) {   // SMILES doesn't support quintuple bonds, thus we use quadruple
					if (builder.length() != startLength)
						builder.append(',');
					builder.append('$');
					}
				if ((bondQFTypes & Molecule.cBondTypeDelocalized) != 0) {
					if (builder.length() != startLength)
						builder.append(',');
					builder.append(':');
					}
				if ((bondQFTypes & Molecule.cBondTypeMetalLigand) != 0) {
					if (builder.length() != startLength)
						builder.append(',');
					builder.append(mMol.isMetalAtom(parentAtom) ? "<-" : "->");
					}
				}
			}
		if (startLength == builder.length() && (!mMol.isAromaticBond(bond) || (mMode & MODE_KEKULIZED_OUTPUT) != 0)) {
			int bondType = mMol.getBondType(bond) & Molecule.cBondTypeMaskSimple;
			if (bondType == Molecule.cBondTypeSingle) {
				if (mMol.isAromaticAtom(mMol.getBondAtom(0, bond))
				 && mMol.isAromaticAtom(mMol.getBondAtom(1, bond))
				 && (mMode & MODE_KEKULIZED_OUTPUT) == 0
				 && mEZHalfParity[bond] == 0)
					builder.append('-');
			} else if (bondType == Molecule.cBondTypeDouble)
				builder.append('=');
			else if (bondType == Molecule.cBondTypeTriple)
				builder.append('#');
			else if (bondType == Molecule.cBondTypeQuadruple)
				builder.append('$');
			else if (bondType == Molecule.cBondTypeQuintuple)	// not supported by SMILES, we use quadruple instead
				builder.append('$');
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
				 && ((mSmilesIndex[neighbour[0]] < mSmilesIndex[neighbour[1]]) ^ (neighbour[0] < neighbour[1])))
					inversion = !inversion;
				}
			}
		else {
			int[] neighborAtom = new int[4];
			int[] neighborRank = new int[4];
			int index = 0;

			if (parent != -1) {
				neighborAtom[index] = parent;
				neighborRank[index++] = 8 * mSmilesIndex[parent];
				}
			if (mMol.getImplicitHydrogens(atom) != 0) {
				neighborAtom[index] = Integer.MAX_VALUE; // index must be higher than any non-hydrogen neighbor
				neighborRank[index++] = 8 * mSmilesIndex[atom]; // smaller than any potential closure neighbor
				}
			else if (mMol.getConnAtoms(atom) == 3) {    // rare 3-neighbour stereo center as with S,P,N
				neighborAtom[index] = Integer.MAX_VALUE; // index must be higher than any non-hydrogen neighbor
				neighborRank[index++] = 8 * mSmilesIndex[atom]; // smaller than any potential closure neighbor
				}
			for (int i=0; i<mMol.getConnAtoms(atom); i++) {
				int connAtom = mMol.getConnAtom(atom, i);
				if (connAtom != parent) {
					neighborAtom[index] = connAtom;
					neighborRank[index++] = getSmilesRank(atom, i);
					}
				}

			inversion = isInverseOrder(neighborAtom, neighborRank);
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

	private boolean isInverseOrder(int[] neighborAtom, int[] neighborRank) {
		boolean inversion = false;
		for (int i=1; i<4; i++) {
			for (int j=0; j<i; j++) {
				if (neighborAtom[j] > neighborAtom[i])
					inversion = !inversion;
				if (neighborRank[j] > neighborRank[i])
					inversion = !inversion;
				}
			}
		return inversion;
		}

	/**
	 * @param atom for which to return a neighbor's smiles rank
	 * @param neighbourIndex index for getConnAtoms() to get neighbor atom and bond
	 * @return neighbor's position rank in smiles from a perspective of atom using closure digit positions for closure neighbors
	 */
	private int getSmilesRank(int atom, int neighbourIndex) {
		int bond = mMol.getConnBond(atom, neighbourIndex);
		int neighbour = mMol.getConnAtom(atom, neighbourIndex);
		if (mClosureNumber[bond] != 0) {
			// if neighbour is attached via a closure digit, then the rank is based primarily on atom's position
			// in the smiles and secondary on the count of other closures at atom that precede this closure
			int rank = 8 * mSmilesIndex[atom] + 1;
			int[] closureNeighbour = mGraphAtomList.get(mSmilesIndex[atom]).closureNeighbour;
			for (int i=0; i<closureNeighbour.length && neighbour != closureNeighbour[i]; i++)
				rank++;
			return rank;
			}

		// if the neighbour is not a closure return a rank based on its atom index in the smiles
		return 8 * mSmilesIndex[neighbour];
	}

	private boolean isOrganic(int atomicNo) {
		return (atomicNo >= 5 && atomicNo <= 9)     // B,C,N,O,F
			|| (atomicNo >= 15 && atomicNo <= 17)   // P,S,Cl
			|| atomicNo == 35                       // Br
			|| atomicNo == 53;                      // I
	}
}

class SmilesAtom {
	public int atom,bond,parent;
	public boolean isSideChainStart,isSideChainEnd;
	public int[] closureNeighbour;
	public boolean[] closureOpens;

	public SmilesAtom(int atom, int bond, int parent, boolean isSideChainStart, boolean isSideChainEnd) {
		this.atom = atom;
		this.bond = bond;
		this.parent = parent;
		this.isSideChainStart = isSideChainStart;
		this.isSideChainEnd = isSideChainEnd;
	}

	public boolean isOpeningClosureTo(int atom) {
		if (closureNeighbour != null)
			for (int i=0; i<closureNeighbour.length; i++)
				if (atom == closureNeighbour[i] && closureOpens[i])
					return true;

		return false;
	}
}