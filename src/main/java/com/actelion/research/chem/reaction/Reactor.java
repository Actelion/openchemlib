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

package com.actelion.research.chem.reaction;

import com.actelion.research.chem.*;
import com.actelion.research.chem.coords.CoordinateInventor;
import com.actelion.research.util.IntArrayComparator;

import java.util.ArrayList;
import java.util.Arrays;

public class Reactor {
	public static final int MODE_RETAIN_COORDINATES = 1;
	public static final int MODE_FULLY_MAP_REACTIONS = 2;
	public static final int MODE_REMOVE_DUPLICATE_PRODUCTS = 4;
	public static final int MODE_ALLOW_CHARGE_CORRECTIONS = 8;

	private final Reaction		mGenericReaction;
	private final SSSearcher	mSSSearcher;
	private	final StereoMolecule[] mReactant;
	private final int[][]		mMinFreeValence;	// minimum required free valence on reactant atoms
	private final boolean[][]	mIsReactionCenter;	// reaction center flags on product atoms
	private final boolean		mRetainCoordinates,mFullyMapReactions,mUniqueProductsOnly,mAllowChargeCorrections;
	private final int			mMaxCount;
	private int					mMaxGenericMapNo,mReactantMatchCombinationCount;
	private final ArrayList<int[]>[] mMatchList;
	private int[][][]			mReactantMapNo;	// Reactant mapNos of the real reactant change with every reactant substructure match
	private StereoMolecule[][]	mProduct;
	private String[][]			mIDCode;
	private boolean[]			mIsDuplicate;

	/**
	 * Constructs a Reactor that is able to construct products from a generic reaction (transformation)
	 * and a list of real world reactants. Reactors built with this constructor will not use product
	 * atom coordinates from the generic products, nor will they be able to generate fully mapped reactions.
	 * They will produce and return all unique products or reactions.
	 * @param reaction generic reactions consisting of reactant substructures with optional query features
	 */
	public Reactor(Reaction reaction) {
		this(reaction, false, false, Integer.MAX_VALUE, true);
		}

	/**
	 * Constructs a Reactor that is able to construct products from a generic reaction (transformation)
	 * and a list of real world reactants. Reactors built with this constructor will be able to generate
	 * fully mapped reactions, They will produce and return all unique products or reactions.
	 * @param reaction generic reactions consisting of reactant substructures with optional query features
	 * @param retainCoordinates if true, then atom coordinates from the generic products are taken into the real world products
	 */
	public Reactor(Reaction reaction, boolean retainCoordinates) {
		this(reaction, (retainCoordinates ? MODE_RETAIN_COORDINATES : 0) + MODE_REMOVE_DUPLICATE_PRODUCTS, Integer.MAX_VALUE);
		}

	/**
	 * Constructs a Reactor that is able to construct products from a generic reaction (transformation)
	 * and a list of real world reactants. The reactor may or may not use atoms coordinates from the generic
	 * products, when generating clean product atom coordinates.
	 * @param reaction generic reactions consisting of reactant substructures with optional query features
	 * @param retainCoordinates if true, then atom coordinates from the generic products are taken into the real world products
	 * @param fullyMapReactions if true, then real world reactants and products will have valid mapping numbers after product generation
	 * @param maxProducts maximum number of products/reactions to be enumerated
	 * @param uniqueOnly whether to skip duplicate products/reactions because of starting material symmetry
	 */
	@Deprecated // Use: Reactor(Reaction reaction, int mode, int maxProducts)
	public Reactor(Reaction reaction, boolean retainCoordinates, boolean fullyMapReactions, int maxProducts, boolean uniqueOnly) {
		this(reaction, (retainCoordinates ? MODE_RETAIN_COORDINATES : 0)
						+ (fullyMapReactions ? MODE_FULLY_MAP_REACTIONS : 0)
						+ (uniqueOnly ? MODE_REMOVE_DUPLICATE_PRODUCTS : 0), maxProducts);
	}

	/**
	 * Constructs a Reactor that is able to build products from a generic reaction (transformation)
	 * and a list of real world reactants. These MODE_xxx flags can be used to finetune the Reactor's
	 * behaviour:<br>
	 * MODE_RETAIN_COORDINATES: If set, then atom coordinates from the generic products are taken into the
	 * real world products as fixed atom coordinates set when creating atom coordinates for all product atoms.<br>
	 * MODE_FULLY_MAP_REACTIONS: If true, then real world reactants and products will have valid mapping numbers
	 * after product generation. getFullyMappedReactions() can be used after each setReactant() call to get
	 * complete and mapped Reaction objects, one for every potential product.<br>
	 * MODE_REMOVE_DUPLICATE_PRODUCTS: Whether to skip duplicate products/reactions because of starting material
	 * symmetry.<br>
	 * MODE_ALLOW_CHARGE_CORRECTIONS: If set, then the Reactor will try to neutralize real world reactant atoms,
	 * if the matching generic reactant atom has no charge. E.g. this would allow carboxylate anions to represent
	 * a carboxylic acid in an esterification despite formally the anion oxygen has no free valence to form a new bond.
	 * WARNING: This may produce products with unbalanced atom charges.<br>
	 * @param reaction generic reaction with reactant & product substructures with optional query features (reactant)
	 * @param mode 0 or any combination of MODE_xxx flags
	 * @param maxProducts limitation of products/reactions to be enumerated in case of heavy symmetries
	 */
	public Reactor(Reaction reaction, int mode, int maxProducts) {
		mRetainCoordinates = (mode & MODE_RETAIN_COORDINATES) != 0;
		mFullyMapReactions = (mode & MODE_FULLY_MAP_REACTIONS) != 0;
		mUniqueProductsOnly = (mode & MODE_REMOVE_DUPLICATE_PRODUCTS) != 0;
		mAllowChargeCorrections = (mode & MODE_ALLOW_CHARGE_CORRECTIONS) != 0;
		mMaxCount = maxProducts;

		mGenericReaction = reaction;
		mReactant = new StereoMolecule[reaction.getReactants()];

		// for sub-structure-search all generic reactants must be fragments
		for (int i=0; i<reaction.getReactants(); i++) {
			reaction.getReactant(i).setFragment(true);
			reaction.getReactant(i).ensureHelperArrays(Molecule.cHelperParities);
			}

		for (int i=0; i<reaction.getProducts(); i++)
			reaction.getProduct(i).ensureHelperArrays(Molecule.cHelperParities);

		removeMappingErrors(reaction);

					// calculate minimum free valence of reactant atoms
		mMinFreeValence = new int[reaction.getReactants()][];
		for (int i=0; i<reaction.getReactants(); i++) {
			StereoMolecule reactant = reaction.getReactant(i);
			mMinFreeValence[i] = new int[reactant.getAtoms()];
			for (int j=0; j<reactant.getAtoms(); j++) {
				int mapNo = reactant.getAtomMapNo(j);
				if (mapNo != 0) {
					for (int k=0; k<reaction.getProducts(); k++) {
						StereoMolecule product = reaction.getProduct(k);
						for (int l=0; l<product.getAtoms(); l++) {
							if (product.getAtomMapNo(l) == mapNo) {
								int dif = reactant.getFreeValence(j) - product.getFreeValence(l);
								mMinFreeValence[i][j] = Math.max(dif, 0);
								}
							}
						}
					}
				}
			}

		// Find and mark all reaction center atoms in product.
		// These are mapped atoms that have a different number of neighbours in reactant or product
		// or different neighbours when judging based on mapping numbers or atomic numbers.
		mIsReactionCenter = new boolean[reaction.getProducts()][];
		for (int i=0; i<reaction.getProducts(); i++) {
			StereoMolecule product = reaction.getProduct(i);
			mIsReactionCenter[i] = new boolean[product.getAtoms()];
			for (int j=0; j<product.getAtoms(); j++) {
				int mapNo = product.getAtomMapNo(j);
				if (mapNo != 0) {
					for (int k=0; k<reaction.getReactants(); k++) {
						StereoMolecule reactant = reaction.getReactant(k);
						for (int l=0; l<reactant.getAtoms(); l++) {
							if (reactant.getAtomMapNo(l) == mapNo) {
								if (new IntArrayComparator().compare(
										getSortedConnectedMapNos(product, j),
										getSortedConnectedMapNos(reactant, l)) != 0
								 || new IntArrayComparator().compare(
										getSortedConnectedAtomicNos(product, j),
										getSortedConnectedAtomicNos(reactant, l)) != 0)
									mIsReactionCenter[i][j] = true;
								}
							}
						}
					}
				}
			}

		mSSSearcher = new SSSearcher();
		mMatchList = new ArrayList[mReactant.length];
		}

	private int[] getSortedConnectedAtomicNos(StereoMolecule mol, int atom) {
		int[] atomicNo = new int[mol.getConnAtoms(atom) - mol.getExcludedNeighbourCount(atom)];
		int index = 0;
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int connAtom = mol.getConnAtom(atom, i);
			if ((mol.getAtomQueryFeatures(connAtom) & Molecule.cAtomQFExcludeGroup) == 0)
				atomicNo[index++] = mol.getAtomicNo(connAtom);
			}
		Arrays.sort(atomicNo);
		return atomicNo;
		}

	private int[] getSortedConnectedMapNos(StereoMolecule mol, int atom) {
		int[] mapNo = new int[mol.getConnAtoms(atom) - mol.getExcludedNeighbourCount(atom)];
		int index = 0;
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int connAtom = mol.getConnAtom(atom, i);
			if ((mol.getAtomQueryFeatures(connAtom) & Molecule.cAtomQFExcludeGroup) == 0)
				mapNo[index++] = mol.getAtomMapNo(connAtom);
			}
		Arrays.sort(mapNo);
		return mapNo;
		}

	/**
	 * Removes unmatched mapping numbers and mapping numbers on exclude atoms
	 * @param rxn
	 */
	private void removeMappingErrors(Reaction rxn) {
		mMaxGenericMapNo = 0;
		for (int i=0; i<rxn.getMolecules(); i++) {
			StereoMolecule mol = mGenericReaction.getMolecule(i);
			for (int atom=0; atom<mol.getAtoms(); atom++) {
				if ((mol.getAtomQueryFeatures(atom) & Molecule.cAtomQFExcludeGroup) != 0)
					mol.setAtomMapNo(atom, 0, false);
				else if (mMaxGenericMapNo < mol.getAtomMapNo(atom))
					mMaxGenericMapNo = mol.getAtomMapNo(atom);
				}
			}

		boolean[] used = new boolean[mMaxGenericMapNo+1];
		for (int i=0; i<rxn.getReactants(); i++) {
			StereoMolecule reactant = rxn.getReactant(i);
			for (int atom=0; atom<reactant.getAtoms(); atom++) {
				int mapNo = reactant.getAtomMapNo(atom);
				if (mapNo != 0)
					used[mapNo] = true;
				}
			}

		for (int i=0; i<rxn.getProducts(); i++) {
			StereoMolecule product = rxn.getProduct(i);
			for (int atom=0; atom<product.getAtoms(); atom++) {
				int mapNo = product.getAtomMapNo(atom);
				if (mapNo != 0) {
					if (!used[mapNo])
						product.setAtomMapNo(atom, 0, false);
					else
						used[mapNo] = false;
					}
				}
			}

		for (int i=0; i<rxn.getReactants(); i++) {
			StereoMolecule reactant = rxn.getReactant(i);
			for (int atom=0; atom<reactant.getAtoms(); atom++) {
				int mapNo = reactant.getAtomMapNo(atom);
				if (mapNo != 0 && used[mapNo])
					reactant.setAtomMapNo(atom, 0, false);
				}
			}
		}

	/**
	 * Defines a real world reactant. This triggers a sub-structure search of the respective
	 * generic reactant in the real reactant resulting in potentially multiple matches.
	 * Matches are removed, if free valences are not compatible with the formation of new bonds.
	 * @param no
	 * @param reactant with correctly set parity flags
	 * @return whether the generic reactant was found in this real world reactant
	 */
	@SuppressWarnings("unchecked")
	public boolean setReactant(int no, StereoMolecule reactant) {
		mReactant[no] = reactant;
		mProduct = null;
		StereoMolecule genericReactant = mGenericReaction.getReactant(no);

		mSSSearcher.setMol(genericReactant, mReactant[no]);
		int matchMode = SSSearcher.cDefaultMatchMode + (mAllowChargeCorrections ? 0 : SSSearcher.cMatchAtomCharge);
		if (mSSSearcher.findFragmentInMolecule(SSSearcher.cCountModeRigorous, matchMode) == 0) {
			mMatchList[no] = new ArrayList<>();
			mReactantMatchCombinationCount = 0;
			return false;
			}

		// eliminate matches where reaction would exceed an atom valence
		mMatchList[no] = mSSSearcher.getMatchList();
		for (int j=mMatchList[no].size()-1; j>=0; j--) {
			int[] matchingAtom = mMatchList[no].get(j);
			for (int k=0; k<matchingAtom.length; k++) {
				if (matchingAtom[k] != -1) {
					// If mAllowChargeCorrections==true, then we could neutralize atom charges, to be more forgiving:
					// We don't use the restrictive SSSearcher.cMatchAtomCharge in SSS matching above and we remove
					// the charge correction for real world reactant atoms when performing the free valence check.
					int correction = 0;
					if (mAllowChargeCorrections
					 && genericReactant.getAtomCharge(k) == 0
					 && mReactant[no].getAtomCharge(matchingAtom[k]) != 0) {
						int valence = mReactant[no].getMaxValenceUncharged(matchingAtom[k]);
						correction = mReactant[no].getElectronValenceCorrection(matchingAtom[k], valence);
						}

					if (mMinFreeValence[no][k] > 0
					 && mMinFreeValence[no][k] > mReactant[no].getFreeValence(matchingAtom[k]) - correction) {
						mMatchList[no].remove(j);
						break;
						}
					}
				}
			}

		mReactantMatchCombinationCount = 1;
		for (int i=0; i<mReactant.length; i++)
			mReactantMatchCombinationCount *= (mMatchList[i] == null) ? 0 : mMatchList[i].size();

		return true;
		}


	public StereoMolecule[][] getProducts() {
		ensureProducts();

		StereoMolecule[][] uniqueProduct = new StereoMolecule[getProductCount()][];
		int count = 0;
		for (int i=0; i<mReactantMatchCombinationCount && count<mMaxCount; i++)
			if (!mUniqueProductsOnly || !mIsDuplicate[i])
				uniqueProduct[count++] = mProduct[i];

		return uniqueProduct;
		}


	public String[][] getProductIDCodes() {
		ensureProducts();

		String[][] uniqueIDCode = new String[getProductCount()][];
		int count = 0;
		for (int i=0; i<mReactantMatchCombinationCount && count<mMaxCount; i++)
			if (!mUniqueProductsOnly || !mIsDuplicate[i])
				uniqueIDCode[count++] = mIDCode[i];

		return uniqueIDCode;
		}


	private int getProductCount() {
		int count = mReactantMatchCombinationCount;
		if (mUniqueProductsOnly)
			for (int i=0; i<mReactantMatchCombinationCount; i++)
				if (mIsDuplicate[i])
					count--;

		if (mMaxCount > 0 && count > mMaxCount)
			count = mMaxCount;

		return count;
		}


	/**
	 * After instantiating the Reactor with MODE_FULLY_MAP_REACTIONS and after supplying real reactants,
	 * this method may be used to construct the fully mapped reaction. If one or more reactants have multiple
	 * matches of their generic reactants, then multiple reactions leading to isomeric products are returned.
	 * After calling this method, one or more new real reactants may provided with setReactant()
	 * before calling this method again.
	 * @return array with one or more reactions
	 */
	public Reaction[] getFullyMappedReactions() {
		ensureProducts();

		Reaction[] reaction = new Reaction[getProductCount()];

		int count = 0;
		int index = 0;
		int[] matchIndex = new int[mReactant.length];
		do {
			if (!mIsDuplicate[index] || !mUniqueProductsOnly) {
				reaction[count] = new Reaction();
				for (int i=0; i<mReactant.length; i++) {
					StereoMolecule reactant = mReactant[i].getCompactCopy();
					for (int atom=0; atom<mReactantMapNo[i][matchIndex[i]].length; atom++)
						reactant.setAtomMapNo(atom, mReactantMapNo[i][matchIndex[i]][atom], false);
					reaction[count].addReactant(reactant);
					}

				for (StereoMolecule product: mProduct[index])
					reaction[count].addProduct(product);

				count++;
				}
			index++;
			} while (nextMatchListCombination(matchIndex) && count < mMaxCount);

		return reaction;
		}


	private void ensureProducts() {
		if (mProduct != null)
			return;

		if (mFullyMapReactions) {
			int firstFreeMapNo = mMaxGenericMapNo+1;
			mReactantMapNo = new int[mReactant.length][][];

			for (int i=0; i<mReactant.length; i++) {
				mReactantMapNo[i] = new int[mMatchList[i].size()][];
				for (int j=0; j<mMatchList[i].size(); j++)
					mReactantMapNo[i][j] = getReactantMapNos(i, mMatchList[i].get(j), firstFreeMapNo);

				firstFreeMapNo += mReactant[i].getAtoms();

				// atoms, which match any (non-exclude) atom of generic reactant, don't need a new mapping number
				StereoMolecule gr = mGenericReaction.getReactant(i);
				for (int atom=0; atom<gr.getAtoms(); atom++)
					if ((gr.getAtomQueryFeatures(atom) & Molecule.cAtomQFExcludeGroup) == 0)
						firstFreeMapNo--;
				}
			}

		mProduct = new StereoMolecule[mReactantMatchCombinationCount][mGenericReaction.getProducts()];
		mIDCode = new String[mReactantMatchCombinationCount][mGenericReaction.getProducts()];
		mIsDuplicate = new boolean[mReactantMatchCombinationCount];

		SortedStringList uniqueList = new SortedStringList();
		StringBuilder uniqueBuilder = new StringBuilder();

		int index = 0;
		int count = 0;
		if (mReactantMatchCombinationCount != 0) {
			int[] matchListIndex = new int[mReactant.length];
			do {
				uniqueBuilder.setLength(0);
				for (int i=0; i<mGenericReaction.getProducts(); i++) {
					mProduct[index][i] = generateProduct(mMatchList, matchListIndex, i);
					mIDCode[index][i] = (new Canonizer(mProduct[index][i]).getIDCode());
					uniqueBuilder.append(mIDCode[index][i]);
					}

				mIsDuplicate[index] = (uniqueList.addString(uniqueBuilder.toString()) == -1);

				if (!mUniqueProductsOnly || !mIsDuplicate[index])
					count++;

				index++;
				} while (nextMatchListCombination(matchListIndex) && count < mMaxCount);
			}
		}


	private boolean nextMatchListCombination(int[] matchListIndex) {
		for (int i=0; i<matchListIndex.length; i++) {
			if (matchListIndex[i] < mMatchList[i].size() - 1) {
				matchListIndex[i]++;
				return true;
				}
			matchListIndex[i] = 0;
			}
		return false;
		}


	/**
	 * Constructs an array of real reactant atoms mapping numbers matchin atom indexes.
	 * Copies generic reactant's mapping numbers first by considering the matchList.
	 * Then assigns new unused mapping numbers to the remaining reactant atoms.
	 * @param reactant
	 * @param matchList
	 * @param firstMapNo
	 * @return
	 */
	private int[] getReactantMapNos(int reactant, int[] matchList, int firstMapNo) {
		int[] reactantMapNo = new int[mReactant[reactant].getAtoms()];
		Arrays.fill(reactantMapNo, -1);

		StereoMolecule genericReactant = mGenericReaction.getReactant(reactant);
		for (int atom=0; atom<genericReactant.getAtoms(); atom++)
			if (matchList[atom] != -1)
				reactantMapNo[matchList[atom]] = genericReactant.getAtomMapNo(atom);

		for (int i=0; i<reactantMapNo.length; i++)
			if (reactantMapNo[i] == -1)
				reactantMapNo[i] = firstMapNo++;

		return reactantMapNo;
		}


	private StereoMolecule generateProduct(ArrayList<int[]>[] matchList, int[] matchListIndex, int genericProductNo) {
		StereoMolecule genericProduct = mGenericReaction.getProduct(genericProductNo);

		StereoMolecule product = new StereoMolecule();

		int esrGroupCountAND = 0;
		int esrGroupCountOR = 0;
		ArrayList<Integer> racemicAtomList = null;

		for (int i=0; i<mReactant.length; i++) {
			StereoMolecule genericReactant = mGenericReaction.getReactant(i);
			mReactant[i].ensureHelperArrays(Molecule.cHelperNeighbours);
			int[] matchingAtom = matchList[i].get(matchListIndex[i]);
			int[] mapNo = new int[mReactant[i].getAtoms()];
			boolean[] excludeAtom = new boolean[mReactant[i].getAtoms()];
			boolean[] excludeBond = new boolean[mReactant[i].getBonds()];

			// Exclude atoms from real reactants, which exist in generic reactant and are not mapped in generic reaction
			// (including attached bonds)
			for (int j=0; j<genericReactant.getAtoms(); j++) {
				if (matchingAtom[j] != -1) {	// non-exclude-group atoms only
					if (genericReactant.getAtomMapNo(j) == 0) {
						int excludedAtom = matchingAtom[j];
						excludeAtom[excludedAtom] = true;
						for (int k=0; k<mReactant[i].getConnAtoms(excludedAtom); k++)
							excludeBond[mReactant[i].getConnBond(excludedAtom, k)] = true;
						// To cover the rare case where we have metal bonded neighbours in generic and real reactants:
						for (int k=mReactant[i].getAllConnAtoms(excludedAtom); k<mReactant[i].getAllConnAtomsPlusMetalBonds(excludedAtom); k++)
							excludeBond[mReactant[i].getConnBond(excludedAtom, k)] = true;
						}
					else {
						mapNo[matchingAtom[j]] = genericReactant.getAtomMapNo(j);
						}
					}
				}

			// Exclude bonds from real reactants, which connect mapped atoms in generic reaction
			for (int j=0; j<genericReactant.getBonds(); j++) {
				int bondAtom1 = genericReactant.getBondAtom(0, j);
				int bondAtom2 = genericReactant.getBondAtom(1, j);
				if (genericReactant.getAtomMapNo(bondAtom1) != 0
				 && genericReactant.getAtomMapNo(bondAtom2) != 0) {
					int atom1 = matchingAtom[bondAtom1];
					int atom2 = matchingAtom[bondAtom2];
					if (atom1 != -1 && atom2 != -1) {
						for (int k=0; k<mReactant[i].getBonds(); k++) {
							if ((mReactant[i].getBondAtom(0, k) == atom1
							  && mReactant[i].getBondAtom(1, k) == atom2)
							 || (mReactant[i].getBondAtom(0, k) == atom2
							  && mReactant[i].getBondAtom(1, k) == atom1)) {
								excludeBond[k] = true;
								break;
								}
							}
						}
					}
				}

			int[] newAtomNo = new int[mReactant[i].getAtoms()];

			for (int j=0; j<mReactant[i].getAtoms(); j++) {
				if (!excludeAtom[j]) {
					newAtomNo[j] = mReactant[i].copyAtom(product, j, esrGroupCountAND, esrGroupCountOR);

					if (mFullyMapReactions)
						product.setAtomMapNo(newAtomNo[j], mReactantMapNo[i][matchListIndex[i]][j], false);

					if (mapNo[j] != 0) {  // take charge and radical from generic product atoms
						product.setAtomMapNo(newAtomNo[j], mapNo[j], false);
						for (int k=0; k<genericProduct.getAtoms(); k++) {
							if (genericProduct.getAtomMapNo(k) == mapNo[j]) {
								product.setAtomCharge(newAtomNo[j], genericProduct.getAtomCharge(k));
								product.setAtomRadical(newAtomNo[j], genericProduct.getAtomRadical(k));
								if (mIsReactionCenter[genericProductNo][k]) {
									int parity = mReactant[i].getAtomParity(j);
									long productParityHint = genericProduct.getAtomQueryFeatures(k) & Molecule.cAtomQFRxnParityHint;
									if (productParityHint == 0) {
										product.setAtomParity(newAtomNo[j], Molecule.cAtomParityUnknown, false);
										}
									else if (productParityHint == Molecule.cAtomQFRxnParityRacemize) {
										if (racemicAtomList == null)
											racemicAtomList = new ArrayList<>();
										racemicAtomList.add(newAtomNo[j]);
										}
									else {
										if (parity == Molecule.cAtomParity1 || parity == Molecule.cAtomParity2) {
											// Excluded atoms are added later to the product. Therefore, the parity needs
											// to be inverted for every atom index that is between the current and the atom's new index
											// Without further potential correction of final order of here excluded atoms
											// this only works for substitution reactions of one atom for a new one!!!
											boolean inversion = (productParityHint == Molecule.cAtomQFRxnParityInvert);
											for (int l=0; l<mReactant[i].getConnAtoms(j); l++) {
												int connAtom = mReactant[i].getConnAtom(j, l);
												if (excludeAtom[connAtom]) {
													for (int m=0; m<mReactant[i].getConnAtoms(j); m++) {
														int otherConnAtom = mReactant[i].getConnAtom(j, m);
														if (!excludeAtom[otherConnAtom] && otherConnAtom > connAtom)
															inversion = !inversion;
														}
													}
												}
											if (inversion)
												parity = (parity == Molecule.cAtomParity1) ? Molecule.cAtomParity2 : Molecule.cAtomParity1;
											product.setAtomParity(newAtomNo[j], parity, false);
											}
										else {
											product.setAtomParity(newAtomNo[j], 0, false);
											}
										}
									}
								break;
								}
							}
						}
					}
				}

			// Delocalized double bonds are copied as delocalized to avoid cumulated double bonds
			// when delocalized double bonds from generic product and from real reactant touch.
			for (int j=0; j<mReactant[i].getBonds(); j++)
				if (!excludeBond[j])
					mReactant[i].copyBond(product, j, esrGroupCountAND, esrGroupCountOR, newAtomNo, true);

			esrGroupCountAND = product.renumberESRGroups(Molecule.cESRTypeAnd);
			esrGroupCountOR = product.renumberESRGroups(Molecule.cESRTypeOr);
			}

		// copy all unmapped atoms of generic product
		// and setup newAtomNo array for all(!!!) atoms of generic product
		int[] newAtomNo = new int[genericProduct.getAtoms()];

		for (int j=0; j<genericProduct.getAtoms(); j++) {
			int mapNo = genericProduct.getAtomMapNo(j);
			if (mapNo == 0) {
				newAtomNo[j] = genericProduct.copyAtom(product, j, esrGroupCountAND, esrGroupCountOR);
				}
			else {
				for (int k=0; k<product.getAllAtoms(); k++) {
					if (product.getAtomMapNo(k) == mapNo) {
						newAtomNo[j] = k;
						break;
						}
					}
				}
			}

		// mark atoms of generic product to retain coordinates when creating new ones
		if (mRetainCoordinates) {
			for (int j=0; j<genericProduct.getAtoms(); j++) {
				product.setAtomMarker(newAtomNo[j], true);
				product.setAtomX(newAtomNo[j], genericProduct.getAtomX(j));
				product.setAtomY(newAtomNo[j], genericProduct.getAtomY(j));
				}
			}

		// copy all bonds of generic product, except for bridge bonds
		for (int gpBond=0; gpBond<genericProduct.getBonds(); gpBond++) {
			if (genericProduct.isBondBridge(gpBond))
				continue;

			int productBond = genericProduct.copyBond(product, gpBond, esrGroupCountAND, esrGroupCountOR, newAtomNo, false);

			// If we have multiple allowed bond types on generic product bond, then we need to calculate
			// the new bond order from reactant bond order and generic reaction bond order change.
			int productQFBondType = genericProduct.getBondQueryFeatures(gpBond) & Molecule.cBondQFBondTypes;
			if (productQFBondType != 0) {
				boolean found = false;
				product.setBondQueryFeature(productBond, Molecule.cBondQFBondTypes, false);
				int gpMapNo1 = genericProduct.getAtomMapNo(genericProduct.getBondAtom(0, gpBond));
				int gpMapNo2 = genericProduct.getAtomMapNo(genericProduct.getBondAtom(1, gpBond));
				for (int i=0; !found && i<mGenericReaction.getReactants(); i++) {
					StereoMolecule genericReactant = mGenericReaction.getReactant(i);
					for (int grBond=0; !found && grBond<genericReactant.getBonds(); grBond++) {
						int grAtom1 = genericReactant.getBondAtom(0, grBond);
						int grAtom2 = genericReactant.getBondAtom(1, grBond);
						int grMapNo1 = genericReactant.getAtomMapNo(grAtom1);
						int grMapNo2 = genericReactant.getAtomMapNo(grAtom2);
						if (grMapNo1 != 0 && grMapNo2 != 0
						 && ((grMapNo1 == gpMapNo1 && grMapNo2 == gpMapNo2) || (grMapNo1 == gpMapNo2 && grMapNo2 == gpMapNo1))) {
							found = true;

							int reactantQFBondType = genericReactant.getBondQueryFeatures(grBond) & Molecule.cBondQFBondTypes;
							if (reactantQFBondType != 0) {
								int[] matchingAtom = matchList[i].get(matchListIndex[i]);
								int rAtom1 = matchingAtom[grAtom1];
								int rAtom2 = matchingAtom[grAtom2];
								int rBond = mReactant[i].getBond(rAtom1, rAtom2);
								int rBondType = mReactant[i].getBondType(rBond);
								int rBondOrder = mReactant[i].getBondOrder(rBond);

								// only consider simple bond order features
								reactantQFBondType &= Molecule.cBondTypeSingle | Molecule.cBondTypeDouble | Molecule.cBondTypeTriple;
								productQFBondType &= Molecule.cBondTypeSingle | Molecule.cBondTypeDouble | Molecule.cBondTypeTriple;

								// increase in bond order
								if (reactantQFBondType == (Molecule.cBondTypeSingle | Molecule.cBondTypeDouble)
								 && productQFBondType == (Molecule.cBondTypeDouble | Molecule.cBondTypeTriple)) {
									product.setBondType(productBond, rBondOrder <= 1 ? Molecule.cBondTypeDouble : Molecule.cBondTypeTriple);
									}
								// decrease in bond order
								else if (reactantQFBondType == (Molecule.cBondTypeDouble | Molecule.cBondTypeTriple)
								 && productQFBondType == (Molecule.cBondTypeSingle | Molecule.cBondTypeDouble)) {
									product.setBondType(productBond, rBondOrder == 3 ? Molecule.cBondTypeDouble : Molecule.cBondTypeSingle);
									}
								else {
									product.setBondType(productBond, rBondType);
									}
								}
							}
						}
					}
				}
			}

		product.ensureHelperArrays(Molecule.cHelperNeighbours);

		// copy corrected atom parities of generic product reaction center atoms
		boolean esrGroupsChanged = false;
		for (int j=0; j<genericProduct.getAtoms(); j++) {
			// copy atom parities of unmapped generic product atoms and potentially convert them according to new neighbour order
			if (genericProduct.getAtomMapNo(j) == 0) {
				int parity = genericProduct.translateTHParity(j, newAtomNo);
				product.setAtomParity(newAtomNo[j], parity, false);
				}
			// copy corrected atom parities of generic product reaction center atoms
			else if (mIsReactionCenter[genericProductNo][j]
				  && copyConfigurationFromGenericProduct(genericProduct, j, product, newAtomNo[j])) {
				int parity = genericProduct.translateTHParity(j, newAtomNo);
				product.setAtomParity(newAtomNo[j], parity, false);
				if (parity == Molecule.cAtomParity1
						|| parity == Molecule.cAtomParity2) {
					int esrType = genericProduct.getAtomESRType(j);
					int esrGroup = genericProduct.getAtomESRGroup(j);
					if (esrType == Molecule.cESRTypeAnd)
						esrGroup += esrGroupCountAND;
					else if (esrType == Molecule.cESRTypeOr)
						esrGroup += esrGroupCountOR;

					product.setAtomESR(newAtomNo[j], esrType, esrGroup);
					esrGroupsChanged = true;
					}
				}
			}

		if (esrGroupsChanged) {
			esrGroupCountAND = product.renumberESRGroups(Molecule.cESRTypeAnd);
			esrGroupCountOR = product.renumberESRGroups(Molecule.cESRTypeOr);
			}

		// copy pseudo parities of generic product if there are no parity hints
		boolean[] pseudoParityHandled = new boolean[genericProduct.getAllAtoms()];
		for (int j=0; j<genericProduct.getAtoms(); j++)
			if (mIsReactionCenter[genericProductNo][j]
			 && !pseudoParityHandled[j]
			 && genericProduct.isAtomParityPseudo(j)
			 && (genericProduct.getAtomQueryFeatures(j) & Molecule.cAtomQFRxnParityHint) == 0
			 && genericProduct.getConnAtoms(j) == product.getConnAtoms(newAtomNo[j]))
				copyPseudoParitiesFromGenericProduct(genericProduct, j, product, newAtomNo, pseudoParityHandled);

		// delete all fragments from product which are not connected to generic product
		boolean[] includeAtom = new boolean[product.getAllAtoms()];
		for (int j : newAtomNo) includeAtom[j] = true;
		boolean found = true;
		while (found) {
			found = false;
			for (int bond=0; bond<product.getAllBonds(); bond++) {
				int atom1 = product.getBondAtom(0, bond);
				int atom2 = product.getBondAtom(1, bond);
				if (includeAtom[atom1] && !includeAtom[atom2]) {
					includeAtom[atom2] = true;
					found = true;
					}
				else if (includeAtom[atom2] && !includeAtom[atom1]) {
					includeAtom[atom1] = true;
					found = true;
					}
				}
			}
		for (int atom=0; atom<product.getAllAtoms(); atom++)
			product.setAtomSelection(atom, !includeAtom[atom]);
		product.deleteSelectedAtoms();

		if (!mFullyMapReactions)
			for (int atom=0; atom<product.getAllAtoms(); atom++)
				product.setAtomMapNo(atom, 0, false);

		// Here we need to convert delocalized double bonds from real reactant to nicely match
		// alternating delocalized double bonds from generic product.
		new AromaticityResolver(product).locateDelocalizedDoubleBonds(null);

		product.setParitiesValid(0);

		int mode = CoordinateInventor.MODE_REMOVE_HYDROGEN
				 | (mRetainCoordinates ? CoordinateInventor.MODE_PREFER_MARKED_ATOM_COORDS : 0);
		new CoordinateInventor(mode).invent(product);

		if (racemicAtomList != null) {
			product.ensureHelperArrays(Molecule.cHelperParities);
			for (int atom:racemicAtomList) {
				if (product.isAtomStereoCenter(atom)) {
					product.setAtomParity(atom, Molecule.cAtomParity1, false);
					product.setAtomESR(atom, Molecule.cESRTypeAnd, esrGroupCountAND++);
					}
				}
			}

		return product;
		}

//	private boolean reactantQueryFeatureMatchStereoCenter(int reactantMapNo) {
//		for (int i=0; i<mGenericReaction.getReactants(); i++) {
//			StereoMolecule reactant = mGenericReaction.getReactant(i);
//			for (int atom=0; atom<reactant.getAtoms(); atom++)
//				if (reactant.getAtomMapNo(atom) == reactantMapNo)
//					return (reactant.getAtomQueryFeatures(atom) & Molecule.cAtomQFMatchStereo) != 0;
//			}
//		return false;
//		}

	private void copyPseudoParitiesFromGenericProduct(StereoMolecule genericProduct, int firstPseudoParityAtom, StereoMolecule product, int[] newAtomNo, boolean[] pseudoParityHandled) {
		int[] relatedAtom = new int[genericProduct.getAtoms()];
		int relatedAtomCount = findRelatedPseudoParityAtoms(genericProduct, firstPseudoParityAtom, relatedAtom);

		boolean matchingAbsParityFound = false;
		boolean invertedAbsParityFound = false;
		ArrayList<int[]> pseudoParityList = new ArrayList<>();

		for (int i=0; i<relatedAtomCount; i++) {
			if (genericProduct.isAtomParityPseudo(relatedAtom[i])) {
				int pseudoParity = genericProduct.translateTHParity(relatedAtom[i], newAtomNo);
				if (pseudoParity == Molecule.cAtomParity1
				 || pseudoParity == Molecule.cAtomParity2) {
					// Real product parities may not be pseudo if the real reactant has additional features.
					// These are given and cannot be inverted because of pseudo parity requirements in the generic product.
					// We check for compatibility and only apply relative pseudo parities, if possible.
					int parity = product.getAtomParity(newAtomNo[relatedAtom[i]]);
					if (parity == Molecule.cAtomParity1 || parity == Molecule.cAtomParity2) {
						if (pseudoParity == parity)
							matchingAbsParityFound = true;
						else
							invertedAbsParityFound = true;
						}
					}

				// We collect a list of all pseudo parity atoms with their translated pseudo parities.
				// If where we have existing real parities on the same atoms and if these
				// - either match the required pseudo parities, then we can add pseudo parities to atoms without a real parity
				// - or match the all inverted pseudo parities, then we can add inverted pseudo parities to atoms without a real parity
				int[] ps = new int[2];
				ps[0] = relatedAtom[i];
				ps[1] = pseudoParity;
				pseudoParityList.add(ps);
				}
			}

		if (!matchingAbsParityFound || !invertedAbsParityFound) {
			for (int[] pseudoParity : pseudoParityList) {
				int ps = pseudoParity[1];
				if ((ps == Molecule.cAtomParity1 || ps == Molecule.cAtomParity2) && invertedAbsParityFound)
					ps = (ps == Molecule.cAtomParity1) ? Molecule.cAtomParity2 : Molecule.cAtomParity1;
				product.setAtomParity(newAtomNo[pseudoParity[0]], ps, true);
				}
			}
		}

	private int findRelatedPseudoParityAtoms(StereoMolecule genericProduct, int firstAtom, int[] relatedAtom) {
		boolean[] isRelated = new boolean[genericProduct.getAtoms()];
		relatedAtom[0] = firstAtom;
		int relatedAtoms = 1;
		isRelated[firstAtom] = true;
		for (int current=0; current<relatedAtoms; current++) {
			for (int i=0; i<genericProduct.getConnAtoms(relatedAtom[current]); i++) {
				int connBond = genericProduct.getConnBond(relatedAtom[current],i);
				if (genericProduct.isRingBond(connBond) || genericProduct.getBondOrder(connBond) == 2 || genericProduct.isBINAPChiralityBond(connBond)) {
					int connAtom = genericProduct.getConnAtom(relatedAtom[current],i);
					if (!isRelated[connAtom]) {
						relatedAtom[relatedAtoms++] = connAtom;
						isRelated[connAtom] = true;
						}
					}
				}
			}
		return relatedAtoms;
		}

	private boolean copyConfigurationFromGenericProduct(StereoMolecule genericProduct, int genericProductAtom, StereoMolecule product, int productAtom) {
		// If we have a hint on the generic product atom, then we have already taken the generic reactant atom parity.
		// Otherwise, if we have a stereo center on the generic product atom, we may consider copying that...
		if (genericProduct.getAtomParity(genericProductAtom) != Molecule.cAtomParityNone
		 && !genericProduct.isAtomParityPseudo(genericProductAtom)
		 && (genericProduct.getAtomQueryFeatures(genericProductAtom) & Molecule.cAtomQFRxnParityHint) == 0) {
			int mapNo = genericProduct.getAtomMapNo(genericProductAtom);
			if (mapNo != 0) {
				for (int i=0; i<mGenericReaction.getReactants(); i++) {
					StereoMolecule genericReactant = mGenericReaction.getReactant(i);
					for (int reactantAtom=0; reactantAtom<genericReactant.getAtoms(); reactantAtom++) {
						if (genericReactant.getAtomMapNo(reactantAtom) == mapNo) {
							// If the generic reactant atom is not even a pseudo stereo center,
							// then we assume that any product configuration is intentional and correct
							if (genericProduct.getConnAtoms(genericProductAtom) == product.getConnAtoms(productAtom))
								return true;

							// If the reactant atom is explicitly set to match the stereo center,
							// we assume that the product also contains a well specified stereo center to be copied.
							if ((genericReactant.getAtomQueryFeatures(reactantAtom) & Molecule.cAtomQFMatchStereo) != 0)
								return true;

							// If we have an explicit change of the ESR type, we also assume that this is intentional
							// and copy the stereo configuration from the product.
							if (genericReactant.getAtomESRType(reactantAtom) != genericProduct.getAtomESRType(genericProductAtom))
								return true;

							return false;
							}
						}
					}
				}
			}
		return false;
		}
	}
