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

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;
import java.util.Arrays;

import static com.actelion.research.chem.SSSearcher.cCountModeRigorous;
import static com.actelion.research.chem.SSSearcher.cDefaultMatchMode;

/**
 * The SRSearcher class handles reaction-sub-structure searches. Correctly, the class
 * should be named SuperReactionSearcher, because it is rather a search for super reactions
 * for a given query reaction. The query reaction may also be called transformation and
 * may contain atom or bond based query features.
 */
public class SRSearcher {
	private StereoMolecule mQueryReactantBuffer,mQueryProductBuffer,mReactantBuffer,mProductBuffer;
	private StereoMolecule mQueryReactant,mQueryProduct,mReactant,mProduct;
	private SSSearcher mReactantSearcher,mProductSearcher;
	private boolean mQueryIsPreprocessed,mReactionIsPreprocessed;
	private int mQueryMaxMapNo,mMaxMapNo;
	private byte[] mQueryCode,mQueryCoords,mQueryMapping,mReactionCode,mReactionCoords,mReactionMapping;
	private long[] mQueryReactantFFP,mQueryProductFFP,mReactantFFP,mProductFFP;
	private int[] mQueryReactantToProductAtom,mQueryReactantToProductBond,mReactantToProductAtom,mReactantToProductBond;
	private int[] mProductMatch,mQueryNeighborDelta,mNeighborDelta;

	public SRSearcher() {
		mReactantSearcher = new SSSearcher() {
			@Override public boolean areAtomsSimilar(int moleculeAtom, int fragmentAtom) {
				return super.areAtomsSimilar(moleculeAtom, fragmentAtom) && productAtomMatches(moleculeAtom, fragmentAtom);
				}
			@Override public boolean areBondsSimilar(int moleculeBond, int fragmentBond) {
				return super.areBondsSimilar(moleculeBond, fragmentBond) && productBondMatches(moleculeBond, fragmentBond);
				}
			};
		mProductSearcher = new SSSearcher() {
			@Override public boolean areAtomsSimilar(int moleculeAtom, int fragmentAtom) {
				return (mProductMatch == null || mProductMatch[fragmentAtom] == -1 || mProductMatch[fragmentAtom] == moleculeAtom)
					&& super.areAtomsSimilar(moleculeAtom, fragmentAtom);
				}
			};
		}

	public void setQuery(byte[] rxncode, byte[] rxnmapping, byte[] rxncoords, long[] reactantFFP, long[] productFFP) {
		mQueryCode = rxncode;
		mQueryMapping = rxnmapping;
		mQueryCoords = rxncoords;
		mQueryIsPreprocessed = false;
		mQueryReactant = null;
		mQueryReactantFFP = reactantFFP;
		mQueryProduct = null;
		mQueryProductFFP = productFFP;
		}

	public void setReaction(byte[] rxncode, byte[] rxnmapping, byte[] rxncoords, long[] reactantFFP, long[] productFFP) {
		mReactionCode = rxncode;
		mReactionMapping = rxnmapping;
		mReactionCoords = rxncoords;
		mReactionIsPreprocessed = false;
		mReactant = null;
		mReactantFFP = reactantFFP;
		mProduct = null;
		mProductFFP = productFFP;
		}

	/**
	 * This defines the query reaction (or transformation).
	 * Typically, this method is called once, while setReaction() is called many times,
	 * if a reaction collection is searched for hits. For acceleration through ffp based
	 * pre-screening, you should use this method to supply query ffps.
	 * If the query reaction contains multiple reactants or multiple products,
	 * then setQuery() merges these into one reactant and one product molecule.
	 * If you call setQuery() with the same query reactions multiple times, then
	 * for a maximum of performance you should cache merged query reactions and pass these.
	 * Merging can be done by getMergedCopy() of the reaction class.
	 * @param query
	 * @param reactantFFP
	 * @param productFFP
	 */
	public void setQuery(Reaction query, long[] reactantFFP, long[] productFFP) {
		mQueryCode = null;
		mQueryReactantFFP = reactantFFP;
		mQueryProductFFP = productFFP;
		mQueryIsPreprocessed = false;
		if (query == null || query.getReactants() == 0 || query.getProducts() == 0) {
			mQueryReactant = null;
			mQueryProduct = null;
			return;
			}

		splitQuery(query);
		}

	/**
	 * This defines the query reaction (or transformation).
	 * Typically, this method is called once, while setReaction() is called many times,
	 * if a reaction collection is searched for hits. For acceleration through ffp based
	 * pre-screening, you should use this method to supply query ffps.
	 * If the query reaction contains multiple reactants or multiple products,
	 * these are merged into one molecule each.
	 * Thus, for a maximum of performance you may avoid this step by parsing a reaction
	 * that contains one reactant and one product only.
	 * @param reaction
	 * @param reactantFFP
	 * @param productFFP
	 */
	public void setReaction(Reaction reaction, long[] reactantFFP, long[] productFFP) {
		mReactionCode = null;
		mReactantFFP = reactantFFP;
		mProductFFP = productFFP;
		mReactionIsPreprocessed = false;
		if (reaction == null || reaction.getReactants() == 0 || reaction.getProducts() == 0) {
			mReactant = null;
			mProduct = null;
			return;
		}

		splitReaction(reaction);
	}

	/**
	 * This defines the query reaction (or transformation).
	 * Typically, this method is called once, while setReaction() is called many times,
	 * if a reaction collection is searched for hits.
	 * If the query reaction contains multiple reactants or multiple products,
	 * these are merged into one molecule each.
	 * Thus, for a maximum of performance you may avoid this step by parsing a reaction
	 * that contains one reactant and one product only.
	 * @param query
	 */
	public void setQuery(Reaction query) {
		mQueryCode = null;
		mQueryReactantFFP = null;
		mQueryProductFFP = null;
		mQueryIsPreprocessed = false;
		if (query == null || query.getReactants() == 0 || query.getProducts() == 0) {
			mQueryReactant = null;
			mQueryProduct = null;
			return;
			}

//		if (!query.isPerfectlyMapped())
//			return;

		splitQuery(query);
		}

	public void setReaction(Reaction reaction) {
		mReactionCode = null;
		mReactantFFP = null;
		mProductFFP = null;
		mReactionIsPreprocessed = false;
		if (reaction == null || reaction.getReactants() == 0 || reaction.getProducts() == 0) {
			mReactant = null;
			mProduct = null;
			return;
			}

		splitReaction(reaction);
		}

	public void stop() {
		mReactantSearcher.stop();
		mProductSearcher.stop();
		}

	private void preprocessQuery() {
		if (!mQueryIsPreprocessed) {
			mQueryMaxMapNo = getHighestMapNo(mQueryReactant, mQueryProduct);
			mReactantSearcher.setFragment(mQueryReactant);
			mProductSearcher.setFragment(mQueryProduct);

			if (mQueryReactant != null && mQueryProduct != null) {
				mQueryReactantToProductAtom = createReactantToProductAtomMap(mQueryReactant, mQueryProduct, mQueryMaxMapNo);
				mQueryReactantToProductBond = createReactantToProductBondMap(mQueryReactant, mQueryProduct, mQueryReactantToProductAtom);
				mQueryNeighborDelta = createMappedAtomNeighborDeltas(mQueryReactant, mQueryProduct, mQueryReactantToProductAtom);
				}

			mQueryIsPreprocessed = true;
			}
		}

	private void preprocessReaction() {
		if (!mReactionIsPreprocessed) {
			mMaxMapNo = getHighestMapNo(mReactant, mProduct);
			mReactantSearcher.setMolecule(mReactant);
			mProductSearcher.setMolecule(mProduct);

			if (mReactant != null && mProduct != null) {
				mReactantToProductAtom = createReactantToProductAtomMap(mReactant, mProduct, mMaxMapNo);
				mReactantToProductBond = createReactantToProductBondMap(mReactant, mProduct, mReactantToProductAtom);
				mNeighborDelta = createMappedAtomNeighborDeltas(mReactant, mProduct, mReactantToProductAtom);
				}

			mReactionIsPreprocessed = true;
			}
		}

	private int getHighestMapNo(StereoMolecule reactant, StereoMolecule product) {
		int maxMapNo = 0;
		for (int atom=0; atom<reactant.getAllAtoms(); atom++)
			maxMapNo = Math.max(maxMapNo, reactant.getAtomMapNo(atom));
		for (int atom=0; atom<product.getAllAtoms(); atom++)
			maxMapNo = Math.max(maxMapNo, product.getAtomMapNo(atom));
		return maxMapNo;
		}

	/**
	 * Performs a reaction substructure search with the SSSearcher.cDefaultMatchMode
	 * @return
	 */
	public boolean isQueryInReaction() {
		return isQueryInReaction(cDefaultMatchMode);
		}

	/**
	 * @param matchMode cDefaultMatchMode or combination of cMatchAtomCharge, cMatchAtomMass, cMatchDBondToDelocalized, cMatchAromDBondToDelocalized
	 * @return
	 */
	public boolean isQueryInReaction(int matchMode) {
		if (mQueryReactantFFP != null && mQueryProductFFP != null && mReactantFFP != null && mProductFFP != null) {
			for (int i=0; i<mReactantFFP.length; i++)
				if ((mQueryReactantFFP[i] & ~mReactantFFP[i]) != 0)
					return false;
			for (int i=0; i<mProductFFP.length; i++)
				if ((mQueryProductFFP[i] & ~mProductFFP[i]) != 0)
					return false;
			}

		if (mQueryReactant == null && mQueryProduct == null)
			splitQuery(ReactionEncoder.decode(mQueryCode, mQueryMapping, mQueryCoords, null, null, false));

		if (mReactant == null && mProduct == null)
			splitReaction(ReactionEncoder.decode(mReactionCode, mReactionMapping, mReactionCoords, null, null, false));

		if (mQueryReactant == null || mQueryProduct == null || mReactant == null || mProduct == null)
			return false;

		preprocessQuery();
		preprocessReaction();

		mProductMatch = null;	// must be null for reactent side search

		mProductSearcher.setupAtomAndBondFeatures(matchMode);

		int count = mReactantSearcher.findFragmentInMolecule(cCountModeRigorous, matchMode);
		if (count == 0)
			return false;

		mProductMatch = new int[mQueryProduct.getAllAtoms()];

		// For every substructure match on the reactant side, we need to check, whether we find a corresponding substructure match
		// on the product side. For that we create a product match map from the reactant substructure match by applying the atom
		// mapping from query and reaction. The substructure search on the product side only allows query atoms to match, if they
		// are unmapped or correspond to the reactant match via atom mapping.
		ArrayList<int[]> matchList = mReactantSearcher.getMatchList();
		for (int[] match:matchList) {
			Arrays.fill(mProductMatch, -1);
			for (int i=0; i<match.length; i++)
				if (mQueryReactantToProductAtom[i] != -1)
					mProductMatch[mQueryReactantToProductAtom[i]] = mReactantToProductAtom[match[i]];

			if (mProductSearcher.isFragmentInMolecule())
				return true;
			}

		return false;
		}

	private void splitQuery(Reaction query) {
		if (query == null) {
			mQueryReactant = null;
			mQueryProduct = null;
			return;
			}

		if (query.getReactants() == 1) {
			mQueryReactant = query.getReactant(0);
			}
		else if (query.getReactants() > 1) {
			if (mQueryReactantBuffer == null)
				mQueryReactantBuffer = new StereoMolecule();
			mQueryReactant = mQueryReactantBuffer;

			query.getReactant(0).copyMolecule(mQueryReactant);
			for (int i=1; i<query.getReactants(); i++)
				mQueryReactant.addMolecule(query.getReactant(i));
			}

		if (query.getProducts() == 1) {
			mQueryProduct = query.getProduct(0);
			}
		else {
			if (mQueryProductBuffer == null)
				mQueryProductBuffer = new StereoMolecule();
			mQueryProduct = mQueryProductBuffer;

			query.getProduct(0).copyMolecule(mQueryProduct);
			for (int i=1; i<query.getProducts(); i++)
				mQueryProduct.addMolecule(query.getProduct(i));
			}

		mQueryReactant.setFragment(true);
		mQueryProduct.setFragment(true);
		}

	private void splitReaction(Reaction reaction) {
		if (reaction == null) {
			mReactant = null;
			mProduct = null;
			return;
			}

		if (reaction.getReactants() == 1) {
			mReactant = reaction.getReactant(0);
			}
		else if (reaction.getReactants() > 1) {
			if (mReactantBuffer == null)
				mReactantBuffer = new StereoMolecule();
			mReactant = mReactantBuffer;

			reaction.getReactant(0).copyMolecule(mReactant);
			for (int i=1; i<reaction.getReactants(); i++)
				mReactant.addMolecule(reaction.getReactant(i));
			}

		if (reaction.getProducts() == 1) {
			mProduct = reaction.getProduct(0);
			}
		else {
			if (mProductBuffer == null)
				mProductBuffer = new StereoMolecule();
			mProduct = mProductBuffer;

			reaction.getProduct(0).copyMolecule(mProduct);
			for (int i=1; i<reaction.getProducts(); i++)
				mProduct.addMolecule(reaction.getProduct(i));
			}
		}

	/**
	 * Assuming that the product doesn't have duplicate mapping numbers, this method creates a lookup array
	 * to get product atom indexes from reactant atom indexes.
	 */
	private int[] createReactantToProductAtomMap(StereoMolecule reactant, StereoMolecule product, int maxMapNo) {
		int[] mapNoToProduct = new int[maxMapNo + 1];
		for (int atom = 0; atom < product.getAllAtoms(); atom++)
			mapNoToProduct[product.getAtomMapNo(atom)] = atom;

		int[] reactantToProductAtom = new int[reactant.getAllAtoms()];
		for (int atom = 0; atom < reactant.getAllAtoms(); atom++) {
			int mapNo = reactant.getAtomMapNo(atom);
			reactantToProductAtom[atom] = (mapNo != 0) ? mapNoToProduct[mapNo] : -1;
			}

		return reactantToProductAtom;
		}

	/**
	 * Assuming that the product doesn't have duplicate mapping numbers, this method creates a lookup array
	 * to get product bond indexes from reactant bond indexes.
	 */
	private int[] createReactantToProductBondMap(StereoMolecule reactant, StereoMolecule product, int[] reactantToProductAtom) {
		int[] reactantToProductBond = new int[reactant.getAllBonds()];
		product.ensureHelperArrays(Molecule.cHelperNeighbours);
		for (int bond=0; bond<reactant.getAllBonds(); bond++) {
			int atom1 = reactantToProductAtom[reactant.getBondAtom(0, bond)];
			int atom2 = reactantToProductAtom[reactant.getBondAtom(1, bond)];
			reactantToProductBond[bond] = (atom1 != -1 && atom2 != -1) ? product.getBond(atom1, atom2) : -1;
			}
		return reactantToProductBond;
		}


	private int[] createMappedAtomNeighborDeltas(StereoMolecule reactant, StereoMolecule product, int[]reactantToProductAtom) {
		int[] neighborDelta = new int[reactantToProductAtom.length];

		for (int atom=0; atom<reactantToProductAtom.length; atom++)
			if (reactantToProductAtom[atom] != -1)
				neighborDelta[atom] = getNonExcludedNeighbours(product, reactantToProductAtom[atom]) - getNonExcludedNeighbours(reactant, atom);

		return neighborDelta;
		}

	private int getNonExcludedNeighbours(StereoMolecule mol, int atom) {
		int count = 0;
		for (int i=0; i<mol.getConnAtoms(atom); i++)
			if ((mol.getAtomQueryFeatures(mol.getConnAtom(atom, i)) & Molecule.cAtomQFExcludeGroup) == 0)
				count++;
		return count;
		}

	private boolean productAtomMatches(int reactantAtom, int reactantQueryAtom) {
		// if the query atom is mapped than we only allow mapped reaction atoms
		if (mQueryReactant.getAtomMapNo(reactantQueryAtom) != 0
		 && mReactant.getAtomMapNo(reactantAtom) == 0)
			return false;

		if (mNeighborDelta[reactantAtom] != mQueryNeighborDelta[reactantQueryAtom])
			return false;

		int productAtom = mReactantToProductAtom[reactantAtom];
		int productQueryAtom = mQueryReactantToProductAtom[reactantQueryAtom];

		if (productAtom == -1 || productQueryAtom == -1)
			return true;

		return mProductSearcher.areAtomsSimilar(productAtom, productQueryAtom);
		}

	private boolean productBondMatches(int reactantBond, int reactantQueryBond) {
		int productBond = mReactantToProductBond[reactantBond];
		int productQueryBond = mQueryReactantToProductBond[reactantQueryBond];

		if (productBond == -1 && productQueryBond != -1)
			return false;

		if (productBond == -1 || productQueryBond == -1)
			return true;

		return mProductSearcher.areBondsSimilar(productBond, productQueryBond);
		}
	}
