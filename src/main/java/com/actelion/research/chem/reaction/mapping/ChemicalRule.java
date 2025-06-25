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

package com.actelion.research.chem.reaction.mapping;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;
import com.actelion.research.util.SortedList;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * A ChemicalRule is basically a chemical reaction (transformation) defined by a reaction substructure
 * and a product substructure with full stoichiometry and completely mapped atoms.
 * ChemicalRules are used the following way:<br>
 * - a substructure search locates all matches of the rule's reactant structure in the to-be-mapped reactant<br>
 * - for every match the transformation of the ChemicalRule is applied to the to-be-mapped reactant<br>
 * - the modified reactant is similarity-graph-mapped with the original product and scored<br>
 * - the score of the unmodified reaction is compared to all modified reaction scores considering the ChemicalRule's score delta<br>
 * - the best scoring mapping is taken as final mapping<br>
 * In order to facilitate an efficient application of the rule to any query reactant, the ChemicalRule object
 * maintains an array of ChemicalRuleBond objects, which describe those bonds that need to be changed,
 * created, or broken in the query reactant. A
 */
public class ChemicalRule {
	private final String mName;
	private final boolean mIsStoichiometrical;
	private final float mPanalty;
	private final StereoMolecule mReactant,mProduct;
	private final ChemicalRuleBond[] mRuleBonds;
	private final ArrayList<RuleIncomingFragment> mRuleFragmentList;
	private final int[] mReactantAtomsForRemoval,mMapNoToReactantAtom;
	private int[] mInvertedTHParity,mReactantAtomSymmetryConstraint;

	public ChemicalRule(String name, String idcode) {
		mName = name;

		Reaction rxn = ReactionEncoder.decode(idcode, false);
		if (rxn.getReactants() != 1 || rxn.getProducts() != 1)
			System.out.println("ERROR: Rule '"+mName+"' doesn't contain exactly one reactant and one product!");
		if (!rxn.isFragment())
			System.out.println("ERROR: Rule '"+mName+"' reactant and product are not marked as fragment!");
		mReactant = rxn.getReactant(0);
		mProduct = rxn.getProduct(0);
		mReactant.ensureHelperArrays(Molecule.cHelperNeighbours);
		mProduct.ensureHelperArrays(Molecule.cHelperSymmetrySimple);

		// key: lower bondAtom mapNo, higher bondAtom mapNo
		// value: reactantAtom1, reactantAtom2, productBondOrder

		SortedList<ChemicalRuleBond> bondList = new SortedList<>();

		mMapNoToReactantAtom = new int[mReactant.getAtoms()+1];
		mMapNoToReactantAtom[0] = -1;    // to cause exceptions in case of faulty logic
		for (int atom=0; atom<mReactant.getAtoms(); atom++)
			if (mReactant.getAtomMapNo(atom) != 0)
				mMapNoToReactantAtom[mReactant.getAtomMapNo(atom)] = atom;

		calculateReactantAtomSymmetryConstraints(mMapNoToReactantAtom);

		boolean[] reactantBondFoundInProduct = new boolean[mReactant.getBonds()];
		for (int productBond=0; productBond<mProduct.getBonds(); productBond++) {
			int productAtom1 = mProduct.getBondAtom(0, productBond);
			int productAtom2 = mProduct.getBondAtom(1, productBond);
			int mapNo1 = mProduct.getAtomMapNo(productAtom1);
			int mapNo2 = mProduct.getAtomMapNo(productAtom2);
			if (mapNo1 != 0 && mapNo2 != 0) {   // unmapped new atoms and exclude atoms (not mapped) don't go into the bond list
				int atom1 = mMapNoToReactantAtom[mapNo1];
				int atom2 = mMapNoToReactantAtom[mapNo2];
				int productBondType = mProduct.getBondTypeSimple(productBond);
				int reactantBond = mReactant.getBond(atom1, atom2);
				if (reactantBond == -1) {
					bondList.add(new ChemicalRuleBond(atom1, atom2, mapNo1, mapNo2, productBondType));
					}
				else {
					// Note: we don't support quadruple or quintuple bonds here!
					if ((mReactant.getBondQueryFeatures(reactantBond) & Molecule.cBondQFBondTypes) == 0) {
						int reactantBondType = mReactant.getBondTypeSimple(reactantBond);
						if (reactantBondType != productBondType)
							bondList.add(new ChemicalRuleBond(atom1, atom2, mapNo1, mapNo2, productBondType));
						}
					reactantBondFoundInProduct[reactantBond] = true;
					}
				}
			}

		for (int bond=0; bond<mReactant.getBonds(); bond++) {
			if (!reactantBondFoundInProduct[bond]) {
				int atom1 = mReactant.getBondAtom(0, bond);
				int atom2 = mReactant.getBondAtom(1, bond);
				if (!mReactant.isExcludeGroupAtom(atom1)
				 && !mReactant.isExcludeGroupAtom(atom2)) {
					int mapNo1 = mReactant.getAtomMapNo(atom1);
					int mapNo2 = mReactant.getAtomMapNo(atom2);
					if (mapNo1 != 0 || mapNo2 != 0)   // exclude atoms are not mapped and don't go into the bond list
						bondList.add(new ChemicalRuleBond(atom1, atom2, mapNo1, mapNo2, ChemicalRuleBond.BOND_TYPE_DELETE));
					}
				}
			}

		mRuleBonds = bondList.toArray(new ChemicalRuleBond[0]);

		// Collect all non-mapped disappearing atoms in rule's reactant
		int count = 0;
		for (int atom=0; atom<mReactant.getAtoms(); atom++)
			if (mReactant.getAtomMapNo(atom) == 0 && !mReactant.isExcludeGroupAtom(atom))
				count++;
		mReactantAtomsForRemoval = new int[count];
		count = 0;
		for (int atom=0; atom<mReactant.getAtoms(); atom++)
			if (mReactant.getAtomMapNo(atom) == 0 && !mReactant.isExcludeGroupAtom(atom))
				mReactantAtomsForRemoval[count++] = atom;

		// List all non-mapped incoming fragments in rule's product
		mRuleFragmentList = new ArrayList<>();
		boolean[] isHandledFragmentAtom = new boolean[mProduct.getAtoms()];
		for (int atom=0; atom<mProduct.getAtoms(); atom++)
			if (mProduct.getAtomMapNo(atom) == 0
			 && !mProduct.isExcludeGroupAtom(atom)
			 && !isHandledFragmentAtom[atom])
				mRuleFragmentList.add(new RuleIncomingFragment(mProduct, atom, isHandledFragmentAtom));

		mIsStoichiometrical = mReactantAtomsForRemoval.length == 0 && mRuleFragmentList.isEmpty();

		mInvertedTHParity = new int[0];
		for (int productAtom=0; productAtom<mProduct.getAtoms(); productAtom++) {
			int productParity = mProduct.getAtomParity(productAtom);
			if (productParity == Molecule.cAtomParity1
			 || productParity == Molecule.cAtomParity2) {
				int reactantAtom = mMapNoToReactantAtom[mProduct.getAtomMapNo(productAtom)];
				if (reactantAtom != -1) {
					int reactantParity = mReactant.getAtomParity(reactantAtom);
					if (isTHParityInversion(productAtom, mMapNoToReactantAtom) == (reactantParity == productParity))
						addInvertedParityAtom(reactantAtom);
					}
				}
			}

		mPanalty = calculatePenalty();
		}

	private boolean isStoichiometrical() {
		return mIsStoichiometrical;
	}

	private float calculatePenalty() {
		MappingScorer scorer = new MappingScorer(mReactant, mProduct);
		int[] reactantMapNo = new int[mReactant.getAllAtoms()];
		for (int i=0; i<mReactant.getAllAtoms(); i++)
			reactantMapNo[i] = mReactant.getAtomMapNo(i);
		int[] productMapNo = new int[mProduct.getAllAtoms()];
		for (int i=0; i<mProduct.getAllAtoms(); i++)
			productMapNo[i] = mProduct.getAtomMapNo(i);

		return -scorer.scoreMapping(scorer.createReactantToProductAtomMap(reactantMapNo, productMapNo)) * 0.25f;
//		return -scorer.scoreMapping(scorer.createReactantToProductAtomMap(reactantMapNo, productMapNo)) - 1.5f;	// is a positive value

		// The idea is that when a rule is applied, then the score should be better than
		// the simple calculated score from bond changes, because we know that we use
		// reasonable chemistry. How much better, whether a constant of whether rule-dependent
		// remains to be determined...
	}

	private boolean isTHParityInversion(int reactantAtom, int[] mapNoToProduct) {
		boolean inversion = false;
		if (mReactant.getAtomPi(reactantAtom) == 0) {
			for (int i=1; i<mReactant.getConnAtoms(reactantAtom); i++) {
				for (int j=0; j<i; j++) {
					int connAtom1 = mReactant.getConnAtom(reactantAtom, i);
					int connAtom2 = mReactant.getConnAtom(reactantAtom, j);
					int connMapNo1 = mReactant.getAtomMapNo(connAtom1);
					int connMapNo2 = mReactant.getAtomMapNo(connAtom2);
					if ((mapNoToProduct[connMapNo1] > mapNoToProduct[connMapNo2]) ^ (connAtom1 > connAtom2))
						inversion = !inversion;
					}
				}
			}
		return inversion;
		}

	private void addInvertedParityAtom(int atom) {
		mInvertedTHParity = Arrays.copyOf(mInvertedTHParity, mInvertedTHParity.length+1);
		mInvertedTHParity[mInvertedTHParity.length-1] = atom;
		}

	/**
	 * Applies the chemical rule to the reactant to basically convert the given reactant
	 * into the rule's product, which will be the basis then for similarity graph based mapping
	 * in the hope that the mapping is a much simpler one then.
	 * Applying the rule means changing connectivity and possibly removing some atoms that are
	 * part of the rule's reactant but not of the rule's product.
	 * @param reactant
	 * @param match
	 * return
	 */
	public int[] apply(StereoMolecule reactant, int[] match) {
		reactant.ensureHelperArrays(Molecule.cHelperNeighbours);

		// If we have unmapped atoms in the rule's reactant, then mark these atoms for removal.
		// Also mark all connected atoms for removal that either are not represented in the match
		// or that are not mapped in the match.
/*		if (removeLeavingAndAddIncomingAtoms) {
			if (mReactantAtomsForRemoval.length != 0) {
				boolean[] isMappedRuleAtom = new boolean[reactant.getAtoms()];
				for (int i=0; i<match.length; i++)
					if (match[i] != -1 && mReactant.getAtomMapNo(i) != 0)
						isMappedRuleAtom[match[i]] = true;

				for (int atom : mReactantAtomsForRemoval)
					if (!reactant.isAtomMarkedForDeletion(match[atom]))
						markSubstituentForDeletion(reactant, match[atom], isMappedRuleAtom);
				}
			}*/

		for (ChemicalRuleBond ruleBond:mRuleBonds) {
			int reactantAtom1 = match[ruleBond.atom1];
			int reactantAtom2 = match[ruleBond.atom2];
			int reactantBond = reactant.getBond(reactantAtom1, reactantAtom2);
			if (reactantBond == -1)
				reactant.addBond(reactantAtom1, reactantAtom2, ruleBond.newBondType);
			else if (ruleBond.newBondType == ChemicalRuleBond.BOND_TYPE_DELETE)
				reactant.markBondForDeletion(reactantBond);
			else if (ruleBond.newBondType != ChemicalRuleBond.BOND_TYPE_KEEP_UNCHANGED)
				reactant.setBondType(reactantBond, ruleBond.newBondType);
			}

		int[] originalToAppliedAtom = reactant.deleteMarkedAtomsAndBonds();

//		if (removeLeavingAndAddIncomingAtoms)
//			for (RuleIncomingFragment ruleFragment : mRuleFragmentList)
//				ruleFragment.apply(reactant, match, originalToAppliedAtom, mMapNoToReactantAtom);

		if (mInvertedTHParity.length != 0) {
			reactant.ensureHelperArrays(Molecule.cHelperRings);
			for (int atom:mInvertedTHParity) {
				int reactantAtom = (originalToAppliedAtom == null) ? match[atom] : originalToAppliedAtom[match[atom]];
				int reactantParity = reactant.getAtomParity(reactantAtom);
				reactant.setAtomParity(reactantAtom, reactantParity == Molecule.cAtomParity1 ?
						Molecule.cAtomParity2 : Molecule.cAtomParity1, false);
				reactant.setStereoBondFromAtomParity(reactantAtom);
			}
		}

		return originalToAppliedAtom;
	}

	public boolean[][] createVetoMatrix(int reactantAtoms, int[] reactantMatch, int productAtoms, int[] productMatch) {
		if (mReactantAtomsForRemoval.length == 0 || mRuleFragmentList.isEmpty())
			return null;

		boolean[][] vetoMatrix = new boolean[reactantAtoms][];
		boolean[] productMask = new boolean[productAtoms];
		for (int i=0; i<productMatch.length; i++)
			if (productMatch[i] != -1 && mProduct.getAtomMapNo(i) == 0)
				productMask[productMatch[i]] = true;
		for (int i=0; i<reactantMatch.length; i++)
			if (reactantMatch[i] != -1 && mReactant.getAtomMapNo(i) == 0)
				vetoMatrix[reactantMatch[i]] = productMask;

		return vetoMatrix;
		}

	/**
	 * Starting from atom marks all connected atoms for removal, which either are not part of the match
	 * or where the associated atom in the matching rule reactant has no mapping number.
	 * @param reactant of reaction to be mapped
	 * @param atom first atom of potentially more to be removed from reactant
	 * @param isMappedRuleAtom
	 */
	private void markSubstituentForDeletion(StereoMolecule reactant, int atom, boolean[] isMappedRuleAtom) {
		int[] workAtom = new int[reactant.getAtoms()];
		workAtom[0] = atom;
		int highest = 0;
		for (int i=0; i<=highest; i++) {
			int currentAtom = workAtom[highest];
			reactant.markAtomForDeletion(currentAtom);
			for (int j=0; j<reactant.getConnAtoms(currentAtom); j++) {
				int connAtom = reactant.getConnAtom(currentAtom, j);
				if (!isMappedRuleAtom[connAtom] && !reactant.isAtomMarkedForDeletion(connAtom)) {
					reactant.markAtomForDeletion(connAtom);
					workAtom[++highest] = connAtom;
				}
			}
		}
	}

	public StereoMolecule getReactant() {
		return mReactant;
	}

	public StereoMolecule getProduct() {
		return mProduct;
	}

	public String getName() {
		// extention S:stoichiometrical; L:leavingAtoms; I:incomingAtoms; U:unbalenced (leaving and incoming atoms)
		return mName + (mIsStoichiometrical ? "_S" : mReactantAtomsForRemoval.length == 0 ? "_L" : mRuleFragmentList.isEmpty() ? "_I" : "_U");
	}

	public float getPanalty() {
		return mPanalty;
	}

	/**
	 * If the rule's reactant matches the real reactant multiple times,
	 * then some of these matches may be symmetrically equivalent. To avoid building and
	 * scoring redundant mapping graphs, these should be sorted out early. Reasons for
	 * redundant matches may be:<br>
	 * - if the rule reactant is one fragment, this may be symmetrical Cn or Dn<br>
	 * - the rule reactant may contain multiple equivalent fragments, e.g. metathesis<br>
	 * - matching atoms in the real reactant may be symmetrical<br>
	 * Otherwise, there certain causes may exist, that break these symmetries:<br>
	 * - a symmetrical rule fragment must not be considered symmetrical, if it doesn't react
	 *   symmetrically, i.e. if its matching rule product atoms are not equivalent anymore.
	 * - in case of multiple symmetrical fragments in the rule's reactant (e.g. metathesis),
	 *   inverted/rotated individual fragment matches cause different products, that means
	 *   that the relative match orientation of multiple symmetrical fragments breaks symmetry,
	 *   if the real matching atoms are not equivalent.<br>
	 * This method calculates symmetry breaking values using these two reasons to be passed to
	 * the reactant rule substructure searcher.
	 */
	private void calculateReactantAtomSymmetryConstraints(int[] mapNoToReactantAtom) {
		// break symmetries because of un-symmetrical rule products
		mReactantAtomSymmetryConstraint = new int[mReactant.getAtoms()];
		for (int atom=0; atom<mProduct.getAtoms(); atom++)
			if (mProduct.getAtomMapNo(atom) != 0)
				mReactantAtomSymmetryConstraint[mapNoToReactantAtom[mProduct.getAtomMapNo(atom)]] = mProduct.getSymmetryRank(atom);

		int[] fragmentNo = new int[mReactant.getAllAtoms()];
		int fragmentCount = mReactant.getFragmentNumbers(fragmentNo, false, false);
		if (fragmentCount > 1) {
			int[] atomIndex = new int[fragmentCount];
			for (int atom=0; atom<mReactant.getAtoms(); atom++)
				mReactantAtomSymmetryConstraint[atom] |= (atomIndex[fragmentNo[atom]]++) << 5;
			}
		}

	public int[] getReactantAtomSymmetryConstraints() {
		return mReactantAtomSymmetryConstraint;
		}
/*
	public ChemicalRuleBond[] getBondsToModify() {
		return mRuleBonds;
	}*/
}

class ChemicalRuleBond implements Comparable<ChemicalRuleBond> {
	static final int BOND_TYPE_KEEP_UNCHANGED = -2;
	static final int BOND_TYPE_DELETE = -1;

	int atom1,atom2,mapNo1,mapNo2,newBondType;

	public ChemicalRuleBond(int atom1, int atom2, int mapNo1, int mapNo2, int newBondType) {
		if (mapNo1 < mapNo2) {
			this.atom1 = atom1;
			this.atom2 = atom2;
			this.mapNo1 = mapNo1;
			this.mapNo2 = mapNo2;
		}
		else {
			this.atom1 = atom2;
			this.atom2 = atom1;
			this.mapNo1 = mapNo2;
			this.mapNo2 = mapNo1;
		}
		this.newBondType = newBondType;
	}

	@Override
	public int compareTo(ChemicalRuleBond crb) {
		if (mapNo1< crb.mapNo1)
			return -1;
		if (mapNo1> crb.mapNo1)
			return 1;
		if (mapNo2< crb.mapNo2)
			return -1;
		if (mapNo2> crb.mapNo2)
			return 1;
		return 0;
	}
}

class RuleIncomingFragment {
	private final StereoMolecule mRuleProduct;
	private final int[] mGraphAtom;
	private final int mGraphAtoms;

	public RuleIncomingFragment(StereoMolecule ruleProduct, int rootAtom, boolean[] isUsedAtom) {
		mGraphAtom = new int[ruleProduct.getAtoms()];
		mGraphAtom[0] = rootAtom;
		isUsedAtom[rootAtom] = true;
		int current = 0;
		int highest = 0;
		while (current <= highest) {
			int connAtoms = ruleProduct.getConnAtoms(mGraphAtom[current]);
			for (int i=0; i<connAtoms; i++) {
				int candidate = ruleProduct.getConnAtom(mGraphAtom[current], i);
				if (ruleProduct.getAtomMapNo(candidate) == 0
				 && !ruleProduct.isExcludeGroupAtom(candidate)
				 && !isUsedAtom[candidate]) {
					mGraphAtom[++highest] = candidate;
					isUsedAtom[candidate] = true;
				}
			}
			current++;
		}

		mGraphAtoms = current;
		mRuleProduct = ruleProduct;
	}

	public void apply(StereoMolecule reactant, int[] match, int[] originalToAppliedAtom, int[] mapNoToReactantAtom) {
		// copy fragment atoms first
		int[] ruleToReactantAtom = new int[mRuleProduct.getAtoms()];
		for (int i=0; i<mGraphAtoms; i++)
			ruleToReactantAtom[mGraphAtom[i]] = mRuleProduct.copyAtom(reactant, mGraphAtom[i], 0, 0);

		// then copy bonds between fragment atoms and between a fragment atom and a mapped atom
		boolean[] isBondHandled = new boolean[mRuleProduct.getBonds()];
		for (int i=0; i<mGraphAtoms; i++) {
			for (int j=0; j<mRuleProduct.getConnAtoms(mGraphAtom[i]); j++) {
				int connAtom = mRuleProduct.getConnAtom(mGraphAtom[i], j);
				if (!mRuleProduct.isExcludeGroupAtom(connAtom)) {
					int connBond = mRuleProduct.getConnBond(mGraphAtom[i], j);
					if (mRuleProduct.getAtomMapNo(connAtom) != 0) {
						int matchingConnAtom = match[mapNoToReactantAtom[mRuleProduct.getAtomMapNo(connAtom)]];
						mRuleProduct.copyBond(reactant, connBond, 0, 0,
								ruleToReactantAtom[mGraphAtom[i]], originalToAppliedAtom == null ? matchingConnAtom : originalToAppliedAtom[matchingConnAtom], false);
					}
					else if (!isBondHandled[connBond]) {
						mRuleProduct.copyBond(reactant, connBond, 0, 0,
								ruleToReactantAtom[mGraphAtom[i]], ruleToReactantAtom[connAtom], false);
						isBondHandled[connBond] = true;
					}
				}
			}
		}
	}
}
