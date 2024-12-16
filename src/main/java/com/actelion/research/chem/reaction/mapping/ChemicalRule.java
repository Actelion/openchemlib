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
	private String mName,mIDCode;
	private float mPanalty;
	private StereoMolecule mReactant,mProduct;
	private ChemicalRuleBond[] mRuleBonds;
	private int[] mInvertedTHParity;
	private int[] mReactantAtomSymmetryConstraint;

	public ChemicalRule(String name, String idcode, float panalty) {
		mName = name;
		mIDCode = idcode;
		mPanalty = panalty;
	}

	public void initialize() {
		Reaction rxn = ReactionEncoder.decode(mIDCode, false);
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

		int[] mapNoToReactantAtom = new int[mReactant.getAtoms()+1];
		mapNoToReactantAtom[0] = -1;    // to cause exceptions in case of faulty logic
		for (int atom=0; atom<mReactant.getAtoms(); atom++)
			mapNoToReactantAtom[mReactant.getAtomMapNo(atom)] = atom;

		calculateReactantAtomSymmetryConstraints(mapNoToReactantAtom);

		boolean[] reactantBondFoundInProduct = new boolean[mReactant.getBonds()];
		for (int productBond=0; productBond<mProduct.getBonds(); productBond++) {
			int mapNo1 = mProduct.getAtomMapNo(mProduct.getBondAtom(0, productBond));
			int mapNo2 = mProduct.getAtomMapNo(mProduct.getBondAtom(1, productBond));
			if (mapNo1 != 0 && mapNo2 != 0) {   // exclude atoms are not mapped and don't go into the bond list
				int atom1 = mapNoToReactantAtom[mapNo1];
				int atom2 = mapNoToReactantAtom[mapNo2];
				int productBondType = mProduct.getBondTypeSimple(productBond);
				int reactantBond = mReactant.getBond(atom1, atom2);
				if (reactantBond == -1) {
					bondList.add(new ChemicalRuleBond(atom1, atom2, mapNo1, mapNo2, productBondType));
					}
				else {
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
				int mapNo1 = mReactant.getAtomMapNo(atom1);
				int mapNo2 = mReactant.getAtomMapNo(atom2);
				if (mapNo1 != 0 && mapNo2 != 0)   // exclude atoms are not mapped and don't go into the bond list
					bondList.add(new ChemicalRuleBond(atom1, atom2, mapNo1, mapNo2, ChemicalRuleBond.BOND_TYPE_DELETE));
				}
			}

		mRuleBonds = bondList.toArray(new ChemicalRuleBond[0]);

		mInvertedTHParity = new int[0];
		for (int productAtom=0; productAtom<mProduct.getAtoms(); productAtom++) {
			int productParity = mProduct.getAtomParity(productAtom);
			if (productParity == Molecule.cAtomParity1
			 || productParity == Molecule.cAtomParity2) {
				int reactantAtom = mapNoToReactantAtom[mProduct.getAtomMapNo(productAtom)];
				if (reactantAtom != -1) {
					int reactantParity = mReactant.getAtomParity(reactantAtom);
					if (isTHParityInversion(productAtom, mapNoToReactantAtom) == (reactantParity == productParity))
						addInvertedParityAtom(reactantAtom);
					}
				}
			}

		calculatePenalty();
		}

	private void calculatePenalty() {
		MappingScorer scorer = new MappingScorer(mReactant, mProduct);
		int[] reactantMapNo = new int[mReactant.getAllAtoms()];
		for (int i=0; i<mReactant.getAllAtoms(); i++)
			reactantMapNo[i] = mReactant.getAtomMapNo(i);
		int[] productMapNo = new int[mProduct.getAllAtoms()];
		for (int i=0; i<mProduct.getAllAtoms(); i++)
			productMapNo[i] = mProduct.getAtomMapNo(i);

		mPanalty = -scorer.scoreMapping(scorer.createReactantToProductAtomMap(reactantMapNo, productMapNo)) * 0.25f;
//		mPanalty = -scorer.scoreMapping(scorer.createReactantToProductAtomMap(reactantMapNo, productMapNo)) - 1.5f;	// is a positive value

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

	public void apply(StereoMolecule reactant, int[] match) {
		reactant.ensureHelperArrays(Molecule.cHelperNeighbours);
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
		reactant.deleteMarkedAtomsAndBonds();

		if (mInvertedTHParity.length != 0) {
			reactant.ensureHelperArrays(Molecule.cHelperRings);
			for (int atom:mInvertedTHParity) {
				int reactantAtom = match[atom];
				int reactantParity = reactant.getAtomParity(reactantAtom);
				reactant.setAtomParity(reactantAtom, reactantParity == Molecule.cAtomParity1 ?
						Molecule.cAtomParity2 : Molecule.cAtomParity1, false);
				reactant.setStereoBondFromAtomParity(reactantAtom);
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
		return mName;
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

	public void setNewBondType(int newBondType) {
		this.newBondType = newBondType;
	}
}
