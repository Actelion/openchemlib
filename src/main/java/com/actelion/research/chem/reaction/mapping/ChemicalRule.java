package com.actelion.research.chem.reaction.mapping;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.chem.reaction.ReactionEncoder;
import com.actelion.research.util.SortedList;

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

	public ChemicalRule(String name, String idcode, int panalty) {
		mName = name;
		mIDCode = idcode;
		mPanalty = panalty;
	}

	public void initialize() {
		Reaction rxn = ReactionEncoder.decode(mIDCode, false);
		mReactant = rxn.getReactant(0);
		mProduct = rxn.getProduct(0);
		mReactant.ensureHelperArrays(Molecule.cHelperNeighbours);
		mProduct.ensureHelperArrays(Molecule.cHelperNeighbours);

		// key: lower bondAtom mapNo, higher bondAtom mapNo
		// value: reactantAtom1, reactantAtom2, productBondOrder

		SortedList<ChemicalRuleBond> bondList = new SortedList<>();

		int[] mapNoToReactantAtom = new int[mReactant.getAtoms()+1];
		for (int atom=0; atom<mReactant.getAtoms(); atom++)
			mapNoToReactantAtom[mReactant.getAtomMapNo(atom)] = atom;

		for (int bond=0; bond<mProduct.getBonds(); bond++) {
			int mapNo1 = mProduct.getAtomMapNo(mProduct.getBondAtom(0, bond));
			int mapNo2 = mProduct.getAtomMapNo(mProduct.getBondAtom(1, bond));
			int atom1 = mapNoToReactantAtom[mapNo1];
			int atom2 = mapNoToReactantAtom[mapNo2];
			int productBondType = ((mProduct.getBondQueryFeatures(bond) & Molecule.cBondQFBondTypes) != 0) ?
					ChemicalRuleBond.BOND_TYPE_KEEP_UNCHANGED : mProduct.getBondTypeSimple(bond);
			bondList.add(new ChemicalRuleBond(atom1, atom2, mapNo1, mapNo2, productBondType));
		}

		for (int bond=0; bond<mReactant.getBonds(); bond++) {
			int atom1 = mReactant.getBondAtom(0, bond);
			int atom2 = mReactant.getBondAtom(1, bond);
			int mapNo1 = mReactant.getAtomMapNo(atom1);
			int mapNo2 = mReactant.getAtomMapNo(atom2);
			ChemicalRuleBond ruleBond = new ChemicalRuleBond(atom1, atom2, mapNo1, mapNo2, ChemicalRuleBond.BOND_TYPE_DELETE);
			if (!bondList.contains(ruleBond))   // if we don't have a product bond for these two atoms, we need to delete...
				bondList.add(ruleBond);
		}

		mRuleBonds = bondList.toArray(new ChemicalRuleBond[0]);
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
