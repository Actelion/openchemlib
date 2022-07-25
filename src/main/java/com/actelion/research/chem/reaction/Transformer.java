package com.actelion.research.chem.reaction;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.SortedList;

import java.util.ArrayList;

/**
 * Like a Reactor, a Transformer uses a transformation rule (passed as generic Reaction) to apply defined
 * bond changes to one or more real world molecules. It is much faster than a Reactor, but requires the given
 * transformation to consist of one reactant and one product only. It also requires all atoms to be mapped.
 * During Transformer construction the given Reaction is translated into rule list of bond changes.
 * When a molecule is given with setMolecule(), then substructure search is performed to determine all
 * possible matches of the transformation reactant. Matches, which would exceed an atom valence when applying the
 * transformation, are sorted out. applyTransformation() can be called then up to the number of valid matches to
 * construct transformed molecules.
 */
public class Transformer {
	private SortedList<TransformerRule> mRuleList;
	private int[] mMinFreeValence;	// minimum required free valence on reactant atoms
	private StereoMolecule mTargetMolecule;
	private Reaction mReaction;
	private SSSearcher mSSSearcher;
	private ArrayList<int[]> mMatchList;

	public Transformer(Reaction reaction) {
		assert reaction.getReactants() == 1 && reaction.getProducts() == 1 && reaction.isPerfectlyMapped();
		mReaction = reaction;
		StereoMolecule reactant = reaction.getReactant(0);
		StereoMolecule product = reaction.getProduct(0);
		reactant.ensureHelperArrays(Molecule.cHelperRings);
		product.ensureHelperArrays(Molecule.cHelperRings);
		mRuleList = new SortedList<>();
		for (int bond=0; bond<reactant.getBonds(); bond++)
			mRuleList.add(new TransformerRule(reactant, bond));
		for (int bond=0; bond<product.getBonds(); bond++) {
			TransformerRule rule = new TransformerRule(product, bond);
			int index = mRuleList.getIndex(rule);
			if (index == -1) {
				rule.finishProductOnly();
				mRuleList.add(rule);
				}
			else {
				mRuleList.get(index).finishWithProduct(rule);
				}
			}
		for (int i=mRuleList.size()-1; i>=0; i--) {
			TransformerRule rule = mRuleList.get(i);
			if (!rule.isFinished())
				rule.finishNoProduct();
			else if (mRuleList.get(i).isNoChange())
				mRuleList.remove(i);
			}

		// calculate minimum free valence of reactant atoms
		mMinFreeValence = new int[reactant.getAtoms()];
		for (int i=0; i<reactant.getAtoms(); i++) {
			for (int j=0; j<product.getAtoms(); j++) {
				if (product.getAtomMapNo(j) == reactant.getAtomMapNo(i)) {
					int dif = reactant.getFreeValence(i) - product.getFreeValence(j);
					mMinFreeValence[i] = (dif > 0) ? dif : 0;
					break;
					}
				}
			}

		mSSSearcher = new SSSearcher();
		}

	/**
	 * Use this method before calling applyTransformation() to actually perform the bond changes.
	 * This runs a substructure search of the transformation's reactant on the passed molecule
	 * to determine all possible matches and, thus, the number of possible transformations that
	 * can be applied on the molecule without causing a valence problem.
	 * @param molecule
	 * @param countMode one of SSSearch.cCountMode...
	 * @return number of valid transformation to be applied to this molecule
	 */
	public int setMolecule(StereoMolecule molecule, int countMode) {
		mTargetMolecule = molecule;
		StereoMolecule reactant = mReaction.getReactant(0);

		mSSSearcher.setMol(reactant, mTargetMolecule);
		int matchMode = SSSearcher.cDefaultMatchMode;
		if (mSSSearcher.findFragmentInMolecule(countMode, matchMode) == 0) {
			mMatchList = null;
			return 0;
			}

		// eliminate matches where reaction would exceed an atom valence
		mMatchList = mSSSearcher.getMatchList();
		for (int j=mMatchList.size()-1; j>=0; j--) {
			int[] matchingAtom = mMatchList.get(j);
			for (int k=0; k<matchingAtom.length; k++) {
				if (matchingAtom[k] != -1) {
					if (mMinFreeValence[k] > 0
					 && mMinFreeValence[k] > mTargetMolecule.getFreeValence(matchingAtom[k])) {
						mMatchList.remove(j);
						break;
					}
				}
			}
		}

		return mMatchList.size();
		}

	/**
	 * Applies the transformation to the given molecule using the chosen substructure match.
	 * @param mol
	 * @param no must be smaller than the number of valid matches returned by setMolecule()
	 */
	public void applyTransformation(StereoMolecule mol, int no) {
		int[] matchAtom = mMatchList.get(no);
		for (int i=0; i<mRuleList.size(); i++) {
			TransformerRule rule = mRuleList.get(i);
			rule.adaptBondOrder(mol, matchAtom);
			}
		mol.deleteMarkedAtomsAndBonds();
		}
	}

class TransformerRule extends Object implements Comparable<TransformerRule> {
	private enum TYPE { CREATE, CHANGE_DIF, CHANGE_ABS, REMOVE, NO_CHANGE, UNKNOWN };

	private int mAtom1,mAtom2,mMapNo1,mMapNo2,mBondQFTypes,mTargetBondType,mTargetBondDif;
	private TYPE mType;

	public TransformerRule(StereoMolecule mol, int bond) {
		mAtom1 = mol.getBondAtom(0, bond);
		mAtom2 = mol.getBondAtom(1, bond);
		int mapNo1 = mol.getAtomMapNo(mAtom1);
		int mapNo2 = mol.getAtomMapNo(mAtom2);
		if (mapNo1 < mapNo2) {
			mMapNo1 = mapNo1;
			mMapNo2 = mapNo2;
			}
		else {
			mMapNo1 = mapNo2;
			mMapNo2 = mapNo1;
			}

		mBondQFTypes = mol.getBondQueryFeatures(bond) & Molecule.cBondQFBondTypes;

//		if (mol.getBondType(bond) == Molecule.cBondTypeDelocalized || mol.isDelocalizedBond(bond))
//			mBondQFTypes |= Molecule.cBondQFDelocalized;
//		else {
//			int order = mol.getBondOrder(bond);
//			if (order == 0)
//				mBondQFTypes |= Molecule.cBondQFMetalLigand;
//			else if (order == 1)
//				mBondQFTypes |= Molecule.cBondQFSingle;
//			else if (order == 2)
//				mBondQFTypes |= Molecule.cBondQFDouble;
//			else if (order == 3)
//				mBondQFTypes |= Molecule.cBondQFTriple;
//			}

		// we don't use the delocalized type
		int bondType = mol.getBondTypeSimple(bond);
		if (bondType == Molecule.cBondTypeMetalLigand)
			mBondQFTypes |= Molecule.cBondQFMetalLigand;
		else if (bondType == Molecule.cBondTypeDouble)
			mBondQFTypes |= Molecule.cBondQFDouble;
		else if (bondType == Molecule.cBondTypeTriple)
			mBondQFTypes |= Molecule.cBondQFTriple;
		else
			mBondQFTypes |= Molecule.cBondQFSingle;

		mType = TYPE.UNKNOWN;
		}

	@Override
	public boolean equals(Object o) {
		return mMapNo1 == ((TransformerRule)o).mMapNo1 && mMapNo2 == ((TransformerRule)o).mMapNo2;
		}

	@Override
	public int compareTo(TransformerRule o) {
		if (mMapNo1 != o.mMapNo1)
			return mMapNo1 < o.mMapNo1 ? -1 : 1;
		if (mMapNo2 != o.mMapNo2)
			return mMapNo2 < o.mMapNo2 ? -1 : 1;
		return 0;
		}

	public boolean isNoChange() {
		return mType == TYPE.NO_CHANGE;
		}

	public boolean isFinished() {
		return mType != TYPE.UNKNOWN;
		}

	public void finishWithProduct(TransformerRule productRule) {
		if (mBondQFTypes == productRule.mBondQFTypes) {
			mType = TYPE.NO_CHANGE;
			}
		else if (Integer.bitCount(productRule.mBondQFTypes) == 1) {
			mTargetBondType = productRule.mBondQFTypes == Molecule.cBondQFSingle ? Molecule.cBondTypeSingle
							: productRule.mBondQFTypes == Molecule.cBondQFDouble ? Molecule.cBondTypeDouble
							: productRule.mBondQFTypes == Molecule.cBondQFTriple ? Molecule.cBondTypeTriple
							: productRule.mBondQFTypes == Molecule.cBondQFMetalLigand ? Molecule.cBondTypeMetalLigand
							: Molecule.cBondTypeDelocalized;
			mType = TYPE.CHANGE_ABS;
			}
		else {
			int sourceBondOrder = (mBondQFTypes & Molecule.cBondQFMetalLigand) != 0 ? 0
								: (mBondQFTypes & Molecule.cBondQFSingle) != 0 ? 1
								: (mBondQFTypes & Molecule.cBondQFDouble) != 0 ? 2 : 3;
			int targetBondOrder = (productRule.mBondQFTypes & Molecule.cBondQFMetalLigand) != 0 ? 0
								: (productRule.mBondQFTypes & Molecule.cBondQFSingle) != 0 ? 1
								: (productRule.mBondQFTypes & Molecule.cBondQFDouble) != 0 ? 2 : 3;
			if (targetBondOrder == sourceBondOrder)
				mType = TYPE.NO_CHANGE;
			else {
				mTargetBondDif = targetBondOrder - sourceBondOrder;
				mType = TYPE.CHANGE_DIF;
				}
			}
		}

	public void finishProductOnly() {
		mType = TYPE.CREATE;
		}

	public void finishNoProduct() {
		mType = TYPE.REMOVE;
		}

	public boolean adaptBondOrder(StereoMolecule mol, int[] matchAtom) {
		int atom1 = matchAtom[mAtom1];
		int atom2 = matchAtom[mAtom2];

		if (mType == TYPE.CREATE) {
//			int bondOrder = mTargetBondType == Molecule.cBondTypeMetalLigand ? 0
//						  : mTargetBondType == Molecule.cBondTypeDouble ? 2
//						  : mTargetBondType == Molecule.cBondTypeTriple ? 3 : 1;

//			if (mol.getFreeValence(atom1) < bondOrder
//			 || mol.getFreeValence(atom2) < bondOrder)
//				return false;

			mol.addBond(atom1, atom2, mTargetBondType);
			return true;
			}
		else {
			int bond = mol.getBond(atom1, atom2);
			if (mType == TYPE.CHANGE_ABS) {
//				int bondOrderChange = (mTargetBondType == Molecule.cBondTypeMetalLigand ? 0
//									 : mTargetBondType == Molecule.cBondTypeDouble ? 2
//									 : mTargetBondType == Molecule.cBondTypeTriple ? 3 : 1) - mol.getBondOrder(bond);

//				if (mol.getFreeValence(atom1) < bondOrderChange
//				 || mol.getFreeValence(atom2) < bondOrderChange)
//					return false;

				mol.setBondType(bond, mTargetBondType);
				return true;
				}
			else if (mType == TYPE.CHANGE_DIF) {
//				if (mTargetBondDif > 0
//				 && (mol.getFreeValence(atom1) < mTargetBondDif
//				  || mol.getFreeValence(atom2) < mTargetBondDif))
//					return false;

				int targetBondOrder = mol.getBondOrder(bond) + mTargetBondDif;
				if (targetBondOrder >= 0 && targetBondOrder <= 3) {
					int targetBondType = targetBondOrder == 0 ? Molecule.cBondTypeMetalLigand
									   : targetBondOrder == 1 ? Molecule.cBondTypeSingle
									   : targetBondOrder == 2 ? Molecule.cBondTypeDouble : Molecule.cBondTypeTriple;
					mol.setBondType(bond, targetBondType);
					return true;
					}
				}
			else if (mType == TYPE.REMOVE) {
				mol.markBondForDeletion(bond);
				return true;
				}
			}
		return false;
		}
	}