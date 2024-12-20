package com.actelion.research.chem.reaction;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

public class TransformerRule extends Object implements Comparable<TransformerRule> {
	public enum TYPE { CREATE, CHANGE_DIF, CHANGE_ABS, REMOVE, NO_CHANGE, UNKNOWN };

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
			mBondQFTypes |= Molecule.cBondTypeMetalLigand;
		else if (bondType == Molecule.cBondTypeDouble)
			mBondQFTypes |= Molecule.cBondTypeDouble;
		else if (bondType == Molecule.cBondTypeTriple)
			mBondQFTypes |= Molecule.cBondTypeTriple;
		else
			mBondQFTypes |= Molecule.cBondTypeSingle;

		mType = TYPE.UNKNOWN;
	}

	public int getAtom1() {
		return mAtom1;
	}

	public int getAtom2() {
		return mAtom2;
	}

	public TYPE getType() {
		return mType;
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
			mTargetBondType = productRule.mBondQFTypes == Molecule.cBondTypeSingle ? Molecule.cBondTypeSingle
					: productRule.mBondQFTypes == Molecule.cBondTypeDouble ? Molecule.cBondTypeCross
					: productRule.mBondQFTypes == Molecule.cBondTypeTriple ? Molecule.cBondTypeTriple
					: productRule.mBondQFTypes == Molecule.cBondTypeMetalLigand ? Molecule.cBondTypeMetalLigand
					: Molecule.cBondTypeDelocalized;
			mType = TYPE.CHANGE_ABS;
		}
		else {
			int sourceBondOrder = (mBondQFTypes & Molecule.cBondTypeMetalLigand) != 0 ? 0
					: (mBondQFTypes & Molecule.cBondTypeSingle) != 0 ? 1
					: (mBondQFTypes & Molecule.cBondTypeDouble) != 0 ? 2 : 3;
			int targetBondOrder = (productRule.mBondQFTypes & Molecule.cBondTypeMetalLigand) != 0 ? 0
					: (productRule.mBondQFTypes & Molecule.cBondTypeSingle) != 0 ? 1
					: (productRule.mBondQFTypes & Molecule.cBondTypeDouble) != 0 ? 2 : 3;
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

	/**
	 * Make sure that the molecule's helper status is at least cHelperNeighbours
	 * @param mol
	 * @param matchAtom
	 * @return
	 */
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
					mol.setBondType(bond, Molecule.bondOrderToType(targetBondOrder, true));
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

	public int addBondCode(StereoMolecule mol, int[] matchAtom, boolean[] bondHandled) {
		int atom1 = matchAtom[mAtom1];
		int atom2 = matchAtom[mAtom2];

		if (mType == TYPE.CREATE)
			return getBondCode(atom1, atom2, Molecule.bondTypeToOrder(mTargetBondType));

		int bond = mol.getBond(atom1, atom2);
		if (mType == TYPE.CHANGE_ABS) {
			bondHandled[bond] = true;
			return getBondCode(atom1, atom2, Molecule.bondTypeToOrder(mTargetBondType));
		}

		if (mType == TYPE.CHANGE_DIF) {
			int targetBondOrder = mol.getBondOrder(bond) + mTargetBondDif;
			if (targetBondOrder >= 0 && targetBondOrder <= 3) {
				bondHandled[bond] = true;
				return getBondCode(atom1, atom2, targetBondOrder);
			}
			else {
				return -1;
			}
		}

		if (mType == TYPE.REMOVE) {
			bondHandled[bond] = true;
			return -1;
		}

		return -1;
	}

	public static int getBondCode(int atom1, int atom2, int bondOrder) {
		return (atom1 < atom2) ?
				(atom1 << 17) + (atom2 << 2) + bondOrder
				: (atom2 << 17) + (atom1 << 2) + bondOrder;
	}
}
