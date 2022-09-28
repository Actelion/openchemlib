package com.actelion.research.chem.reaction;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.util.SortedList;

import java.util.ArrayList;
import java.util.Arrays;

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
	private int mBondGain;
	private StereoMolecule mReactant,mTargetMolecule;
	private SSSearcher mSSSearcher;
	private ArrayList<int[]> mMatchList;
	private String mName;

	public Transformer(StereoMolecule reactant, StereoMolecule product, String name) {
		mName = name;
		reactant.ensureHelperArrays(Molecule.cHelperRings);
		product.ensureHelperArrays(Molecule.cHelperRings);
		mBondGain = product.getBonds() - reactant.getBonds();
		mReactant = reactant;
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

		mSSSearcher.setMol(mReactant, mTargetMolecule);
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
	 * Creates a canonical representation of the bond state after the transformation
	 * without actually touching the molecule itself.
	 * Together with a history of created bonds states, this can be used to avoid
	 * transformations that lead to products that have been seen before.
	 * @return
	 */
	public int[] getTransformedBondList(StereoMolecule mol, int matchNo) {
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);

		int[] buffer = new int[mol.getBonds()+mBondGain];
		int index = 0;

		int[] matchAtom = mMatchList.get(matchNo);
		boolean[] bondHandled = new boolean[mol.getBonds()];
		for (int i=0; i<mRuleList.size(); i++) {
			int code = mRuleList.get(i).addBondCode(mol, matchAtom, bondHandled);
			if (code != -1)
				buffer[index++] = code;
			}

		for (int bond=0; bond<mol.getBonds(); bond++)
			if (!bondHandled[bond])
				buffer[index++] = TransformerRule.getBondCode(mol.getBondAtom(0, bond),
															  mol.getBondAtom(1, bond),
															  mol.getBondOrder(bond));

		Arrays.sort(buffer);
		return buffer;
		}

	/**
	 * Applies the transformation to the given molecule using the chosen substructure match.
	 * @param mol
	 * @param matchNo must be smaller than the number of valid matches returned by setMolecule()
	 */
	public void applyTransformation(StereoMolecule mol, int matchNo) {
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		int[] matchAtom = mMatchList.get(matchNo);
		for (int i=0; i<mRuleList.size(); i++) {
			TransformerRule rule = mRuleList.get(i);
			rule.adaptBondOrder(mol, matchAtom);
			}
		mol.deleteMarkedAtomsAndBonds();
		}

	public String getName() {
		return mName;
		}
	}
