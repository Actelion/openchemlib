package com.actelion.research.chem.sar;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;

public class SARMolecule {
	private StereoMolecule mMol,mBuffer;
	private SARScaffold mScaffold;
	private String[] mSubstituent;
	private int[] mSubstituentBondOrder;
	private boolean[] mSubstituentConnectsBack;

	protected SARMolecule(StereoMolecule mol, SARScaffold scaffold, int[] molToCoreAtom, StereoMolecule buffer) {
		mMol = mol;
		mScaffold = scaffold;
		mBuffer = buffer;

		buildSubstituents(molToCoreAtom);
	}

	public SARScaffold getScaffold() {
		return mScaffold;
	}

	public String[] getSubstituents() {
		return mSubstituent;
	}

	protected void clear() {
		mScaffold = null;
	}

	protected void checkSubstituents() {
		for (int exitVectorIndex = 0; exitVectorIndex<mScaffold.getExitVectorCount(); exitVectorIndex++)
			mScaffold.checkSubstituent(mSubstituent[exitVectorIndex], exitVectorIndex, mSubstituentBondOrder[exitVectorIndex]);
	}

	protected void removeUnchangingSubstituents() {
		for (int exitVectorIndex = 0; exitVectorIndex<mScaffold.getExitVectorCount(); exitVectorIndex++)
			if (mSubstituent[exitVectorIndex] != null
			 && !mScaffold.getExitVector(exitVectorIndex).substituentVaries())
				mSubstituent[exitVectorIndex] = null;
	}

	protected void buildSubstituents(int[] molToCoreAtom) {
		StereoMolecule core = mScaffold.getCoreStructure();

		int[] coreToMolAtom = new int[core.getAtoms()];
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (molToCoreAtom[atom] != -1)
				coreToMolAtom[molToCoreAtom[atom]] = atom;

		// We change all molecule atoms, which belong to the core, to connection point atoms
		// in order to easily extract/copy substituents including connections points from the molecule.
		for (int i=0; i<core.getAtoms(); i++)
			mMol.setAtomicNo(coreToMolAtom[i], 0);

		mSubstituentBondOrder = new int[mScaffold.getExitVectorCount()];
		mSubstituentConnectsBack = new boolean[mScaffold.getExitVectorCount()];
		mSubstituent = new String[mScaffold.getExitVectorCount()];

		int[] rootAtom = new int[mScaffold.getExitVectorCount()];
		int[] exitAtom = new int[mScaffold.getExitVectorCount()];
		for (int evi=0; evi<mScaffold.getExitVectorCount(); evi++) {
			rootAtom[evi] = coreToMolAtom[mScaffold.getCoreAtom(evi)];
			exitAtom[evi] = mScaffold.getExitVectorAtom(mMol, coreToMolAtom, molToCoreAtom, evi);
		}

		// For all exit vectors that carry substituents in the molecule,
		// create substituent idcodes and assign them to the respective exit vectors.
		for (int evi=0; evi<mScaffold.getExitVectorCount(); evi++)
			if (exitAtom[evi] != -1)
				encodeSubstituent(evi, rootAtom, exitAtom, molToCoreAtom);
	}

	private void encodeSubstituent(int exitVectorIndex, int[] rootAtoms, int[] exitAtoms, int[] molToCoreAtom) {
		int exitAtom = exitAtoms[exitVectorIndex];
		int rootAtom = rootAtoms[exitVectorIndex];
		int rootBond = mMol.getBond(rootAtom, exitAtom);

		int[] workAtom = new int[mMol.getAtoms()];
		boolean[] isSubstituentAtom = new boolean[mMol.getAtoms()];
		boolean[] isSubstituentBond = new boolean[mMol.getBonds()];
		ArrayList<BackConnection> backConnectionList = new ArrayList<>();

		isSubstituentAtom[rootAtom] = true;
		isSubstituentAtom[exitAtom] = true;
		isSubstituentBond[rootBond] = true;
		workAtom[0] = rootAtom;
		workAtom[1] = exitAtom;
		int current = 1;
		int highest = 1;
		while (current <= highest) {
			for (int i=0; i<mMol.getConnAtoms(workAtom[current]); i++) {
				int candidate = mMol.getConnAtom(workAtom[current], i);
				if (molToCoreAtom[candidate] == -1) {  // atom doesn't belong to core
					if (!isSubstituentAtom[candidate]) {
						isSubstituentAtom[candidate] = true;
						isSubstituentBond[mMol.getConnBond(workAtom[current], i)] = true;
						workAtom[++highest] = candidate;
					}
					else {  // ring closure within substituent
						isSubstituentBond[mMol.getConnBond(workAtom[current], i)] = true;
					}
				}
				else if (current != 1 || candidate != rootAtom) {
					mSubstituentConnectsBack[exitVectorIndex] = true;
					backConnectionList.add(new BackConnection(workAtom[current], candidate, mMol.getConnBondOrder(workAtom[current], i)));
				}
			}
			current++;
		}

		mBuffer.clear();

		int[] molToSubstituentAtom = mMol.copyMoleculeByBonds(mBuffer, isSubstituentBond, false, null);
		for (BackConnection backConnection:backConnectionList) {
			int atom = mBuffer.addAtom(0);
			mBuffer.addBond(molToSubstituentAtom[backConnection.substituentAtom], atom, backConnection.bondType);
			for (int evi=0; evi<exitAtoms.length; evi++) {
				if (rootAtoms[evi] == backConnection.backEndAtom
				 && exitAtoms[evi] == backConnection.substituentAtom) {
					mBuffer.setAtomCustomLabel(atom, Integer.toString(evi));
					break;
				}
			}
		}

		mBuffer.setFragment(false);

		if (!CoreBasedSARAnalyzer.DISTINGUISH_STEREO_CENTERS)
			mBuffer.stripStereoInformation();

		// if substituent is a ring forming bridge to the startatom
		for (int bond=mBuffer.getAllBonds()-1; bond>=0; bond--)
			if (mBuffer.getAtomicNo(mBuffer.getBondAtom(0, bond)) == 0
			 && mBuffer.getAtomicNo(mBuffer.getBondAtom(1, bond)) == 0)
				mBuffer.deleteBond(bond);

		mSubstituent[exitVectorIndex] = new Canonizer(mBuffer, Canonizer.ENCODE_ATOM_CUSTOM_LABELS).getIDCode();
		mSubstituentBondOrder[exitVectorIndex] = mMol.getBondOrder(rootBond);
	}

	/**
	 * In case of substituent atom connecting back to the core structure, the exit vector index is encoded as label
	 * in the substituent idcode.
	 * This is needed within the check for varying substituents, because of the label chains with inverted direction
	 * are recognized as different substituents. Also, otherwise equal chains that connect back to different exit
	 * vectors are also recognized as being different.
	 * After the check for varying substituents and once we have a mapping from exit vector index to R-group index,
	 * we need to exchange the label by a new one with the R-Group index, which should be finally displayed to the user.
	 */
	protected void correctSubstituentRingClosureLabels() {
		if (mSubstituentConnectsBack != null) {
			for (int exitVectorIndex = 0; exitVectorIndex<mScaffold.getExitVectorCount(); exitVectorIndex++) {
				ExitVector exitVector = mScaffold.getExitVector(exitVectorIndex);
				if (exitVector.substituentVaries()
				 && mSubstituentConnectsBack[exitVectorIndex]
				 && mSubstituent[exitVectorIndex] != null) {
					String newIDCode = mScaffold.getOldToNewMap().get(mSubstituent[exitVectorIndex]);
					if (newIDCode != null) {
						mSubstituent[exitVectorIndex] = newIDCode;
					}
					else {
						StereoMolecule s = new IDCodeParser().getCompactMolecule(mSubstituent[exitVectorIndex]);
						for (int atom=0; atom<s.getAllAtoms(); atom++) {
							String label = s.getAtomCustomLabel(atom);
							if (label != null) {
								int backConnectionIndex = Integer.parseInt(label);
								s.setAtomCustomLabel(atom, Integer.toString(mScaffold.getExitVector(backConnectionIndex).getRGroupNo()));
							}
						}
						newIDCode = new Canonizer(s, Canonizer.ENCODE_ATOM_CUSTOM_LABELS).getIDCode();
						mScaffold.getOldToNewMap().put(mSubstituent[exitVectorIndex], newIDCode);
						mSubstituent[exitVectorIndex] = newIDCode;
					}
				}
			}
		}
	}
}

class BackConnection {
	public int substituentAtom,backEndAtom,bondType;

	public BackConnection(int substituentAtom, int backEndAtom, int bondOrder) {
		this.substituentAtom = substituentAtom;
		this.backEndAtom = backEndAtom;
		this.bondType = bondOrder == 3 ? Molecule.cBondTypeTriple : bondOrder == 2 ? Molecule.cBondTypeDouble : Molecule.cBondTypeSingle;
	}
}