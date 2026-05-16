package org.openmolecules.chem.interaction;

import com.actelion.research.chem.*;

public abstract class AtomClassifier implements InteractionAtomClassifier {
	public static final int TYPE_UNKNOWN = 0;

	private SSSearcher[] mFragmentSearcher;
	private StereoMolecule[] mFragment;
	private int[][] mFlaggedAtom;

	protected void initializeStaticStuff(String[][] idcodesWithTypeNames, StereoMolecule[] fragment, int[][] flaggedAtom) {
		for (int i=0; i<fragment.length; i++) {
			fragment[i] = new IDCodeParserWithoutCoordinateInvention().getCompactMolecule(idcodesWithTypeNames[i][0]);
			if (fragment[i] == null)
				System.out.println("$$$ ERROR: Parsing ligand idcode: " + idcodesWithTypeNames[i][0]);
			else
				fragment[i].setName(idcodesWithTypeNames[i][1]);

			fragment[i].ensureHelperArrays(Molecule.cHelperRings);
			int count = 0;
			for (int atom=0; atom<fragment[i].getAtoms(); atom++)
				if ("]*".equals(fragment[i].getAtomCustomLabel(atom)))
					count++;
			flaggedAtom[i] = new int[count];
			count = 0;
			for (int atom=0; atom<fragment[i].getAtoms(); atom++)
				if ("]*".equals(fragment[i].getAtomCustomLabel(atom)))
					flaggedAtom[i][count++] = atom;
		}
	}

	protected void initialize(StereoMolecule[] fragment, int[][] flaggedAtom) {
		mFragment = fragment;
		mFlaggedAtom = flaggedAtom;
		mFragmentSearcher = new SSSearcher[mFragment.length];
		for (int i=0; i<mFragment.length; i++) {
			if (mFragment[i] != null) {
				mFragmentSearcher[i] = new SSSearcher();
				mFragmentSearcher[i].setFragment(mFragment[i]);
			}
		}
	}

	public int[] classifyAtoms(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		int[] atomType = new int[mol.getAtoms()];

		for (int i=0; i<mFragmentSearcher.length; i++) {
			if (mFragmentSearcher[i] != null) {
				mFragmentSearcher[i].setMolecule(mol);
				if (mFragmentSearcher[i].findFragmentInMolecule(SSSearcher.cCountModeRigorous, SSSearcher.cDefaultMatchMode) != 0)
					for (int[] match : mFragmentSearcher[i].getMatchList())
						for (int atom : mFlaggedAtom[i])
							atomType[match[atom]] = getAtomTypeIndex(mFragment[i].getName());
			}
		}

		return atomType;
	}

	abstract public int getAtomTypeIndex(String typeName);
}
