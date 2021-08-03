package com.actelion.research.chem;

public class SubstructureFilter implements MoleculeFilter {
	private SSSearcher mSearcher;
	private int mMatchMode;

	public SubstructureFilter(StereoMolecule substructure) {
		this(substructure, SSSearcher.cDefaultMatchMode);
		}

	public SubstructureFilter(StereoMolecule substructure, int matchMode) {
		mSearcher = new SSSearcher();
		mSearcher.setFragment(substructure);
		mMatchMode = matchMode;
		}

	@Override
	public boolean moleculeQualifies(StereoMolecule mol) {
		mSearcher.setMolecule(mol);
		return mSearcher.isFragmentInMolecule(mMatchMode);
		}
	}
