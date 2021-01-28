package com.actelion.research.chem;

public class SubstructureFilter implements MoleculeFilter {
	private SSSearcher mSearcher;

	public SubstructureFilter(StereoMolecule substructure) {
		mSearcher = new SSSearcher();
		mSearcher.setFragment(substructure);
		}

	@Override
	public boolean moleculeQualifies(StereoMolecule mol) {
		mSearcher.setMolecule(mol);
		return mSearcher.isFragmentInMolecule();
		}
	}
