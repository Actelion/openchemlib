package com.actelion.research.chem.coords;

import com.actelion.research.chem.StereoMolecule;

public class InventorTemplate {
	private StereoMolecule mFragment;
	private int[] mFFP;
	private double mAVBL;

	public InventorTemplate(StereoMolecule fragment, int[] ffp) {
		mFragment = fragment;
		mFFP = ffp;
		}

	public void normalizeCoordinates() {
		mAVBL = mFragment.getAverageBondLength();
		}

	public double getNormalizedAtomX(int atom) {
		return mFragment.getAtomX(atom) / mAVBL;
		}

	public double getNormalizedAtomY(int atom) {
		return mFragment.getAtomY(atom) / mAVBL;
		}

	public StereoMolecule getFragment() {
		return  mFragment;
		}

	public int[] getFFP() {
		return mFFP;
		}
	}
