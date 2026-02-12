package com.actelion.research.chem.coords;

import com.actelion.research.chem.StereoMolecule;

public class InventorTemplate {
	private final StereoMolecule mFragment;
	private final long[] mFFP;
	private double mAVBL;
	private final boolean mKeepAbsoluteOrientation;

	public InventorTemplate(StereoMolecule fragment, long[] ffp, boolean keepAbsoluteOrientation) {
		mFragment = fragment;
		mFFP = ffp;
		mKeepAbsoluteOrientation = keepAbsoluteOrientation;
		}

	public boolean keepAbsoluteOrientation() {
		return mKeepAbsoluteOrientation;
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

	public long[] getFFP() {
		return mFFP;
		}
	}
