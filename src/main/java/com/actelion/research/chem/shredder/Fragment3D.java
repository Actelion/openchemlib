package com.actelion.research.chem.shredder;

import com.actelion.research.chem.conf.TorsionDescriptor;

public class Fragment3D implements Comparable<Fragment3D> {
	private final String mIDCode;
	private String mIDCoords;
	private final TorsionDescriptor mTorsions;
	private final int[] mExitAtoms;

	public Fragment3D(String idcode, String coords, TorsionDescriptor td, int[] exitAtoms) {
		this.mIDCode = idcode;
		this.mIDCoords = coords;
		this.mTorsions = td;
		this.mExitAtoms = exitAtoms;
		}

	public String getIDCode() {
		return mIDCode;
		}

	public String getIDCoordinates() {
		return mIDCoords;
		}

	public int[] getExitAtoms() {
		return mExitAtoms;
		}

	public void setCoordinates(String coords) {
		mIDCoords = coords;
		}

	public boolean equals(Fragment3D f) {
		boolean idcodeDiffers = !mIDCode.equals(f.mIDCode);
		return !(idcodeDiffers || (mTorsions != null && !mTorsions.equals(f.mTorsions)));
		}

	@Override public int compareTo(Fragment3D f) {
		int comparison = mIDCode.compareTo(f.mIDCode);
		if (comparison != 0 || mTorsions == null)	// different structure or no rotatable bonds
			return comparison;

		return mTorsions.compareTo(f.mTorsions);
		}
	}
