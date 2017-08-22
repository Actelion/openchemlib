package com.actelion.research.chem.coords;

import com.actelion.research.chem.SSSearcherWithIndex;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.DescriptorHandlerFFP512;

public class InventorTemplate {
	private StereoMolecule mFragment;
	private int[] mFFP;

	public InventorTemplate(StereoMolecule fragment) {
		mFragment = fragment;
		}

	public InventorTemplate(StereoMolecule fragment, int[] ffp) {
		mFragment = fragment;
		mFFP = ffp;
		}

	public StereoMolecule getFragment() {
		return  mFragment;
		}

	public int[] getFFP() {
		return mFFP;
		}
	}
