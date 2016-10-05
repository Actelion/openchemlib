package com.actelion.research.chem.coords;

/**
 * Created by thomas on 9/28/16.
 */
public class InventorCharge {
	InventorFragment fragment;
	public int atom,charge;

	public InventorCharge(InventorFragment fragment, int atom, int charge) {
		this.fragment = fragment;
		this.atom = atom;
		this.charge = charge;
	}
}
