package com.actelion.research.chem.forcefield;

public interface ForceField {
	public double getTotalEnergy(double[] pos);

	/**
	 * updates the gradient of the ForceField and returns the gradient scale
	 * @return
	 */
	public double updateGradient();

	public double[] getCurrentPositions();

	public double[] getCurrentPositionsMapped(); //remaps the order of the coordinates to originally submitted ones, before the reordering by ensureHelperArrays

	public void addListener(ForceFieldChangeListener listener);

	public int minimise();

	public void interrupt();
}
