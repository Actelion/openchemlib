package com.actelion.research.chem.forcefield;

public interface ForceField {
	public double getTotalEnergy(double[] pos);
	
	public double getTotalEnergy();

	/**
	 * updates the gradient of the ForceField and returns the gradient scale
	 * @return
	 */
	public double updateGradient();

	public double[] getCurrentPositions();

	public void addListener(ForceFieldChangeListener listener);
	
	public void setFixedAtoms(int[] atoms);
	
	
	/**
	 * zeroes out gradient for fixed atoms
	 */
	
	public void zeroGradient();

	public int minimise();

	public void interrupt();
}
