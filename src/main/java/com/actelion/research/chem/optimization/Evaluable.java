package com.actelion.research.chem.optimization;

public interface Evaluable {

	double[] getState();

	double getFGValue(double[] grad);

	void setState(double[] transform);

}
