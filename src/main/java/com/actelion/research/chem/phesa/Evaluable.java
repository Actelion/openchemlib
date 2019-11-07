package com.actelion.research.chem.phesa;

public interface Evaluable {

	double[] getState();

	double getFGValue(double[] grad);

	void setState(double[] transform);

}
