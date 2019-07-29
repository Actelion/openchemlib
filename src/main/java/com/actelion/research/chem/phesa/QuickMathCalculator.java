package com.actelion.research.chem.phesa;

import com.actelion.research.chem.PeriodicTable;

public class QuickMathCalculator {
	
	private double[] precalcExp;
	private double[][] precalcPrefactors;
	private static QuickMathCalculator sInstance;
	
	
	private QuickMathCalculator() {
		precalculatePrefactors();
		precalculateExp();
		
	}
	
	public static QuickMathCalculator getInstance() {
		if(sInstance==null) {
			synchronized(QuickMathCalculator.class) {
				if(sInstance==null) {
					sInstance = new QuickMathCalculator();
				}
			}
			
		}
		return sInstance;
	}
	
	private void precalculatePrefactors() {
		precalcPrefactors = new double[54][54];
		for(int i=1;i<54;i++) { //last element taken into account is Iodine
			for(int j=1;j<54;j++) {
				double vdwR1 = PeriodicTable.getElement(i).getVDWRadius();
				double vdwR2 = PeriodicTable.getElement(j).getVDWRadius();
				double alphaSum = MolecularVolume.alpha_pref/(vdwR1*vdwR1) + MolecularVolume.alpha_pref/(vdwR2*vdwR2);
				precalcPrefactors[i][j] = Math.pow((Math.PI/alphaSum), 1.5); 		
			}
		}
	}
	
	private void precalculateExp() {
		precalcExp = new double[1000];
		for(int i=0;i<1000;i++) {
			precalcExp[i] = Math.exp(-i*0.01);
			
		}			
	}
	
	public double quickExp(double c) { //fast approximation of Math.exp using a lookup table and linear interpolation 
		int index = -1*(int)(c*100);
		double exp1 = precalcExp[index];
		double exp2 = precalcExp[index+1];
		double f = -c*100-index;
		double exp = exp1 + f*(exp2-exp1);
		return exp;
	}
	
	public double getPrefactor(int a1, int a2) {
		return precalcPrefactors[a1][a2];
	}

}
