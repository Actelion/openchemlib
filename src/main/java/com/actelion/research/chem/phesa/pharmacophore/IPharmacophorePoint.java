package com.actelion.research.chem.phesa.pharmacophore;


import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

public interface IPharmacophorePoint {
	
	
	public enum Functionality {ACCEPTOR(PharmacophoreCalculator.ACCEPTOR_ID), DONOR(PharmacophoreCalculator.DONOR_ID), 
		NEG_CHARGE(PharmacophoreCalculator.CHARGE_NEG_ID), POS_CHARGE(PharmacophoreCalculator.CHARGE_POS_ID),
		AROM_RING(PharmacophoreCalculator.AROM_RING_ID),EXIT_VECTOR(PharmacophoreCalculator.EXIT_VECTOR_ID);
		private final int index;
		Functionality(int index) {
			this.index = index;
		}
		public int getIndex() {
			return this.index;
		}
		
	}
	
	
	public Coordinates getCenter();
	
	public void updateCoordinates(StereoMolecule mol);
	
	public void updateCoordinates(Conformer conf);
	
	public Coordinates getDirectionality();
	
	public String encode();
	
	public double getSimilarity (IPharmacophorePoint pp);
		
	public int getCenterID();
	
	public void setDirectionality(Coordinates directionality);
	
	public void updateAtomIndeces(int[] map);
	
	public IPharmacophorePoint copyPharmacophorePoint();
	
	public void getDirectionalityDerivativeCartesian(double[] grad, double[] v, Coordinates di, double sim);
	
	public int getFunctionalityIndex();
	
	
	default public double getVectorSimilarity(IPharmacophorePoint pp2,Coordinates directionalityMod) {
		double vectorSim = 0.0;
		vectorSim = getDirectionality().dot(directionalityMod);
		if (vectorSim<0.0) {
			vectorSim = 0.0;
		}
		return vectorSim;
	}
		

	default public double getVectorSimilarity(IPharmacophorePoint pp2) {
		return getVectorSimilarity(pp2, pp2.getDirectionality());
	}
	
	
	
	
}
