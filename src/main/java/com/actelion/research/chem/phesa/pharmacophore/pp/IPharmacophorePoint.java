package com.actelion.research.chem.phesa.pharmacophore.pp;


import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.alignment3d.transformation.Rotation;
import com.actelion.research.chem.alignment3d.transformation.Scaling;
import com.actelion.research.chem.alignment3d.transformation.Transformation;
import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;

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
	
	public void updateCoordinates(Coordinates[] coords);
	
	public Coordinates getDirectionality();
	
	public Coordinates getRotatedDirectionality(double[][] m, double scaleFactor);
	
	public String encode();
	
	public double getSimilarity (IPharmacophorePoint pp);
		
	public int getCenterID();
	
	public void setCenterID(int id);
	
	public void setDirectionality(Coordinates directionality);
	
	public void updateAtomIndices(int[] map);

	int [] getAtomIndices();
	
	public IPharmacophorePoint copyPharmacophorePoint();
	
	public void getDirectionalityDerivativeCartesian(double[] grad, double[] v, Coordinates di, double sim);
	
	public int getFunctionalityIndex();
	
	default public void applyTransformation(Transformation transform) {
		transform.apply(getCenter());
		if(transform instanceof Rotation) {
			Rotation rot = (Rotation) transform;
			Coordinates direc = getRotatedDirectionality(rot.getRotation(),1.0);
			setDirectionality(direc);
		}
		else if(transform instanceof Scaling) {
			Scaling scaling = (Scaling) transform;
			double[][] rot = new double[][] {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
			Coordinates direc = getRotatedDirectionality(rot,scaling.getScalingFactor());
			setDirectionality(direc);
		}
	}
	
	
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
