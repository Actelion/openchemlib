package com.actelion.research.chem.phesa.pharmacophore;


import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;

public interface IPharmacophorePoint {
	
	
	public Coordinates getCenter();
	
	public void updateCoordinates(StereoMolecule mol);
	
	public Coordinates getDirectionality();
	
	public String encode();
	
	public double getSimilarity (IPharmacophorePoint pp);
		
	public int getCenterID();
	
	public void setDirectionality(Coordinates directionality);
	
	

}
