package com.actelion.research.chem.phesa;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

public interface IPharmacophorePoint {
	
	
	public Coordinates getCenter();
	
	public void updateCoordinates(StereoMolecule mol);
	
	public void updateCoordinates(Conformer conf);
	
	public Coordinates getDirectionality();
	
	public String encode();
	
	public double getSimilarity (IPharmacophorePoint pp);
	
	public int getInteractionClass();
	
	public int getCenterID();
	
	public void setDirectionality(Coordinates directionality);
	

	
	


}
