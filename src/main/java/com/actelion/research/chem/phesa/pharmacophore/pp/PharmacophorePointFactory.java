package com.actelion.research.chem.phesa.pharmacophore.pp;

import com.actelion.research.chem.StereoMolecule;

public class PharmacophorePointFactory {
	
	public static IPharmacophorePoint fromString(String ppString, StereoMolecule mol) {
		String type = ppString.split(" ")[0];
		if(type.equals("a"))
			return AcceptorPoint.fromString(ppString, mol);
		else if(type.equals("d"))
			return DonorPoint.fromString(ppString, mol);
		else if(type.equals("i"))
			return ChargePoint.fromString(ppString, mol);
		else if(type.equals("r"))
			return AromRingPoint.fromString(ppString, mol);
		else if(type.equals("e"))
			return ExitVectorPoint.fromString(ppString, mol);
		else if(type.equals("s"))
			return SimplePharmacophorePoint.fromString(ppString, mol);
		
		else 
			return null;
		
	}

}
