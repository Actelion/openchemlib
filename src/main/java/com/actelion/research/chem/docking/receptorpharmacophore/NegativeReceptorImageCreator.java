package com.actelion.research.chem.docking.receptorpharmacophore;

import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.phesa.BindingSiteVolume;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.PheSAAlignment;

public class NegativeReceptorImageCreator {
	
	
	
	public static  BindingSiteVolume create(StereoMolecule lig, StereoMolecule rec)  {
		NegativeReceptorImage negImg;
		Matrix rotation;
		Coordinates origCOM;
		Molecule3D ligand = new Molecule3D(lig);
		ligand.ensureHelperArrays(Molecule.cHelperCIP);
		Molecule3D receptor = new Molecule3D(rec);
		receptor.ensureHelperArrays(Molecule.cHelperCIP);
		MolecularVolume molVol = new MolecularVolume(ligand);
		origCOM  = new Coordinates(molVol.getCOM());
		Conformer conf = new Conformer(ligand);
		//align PMI of ligand with the axes of the coordinate system for efficient grid creation
		rotation = PheSAAlignment.preProcess(conf, molVol);
		rotateMols(receptor,ligand,rotation,origCOM);
		negImg = new NegativeReceptorImage(ligand,receptor, 0.4, new Coordinates(4.0,4.0,4.0));
		BindingSiteVolume bsVolume = negImg.calculate();
		double[][] rot = rotation.getTranspose().getArray();
		bsVolume.getPPGaussians().stream().forEach(e -> {
			e.getCenter().rotate(rot);
			e.getCenter().add(origCOM);
		});
		bsVolume.getAtomicGaussians().stream().forEach(e -> {
			e.getCenter().rotate(rot);
			e.getCenter().add(origCOM);
		});

		return bsVolume;
	}
	
	private static void rotateMols(Molecule3D receptor, Molecule3D ligand, Matrix rotation, Coordinates origCOM) {
		double[][] rot = rotation.getArray();
		for(int a=0;a<ligand.getAllAtoms();a++) {
			Coordinates c = ligand.getCoordinates(a);
			c.sub(origCOM);
			c.rotate(rot);
		}
		for(int a=0;a<receptor.getAllAtoms();a++) {
			Coordinates c = receptor.getCoordinates(a);
			c.sub(origCOM);
			c.rotate(rot);
		}
		
		
	}

}
