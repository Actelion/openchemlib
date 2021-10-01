package com.actelion.research.chem.docking.receptorpharmacophore;

import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.transformation.Rotation;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;
import com.actelion.research.chem.alignment3d.transformation.Translation;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.phesa.ShapeVolume;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.PheSAAlignment;

public class NegativeReceptorImageCreator {
	
	
	
	public static  ShapeVolume create(StereoMolecule lig, StereoMolecule rec,TransformationSequence transformation)  {
		NegativeReceptorImage negImg;
		Molecule3D ligand = new Molecule3D(lig);
		ligand.ensureHelperArrays(Molecule.cHelperCIP);
		Molecule3D receptor = new Molecule3D(rec);
		receptor.ensureHelperArrays(Molecule.cHelperCIP);
		MolecularVolume molVol = new MolecularVolume(ligand);
		Coordinates origCOM = new Coordinates(molVol.getCOM());
		Conformer conf = new Conformer(ligand);
		//align PMI of ligand with the axes of the coordinate system for efficient grid creation
		Rotation rotation = molVol.preProcess(conf);
		rotateMols(receptor,ligand,rotation,origCOM);
		negImg = new NegativeReceptorImage(ligand,receptor, 0.4, new Coordinates(4.0,4.0,4.0));
		ShapeVolume bsVolume = negImg.calculate();
		Rotation rot = rotation.getInvert();
		Translation trans = new Translation(new double[] {origCOM.x,origCOM.y,origCOM.z});
	    transformation.addTransformation(rot);
		transformation.addTransformation(trans);



		return bsVolume;
	}
	
	private static void rotateMols(Molecule3D receptor, Molecule3D ligand, Rotation rotation, Coordinates origCOM) {
		for(int a=0;a<ligand.getAllAtoms();a++) {
			Coordinates c = ligand.getCoordinates(a);
			c.sub(origCOM);
			rotation.apply(c);
		}
		for(int a=0;a<receptor.getAllAtoms();a++) {
			Coordinates c = receptor.getCoordinates(a);
			c.sub(origCOM);
			rotation.apply(c);
		}
		
		
	}

}
