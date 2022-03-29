package com.actelion.research.chem.docking;

import java.util.Random;


import com.actelion.research.calc.Matrix;
import com.actelion.research.calc.SingularValueDecomposition;
import com.actelion.research.chem.AtomFunctionAnalyzer;
import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.AtomAssembler;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.phesa.AtomicGaussian;
import com.actelion.research.chem.phesa.PheSAAlignment;

public class DockingUtils {
	
	private DockingUtils() {};
	
	public static Coordinates getCOM(Conformer conf) {
		int counter = 0;
		Coordinates com = new Coordinates();
		for(int a=0;a<conf.getMolecule().getAtoms();a++) {
			Coordinates coords = conf.getCoordinates(a);
			com.add(coords);
			counter++;
		}
		com.scale(1.0/counter);
		return com;
	}
	
	public static Coordinates getCOM(StereoMolecule conf) {
		int counter = 0;
		Coordinates com = new Coordinates();
		for(int a=0;a<conf.getAtoms();a++) {
			Coordinates coords = conf.getCoordinates(a);
			com.add(coords);
			counter++;
		}
		com.scale(1.0/counter);
		return com;
	}
	
	public static Matrix createInitialOrientation(Conformer conf) {
		Matrix m = calculateMassCovarianceMatrix(conf);
		SingularValueDecomposition svd = new SingularValueDecomposition(m.getArray(),null,null);
		Matrix u = new Matrix(svd.getU());
		double det = u.det();
		if(det<0) {
			u.set(0,1,-u.get(0, 1));
			u.set(1,1,-u.get(1, 1));
			u.set(2,1,-u.get(2, 1));
		}
		PheSAAlignment.rotateMol(conf,u.getArray());
		return u;
		
	}
	
	public static Matrix calculateMassCovarianceMatrix(Conformer conf) {
		Matrix massMatrix = new Matrix(3,3); 
		int counter = 0;
		for (int a=0;a<conf.getMolecule().getAllAtoms();a++){
			Coordinates coords = conf.getCoordinates(a);
			counter++;
			double value = coords.x*coords.x;
			massMatrix.addToElement(0,0,value);
			value = coords.x*coords.y;
			massMatrix.addToElement(0,1,value);
			value = coords.x*coords.z;
			massMatrix.addToElement(0,2,value);
			value = coords.y*coords.y;
			massMatrix.addToElement(1,1,value);
			value = coords.y*coords.z;
			massMatrix.addToElement(1,2,value);
			value = coords.z*coords.z;
			massMatrix.addToElement(2,2,value);	
		}
		massMatrix.set(0,0,massMatrix.get(0,0)/counter);
		massMatrix.set(0,1,massMatrix.get(0,1)/counter);
		massMatrix.set(0,2,massMatrix.get(0,2)/counter);
		massMatrix.set(1,1,massMatrix.get(1,1)/counter);
		massMatrix.set(1,2,massMatrix.get(1,2)/counter);
		massMatrix.set(2,2,massMatrix.get(2,2)/counter);
		massMatrix.set(1,0,massMatrix.get(0,1));
		massMatrix.set(2,0,massMatrix.get(0,2));
		massMatrix.set(2,1,massMatrix.get(1,2));
		
		return massMatrix;
	}
	
	public static Coordinates randomVectorInSphere(Random r) {
		double r1 = 1.0;
		double r2 = 0.0;
		double r3 = 0.0;
		boolean notDone = true;
		while(notDone) {
			r1 = -1+2*r.nextDouble();
			r2 = -1+2*r.nextDouble();
			r3 = -1+2*r.nextDouble();
			if(r1*r1+r2*r2+r3*r3<1.0)
				notDone = false;
		}
		return new Coordinates(r1,r2,r3);
	}
	
	public static void repairMolecule3D(StereoMolecule lig) {
		new Canonizer(lig);
		lig.normalizeAmbiguousBonds();
		repairQuaternaryNitrogen(lig);
		repairCarboxylate(lig);
		addImplicitHydrogens(lig);
	}
	
	private static void repairQuaternaryNitrogen(StereoMolecule mol){
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i) == 7) {
				if(mol.getOccupiedValence(i)==4){
					if(mol.getAtomCharge(i)==0) {
						mol.setAtomCharge(i, 1);
					}
				}
			}
		}
		mol.ensureHelperArrays(Molecule.cHelperRings);
	}
	
	public static void assignLikelyProtonationStates(StereoMolecule mol) {
		for(int a=0;a<mol.getAtoms();a++) {
			if(mol.getAtomicNo(a)==7) {
				if (AtomFunctionAnalyzer.isBasicNitrogen(mol, a))
					mol.setAtomCharge(a, +1);
					
			}
			if(mol.getAtomicNo(a)==8) {
				if (AtomFunctionAnalyzer.isAcidicOxygen(mol, a))
					mol.setAtomCharge(a, -1);
					
			}
		}
		addImplicitHydrogens(mol);
	}
	
	
	public static void addImplicitHydrogens(StereoMolecule mol) {

		 new AtomAssembler(mol).addImplicitHydrogens();


		mol.ensureHelperArrays(Molecule.cHelperNeighbours);

	}
	
	private static void repairCarboxylate(StereoMolecule mol){
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i) == 8) {
				if(mol.getOccupiedValence(i)==1){
					if(mol.getAtomCharge(i)==0) {
						mol.setAtomCharge(i, -1);
					}
				}
			}
		}
		mol.ensureHelperArrays(Molecule.cHelperRings);
	}
	

}
