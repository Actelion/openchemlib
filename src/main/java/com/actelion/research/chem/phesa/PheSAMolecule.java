package com.actelion.research.chem.phesa;

import com.actelion.research.chem.StereoMolecule;
import java.util.ArrayList;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;

/** 
 * @version: 1.0, February 2018
 * Author: J. Wahl
 * contains all information about a molecule (atoms, connectivity) and its shape including flexibilty (Molecular Volumes of a conformational ensemble)
 * the Coordinates of the hydrogens are contained and therefore any molecular conformer can be back-calculated on request
 * the heavy atom coordinates are implicitly given as the centers of the atomic gaussians 
 * 
*/
public class PheSAMolecule {
	private StereoMolecule mol;

	// As many objects as conformers.
	private ArrayList<MolecularVolume> shape;

	
	
	public PheSAMolecule() {
		this.mol = new StereoMolecule();
		this.shape = new ArrayList<MolecularVolume>();
	}
	
	public PheSAMolecule(StereoMolecule mol, MolecularVolume shape) {
		this.mol = mol;
		this.shape = new ArrayList<MolecularVolume>();
		this.shape.add(shape);
	}
	
	public PheSAMolecule(StereoMolecule mol, ArrayList<MolecularVolume> shape) {
		this.mol = mol;
		this.shape = shape;


	}

	/**
	 * Returns the corresponding conformer of a molecular volume
	 * @param molVol
	 * @return
	 */
	public StereoMolecule getConformer(MolecularVolume molVol) {
		int nrOfAtoms = mol.getAllAtoms();
		StereoMolecule conformer = new StereoMolecule(mol);
		int hydrogenCounter = 0;
		ArrayList<Coordinates> hydrogens = molVol.getHydrogens();
		for(int i=0;i<nrOfAtoms;i++) {
			if(mol.getAtomicNo(i)==1){
				conformer.getCoordinates(i).set(hydrogens.get(hydrogenCounter));
				hydrogenCounter+=1;
			} 

		for(int j=0;j<molVol.getAtomicGaussians().size();j++) {
			int atomId =molVol.getAtomicGaussians().get(j).getAtomId();
			conformer.getCoordinates(atomId).set(molVol.getAtomicGaussians().get(j).getCenter());
		}
		}
		conformer.ensureHelperArrays(Molecule.cHelperNeighbours);
		return conformer;
	}
	
	public StereoMolecule getConformer(int index) {
		return getConformer(shape.get(index));
	}
	
	public StereoMolecule getMolecule() {
		return this.mol;
	}
	
	public ArrayList<MolecularVolume> getVolumes() {
		return this.shape;
	}
	
	
	
}
