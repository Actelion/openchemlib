package com.actelion.research.chem.descriptor.pharmacophoretree;


import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;

/**
 * calculates features for all atoms in a molecule: 
 * H-bond acceptors/donors, negatively/positively ionizable, aromatic, lipophilic
 * @author joel
 *
 */


public class FeatureCalculator {
	StereoMolecule mol;
	private Set<Integer> donorAtoms;
	private Set<Integer> acceptorAtoms;
	private Set<Integer> negChargeAtoms;
	private Set<Integer> posChargeAtoms;
	private Set<Integer> aromAtoms;
	private Set<Integer> lipoAtoms;
	
	
	public FeatureCalculator(StereoMolecule mol) {
		this.mol = mol;
		donorAtoms = new HashSet<Integer>();
		acceptorAtoms = new HashSet<Integer>();
		negChargeAtoms = new HashSet<Integer>();
		posChargeAtoms = new HashSet<Integer>();
		aromAtoms = new HashSet<Integer>();
		lipoAtoms = new HashSet<Integer>();
		
	}
		
	public void calculate() {
		IonizableGroupDetector2D detector = new IonizableGroupDetector2D(mol);
		detector.detect();
		negChargeAtoms = detector.getNegIonizableAtoms();
		posChargeAtoms = detector.getPosIonizableAtoms();
		getAtomFeatures();
	}
	
	private void getAtomFeatures() {

		for(int i=0;i<mol.getAtoms();i++) {
			if(mol.isAromaticAtom(i))
				aromAtoms.add(i);
			else {
				for(int n=0;n<mol.getConnAtoms(i);n++)
					if(mol.getConnBondOrder(i, n)==2)
						aromAtoms.add(i);
			}

			if (mol.getAtomicNo(i)==7 || mol.getAtomicNo(i)==8) {
				if(PharmacophoreCalculator.isAcceptor(mol,i)) {
					acceptorAtoms.add(i);
				}
				if(PharmacophoreCalculator.isDonorHeavyAtom(mol,i)) {
					donorAtoms.add(i);
				}
				else 
					lipoAtoms.add(i);
			}
			else {
				lipoAtoms.add(i);
			}
	}
	}
	
	public int[][] getAtomFunctionalities() {
		int[][] features = new int[mol.getAtoms()][PharmacophoreCalculator.LIPO_ID+1];
		for(int a=0;a<mol.getAtoms();a++) {
			int[] atomFeatures = features[a];
			if(donorAtoms.contains(a))
				atomFeatures[PharmacophoreCalculator.DONOR_ID]++;
			if (acceptorAtoms.contains(a))
					atomFeatures[PharmacophoreCalculator.ACCEPTOR_ID]++;
			if (negChargeAtoms.contains(a))
				atomFeatures[PharmacophoreCalculator.CHARGE_NEG_ID]++;
			if (posChargeAtoms.contains(a))
				atomFeatures[PharmacophoreCalculator.CHARGE_POS_ID]++;
			if (lipoAtoms.contains(a))
				atomFeatures[PharmacophoreCalculator.LIPO_ID]++;
			if (aromAtoms.contains(a))
				atomFeatures[PharmacophoreCalculator.AROM_ID]++;
		}
		return features;
	}

	
}
