package com.actelion.research.chem.potentialenergy;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.AbstractMap.SimpleEntry;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.torsionstrain.StatisticalTorsionPotential;
import com.actelion.research.chem.conf.torsionstrain.StatisticalTorsionTerm;
import com.actelion.research.chem.docking.scoring.idoscore.InteractionTerm;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.chem.interactionstatistics.InteractionDistanceStatistics;


public class EmpiricalLigandStrain implements PotentialEnergyTerm {
	
	private Conformer conf;
	private List<PotentialEnergyTerm> strain;
	private BondRotationHelper torsionHelper;
	private List<int[]> ligAtomPairs; 
	private int[] atomTypes;

	public EmpiricalLigandStrain(Conformer conf, int[] atomTypes, BondRotationHelper torsionHelper) {
		this.conf = conf;
		this.atomTypes = atomTypes;
		this.torsionHelper = torsionHelper;
		
		init();
	}
	
	private void init() {
		StereoMolecule mol = conf.getMolecule();
		strain = new ArrayList<PotentialEnergyTerm>();
		ligAtomPairs = new ArrayList<int[]>();
		List<Integer> ligAtomTypesList = new ArrayList<>();
		for(int a=0;a<mol.getAtoms();a++) 
			ligAtomTypesList.add(InteractionAtomTypeCalculator.getAtomType(mol, a));

		findLigAtomPairs();
		for(int[] pair : ligAtomPairs) {
			int at1 = pair[0];
			int at2 = pair[1];
			InteractionTerm term = InteractionTerm.create(conf,conf,at1,at2,atomTypes,atomTypes);
			if(term==null)
				continue;
			strain.add(term);
		}
			//add torsions around rotatable bonds to term list and restrain the others to the current value, to prevent distortions
			for(int b=0;b<conf.getMolecule().getBonds();b++) {
				for(int rotBondIndex=0;rotBondIndex<torsionHelper.getRotatableBonds().length;rotBondIndex++) {
					int rotBond = torsionHelper.getRotatableBonds()[rotBondIndex];
					if(b==rotBond) {
						int[] torsionAtoms = torsionHelper.getTorsionAtoms()[rotBondIndex];
						String torsionID = torsionHelper.getTorsionIDs()[rotBondIndex];
						StatisticalTorsionTerm term = StatisticalTorsionTerm.create(conf,torsionAtoms,torsionID);
						if(term==null)
							continue;
						strain.add(term);
						break;
					}
				}
			}
		
	}
	
	private void findLigAtomPairs() {
		Set<SimpleEntry<Integer,Integer>> invalidPairs  = new HashSet<SimpleEntry<Integer,Integer>>();
		StereoMolecule mol = conf.getMolecule();
		for(int a=0;a<mol.getAtoms();a++) {
			for(int i=0;i<mol.getConnAtoms(a);i++) {
				int aa = mol.getConnAtom(a, i);
				SimpleEntry<Integer,Integer> entry = a<aa ? new SimpleEntry<>(a,aa) : new SimpleEntry<>(a,aa);
				invalidPairs.add(entry);
				for(int j=0;j<mol.getConnAtoms(aa);j++) {
					int aaa = mol.getConnAtom(aa, j);
					entry = a<aaa ? new SimpleEntry<>(a,aaa) : new SimpleEntry<>(a,aaa);
					invalidPairs.add(entry);
					for(int k=0;k<mol.getConnAtoms(aaa);k++) {
						int aaaa = mol.getConnAtom(aaa, k);
						entry = a<aaaa ? new SimpleEntry<>(a,aaaa) : new SimpleEntry<>(a,aaaa);
						invalidPairs.add(entry);
					}
				}
			}
		}
		for(int i=0;i<mol.getAtoms();i++) {
			for(int j=i+1;j<mol.getAtoms();j++) {
				SimpleEntry<Integer,Integer> entry =new SimpleEntry<>(i,j);
				if(invalidPairs.contains(entry))
					continue;
				else 
					ligAtomPairs.add(new int[] {i,j});
					
			}
		}
	}
	

	@Override
	public double getFGValue(double[] grad) {
		double energy = 0.0;
		for(PotentialEnergyTerm term : strain) {
			energy+=term.getFGValue(grad);
		}
		return energy;
	}


	
	

}
