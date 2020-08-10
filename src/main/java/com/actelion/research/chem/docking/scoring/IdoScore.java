package com.actelion.research.chem.docking.scoring;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.AbstractMap.SimpleEntry;
import java.util.stream.IntStream;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.torsionstrain.StatisticalTorsionPotential;
import com.actelion.research.chem.conf.torsionstrain.StatisticalTorsionTerm;
import com.actelion.research.chem.docking.scoring.idoscore.InteractionTerm;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.chem.interactionstatistics.InteractionDistanceStatistics;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.potentialenergy.AngleConstraint;
import com.actelion.research.chem.potentialenergy.BondConstraint;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;
import com.actelion.research.chem.potentialenergy.TorsionConstraint;

public class IdoScore extends AbstractScoringEngine {
	
	private List<PotentialEnergyTerm> ligStrain;
	private List<PotentialEnergyTerm> interactionEnergy;	
	private BondRotationHelper torsionHelper;
	private int[] receptorAtomTypes;
	private int[] ligAtomTypes;
	private List<int[]> ligAtomPairs; //separated by more than 3 bonds for internal strain

	public IdoScore(Molecule3D receptor, Set<Integer> bindingSiteAtoms, int[] receptorAtomTypes, MoleculeGrid grid) {
		super(receptor, bindingSiteAtoms, grid);
		this.receptorAtomTypes = receptorAtomTypes;
	}
	
	public void init(Conformer candidatePose) {
		this.candidatePose = candidatePose;
		InteractionDistanceStatistics.getInstance().initialize();
		StatisticalTorsionPotential.getInstance().initialize();
		StereoMolecule mol = candidatePose.getMolecule();
		state = new double[3*mol.getAllAtoms()];
		ligStrain = new ArrayList<PotentialEnergyTerm>();
		List<Integer> ligAtomTypesList = new ArrayList<>();
		for(int a=0;a<mol.getAtoms();a++) {
			ligAtomTypesList.add(InteractionAtomTypeCalculator.getAtomType(mol, a));
		}
		ligAtomTypes = new int[ligAtomTypesList.size()];
		IntStream.range(0, ligAtomTypes.length).forEach(e -> ligAtomTypes[e] = ligAtomTypesList.get(e));
		updateState();
		torsionHelper = new BondRotationHelper(mol);
		ligAtomPairs = new ArrayList<int[]>();
		createStrainFunction();
		initiateInteractionTerms();
		
		
	}
	
	private void initiateInteractionTerms() {

		interactionEnergy = new ArrayList<PotentialEnergyTerm>();

		for(int p : bindingSiteAtoms) {
			for(int l=0;l<candidatePose.getMolecule().getAtoms();l++) {
				interactionEnergy.add(InteractionTerm.create(receptorConf, candidatePose, p,l, receptorAtomTypes, ligAtomTypes));
			}
		}
		

	}
	
	private void createStrainFunction() {
		findLigAtomPairs();
		for(int[] pair : ligAtomPairs) {
			int at1 = pair[0];
			int at2 = pair[1];
			InteractionTerm term = InteractionTerm.create(candidatePose,candidatePose,at1,at2,ligAtomTypes,ligAtomTypes);
			if(term==null)
				continue;
			ligStrain.add(term);
		}
		//add torsions around rotatable bonds to term list and restrain the others to the current value, to prevent distortions
		for(int b=0;b<candidatePose.getMolecule().getBonds();b++) {
			boolean isRotBond = false;
			for(int rotBondIndex=0;rotBondIndex<torsionHelper.getRotatableBonds().length;rotBondIndex++) {
				int rotBond = torsionHelper.getRotatableBonds()[rotBondIndex];
				if(b==rotBond) {
					isRotBond = true;
					int[] torsionAtoms = torsionHelper.getTorsionAtoms()[rotBondIndex];
					String torsionID = torsionHelper.getTorsionIDs()[rotBondIndex];
					StatisticalTorsionTerm term = StatisticalTorsionTerm.create(candidatePose,torsionAtoms,torsionID);
					if(term==null)
						continue;
					ligStrain.add(term);
					break;
				}
			}
			if(!isRotBond) {
				int at2 = candidatePose.getMolecule().getBondAtom(0, b);
				int at3 = candidatePose.getMolecule().getBondAtom(1, b);
				if(candidatePose.getMolecule().getConnAtoms(at2)==1 || candidatePose.getMolecule().getConnAtoms(at3)==1) //terminal bond
					continue;
				int at1 = candidatePose.getMolecule().getConnAtom(at2, 0) == at3 ? candidatePose.getMolecule().getConnAtom(at2, 1) 
						: candidatePose.getMolecule().getConnAtom(at2, 0);
				int at4 = candidatePose.getMolecule().getConnAtom(at3, 0) == at2 ? candidatePose.getMolecule().getConnAtom(at3, 1) :
					candidatePose.getMolecule().getConnAtom(at3, 0);
				int[] torsionAtoms = new int[] {at1,at2,at3,at4};
				Coordinates c1 = candidatePose.getCoordinates(at1);
				Coordinates c2 = candidatePose.getCoordinates(at2);
				Coordinates c3 = candidatePose.getCoordinates(at3);
				Coordinates c4 = candidatePose.getCoordinates(at4);
				double dihedral = 360.0*Coordinates.getDihedral(c1, c2, c3, c4)/(2*Math.PI);
				TorsionConstraint constraint = new TorsionConstraint(candidatePose,torsionAtoms,dihedral,5.0);
				ligStrain.add(constraint);
			}
			int at1 = candidatePose.getMolecule().getBondAtom(0, b);
			int at2 = candidatePose.getMolecule().getBondAtom(0, 1);
			Coordinates c1 = candidatePose.getMolecule().getCoordinates(at1);
			Coordinates c2 = candidatePose.getMolecule().getCoordinates(at2);
			double dist = c1.distance(c2);
			BondConstraint constraint = new BondConstraint(candidatePose,new int[] {at1,at2},dist);
			ligStrain.add(constraint);
		}
		// create angle constraints
        for (int at1=0; at1<candidatePose.getMolecule().getAtoms(); at1++) {
            if (candidatePose.getMolecule().getConnAtoms(at1) > 1) {
                for (int i=0; i<candidatePose.getMolecule().getConnAtoms(at1); i++) {
                    int at2 = candidatePose.getMolecule().getConnAtom(at1, i);
                    for (int k=i+1; k<candidatePose.getMolecule().getConnAtoms(at1); k++) {
                        int at3 = candidatePose.getMolecule().getConnAtom(at1, k);
                        Coordinates c1 = candidatePose.getMolecule().getCoordinates(at1);
                        Coordinates c2 = candidatePose.getMolecule().getCoordinates(at2);
                        Coordinates c3 = candidatePose.getMolecule().getCoordinates(at3);
                        Coordinates v1 = c2.subC(c1);
                        Coordinates v2 = c3.subC(c1);
                        v1.unit();
                        v2.unit();
                        double alpha = Math.acos(v1.dot(v2));
                        alpha = 180.0*alpha/Math.PI;
                        AngleConstraint constraint = new AngleConstraint(candidatePose,new int[] {at2,at1,at3},alpha);
                        ligStrain.add(constraint);
                    }
                }
            }
        }

	}
	
	private void findLigAtomPairs() {
		Set<SimpleEntry<Integer,Integer>> invalidPairs  = new HashSet<SimpleEntry<Integer,Integer>>();
		StereoMolecule mol = candidatePose.getMolecule();
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
	public double getFGValue(double[] gradient) {
		double energy = 0.0;
		for(PotentialEnergyTerm term : ligStrain) {
			energy+=term.getFGValue(gradient);
		}
		for(PotentialEnergyTerm term : interactionEnergy) {
			energy+=term.getFGValue(gradient);
		}
		energy+=getBumpTerm();
		return energy;
			
	}





}
