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
import com.actelion.research.chem.potentialenergy.EmpiricalLigandStrain;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;
import com.actelion.research.chem.potentialenergy.TorsionConstraint;

public class IdoScore extends AbstractScoringEngine {
	
	private List<PotentialEnergyTerm> ligStrain;
	private List<PotentialEnergyTerm> constraint;
	private List<PotentialEnergyTerm> interactionEnergy;	
	private BondRotationHelper torsionHelper;
	private int[] receptorAtomTypes;
	private int[] ligAtomTypes;

	public IdoScore(StereoMolecule receptor, Set<Integer> bindingSiteAtoms, int[] receptorAtomTypes, MoleculeGrid grid) {
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
		constraint = new ArrayList<PotentialEnergyTerm>();
		interactionEnergy = new ArrayList<PotentialEnergyTerm>();
		List<Integer> ligAtomTypesList = new ArrayList<>();
		for(int a=0;a<mol.getAtoms();a++) {
			ligAtomTypesList.add(InteractionAtomTypeCalculator.getAtomType(mol, a));
		}
		ligAtomTypes = new int[ligAtomTypesList.size()];
		IntStream.range(0, ligAtomTypes.length).forEach(e -> ligAtomTypes[e] = ligAtomTypesList.get(e));
		updateState();
		torsionHelper = new BondRotationHelper(mol);
		createConstraintFunction();
		initiateInteractionTerms();
		
		
	}
	
	private void initiateInteractionTerms() {
		for(int p : bindingSiteAtoms) {
			for(int l=0;l<candidatePose.getMolecule().getAtoms();l++) {
				PotentialEnergyTerm term = InteractionTerm.create(receptorConf, candidatePose, p,l, receptorAtomTypes, ligAtomTypes);
				if(term!=null)
					interactionEnergy.add(term);
			}
		}
		

	}
	
	private void createConstraintFunction() {
		EmpiricalLigandStrain strain = new EmpiricalLigandStrain(candidatePose,ligAtomTypes,torsionHelper);
		ligStrain.add(strain);
		//add torsions around rotatable bonds to term list and restrain the others to the current value, to prevent distortions
		for(int b=0;b<candidatePose.getMolecule().getBonds();b++) {
			boolean isRotBond = false;
			for(int rotBondIndex=0;rotBondIndex<torsionHelper.getRotatableBonds().length;rotBondIndex++) {
				int rotBond = torsionHelper.getRotatableBonds()[rotBondIndex];
				if(b==rotBond) {
					isRotBond = true;
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
				TorsionConstraint tConstraint = new TorsionConstraint(candidatePose,torsionAtoms,dihedral,5.0);
				constraint.add(tConstraint);
			}
			int at1 = candidatePose.getMolecule().getBondAtom(0, b);
			int at2 = candidatePose.getMolecule().getBondAtom(0, 1);
			Coordinates c1 = candidatePose.getMolecule().getCoordinates(at1);
			Coordinates c2 = candidatePose.getMolecule().getCoordinates(at2);
			double dist = c1.distance(c2);
			BondConstraint bConstraint = new BondConstraint(candidatePose,new int[] {at1,at2},dist);
			constraint.add(bConstraint);
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
                        AngleConstraint aConstraint = new AngleConstraint(candidatePose,new int[] {at2,at1,at3},alpha);
                        constraint.add(aConstraint);
                    }
                }
            }
        }

	}
	

	
	public double getStrain(double[] gradient) {
		double energy = 0.0;
		for(PotentialEnergyTerm term : ligStrain) {
			energy+=term.getFGValue(gradient);
		} 
		return energy;
	}


	@Override
	public double getFGValue(double[] gradient) {
		double energy = 0.0;
		for(PotentialEnergyTerm term : ligStrain) {
			energy+=term.getFGValue(gradient);
		} 
		for(PotentialEnergyTerm term : constraint) {
			energy+=term.getFGValue(gradient);
		} 
		
		for(PotentialEnergyTerm term : interactionEnergy) {
			energy+=term.getFGValue(gradient);
		}
		energy+=getBumpTerm();
		return energy;
			
	}
	
	public double score() {
		double[] gradient = new double[candidatePose.getMolecule().getAllAtoms()*3];
		double energy = 0.0;

		
		for(PotentialEnergyTerm term : interactionEnergy) {
			energy+=term.getFGValue(gradient);
		}

		return energy;
			
	}
	






}
