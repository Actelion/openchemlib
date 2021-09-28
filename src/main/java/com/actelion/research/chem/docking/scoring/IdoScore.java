package com.actelion.research.chem.docking.scoring;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.AbstractMap.SimpleEntry;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.torsionstrain.StatisticalTorsionPotential;
import com.actelion.research.chem.conf.torsionstrain.StatisticalTorsionTerm;
import com.actelion.research.chem.docking.LigandPose;
import com.actelion.research.chem.docking.scoring.idoscore.InteractionTerm;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.chem.interactionstatistics.InteractionDistanceStatistics;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.potentialenergy.AngleConstraint;
import com.actelion.research.chem.potentialenergy.BondConstraint;
import com.actelion.research.chem.potentialenergy.EmpiricalLigandStrain;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;
import com.actelion.research.chem.potentialenergy.TorsionConstraint;

public class IdoScore extends AbstractScoringEngine {
	
	
	private static final double STRAIN_CUTOFF = 20.0;
	private List<PotentialEnergyTerm> ligStrain;
	private List<PotentialEnergyTerm> constraint;
	private List<PotentialEnergyTerm> interactionEnergy;	
	private BondRotationHelper torsionHelper;
	private int[] receptorAtomTypes;
	private int[] ligAtomTypes;
	private ForceFieldMMFF94 ff;
	private double e0;

	public IdoScore(StereoMolecule receptor, Set<Integer> bindingSiteAtoms, int[] receptorAtomTypes, MoleculeGrid grid) {
		super(receptor, bindingSiteAtoms, grid);
		this.receptorAtomTypes = receptorAtomTypes;
	}
	
	@Override
	public void init(LigandPose candidatePose, double e0) {
		this.candidatePose = candidatePose;
		this.e0 = e0;
		ligStrain = new ArrayList<PotentialEnergyTerm>();
		constraint = new ArrayList<PotentialEnergyTerm>();
		interactionEnergy = new ArrayList<PotentialEnergyTerm>();
		List<Integer> ligAtomTypesList = new ArrayList<>();
		StereoMolecule mol = candidatePose.getLigConf().getMolecule();
		for(int a=0;a<mol.getAtoms();a++) {
			ligAtomTypesList.add(InteractionAtomTypeCalculator.getAtomType(mol, a));
		}
		ligAtomTypes = new int[ligAtomTypesList.size()];

		IntStream.range(0, ligAtomTypes.length).forEach(e -> ligAtomTypes[e] = ligAtomTypesList.get(e));

		torsionHelper = new BondRotationHelper(mol);
		
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		ffOptions.put("dielectric constant", 80.0);
		
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		ff = new ForceFieldMMFF94(mol, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);


		initiateInteractionTerms();
		
		
	}
	
	private void initiateInteractionTerms() {
		for(int p : bindingSiteAtoms) {
			for(int l=0;l<candidatePose.getLigConf().getMolecule().getAtoms();l++) {
				PotentialEnergyTerm term = InteractionTerm.create(receptorConf, candidatePose.getLigConf(), p,l, receptorAtomTypes, ligAtomTypes);
				if(term!=null)
					interactionEnergy.add(term);
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
		for(PotentialEnergyTerm term : interactionEnergy) {
			energy+=term.getFGValue(gradient);
		}

		
		energy+=getBumpTerm();
		ff.setState(candidatePose.getCartState());
		double ffEnergy = ff.getTotalEnergy();
		if(ffEnergy-e0>STRAIN_CUTOFF) {
			ff.addGradient(gradient);
			energy+=(ffEnergy-e0);
		}


		return energy;
	}
	
	@Override 
	public void updateState() {
		ff.setState(candidatePose.getState());
	}
		
	@Override
	public double getScore() {
		double[] gradient = new double[candidatePose.getLigConf().getMolecule().getAllAtoms()*3];
		double energy = getBumpTerm();

		
		for(PotentialEnergyTerm term : interactionEnergy) {
			energy+=term.getFGValue(gradient);
		}

		return energy;
			
	}

	@Override
	public Map<String, Double> getContributions() {
		// TODO Auto-generated method stub
		return null;
	}









}
