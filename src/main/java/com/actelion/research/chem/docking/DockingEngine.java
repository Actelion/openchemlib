package com.actelion.research.chem.docking;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.openmolecules.chem.conf.gen.ConformerGenerator;

import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.docking.scoring.AbstractScoringEngine;
import com.actelion.research.chem.docking.scoring.ChemPLP;
import com.actelion.research.chem.docking.scoring.IdoScore;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.optimization.OptimizerLBFGS;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.PheSAAlignment;
import com.actelion.research.chem.phesa.PheSAMolecule;

public class DockingEngine {
	
	public enum ScoringFunction {CHEMPLP,IDOSCORE;}
	private static final int DEFAULT_NR_MC_STEPS = 100;
	private static final int DEFAULT_START_POSITIONS = 5;
	private static final double BOLTZMANN_FACTOR = 1.2; //as for AutoDock Vina
	public static final double GRID_DIMENSION = 6.0;
	public static final double GRID_RESOLUTION = 0.5;
	

	private Matrix rotation; //for initial prealignment to native ligand
	private Coordinates origCOM;
	private int mcSteps;
	private Random random;
	private AbstractScoringEngine engine;
	private int startPositions;
	
	public DockingEngine(Molecule3D receptor, Molecule3D nativeLig, int mcSteps, int startPositions,
			ScoringFunction scoringFunction) {
		Molecule3D nativeLigand = new Molecule3D(nativeLig);
		nativeLigand.ensureHelperArrays(Molecule.cHelperCIP);
		MolecularVolume molVol = new MolecularVolume(nativeLigand);
		origCOM  = new Coordinates(molVol.getCOM());
		Conformer conf = new Conformer(nativeLigand);
		rotation = PheSAAlignment.preProcess(conf, molVol);
		this.startPositions = startPositions;
		preprocess(receptor,nativeLigand);
		
		MoleculeGrid grid = new MoleculeGrid(nativeLigand,DockingEngine.GRID_RESOLUTION,
				new Coordinates(DockingEngine.GRID_DIMENSION,DockingEngine.GRID_DIMENSION,
						DockingEngine.GRID_DIMENSION));
		
		Set<Integer> bindingSiteAtoms = new HashSet<Integer>();
		
		if(scoringFunction==ScoringFunction.CHEMPLP) {
			DockingEngine.getBindingSiteAtoms(receptor, bindingSiteAtoms, grid, true);
			engine = new ChemPLP(receptor,bindingSiteAtoms,grid);
		}
		else if(scoringFunction==ScoringFunction.IDOSCORE) {
			DockingEngine.getBindingSiteAtoms(receptor, bindingSiteAtoms, grid, false);
			engine = new IdoScore(receptor,bindingSiteAtoms, getReceptorAtomTypes(receptor),grid);
		}
		
		this.mcSteps = mcSteps;
		this.random = new Random(12354L);


	}
	
	public DockingEngine(Molecule3D receptor, Molecule3D nativeLigand) {
		this(receptor,nativeLigand,DEFAULT_NR_MC_STEPS,DEFAULT_START_POSITIONS,ScoringFunction.IDOSCORE);
	}
	
	
	
	public StereoMolecule dockMolecule(StereoMolecule mol) {
		/*
		DescriptorHandlerShape dhs = new DescriptorHandlerShape();
		PheSAMolecule pheSAMol = dhs.createDescriptor(mol);
		dhs.getSimilarity(refVol,pheSAMol);
		StereoMolecule alignedMol = dhs.getPreviousAlignment()[1];
		for(int a=0;a<mol.getAllAtoms();a++) {
			mol.setAtomX(a, alignedMol.getAtomX(a));
			mol.setAtomY(a, alignedMol.getAtomY(a));
			mol.setAtomZ(a, alignedMol.getAtomZ(a));
		}
		PheSAAlignment.rotateMol(mol, rotation);
		mol.translate(origCOM.x,origCOM.y,origCOM.z);
		*/

		ConformerGenerator confGen = new ConformerGenerator(LigandPose.SEED,true);
		confGen.initializeConformers(mol);
		Conformer bestPose = null;
		double bestEnergy = Double.MAX_VALUE;
		
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		ffOptions.put("dielectric constant", 10.0);
		
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		for(int i=0;i<startPositions;i++) {
			StereoMolecule conf = new StereoMolecule(mol);
			conf.ensureHelperArrays(Molecule.cHelperCIP);
			confGen.getNextConformerAsMolecule(conf);
			ForceFieldMMFF94 ff = new ForceFieldMMFF94(conf, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
			ff.minimise();
			Conformer ligConf = new Conformer(conf);
			Coordinates com = DockingUtils.getCOM(ligConf);
			Coordinates translate = com.scale(-1.0);
			for(int a=0;a<ligConf.getMolecule().getAllAtoms();a++) {
				Coordinates c = ligConf.getCoordinates(a);
				c.add(translate);
			}
			for(double[] transform : PheSAAlignment.initialTransform(1)) {
				Conformer newLigConf = new Conformer(ligConf);
				PheSAAlignment.rotateMol(newLigConf, transform);
				
				LigandPose pose = initiate(newLigConf);
				double energy = mcSearch(pose);
				if(energy<bestEnergy) {
					bestEnergy = energy;
					bestPose = pose.getLigConf();
				}
			}
		}
		System.out.println(bestEnergy);
		StereoMolecule best = bestPose.toMolecule();
		double[][] rot = rotation.getTranspose().getArray();
		PheSAAlignment.rotateMol(best, rot);
		PheSAAlignment.translateMol(best, new double[] {origCOM.x, origCOM.y, origCOM.z} );
		return best;
		
	}
	
	private double mcSearch(LigandPose pose) {
		List<Conformer> allPoses = new ArrayList<>();
		double[] bestState = new double[pose.getState().length];
		double[] oldState = new double[pose.getState().length];
		double[] state = new double[pose.getState().length];
		double bestEnergy = -Float.MAX_VALUE;
		double oldEnergy = -Float.MAX_VALUE;
		double energy = -Float.MAX_VALUE;
		OptimizerLBFGS optimizer = new OptimizerLBFGS(200,0.001);
		oldState = optimizer.optimize(pose);

		bestState = oldState;
		oldEnergy = pose.getFGValue(new double[bestState.length]);

		bestEnergy = oldEnergy;
		for(int i=0;i<mcSteps;i++) {
			pose.randomPerturbation();
			state = optimizer.optimize(pose);
			energy = pose.getFGValue(new double[bestState.length]);

			if(energy<bestEnergy) {
				bestEnergy = energy;
				bestState = state;
			}
			if(energy<oldEnergy) { //accept new pose
				oldEnergy = energy;
				oldState = state;
			}
			else {
				double dE = energy-oldEnergy;
				double randNr = random.nextDouble();
				double probability = Math.exp(-dE/BOLTZMANN_FACTOR); // Metropolis-Hastings
				if(randNr<probability) { //accept
					oldEnergy = energy;
					state = oldState;
				}
				else { //reject
					pose.setState(oldState);
				}
				
			}
			
		}

		pose.setState(bestState);
		return bestEnergy;
		
	}
	
	private LigandPose initiate(Conformer ligConf) {
		LigandPose pose = new LigandPose(ligConf, engine);
		
		return pose;
		
		
	}
	
	public static void getBindingSiteAtoms(StereoMolecule receptor, Set<Integer> bindingSiteAtoms, MoleculeGrid grid,
			boolean includeHydrogens) {
		int[] gridSize = grid.getGridSize();
		int atoms = receptor.getAtoms();
		if(includeHydrogens)
			atoms = receptor.getAllAtoms();
		
		for(int i=0;i<atoms;i++) {
			int[] gridC = grid.getGridCoordinates(receptor.getCoordinates(i));
			int x = gridC[0];
			int y = gridC[1];
			int z = gridC[2];	
			if(x>0 && x<gridSize[0]) {
				if(y>0 && y<gridSize[1]) {
					if(z>0 && z<gridSize[2]) {
						bindingSiteAtoms.add(i);
					}
			}
		}
		}
		
	}
	
	public static int[] getReceptorAtomTypes(StereoMolecule receptor) {
		int[] receptorAtomTypes = new int[receptor.getAtoms()];
		for(int i=0;i<receptor.getAtoms();i++) {
			receptorAtomTypes[i] = InteractionAtomTypeCalculator.getAtomType(receptor, i);
		}
		return receptorAtomTypes;
		
		
	}
	/**
	 * rotates receptor and ligand to principal moments of inertia of ligand, for efficient grid creation
	 * @param receptor
	 * @param ligand
	 * @return
	 */
	private void preprocess(Molecule3D receptor, Molecule3D ligand) {
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
