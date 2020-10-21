package com.actelion.research.chem.docking;


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
import com.actelion.research.chem.conf.ConformerSet;
import com.actelion.research.chem.conf.ConformerSetGenerator;
import com.actelion.research.chem.docking.scoring.AbstractScoringEngine;
import com.actelion.research.chem.docking.scoring.ChemPLP;
import com.actelion.research.chem.docking.scoring.IdoScore;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.optimization.OptimizerLBFGS;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.PheSAAlignment;

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
		this(receptor,nativeLigand,DEFAULT_NR_MC_STEPS,DEFAULT_START_POSITIONS,ScoringFunction.CHEMPLP);
	}
	
	
	
	public StereoMolecule dockMolecule(StereoMolecule mol) throws DockingFailedException {

		Conformer bestPose = null;
		double bestEnergy = Double.MAX_VALUE;
		
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		ffOptions.put("dielectric constant", 80.0);
		
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);

		ConformerSetGenerator confSetGen = new ConformerSetGenerator(100,ConformerGenerator.STRATEGY_LIKELY_RANDOM, false,
				LigandPose.SEED);
		ConformerSet confSet = confSetGen.generateConformerSet(mol);
		ConformerSet initialPos = new ConformerSet();
		double eMin = Double.MAX_VALUE;
		for(Conformer conformer : confSet) {
			if(conformer!=null) {
				StereoMolecule conf = conformer.toMolecule(new StereoMolecule(mol));
				conf.ensureHelperArrays(Molecule.cHelperParities);
				ForceFieldMMFF94 mmff = new ForceFieldMMFF94(conf, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
				mmff.minimise();
				Conformer ligConf = new Conformer(conf);
				initialPos.add(ligConf);
				if(initialPos.size()>=startPositions)
					break;
				double e = mmff.getTotalEnergy();
				if(e<eMin)
					eMin = e;
			}
		}

		for(Conformer ligConf : initialPos) {
			PheSAAlignment.preProcess(ligConf, new MolecularVolume(ligConf.toMolecule()));
			for(double[] transform : PheSAAlignment.initialTransform(1)) {
				Conformer newLigConf = new Conformer(ligConf);
				PheSAAlignment.rotateMol(newLigConf, transform);
				
				LigandPose pose = initiate(newLigConf,eMin);
				double energy = mcSearch(pose);
				if(energy<bestEnergy) {
					bestEnergy = energy;
					bestPose = pose.getLigConf();
					}
				}
			}
		if(bestPose!=null) {
			StereoMolecule best = bestPose.toMolecule();
			double[][] rot = rotation.getTranspose().getArray();
			PheSAAlignment.rotateMol(best, rot);
			PheSAAlignment.translateMol(best, new double[] {origCOM.x, origCOM.y, origCOM.z} );
			System.out.println("bestEnergy");
			System.out.println(bestEnergy);
			return best;
		}
		else {
			throw new DockingFailedException("docking failed");
		}
		
	}
	
	private double mcSearch(LigandPose pose) {
		double[] bestState = new double[pose.getState().length];
		double[] oldState = new double[pose.getState().length];
		double[] state = new double[pose.getState().length];
		double bestEnergy = -Float.MAX_VALUE;
		double oldEnergy = -Float.MAX_VALUE;
		double energy = -Float.MAX_VALUE;
		OptimizerLBFGS optimizer = new OptimizerLBFGS(200,0.001);
		//oldState = pose.getState();
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
	
	private LigandPose initiate(Conformer ligConf, double e0) {
		LigandPose pose = new LigandPose(ligConf, engine, e0);
		
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
