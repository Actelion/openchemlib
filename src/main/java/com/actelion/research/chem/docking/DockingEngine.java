package com.actelion.research.chem.docking;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.optimization.OptimizerLBFGS;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.PheSAAlignment;
import com.actelion.research.chem.phesa.PheSAMolecule;

public class DockingEngine {
	
	private static final int DEFAULT_NR_MC_STEPS = 200;
	private static final int DEFAULT_START_POSITIONS = 5;
	private static final double BOLTZMANN_FACTOR = 0.593;
	
	private StereoMolecule receptor;
	private MoleculeGrid grid;
	private PheSAMolecule refVol;
	private Matrix rotation; //for initial prealignment to native ligand
	private Coordinates origCOM;
	private int mcSteps;
	private int startPositions;
	private Random random;
	
	public DockingEngine(StereoMolecule receptor, StereoMolecule nativeLigand, int mcSteps, int startPositions) {
		grid = new MoleculeGrid(nativeLigand, 0.5, new Coordinates(8.0,8.0,8.0),false);
		MolecularVolume molVol = new MolecularVolume(nativeLigand);
		origCOM  = new Coordinates(molVol.getCOM());
		Conformer conf = new Conformer(nativeLigand);
		rotation = PheSAAlignment.preProcess(conf, molVol).getTranspose();
		refVol = new PheSAMolecule(nativeLigand,molVol);
		this.mcSteps = mcSteps;
		this.startPositions = startPositions;
		this.random = new Random(12354L);
		this.receptor = receptor;
	}
	
	public DockingEngine(StereoMolecule receptor, StereoMolecule nativeLigand) {
		this(receptor,nativeLigand,DEFAULT_NR_MC_STEPS,DEFAULT_START_POSITIONS);
	}
	
	
	
	public StereoMolecule dockMolecule(StereoMolecule mol) {
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
		Conformer ligConf = new Conformer(mol);
		LigandPose pose = initiate(ligConf);
		mcSearch(pose);
		return ligConf.toMolecule();
		
	}
	
	private void mcSearch(LigandPose pose) {
		System.out.println("start docking");
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
			System.out.println(i);
			pose.randomPerturbation();
			state = optimizer.optimize(pose);
			energy = pose.getFGValue(new double[bestState.length]);
			System.out.println(energy);
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
		
	}
	
	private LigandPose initiate(Conformer ligConf) {
		int[] gridSize = grid.getGridSize();
		Conformer recConf = new Conformer(receptor);
		Set<Integer> receptorAtoms = new HashSet<Integer>();
		int[] receptorAtomTypes = new int[receptor.getAtoms()];
		for(int i=0;i<receptor.getAtoms();i++) {
			receptorAtomTypes[i] = InteractionAtomTypeCalculator.getAtomType(receptor, i);
			int[] gridC = grid.getGridCoordinates(receptor.getCoordinates(i));
			int x = gridC[0];
			int y = gridC[1];
			int z = gridC[2];	
			if(x>0 && x<gridSize[0]) {
				if(y>0 && y<gridSize[1]) {
					if(z>0 && z<gridSize[2]) {
						receptorAtoms.add(i);
					}
			}
		}
		}
		LigandPose pose = new LigandPose(ligConf,recConf,receptorAtoms,receptorAtomTypes, grid);
		
		return pose;
		
		
	}

	
	
	
	

}
