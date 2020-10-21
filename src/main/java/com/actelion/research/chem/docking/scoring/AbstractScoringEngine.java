package com.actelion.research.chem.docking.scoring;

import java.util.Set;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.docking.LigandPose;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.optimization.Evaluable;


public abstract class AbstractScoringEngine  {
	
	private double BUMP_PENALTY = 500;
	private int BUMP_RADIUS = 3;
	
	protected Conformer receptorConf;
	protected Set<Integer> bindingSiteAtoms;
	protected LigandPose candidatePose;
	protected MoleculeGrid grid;
	
	public AbstractScoringEngine(StereoMolecule receptor, Set<Integer> bindingSiteAtoms, MoleculeGrid grid) {
		this.receptorConf = new Conformer(receptor);
		this.bindingSiteAtoms = bindingSiteAtoms;
		this.grid = grid;
	}
	
	public LigandPose getCandidatePose() {
		return candidatePose; 
	}
	

	
	public double getBumpTerm() {
		double bumpTerm = 0.0;
		int[] gridSize = grid.getGridSize();
		for(int a=0;a<candidatePose.getLigConf().getMolecule().getAllAtoms();a++) {
			Coordinates c = candidatePose.getLigConf().getCoordinates(a);
			int[] gridC = grid.getGridCoordinates(c);
			int x = gridC[0];
			int y = gridC[1];
			int z = gridC[2];	
			if(x<BUMP_RADIUS || x>(gridSize[0]-BUMP_RADIUS)) {
				bumpTerm = BUMP_PENALTY;
				break;
			}
			else if(y<BUMP_RADIUS || y>(gridSize[1]-BUMP_RADIUS)) {
				bumpTerm = BUMP_PENALTY;
				break;
			}
			else if(z<BUMP_RADIUS || z>(gridSize[2]-BUMP_RADIUS)) {
				bumpTerm = BUMP_PENALTY;
				break;	
				}
			}
		return bumpTerm;
		}

	public abstract void init(LigandPose candidatePose, double e0);
	
	public abstract void updateState();

	public abstract double getFGValue(double[] grad);
	
	public abstract double getScore();
	

}
