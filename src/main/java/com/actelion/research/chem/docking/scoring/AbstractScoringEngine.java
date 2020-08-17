package com.actelion.research.chem.docking.scoring;

import java.util.Set;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule3D;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.chem.optimization.Evaluable;


public abstract class AbstractScoringEngine implements Evaluable  {
	
	private double BUMP_PENALTY = 500;
	private int BUMP_RADIUS = 3;
	
	protected Conformer receptorConf;
	protected Set<Integer> bindingSiteAtoms;
	protected Conformer candidatePose;
	protected double[] state;
	protected MoleculeGrid grid;
	
	public AbstractScoringEngine(StereoMolecule receptor, Set<Integer> bindingSiteAtoms, MoleculeGrid grid) {
		this.receptorConf = new Conformer(receptor);
		this.bindingSiteAtoms = bindingSiteAtoms;
		this.grid = grid;
	}
	
	public Conformer getCandidatePose() {
		return candidatePose;
	}
	
	public void updateState() {
		for(int a=0;a<candidatePose.getMolecule().getAllAtoms();a++) {
			Coordinates c = candidatePose.getCoordinates(a);
			state[3*a] = c.x;
			state[3*a+1] = c.y;
			state[3*a+2] = c.z;
		}
	}
	

	@Override
	public void setState(double[] state){
		assert this.state.length==state.length;
		for(int i=0;i<state.length;i++) {
			this.state[i] = state[i];
		}
		for(int a=0;a<candidatePose.getMolecule().getAllAtoms();a++) {
			Coordinates c = new Coordinates(state[3*a],state[3*a+1],state[3*a+2]);
			candidatePose.setCoordinates(a, c);
		}
	}
	
	@Override
	public double[] getState() {
		return this.getState(new double[state.length]);
	}
	
	public double[] getState(double[] v){
		for(int i=0;i<this.state.length;i++) {
			v[i] = state[i];
			
		}
		return v;
	}
	
	public double getBumpTerm() {
		double bumpTerm = 0.0;
		int[] gridSize = grid.getGridSize();
		for(Coordinates c : candidatePose.getCoordinates()) {
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

	public abstract void init(Conformer candidatePose);

	

}
