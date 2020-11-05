package com.actelion.research.chem.docking;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;


import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.docking.scoring.AbstractScoringEngine;
import com.actelion.research.chem.optimization.Evaluable;
import com.actelion.research.chem.phesa.PheSAAlignment;
import com.actelion.research.chem.phesa.Quaternion;


public class LigandPose implements Evaluable{
	
	private static double MOVE_AMPLITUDE = 2.0;
	private BondRotationHelper torsionHelper;
	//private double[] torsionValues;
	//private Quaternion rotation; //quaternion
	//private double[] translation;
	private double[] state; //coordinates
	private Conformer ligConf;
	private StereoMolecule mol;
	private AbstractScoringEngine engine;
	private Map<Integer,int[][]> rearAtoms;
	public static long SEED = 12345L;
	
	public LigandPose(Conformer ligConf, AbstractScoringEngine engine, double e0) {
		
		this.engine = engine;
		this.ligConf = ligConf;
		init(e0);
	}
	
	private void init(double e0) {
		mol = ligConf.getMolecule();
		state = new double[3*mol.getAllAtoms()];
		torsionHelper = new BondRotationHelper(mol);
		engine.init(this,e0);
		updateState();
		rearAtoms = new HashMap<Integer,int[][]>();
		assessRearAtoms();
		
	}

	//constrain bonds that are not rotatable, constrain bond lengths and angles
	public double getFGValue(double[] gradient) {
		for(int i=0;i<gradient.length;i++) {
			gradient[i] = 0.0;
		}
		double energy = engine.getFGValue(gradient);
		return energy;
			
	}
	

	
	public void updateState() {
		for(int a=0;a<mol.getAllAtoms();a++) {
			Coordinates c = ligConf.getCoordinates(a);
			state[3*a] = c.x;
			state[3*a+1] = c.y;
			state[3*a+2] = c.z;
		}
		engine.updateState();

	}
	

	@Override
	public void setState(double[] state){
		assert this.state.length==state.length;
		for(int i=0;i<state.length;i++) {
			this.state[i] = state[i];
		}
		for(int a=0;a<mol.getAllAtoms();a++) {
			Coordinates c = new Coordinates(state[3*a],state[3*a+1],state[3*a+2]);
			ligConf.setCoordinates(a, c);
		}
	}

	public double[] getState(double[] v){
		for(int i=0;i<this.state.length;i++) {
			v[i] = state[i];
			
		}
		return v;
	}
	
	public double getGyrationRadius() {
		Coordinates com = DockingUtils.getCOM(ligConf);
		double r = 0.0;
		int counter = 0;
		for(int a=0;a<ligConf.getMolecule().getAtoms();a++) {
			Coordinates c = ligConf.getCoordinates(a);
			r+= c.distanceSquared(com);
			counter++;
		}
		r/=counter;
		return Math.sqrt(r);
	}
	
	public double[] getState() {
		return this.getState(new double[state.length]);
	}
	
	public void randomPerturbation(Random random) {
		int num = (int) (3*random.nextDouble());
		if(num==0) { //translation
			Coordinates shift = DockingUtils.randomVectorInSphere(random).scale(MOVE_AMPLITUDE);
			for(int a=0;a<ligConf.getMolecule().getAllAtoms();a++) { 
				Coordinates c = ligConf.getCoordinates(a);
				c.add(shift);
			}
		}
		else if(num==1) {
			double r = getGyrationRadius();
			Coordinates rot = DockingUtils.randomVectorInSphere(random).scale(MOVE_AMPLITUDE/r);
			double angle = rot.dist();
			Quaternion q = new Quaternion(1.0,0.0,0.0,0.0);
			if(angle>0.0001) {
				Coordinates axis = rot.scale(1.0/angle);
				q = new Quaternion(axis,angle);
			}
			Matrix m = q.getRotMatrix();
			Coordinates com = DockingUtils.getCOM(ligConf);
			ligConf.translate(-com.x, -com.y, -com.z);
			PheSAAlignment.rotateMol(ligConf, m);
			ligConf.translate(com.x, com.y, com.z);

			
		}
		else  {
			if(torsionHelper.getRotatableBonds().length==0)
				return;
			double rnd = random.nextDouble();
			rnd*=180.0;
			double rotateBy = random.nextBoolean() ? rnd : -rnd;
			rotateBy = rotateBy*Math.PI/180.0;
			int bond = random.nextInt(torsionHelper.getRotatableBonds().length);
			bond = torsionHelper.getRotatableBonds()[bond];
			
			torsionHelper.rotateAroundBond(bond, rotateBy,ligConf,random.nextBoolean());

			
		}
		updateState();
		
	}
	
	public Conformer getLigConf() {
		return ligConf;
	}
	
	private void assessRearAtoms() {
		StereoMolecule mol = ligConf.getMolecule();
		for(int b=0;b<torsionHelper.getRotatableBonds().length;b++) {
			int[][] rAtoms = new int[2][];
			int[] atoms = torsionHelper.getTorsionAtoms()[b];
			rearAtoms.put(b,rAtoms);
			if(atoms[0]==-1) {
				rAtoms[0] = new int[2];
				int index = 0;
				for(int i=0;i<mol.getConnAtoms(atoms[1]);i++) {
					if(mol.getConnAtom(atoms[1], i)==atoms[2])
						continue;
					else {
						rAtoms[0][index] = mol.getConnAtom(atoms[1], i);
						index++;
					}
				}
			}
			
			if(atoms[3]==-1) {
				rAtoms[1] = new int[2];
				int index = 0;
				for(int i=0;i<mol.getConnAtoms(atoms[2]);i++) {
					if(mol.getConnAtom(atoms[2], i)==atoms[1])
						continue;
					else {
						rAtoms[1][index] = mol.getConnAtom(atoms[2], i);
						index++;
					}
				}
			}
			}
			
		}
		
}
	
	
	


