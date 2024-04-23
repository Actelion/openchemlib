package com.actelion.research.chem.docking;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;


import com.actelion.research.calc.Matrix;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.transformation.ExponentialMap;
import com.actelion.research.chem.alignment3d.transformation.Quaternion;
import com.actelion.research.chem.alignment3d.transformation.Rotation;
import com.actelion.research.chem.alignment3d.transformation.RotationDerivatives;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;
import com.actelion.research.chem.alignment3d.transformation.Translation;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.TorsionDB;
import com.actelion.research.chem.docking.scoring.AbstractScoringEngine;
import com.actelion.research.chem.optimization.Evaluable;
import com.actelion.research.chem.optimization.MCHelper;
import com.actelion.research.chem.potentialenergy.PositionConstraint;
import com.actelion.research.chem.potentialenergy.PotentialEnergyTerm;


public class LigandPose implements Evaluable{
	
	private static double MOVE_AMPLITUDE = 2.0;
	private BondRotationHelper torsionHelper;
	//private Coordinates rotationCenter;
	private double[] state; //internal coordinates: translation,rotation,dihedral
	private Conformer ligConf;
	private Coordinates[] origCoords;
	private Coordinates[] cachedCoords; //used for gradient calculation: ligand coordinates with adjusted dihedral angles, but before rotation and translation
	private double[][] dRdvi1;
	private double[][] dRdvi2;
	private double[][] dRdvi3;
	private StereoMolecule mol;
	private AbstractScoringEngine engine;
	public static long SEED = 12345L;
	private Coordinates origCOM;
	private MCHelper mcHelper;
	//private int[] mcsRotBondIndeces; //for MCS docking, only bonds not part of the MCS are sampled
	
	
	public LigandPose(Conformer ligConf, AbstractScoringEngine engine, double e0) {
		
		this.engine = engine;
		this.ligConf = ligConf;
		init(e0);
	}
	
	/**
	 * for MCS docking:create array of bond indices that are allowed to be permuted
	 * @param constraints
	 */
	public void setMCSBondConstraints(List<Integer> constraints) {
		int[] rotBonds = torsionHelper.getRotatableBonds();
		List<Integer> allowedIndices = new ArrayList<>();
		for(int rbIndex=0;rbIndex<rotBonds.length;rbIndex++) {
			int rb = rotBonds[rbIndex];
			if(!constraints.contains(rb))
				allowedIndices.add(rbIndex);
		}
		int [] mcsRotBondIndeces = new int[allowedIndices.size()];
		for(int i=0;i<allowedIndices.size();i++) {
			mcsRotBondIndeces[i]=allowedIndices.get(i);
		}
		mcHelper.setMcsRotBondIndeces(mcsRotBondIndeces);
		
	}
	
	private void init(double e0) {
		mol = ligConf.getMolecule();
		torsionHelper = new BondRotationHelper(mol,true);
		engine.init(this,e0);
		setInitialState();
		origCoords = new Coordinates[ligConf.getMolecule().getAllAtoms()];
		cachedCoords = new Coordinates[ligConf.getMolecule().getAllAtoms()];
		origCOM = new Coordinates();
		for(int a=0;a<ligConf.getMolecule().getAllAtoms();a++) {
			origCoords[a] = new Coordinates(ligConf.getCoordinates(a));
			cachedCoords[a] = new Coordinates(ligConf.getCoordinates(a));
			origCOM.add(cachedCoords[a]);
		}
		origCOM.scale(1.0/cachedCoords.length);
		dRdvi1 = new double[3][3];
		dRdvi2 = new double[3][3];
		dRdvi3 = new double[3][3];
		mcHelper = new MCHelper(torsionHelper,null,new Random(SEED));
		
	}
	
	private void resetLigCoordinates() {
		for(int a=0;a<ligConf.getMolecule().getAllAtoms();a++) {
			ligConf.setX(a, origCoords[a].x);
			ligConf.setY(a, origCoords[a].y);
			ligConf.setZ(a, origCoords[a].z);
		}
	}
	/**
	 * Fuhrmann J, Rurainski A, Lenhof HP, Neumann D. A new method for the gradient-based optimization of molecular complexes. 
	 * J Comput Chem. 2009 Jul 15;30(9):1371-8. doi: 10.1002/jcc.21159. PMID: 19031415.
	 */
	//constrain bonds that are not rotatable, constrain bond lengths and angles
	public double getFGValue(double[] gradient) {
		double[] coordGrad = new double[ligConf.getMolecule().getAllAtoms()*3];
		for(int i=0;i<gradient.length;i++) {
			gradient[i] = 0.0;
		}
		double energy = engine.getFGValue(coordGrad);
		//to inner coordinates
		//1. with respect to translational DOG
		for(int a=0;a<ligConf.getMolecule().getAllAtoms();a++) {
			gradient[0] += coordGrad[3*a]; 
			gradient[1] += coordGrad[3*a+1]; 
			gradient[2] += coordGrad[3*a+2]; 
		}
		//2. orientational 
		//with respect to vector of exponential mapping p
		// dE/dpj = Tj*vi'*dE/dx
		//vi': atomic position (after adjustment of torsion values)
		double[] p = new double[] {state[3],state[4],state[5]};
		RotationDerivatives transformDerivatives = new RotationDerivatives(p);
		transformDerivatives.dRdv(0, dRdvi1);
		transformDerivatives.dRdv(1, dRdvi2);
		transformDerivatives.dRdv(2, dRdvi3);
		for(int a=0;a<ligConf.getMolecule().getAllAtoms();a++) {
			Coordinates vi = cachedCoords[a];
			Coordinates Tj_vi = vi.rotateC(dRdvi1);
			gradient[3] += coordGrad[3*a]*Tj_vi.x+coordGrad[3*a+1]*Tj_vi.y+coordGrad[3*a+2]*Tj_vi.z;
			Tj_vi = vi.rotateC(dRdvi2);
			gradient[4] += coordGrad[3*a]*Tj_vi.x+coordGrad[3*a+1]*Tj_vi.y+coordGrad[3*a+2]*Tj_vi.z;
			Tj_vi = vi.rotateC(dRdvi3);
			gradient[5] += coordGrad[3*a]*Tj_vi.x+coordGrad[3*a+1]*Tj_vi.y+coordGrad[3*a+2]*Tj_vi.z;
		}
		//3. torsional gradient
		for(int b=0;b<torsionHelper.getRotatableBonds().length;b++) {
			int[] rotatedAtoms = torsionHelper.getSmallerSideAtomLists()[b];
			int j = torsionHelper.getRotationCenters()[b];
			int k = torsionHelper.getTorsionAtoms()[b][1] == j ? torsionHelper.getTorsionAtoms()[b][2] : torsionHelper.getTorsionAtoms()[b][1];
			Coordinates v1 = ligConf.getCoordinates(k).subC(ligConf.getCoordinates(j));

			for(int i : rotatedAtoms) {
				Coordinates v2 = 
						ligConf.getCoordinates(i).subC(ligConf.getCoordinates(j));
				Coordinates dx_dphi = v1.cross(v2);
				gradient[6+b] += dx_dphi.x*coordGrad[3*i] + dx_dphi.y*coordGrad[3*i+1] + 
						dx_dphi.z*coordGrad[3*i+2];
			}
			
		}
		
		return energy;
			
	}
	
	public Map<String,Double> getContributions() {
		return engine.getContributions();
	}
	
	public void setInitialState() {
		int elements = 3+3+torsionHelper.getRotatableBonds().length; //3 translational, 3 rotational, 3 torsion
		state = new double[elements];
		state[0] = 0.0;
		state[1] = 0.0;
		state[2] = 0.0;
		Quaternion quat = new Quaternion(1.0,0.0,0.0,0.0);
		ExponentialMap emap = new ExponentialMap(quat);
		state[3] = emap.getP().x;
		state[4] = emap.getP().y;
		state[5] = emap.getP().z;
		for(int b=0;b<torsionHelper.getRotatableBonds().length;b++) {
			int[] atoms = torsionHelper.getTorsionAtoms()[b];
			state[6+b] = TorsionDB.calculateTorsionExtended(ligConf, atoms);
		}
	}
	
	public void updateLigandCoordinates() {
		resetLigCoordinates();
		//1. update dihedral angles
		//2. translate COM of ligand to rotation center
		//3. apply rotation 
		//4. translate back
		//5. apply translation 
		updateDihedralAngles();
		for(int a=0;a<ligConf.getMolecule().getAllAtoms();a++) {
			cachedCoords[a] = new Coordinates(ligConf.getCoordinates(a));
		}
		ExponentialMap eMap = new ExponentialMap(state[3],state[4],state[5]);
		Quaternion q = eMap.toQuaternion();
		Translation trans1 = new Translation(origCOM.scaleC(-1.0));
		Translation trans2 = new Translation(origCOM);
		Rotation rot = new Rotation(q.getRotMatrix().getArray());
		Translation t = new Translation(state[0],state[1],state[2]);
		TransformationSequence transformation = new TransformationSequence();
		transformation.addTransformation(trans1);
		transformation.addTransformation(rot);
		transformation.addTransformation(trans2);
		transformation.addTransformation(t);
		transformation.apply(ligConf);
	}
	
	private void updateDihedralAngles() {
		for(int b=0;b<torsionHelper.getRotatableBonds().length;b++) {
			double targetTorsion = state[6+b];
			int[] atoms = torsionHelper.getTorsionAtoms()[b];
			double currentTorsion = TorsionDB.calculateTorsionExtended(ligConf, atoms);
			double deltaTorsion = targetTorsion - currentTorsion;
			torsionHelper.rotateAroundBond(b, deltaTorsion,ligConf,false);
		}
	}


	@Override
	public void setState(double[] state){
		assert this.state.length==state.length;
		for(int i=0;i<state.length;i++) {
			if(i>5) { //torsions
				if(state[i]>Math.PI) {
					state[i] -= 2*Math.PI;
				}
			}
			this.state[i] = state[i];
	
		}
		updateLigandCoordinates();
	}

	public double[] getState(double[] v){
		for(int i=0;i<this.state.length;i++) {
			v[i] = state[i];
			
		}
		return v;
	}
	
	public double[] getCartState(){
		double[] cartState = new double[3*ligConf.getMolecule().getAllAtoms()];
		for(int a=0;a<ligConf.getMolecule().getAllAtoms();a++) {
			cartState[3*a] = ligConf.getCoordinates(a).x;
			cartState[3*a+1] = ligConf.getCoordinates(a).y;
			cartState[3*a+2] = ligConf.getCoordinates(a).z;
			
		}
		return cartState;
	}

	
	public double getScore() {
		return this.engine.getScore();
	}
	
	public double[] getState() {
		return this.getState(new double[state.length]);
	}
	
	public void randomPerturbation() {
		mcHelper.randomPerturbation(ligConf, state);
		updateLigandCoordinates();
		
	}
	
	public void addPositionalConstraints(double d) {
		for(int a=0;a<ligConf.getMolecule().getAtoms();a++) {
			PositionConstraint constraint = new PositionConstraint(ligConf,a,50.0, d);
			engine.addConstraint(constraint);
		}
	}
	
	public void addConstraint(PositionConstraint constraint) {
		engine.addConstraint(constraint);
	}
	
	public void removeConstraints() {
		engine.removeConstraints();
	}
	
	
	public Conformer getLigConf() {
		return ligConf;
	}
	

	

}
	
	
	


