package com.actelion.research.chem.phesaflex;

import java.util.Arrays;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.OptimizerLBFGS;
import com.actelion.research.chem.phesa.PheSAAlignment;


/**
 * Performs flexible Alignment of two Molecules that are prealigned. A base molecule is thereby
 * aligned to a rigid template. The internal degrees of freedom of the base molecule are flexible
 * in order to not create molecular conformations with too high strain energy, the force field potential
 * energy is part of the objective function of the alignment optimization. The procedure is inspired by 
 * doi: 10.1021/acs.jcim.7b00618, but uses analytical gradients and optimizes the directionality overlap of
 * pharmacophore features. Furthermore, in addition to local optimization, Monte Carlo steps are performed to randomly
 * perturb dihedral angles of the base molecule.
 * @author JW
 *
 */


public class FlexibleShapeAlignment {
	private static final int MC_STEPS = 50;
	private StereoMolecule refMol;
	private StereoMolecule fitMol;
	private MolecularVolume refVol;
	private MolecularVolume fitVol;

	
	public FlexibleShapeAlignment(StereoMolecule refMol, StereoMolecule fitMol) {
		this(refMol, fitMol, new MolecularVolume(refMol), new MolecularVolume(fitMol));
	}
	
	public FlexibleShapeAlignment(StereoMolecule refMol,StereoMolecule fitMol, MolecularVolume refVol, MolecularVolume fitVol) {
		this.refMol = refMol;
		this.fitMol = fitMol;
		this.refVol = refVol;
		this.fitVol = fitVol;
	}
	
	public double align() {
		
		PheSAAlignment shapeAlign = new PheSAAlignment(refVol,fitVol);
		
		double e0 = calcMin(fitMol);
		
		double[] v = new double[3*fitMol.getAllAtoms()];

		boolean[] isHydrogen = new boolean[fitMol.getAllAtoms()];
		for(int at=0;at<fitMol.getAllAtoms();at++) {
			
			isHydrogen[at] = fitMol.getAtomicNo(at)==1 ? true : false;
		}
		EvaluableFlexibleOverlap eval = new EvaluableFlexibleOverlap(shapeAlign, refMol, fitMol, isHydrogen, v);
		eval.setE0(e0);
		OptimizerLBFGS opt = new OptimizerLBFGS(200,0.001);
		opt.optimize(eval);
		eval.getState(v);
		double t0 = getTanimoto(eval,shapeAlign);

		MetropolisMonteCarloHelper mcHelper = new MetropolisMonteCarloHelper(fitMol);
		double told = t0;
		for(int i=0;i<MC_STEPS;i++) {
			double [] vold = Arrays.stream(v).toArray();
			// now copy v
			mcHelper.step();
			eval = new EvaluableFlexibleOverlap(shapeAlign, refMol, fitMol, isHydrogen, v);
			eval.setE0(e0);
			opt = new OptimizerLBFGS(200,0.001);
			opt.optimize(eval);
			double tnew = getTanimoto(eval,shapeAlign);
			if(!mcHelper.accept(told, tnew)) {
				v = vold;
				eval.setState(v);
			}
			else {
				eval.getState(v);
				told = tnew;
			}
		}
		return getTanimoto(eval,shapeAlign);
	}
	
	private double getTanimoto(EvaluableFlexibleOverlap eval, PheSAAlignment shapeAlign) { 
		double Obb = eval.getFGValueShapeSelf(new double[3*fitMol.getAllAtoms()], shapeAlign.getMolGauss(),false);
		double Oaa = eval.getFGValueShapeSelf(new double[3*fitMol.getAllAtoms()], shapeAlign.getRefMolGauss(),true);
		double Oab = eval.getFGValueShape(new double[3*fitMol.getAllAtoms()]);
		return (Oab/(Oaa+Obb-Oab));
	}
	
	public double calcMin(StereoMolecule fitMol) {
		
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		ForceFieldMMFF94 forceField = new ForceFieldMMFF94(new StereoMolecule(fitMol), ForceFieldMMFF94.MMFF94SPLUS);
		forceField.minimise();
		double e0 = forceField.getTotalEnergy();
		return e0;
		
		
	}
	

}
