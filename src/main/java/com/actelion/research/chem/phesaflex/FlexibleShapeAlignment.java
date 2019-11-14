package com.actelion.research.chem.phesaflex;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.forcefield.mmff.PositionConstraint;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.OptimizerLBFGS;
import com.actelion.research.chem.phesa.PheSAAlignment;


/**
 * Performs flexible Alignment of two Molecules that are prealigned. A fit molecule is thereby
 * aligned to a rigid template. The internal degrees of freedom of the fit molecule are flexible
 * in order to not create molecular conformations with too high strain energy, the force field potential
 * energy is part of the objective function of the alignment optimization. The procedure is inspired by 
 * doi: 10.1021/acs.jcim.7b00618, but uses analytical gradients and optimizes the directionality overlap of
 * pharmacophore features. Furthermore, in addition to local optimization, Monte Carlo steps are performed to randomly
 * perturb dihedral angles of the base molecule.
 * The first step is to perform constrained minimizations of the rigidly prealigned fit molecule. The well-depth of the positional
 * constraint is gradually increased until to strain energy is lower than 10 kcal/mol.
 * This yields a good initial guess of a energetically low-lying conformation of the fit molecule that still
 * has a reasonable PheSA overlap. From this first conformation, a Monte Carlo optimization with local minimization is started. 
 * The objective function is:
 * O = -Tphesa + lambda*Estrain*Estrain 
 * Estrain = E-E0, E0=10 kcal/mol
 * the strain penalty is only applied if the strain energy is higher than 10 kcal/mol, otherwise the 
 * objective function is solely depending on the PheSA similarity Tphesa
 * @author JW
 *
 */


public class FlexibleShapeAlignment {
	private static final int MC_STEPS = 50;
	public static final double ENERGY_CUTOFF = 10.0;
	private StereoMolecule refMol;
	private StereoMolecule fitMol;
	private MolecularVolume refVol;
	private MolecularVolume fitVol;
	Map<String, Object> ffOptions;
	
	public FlexibleShapeAlignment(StereoMolecule refMol, StereoMolecule fitMol) {
		this(refMol, fitMol, new MolecularVolume(refMol), new MolecularVolume(fitMol));
	}
	
	public FlexibleShapeAlignment(StereoMolecule refMol,StereoMolecule fitMol, MolecularVolume refVol, MolecularVolume fitVol) {
		this.refMol = refMol;
		this.fitMol = fitMol;
		this.refVol = refVol;
		this.fitVol = fitVol;
		ffOptions = new HashMap<String, Object>();
		ffOptions.put("dielectric constant", 4.0);
	}
	
	public double align() {
		
		PheSAAlignment shapeAlign = new PheSAAlignment(refVol,fitVol);
		
		double e0 = calcMin(fitMol);
		if(Double.isNaN(e0)) {
			System.err.print("no force field parameters for this structure");
			return 0.0;
		}
		
		restrainedRelaxation(fitMol,e0);
		
		double[] v = new double[3*fitMol.getAllAtoms()];

		boolean[] isHydrogen = new boolean[fitMol.getAllAtoms()];
		for(int at=0;at<fitMol.getAllAtoms();at++) {
			
			isHydrogen[at] = fitMol.getAtomicNo(at)==1 ? true : false;
		}
		EvaluableFlexibleOverlap eval = new EvaluableFlexibleOverlap(shapeAlign, refMol, fitMol, isHydrogen, v, ffOptions);
		eval.setE0(e0);
		OptimizerLBFGS opt = new OptimizerLBFGS(200,0.001);
		opt.optimize(eval);
		eval.getState(v);
		double t0 = getTanimoto(eval,shapeAlign);

		MetropolisMonteCarloHelper mcHelper = new MetropolisMonteCarloHelper(fitMol);
		double told = t0;
		for(int i=0;i<MC_STEPS;i++) {
			double [] vold = Arrays.stream(v).toArray(); // now copy v
			mcHelper.step();
			eval = new EvaluableFlexibleOverlap(shapeAlign, refMol, fitMol, isHydrogen, v, ffOptions);
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
			ForceFieldMMFF94 forceField = new ForceFieldMMFF94(fitMol, ForceFieldMMFF94.MMFF94SPLUS,ffOptions);
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
		ForceFieldMMFF94 forceField = new ForceFieldMMFF94(new StereoMolecule(fitMol), ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
		forceField.minimise();
		double e0 = forceField.getTotalEnergy();
		return e0;
		
		
	}
	
	public void restrainedRelaxation(StereoMolecule fitMol, double e0) {
		
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		double init = 0.2;
		boolean notRelaxed = true;
		while(notRelaxed) {
			ForceFieldMMFF94 forceField = new ForceFieldMMFF94(fitMol, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
			PositionConstraint constraint = new PositionConstraint(fitMol,50,init);
			forceField.addEnergyTerm(constraint);
			forceField.minimise();
			double e = forceField.getTotalEnergy();
			notRelaxed = e-e0>ENERGY_CUTOFF;
			init += 0.2;
		}
		
	}
	
	
	

}
