package com.actelion.research.chem.phesaflex;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer.PheSASetting;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer.SimilarityMode;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;
import com.actelion.research.chem.conf.BondRotationHelper;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import com.actelion.research.chem.forcefield.mmff.MMFFPositionConstraint;
import com.actelion.research.chem.optimization.MCHelper;
import com.actelion.research.chem.optimization.OptimizerLBFGS;
import com.actelion.research.chem.phesa.MolecularVolume;
import com.actelion.research.chem.phesa.PheSAAlignment;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;


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
	private final StereoMolecule refMol;
	private final StereoMolecule fitMol;
	private final MolecularVolume refVol;
	private final MolecularVolume fitVol;
	private final Map<String, Object> ffOptions;
	private PheSASetting settings;
	
	public PheSASetting getSettings() {
		return settings;
	}

	public void setSettings(PheSASetting settings) {
		this.settings = settings;
	}

	public FlexibleShapeAlignment(StereoMolecule refMol, StereoMolecule fitMol) {
		this(refMol, fitMol, new MolecularVolume(refMol), new MolecularVolume(fitMol));
	}
	
	
	public FlexibleShapeAlignment(StereoMolecule refMol,StereoMolecule fitMol, MolecularVolume refVol, MolecularVolume fitVol) {
		this.refMol = refMol;
		this.fitMol = fitMol;
		this.refVol = refVol;
		this.fitVol = fitVol;
		ffOptions = new HashMap<>();
		ffOptions.put("dielectric constant", 4.0);
		settings = new PheSASetting();
	}
	
	public double[] align() {
		double[] result = new double[3];
		PheSAAlignment shapeAlign = new PheSAAlignment(refVol,fitVol);
		
		double e0 = calcMin(fitMol);
		if(Double.isNaN(e0)) {
			System.err.println("no force field parameters for this structure");
			return result;
		}
		
		restrainedRelaxation(fitMol,e0);

		boolean[] isHydrogen = new boolean[fitMol.getAllAtoms()];
		for(int at=0;at<fitMol.getAllAtoms();at++) {
			
			isHydrogen[at] = fitMol.getAtomicNo(at) == 1;
		}
		Conformer fitConf = new Conformer(fitMol);
		BondRotationHelper  torsionHelper = new BondRotationHelper(fitConf.getMolecule(),true);
		EvaluableFlexibleOverlap eval = new EvaluableFlexibleOverlap(shapeAlign, refMol, fitConf, torsionHelper, 
				settings, isHydrogen, ffOptions);
		eval.setState(eval.getState());
		eval.setE0(e0);
		OptimizerLBFGS opt = new OptimizerLBFGS(200,0.001);
		double t0 = opt.optimize(eval)[0];
		//eval.getState(v);
		Random random = new Random(12345L);
		MCHelper mcHelper = new MCHelper(torsionHelper, null, random);
		//MetropolisMonteCarloHelper mcHelper = new MetropolisMonteCarloHelper(fitMol);
		boolean success = torsionHelper.getRotatableBonds().length!=0;
		double[] v = eval.getState();
		double told = t0;
		if(success) {
			for(int i=0;i<MC_STEPS;i++) {
				double [] vold = Arrays.stream(v).toArray(); // now copy v
				mcHelper.torsionPerturbation(fitConf, v);
				eval.setState(v);
				eval.setE0(e0);
				opt = new OptimizerLBFGS(200,0.001);
				opt.optimize(eval);
				double tnew = getSimilarity(eval,shapeAlign);
				if(!mcHelper.accept(told, tnew)) {
					v = vold;
					eval.setState(v);
				}
				else {
					eval.getState(v);
					told = tnew;
				}
	
			}
		}
		
		Conformer alignedConf = eval.getFitConf();
		for(int a=0;a<alignedConf.getMolecule().getAllAtoms();a++) {
			fitMol.setAtomX(a, alignedConf.getX(a));
			fitMol.setAtomY(a, alignedConf.getY(a));
			fitMol.setAtomZ(a, alignedConf.getZ(a));
		}
		result = getResult();
		return result;
	}
	
	private double getSimilarity(EvaluableFlexibleOverlap eval, PheSAAlignment shapeAlign) { 
		boolean tversky = settings.getSimMode() != SimilarityMode.TANIMOTO;
		double ppWeight = settings.getPpWeight();
		double tverskyCoeff = settings.getSimMode()==SimilarityMode.TVERSKY ? PheSAAlignment.TVERSKY_COEFFICIENT : 1.0-PheSAAlignment.TVERSKY_COEFFICIENT;
		double Obb = eval.getFGValueShapeSelf(new double[3*fitMol.getAllAtoms()], shapeAlign.getMolGauss(),false);
		double Oaa = eval.getFGValueShapeSelf(new double[3*refMol.getAllAtoms()], shapeAlign.getRefMolGauss(),true);
		double Oab = eval.getFGValueShape(new double[3*fitMol.getAllAtoms()]);
		double Tshape = 0.0;
		if(tversky)
			Tshape = Oab/(tverskyCoeff*Obb+(1.0-tverskyCoeff)*Oaa);
		else
			Tshape = Oab/(Oaa+Obb-Oab);
		if(!tversky && Tshape>1.0) //can happen because of weights
			Tshape = 1.0f;	
		double correctionFactor = shapeAlign.getRefMolGauss().getPPGaussians().size()/shapeAlign.getRefMolGauss().getPPGaussians().stream().mapToDouble(g -> g.getWeight()).sum();
		double ObbPP = eval.getFGValueSelfPP(shapeAlign.getMolGauss(),false);
		double OaaPP = eval.getFGValueSelfPP(shapeAlign.getRefMolGauss(),true);
		double OabPP = eval.getFGValuePP();
		double Tpp = 0.0;
		if(shapeAlign.getRefMolGauss().getPPGaussians().isEmpty() && shapeAlign.getMolGauss().getPPGaussians().isEmpty())
			Tpp = 1.0;
		else {
			if(tversky)
				Tpp = OabPP/(tverskyCoeff*ObbPP+(1.0-tverskyCoeff)*OaaPP);
			else
				Tpp = (OabPP/(OaaPP+ObbPP-OabPP));
		}
		Tpp*=correctionFactor;
		if(!tversky && Tshape>1.0) //can happen because of weights
			Tshape = 1.0f;
		if(!tversky && Tpp>1.0) //can happen because of weights
			Tpp = 1.0f;

		return (1.0f-(float)ppWeight)*Tshape + (float)ppWeight*Tpp;
	}
	
	
	
	private double[] getResult() { 
		TransformationSequence sequence = new TransformationSequence();
		PheSAAlignment pa = new PheSAAlignment(fitMol,refMol, settings.getPpWeight());
		double[] r = pa.findAlignment(new double[][] {{0.00, 0.00, 0.00,0.0,0.0,0.0}},sequence,false,settings.getSimMode());
		return new double[] {r[0],r[1],r[2], r[3]};
	}
	
	public double calcMin(StereoMolecule fitMol) {
		try {
			ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
			ForceFieldMMFF94 forceField = new ForceFieldMMFF94(new StereoMolecule(fitMol), ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
			forceField.minimise();
			return forceField.getTotalEnergy();
		} catch (RuntimeException rte) {
			return Double.NaN;
		}
	}
	
	public void restrainedRelaxation(StereoMolecule fitMol, double e0) {
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		double init = 0.2;
		boolean notRelaxed = true;
		int maxCycles = 10;
		int cycles = 0;
		while(notRelaxed && cycles<maxCycles) {
			ForceFieldMMFF94 forceField = new ForceFieldMMFF94(fitMol, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
			MMFFPositionConstraint constraint = new MMFFPositionConstraint(fitMol,50,init);
			forceField.addEnergyTerm(constraint);
			forceField.minimise();
			double e = forceField.getTotalEnergy();
			notRelaxed = (e>e0) && (e-e0>ENERGY_CUTOFF);
			init += 0.2;
			cycles++;
		}
	}
}
