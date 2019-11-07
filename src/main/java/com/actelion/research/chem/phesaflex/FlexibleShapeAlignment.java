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
 * energy is part of the objective function of the alignment optimization
 * @author joelwahl
 *
 */


public class FlexibleShapeAlignment {
	private StereoMolecule templateMol;
	private StereoMolecule fitMol;
	private StereoMolecule alignedMol;

	
	public FlexibleShapeAlignment(StereoMolecule templateMol, StereoMolecule fitMol) {
		this.templateMol = templateMol;
		this.fitMol = fitMol;
		this.alignedMol = new StereoMolecule();
	}
	
	public double align() {
		

		MolecularVolume templateVol = new MolecularVolume(templateMol); 
		MolecularVolume fitVol = new MolecularVolume(fitMol);
		
		PheSAAlignment shapeAlign = new PheSAAlignment(templateVol,fitVol);
		
		double e0 = calcMin(fitMol);
		
		double[] v = new double[3*fitMol.getAllAtoms()];

		boolean[] isHydrogen = new boolean[fitMol.getAllAtoms()];
		for(int at=0;at<fitMol.getAllAtoms();at++) {
			
			isHydrogen[at] = fitMol.getAtomicNo(at)==1 ? true : false;
		}
		EvaluableFlexibleOverlap eval = new EvaluableFlexibleOverlap(shapeAlign, fitMol, isHydrogen, v);
		eval.setE0(e0);
		OptimizerLBFGS opt = new OptimizerLBFGS(200,0.001);
		eval.setState(v);
		opt.optimize(eval);
		eval.getState(v);
		double Obb = eval.getFGValueShapeSelf(new double[v.length], shapeAlign.getMolGauss(),false);
		double Oaa = eval.getFGValueShapeSelf(new double[v.length], shapeAlign.getRefMolGauss(),true);
		double Oab = eval.getFGValueShape(new double[v.length]);

		for(int a=0,i=0;i<fitMol.getAllAtoms();i++) {
			fitMol.setAtomX(i,v[a++]);
			fitMol.setAtomY(i,v[a++]);
			fitMol.setAtomZ(i,v[a++]);
		}
		
		return (Oab/(Oaa+Obb-Oab));
		
		
		
		
	
		
			
	}
	
	public double calcMin(StereoMolecule fitMol) {
		
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		ForceFieldMMFF94 forceField = new ForceFieldMMFF94(new StereoMolecule(fitMol), ForceFieldMMFF94.MMFF94SPLUS);
		forceField.minimise();
		double e0 = forceField.getTotalEnergy();
		return e0;
		
		
	}
	
	public StereoMolecule getAlignedMol() {
		return this.alignedMol;
	}

}
