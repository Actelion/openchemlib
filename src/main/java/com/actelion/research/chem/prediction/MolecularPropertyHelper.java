/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 *
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.chem.prediction;

import com.actelion.research.chem.AtomFunctionAnalyzer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.MolecularFlexibilityCalculator;

public class MolecularPropertyHelper {
	public static final int MOLECULAR_PROPERTY_MOLWEIGHT = 0;
	public static final int MOLECULAR_PROPERTY_CLOGP = 1;
	public static final int MOLECULAR_PROPERTY_CLOGS = 2;
	public static final int MOLECULAR_PROPERTY_TPSA = 3;
	public static final int MOLECULAR_PROPERTY_HDONORS = 4;
	public static final int MOLECULAR_PROPERTY_HACCEPTORS = 5;
	public static final int MOLECULAR_PROPERTY_FLEXIBILITY = 6;
	public static final int MOLECULAR_PROPERTY_COMPLEXITY = 7;
	public static final int MOLECULAR_PROPERTY_SHAPE = 8;
	public static final int MOLECULAR_PROPERTY_ROTATABLEBONDS = 9;
	public static final int MOLECULAR_PROPERTY_STEREOCENTERS = 10;
	public static final int MOLECULAR_PROPERTY_SMALLRINGCOUNT = 11;
	public static final int MOLECULAR_PROPERTY_AROMRINGCOUNT = 12;
	public static final int MOLECULAR_PROPERTY_BASIC_NITROGENS = 13;
	public static final int MOLECULAR_PROPERTY_ACIDIC_OXYGENS = 14;

	public static final PropertySpecification[] SPEC = {
			new PropertySpecification("molweight", "Molecular weight", "", "400", 50f, 0f, 800f),
			new PropertySpecification("cLogP", "cLogP", "", "4", 0.5f, 0f, 8f),
			new PropertySpecification("cLogS", "cLogS", "-4", "", 0.5f, -8f, 2f),
			new PropertySpecification("tpsa", "Polar surface area", "", "120", 20f, 0f, 250f),
			new PropertySpecification("donors", "H-Donors", "", "5", 0.5f, 0, 8f),
			new PropertySpecification("acceptors", "H-Acceptors", "", "10", 1.0f, 0, 16f),
			new PropertySpecification("flexibility", "Molecular flexibility", "0.3", "0.7", 0.05f, 0f, 1f),
			new PropertySpecification("complexity", "Molecular complexity", "0.8", "", 0.05f, 0f, 1f),
			new PropertySpecification("shape", "Molecular shape", "", "0.5", 0.05f, 0f, 1f),
			new PropertySpecification("rotatableBonds", "Rotatable bond count", "", "4", 5f, 0f, 20f),
			new PropertySpecification("stereoCenters", "Stereo center count", "1", "3", 2f, 0f, 8f),
			new PropertySpecification("smallRings", "Ring count (<= 7 atoms)", "2", "", 1f, 0f, 10f),
			new PropertySpecification("aromaticRings", "Aromatic ring count", "", "2", 1f, 0f, 6f),
			new PropertySpecification("basicN", "Basic nitrogen count", "1", "", 0.5f, 0f, 6f),
			new PropertySpecification("acidicO", "Acidic oxygen count", "1", "", 0.5f, 0, 6f),
	};

	public static float calculateProperty(StereoMolecule mol, int type) {
		return (type == MOLECULAR_PROPERTY_MOLWEIGHT) ? mol.getMolweight()
				: (type == MOLECULAR_PROPERTY_CLOGP) ? new CLogPPredictor().assessCLogP(mol)
				: (type == MOLECULAR_PROPERTY_CLOGS) ? new SolubilityPredictor().assessSolubility(mol)
				: (type == MOLECULAR_PROPERTY_TPSA) ? new PolarSurfaceAreaPredictor().assessPSA(mol)
				: (type == MOLECULAR_PROPERTY_HDONORS) ? getHDonorCount(mol)
				: (type == MOLECULAR_PROPERTY_HACCEPTORS) ? getHAcceptorCount(mol)
				: (type == MOLECULAR_PROPERTY_ROTATABLEBONDS) ? mol.getRotatableBondCount()
				: (type == MOLECULAR_PROPERTY_FLEXIBILITY) ? new MolecularFlexibilityCalculator().calculateMolecularFlexibility(mol)
				: (type == MOLECULAR_PROPERTY_COMPLEXITY) ? FastMolecularComplexityCalculator.assessComplexity(mol)
				: (type == MOLECULAR_PROPERTY_SHAPE) ? MolecularShapeCalculator.assessShape(mol)
				: (type == MOLECULAR_PROPERTY_STEREOCENTERS) ? mol.getStereoCenterCount()
				: (type == MOLECULAR_PROPERTY_SMALLRINGCOUNT) ? mol.getRingSet().getSize()
				: (type == MOLECULAR_PROPERTY_AROMRINGCOUNT) ? mol.getAromaticRingCount()
				: (type == MOLECULAR_PROPERTY_BASIC_NITROGENS) ? getBasicNitrogenCount(mol)
				: (type == MOLECULAR_PROPERTY_ACIDIC_OXYGENS) ? getAcidicOxygenCount(mol)
				: Float.NaN;
	}

	public static int getPropertyCount() {
		return SPEC.length;
	}

	public static String getPropertyName(int type) {
		return SPEC[type].name;
	}

	public static String getPropertyCode(int type) {
		return SPEC[type].code;
	}

	public static String getPreferredMin(int type) {
		return SPEC[type].min;
	}

	public static String getPreferredMax(int type) {
		return SPEC[type].max;
	}

	public static float getRangeMin(int type) {
		return SPEC[type].rangeMin;
	}

	public static float getRangeMax(int type) {
		return SPEC[type].rangeMax;
	}

	public static double getValuation(double value, double min, double max, double halfWidth) {
		if (Double.isNaN(value))
			return Double.NaN;

		double v = 1.0;
		if (!Double.isNaN(min))
			v *= 1.0 / (1.0 + Math.exp((min - value) / halfWidth));
		if (!Double.isNaN(max))
			v *= 1.0 / (1.0 + Math.exp((value - max) / halfWidth));
		return v;
	}

	public static int getTypeFromCode(String code) {
		for (int i=0; i<SPEC.length; i++)
			if (SPEC[i].code.equals(code))
				return i;
		return -1;
	}

	public static int getTypeFromName(String name) {
		for (int i=0; i<SPEC.length; i++)
			if (SPEC[i].name.equals(name))
				return i;
		return -1;
	}

	public static String getMinText(int type) {
		return SPEC[type].min;
	}

	public static String getMaxText(int type) {
		return SPEC[type].max;
	}

	public static float getHalfFitnessWidth(int type) {
		return SPEC[type].halfWidth;
	}

	private static int getHAcceptorCount(StereoMolecule mol) {
		int count = 0;
		for (int atom = 0; atom < mol.getAllAtoms(); atom++)
			if (mol.getAtomicNo(atom) == 7 || mol.getAtomicNo(atom) == 8)
				count++;
		return count;
	}

	private static int getHDonorCount(StereoMolecule mol) {
		int count = 0;
		for (int atom = 0; atom < mol.getAllAtoms(); atom++)
			if ((mol.getAtomicNo(atom) == 7 || mol.getAtomicNo(atom) == 8) && mol.getAllHydrogens(atom) > 0)
				count++;
		return count;
	}

	private static int getBasicNitrogenCount(StereoMolecule mol) {
		int count = 0;
		mol.ensureHelperArrays(Molecule.cHelperRings);
		for (int atom=0; atom<mol.getAtoms(); atom++)
			if (AtomFunctionAnalyzer.isBasicNitrogen(mol, atom))
				count++;
		return count;
	}

	private static int getAcidicOxygenCount(StereoMolecule mol) {
		int count = 0;
		mol.ensureHelperArrays(Molecule.cHelperRings);
		for (int atom=0; atom<mol.getAtoms(); atom++)
			if (AtomFunctionAnalyzer.isAcidicOxygen(mol, atom))
				count++;
		return count;
	}

}

class PropertySpecification {
	float halfWidth,rangeMin,rangeMax;	// rangeMin and rangeMax should define a typical range to display the distribution of this property
	String code,name,min,max;

	public PropertySpecification(String code, String name, String min, String max, float halfWidth, float rangeMin, float rangeMax) {
		this.code = code;
		this.name = name;
		this.min = min;
		this.max = max;
		this.halfWidth = halfWidth;
		this.rangeMin = rangeMin;
		this.rangeMax = rangeMax;
	}
}
