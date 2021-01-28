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

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;


public class TotalSurfaceAreaPredictor extends PolarSurfaceAreaPredictor {
	public static final float cPSAUnknown = -1.0f;

	public static final String[] cNonPolarAtomTypeName = {
		"[B](-*)(-*)-*",
		"[BH2]-*",
		"[B-](-*)(-*)(-*)-*",	//x
		"[C](-*)(-*)(-*)-*",
		"[C](-*)(-*)=*",
		"[C](=*)=*",
		"[C](-*)#*",
		"[CH](-*)(-*)-*",
		"[CH](-*)=*",
		"[CH]#*",
		"[CH2](-*)-*",
		"[CH2]=*",
		"[CH3]-*",
		"[c](:*)(:*):*",
		"[c](-*)(:*):*",
		"[cH](:*):*",
		"[F]-*",
		"[Cl]-*",
		"[Br]-*",
		"[I]-*",
		};

	/* These are the increments from Actelion3D surface calculation (VdW radii, 1.4A probe),
	 * which creates somewhat smaller values than the Schroedinger method.
	 * The average error of the PLS prediction (ChEMBL training set) was about 7 square angstrom. */
	private static final float[] cPolarIncrement = {
		 3.55f, 11.51f, 13.56f,  4.34f, 25.88f,  1.55f, 11.46f, 21.62f, 15.60f, 15.27f,
		-7.62f, -5.20f, 17.76f, 17.59f,  0.00f,  7.74f, -1.15f, 14.48f, 10.97f,  7.95f,
		 6.87f,  1.22f, 13.99f,-18.78f,  2.76f,  0.00f, 10.00f, 17.65f, 13.04f, 13.10f,
		14.28f, 14.11f, 17.49f, 28.00f,  9.68f,  2.50f, 22.61f, 20.36f,  0.00f,  0.00f,
		 0.00f,  3.98f,  8.45f }; 

	private static final float[] cNonPolarIncrement = {
		 7.85f, 20.62f,  0.04f, -3.09f,  4.71f, 12.52f, 14.90f,  5.75f, 13.25f, 15.66f,
		13.76f, 18.04f, 19.01f,  5.09f,  5.46f, 12.21f, 13.10f, 22.17f, 25.38f, 33.03f };

	/* These are the increments approximating the Schroedinger method (VdW radii, 1.4A probe),
	 * which creates somewhat smaller values than the Schroedinger method.
	 * The average error of the PLS prediction (ChEMBL training set) was about 10 square angstrom.
	private static final float[] cPolarIncrement = {
		 3.96f, 12.12f,  7.11f,  3.65f, 18.40f, -0.46f, 12.34f, 18.64f, 13.41f, 13.51f,
	   -10.84f, -0.86f,  1.55f,  6.00f,  0.00f,  2.10f,  0.10f, 17.22f, 10.85f,  7.53f,
	     6.85f,  0.72f, 13.09f, -0.20f,  4.69f,  0.00f, 10.40f, 16.13f, 11.87f, 12.49f,
	    11.65f, 13.61f, 14.99f, 24.11f,  7.06f,  1.74f, 20.31f, 17.97f,  0.00f,  0.00f,
	     0.00f,  4.48f,  5.87f };

	private static final float[] cNonPolarIncrement = {
		 3.11f, -0.00f, -2.45f, -3.22f,  5.06f, 10.02f, 16.42f,  5.82f, 14.71f, 13.37f,
		15.93f, 17.64f, 20.61f,  6.15f,  6.78f, 12.69f, 11.75f, 17.60f, 21.76f, 28.12f }; 
	*/

	public TotalSurfaceAreaPredictor() {
		}

	public static int getNonPolarAtomTypeCount() {
		return cNonPolarAtomTypeName.length;
		}

	/**
	 * Calculates the topological total surface area (TPSA) of a molecule as a sum of
	 * contributions of its polar and non-atom-types. This method the contributions of all polar
	 * atoms. Contributions of polar and non-polar atom types were calculated as follows:
	 * 92879 unique reliable active structures with molecular weights between 200 and 700
	 * were selected from the ChEMBL_20 database. For every structure one conformer was
	 * generated and minimized by the MMFF94 forcefield. Total surface areas were calculated
	 * as solvent accessible surface using Van der Waals radii and a probe of 1.4 angstrom
	 * radius. Surface contribution values were for 43 polar and 20 non polar atom types were
	 * calculated as a multi-linear regression with a PLS approach.
	 * @param mol
	 * @return
	 */
	public float assessTotalSurfaceArea(StereoMolecule mol) {
		return assessNonPolarSurfaceArea(mol) + assessPolarSurfaceArea(mol);
		}

	/**
	 * Calculates the topological polar surface area (TPSA) of a molecule as a sum of
	 * contributions of its polar atom-types. This method the contributions of all polar
	 * atoms. Contributions of polar and non-polar atom types were calculated as follows:
	 * 92879 unique reliable active structures with molecular weights between 200 and 700
	 * were selected from the ChEMBL_20 database. For every structure one conformer was
	 * generated and minimized by the MMFF94 forcefield. Total surface areas were calculated
	 * as solvent accessible surface using Van der Waals radii and a probe of 1.4 angstrom
	 * radius. Surface contribution values were for 43 polar and 20 non polar atom types were
	 * calculated as a multi-linear regression with a PLS approach.
	 * @param mol
	 * @return
	 */
	public float assessPolarSurfaceArea(StereoMolecule mol) {
		int[] count = getPolarAtomTypeCounts(mol);

		float psa = 0.0f;
		for (int i=0; i<cPolarIncrement.length; i++)
			psa += count[i] * cPolarIncrement[i];

		return psa;
		}

	/**
	 * Calculates the topological non-polar surface area (TPSA) of a molecule as a sum of
	 * contributions of its non-polar atom-types. This method the contributions of all polar
	 * atoms. Contributions of polar and non-polar atom types were calculated as follows:
	 * 92879 unique reliable active structures with molecular weights between 200 and 700
	 * were selected from the ChEMBL_20 database. For every structure one conformer was
	 * generated and minimized by the MMFF94 forcefield. Total surface areas were calculated
	 * as solvent accessible surface using Van der Waals radii and a probe of 1.4 angstrom
	 * radius. Surface contribution values were for 43 polar and 20 non polar atom types were
	 * calculated as a multi-linear regression with a PLS approach.
	 * @param mol
	 * @return
	 */
	public float assessNonPolarSurfaceArea(StereoMolecule mol) {
		int[] count = getNonPolarAtomTypeCounts(mol);

		float npsa = 0.0f;
		for (int i=0; i<cNonPolarIncrement.length; i++)
			npsa += count[i] * cNonPolarIncrement[i];

		return npsa;
		}

	/**
	 * Calculates the relative (fractional) polar surface area from polar and non-polar atom
	 * contributions. This method does not use the Peter Ertl increments.
	 * @param mol
	 * @return
	 */
	public float assessRelativePolarSurfaceArea(StereoMolecule mol) {
		float psa = assessPolarSurfaceArea(mol);
		float npsa = assessNonPolarSurfaceArea(mol);
		return psa / (psa + npsa);
		}

	public ParameterizedStringList getDetail(StereoMolecule mol) {
		ParameterizedStringList detail = new ParameterizedStringList();
		detail.add("The total surface area prediction is based on an atom-type based",
							ParameterizedStringList.cStringTypeText);
		detail.add("increment system. Recognized atom types and their contributions are:",
							ParameterizedStringList.cStringTypeText);

		addNonPolarAtomTypeIncrements(mol, detail);
		addPolarAtomTypeIncrements(mol, detail);

		return detail;
		}

	private void addPolarAtomTypeIncrements(StereoMolecule mol, ParameterizedStringList detail) {
		int[] count = getPolarAtomTypeCounts(mol);

		for (int i=0; i<cPolarIncrement.length; i++)
			if (count[i] != 0)
				detail.add(""+count[i]+" * "+cPolarIncrement[i]+"   AtomType: "
								  +cPolarAtomTypeName[i],ParameterizedStringList.cStringTypeText);
		}

	private void addNonPolarAtomTypeIncrements(StereoMolecule mol, ParameterizedStringList detail) {
		int[] count = getNonPolarAtomTypeCounts(mol);

		for (int i=0; i<cNonPolarIncrement.length; i++)
			if (count[i] != 0)
				detail.add(""+count[i]+" * "+cNonPolarIncrement[i]+"   AtomType: "
								  +cNonPolarAtomTypeName[i],ParameterizedStringList.cStringTypeText);
		}

	public int[] getNonPolarAtomTypeCounts(StereoMolecule mol) {
		int[] count = new int[cNonPolarIncrement.length+2];

		mol.ensureHelperArrays(Molecule.cHelperRings);

		for (int atom=0; atom<mol.getAtoms(); atom++)
			count[getNonPolarAtomType(mol, atom)]++;

		return count;
		}

	private int getNonPolarAtomType(StereoMolecule mol, int atom) {
		switch (mol.getAtomicNo(atom)) {
		case 5:
			if (mol.getConnAtoms(atom) == 3
			 && mol.getAtomCharge(atom) == 0
			 && mol.getAtomPi(atom) == 0)
				return 0;
			if (mol.getConnAtoms(atom) == 1
			 && mol.getAtomCharge(atom) == 0
			 && mol.getAtomPi(atom) == 0)
				return 1;
			if (mol.getConnAtoms(atom) == 4
			 && mol.getAtomCharge(atom) == -1
			 && mol.getAtomPi(atom) == 0)
				return 2;
			return cNonPolarIncrement.length+1;	// unrecognized B
		case 6:
			if (mol.isAromaticAtom(atom)) {
				if (mol.getAllHydrogens(atom) != 0)
					return 15;
				for (int i=0; i<mol.getConnAtoms(atom); i++)
					if (!mol.isAromaticBond(mol.getConnBond(atom, i)))
						return 14;
				return 13;
				}
			else {
				if (mol.getAtomCharge(atom) == 0) {
					switch (mol.getAllHydrogens(atom)) {
					case 0:	// hydrogens
						switch (mol.getAtomPi(atom)) {
						case 0:	// pi
							return 3;
						case 1:	// pi
							return 4;
						case 2:	// pi
							if (mol.getConnBondOrder(atom, 0) == 2)
								return 5;
							else
								return 6;
							}
						break;
					case 1:	// hydrogens
						switch (mol.getAtomPi(atom)) {
						case 0:	// pi
							return 7;
						case 1:	// pi
							return 8;
						case 2:	// pi
							return 9;
							}
						break;
					case 2:	// hydrogens
						switch (mol.getAtomPi(atom)) {
						case 0:	// pi
							return 10;
						case 1:	// pi
							return 11;
							}
					case 3:	// hydrogens
						return 12;
						}
					}
				}
			return cNonPolarIncrement.length+1;	// unrecognized C
		case 9:	// F
			return 16;
		case 17:	// Cl
			return 17;
		case 35:	// Br
			return 18;
		case 53:	// I
			return 19;
		default:
			return cNonPolarIncrement.length;	// undefined type
			}
		}
	}
