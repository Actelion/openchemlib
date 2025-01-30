package com.actelion.research.chem;

import com.actelion.research.chem.coords.CoordinateInventor;

/**
 * MoleculeStandardizer
 * @author Modest von Korff, Thomas Sander
 * @version 1.0
 * Apr 5, 2012 MvK: Start implementation
 * Oct 2020 MvK,TLS: improved performance and adapted standardization according to following publication:
 * Bento, A. P., Hersey, A., Félix, E., Landrum, G., Gaulton, A., Atkinson, F., ... & Leach, A. R. (2020).
 * An open source chemical structure curation pipeline using RDKit. Journal of Cheminformatics, 12(1), 1-16.<br>
 * Exceptions: - S=O is not transformed into S(+)-O(-)<br>
 *             - If charges (e.g. quaternary nitrogen) cannot be balanced, then Na(+) or Cl(-) are added to neutralize as last resort<br>
 */
public class MoleculeStandardizer {
	public static final int MODE_STANDARDIZE_ONLY = 0;
	public static final int MODE_LARGEST_FRAGMENT = 1;
	public static final int MODE_REMOVE_ISOTOPS = 2;
	public static final int MODE_ADD_NA_AND_CL = 4;
	public static final int MODE_PROHIBIT_REMAINING_CHARGE = 8;

	public static final int MODE_GET_PARENT = MODE_LARGEST_FRAGMENT + MODE_REMOVE_ISOTOPS;	// for compatibility

	/**
	 * Under normal circumstances, one should never need to standardize a molecule from an idcode,
	 * because molecules should be standardized before generating the canonical encoding.
	 * An exception is when generating the parent structure using mode MODE_GET_PARENT and
	 * potentially MODE_ADD_NA_AND_CL.
	 * @param idcode
	 * @param coordinates if null the result may change.
	 * @param mode 0 or any combination of MODE_LARGEST_FRAGMENT, MODE_REMOVE_ISOTOPS, and MODE_ADD_NA_AND_CL
	 * @return
	 * @throws Exception
	 */
	public static StereoMolecule getStandardized(String idcode, String coordinates, int mode) throws Exception {
		StereoMolecule mol = new IDCodeParser().getCompactMolecule(idcode, coordinates);
		standardize(mol, mode);
		return mol;
	}

	/**
	 * Standardises a molecule and fixes some structural errors.
	 * Typically, this is done before canonicalization.
	 * It includes the following changes:<br>
	 * - different forms of functional groups (e.g. nitro) are normalized to a preferred one<br>
	 * - charged acidic or basic atoms are (de-)protonated to remove charges and neutralize the molecule, if possible.<br>
	 * - alkali/earthalkali/halogene atoms, if charged despite being covalently bound, get uncharged<br>
	 * - trivalent, uncharged oxygens get a positive charge<br>
	 * - unusual amide tautomeric structures, if not in a ring, are inverted<br>
	 * - uncharged isocyano groups get proper charges to validate valences<br>
	 * - wrongly charged azido groups get proper charges to validate valences<br>
	 * - uncharged, quarternary nitrogens get a positive charge<br>
	 * If mode includes MODE_GET_PARENT, then only the largest, normalized fragment is kept.
	 * If mode includes MODE_ADD_NA_AND_CL, then molecules, that are still charged after normalization,
	 * e.g. quarternary ammonium, are neutralized by adding the right amount of Na+ or Cl- ions.
	 * @param mol
	 * @param mode 0 or any combination of MODE_LARGEST_FRAGMENT, MODE_REMOVE_ISOTOPS, and MODE_ADD_NA_AND_CL, MODE_PROHIBIT_REMAINING_CHARGE
	 * @throws Exception
	 */
	public static void standardize(StereoMolecule mol, int mode) throws Exception {
		if ((mode & MODE_LARGEST_FRAGMENT) != 0)
			mol.stripSmallFragments();

		if ((mode & MODE_REMOVE_ISOTOPS) != 0)
			mol.stripIsotopInfo();

		repairAndUnify(mol);

		mol.normalizeAmbiguousBonds();

		int remainingCharge = mol.canonizeCharge(true, true);

		if (remainingCharge != 0)
			remainingCharge = neutralizeCharges(mol, mode, remainingCharge);

		if ((mode & MODE_PROHIBIT_REMAINING_CHARGE) != 0 && remainingCharge != 0)
			throw new Exception("Couldn't neutralize molecule.");

		mol.validateAtomQueryFeatures();
		mol.validateBondQueryFeatures();
	}

	/**
	 * Repairs wrongly uncharged quaternary nitrogen. Unifies carbonyl acid groups,
	 * sulfonic acid, phosphoric acid, phenolic oxygen. Means: negative charges are removed.
	 * Adds Na+ or Cl- for final charge equilibration.
	 * @param mol
	 */
	private static void repairAndUnify(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);

		repairCovalentBoundChargedAlkaliAndHalogen(mol);
		chargeTrivalentOxygen(mol);
		repairBadAmideTautomer(mol);
		repairQuaternaryNitrogen(mol);
		unifyIsoCyano(mol);
		unifyAzido(mol);
	}

	private static int neutralizeCharges(StereoMolecule mol, int mode, int totalCharge) {
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);

		for (int atom=0; atom<mol.getAllAtoms() && totalCharge>0; atom++) {
			if (AtomFunctionAnalyzer.isAcidicOxygen(mol, atom)) {
				mol.setAtomCharge(atom, -1);
				totalCharge--;
			}
		}

		for (int atom=0; atom<mol.getAllAtoms() && totalCharge<0; atom++) {
			if (AtomFunctionAnalyzer.isBasicNitrogen(mol, atom)) {
				mol.setAtomCharge(atom, 1);
				totalCharge++;
			}
		}

		if (totalCharge != 0 && (mode & MODE_ADD_NA_AND_CL) != 0) {
			for (int atom=0; atom<mol.getAllAtoms(); atom++)
				mol.setAtomMarker(atom, true);

			// add Cl-
			while (totalCharge > 0) {
				int ind = mol.addAtom(17);
				mol.setAtomCharge(ind, -1);
				totalCharge--;
			}

			// add Na+
			while (totalCharge < 0) {
				int ind = mol.addAtom(11);
				mol.setAtomCharge(ind, 1);
				totalCharge++;
			}

			new CoordinateInventor(CoordinateInventor.MODE_KEEP_MARKED_ATOM_COORDS
					+ CoordinateInventor.MODE_REMOVE_HYDROGEN).invent(mol);
		}

		return totalCharge;
	}

	/**
	 * Remove wrong charges on halogen and (earth)alkali atoms, if they are
	 * covalently bound.
	 * @param mol
	 */
	private static void repairCovalentBoundChargedAlkaliAndHalogen(StereoMolecule mol) {
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mol.isHalogene(atom)) {
				if (mol.getOccupiedValence(atom) == 1
				 && mol.getAtomCharge(atom) == -1) {
					mol.setAtomCharge(atom, 0);
					mol.setAtomAbnormalValence(atom, -1);
					}
				continue;
				}

			if (mol.isAlkaliMetal(atom)) {
				if (mol.getOccupiedValence(atom) == 1
				 && mol.getAtomCharge(atom) == 1) {
					mol.setAtomCharge(atom, 0);
					mol.setAtomAbnormalValence(atom, -1);
					}
				continue;
				}

			if (mol.isEarthAlkaliMetal(atom)) {
				if (mol.getOccupiedValence(atom) == 2
				 && mol.getAtomCharge(atom) == 2) {
					mol.setAtomCharge(atom, 0);
					mol.setAtomAbnormalValence(atom, -1);
					}
				continue;
				}
			}
		}

	private static void chargeTrivalentOxygen(StereoMolecule mol) {
		for (int atom=0; atom<mol.getAtoms(); atom++)
			if ((mol.getAtomicNo(atom) == 8
			  || mol.getAtomicNo(atom) == 16)
			 && mol.getOccupiedValence(atom) == 3
			 && mol.getAtomCharge(atom) != 1)
				mol.setAtomCharge(atom,1);
		}

	private static void repairBadAmideTautomer(StereoMolecule mol) {
		for (int oxygen=0; oxygen<mol.getAtoms(); oxygen++) {
			if (mol.getAtomicNo(oxygen) == 8
			 && mol.getConnAtoms(oxygen) == 1
			 && mol.getConnBondOrder(oxygen, 0) == 1) {
				int carbon = mol.getConnAtom(oxygen, 0);
				if (mol.getAtomicNo(carbon) == 6
				 && mol.getAtomPi(carbon) == 1) {
					for (int i=0; i<mol.getConnAtoms(carbon); i++) {
						if (mol.getConnBondOrder(carbon, i) == 2) {
							int nitrogen = mol.getConnAtom(carbon, i);
							if (mol.getAtomicNo(nitrogen) == 7
							 && !mol.isRingAtom(nitrogen)) {
								boolean hasResonance = false;
								for (int j=0; j<mol.getConnAtoms(nitrogen); j++) {
									int connAtom = mol.getConnAtom(nitrogen, j);
									if (connAtom != carbon
									 && mol.getAtomPi(connAtom) != 0) {
										hasResonance = true;
										break;
									}
								}
								if (!hasResonance) {
									mol.setBondType(mol.getConnBond(oxygen, 0), Molecule.cBondTypeDouble);
									mol.setBondType(mol.getConnBond(carbon, i), Molecule.cBondTypeSingle);
									break;
								}
							}
						}
					}
				}
			}
		}
	}

	/**
	 * -N#C, N positive charged, C negative charged.
	 * @param mol
	 */
	private static void unifyIsoCyano(StereoMolecule mol) {
		for (int bond=0; bond<mol.getBonds(); bond++) {
			if (mol.getBondType(bond) == Molecule.cBondTypeTriple) {
				for (int i=0; i<2; i++) {
					int atom1 = mol.getBondAtom(i, bond);
					int atom2 = mol.getBondAtom(1-i, bond);
					if (mol.getAtomicNo(atom1) == 7
					 && mol.getConnAtoms(atom1) == 2
					 && mol.getAtomicNo(atom2) == 6
					 && mol.getConnAtoms(atom2) == 1) {
						if (mol.getAtomCharge(atom1) != 1)
							mol.setAtomCharge(atom1, 1);
						if (mol.getAtomCharge(atom2) != -1)
							mol.setAtomCharge(atom2, -1);
						break;
					}
				}
			}
		}
	}
	
	private static void unifyAzido(StereoMolecule mol) {
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mol.getAtomicNo(atom) == 7
			 && mol.getConnAtoms(atom) == 2
			 && mol.getConnBondOrder(atom, 0) == 2
			 && mol.getConnBondOrder(atom, 1) == 2
			 && mol.getAtomicNo(mol.getConnAtom(atom, 0)) == 7
			 && mol.getAtomicNo(mol.getConnAtom(atom, 1)) == 7) {
				for (int i=0; i<2; i++) {
					int atom1 = mol.getConnAtom(atom, i);
					int atom2 = mol.getConnAtom(atom, 1-i);
					if (mol.getConnAtoms(atom1) == 1
					 && mol.getConnAtoms(atom2) == 2) {
						if (mol.getAtomCharge(atom) != 1)
							mol.setAtomCharge(atom, 1);
						if (mol.getAtomCharge(atom1) != -1)
							mol.setAtomCharge(atom1, -1);
						if (mol.getAtomCharge(atom2) != 0)
							mol.setAtomCharge(atom2, 0);
						break;
					}
				}
			}
		}
	}
	
	/**
	 * Sets the charge to +1.
	 * @param mol
	 */
	private static void repairQuaternaryNitrogen(StereoMolecule mol){
		for (int i=0; i<mol.getAllAtoms(); i++)
			if (mol.getAtomicNo(i) == 7
			 && mol.getOccupiedValence(i) == 4
			 && mol.getAtomCharge(i) != 1)
				mol.setAtomCharge(i, 1);
	}
}
