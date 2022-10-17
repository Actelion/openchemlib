package com.actelion.research.chem.contrib;


import com.actelion.research.chem.*;

public class DiastereotopicAtomID {

	private static int xAtomicNumber = Molecule.getAtomicNoFromLabel("X", Molecule.cPseudoAtomX);

	public static String[] getAtomIds(StereoMolecule molecule) {
		addMissingChirality(molecule);
		
		int numberAtoms=molecule.getAllAtoms();
		String[] ids=new String[numberAtoms];
		StereoMolecule tempMolecule;
		for (int iAtom=0; iAtom<numberAtoms; iAtom++) {
			tempMolecule=molecule.getCompactCopy();
			changeAtom(tempMolecule, iAtom);
			makeRacemic(tempMolecule);
			// We need to ensure the helper array in order to get correctly the result of racemisation
			ids[iAtom]=(new Canonizer(tempMolecule, Canonizer.ENCODE_ATOM_CUSTOM_LABELS)).getIDCode();
		}
		return ids;
	}

	/**
	 * The problem is that sometimes we need to add chiral bond that was not planned because it is the same group
	 * This is the case for example for the valine where the 2 C of the methyl groups are diastereotopic
	 * @param molecule
	 */
	public static void addMissingChirality(StereoMolecule molecule, int esrType) {
		for (int iAtom=0; iAtom<molecule.getAllAtoms(); iAtom++) {
			StereoMolecule tempMolecule=molecule.getCompactCopy();
			changeAtomForStereo(tempMolecule, iAtom);
			// After copy, helpers must be recalculated
			tempMolecule.ensureHelperArrays(Molecule.cHelperParities);
			// We need to have >0 and not >1 because there could be unspecified chirality in racemate
			for (int i=0; i<tempMolecule.getAtoms(); i++) {
				// changed from from handling below; TLS 9.Nov.2015
				if (tempMolecule.isAtomStereoCenter(i) && tempMolecule.getStereoBond(i)==-1) {
					int stereoBond = tempMolecule.getAtomPreferredStereoBond(i);
					if (stereoBond != -1) {
						molecule.setBondType(stereoBond, Molecule.cBondTypeUp);
						if (molecule.getBondAtom(1, stereoBond) == i) {
							int connAtom = molecule.getBondAtom(0, stereoBond);
							molecule.setBondAtom(0, stereoBond, i);
							molecule.setBondAtom(1, stereoBond, connAtom);
						}
					// To me it seems that we have to add all stereo centers into AND group 0.	TLS 9.Nov.2015
					molecule.setAtomESR(i, esrType, 0);
					}
				}
			}
		}
	}

	public static void addMissingChirality(StereoMolecule molecule) {
		addMissingChirality(molecule, Molecule.cESRTypeAnd);
	}

	private static void changeAtomForStereo(StereoMolecule molecule, int iAtom) {
		// need to force the change to an excotic atom to check if it is really chiral
		molecule.setAtomicNo(iAtom, xAtomicNumber);
	}

	private static void changeAtom(StereoMolecule molecule, int iAtom) {
		molecule.setAtomCustomLabel(iAtom, molecule.getAtomLabel(iAtom)+"*");
		if (molecule.getAtomicNo(iAtom)==1) {
			molecule.setAtomicNo(iAtom, xAtomicNumber);
		} else {
			// we can not use X because we would have problems with valencies if it is
			// expanded hydrogens or not
			// we can not only use a custom label because it does not cound for the canonisation
			molecule.setAtomMass( iAtom,molecule.getAtomMass(iAtom)+5);
		}
	}


	private static void makeRacemic(StereoMolecule molecule) {
		// if we don't calculate this we have 2 epimers
		molecule.ensureHelperArrays(Molecule.cHelperParities);
		// we need to make one group "AND" for chiral (to force to racemic, this means diastereotopic and not enantiotopic)
		for (int i=0; i<molecule.getAllAtoms(); i++) {
			if (molecule.getAtomParity(i)!=Molecule.cAtomParityNone) {
				molecule.setAtomESR(i, Molecule.cESRTypeAnd, 0);	// changed to group 0; TLS 9.Nov.2015
			}
		}
	}

	public static String test(StereoMolecule molecule) {		// changed to use StringBuilder; TLS 9.Nov.2015
		String[] ids=getAtomIds(molecule);
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<ids.length; i++)
			sb.append(Integer.toString(i)).append(" - ").append(ids[i]).append("\n");
		return sb.toString();
	}
	
	/**
	 * In order to debug we could number the group of diastereotopic atoms
	 * 
	 * @param molecule
	 */
	public static void markDiastereotopicAtoms(StereoMolecule molecule) {		// changed from markDiastereo(); TLS 9.Nov.2015
		String[] ids = getAtomIds(molecule);
		SortedStringList analyzed = new SortedStringList();
		int group=0;
		for (String id:ids) {
			System.out.println(id+" - "+group);
			if (!analyzed.contains(id)) {
				analyzed.addString(id);
				for (int iAtom=0; iAtom<ids.length; iAtom++) {
					if (id.equals(ids[iAtom])) {
						molecule.setAtomCustomLabel(iAtom, Integer.toString(group));
					}
				}
				group++;
			}
		}
	}
}
