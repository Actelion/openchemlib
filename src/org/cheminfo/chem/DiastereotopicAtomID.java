package org.cheminfo.chem;


import java.util.Vector;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.ExtendedMolecule;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

public class DiastereotopicAtomID {
	
	public static String[] getAtomIds(StereoMolecule molecule) {		
		addMissingChirality(molecule);
		
		int numberAtoms=molecule.getAllAtoms();
		String[] ids=new String[numberAtoms];
		StereoMolecule tempMolecule;
		for (int iAtom=0; iAtom<numberAtoms; iAtom++) {
			tempMolecule=molecule.getCompactCopy();
			
			// Temporary bug fix
			tempMolecule.ensureHelperArrays(ExtendedMolecule.cHelperCIP);

			changeAtom(tempMolecule, iAtom);
			makeRacemic(tempMolecule);
			ids[iAtom]=(new Canonizer(tempMolecule, Canonizer.ENCODE_ATOM_CUSTOM_LABELS)).getIDCode();
		}
		return ids;
	}

	/**
	 * The problem is that sometimes we need to add chiral bond that was not planned because it is the same group
	 * This is the case for example for the valine where the 2 C of the methyl groups are diastereotopic
	 * @param molecule
	 */
	public static void addMissingChirality(StereoMolecule molecule) {

		
		int numberAtoms=molecule.getAllAtoms();
		StereoMolecule tempMolecule;
		

		for (int iAtom=0; iAtom<numberAtoms; iAtom++) {
			tempMolecule=molecule.getCompactCopy();
			// after copy we need to recalculate the helpers ...
			tempMolecule.ensureHelperArrays(Molecule.cHelperCIP);
			
			changeAtom(tempMolecule, iAtom);
			// we need to have >0 and not >1 because there could be unspecified chirality in racemate
			if (tempMolecule.getStereoCenterCount()>0) {
				for (int i=0; i<tempMolecule.getAtoms(); i++) {
					if ((tempMolecule.getStereoProblem(i)) && (tempMolecule.getStereoBond(i)==-1)) {
						// we need to add a stereobond
						// we find a single bond that we could convert ...
						tempMolecule.ensureHelperArrays(Molecule.cHelperNeighbours);
						int bondToChange=-1;	// We will prefer to change a bond that is not in the ring !
				        for (int j=0; j<tempMolecule.getAllConnAtoms(i); j++) {
				            int bond = tempMolecule.getConnBond(i, j);
				            if ((tempMolecule.getBondType(bond)==Molecule.cBondTypeSingle) && (tempMolecule.getBondAtom(0,bond)==i)) {
				            	bondToChange=bond;
				// TODO: reactivate the change
				            	if (tempMolecule.getBondRingSize(bond)==0) {
					            	molecule.changeBond(bond, Molecule.cBondTypeUp);
					            	bondToChange=-1;
									break;
				            	}
				            }
				        }
				        if (bondToChange>-1) {
				 //       	molecule.changeBond(bondToChange, Molecule.cBondTypeUp);
				        }
					}
				}
			}
		}
	}
	
	
	private static void changeAtom(StereoMolecule molecule, int iAtom) {
			molecule.setAtomCustomLabel(iAtom, molecule.getAtomLabel(iAtom)+"*");
			molecule.setAtomicNo(iAtom, Molecule.getAtomicNoFromLabel("X"));
	}
	
	private static void makeRacemic(StereoMolecule molecule) {
		// if we don't calculate this we have 2 epimers
		molecule.ensureHelperArrays(Molecule.cHelperCIP);
		// we need to make one group "AND" for chiral (to force to racemic, this means diastereotopic and not enantiotopic)
		for (int i=0; i<molecule.getAllAtoms(); i++) {
			if (molecule.getAtomParity(i)!=Molecule.cAtomParityNone) {
				molecule.setAtomESR(i, Molecule.cESRTypeAnd, 1);
			}
		}
	}
	
	
	
	public static String test(StereoMolecule molecule) {
		String[] ids=getAtomIds(molecule);
		String toReturn="";
		for (int i=0; i<ids.length; i++) {
			toReturn+=i+" - "+ids[i]+"\r\n";
		}
		return toReturn;
	}
	
	/**
	 * In order to debug we could number the group of diastereotopic atoms
	 * 
	 * @param molecule
	 */
	public static void markDiastereo(StereoMolecule molecule) {
		String[] ids=getAtomIds(molecule);
		Vector<String> analyzed=new Vector<String> ();
		int group=0;
		for (String id : ids) {
			System.out.println(id+" - "+group);
			if (!analyzed.contains(id)) {
				analyzed.add(id);
				for (int iAtom=0; iAtom<ids.length; iAtom++) {
					if (id.equals(ids[iAtom])) {
						molecule.setAtomCustomLabel(iAtom, group+"");
					}
				}
				group++;
			}
		}
	}
}
