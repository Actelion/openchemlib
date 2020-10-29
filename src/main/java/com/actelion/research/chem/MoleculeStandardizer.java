package com.actelion.research.chem;

import java.util.List;

/**
 * 
 * 
 * MoleculeStandardizer
 * <p>Copyright: Actelion Ltd., Inc. All Rights Reserved
 * This software is the proprietary information of Actelion Pharmaceuticals, Ltd.
 * Use is subject to license terms.</p>
 * @author Modest von Korff
 * @version 1.0
 * Apr 5, 2012 MvK: Start implementation
 * 22.10.2020 MvK: added standardization loosely following publication:
 * Bento, A. P., Hersey, A., Félix, E., Landrum, G., Gaulton, A., Atkinson, F., ... & Leach, A. R. (2020).
 * An open source chemical structure curation pipeline using RDKit. Journal of Cheminformatics, 12(1), 1-16.
 *
 * Considered was also the Python RDKit code from Greg on Github.
 *
 * Except: S=O is not transformed into S+-O-
 */
public class MoleculeStandardizer {

	// Nothing happens
	public static final int MODE_0 = 0;

	public static final int MODE_REPAIR_AND_BALANCE = 1;

	public static final int MODE_STRIP_SMALL = 2;

	private static final String idCodeNegNAtCarbon = "eF``zLFD@";
	
	private static final String idCodeNegNAt2Carbons = "eM``zN``";

	private static final String idCodeNegCarbon = "fHW@P";
	
	private static final String idCodeNegSulfur = "fHpPxB";
	
	private static final String idCodeNegOxygenAtCarbon = "eFH`zLD";
	
	private static final String idCodeNegOxygenAtP = "eFJh|KhpP";
	
	private static final String idCodeNegOxygenAtS = "eFJhBKhpP";

	private static final String idCodeIsoCyanoUncharged = "eM`BO``";

	private static final String idCodeAzidoUncharged = "gCl@ADeJpD";

	private static final String idCodeSulfOxideCharged = "eFJXBB[hpP";

	// Wrong drawn nitro group (no charges to N and O assigned)
	private static final String idCodeBadNitro = "gCh`HDdsPFDP";

	// Charged at N
	private static final String idCodeBadNitroCharged = "gChaHEIRYhCBH";

	// Wrong tautomer, oxygen not charged
	// OC=NC
	private static final String idCodeBadAmide = "gCi@DDefDeeD\u007FhK@";

	// Carboxyl group to detect acidic H.
	private static final String idCodeCarboxylGroup = "gC``@dfZ@~bnLqLG@HA@";

	private static MoleculeStandardizer instance;
	
	
	private IDCodeParser parser;



	private StereoMolecule molSulfOxideCharged;
	private StereoMolecule molBadAmide;
	private StereoMolecule molCarboxylicGroup;

	private StereoMolecule [] arrNeutralizeNegative;

	private StereoMolecule [] arrNeutralizeNegativeNitrogen;
	
	private StereoMolecule [] arrRepairNitro;
	
	private StereoMolecule [] arrUnifyAzido;
	
	private StereoMolecule [] arrUnifyIsoCyano;


	
	public MoleculeStandardizer() {
		
		parser = new IDCodeParser(false);
		
		arrNeutralizeNegativeNitrogen = initNeutralizeNegativeNitrogen();
		
		arrNeutralizeNegative = initNeutralizeNegative();
		
		arrRepairNitro = initRepairNitro();
		
		arrUnifyAzido = initUnifyAzido();
		
		arrUnifyIsoCyano = initUnifyIsoCyano();

		molSulfOxideCharged = parser.getCompactMolecule(idCodeSulfOxideCharged);
		molSulfOxideCharged.ensureHelperArrays(Molecule.cHelperRings);

		molBadAmide = parser.getCompactMolecule(idCodeBadAmide);
		molBadAmide.ensureHelperArrays(Molecule.cHelperRings);

		molCarboxylicGroup = parser.getCompactMolecule(idCodeCarboxylGroup);
		molCarboxylicGroup.ensureHelperArrays(Molecule.cHelperRings);

	}
	
	
	private StereoMolecule [] initNeutralizeNegativeNitrogen(){
		
		StereoMolecule [] arr = new StereoMolecule [2];
		
		int cc=0;
		
		arr[cc++]=parser.getCompactMolecule(idCodeNegNAtCarbon);
		
		arr[cc++]=parser.getCompactMolecule(idCodeNegNAt2Carbons);
		
		for (int i = 0; i < arr.length; i++) {
			arr[i].ensureHelperArrays(Molecule.cHelperRings);
		}
		
		return arr;
	}
	
	private StereoMolecule [] initNeutralizeNegative(){
		
		StereoMolecule [] arr = new StereoMolecule [5];
		
		int cc=0;
		
		arr[cc++]=parser.getCompactMolecule(idCodeNegCarbon);
		
		arr[cc++]=parser.getCompactMolecule(idCodeNegSulfur);
		
		arr[cc++]=parser.getCompactMolecule(idCodeNegOxygenAtCarbon);
		
		arr[cc++]=parser.getCompactMolecule(idCodeNegOxygenAtP);
		
		arr[cc++]=parser.getCompactMolecule(idCodeNegOxygenAtS);
		
		for (int i = 0; i < arr.length; i++) {
			arr[i].ensureHelperArrays(Molecule.cHelperRings);
		}
		
		return arr;
	}
	
	private StereoMolecule [] initRepairNitro(){
		
		StereoMolecule [] arr = new StereoMolecule [2];
		
		int cc=0;
		
		arr[cc++]=parser.getCompactMolecule(idCodeBadNitro);
		
		arr[cc++]=parser.getCompactMolecule(idCodeBadNitroCharged);
		
		for (int i = 0; i < arr.length; i++) {
			arr[i].ensureHelperArrays(Molecule.cHelperRings);
		}
		
		return arr;
	}
	
	private StereoMolecule [] initUnifyAzido(){
		
		StereoMolecule [] arr = new StereoMolecule [1];
		
		int cc=0;
		
		arr[cc++]=parser.getCompactMolecule(idCodeAzidoUncharged);
		
		for (int i = 0; i < arr.length; i++) {
			arr[i].ensureHelperArrays(Molecule.cHelperRings);
		}
		
		return arr;
	}
	
	private StereoMolecule [] initUnifyIsoCyano(){
		
		StereoMolecule [] arr = new StereoMolecule [1];
		
		int cc=0;
		
		arr[cc++]=parser.getCompactMolecule(idCodeIsoCyanoUncharged);
		
		for (int i = 0; i < arr.length; i++) {
			arr[i].ensureHelperArrays(Molecule.cHelperRings);
		}
		
		return arr;
	}
	
	/**
	 * Returns a fragment, no counter ion is given.
	 * @param idcode
	 * @param coordinates if null the result may change.
	 * @return
	 * @throws Exception
	 */
	public StereoMolecule getStandardized(String idcode, String coordinates) throws Exception {
		
		StereoMolecule mol = parser.getCompactMolecule(idcode, coordinates);
		
		return getStandardized(mol);

	}
	
	/**
	 * Returns a fragment, no counter ion is given.
	 * @param mol
	 * @return
	 * @throws Exception
	 */
	public StereoMolecule getStandardized(StereoMolecule mol) throws Exception {
		return getStandardized(mol, MODE_STRIP_SMALL);
	}

	public StereoMolecule getStandardized(StereoMolecule mol, int mode) throws Exception {

		StereoMolecule molStandard = new StereoMolecule(mol);

		molStandard.ensureHelperArrays(Molecule.cHelperRings);

		if(mode == MODE_0){
			return molStandard;
		}

		if(mode> MODE_REPAIR_AND_BALANCE) {
			// Strip small frags before repair
			molStandard.stripSmallFragments();
		}

		// Repair
		repairAndUnify(molStandard, mode);

		molStandard.normalizeAmbiguousBonds();

		// molStandard.canonizeCharge(false);


		molStandard.ensureHelperArrays(Molecule.cHelperRings);
		molStandard.canonizeCharge(true);

		molStandard.ensureHelperArrays(Molecule.cHelperRings);

		if(mode> MODE_REPAIR_AND_BALANCE) {
			// Remove everything we added for charge balancing.
			molStandard.stripSmallFragments();
		}

		return molStandard;
	}

	/**
	 * Repairs wrong charged quaternary nitrogen. Unifies carbonyl acid groups, 
	 * sulfonic acid, phosphoric acid, phenolic oxygen. Means: negative charges are removed.
	 * Adds Na+ or Cl- for final charge equilibration. Calls ensure helper arrays at the end.
	 * @param mol
	 * @return true if an atom was added.
	 */


	public void repairAndUnify(StereoMolecule mol, int mode) throws Exception {

		// ChEMBL creates charged sulfoxides S+-O-
		// Here the neutral form is generated S=O
		neutralizeSulfoxideCharged(mol);

		chargeTrivalentOxygen(mol);

		// Handelt by resolve ambigous bonds in Molecule
		// removeCovalentAlkalineBonds(mol);

		repairBadAmideTautomer(mol);

		// Has to be here.
		// Prepares molecule for canonize charges
		neutralizeNegative(mol);

		neutralizeNegativeNitrogen(mol);

		repairQuaternaryNitrogen(mol);

		repairTertiaryNitrogen(mol);

		repairNitro(mol);

		unifyIsoCyano(mol);

		unifyAzido(mol);

		// Prepares molecule for canonize charges
		// balanceCharges(mol, mode);
	}


	private boolean balanceCharges(StereoMolecule mol, int mode) throws Exception {

		boolean added = false;

		//
		// Balance charges
		//
		int totalCharge = 0;
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			totalCharge += mol.getAtomCharge(i);
		}


		//
		// Remove single atom positive ions
		//
		if(mode>MODE_REPAIR_AND_BALANCE) {
			while (totalCharge > 0) {
				boolean removed = false;
				for (int i = 0; i < mol.getAtoms(); i++) {
					if (totalCharge >= mol.getAtomCharge(i)) {
						if (mol.getConnAtoms(i) == 0) {
							mol.deleteAtom(i);
							mol.ensureHelperArrays(Molecule.cHelperRings);
							totalCharge--;
							removed = true;
						}
					}
				}
				if (!removed) {
					break;
				}
			}
		}

//		if(removeAcidicHydrogen(mol)){
//			totalCharge--;
//		}

		totalCharge = mol.canonizeCharge(true);

		//
		// add Cl-
		//
		while (totalCharge > 0) {
			int ind = mol.addAtom(17);
			mol.setAtomCharge(ind, -1);
			totalCharge--;
			added = true;
		}

		//
		// add Na+
		//
		while (totalCharge < 0) {
			int ind = mol.addAtom(11);
			mol.setAtomCharge(ind, 1);
			totalCharge++;
			added = true;
		}

		mol.ensureHelperArrays(Molecule.cHelperRings);

		return added;

	}

	private boolean removeAcidicHydrogen(StereoMolecule mol) {

		boolean removed = false;

		SSSearcher sss = new SSSearcher();
		sss.setMol(molCarboxylicGroup, mol);
		if(sss.findFragmentInMolecule()>0) {
			List<int[]> liArrMatch = sss.getMatchList();

			int [] arr = liArrMatch.get(0);

			for (int atIndex : arr) {

				if(mol.getAtomicNo(atIndex)==8){

					if(mol.getImplicitHydrogens(atIndex)==1){
						mol.setAtomCharge(atIndex, -1);
						removed = true;
						break;
					}
				}
			}
		}

		return removed;
	}


	/**
	 * ChEMBL publication
	 * @param mol
	 */
	private void chargeTrivalentOxygen(StereoMolecule mol) {
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i) == 8){
				int v = mol.getOccupiedValence(i);
				if(v==3){
					if(mol.getAtomCharge(i)!=1){
						mol.setAtomCharge(i,1);
					}
				}
			}
		}
		mol.ensureHelperArrays(Molecule.cHelperRings);
	}

	private void removeCovalentAlkalineBonds(StereoMolecule mol) {
		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(mol.getAtomicNo(i) == 7 || mol.getAtomicNo(i) == 8){

				int connAts = mol.getConnAtoms(i);

				for (int j = 0; j < connAts; j++) {

					int atConn = mol.getConnAtom(i, j);
					int atNoConn = mol.getAtomicNo(atConn);

					if(PeriodicTable.isAlkaline(atNoConn)){
						int bnd = mol.getBond(i, atConn);
						mol.deleteBond(bnd);

						int charge = mol.getAtomCharge(i);
						charge--;
						mol.setAtomCharge(i, charge);

						mol.setAtomCharge(atConn, 1);
					}
				}

			}
		}
		mol.ensureHelperArrays(Molecule.cHelperRings);
	}

	private void neutralizeNegativeNitrogen(StereoMolecule mol){
		
		for (int i = 0; i < arrNeutralizeNegativeNitrogen.length; i++) {

			StereoMolecule frag = arrNeutralizeNegativeNitrogen[i];
			
			SSSearcher sss = new SSSearcher();
			
			sss.setMol(frag, mol);
			
			if(sss.findFragmentInMolecule()>0){
				List<int []> liArrMatch = sss.getMatchList();
				
				if(liArrMatch == null){
					return;
				}
				
				for (int[] arrMatch : liArrMatch) {
					
					for (int j = 0; j < arrMatch.length; j++) {
						int atom = arrMatch[j];
						
						if(mol.getAtomicNo(atom)==7){
							mol.setAtomCharge(atom, 0);
						} 
					}
				}
			}
		}
		mol.ensureHelperArrays(Molecule.cHelperRings);
	}

	private void neutralizeSulfoxideCharged(StereoMolecule mol){

		SSSearcher sss = new SSSearcher();

		sss.setMol(molSulfOxideCharged, mol);

		if(sss.findFragmentInMolecule()>0){
			List<int []> liArrMatch = sss.getMatchList();

			if(liArrMatch == null){
				return;
			}

			for (int[] arrMatch : liArrMatch) {
				if(arrMatch.length>2){
					throw new RuntimeException("Error in logic of algorithm!");
				}

				int [] arrAtIndex = new int[2];

				for (int i = 0; i < arrMatch.length; i++) {
					int atom = arrMatch[i];
					if(mol.getAtomicNo(atom)==8){
						mol.setAtomCharge(atom, 0);
						arrAtIndex[i]=atom;
					} else if(mol.getAtomicNo(atom)==16){
						mol.setAtomCharge(atom, 0);
						arrAtIndex[i]=atom;
					} else {
						throw new RuntimeException("Error in logic of algorithm!");
					}
				}

				int bnd = mol.getBond(arrAtIndex[0],arrAtIndex[1]);

				mol.changeBond(bnd, Molecule.cBondTypeDouble);
			}
		}
		mol.ensureHelperArrays(Molecule.cHelperRings);
	}

	private void repairBadAmideTautomer(StereoMolecule mol){



		SSSearcher sss = new SSSearcher();

		boolean changed = true;

		while (changed) {
			changed = false;
			sss.setMol(molBadAmide, mol);

			if (sss.findFragmentInMolecule() > 0) {
				List<int[]> liArrMatch = sss.getMatchList();

				if (liArrMatch == null) {
					return;
				}

				for (int[] arrMatch : liArrMatch) {

					for (int i = 0; i < arrMatch.length; i++) {
						int atom = arrMatch[i];
						if (mol.getAtomicNo(atom) == 8) {
							mol.setAtomCharge(atom, 0);

							int atC = mol.getConnAtom(atom, 0);

							if (mol.getAtomicNo(atC) != 6) {
								throw new RuntimeException("Error in logic of algorithm!");
							}

							int connAtsC = mol.getConnAtoms(atC);

							for (int j = 0; j < connAtsC; j++) {
								int atConnConn = mol.getConnAtom(atC, j);

								if (mol.getAtomicNo(atConnConn) == 7) {

									int bndCToN = mol.getBond(atC, atConnConn);
									if (mol.getBondOrder(bndCToN) != 2) {
										throw new RuntimeException("Error in logic of algorithm!");
									}
									mol.changeBond(bndCToN, Molecule.cBondTypeSingle);
									int bndCToO = mol.getBond(atC, atom);
									mol.changeBond(bndCToO, Molecule.cBondTypeDouble);

								}
							}

						}
					}
					changed = true;
					break;
				}
			}

			mol.ensureHelperArrays(Molecule.cHelperRings);
		}
	}

	private void neutralizeNegative(StereoMolecule mol){
		
		for (int i = 0; i < arrNeutralizeNegative.length; i++) {

			StereoMolecule frag = arrNeutralizeNegative[i];
			
			SSSearcher sss = new SSSearcher();
			
			sss.setMol(frag, mol);
			
			if(sss.findFragmentInMolecule()>0){
				List<int []> liArrMatch = sss.getMatchList();
				
				if(liArrMatch == null){
					return;
				}
				
				for (int[] arrMatch : liArrMatch) {
					
					for (int j = 0; j < arrMatch.length; j++) {
						int atom = arrMatch[j];
						
						if(mol.getAtomicNo(atom)==6){
							mol.setAtomCharge(atom, 0);
						} else if(mol.getAtomicNo(atom)==8){
							mol.setAtomCharge(atom, 0);
						} else if(mol.getAtomicNo(atom)==16){
							mol.setAtomCharge(atom, 0);
						}
					}
				}
			}

		}
		mol.ensureHelperArrays(Molecule.cHelperRings);
	}
	
	/**
	 * -N#C, N positive charged, C negative charged.
	 * @param mol
	 */
	private void unifyIsoCyano(StereoMolecule mol){
		
		int tripleBondOrder = 3;
		
		StereoMolecule fragIsoCyano = arrUnifyIsoCyano[0];
		
		SSSearcher sss = new SSSearcher();
		
		sss.setMol(fragIsoCyano, mol);
		
		if(sss.findFragmentInMolecule()>0){
			List<int []> liArrMatch = sss.getMatchList();
			
			if(liArrMatch == null){
				return;
			}
			
			for (int[] arrMatch : liArrMatch) {
				
				for (int i = 0; i < arrMatch.length; i++) {
					int atom = arrMatch[i];
					
					if(mol.getAtomicNo(atom)==7){
						mol.setAtomCharge(atom, 1);
					} else if(mol.getAtomicNo(atom)==6){
						boolean tripleBond=false;
						
						for (int j = 0; j < mol.getConnAtoms(atom); j++) {
							if(mol.getConnBondOrder(atom, j)==tripleBondOrder){
								tripleBond = true;
								break;
							}
						}
						if(tripleBond && (mol.getConnAtoms(atom)==1))
							mol.setAtomCharge(atom, -1);
					}

				}
			}
			mol.ensureHelperArrays(Molecule.cHelperRings);
		}
	}
	
	private void unifyAzido(StereoMolecule mol){
		
		StereoMolecule fragAzido = arrUnifyAzido[0];
		
		SSSearcher sss = new SSSearcher();
		
		sss.setMol(fragAzido, mol);
		
		if(sss.findFragmentInMolecule()>0){
			List<int []> liArrMatch = sss.getMatchList();
			
			if(liArrMatch == null){
				return;
			}
			
			for (int[] arrMatch : liArrMatch) {
				
				for (int i = 0; i < arrMatch.length; i++) {
					int atom = arrMatch[i];
					
					if(mol.getAtomicNo(atom)==7){
						
						// Terminal N, negative charged.
						if(mol.getConnAtoms(atom)==1){
							mol.setAtomCharge(atom, -1);
							break;
						}
						
						boolean centerN=true;
						// Two connected N?-->Center N, positive charged.
						for (int j = 0; j < mol.getConnAtoms(atom); j++) {
							if(mol.getConnAtom(atom, j)!=7){
								centerN = false;
								break;
							}
						}
						if(centerN)
							mol.setAtomCharge(atom, 1);
					} 
				}
			}
			mol.ensureHelperArrays(Molecule.cHelperRings);
		}
	}
	
	/**
	 * Sets the charge to 0.
	 * @param mol
	 */
	private static void repairTertiaryNitrogen(StereoMolecule mol){
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i) == 7) {
				if(mol.getOccupiedValence(i)==3){
					mol.setAtomCharge(i, 0);
				}
			}
		}
		mol.ensureHelperArrays(Molecule.cHelperRings);
	}
	
	/**
	 * Sets the charge to +1.
	 * @param mol
	 */
	private static void repairQuaternaryNitrogen(StereoMolecule mol){
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i) == 7) {
				if(mol.getOccupiedValence(i)==4){
					if(mol.getAtomCharge(i)==0) {
						mol.setAtomCharge(i, 1);
					}
				}
			}
		}
		mol.ensureHelperArrays(Molecule.cHelperRings);
	}
	
	private void repairNitro(StereoMolecule mol){
		
		for (int i = 0; i < arrRepairNitro.length; i++) {
			
			SSSearcher sss = new SSSearcher();
			
			StereoMolecule frag = arrRepairNitro[i];
			
			sss.setMol(frag, mol);
			
			if(sss.findFragmentInMolecule()>0){
				List<int []> liArrMatch = sss.getMatchList();
				
				if(liArrMatch == null){
					return;
				}
				
				for (int[] arrMatch : liArrMatch) {
					for (int j = 0; j < arrMatch.length; j++) {
						int atom = arrMatch[j];
						
						if(mol.getAtomicNo(atom)==7){
							mol.setAtomCharge(atom, 1);
						}else if(mol.getAtomicNo(atom)==8){
							boolean doubleBondedOxy = false;
							for (int k = 0; k < mol.getConnAtoms(atom); k++) {
								if(mol.getConnBondOrder(atom, k)==Molecule.cBondTypeDouble){
									doubleBondedOxy = true;
									break;
								}
							}
							
							if(!doubleBondedOxy){
								mol.setAtomCharge(atom, -1);
							}
						}
					}
				}
			}
			mol.ensureHelperArrays(Molecule.cHelperRings);
		}
		
	}
	
	public static MoleculeStandardizer getInstance(){
		
		if(instance==null){
			instance = new MoleculeStandardizer();
		}
		
		return instance;
	}

}
