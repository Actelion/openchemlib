package com.actelion.research.chem.phesa.pharmacophore;

import com.actelion.research.chem.AtomFunctionAnalyzer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.phesa.pharmacophore.pp.ChargePoint;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;


/**
 * derives potentially ionizable Groups, independent of the chosen explicit protonation state
 * @author wahljo1
 *
 */

public class IonizableGroupDetector {
	
	private StereoMolecule mol;
	private List<ArrayList<Integer>> ionizableGroups;
	private RingCollection ringCollection;
	
	public IonizableGroupDetector(StereoMolecule mol) {
		this.mol = mol;
		ionizableGroups = new ArrayList<ArrayList<Integer>>();
		ringCollection = mol.getRingSet();
	}
	
	public ArrayList<ChargePoint> detect() {
		ArrayList<ChargePoint> chargePoints = new ArrayList<ChargePoint>();
		ArrayList<Integer> ionizableGroup;
		//detect tetrazoles
		for(int r=0;r<ringCollection.getSize();r++) {
			ArrayList<Integer> tetrazole = new ArrayList<Integer>();
			int[] ringAtoms = ringCollection.getRingAtoms(r);
			for(Integer atom : ringAtoms) {
				if(alreadyDetected(atom)) continue;
				if(mol.getAtomicNo(atom)==7 && mol.isAromaticAtom(atom) && mol.getConnAtoms(atom)<=2)
					tetrazole.add(atom);
			}
			if(tetrazole.size()==4) {
				ionizableGroups.add(tetrazole);
				ChargePoint cp = new ChargePoint(mol,tetrazole.get(0),Arrays.asList(tetrazole.get(1),
						tetrazole.get(2),tetrazole.get(3)),-1);
				chargePoints.add(cp);
			}
		}
		for(int a=0;a<mol.getAtoms();a++) {
			if(alreadyDetected(a)) continue;
			if(mol.getAtomicNo(a)==8) { //oxygen
				if(mol.getConnAtoms(a)==0)
					continue;
				int aa = mol.getConnAtom(a,0);
				if(alreadyDetected(aa)) continue;
				if(AtomFunctionAnalyzer.isAcidicOxygen(mol, a, false)) { //COOH,SO3H,PO3H2, N(+)-OH
					if(mol.getAtomicNo(aa)==6) { //COOH
						ionizableGroup = new ArrayList<Integer>();
						ionizableGroup.add(a);
						ionizableGroup.add(aa);
						int aaa1 = mol.getConnAtom(aa, 0);
						int aaa2 = mol.getConnAtom(aa, 1);
						int aaa3 = mol.getConnAtom(aa, 2);
						int aaa = (aaa1!=a && mol.getAtomicNo(aaa1)==8) ? aaa1 : 
							(aaa2!=a && mol.getAtomicNo(aaa2)==8) ? aaa2 : 
								aaa3;
						if(alreadyDetected(aaa)) continue;
						ionizableGroup.add(aaa);
						ionizableGroups.add(ionizableGroup);
						ChargePoint cp = new ChargePoint(mol,aa,new ArrayList<Integer>(),-1);
						chargePoints.add(cp);
						continue;
					}
					else if (mol.getAtomicNo(aa)==15) {//POO3H2
						ionizableGroup = new ArrayList<Integer>();
						ionizableGroup.add(a);
						ionizableGroup.add(aa);
						for(int i=0;i<mol.getConnAtoms(aa);i++) {
							int aaa = mol.getConnAtom(aa, i);
							if(mol.getAtomicNo(aaa) ==8 && aaa!=a) {
								if(alreadyDetected(aaa)) continue;
								ionizableGroup.add(aaa);
							}
						}
						ionizableGroups.add(ionizableGroup);
						ChargePoint cp = new ChargePoint(mol,aa,new ArrayList<Integer>(),-1);
						chargePoints.add(cp);
						continue;					
					}
					else if (mol.getAtomicNo(aa)==16) {//SOO3H
						ionizableGroup = new ArrayList<Integer>();
						ionizableGroup.add(a);
						ionizableGroup.add(aa);
						for(int i=0;i<mol.getConnAtoms(aa);i++) {
							int aaa = mol.getConnAtom(aa, i);
							if(mol.getAtomicNo(aaa) ==8 && aaa!=a) {
								if(alreadyDetected(aaa)) continue;
								ionizableGroup.add(aaa);
							}
							
						}
						ionizableGroups.add(ionizableGroup);
						ChargePoint cp = new ChargePoint(mol,aa,new ArrayList<Integer>(),-1);
						chargePoints.add(cp);
						continue;					
					}
				}
				}
				else if(mol.getAtomicNo(a)==7) {
					if(!mol.isAromaticAtom(a) && mol.getConnAtoms(a)<=2) { //HNR2 or H2NR
							boolean found=false;
							for(int i=0;i<mol.getConnAtoms(a) && !found;i++) { //search for amidine
								int nDBs = 0;
								int aa = mol.getConnAtom(a, i);
								if(alreadyDetected(aa)) continue;
								if(mol.getAtomicNo(aa)==6) {
									if(mol.getBondOrder(mol.getBond(a, aa))==2)nDBs++;
									for(int j=0;j<mol.getConnAtoms(aa) && !found;j++) {
										int aaa = mol.getConnAtom(aa, j);
										if(mol.isAromaticAtom(aaa))
											continue;
										if(aaa==a) continue;
										if(alreadyDetected(aaa)) continue;
										if(mol.getAtomicNo(aaa)==7 && mol.getConnAtoms(aaa)<=2) {;
											if(mol.getBondOrder(mol.getBond(aa, aaa))==2)nDBs++;
											if(nDBs==1) { //Amidine
												ionizableGroup = new ArrayList<Integer>();
												ionizableGroup.add(a);
												ionizableGroup.add(aa);
												ionizableGroup.add(aaa);
												ionizableGroups.add(ionizableGroup);
												ChargePoint cp = new ChargePoint(mol,aa,new ArrayList<Integer>(),1);
												chargePoints.add(cp);
												found = true;

											}
												
										}
									
								}
							}
							
					}
					
				}
					if(alreadyDetected(a))continue;
					if(AtomFunctionAnalyzer.isBasicNitrogen(mol, a)) {
						ionizableGroup = new ArrayList<Integer>();
						ionizableGroup.add(a);
						ionizableGroups.add(ionizableGroup);
						ChargePoint cp = new ChargePoint(mol,a,new ArrayList<Integer>(),1);
						chargePoints.add(cp);
						continue;
					}
			}
			if(alreadyDetected(a))continue;
			else {
				int charge = mol.getAtomCharge(a);
				if(charge!=0 && !hasCounterChargedNeighbour(a)) {
					charge = charge>0 ? 1 : -1;
					ChargePoint cp = new ChargePoint(mol,a,new ArrayList<Integer>(),charge);
					chargePoints.add(cp);
				}
			}
		}
		return chargePoints;
	}

	public List<ArrayList<Integer>> getIonizableGroups() {
		return ionizableGroups;
	}

	private boolean hasCounterChargedNeighbour(int a) {
		for(int aa=0;aa<mol.getConnAtoms(a);aa++) {
			if(mol.getAtomCharge(a)*mol.getAtomCharge(mol.getConnAtom(a,aa))<0)
				return true;
		}
		return false;
	}
	
	private boolean alreadyDetected(int a) {
		return ionizableGroups.stream().flatMap(List::stream).collect(Collectors.toList()).contains(a) ? true : false;
	}
	
	/**
	 * independent if acid is protonated or not
	 * @param mol
	 * @param atom
	 * @return
	 */
	private static boolean isPartOfAcid(StereoMolecule mol, int atom) {
		if (mol.getAtomicNo(atom) != 8
		 || mol.getConnAtoms(atom) != 1
		 || mol.getConnBondOrder(atom, 0) != 1)
			return false;

		int connAtom = mol.getConnAtom(atom, 0);

		// COOH
		if(mol.getAtomicNo(connAtom)==6){
			int nConnected2C = mol.getConnAtoms(connAtom);
			for (int i = 0; i < nConnected2C; i++) {
				int indexAtom = mol.getConnAtom(connAtom, i);
				
				if(indexAtom==atom){
					continue;
				}
				
				if(mol.getAtomicNo(indexAtom) != 8){
					continue;
				}
				
				int indexBond = mol.getBond(connAtom, indexAtom);
				
				if(mol.getBondType(indexBond)==Molecule.cBondTypeDouble)
					return true;
			}
		} else if (mol.getAtomicNo(connAtom) == 7) {
			if (mol.getAtomCharge(connAtom) == 1) // (N+)-OH
				return true;
		} else if (mol.getAtomicNo(connAtom) == 16) { // CSOOOH
			int nConnected2S = mol.getConnAtoms(connAtom);
			
			int nDoubleBondedO2S=0;
			for (int i = 0; i < nConnected2S; i++) {
				int indexAtom = mol.getConnAtom(connAtom, i);
				
				if(indexAtom==atom)
					continue;

				if(mol.getAtomicNo(indexAtom) != 8)
					continue;

				int indexBond = mol.getBond(connAtom, indexAtom);
				
				if(mol.getBondType(indexBond)==Molecule.cBondTypeDouble)
					nDoubleBondedO2S++;
			}
			
			if(nDoubleBondedO2S == 2)
				return true;
		} else if(AtomFunctionAnalyzer.isAcidicOxygenAtPhosphoricAcid(mol, atom)) // CP=O(OH)(OH)
			return true;

		return false;
	}


}
