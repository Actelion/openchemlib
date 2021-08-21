package com.actelion.research.chem.phesa.pharmacophore;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.phesa.pharmacophore.pp.ChargePoint;


/**
 * derives charged groups
 * @author wahljo1
 *
 */

public class ChargedGroupDetector {
	
	private StereoMolecule mol;
	private List<ArrayList<Integer>> chargedGroups;
	private RingCollection ringCollection;
	
	public ChargedGroupDetector(StereoMolecule mol) {
		this.mol = mol;
		chargedGroups = new ArrayList<ArrayList<Integer>>();
		ringCollection = mol.getRingSet();
	}
	
	public ArrayList<ChargePoint> detect() {
		ArrayList<ChargePoint> chargePoints = new ArrayList<ChargePoint>();
		//detect tetrazoles
		for(int r=0;r<ringCollection.getSize();r++) {
			ArrayList<Integer> tetrazole = new ArrayList<Integer>();
			int totCharge = 0;
			int[] ringAtoms = ringCollection.getRingAtoms(r);
			for(Integer atom : ringAtoms) {
				if(alreadyDetected(atom)) continue;
				if(mol.getAtomicNo(atom)==7 && mol.isAromaticAtom(atom) && mol.getConnAtoms(atom)<=2) {
					tetrazole.add(atom);
					totCharge += mol.getAtomCharge(atom);
				}
			}
			if(tetrazole.size()==4 && totCharge<0) {
				chargedGroups.add(tetrazole);
				ChargePoint cp = new ChargePoint(mol,tetrazole.get(0),Arrays.asList(tetrazole.get(1),
						tetrazole.get(2),tetrazole.get(3)),totCharge);
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
				ChargePoint cp = isPartOfChargedAcid(mol,a);
				if(cp==null)
					continue;
				chargePoints.add(cp);
				ArrayList<Integer> group = new ArrayList<>();
				group.add(cp.getChargeAtom());
				group.addAll(cp.getNeighbours());
				chargedGroups.add(group);
				}
			else if(mol.getAtomicNo(a)==7) {
				int carbonCenter = -1;
				int charge = 0;
				List<Integer> neighbours = new ArrayList<>();
				if(!mol.isAromaticAtom(a) && mol.getConnAtoms(a)<=2) { //HNR2 or H2NR
					neighbours.add(a);
					charge+=mol.getAtomCharge(a);
					int nDBs = 0;
					for(int i=0;i<mol.getConnAtoms(a);i++) { //search for amidine
						int aa = mol.getConnAtom(a, i);
						if(alreadyDetected(aa)) continue;
						if(mol.getAtomicNo(aa)==6) {
							carbonCenter = aa;
							charge+=mol.getAtomCharge(aa);
							if(mol.getBondOrder(mol.getBond(a, aa))==2)nDBs++;
							for(int j=0;j<mol.getConnAtoms(aa);j++) {
								int aaa = mol.getConnAtom(aa, j);
								if(mol.isAromaticAtom(aaa))
									continue;
								if(aaa==a) continue;
								if(alreadyDetected(aaa)) continue;
								if(mol.getAtomicNo(aaa)==7 && mol.getConnAtoms(aaa)<=2) {
									neighbours.add(aaa);
									charge+=mol.getAtomCharge(aaa);
									if(mol.getBondOrder(mol.getBond(aa, aaa))==2)
										nDBs++;

									}
										
								}
							
						}
					}
					if(nDBs>1 && charge>0) {
						ChargePoint cp = new ChargePoint(mol,carbonCenter,neighbours,charge);
						chargePoints.add(cp);
						ArrayList<Integer> group = new ArrayList<>();
						group.add(cp.getChargeAtom());
						group.addAll(cp.getNeighbours());
						chargedGroups.add(group);
						
				}
				
			}

		}
		if(alreadyDetected(a))continue;
		else {
			int charge = mol.getAtomCharge(a);
			if(charge!=0 && !hasCounterChargedNeighbour(a)) {
				charge = charge>0 ? 1 : -1;
				ChargePoint cp = new ChargePoint(mol,a,new ArrayList<Integer>(),charge);
				chargePoints.add(cp);
				ArrayList<Integer> group = new ArrayList<>();
				group.add(a);
				chargedGroups.add(group);
			}
			}
		}
		return chargePoints;
	}

	public List<ArrayList<Integer>> getChargedGroups() {
		return chargedGroups;
	}

	private boolean hasCounterChargedNeighbour(int a) {
		for(int aa=0;aa<mol.getConnAtoms(a);aa++) {
			if(mol.getAtomCharge(a)*mol.getAtomCharge(mol.getConnAtom(a,aa))<0)
				return true;
		}
		return false;
	}
	
	private boolean alreadyDetected(int a) {

		boolean isDetected = chargedGroups.stream().flatMap(List::stream).collect(Collectors.toList()).contains(a) ? true : false;
		return isDetected;
	}
	
	/**
	 * independent if acid is protonated or not
	 * @param mol
	 * @param atom
	 * @return
	 */
	private static ChargePoint isPartOfChargedAcid(StereoMolecule mol, int atom) {
		if (mol.getAtomicNo(atom) != 8
		 || mol.getConnAtoms(atom) != 1
		 || mol.getConnBondOrder(atom, 0) != 1)
			return null;

		int connAtom = mol.getConnAtom(atom, 0);
		
		// COO-
		if(mol.getAtomicNo(connAtom)==6){
			List<Integer> neighbours = new ArrayList<>();
			neighbours.add(atom);
			int charge = 0;
			charge+=mol.getAtomCharge(atom);
			charge+=mol.getAtomCharge(connAtom);
			int nConnected2C = mol.getConnAtoms(connAtom);
			for (int i = 0; i < nConnected2C; i++) {
				int indexAtom = mol.getConnAtom(connAtom, i);
				
				if(indexAtom==atom){
					continue;
				}
				
				if(mol.getAtomicNo(indexAtom) != 8){
					continue;
				}
				
				charge+=mol.getAtomCharge(indexAtom);
				neighbours.add(indexAtom);
			}
			if(charge<0) {
				ChargePoint cp = new ChargePoint(mol,connAtom,neighbours,charge);
				return cp;
			
			}

		} else if (mol.getAtomicNo(connAtom) == 16) { // CSOOOH
			int nConnected2S = mol.getConnAtoms(connAtom);
			List<Integer> neighbours = new ArrayList<>();
			neighbours.add(atom);
			int charge = 0;
			charge+=mol.getAtomCharge(atom);
			charge+=mol.getAtomCharge(connAtom);
			for (int i = 0; i < nConnected2S; i++) {
				int indexAtom = mol.getConnAtom(connAtom, i);
				
				if(indexAtom==atom)
					continue;

				if(mol.getAtomicNo(indexAtom) != 8)
					continue;
				charge+=mol.getAtomCharge(indexAtom);
				neighbours.add(indexAtom);
			}
			if(charge<0) {
				ChargePoint cp = new ChargePoint(mol,connAtom,neighbours,charge);
				return cp;
			}
			
		
		} else if(mol.getAtomicNo(connAtom)==15){ // CP=O(OH)(OH)
			int nConnected2P = mol.getConnAtoms(connAtom);
			List<Integer> neighbours = new ArrayList<>();
			neighbours.add(atom);
			int charge = 0;
			charge+=mol.getAtomCharge(atom);
			charge+=mol.getAtomCharge(connAtom);
			for (int i = 0; i < nConnected2P; i++) {
				int indexAtom = mol.getConnAtom(connAtom, i);
				if(indexAtom==atom)
					continue;

				if(mol.getAtomicNo(indexAtom) != 8)
					continue;
				charge+=mol.getAtomCharge(indexAtom);
				neighbours.add(indexAtom);
			}
			if(charge<0) {
				ChargePoint cp = new ChargePoint(mol,connAtom,neighbours,charge);
				return cp;
			}
		}
		return null;

	}


}
