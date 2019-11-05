package com.actelion.research.chem.phesa.pharmacophore;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.actelion.research.chem.AtomFunctionAnalyzer;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;


public class IonizableGroupDetector {
	
	private StereoMolecule mol;
	private List<ArrayList<Integer>> ionizableGroups;
	private RingCollection ringCollection;
	
	public IonizableGroupDetector(StereoMolecule mol) {
		this.mol = mol;
		ionizableGroups = new ArrayList<ArrayList<Integer>>();
		ringCollection = mol.getRingSet();
	}
	
	public ArrayList<PPGaussian> detect() {
		ArrayList<PPGaussian> chargePoints = new ArrayList<PPGaussian>();
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
				chargePoints.add(new PPGaussian(7,cp));
			}
			
			
		}
		for(int a=0;a<mol.getAtoms();a++) {
			if(alreadyDetected(a)) continue;
			if(mol.getAtomicNo(a)==8) { //oxygen
				int aa = mol.getConnAtom(a,0);
				if(alreadyDetected(aa)) continue;
				if(AtomFunctionAnalyzer.isAcidicOxygen(mol, a)) { //COOH,SO3H,PO3H2, N(+)-OH
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
						chargePoints.add(new PPGaussian(6,cp));
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
						chargePoints.add(new PPGaussian(15,cp));
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
						chargePoints.add(new PPGaussian(16,cp));
						continue;					
					}
				}
				}
				else if(mol.getAtomicNo(a)==7) {
					
					if(mol.getConnAtoms(a)<=2) { //HNR2 or H2NR
							int nDBs = 0;
							boolean found=false;
							for(int i=0;i<mol.getConnAtoms(a) && !found;i++) { //search for amidine
								int aa = mol.getConnAtom(a, i);
								if(alreadyDetected(aa)) continue;
								if(mol.getAtomicNo(aa)==6) {
									if(mol.getBondOrder(mol.getBond(a, aa))==2)nDBs++;
									for(int j=0;j<mol.getConnAtoms(aa) && !found;j++) {
										int aaa = mol.getConnAtom(aa, j);
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
												chargePoints.add(new PPGaussian(6,cp));
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
						chargePoints.add(new PPGaussian(7,cp));
						continue;
					}
			}
			if(alreadyDetected(a))continue;
			else {
				int charge = mol.getAtomCharge(a);
				if(charge!=0 && !hasCounterChargedNeighbour(a)) {
					charge = charge>0 ? 1 : -1;
					ChargePoint cp = new ChargePoint(mol,a,new ArrayList<Integer>(),charge);
					chargePoints.add(new PPGaussian(mol.getAtomicNo(a),cp));
				}
			}
		}
		return chargePoints;
	}
		

	
	private boolean hasCounterChargedNeighbour(int a) {
		for(int aa=0;aa<mol.getConnAtoms(a);aa++) {
			if(mol.getAtomCharge(a)*mol.getAtomCharge(mol.getConnAtom(a,aa))<0)
				return true;
		}
		return false;
	}
	
	private boolean alreadyDetected(int a) {

		boolean isDetected = ionizableGroups.stream().flatMap(List::stream).collect(Collectors.toList()).contains(a) ? true : false;
		return isDetected;
	}


}
