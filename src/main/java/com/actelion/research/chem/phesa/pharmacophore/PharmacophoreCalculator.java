package com.actelion.research.chem.phesa.pharmacophore;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.AtomFunctionAnalyzer;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;

public class PharmacophoreCalculator {
	
	private PharmacophoreCalculator() {}
	
	public static List<IPharmacophorePoint> getPharmacophorePoints(StereoMolecule mol) {
		List<IPharmacophorePoint> ppPoints = new ArrayList<IPharmacophorePoint>();
		IonizableGroupDetector detector = new IonizableGroupDetector(mol);
		ppPoints.addAll(detector.detect());
		for(int i=0;i<mol.getAllAtoms();i++) {
			if (mol.getAtomicNo(i)==1) {
				if(isDonorHydrogen(mol,i)) {
					int d = mol.getConnAtom(i,0);
					int interactionClass = InteractionAtomTypeCalculator.getAtomType(mol, d);
					if(interactionClass<0) {
						continue;
					}
					DonorPoint dp = new DonorPoint(mol,d,i,interactionClass);
					ppPoints.add(dp);
				}
			}
			else if (mol.getAtomicNo(i)==7 || mol.getAtomicNo(i)==8) {
				if(isAcceptor(mol,i)) {
					int neighbours = mol.getAllConnAtoms(i);
					List<Integer> neighbourList= new ArrayList<Integer>();
					for(int j=0;j<neighbours;j++) 
						neighbourList.add(mol.getConnAtom(i,j));
					
					int interactionClass = InteractionAtomTypeCalculator.getAtomType(mol, i);
					if(interactionClass<0) {
						continue;
					}
					if(mol.getAtomicNo(i)==8 && neighbours==1 && (mol.getConnBondOrder(i, 0)==2 || AtomFunctionAnalyzer.isAcidicOxygen(mol, i) )) {
						int a1 = mol.getConnAtom(i,0);
						if(!(mol.getAtomicNo(a1)==16 || mol.getAtomicNo(a1)==15)) {		
							int aa1 = mol.getConnAtom(a1,0);
							if(aa1==i) 
								aa1 = mol.getConnAtom(a1,1);
							neighbourList.add(aa1);
							AcceptorPoint ap = new AcceptorPoint(mol,i,neighbourList,interactionClass,1);
							ppPoints.add(ap);
							List<Integer> neighbourList2 = new ArrayList<Integer>();
							for(int neighbour : neighbourList) {
								neighbourList2.add(neighbour);
							}
							ap = new AcceptorPoint(mol,i,neighbourList2,interactionClass,2);
							ppPoints.add(ap);
							}
						else { 
							AcceptorPoint ap = new AcceptorPoint(mol,i,neighbourList,interactionClass);
							ppPoints.add(ap);
						}
	
					}
					else {
					AcceptorPoint ap = new AcceptorPoint(mol,i,neighbourList,interactionClass);
					ppPoints.add(ap);
					}
			}
		}
		}
		return ppPoints;
	}
	
	private static boolean isAcceptor(StereoMolecule mol, int a) {
		if (mol.getAtomCharge(a)<=0) { //charge is not positive -> no acceptor
			if (mol.isAromaticAtom(a)) { 
				if (mol.getAllConnAtoms(a)<3) {
					return true;
				}
				else {
					return false; //is in aromatic ring and has at least 3 bonded atoms -> no acceptor 
				}
			}
			else if (mol.getAtomicNo(a)==7){ // atom is not aromatic
				if(mol.getConnBondOrder(a, 0)==3) //nitrile
					return true;
				if (mol.isFlatNitrogen(a)) 
					return false;
				for(int i=0;i<mol.getAllConnAtoms(a);i++) {
					int aa = mol.getConnAtom(a, i);
					if (mol.getAtomicNo(aa)==6) {
						for(int j=0;j<mol.getAllConnAtoms(aa);j++) {
							int aaa = mol.getConnAtom(aa,j);
							if(a==aaa) continue;
							if (mol.getBondOrder(mol.getBond(aa,aaa))==2) { 
								if (mol.getAtomicNo(aaa)==7) return false; //amide structure
								if (mol.getAtomicNo(aaa)==8) return false;
								if (mol.getAtomicNo(aaa)==16) return false;
							}
						}
					}

				}
			}

		return true;
		}
		else return false;
	}
	
	private static boolean isDonorHydrogen(StereoMolecule mol, int h) {
		int dh = mol.getConnAtom(h, 0);
		if (mol.getAtomCharge(dh)>=0 && (mol.getAtomicNo(dh)==7 || mol.getAtomicNo(dh)==8) ) { //charge is not positive -> no acceptor
			return true;
		}
		else return false;
	}
	

}
