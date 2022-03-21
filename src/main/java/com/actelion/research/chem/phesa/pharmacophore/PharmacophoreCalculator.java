package com.actelion.research.chem.phesa.pharmacophore;

import com.actelion.research.chem.AtomFunctionAnalyzer;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.interactionstatistics.InteractionAtomTypeCalculator;
import com.actelion.research.chem.phesa.pharmacophore.pp.AcceptorPoint;
import com.actelion.research.chem.phesa.pharmacophore.pp.AromRingPoint;
import com.actelion.research.chem.phesa.pharmacophore.pp.DonorPoint;
import com.actelion.research.chem.phesa.pharmacophore.pp.IPharmacophorePoint;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class PharmacophoreCalculator {
	
	public static final int ACCEPTOR_ID = 0;
	public static final int DONOR_ID = 1;
	public static final int CHARGE_NEG_ID = 2;
	public static final int CHARGE_POS_ID = 3;
	public static final int AROM_ID = 4;
	public static final int LIPO_ID = 5;
	public static final int AROM_RING_ID = 6;
	public static final int EXIT_VECTOR_ID = 7;
	public static final int MAX_ID = 7;
	
	private PharmacophoreCalculator() {}
	
	
	public static List<IPharmacophorePoint> getPharmacophorePoints(StereoMolecule mol) {
		List<IPharmacophorePoint> ppPoints = new ArrayList<IPharmacophorePoint>();
		RingCollection rc = mol.getRingSet();
		for(int r=0;r<rc.getSize();r++) {
			if(!rc.isAromatic(r))
				continue;
			int[] ringAtoms = rc.getRingAtoms(r);
			AromRingPoint arp = new AromRingPoint(mol,ringAtoms[0],Arrays.stream(ringAtoms)
				      .boxed()
				      .collect(Collectors.toList()));
			ppPoints.add(arp);
		}
		for(int i=0;i<mol.getAllAtoms();i++) {
			/*
			if(mol.getAtomicNo(i)==0) { //end point of exit vector
				if(mol.getConnAtoms(i)==0)
					continue;
				for(int c=0;c<mol.getConnAtoms(i);c++) {
					int j = mol.getConnAtom(i, 0);
					ExitVectorPoint evp = new ExitVectorPoint(mol,j,i);
					ppPoints.add(evp);
					}
				}
				*/
			
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
					if(mol.getAtomicNo(i)==8 && neighbours==1 && (mol.getConnBondOrder(i, 0)==2 || AtomFunctionAnalyzer.isAcidicOxygen(mol, i) || mol.getAtomCharge(i)==-1 )) {
						int a1 = mol.getConnAtom(i,0);
						if (mol.getConnAtoms(a1) > 1) { // added this check to prevent OOB exceptions with cropped proteins; TLS 17Oct2021
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
					}
					else if (neighbourList.size() != 0) { // added this check to prevent OOB exceptions with cropped proteins; TLS 17Oct2021
						AcceptorPoint ap = new AcceptorPoint(mol,i,neighbourList,interactionClass);
						ppPoints.add(ap);
					}
				}
			}
		}
		return ppPoints;
	}
	
	public static boolean isAcceptor(StereoMolecule mol, int a) {
		if (mol.getAtomicNo(a)==7 || mol.getAtomicNo(a)==8) {
		if (mol.getAtomCharge(a)<=0) { //charge is positive -> no acceptor
			if (mol.isAromaticAtom(a)) { 
				if (mol.getAllConnAtoms(a)<3) {
					return true;
				}
				else {
					return false; //is in aromatic ring and has at least 3 bonded atoms -> no acceptor 
				}
			}
			else if (mol.getAtomicNo(a)==7){ // atom is not aromatic
				if(mol.getConnAtoms(a) == 1 && mol.getConnBondOrder(a, 0)==3) //nitrile
					return true;
				if (mol.isFlatNitrogen(a)) { 
					for(int b=0;b<mol.getConnAtoms(a);b++) {
						if(mol.getConnBondOrder(a, b)>1)
							return true;
					}
					return false;
					
				}
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
		else {
			return false;
		}
	}
	
	public static boolean isDonorHydrogen(StereoMolecule mol, int h) {
		if(mol.getAtomicNo(h)==1) {
			int dh = mol.getConnAtom(h, 0);
			if (mol.getAtomCharge(dh)>=0 && (mol.getAtomicNo(dh)==7 || mol.getAtomicNo(dh)==8) ) { 
				return true;
			}
			else return false;
		}
		else
			return false;
	}
	
	public static boolean isDonorHeavyAtom(StereoMolecule mol, int d) {
		boolean isDonor = false;
		if (mol.getAtomicNo(d)==7 || mol.getAtomicNo(d)==8) {
			if (mol.getAtomCharge(d)>=0 ) {
				if(mol.getAllHydrogens(d)>0)
				isDonor = true;
			}
		}
		return isDonor;
	}
	

}
