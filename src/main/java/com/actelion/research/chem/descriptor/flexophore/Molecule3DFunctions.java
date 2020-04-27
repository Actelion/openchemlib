/*
 * Copyright (c) 2020.
 * Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
 *
 *  This file is part of DataWarrior.
 *
 *  DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version.
 *
 *  DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with DataWarrior.
 *  If not, see http://www.gnu.org/licenses/.
 *
 *  @author Modest v. Korff
 *
 */

package com.actelion.research.chem.descriptor.flexophore;

import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.chem.*;
import com.actelion.research.chem.descriptor.flexophore.calculator.StructureCalculator;
import com.actelion.research.chem.mcs.ListWithIntVec;
import com.actelion.research.util.datamodel.IntegerDouble;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeNode;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

public class Molecule3DFunctions {

	public static final String LEADING_SPACE_DESCRIPTION = "     ";

	public static final int FLAG_CENTER_ATOM = 1<<4;
	
	public static final int FLAG_AROMATIC_ATOM = 1<<5;

	
	private static final NumberFormat NF = new DecimalFormat("0.000");
	
	private static final int STERIC_HINDRANCE_STRONG = 1000;
	
	private static final int STERIC_HINDRANCE_MEDIUM = 100;
	
	private static final int STERIC_HINDRANCE_WEAK = 10;

	private static final double STERIC_HINDRANCE_FACTOR_STRONG = 1.0;
	
	private static final double STERIC_HINDRANCE_FACTOR_MEDIUM = 1.5;
	
	private static final double STERIC_HINDRANCE_FACTOR_WEAK = 2.0;


	
	public static double [][] getAdjacencyMatrix(Molecule3D mol) {
		
		double arr[][]= new double[mol.getAllAtoms()][mol.getAllAtoms()];
		
		for (int i = 0; i < mol.getAllBonds(); i++) {
			
			int indexAtom1 = mol.getBondAtom(0, i);
			
			int indexAtom2 = mol.getBondAtom(1, i);
			
			arr[indexAtom1][indexAtom2]=1;
			arr[indexAtom2][indexAtom1]=1;
			
		}
		
		return arr;
	}

	
	/**
	 * Protons attached to an O.
	 * @param mol
	 * @return
	 */
	public static List<Integer> getAcidicProtons(Molecule3D mol){
		
		List<Integer> liIndexAtomAcidicProt = new ArrayList<Integer>();
		
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i)==8){
				
				int nConnected = mol.getAllConnAtoms(i);
				
				for (int j = 0; j < nConnected; j++) {
					
					int indexConnAt = mol.getConnAtom(i, j);
					
					if(mol.getAtomicNo(indexConnAt)==1){
						liIndexAtomAcidicProt.add(indexConnAt);
					}
				}
			}
		}
		
		return liIndexAtomAcidicProt;

	}
	
	/**
	 * Donor atoms. Oxygen with a connected hydrogen.
	 * @param mol
	 * @return
	 */
	public static List<Integer> getDeprotonableAtoms(Molecule3D mol){
		
		List<Integer> liIndexAtomAcidic = new ArrayList<Integer>();
		
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i)==8){
				
				int nConnected = mol.getAllConnAtoms(i);
				
				for (int j = 0; j < nConnected; j++) {
					
					int indexConnAt = mol.getConnAtom(i, j);
					
					if(mol.getAtomicNo(indexConnAt)==1){
						liIndexAtomAcidic.add(i);
						break;
					}
				}
			}
		}
		
		return liIndexAtomAcidic;

	}
	
	public static List<Integer> getProtonableAtoms(Molecule3D mol){
		
		List<Integer> liIndexAtomBasic = new ArrayList<Integer>();
		
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			
			if(mol.getAtomicNo(i)==7){
				
				if(isProtonableNitrogen(mol, i)) {
					liIndexAtomBasic.add(i);
					break;
				}
			}
		}
		
		return liIndexAtomBasic;

	}
	
	public static boolean isProtonableNitrogen(Molecule3D mol, int indexAtom){
		boolean protonableN=true;
		
		if(mol.getAtomicNo(indexAtom)!=7){
			return false;
		}
		
		int nConnected = mol.getAllConnAtoms(indexAtom);
		
			
		// Already protonated N?
		boolean hydrogen=false;
		for (int i = 0; i < nConnected; i++) {
			int indexAtomConn = mol.getConnAtom(indexAtom, i);
			
			if(mol.getAtomicNo(indexAtomConn)==1){
				hydrogen=true;
				return true;
			}
		}
		
		if((!hydrogen) && (nConnected==4)){
			return false;
		}
		
		
		int sumBondOrder=0;
		for (int i = 0; i < nConnected; i++) {
			sumBondOrder += mol.getConnBondOrder(indexAtom, i);
		}
		
		if(sumBondOrder==4){
			return false;
		}
		
		return protonableN;
	}
	

	public static boolean has2DCoordinates(Molecule3D mol){
		
		double sumX=0;
		double sumY=0;
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			Coordinates c = mol.getCoordinates(i);
			
			if(c==null){
				continue;
			}
			
			sumX += Math.abs(c.x);
			sumY += Math.abs(c.y);
		}
		
		if(sumX!=0 && sumY!=0){
			return true;
		}
		
		return false;
	}
	
	public static boolean has3DCoordinates(Molecule3D mol){
		
		double sumX=0;
		double sumY=0;
		double sumZ=0;
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			Coordinates c = mol.getCoordinates(i);
			
			if(c==null){
				continue;
			}
			
			sumX += Math.abs(c.x);
			sumY += Math.abs(c.y);
			sumZ += Math.abs(c.z);
		}
		
		if(sumX!=0 && sumY!=0){
			return true;
		}
		
		return false;
	}
	
	public static boolean hasChargedAtom(Molecule3D ff){
		
		boolean charged = false;
		
		for (int i = 0; i < ff.getAllAtoms(); i++) {
			
			if(ff.getAtomCharge(i)!=0){
				charged = true;
				break;
			}
		}
		
		return charged;
	}
	
	/**
	 * Calculates a value for the steric hindrance of the atoms from liIndexAtoms to all other atoms.
	 * @param molecule3D
	 * @param liIndexAtoms
	 * @return
	 */
	public static double getStericHindrance(Molecule3D molecule3D, List<Integer> liIndexAtoms) {
		double hindrance = 0;
		
		HashSet<Integer> hsIndexAtoms = new HashSet<Integer>(liIndexAtoms);
		
		List<Integer> liIndexAtoms2 = new ArrayList<Integer>();
		for (int i = 0; i < molecule3D.getAllAtoms(); i++) {
			if(!hsIndexAtoms.contains(i)){
				liIndexAtoms2.add(i);
			}
		}
		
    	double [][] arrDist = Molecule3DFunctions.getDistanceArray(molecule3D);

		
		for (int i = 0; i < liIndexAtoms.size(); i++) {
			int indexAt1 = liIndexAtoms.get(i);
			
			for (int j = 0; j < liIndexAtoms2.size(); j++) {
				int indexAt2 = liIndexAtoms2.get(i);
				
				double dist = arrDist[indexAt1][indexAt2];
				
				Element el1 = PeriodicTable.getElement(molecule3D.getAtomicNo(indexAt1));
				Element el2 = PeriodicTable.getElement(molecule3D.getAtomicNo(indexAt2));
				
				double r1 = el1.getVDWRadius();
				
				double r2 = el2.getVDWRadius();
				
				double sumRadii = r1+r2;
				
				double distHindranceStrong = dist * STERIC_HINDRANCE_FACTOR_STRONG;
				
				double distHindranceMedium = dist * STERIC_HINDRANCE_FACTOR_MEDIUM;
				
				double distHindranceWeak = dist * STERIC_HINDRANCE_FACTOR_WEAK;
				
				if(sumRadii>distHindranceStrong ) {
					hindrance += STERIC_HINDRANCE_STRONG;
				} else if(sumRadii>distHindranceMedium ) {
					hindrance += STERIC_HINDRANCE_MEDIUM;
				} else if(sumRadii>distHindranceWeak ) {
					hindrance += STERIC_HINDRANCE_WEAK;
				}
			}
		}
		
		
		return hindrance;
	}

	public static Molecule3D add(Molecule3D mol1, Molecule3D mol2) {
		Molecule3D ff = new Molecule3D(mol1);
		
		
		int [] arrMap = new int [mol2.getAllAtoms()];
		for (int i = 0; i < mol2.getAllAtoms(); i++) {
			int indexAt = ff.addAtom(mol2, i);
			
			arrMap[i]=indexAt;
		}
		
		for (int i = 0; i < mol2.getAllBonds(); i++) {
			int indexAt1 = mol2.getBondAtom(0, i);
			int indexAt2 = mol2.getBondAtom(1, i);
			
			int order = mol2.getBondOrder(i);
			
			int indexAt1New = arrMap[indexAt1];
			int indexAt2New = arrMap[indexAt2];
			
			ff.addBond(indexAt1New, indexAt2New, order);
			
		}
		return ff;
	}

	
	public static boolean isPyramedal(Molecule3D molecule3D){
		
		boolean pyramedal=false;
		if(molecule3D.getAllAtoms()==4) {
			for (int i = 0; i < molecule3D.getAllAtoms(); i++) {
				if(molecule3D.getAllConnAtoms(i)==3){
					pyramedal=true;
					break;
				}
			}
		}
		return pyramedal;
	}
	
	public static boolean isTripleBond(Molecule3D molecule3D){
		
		boolean tripleBond=false;
		for (int i = 0; i < molecule3D.getAllBonds(); i++) {
			
			if(molecule3D.getBondOrder(i)==3){
				tripleBond=true;
				break;
			}
		}
		
		return tripleBond;
	}
	
	public static void addMissingHydrogens(Molecule3D molecule3D) {
		for (int i = 0; i < molecule3D.getAtoms(); i++) {
			if(molecule3D.getAtomicNo(i)==8){
				int indNew = molecule3D.addAtom(1);
				molecule3D.addBond(i, indNew, 1);
			}
		}
	}

	public static void writeCharge2Label(Molecule3D molecule3D){
		
		for (int i = 0; i < molecule3D.getAllAtoms(); i++) {
			int charge = molecule3D.getAtomCharge(i);
			if(charge!=0){
				String l = "   " + charge;
				molecule3D.setAtomDescription(i, l);
			} else {
				molecule3D.setAtomDescription(i, "");
			}
		}
		
	}

	private static List<int []> getListFromSSSearcher(StereoMolecule mol, StereoMolecule fragment){
		SSSearcher sss = new SSSearcher();
		sss.setMol(fragment,mol);
		
		int numFrags = sss.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cMatchAromDBondToDelocalized);
		List<int []> liAtomLists = new ArrayList<int []>();
		if(numFrags > 0)  {
			liAtomLists = sss.getMatchList();
			// Remove double atomlists
			ArrayUtilsCalc.removeDoubletsIntOrderIndepend(liAtomLists);
		}

		return liAtomLists;
	}



	/**
	 * 
	 * @param molecule3D Molecule
	 * @param idCodeQ1 The idcode for the substructure (ode fragment)
	 * @param arrIndexAts2Flag the index for atoms in the idCode1q. The 
	 * corresponding atoms in the molecule will be deleted. 
	 */
	public static void flagSubstructure(Molecule3D molecule3D, String idCodeQ1, int [] arrIndexAts2Flag) {

		StereoMolecule query = new StereoMolecule();
		
		query.setFragment(true);
		
		new IDCodeParser().parse(query,idCodeQ1);
			
		SSSearcher sss = new SSSearcher();
		sss.setMol(query, molecule3D);
		
		int numFrags = sss.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cMatchAromDBondToDelocalized);
		List<int[]> liAtomLists = new ArrayList<int[]>();
		if(numFrags > 0)  {
			List<int []> vecMatchList = sss.getMatchList();
			for (Iterator<int []> iter = vecMatchList.iterator(); iter.hasNext();) {
				int [] atomlist = iter.next();
				liAtomLists.add(atomlist);
			}
		}
		
		if(liAtomLists.size() == 0) {
			// Exception ex = new Exception("Substructure not found.");
			// ex.printStackTrace();
			return;
		}
		
		// System.out.println("Found: " + liAtomLists.size());

		
		int [] arrAts2DelInMol = new int [arrIndexAts2Flag.length * liAtomLists.size()];
		int cc = 0;
		// Extract substructure from molecule
		for (Iterator<int []> iter = liAtomLists.iterator(); iter.hasNext();) {
			int [] arrAtomList = iter.next();
			
			for (int i = 0; i < arrIndexAts2Flag.length; i++) {
				arrAts2DelInMol[cc++] =  arrAtomList[arrIndexAts2Flag[i]];
			}
		}

		
		for (int i = 0; i < arrAts2DelInMol.length; i++) {
			molecule3D.setAtomFlag(arrAts2DelInMol[i], Molecule3D.FLAG1, true);
		}
		
	}
	
	/**
	 * Protonates or deprotonates the specified atom, if necessary and possible, to neutralize a charge.
	 * @param molecule3D
	 * @param
	 * @return
	 */
	public static Molecule3D getNeutralized(Molecule3D molecule3D, int indexAtomReactiveCenter) {
		Molecule3D ffNeutral = null;
		
		int charge = molecule3D.getAtomCharge(indexAtomReactiveCenter);
		
		if(charge==0){
			ffNeutral = new Molecule3D(molecule3D);
			
		} else if (charge==1){
		
			ffNeutral = getDeprotonated(molecule3D, indexAtomReactiveCenter);
			
		} else if (charge==-1){
		
			ffNeutral = getProtonated(molecule3D, indexAtomReactiveCenter);
			
		}
		
		return ffNeutral;
	}


	
	public static void copy(Molecule3D molSource, Molecule3D molDestination) {
		molDestination.clear();
		
		for (int i = 0; i < molSource.getAllAtoms(); i++) {
			molDestination.addAtom(molSource.getAtomicNo(i));
			molDestination.setAtomCharge(i, molSource.getAtomCharge(i));
			molDestination.setAtomX(i, molSource.getAtomX(i));
			molDestination.setAtomY(i, molSource.getAtomY(i));
			molDestination.setAtomZ(i, molSource.getAtomZ(i));
		}
		for (int i = 0; i < molSource.getAllBonds(); i++) {
			molDestination.addBond(molSource.getBondAtom(0, i), molSource.getBondAtom(1, i), molSource.getBondOrder(i));
		}
	}

	/**
	 * Changes the coordinates with a random value. Gaussian distribution.
	 * @param molecule3D
	 * @param maxLengthDistortion if 0 just a copy of ff will be returned.
	 * @return
	 */
	public static Molecule3D distortCoordinates(Molecule3D molecule3D, double maxLengthDistortion){
		
		Molecule3D ffShaked = new Molecule3D(molecule3D);
		
		if(maxLengthDistortion==0){
			return ffShaked;
		}
		
		Random rnd = new Random();
		
		for (int i = 0; i < ffShaked.getAllAtoms(); i++) {
			
			double [] arrRND = new double [3]; 
			
			for (int j = 0; j < arrRND.length; j++) {
				
				double d = rnd.nextGaussian() * maxLengthDistortion;
				
				if(Math.abs(d)>maxLengthDistortion){
					
					if(d>0)
						d = maxLengthDistortion;
					else
						d = maxLengthDistortion*(-1);
					
				}
				
				arrRND[j]=d;
				
			}
			
			Coordinates c = ffShaked.getCoordinates(i);
			c.x += arrRND[0];
			c.y += arrRND[1];
			c.z += arrRND[2];
						
			ffShaked.setCoordinates(i, c);
		}
		
		
		return ffShaked;
	}
	
	public static void addRNDCoordinates(Molecule3D mol, double maxNoise){
		
		Random rnd = new Random();
		
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			
			double x = mol.getAtomX(i) + (rnd.nextGaussian() * maxNoise);
			
			double y = mol.getAtomY(i) + (rnd.nextGaussian() * maxNoise);
			
			double z = mol.getAtomZ(i) + (rnd.nextGaussian() * maxNoise);
			
			mol.setAtomX(i, x);
			mol.setAtomY(i, y);
			mol.setAtomZ(i, z);
		}
	}
	

	
	private static void setRNDCoordinates(Molecule3D mol){
		double fac = 100;
		Random rnd = new Random();
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			mol.setAtomX(i,rnd.nextDouble() * fac);
			mol.setAtomY(i,rnd.nextDouble() * fac);
			mol.setAtomZ(i,rnd.nextDouble() * fac);
		}
	}
	
	public static void calculateCoordinatesIterative(Molecule3D mol, double [][] arrDist){
		
		
		int maxcycles = 10000;
		int mSize = mol.getAllAtoms();
		
		setRNDCoordinates(mol);
		
		boolean bFinished=false;
		
		double cycleFactor = 1.0;
		
		
		double start = 0;
		int cc=0;
		
		while(!bFinished){
			handleDistanceConstraint(mol, cycleFactor, arrDist);
			
			if(cc==1) {
				start = 0;
				for (int i = 0; i < mSize; i++) {
					for (int j = i+1; j < mSize; j++) {
						start += Math.abs(arrDist[i][j] - getDistanceFromCoord(mol, i,j));
					}
				}
			}
				
			if(cycleFactor <= 0) {
				bFinished=true;
			}else if(cc == maxcycles){
				bFinished = true;
			}
			
			
			cycleFactor -= 1.0/maxcycles;
			
			cc++;
		}
		
		double end = 0;
		for (int i = 0; i < mSize; i++) {
			for (int j = i+1; j < mSize; j++) {
				end += Math.abs(arrDist[i][j] - getDistanceFromCoord(mol, i,j));
			}
		}
		System.out.println("Start: " + start +" end: " + end);
		
		
	}
	
	public static double getDistanceFromCoord(Molecule3D mol, int atom1, int atom2){
		

		double dx = mol.getAtomX(atom2) - mol.getAtomX(atom1);
		double dy = mol.getAtomY(atom2) - mol.getAtomY(atom1);
		double dz = mol.getAtomZ(atom2) - mol.getAtomZ(atom1);
		double distance = Math.sqrt(dx*dx+dy*dy+dz*dz);
		
		
		return distance;
	}
	
	
	private static void handleDistanceConstraint(Molecule3D mol, double cycleFactor, double[][] mDistance) {
		Random random = new Random();
		
		int mSize = mol.getAllAtoms();


		int atom1 = (int)(random.nextDouble() * mSize);
		int atom2 = (int)(random.nextDouble() * mSize);
		while (atom2 == atom1)
			atom2 = (int)(random.nextDouble() * mSize);

		if (atom1 < atom2) {
			int temp = atom1;
			atom1 = atom2;
			atom2 = temp;
		}

		double dx = mol.getAtomX(atom2) - mol.getAtomX(atom1);
		double dy = mol.getAtomY(atom2) - mol.getAtomY(atom1);
		double dz = mol.getAtomZ(atom2) - mol.getAtomZ(atom1);
		double distance = Math.sqrt(dx*dx+dy*dy+dz*dz);

		double distanceFactor = 0.0;
		
		if (distance < mDistance[atom1][atom2]) {
			distanceFactor = (distance - mDistance[atom1][atom2])
					/ (2 * mDistance[atom1][atom2]);
			if (cycleFactor > 1.0)
				cycleFactor = 1.0;
		} else if (distance > mDistance[atom1][atom2]) {
			distanceFactor = (distance - mDistance[atom1][atom2])
					/ (2 * distance);
		}
		
		if (Math.abs(distanceFactor) > 0.0001) {
			double factor = cycleFactor * distanceFactor;
			mol.setAtomX(atom1, mol.getAtomX(atom1)+dx*factor);
			mol.setAtomX(atom2, mol.getAtomX(atom2)-dx*factor);
			mol.setAtomY(atom1, mol.getAtomY(atom1)+dy*factor);
			mol.setAtomY(atom2, mol.getAtomY(atom2)-dy*factor);
			mol.setAtomZ(atom1, mol.getAtomZ(atom1)+dz*factor);
			mol.setAtomZ(atom2, mol.getAtomZ(atom2)-dz*factor);
		}
	}

	public static double [][] getDistanceArray(Molecule3D mol) {
		
		double arr[][]= new double[mol.getAllAtoms()][mol.getAllAtoms()];
		
		for (int i = 0; i < arr.length; i++) {
			for (int j = i+1; j < arr.length; j++) {
				double dx = mol.getAtomX(i) - mol.getAtomX(j);
				double dy = mol.getAtomY(i) - mol.getAtomY(j);
				double dz = mol.getAtomZ(i) - mol.getAtomZ(j);
				double v= Math.sqrt(dx*dx+dy*dy+dz*dz);
				arr[i][j] = v;
				arr[j][i] = v;
			}
		}
		
		return arr;
	}
	
	public static boolean isParametrized(Molecule3D mol) {
		boolean bOk = true;

		for (int i = 0; i < mol.getAllAtoms(); i++) {
			if (mol.getAtomicNo(i) > 1 && mol.getInteractionAtomType(i) == -1) {
				System.out.println("Atomic no: " + mol.getAtomicNo(i));
				bOk = false;
				break;
			}
		}

		return bOk;
	}
	

	
	public static boolean isConnectedAtoms(Molecule3D ff, List<Integer> liIndexAtoms, int atm2) {
		boolean connected = false;
		
		HashSet<Integer> hsAtomIndex = new HashSet<Integer>(liIndexAtoms);

		int nConnected = ff.getConnAtoms(atm2);

		for (int i = 0; i < nConnected; i++) {
			int indexAtomConnected = ff.getConnAtom(atm2, i);
			
			if(hsAtomIndex.contains(indexAtomConnected)){
				connected = true;
				break;
				
			}
			
		}
	
		return connected;
	}
	
	
	
	/**
	 * 
	 * @param mol
	 * @param atm index of atom
	 * @return false if atom is not C or if one of the neighbours is N, O, F, S, P or Cl.
	 */
	public static boolean isAliphaticAtom(Molecule3D mol, int atm) {
		boolean aliphatic = true;
		
		if(mol.getAtomicNo(atm)!=6)
			return false;
		
		int nConn = mol.getAllConnAtoms(atm);
		
		for (int i = 0; i < nConn; i++) {
			int atmConn = mol.getConnAtom(atm, i);
			
			int atomicNoConn = mol.getAtomicNo(atmConn);
			
			if((atomicNoConn == 7) || 
			   (atomicNoConn == 8) ||
			   (atomicNoConn == 9) ||
			   (atomicNoConn == 15) ||
			   (atomicNoConn == 16) ||
			   (atomicNoConn == 17)) {
				aliphatic = false;
				break;
			}
		}
		
		return aliphatic;
	}

	public static void removeElectronPairs(Molecule3D mol) {
		for (int i = mol.getAllAtoms() - 1; i >= 0; i--) {
			if (mol.getAtomicNo(i) < 1)
				mol.deleteAtom(i);
		}
	}
	
	public static void removeHydrogensAndElectronPairs(Molecule3D mol) {
		// Remove hydrogens
		for (int i = mol.getAllAtoms() - 1; i >= 0; i--) {
			if (mol.getAtomicNo(i) <= 1)
				mol.deleteAtom(i);
		}
	}

	public static void removeElement(Molecule3D mol, int atomicNo) {
		for (int i = mol.getAllAtoms() - 1; i >= 0; i--) {
			if (mol.getAtomicNo(i) == atomicNo)
				mol.deleteAtom(i);
		}
	}

	public static void removeNonAromaticC(Molecule3D mol) {
		for (int i = mol.getAllAtoms() - 1; i >= 0; i--) {
			if (mol.getAtomicNo(i) == 6 && !mol.isAromaticAtom(i)
					&& !mol.isAtomFlag(i, FLAG_CENTER_ATOM))
				mol.deleteAtom(i);
		}

	}

	/**
	 * Removes all Carbon atoms, except the FLAG_CENTER_ATOM is set.
	 * 
	 * @param mol
	 *            Molecule
	 */
	public static void removeCarbon(Molecule3D mol) {
		for (int i = mol.getAllAtoms() - 1; i >= 0; i--) {
			if (mol.getAtomicNo(i) == 6 && !mol.isAtomFlag(i, FLAG_CENTER_ATOM))
				mol.deleteAtom(i);
		}

	}
	public static void flagCarbon(Molecule3D mol) {
		for (int i = mol.getAllAtoms() - 1; i >= 0; i--) {
			if (mol.getAtomicNo(i) == 6 && !mol.isAtomFlag(i, FLAG_CENTER_ATOM))
				mol.setAtomFlag(i, Molecule3D.FLAG1, true);
		}

	}
	
	public static void flagUnflaggedWithFlagCenter(Molecule3D mol) {
		for (int i = mol.getAllAtoms() - 1; i >= 0; i--) {
			if (!mol.isAtomFlag(i, Molecule3D.FLAG1))
				mol.setAtomFlag(i, FLAG_CENTER_ATOM, true);
		}
	}
	
	public static int[][] getConnectionTable(Molecule3D mol){
		int [][] arrConnTable = new int [mol.getAllAtoms()][mol.getAllAtoms()];
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			for (int j = 0; j < mol.getAllConnAtoms(i); j++) {
				int atConn = mol.getConnAtom(i,j);
				int bondIndex = mol.getBond(i, atConn);
				int bondOrder = mol.getBondOrder(bondIndex);
				arrConnTable[i][atConn] = bondOrder;
				arrConnTable[atConn][i] = bondOrder;
			}
		}
		return arrConnTable;
	}

	
	

	/**
	 * 
	 * @param mol
	 *            the atoms in this molecule will be deleteted.
	 * @param liIndices2Delete
	 *            list with Integer as indices for the atoms we will delete.
	 */
	public static final void removeAtomsIfCarbon(Molecule3D mol, List<Integer> liIndices2Delete) {
		Collections.sort(liIndices2Delete);
		Collections.reverse(liIndices2Delete);

		// Unify atom index list.
		int ii = 0;
		while (ii < liIndices2Delete.size() - 1) {
			if (liIndices2Delete.get(ii).equals(liIndices2Delete.get(ii + 1))) {
				liIndices2Delete.remove(ii + 1);
			} else {
				ii++;
			}
		}

		// System.out.println(liIndices2Del);

		// Order is correct, because we made a sort and a reverse before.
		for (ii = 0; ii < liIndices2Delete.size(); ii++) {
			int ind = liIndices2Delete.get(ii).intValue();
			if (mol.getAtomicNo(ind) == 6)
				mol.deleteAtom(ind);
		}

	}
	
	public static final Molecule3D removeAllAtomsWithoutNeighbours(Molecule3D ffMol) {
		Molecule3D ff = new Molecule3D(ffMol);
		
		HashSet<Integer> hsAt2Del = new HashSet<Integer>();
		for (int i = 0; i < ff.getAllAtoms(); i++) {
			if(ff.getConnAtoms(i)==0)
				hsAt2Del.add(i);
		}
		
		List<Integer> liAt2Del = new ArrayList<Integer>(hsAt2Del);
		Collections.sort(liAt2Del);
		Collections.reverse(liAt2Del);
		
		for (Integer at : liAt2Del) {
			ff.deleteAtom(at);
		}
		
		return ff;
	}

	public static final void removeAtoms(Molecule3D mol, List<Integer> liIndices2Delete) {
		Collections.sort(liIndices2Delete);
		Collections.reverse(liIndices2Delete);

		// Unify atom index list.
		int ii = 0;
		while (ii < liIndices2Delete.size() - 1) {
			if (liIndices2Delete.get(ii).equals(liIndices2Delete.get(ii + 1))) {
				liIndices2Delete.remove(ii + 1);
			} else {
				ii++;
			}
		}

		// System.out.println(liIndices2Del);

		// Order is correct, because we made a sort and a reverse before.
		for (ii = 0; ii < liIndices2Delete.size(); ii++) {
			int ind = liIndices2Delete.get(ii).intValue();
			mol.deleteAtom(ind);
		}

	}

	public static final void removeFlaggedAtoms(Molecule3D mol) {
		List<Integer> liIndices2Delete = new ArrayList<Integer>();
		
		for(int i=0; i < mol.getAllAtoms(); i++){
			if(mol.isAtomFlag(i, Molecule3D.FLAG1)) {
				liIndices2Delete.add(new Integer(i));
			}
		}
		removeAtoms(mol, liIndices2Delete);
	}

	/**
	 * 09.03.2020
	 * For what is that good?
	 * @param ffMol
	 * @param idcode2Replace
	 * @param idcodeNewSubstructure
	 * @param arrMapAtoms
	 */
	public static final void replaceAtoms(Molecule3D ffMol, String idcode2Replace, String idcodeNewSubstructure, int [][] arrMapAtoms) {
		
		StereoMolecule mol = new Molecule3D(ffMol);
		
		mol.ensureHelperArrays(Molecule.cHelperRings);
		
		
		StereoMolecule query = new StereoMolecule();
		
		query.setFragment(true);
		
		
		new IDCodeParser().parse(query, idcode2Replace);
		
			
		List<int []> liAtomLists = getListFromSSSearcher(mol, query);
		
		if(liAtomLists.size() == 0) {
			Exception ex = new Exception("Substructure not found.");
			ex.printStackTrace();
			return;
		}
		System.out.println("Found: " + liAtomLists.size());

		
		// Extract substructure from molecule
		for (Iterator<int []> iter = liAtomLists.iterator(); iter.hasNext();) {
			StereoMolecule subStruc = new StereoMolecule(mol);	
			
			int [] arrAtomList = iter.next();
			
			for (int i = mol.getAllAtoms() - 1; i >= 0; i--) {
				
				boolean bInList = false;
				for (int j = 0; j < arrAtomList.length; j++) {
					if(i == arrAtomList[j]){
						subStruc.setAtomCharge(i, i);
						bInList = true;
						break;
					}
				}
				if(!bInList){
					subStruc.setAtomSelection(i, true);
				}
			}
			
			subStruc.deleteSelectedAtoms();
			
			
		}
			
		mol.deleteSelectedAtoms();
		
		// return new Molecule3D(mol);
		return;
		
	}
	public static final void formatAtomDescriptionsForViewer(Molecule3D mol, boolean showIndex) {
		
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			String sDescription = mol.getAtomDescription(i).trim();
			
			if(sDescription==null){
				sDescription="";
			}
			
			if(showIndex)
				sDescription = LEADING_SPACE_DESCRIPTION + i + ";" + sDescription;
			else
				sDescription = LEADING_SPACE_DESCRIPTION + sDescription;
			mol.setAtomDescription(i,sDescription);
		}
	}

	public static final void formatAtomDescriptionsForViewer(Molecule3D mol) {
		
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			String sDescription = mol.getAtomDescription(i);
			
			if(sDescription!=null){
				sDescription = mol.getAtomDescription(i);
				sDescription = LEADING_SPACE_DESCRIPTION + sDescription;
				mol.setAtomDescription(i,sDescription);
			}
		}
	}

	/**
	 * Generates a center atom for each atom of a terminal alkyl group -CC, -C(C)C and
	 * -C(C)(C)C. The interaction type is the core atom of the substituents.
	 * All center atoms have identical coordinates.
	 * The original index of each considered atom is written to Molecule3D.INFO_ATOMGROUP.
	 * @param mol Molecule, the atoms used for calculating the center are flagged with FLAG1
	 * @return number of terminal alkyl groups (Ethyl, Propyl, Butyl).
	 */
	public final static int calcTerminalAlkylGroupsCenter(Molecule3D mol) {

		List<List<Integer>> liEndingAlkylGroups = findTerminalAlkylGroups(mol);

		List<Integer> liAtoms2Delete = new ArrayList<Integer>();

		for (Iterator<List<Integer>> iter = liEndingAlkylGroups.iterator(); iter.hasNext();) {
			List<Integer> liIndices = (ArrayList<Integer>) iter.next();

			int[] arrIndices = ArrayUtilsCalc.toIntArray(liIndices);
			Coordinates coord = getCenterGravity(mol, arrIndices);

			for (int at = 0; at < arrIndices.length; at++) {
				// The interaction type is the core atom of the substituents.
				int iInteractionType = mol.getInteractionAtomType(arrIndices[at]);
				String sDescription = mol.getAtomDescription(arrIndices[at]);
				// Major MM2 interaction type
				// int iMM2Type = mol.getMM2AtomType(arrIndices[at]);

				int index = mol.getAllAtoms();
				
				
				mol.addAtom(6);
				mol.setInteractionAtomType(index, iInteractionType);
				// mol.setMM2AtomType(index, iMM2Type);
				
				String sOrigIndex = Integer.toString(arrIndices[at]);
				mol.setAtomChainId(index, sOrigIndex);

				
				if(sDescription!=null)
					mol.setAtomDescription(index, sDescription);

				mol.setCoordinates(index, coord);
				
				// Set original index
				mol.setAtomChainId(index, arrIndices[at]+"");
				
				mol.setAtomFlag(index, FLAG_CENTER_ATOM, true);
			}

			liAtoms2Delete.addAll(liIndices);
		}
		
		for (int i = 0; i < liAtoms2Delete.size(); i++) {
			int index = liAtoms2Delete.get(i).intValue();
			mol.setAtomFlag(index, Molecule3D.FLAG1, true);
		}

		return liEndingAlkylGroups.size();

	}

	/**
	 * Finds terminal alkyl groups
	 * 
	 * @param mol
	 * @return List of Lists with Integer containing the atom indices. At pos 0
	 *         is the index for the core atom given.
	 */
	private final static List<List<Integer>> findTerminalAlkylGroups(Molecule3D mol) {
		List<List<Integer>> liAlkylGroups = new ArrayList<List<Integer>>();

		Hashtable<Integer, List<Integer>> ht = new Hashtable<Integer, List<Integer>>();

		int[] arrMethylGroups = findTerminalMethylAtCarb(mol);
		for (int i = 0; i < arrMethylGroups.length; i++) {
			int nConnAts = mol.getAllConnAtoms(arrMethylGroups[i]);
			for (int j = 0; j < nConnAts; j++) {
				int indNeighborAt = mol.getConnAtom(arrMethylGroups[i], j);
				int atomicNumberNeighbor = mol.getAtomicNo(indNeighborAt);
				if (atomicNumberNeighbor > 1) {
					Integer intIndNeighborAt = new Integer(indNeighborAt);
					if (ht.containsKey(intIndNeighborAt)) {
						List<Integer> li = ht.get(intIndNeighborAt);
						li.add(new Integer(arrMethylGroups[i]));
					} else {
						List<Integer> li = new ArrayList<Integer>();
						li.add(intIndNeighborAt);
						li.add(new Integer(arrMethylGroups[i]));
						ht.put(intIndNeighborAt, li);
					}
				}
			}
		}

		for (Enumeration<List<Integer>> en = ht.elements(); en.hasMoreElements();) {
			List<Integer> element = en.nextElement();
			liAlkylGroups.add(element);
		}

		return liAlkylGroups;

	}

	/**
	 * 
	 * @param mol
	 * @return indices for all terminal methyl CH3 groups which are connected to
	 *         a carbon.
	 */
	public final static int[] findTerminalMethylAtCarb(Molecule3D mol) {

		int numAts = mol.getAllAtoms();
		int[] arrIndexHeavyNeighbor = new int[numAts];
		for (int i = 0; i < arrIndexHeavyNeighbor.length; i++) {
			arrIndexHeavyNeighbor[i] = -1;
		}

		int numMethyl = 0;
		for (int i = 0; i < numAts; i++) {
			if (mol.getAtomicNo(i) == 6) {
				int nConnAts = mol.getAllConnAtoms(i);
				int nHeavyNeighAts = 0;
				int indexNeighborAt = -1;
				for (int j = 0; j < nConnAts; j++) {
					int indAt = mol.getConnAtom(i, j);
					int atomicNumber = mol.getAtomicNo(indAt);
					if (atomicNumber > 1) {
						nHeavyNeighAts++;
						indexNeighborAt = indAt;
					}
				}

				if (nHeavyNeighAts == 1
						&& mol.getAtomicNo(indexNeighborAt) == 6) {

					arrIndexHeavyNeighbor[i] = indexNeighborAt;
					numMethyl++;
				}
			}
		}

		int[] arrIndexMethylGroups = new int[numMethyl];

		int cc = 0;
		for (int i = 0; i < arrIndexHeavyNeighbor.length; i++) {
			if (arrIndexHeavyNeighbor[i] > -1
					&& mol.getAtomicNo(arrIndexHeavyNeighbor[i]) == 6) {
				arrIndexMethylGroups[cc++] = i;
			}
		}

		return arrIndexMethylGroups;
	}

	public final static Coordinates getCenterGravity(Molecule3D mol) {
		
		int n = mol.getAllAtoms();
		
		int [] indices = new int [n];
		
		for (int i = 0; i < indices.length; i++) {
			indices[i]=i;
		}

		return getCenterGravity(mol, indices);
	}

	public final static Coordinates getCenterGravity(Molecule3D mol, int[] indices) {
		
		Coordinates c = new Coordinates();
		for (int i = 0; i < indices.length; i++) {
			c.x += mol.getAtomX(indices[i]);
			c.y += mol.getAtomY(indices[i]);
			c.z += mol.getAtomZ(indices[i]);
		}
		c.x /= indices.length;
		c.y /= indices.length;
		c.z /= indices.length;

		return c;
	}
	
	public final static Molecule3D getCentered(Molecule3D ff) {
		Molecule3D ffCent = new Molecule3D(ff);
		
		final int n = ff.getAllAtoms();
		
		double meanX = 0, meanY = 0, meanZ = 0;
		
		for (int i = 0; i < n; i++) {
			Coordinates c = ff.getCoordinates(i);
			
			meanX += c.x;
			meanY += c.y;
			meanZ += c.z;
			
		}
		
		meanX /= n;
		meanY /= n;
		meanZ /= n;
		
		for (int i = 0; i < n; i++) {
			Coordinates c = ffCent.getCoordinates(i);
			
			c.x -= meanX;
			c.y -= meanY;
			c.z -= meanZ;
		}
		
		return ffCent;
	}
	
	public final static Molecule3D getCentered(Molecule3D molecule3D, List<Integer> liAtomIndexCenter) {
		Molecule3D ffCent = new Molecule3D(molecule3D);
		
		double meanX = 0, meanY = 0, meanZ = 0;
		
		for (int i = 0; i < liAtomIndexCenter.size(); i++) {
			Coordinates c = ffCent.getCoordinates(liAtomIndexCenter.get(i));
			
			meanX += c.x;
			meanY += c.y;
			meanZ += c.z;
			
		}
		
		meanX /= liAtomIndexCenter.size();
		meanY /= liAtomIndexCenter.size();
		meanZ /= liAtomIndexCenter.size();
		
		for (int i = 0; i < ffCent.getAllAtoms(); i++) {
			Coordinates c = ffCent.getCoordinates(i);
			
			c.x -= meanX;
			c.y -= meanY;
			c.z -= meanZ;
		}
		
		return ffCent;
	}

	/**
	 * Topological centaer atoms are the atoms with the lowest squared sum of topological distances to all atoms.
	 * @param
	 * @return
	 */
	public final static List<Integer> getTopologicalCenter(int [][] arrTopoDist) {
		List<Integer> liTopoCenterAtoms = new ArrayList<Integer>();

		List<IntegerDouble> li = new ArrayList<IntegerDouble>();
		for (int i = 0; i < arrTopoDist.length; i++) {
			double sum=0;
			for (int j = 0; j < arrTopoDist.length; j++) {
				sum += arrTopoDist[i][j]*arrTopoDist[i][j];
			}
			
			li.add(new IntegerDouble(i,sum));
		}
		
		Collections.sort(li, IntegerDouble.getComparatorDouble());
		
		liTopoCenterAtoms.add(li.get(0).getInt());
		
		for (int i = 1; i < li.size(); i++) {
			if(li.get(i).getDouble()==li.get(0).getDouble()){
				liTopoCenterAtoms.add(li.get(i).getInt());
			}
		}
		
		return liTopoCenterAtoms;
	}
	
	public final static List<Integer> getTopologicalCenter(int [][] arrTopoDist, ListWithIntVec ilIndexAtoms) {
		List<Integer> liTopoCenterAtoms = new ArrayList<Integer>();

		List<IntegerDouble> li = new ArrayList<IntegerDouble>();
		
		for (int i = 0; i < ilIndexAtoms.size(); i++) {
			
			int indexAt1 = ilIndexAtoms.get(i);
			
			double sum=0;
			for (int j = 0; j < ilIndexAtoms.size(); j++) {
				
				int indexAt2 = ilIndexAtoms.get(j);
				
				sum += arrTopoDist[indexAt1][indexAt2]*arrTopoDist[indexAt1][indexAt2];
			}
			
			li.add(new IntegerDouble(indexAt1,sum));
		}
		
		Collections.sort(li, IntegerDouble.getComparatorDouble());
		
		liTopoCenterAtoms.add(li.get(0).getInt());
		
		for (int i = 1; i < li.size(); i++) {
			if(li.get(i).getDouble()==li.get(0).getDouble()){
				liTopoCenterAtoms.add(li.get(i).getInt());
			}
		}
		
		return liTopoCenterAtoms;
	}

	/**
	 * Gets the points with the maximum sum of squared topological distances to all atoms.

	 * @param arrTopoDist
	 * @return
	 */
	public final static List<Integer> getPheriphericPoints(int [][] arrTopoDist) {
		
		List<Integer> liTopoCenterAtoms = new ArrayList<Integer>();

		List<IntegerDouble> li = new ArrayList<IntegerDouble>();
		for (int i = 0; i < arrTopoDist.length; i++) {
			double sum=0;
			for (int j = 0; j < arrTopoDist.length; j++) {
				sum += arrTopoDist[i][j]*arrTopoDist[i][j];
			}
			
			li.add(new IntegerDouble(i,sum));
		}
		
		Collections.sort(li, IntegerDouble.getComparatorDouble());
		
		Collections.reverse(li);
		
		liTopoCenterAtoms.add(li.get(0).getInt());
		
		for (int i = 1; i < li.size(); i++) {
			if(li.get(i).getDouble()==li.get(0).getDouble()){
				liTopoCenterAtoms.add(li.get(i).getInt());
			}
		}
		
		return liTopoCenterAtoms;
	}

	public final static int [][] getTopologicalDistanceMatrix(Molecule3D mol) {
		
		return StructureCalculator.getNumberOfBondsBetweenAtoms(mol, mol.getAllBonds(), null);
		
	}
	
	public final static List<Integer> getLongestChain(Molecule3D ff) {
		
		int [][] arrTopoDist = Molecule3DFunctions.getTopologicalDistanceMatrix(ff);
		
		int maxDist=0;
		
		int indexAt=-1;
		for (int i = 0; i < arrTopoDist.length; i++) {
			// We will not start at Hydrogen or electron pairs.
			if(ff.getAtomicNo(i)<2)
				continue;
			
			for (int j = i+1; j < arrTopoDist.length; j++) {
				if(arrTopoDist[i][j]>maxDist){
					maxDist=arrTopoDist[i][j];
					indexAt=i;
				}
			}
		}
		
		List<Integer> liLongestChain = getLongestChain(ff, indexAt);

		return liLongestChain;
	}
	
	/**
	 * User object in DefaultMutableTreeNode is the index of the atom. Root contains indexAtomStart.
	 * @param molecule3D
	 * @param indexAtomStart
	 * @return
	 */
	public final static DefaultMutableTreeNode getTreeFromBroadFirstSearch(Molecule3D molecule3D, int indexAtomStart){
		
		HashSet<Integer> hsVisited = new HashSet<Integer>();
		
		hsVisited.add(indexAtomStart);
		List<DefaultMutableTreeNode> liQueue = new ArrayList<DefaultMutableTreeNode>();
		
		DefaultMutableTreeNode root = new DefaultMutableTreeNode();
		
		root.setUserObject(indexAtomStart);
		
		liQueue.add(root);
		
		while(!liQueue.isEmpty()){
			DefaultMutableTreeNode parent = liQueue.remove(0);
			
			int indexAtom = (Integer)parent.getUserObject();
			
			int nConnceted = molecule3D.getConnAtoms(indexAtom);
			
			for (int i = 0; i < nConnceted; i++) {
				int indexAtomConnceted = molecule3D.getConnAtom(indexAtom, i);
				
				if(molecule3D.getAtomicNo(indexAtomConnceted)<2)
					continue;
				
				if(!hsVisited.contains(indexAtomConnceted)) {
					
					hsVisited.add(indexAtomConnceted);
					
					DefaultMutableTreeNode child = new DefaultMutableTreeNode();
					
					child.setUserObject(indexAtomConnceted);
					
					liQueue.add(child);
					
					parent.add(child);
					
				}
			}
		}
		
		return root;
	}

	/**
	 * Get all possible paths. Only crossing is not allowed.
	 * @param molecule3D
	 * @param indexAtomStart
	 * @return
	 */
	public final static DefaultMutableTreeNode getTreeFromComplete(Molecule3D molecule3D, int indexAtomStart){
		
		List<DefaultMutableTreeNode> liQueue = new ArrayList<DefaultMutableTreeNode>();
		
		DefaultMutableTreeNode root = new DefaultMutableTreeNode();
		
		root.setUserObject(indexAtomStart);
		
		liQueue.add(root);
		
		while(!liQueue.isEmpty()){
			DefaultMutableTreeNode parent = liQueue.remove(0);
			
			int indexAtom = (Integer)parent.getUserObject();
			
			int nConnceted = molecule3D.getConnAtoms(indexAtom);
			
			for (int i = 0; i < nConnceted; i++) {
				int indexAtomConnceted = molecule3D.getConnAtom(indexAtom, i);
				
				if(molecule3D.getAtomicNo(indexAtomConnceted)<2)
					continue;
				
				// Check for crossing
				if(!isInPath(parent, indexAtomConnceted)) {
					
					DefaultMutableTreeNode child = new DefaultMutableTreeNode();
					
					child.setUserObject(indexAtomConnceted);
					
					liQueue.add(child);
					
					parent.add(child);
					
				}
			}
		}
		
		return root;
	}
	
	private static final boolean isInPath(DefaultMutableTreeNode parent, int indexAtom){
		boolean inPath=false;
		
		TreeNode [] arrPath = parent.getPath();
		
		for (int i = 0; i < arrPath.length; i++) {
			DefaultMutableTreeNode node = (DefaultMutableTreeNode)arrPath[i];
			
			int indexAtomNode = (Integer)node.getUserObject();
			
			if(indexAtomNode==indexAtom){
				inPath=true;
				break;
			}
			
		}
		
		return inPath;
	}
	
	
	@SuppressWarnings("unchecked")
	public final static int [] getPath(Molecule3D ff, int indexAt1, int indexAt2){
		
		int [] arrPath = null;
		
		DefaultMutableTreeNode root = getTreeFromBroadFirstSearch(ff, indexAt1);
		
		Enumeration<TreeNode> en = root.breadthFirstEnumeration();
		
		for(;en.hasMoreElements();){
			DefaultMutableTreeNode node = (DefaultMutableTreeNode)en.nextElement();
			
			int indexAtomNode = (Integer)node.getUserObject();
			
			if(indexAt2==indexAtomNode) {
				TreeNode [] arrPathNodes = node.getPath();
				
				arrPath = new int [arrPathNodes.length];
				
				for (int i = 0; i < arrPathNodes.length; i++) {
					arrPath[i]=(Integer)((DefaultMutableTreeNode)arrPathNodes[i]).getUserObject();
				}
			}
			
		}
		
		return arrPath;
	}

	
	/**
	 * Returns the atom indices as layers of the start atom. First layer is indexAtomStart.
	 * @param molecule3D
	 * @param indexAtomStart
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public final static List<List<Integer>> getLayersFromBroadFirstSearch(Molecule3D molecule3D, int indexAtomStart){
		
		DefaultMutableTreeNode root = Molecule3DFunctions.getTreeFromBroadFirstSearch(molecule3D, indexAtomStart);
		
		Enumeration<TreeNode> en = root.breadthFirstEnumeration();
		
		List<List<Integer>> liliIndexAtomLayer = new ArrayList<List<Integer>>();
		
		List<Integer> liIndexAtomLayer = new ArrayList<Integer>();
		
		int level = 0;
		while(en.hasMoreElements()){
			DefaultMutableTreeNode node = (DefaultMutableTreeNode)en.nextElement();
			
			if(node.getLevel()>level){
				liliIndexAtomLayer.add(liIndexAtomLayer);
				liIndexAtomLayer = new ArrayList<Integer>();
				level++;
			}
			liIndexAtomLayer.add((Integer)node.getUserObject());
		}
		
		if(liIndexAtomLayer.size()>0){
			liliIndexAtomLayer.add(liIndexAtomLayer);
		}
		
		return liliIndexAtomLayer;
	}

	public final static List<Integer> getLongestChain(Molecule3D molecule3D, int indexAtomStart) {
		
		HashSet<Integer> hsVisited = new HashSet<Integer>();
		
		hsVisited.add(indexAtomStart);
		
		List<DefaultMutableTreeNode> liQueue = new ArrayList<DefaultMutableTreeNode>();
		
		DefaultMutableTreeNode root = new DefaultMutableTreeNode();
		
		root.setUserObject(indexAtomStart);
		
		liQueue.add(root);
		
		DefaultMutableTreeNode nodeDeepest =  root;
		
		while(!liQueue.isEmpty()){
			DefaultMutableTreeNode parent = liQueue.remove(0);
			
			int indexAtom = (Integer)parent.getUserObject();
			
			int nConnceted = molecule3D.getConnAtoms(indexAtom);
			
			for (int i = 0; i < nConnceted; i++) {
				int indexAtomConnceted = molecule3D.getConnAtom(indexAtom, i);
				
				if(molecule3D.getAtomicNo(indexAtomConnceted)<2)
					continue;
				
				if(!hsVisited.contains(indexAtomConnceted)) {
					
					hsVisited.add(indexAtomConnceted);
					
					DefaultMutableTreeNode child = new DefaultMutableTreeNode();
					
					child.setUserObject(indexAtomConnceted);
					
					liQueue.add(child);
					
					parent.add(child);
					
					if(child.getLevel()>nodeDeepest.getLevel()){
						
						nodeDeepest = child;
					}
				}
			}
		}
		
		TreeNode [] tnPath = nodeDeepest.getPath();
		
		List<Integer> liChain = new ArrayList<Integer>();
		for (int i = 0; i < tnPath.length; i++) {
			DefaultMutableTreeNode node = (DefaultMutableTreeNode)tnPath[i];
			
			liChain.add((Integer)node.getUserObject());
		}
		
		return liChain;
		
	}
	
	
	
	public final static String toString(Molecule3D mol){
		StringBuilder sb = new StringBuilder();
		
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			sb.append(mol.getAtomicNo(i));
			sb.append(":");
			sb.append(NF.format(mol.getAtomX(i)));
			sb.append(",");
			sb.append(NF.format(mol.getAtomY(i)));
			sb.append(",");
			sb.append(NF.format(mol.getAtomZ(i)));
			sb.append("\n");
		}
		
		
		return sb.toString();
	}
	
	/**
	 * Deprotonates and decreases charge by 1.
	 * @param molecule3D
	 * @param indexAtomReactiveCenter
	 * @return null if no hydrogen is attached to indexAtomReactiveCenter.
	 */
	public static Molecule3D getDeprotonated(Molecule3D molecule3D, int indexAtomReactiveCenter){
		
		int indexProtonAcid = -1;
		
		int nConnected2Donor = molecule3D.getAllConnAtoms(indexAtomReactiveCenter);
		
		for (int i = 0; i < nConnected2Donor; i++) {
			int indexAtomConn = molecule3D.getConnAtom(indexAtomReactiveCenter, i);
			
			if(molecule3D.getAtomicNo(indexAtomConn)==1){
				indexProtonAcid = indexAtomConn;
				break;
			}
		}
		
		if(indexProtonAcid==-1){
			return null;
		}
		
		Molecule3D ffDeprotonated = new Molecule3D(molecule3D);
		
		ffDeprotonated.deleteAtom(indexProtonAcid);
		
		int charge = ffDeprotonated.getAtomCharge(indexAtomReactiveCenter) - 1;
		
		ffDeprotonated.setAtomCharge(indexAtomReactiveCenter, charge);
		
		return ffDeprotonated;
	}
	
	/**
	 * Protonates and increases charge by 1.
	 * @param molecule3D
	 * @param
	 * @return
	 */
	public static Molecule3D getProtonated(Molecule3D molecule3D, int indexAtomReactiveCenter){
		
		Molecule3D ffProtonated = new Molecule3D(molecule3D);
		
		int atomicNo = molecule3D.getAtomicNo(indexAtomReactiveCenter);
		
		int connAtoms = molecule3D.getConnAtoms(indexAtomReactiveCenter);
		
		int sumBondOrder = 0;
		for (int i = 0; i < connAtoms; i++) {
			sumBondOrder += molecule3D.getConnBondOrder(indexAtomReactiveCenter, i);
		}
		
		if(atomicNo == 6){
			if(sumBondOrder == 4){ // C
				return null;
			}
		} else if(atomicNo == 7){ // N
			if(sumBondOrder == 4){
				return null;
			}
		} else if(atomicNo == 8){ // O
			if(sumBondOrder == 3){
				return null;
			}
		}
		
		int indexProton = ffProtonated.addAtom(1);
		
		int charge = ffProtonated.getAtomCharge(indexAtomReactiveCenter) + 1;
		
		ffProtonated.setAtomCharge(indexAtomReactiveCenter, charge);
		
		ffProtonated.addBond(indexAtomReactiveCenter, indexProton, 1);
		
		return ffProtonated;
	}

	public static Molecule3D randomizeAtoms(Molecule3D molecule3D){
		
		
		Molecule3D ffPure = new Molecule3D(molecule3D);
		
		removeElectronPairs(ffPure);
		
		List<Integer> liAtomIndex = new ArrayList<Integer>();
		
		for (int i = 0; i < ffPure.getAllAtoms(); i++) {
			liAtomIndex.add(i);
		}
		
		Collections.shuffle(liAtomIndex);
		
		int [] atomMapRND2Mol = ArrayUtilsCalc.toIntArray(liAtomIndex);
		int [] atomMapMol2RNDl = new int [ffPure.getAllAtoms()];
		
		Molecule3D ffRND = new Molecule3D(ffPure.getAllAtoms(), ffPure.getAllBonds());
		
		
		for (int i = 0; i < atomMapRND2Mol.length; i++) {
			
			int indexAtSource = atomMapRND2Mol[i];
			
			int indexAtRND = ffRND.addAtom(ffPure.getAtomicNo(indexAtSource));
			
			Coordinates c = ffPure.getCoordinates(indexAtSource);
			
			ffRND.setCoordinates(indexAtRND, c);
			
			atomMapMol2RNDl[indexAtSource]=indexAtRND;
		}
		
		int nBonds = ffPure.getAllBonds();
		
		for (int i = 0; i < nBonds; i++) {
			
			int indexAt1 = ffPure.getBondAtom(0, i);
			int indexAt2 = ffPure.getBondAtom(1, i);
			
			int indexAtRND1 = atomMapMol2RNDl[indexAt1];
			
			int indexAtRND2 = atomMapMol2RNDl[indexAt2];
			
			int order = ffPure.getBondOrder(i);
			
			ffRND.setBondAtom(0, i, indexAtRND1);
			ffRND.setBondAtom(1, i, indexAtRND2);
			
			ffRND.addBond(indexAtRND1, indexAtRND2, order);
			
		}
		
		return ffRND;
	}

	
}
