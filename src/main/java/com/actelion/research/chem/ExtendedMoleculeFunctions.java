/*
 * Copyright (c) 1997 - 2016
 * Actelion Pharmaceuticals Ltd.
 * Gewerbestrasse 16
 * CH-4123 Allschwil, Switzerland
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the the copyright holder nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
* Copyright (c) 1997 - 2016
* Actelion Pharmaceuticals Ltd.
* Gewerbestrasse 16
* CH-4123 Allschwil, Switzerland
*
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 3. Neither the name of the the copyright holder nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

package com.actelion.research.chem;

import com.actelion.research.calc.ArrayUtilsCalc;
import com.actelion.research.chem.descriptor.DescriptorEncoder;
import com.actelion.research.chem.descriptor.DescriptorHandler;
import com.actelion.research.util.BurtleHasher;

import java.util.*;

/**
 *
 *
 * ExtendedMoleculeFunctions
 * @author Modest von Korff
 * @version 1.0
 * 2005 MvK: Start implementation
 */
public class ExtendedMoleculeFunctions {


	public static final int NUMBER_ELEMENTS = 118;

	public static final int COLOR_FOR_CENTER_SELECTION = Molecule.cAtomColorRed;


	public static final int [] arrRGroupsAtomicNo = {142,143,144,129,130,131,132,133,134,135,136,137,138,139,140,141};

	public static final String [] arrRGroupsSymbol = {"R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11","R12","R13","R14","R15","R16"};


	public static void makeSkeleton(StereoMolecule mol) {
		for (int bond=0; bond<mol.getAllBonds(); bond++)
			mol.setBondType(bond, Molecule.cBondTypeSingle);
		for (int atom=0; atom<mol.getAllAtoms(); atom++)
			mol.setAtomicNo(atom, 6);
	}

	public static void analyzeMolecule(StereoMolecule mol){
		int nAtoms = mol.getAtoms();

		for (int atom = 0; atom < nAtoms; atom++) {
			if(AtomFunctionAnalyzer.isAcidicOxygen(mol, atom)){
				System.out.println("Acidic oxygen " + atom);
			}
			if(AtomFunctionAnalyzer.isNitroGroupN(mol, atom)){
				System.out.println("Nitro group N " + atom);
			}
			if(AtomFunctionAnalyzer.isMemberOfNitroGroup(mol, atom)){
				System.out.println("Member of group " + atom);
			}
			if(AtomFunctionAnalyzer.isBasicNitrogen(mol, atom)){
				System.out.println("Basic nitrogen " + atom);
			}



		}

		System.out.println("Index\tAtomic no\tAtom type");
		for (int indexAtom = 0; indexAtom < mol.getAllAtoms(); indexAtom++) {

			long atomType=-1;
			try {
				atomType = AtomTypeCalculator.getAtomType(mol, indexAtom, AtomTypeCalculator.cPropertiesForCLogPCharges);
			} catch (Exception e) {
				e.printStackTrace();
			}

			System.out.println(indexAtom + "\t" + mol.getAtomicNo(indexAtom) + "\t" + atomType + "\t" + AtomTypeCalculator.toString(atomType));
		}

		System.out.println();
	}




	public static String getColorVal2String(Molecule mol, int indexAtom){

		int color = mol.getAtomColor(indexAtom);

		return getcAtomColor2String(color);
	}

	/**
	 *
	 * @param cAtomColor atom color value from Molecule.
	 * @return
	 */
	public static String getcAtomColor2String(int cAtomColor){
		String sColor = "";

		switch (cAtomColor) {

			case Molecule.cAtomColorNone:
				sColor = "black";
				break;
			case Molecule.cAtomColorBlue:
				sColor = "blue";
				break;
			case Molecule.cAtomColorRed:
				sColor = "red";
				break;
			case Molecule.cAtomColorGreen:
				sColor = "green";
				break;
			case Molecule.cAtomColorMagenta:
				sColor = "magenta";
				break;
			case Molecule.cAtomColorOrange:
				sColor = "orange";
				break;
			case Molecule.cAtomColorDarkGreen:
				sColor = "darkGreen";
				break;
			case Molecule.cAtomColorDarkRed:
				sColor = "darkRed";
				break;

			default:
				break;
		}


		return sColor;
	}

	/**
	 *
	 * @param mol
	 * @param liIndexAtom
	 * @param cAtomColor atom color value from Molecule.
	 * @return
	 */
	public static String getColorRecord(Molecule mol, Collection<Integer> liIndexAtom, int cAtomColor){

		StringBuilder sb = new StringBuilder();

		int cc=0;
		for (int indexAtom : liIndexAtom) {

			String sColor = getcAtomColor2String(cAtomColor);

			sColor += ":"+indexAtom;

			sb.append(sColor);

			if(cc<liIndexAtom.size()-1){
				sb.append(";");
			}
			cc++;
		}

		return sb.toString();
	}

	public static int getNumQueryAtoms(ExtendedMolecule mol, int [] arrAtomicNoQuery) {
		int hetero = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {
			int atomicNo = mol.getAtomicNo(i);
			for (int j = 0; j < arrAtomicNoQuery.length; j++) {
				if(arrAtomicNoQuery[j]==atomicNo){
					hetero++;
					break;
				}
			}
		}

		return hetero;
	}

	public static int getNumBondsNoHydrogen(ExtendedMolecule mol) {
		int bnds = 0;
		for (int i = 0; i < mol.getBonds(); i++) {
			int at1 = mol.getBondAtom(0,i);
			int atNo1 = mol.getAtomicNo(at1);
			int at2 = mol.getBondAtom(1,i);
			int atNo2 = mol.getAtomicNo(at2);
			if(atNo1!=1 && atNo2!=1) {
				bnds++;
			}
		}
		return bnds;
	}
	public static int getNumNonHydrogenAtoms(ExtendedMolecule mol) {
		int non = 0;

		for (int i = 0; i < mol.getAtoms(); i++) {
			int atomicNo = mol.getAtomicNo(i);
			if(atomicNo!=1) {
				non++;
			}
		}

		return non;
	}
	public static int getNumCarbonAtoms(ExtendedMolecule mol) {
		int carbon = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {
			int atomicNo = mol.getAtomicNo(i);
			if(atomicNo==6) {
				carbon++;
			}
		}

		return carbon;
	}

	public static int getNumHeteroAtoms(ExtendedMolecule mol) {
		int hetero = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {
			int atomicNo = mol.getAtomicNo(i);
			if((atomicNo>1) &&
					(atomicNo != 6) &&
					(atomicNo <= NUMBER_ELEMENTS)) {
				hetero++;
			}
		}

		return hetero;
	}

	public static boolean isHetero(ExtendedMolecule mol, int indexAtom){

		boolean hetero = false;

		int atomicNo = mol.getAtomicNo(indexAtom);

		if((atomicNo>1) &&
				(atomicNo != 6) &&
				(atomicNo <= NUMBER_ELEMENTS)) {
			hetero = true;
		}

		return hetero;
	}

	public static int getNumNitroGroupN(StereoMolecule mol) {
		int n = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(AtomFunctionAnalyzer.isNitroGroupN(mol, i)){
				n++;
			}
		}

		return n;
	}

	public static int getNumAmide(StereoMolecule mol) {
		int n = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(AtomFunctionAnalyzer.isAmide(mol, i)){
				n++;
			}
		}

		return n;
	}

	public static int getNumCarboxy(StereoMolecule mol) {
		int n = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(isCarboxyC(mol, i)){
				n++;
			}
		}

		return n;
	}


	public static int getNumAcidicOxygen(StereoMolecule mol) {
		int n = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(AtomFunctionAnalyzer.isAcidicOxygen(mol, i)){
				n++;
			}
		}

		return n;
	}

	public static int getNumBasicNitrogen(StereoMolecule mol) {
		int n = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(AtomFunctionAnalyzer.isBasicNitrogen(mol, i)){
				n++;
			}
		}

		return n;
	}
	public static int getNumAliphaticRingAtoms(ExtendedMolecule mol, int atomicNoQuery) {
		int n = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {
			int atomicNo = mol.getAtomicNo(i);

			if(mol.isRingAtom(i) && atomicNo==atomicNoQuery){
				n++;
			}
		}

		return n;
	}
	/**
	 *
	 * @param mol
	 * @return number of atoms which are not hydrogen.
	 */
	public static int getNumHeavyAtoms(ExtendedMolecule mol) {
		int heavy = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {
			int atomicNo = mol.getAtomicNo(i);
			if(atomicNo>1 && atomicNo <= NUMBER_ELEMENTS) {
				heavy++;
			}
		}

		return heavy;
	}

	public static int getNumAromaticAtoms(ExtendedMolecule mol) {

		int n=0;
		for (int i = 0; i < mol.getAtoms(); i++) {
			if(mol.isAromaticAtom(i)){
				n++;
			}
		}
		return n;
	}

	public static int getNumArylAmine(StereoMolecule mol) {

		int n=0;
		for (int i = 0; i < mol.getAtoms(); i++) {
			if(AtomFunctionAnalyzer.isArylAmine(mol, i)){
				n++;
			}
		}
		return n;
	}

	public static int getNumHeteroAromaticAtoms(ExtendedMolecule mol) {

		int n=0;
		for (int i = 0; i < mol.getAtoms(); i++) {
			if(mol.isRingAtom(i) && mol.isAromaticAtom(i) && mol.getAtomicNo(i)!= 6){
				n++;
			}
		}
		return n;
	}

	public static int getNumIsopropyl(ExtendedMolecule mol) {

		int n=0;
		for (int i = 0; i < mol.getAtoms(); i++) {
			if(mol.isRingAtom(i) && mol.isAromaticAtom(i) && mol.getAtomicNo(i)!= 6){
				n++;
			}
		}
		return n;
	}

	public static int getNumSubstructure(StereoMolecule mol, String idcodeFragment) {

		IDCodeParser parser = new IDCodeParser();

		StereoMolecule frag = parser.getCompactMolecule(idcodeFragment);

		frag.ensureHelperArrays(Molecule.cHelperRings);

		SSSearcher ssSearcher = new SSSearcher();

		ssSearcher.setMol(frag, mol);

		int nMatches = ssSearcher.findFragmentInMolecule();

		return nMatches;
	}



	public static int getBondNo(ExtendedMolecule mol, int atm1,int atm2) {
		for (int bnd=0; bnd<mol.getAllBonds(); bnd++)
			if (((mol.getBondAtom(0,bnd) == atm1)
					&& (mol.getBondAtom(1,bnd) == atm2))
					|| ((mol.getBondAtom(0,bnd) == atm2)
					&& (mol.getBondAtom(1,bnd) == atm1)))
				return bnd;
		return -1;
	}

	public static int getBondOrder(ExtendedMolecule mol, int atm1,int atm2) {

		int bndno = getBondNo(mol, atm1,atm2);
		if(bndno==-1)
			return -1;

		int bndord = mol.getBondOrder(bndno);

		return bndord;
	}

	public static int getBondType(ExtendedMolecule mol, int atm1,int atm2) {

		int bndno = getBondNo(mol, atm1,atm2);
		if(bndno==-1)
			return -1;

		int bndtype = mol.getBondType(bndno);

		return bndtype;
	}

	public static int getBondParity(ExtendedMolecule mol, int atm1,int atm2) {

		int bndno = getBondNo(mol, atm1,atm2);
		if(bndno==-1)
			return -1;

		int bndtype = mol.getBondParity(bndno);

		return bndtype;
	}

	public static boolean deleteBond(ExtendedMolecule mol, int atm1,int atm2) {
		int bndno = getBondNo(mol, atm1,atm2);
		if(bndno==-1)
			return false;

		mol.deleteBond(bndno);

		return true;
	}

	public static String getBiggestFragmentIDCode(String idCode) {

		StereoMolecule mol = new StereoMolecule();
		new IDCodeParser(false).parse(mol, idCode);

		StereoMolecule [] frags = mol.getFragments();

		String idBiggest = idCode;

		if(frags.length > 1) {

			StereoMolecule biggest = mol;
			int maxAtoms = 0;
			for (int ii = 0; ii < frags.length; ii++) {
				if (frags[ii].getAllAtoms() > maxAtoms) {
					maxAtoms = frags[ii].getAllAtoms();
					biggest = frags[ii];
				}
			}
			Canonizer canFrag = new Canonizer(biggest);
			idBiggest = canFrag.getIDCode();
		}

		return idBiggest;
	}

	public static StereoMolecule getBiggestFragment(StereoMolecule mol) {

		StereoMolecule [] frags = mol.getFragments();

		StereoMolecule biggest = mol;

		if(frags.length > 1) {

			int maxAtoms = 0;
			for (int ii = 0; ii < frags.length; ii++) {
				if (frags[ii].getAllAtoms() > maxAtoms) {
					maxAtoms = frags[ii].getAllAtoms();
					biggest = frags[ii];
				}
			}
		}
		return biggest;
	}

	/**
	 * Replaces all hetero-atoms, except hydrogen, with carbon.
	 * @param m
	 * @return
	 */
	public static StereoMolecule getConverted2CarbonSkeleton(StereoMolecule m){

		StereoMolecule skel = new StereoMolecule(m);

		skel.ensureHelperArrays(Molecule.cHelperRings);

		for (int i = 0; i < skel.getAtoms(); i++) {

			if(skel.getAtomicNo(i)>1) {
				skel.setAtomicNo(i, 6);
			}
		}

		skel.ensureHelperArrays(Molecule.cHelperRings);

		return skel;

	}

	public static Comparator<StereoMolecule> getComparatorAtomsBonds(){

		return new Comparator<StereoMolecule>() {

			@Override
			public int compare(StereoMolecule m1, StereoMolecule m2) {

				if(m1.getAllAtoms() > m2.getAllAtoms()){
					return 1;
				}else if(m1.getAllAtoms() < m2.getAllAtoms()){
					return -1;
				} else if(m1.getAllBonds() > m2.getAllBonds()){
					return 1;
				} else if(m1.getAllBonds() < m2.getAllBonds()){
					return -1;
				}

				return 0;
			}
		};

	}

	public static boolean checkBiggestFragmentForUnwanted(StereoMolecule mol, List<Integer> liAtomicNo) {
		boolean bOk=true;

		ExtendedMolecule [] frags = mol.getFragments();
		int indexBiggestFrag = 0;
		if(frags.length > 1) {
			int maxAtoms = 0;
			for (int ii = 0; ii < frags.length; ii++) {
				if (frags[ii].getAllAtoms() > maxAtoms) {
					indexBiggestFrag = ii;
					maxAtoms = frags[ii].getAllAtoms();
				}
			}
		}

		ExtendedMolecule frag = frags[indexBiggestFrag];

		for (int i = 0; i < frag.getAllAtoms(); i++) {
			if(liAtomicNo.contains(new Integer(frag.getAtomicNo(i)))) {
				bOk=false;
				break;
			}
		}
		return bOk;
	}
	/**
	 *
	 * @param mol
	 * @param hsAtomicNo
	 * @return true if an atomic number from the hash set is found.
	 */
	public static boolean containsAtLeastOneAtomicNumbersFromHashSet(ExtendedMolecule mol, HashSet<Integer> hsAtomicNo) {
		boolean bOk=false;

		for (int i = 0; i < mol.getAllAtoms(); i++) {
			if(hsAtomicNo.contains(mol.getAtomicNo(i))) {
				bOk=true;
				break;
			}
		}
		return bOk;
	}

	public static boolean containsHeteroAtom(ExtendedMolecule mol, int [] arrIndexAt) {
		boolean hetero=false;

		for (int indexAt : arrIndexAt) {
			if(mol.getAtomicNo(indexAt)!=6 && mol.getAtomicNo(indexAt)!=1){
				hetero=true;
				break;
			}
		}
		return hetero;
	}

	/**
	 *
	 * @param mol
	 * @param hsAtomicNo
	 * @return true if the molecule contains an atomic number that is not in the hash set.
	 */
	public static boolean containsSolelyAtomicNumbersFromHashSet(ExtendedMolecule mol, HashSet<Integer> hsAtomicNo) {

		boolean allAtomicNosInHashset=true;

		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(!hsAtomicNo.contains(mol.getAtomicNo(i))) {

				allAtomicNosInHashset=false;

				break;
			}
		}

		return allAtomicNosInHashset;
	}

	/**
	 *
	 * @param molecule
	 * @param at index to specify atom in molecule
	 * @param fragment one atom must have the color: COLOR_FOR_CENTER_SELECTION.
	 * @return true if the specified atom in molecule equals the colored atom in fragment.
	 */
	public final static boolean atomAtomSubStrucMatch(StereoMolecule molecule, int at, StereoMolecule fragment) {
		boolean bMatch = false;

		SSSearcher sss = new SSSearcher();
		sss.setMol(fragment, molecule);
		int numFrags = sss.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cMatchAromDBondToDelocalized);
		if(numFrags == 0)
			return false;

		List<int []> liAtomLists = sss.getMatchList();
		ArrayUtilsCalc.removeDoubletsIntOrderIndepend(liAtomLists);

		int atomMarkedInQuery = -1;
		for (int i = 0; i < fragment.getAllAtoms(); i++) {
			if(fragment.getAtomColor(i)== COLOR_FOR_CENTER_SELECTION){
				atomMarkedInQuery = i;
			}
		}

		for (Iterator<int []> iter = liAtomLists.iterator(); iter.hasNext();) {
			int[] arrInd = (int[]) iter.next();
			if(arrInd[atomMarkedInQuery] == at){
				bMatch = true;
				break;
			}
		}

		return bMatch;

	}

	public static double [][] getDistanceArray(ExtendedMolecule mol) {

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


	public final static int getTopologicalDistance(ExtendedMolecule mol, int at1, int at2) {
		int dist = 0;

		if(at1 == at2)
			return 0;

		List<Integer> liExamine = new ArrayList<Integer>();
		List<Integer> liSphere = new ArrayList<Integer>();
		List<Integer> liVisited = new ArrayList<Integer>();

		liExamine.add(new Integer(at1));

		liSphere.add(new Integer(1));

		boolean bFound = false;
		while(!liExamine.isEmpty()) {
			dist = ((Integer)liSphere.remove(0)).intValue();
			int indAtCenter = ((Integer)liExamine.remove(0)).intValue();
			liVisited.add(new Integer(indAtCenter));
			int numJNeighbors = mol.getAllConnAtoms(indAtCenter);
			for (int i = 0; i < numJNeighbors; i++) {
				int indAtNeighb = mol.getConnAtom(indAtCenter, i);
				if(indAtNeighb == at2) {
					bFound = true;
					break;
				}
				if(!liVisited.contains(new Integer(indAtNeighb))) {
					liExamine.add(new Integer(indAtNeighb));
					liSphere.add(new Integer(dist + 1));
				}
			}
			if(bFound)
				break;
		}

		if(!bFound)
			dist = -1;

		return dist;
	}

	public final static int [][] getTopologicalDistanceMatrix(StereoMolecule mol) {
		return getNumberOfBondsBetweenAtoms(mol, mol.getBonds(), null);
	}

	/**
	 * From Joel Freyss developed for the FFMolecule
	 * Computes a matrix of distances between all the atoms in the graph.
	 * Complexity: O(m*n*maxBonds) m = number of bonds, n = number of atoms
	 * @param mol
	 * @param maxBonds
	 * @return an array A[i][j] = nBonds if i and j are connected by less than maxBonds
	 * 						      or -1  otherwise
	 */
	public static int[][] getNumberOfBondsBetweenAtoms(StereoMolecule mol, int maxBonds, int[][] dist) {
		//Initialization of the data
		if(dist==null)  dist = new int[mol.getAtoms()][mol.getAtoms()];
		int N = dist.length;
		for(int i=0; i<N; i++) {
			dist[i][i] = 0;
			for(int j=i+1; j<N; j++) {
				dist[i][j] = dist[j][i] = -1;
			}
		}

		//Algo: Dynamic Programming
		for(int j=0; j<maxBonds; j++) { //Maximum of nBonds bonds
			for(int i=0; i<mol.getBonds(); i++) {
				int a1 = mol.getBondAtom(0, i);
				int a2 = mol.getBondAtom(1, i);
				if(a1>=N || a2>=N) continue;
				for(int a0=0; a0<N; a0++) {
					// Dynamic Programming: dist(a0,a2) = min(dist(a0, a2), dist(a0,a1)+1) [if dist(a0,a1)>0]
					if(dist[a0][a1]>=0 && (dist[a0][a2]==-1 || dist[a0][a1]+1<dist[a0][a2]) && dist[a0][a1]<maxBonds) {
						dist[a2][a0] = (dist[a0][a2] = (dist[a0][a1] + 1));
					}
					if(dist[a0][a2]>=0 && (dist[a0][a1]==-1 || dist[a0][a2]+1<dist[a0][a1]) && dist[a0][a2]<maxBonds) {
						dist[a1][a0] = (dist[a0][a1] = (dist[a0][a2] + 1));
					}
				}
			}
		}
		return dist;
	}

	/**
	 *
	 * @param mol
	 * @param atm index of atom
	 * @return false if atom is not C or if one of the neighbours is N, O, F, S, P or Cl.
	 */
	public static boolean isAliphaticAtom(StereoMolecule mol, int atm) {
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

	public static boolean isAcceptor(StereoMolecule mol, int atom) {
		boolean acceptor = true;

		if(mol.getAtomicNo(atom)!=8 && mol.getAtomicNo(atom)!=7)
			return false;

		return acceptor;
	}

	public static boolean isDonor(StereoMolecule mol, int atom) {
		boolean donor = true;

		if(mol.getAtomicNo(atom)!=8 && mol.getAtomicNo(atom)!=7)
			return false;

		if(mol.getAllHydrogens(atom)==0){
			return false;
		}

		return donor;
	}


	public static boolean isCarbonTwoValencesMinimum(StereoMolecule mol) {

		boolean twoValencesMin = true;

		int atoms = mol.getAtoms();


		for (int i = 0; i < atoms; i++) {

			int atomicNo = mol.getAtomicNo(i);

			if(atomicNo==6){

				int connected = mol.getConnAtoms(i);

				if(connected < 2){
					twoValencesMin = false;
					break;
				}
			}
		}

		return twoValencesMin;
	}

	public static boolean isCarbonOnlyConnected2Hetero(StereoMolecule mol) {

		boolean hetero = true;

		int atoms = mol.getAtoms();

		outer:
		for (int i = 0; i < atoms; i++) {

			int atomicNo = mol.getAtomicNo(i);

			if(atomicNo==6){

				int connected = mol.getConnAtoms(i);

				for (int j = 0; j < connected; j++) {

					int indexConnecetd = mol.getConnAtom(i, j);

					if(mol.getAtomicNo(indexConnecetd)==6){
						hetero = false;
						break outer;
					}

				}

				if(mol.isRingAtom(i)){
					hetero = true;
					break;
				}
			}
		}

		return hetero;
	}

	/**
	 *
	 * @param mol
	 * @param atomIndex
	 * @return true if at least one neighbour is a hetero atom.
	 */
	public static boolean isCarbonConnected2Hetero(StereoMolecule mol, int atomIndex) {

		boolean hetero = false;

		int atomicNo = mol.getAtomicNo(atomIndex);

		if(atomicNo==6){

			int connected = mol.getConnAtoms(atomIndex);

			for (int j = 0; j < connected; j++) {

				int indexConnected = mol.getConnAtom(atomIndex, j);

				if(mol.getAtomicNo(indexConnected)!=6 && mol.getAtomicNo(indexConnected)!=1){
					hetero = true;
					break;
				}
			}
		}

		return hetero;
	}

	/**
	 *
	 * @param mol
	 * @param arrAtomIndex
	 * @return true if one of the atoms from arrAtomIndex is connected to a hetero atom.
	 */
	public static boolean isConnected2Hetero(StereoMolecule mol, int [] arrAtomIndex) {

		boolean hetero = false;

		for (int atomIndex : arrAtomIndex) {

			int connected = mol.getConnAtoms(atomIndex);

			for (int j = 0; j < connected; j++) {

				int indexConnected = mol.getConnAtom(atomIndex, j);

				if(mol.getAtomicNo(indexConnected)!=6 && mol.getAtomicNo(indexConnected)!=1){
					hetero = true;
					break;
				}
			}
		}
		return hetero;
	}

	public static boolean isRingInMolecule(StereoMolecule mol) {

		boolean ring = false;

		int atoms = mol.getAtoms();

		for (int i = 0; i < atoms; i++) {

			if(mol.isRingAtom(i)){
				ring = true;
				break;
			}
		}

		return ring;
	}

	/**
	 *
	 * @param mol
	 * @return true if all atoms are in a ring.
	 */
	public static boolean isRingExclusively(StereoMolecule mol) {

		boolean ring = true;

		int atoms = mol.getAtoms();

		for (int i = 0; i < atoms; i++) {

			if(!mol.isRingAtom(i)){
				ring = false;
				break;
			}
		}

		return ring;
	}

	public static boolean containsFiveBindingCarbon(StereoMolecule mol) {

		boolean five=false;
		for (int i = 0; i < mol.getAtoms(); i++) {
			if(mol.getAtomicNo(i)==6){
				int sumBO=0;
				for (int j = 0; j < mol.getConnAtoms(i); j++) {
					int bo = mol.getConnBondOrder(i,j);
					sumBO+=bo;
				}
				five=(sumBO>4)?true:false;
			}
		}

		return five;
	}

	/**
	 *
	 * @param mol
	 * @param atom
	 * @return true for cyano and iso-cyano.
	 */
	public static boolean isCyanoN(StereoMolecule mol, int atom) {

		if(mol.getAtomicNo(atom)!= 7){
			return false;
		}


		int nConn = mol.getConnAtoms(atom);

		if(nConn != 1){
			return false;
		}

		int indexConn = mol.getConnAtom(atom, 0);

		int bond = mol.getBond(atom, indexConn);

		if(mol.getBondOrder(bond) != 3){
			return false;
		}

		return true;
	}

	public static boolean isThioEther(StereoMolecule mol, int atom) {

		if(mol.getAtomicNo(atom)!= 16){
			return false;
		}


		int nConn = mol.getConnAtoms(atom);

		if(nConn != 2){
			return false;
		}

		boolean thio = true;
		for (int i = 0; i < nConn; i++) {

			int indexConn = mol.getConnAtom(atom, i);

			if(mol.getAtomicNo(indexConn) != 6){
				thio=false;
				break;
			}

		}

		return thio;
	}

	public static boolean isWildcard(StereoMolecule mol, int atom) {

		if(mol.getAtomicNo(atom) == PeriodicTable.ConnectionPoint){
			return true;
		}

		return false;
	}

	/**
	 *
	 * @param mol
	 * @param atom
	 * @return true if atom is S and at let one attached atom is O.
	 */
	public static boolean isSulfoxyGroup(StereoMolecule mol, int atom) {

		if(mol.getAtomicNo(atom)!= 16){
			return false;
		}

		int nConn = mol.getConnAtoms(atom);

		boolean oxy = false;
		for (int i = 0; i < nConn; i++) {

			int indexConn = mol.getConnAtom(atom, i);

			if(mol.getAtomicNo(indexConn) == 8){
				oxy=true;
				break;
			}

		}

		return oxy;
	}

	/**
	 *
	 * @param mol
	 * @param indexAtCentral
	 * @param arrIndexAt
	 * @return true if no connected carbon atom is in arrIndexAt
	 */
	public static boolean isIsolatedCarbon(StereoMolecule mol, int indexAtCentral, int [] arrIndexAt){

		boolean isolated=true;

		int nConnected = mol.getConnAtoms(indexAtCentral);

		boolean [] arrConnected = new boolean[mol.getAtoms()];

		for (int i = 0; i < nConnected; i++) {
			arrConnected[mol.getConnAtom(indexAtCentral,i)]=true;
		}

		for (int indexAt : arrIndexAt) {
			if(!arrConnected[indexAt]){
				continue;
			}

			if(mol.getAtomicNo(indexAt)==6){
				isolated=false;
				break;
			}
		}

		return isolated;
	}

	public static int [] extractAromaticRing(StereoMolecule mol, int [] arrIndexAt){

		RingCollection rc = mol.getRingSet();

		boolean [] arrRingMemberMarker = new boolean[mol.getAtoms()];

		int [] arrIndexAromaticRing = null;

		for (int i = 0; i < rc.getSize(); i++) {

			if(!rc.isAromatic(i)){
				continue;
			}

			Arrays.fill(arrRingMemberMarker, false);

			int [] arrRingAtoms = rc.getRingAtoms(i);

			for (int indexRingAtom : arrRingAtoms) {
				arrRingMemberMarker[indexRingAtom]=true;
			}

			int sum=0;
			for (int indexAt : arrIndexAt) {
				if(arrRingMemberMarker[indexAt]){
					sum++;
				}
			}

			if(sum==arrRingAtoms.length){
				arrIndexAromaticRing=arrRingAtoms;
				break;
			}
		}

		return arrIndexAromaticRing;
	}

	/**
	 * Counts cyano and iso-cyano
	 * @param mol
	 * @return
	 */
	public static int getNumCyanoGroups(StereoMolecule mol){

		int n = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(isCyanoN(mol, i)){
				n++;
			}
		}

		return n;

	}


	public static int getNumAlcoholicOxygen(StereoMolecule mol){

		int n = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(isAlcoholicOxygen(mol, i)){
				n++;
			}
		}

		return n;

	}

	public static int getNumThioEther(StereoMolecule mol){

		int n = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(isThioEther(mol, i)){
				n++;
			}
		}

		return n;

	}

	public static int getNumSulfOxyGroups(StereoMolecule mol){

		int n = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(isSulfoxyGroup(mol, i)){
				n++;
			}
		}

		return n;

	}

	public static int getNumWildcards(StereoMolecule mol){

		int n = 0;

		for (int i = 0; i < mol.getAllAtoms(); i++) {

			if(isWildcard(mol, i)){
				n++;
			}
		}

		return n;

	}

	public static boolean isAlcoholicOxygen(StereoMolecule mol, int atom) {

		if (mol.getAtomicNo(atom) != 8 || mol.isAromaticAtom(atom))
			return false;

		int nConnected = mol.getConnAtoms(atom);

		if(nConnected != 1){

			return false;
		}

		int indConn = mol.getConnAtom(atom, 0);

		if(mol.getAtomicNo(indConn) != 6){
			return false;
		}

		if(AtomFunctionAnalyzer.getNegativeNeighbourCount(mol, indConn) > 1){
			return false;
		}

		return true;
	}

	public static boolean isEtherOxygenAtAromatic(StereoMolecule mol, int atom) {

		boolean aromaticEtherO = false;

		if (mol.getAtomicNo(atom) != 8 || mol.isAromaticAtom(atom))
			return false;

		int nConnected = mol.getConnAtoms(atom);

		if(nConnected != 2){

			return false;
		}

		int indConn1 = mol.getConnAtom(atom, 0);

		int indConn2 = mol.getConnAtom(atom, 1);

		if(mol.isAromaticAtom(indConn1)){

			aromaticEtherO = true;

		} else if(mol.isAromaticAtom(indConn2)){

			aromaticEtherO = true;

		}

		return aromaticEtherO;
	}

	public static boolean isCarboxyC(StereoMolecule mol, int atom) {

		if (mol.getAtomicNo(atom) != 6 || mol.getAtomPi(atom) != 0)
			return false;

		boolean carboxyO = false;

		boolean etherO = false;

		for (int i=0; i<mol.getConnAtoms(atom); i++) {

			int conn = mol.getConnAtom(atom, i);

			if (mol.getAtomicNo(conn) == 6) {

				int bond = mol.getBond(atom, conn);

				int bo = mol.getBondOrder(bond);

				if(bo==1){
					etherO=true;
				} else if(bo==2){
					carboxyO=true;
				}

			}

		}

		boolean carboxyC=false;

		if(carboxyO && etherO){
			carboxyC=true;
		}

		return carboxyC;
	}

	public static int replaceAtoms(ExtendedMolecule [] arr, int atnoOrig,int atnoRpl) {
		int cc=0;
		for (int i = 0; i < arr.length; i++) {
			cc+=replaceAtoms(arr[i], atnoOrig,atnoRpl);
		}
		return cc;
	}

	public static int replaceAtoms(ExtendedMolecule mol, int atnoOrig,int atnoRpl) {
		int cc=0;
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i)==atnoOrig) {
				mol.setAtomicNo(i, atnoRpl);
				cc++;
			}
		}
		return cc;
	}





	/**
	 * Removes all molecules that are a substructure of one of the molecules in the input list.
	 * @param liInput
	 * @return
	 */
	public static final LinkedList<StereoMolecule> removeSubStructures(List<StereoMolecule> liInput){
		SSSearcher sss = new SSSearcher();

		LinkedList<StereoMolecule> li = new LinkedList<StereoMolecule>(liInput);

		for (int i = li.size()-1; i >= 0; i--) {

			StereoMolecule mol1 = li.get(i);

			for (int j = 0; j < li.size(); j++) {
				if(i!=j){

					StereoMolecule mol2 = li.get(j);

					if(mol1.getAtoms()<mol2.getAtoms()){
						sss.setMol(mol1, mol2);

						if(sss.isFragmentInMolecule()){
							li.remove(i);
							break;
						}
					}
				}
			}
		}

		return li;
	}

	/**
	 * Deletes the substructures, if found, in the molecule.
	 * Starts with the first structure in the list.
	 * Helper arrays are ensured.
	 * @param mol
	 * @param liFragment
	 * @return
	 */
	public static StereoMolecule removeSubstructuresFromMolecule(StereoMolecule mol, List<StereoMolecule> liFragment){

		StereoMolecule molReduced = new StereoMolecule(mol);

		molReduced.ensureHelperArrays(Molecule.cHelperRings);

		SSSearcher sss = new SSSearcher();

		HashSet<Integer> hsIndexMatchingAtoms = new HashSet<Integer>();

		for (StereoMolecule frag : liFragment) {

			sss.setMol(frag, molReduced);

			if(sss.findFragmentInMolecule()>0){

				List<int[]> liMatch = sss.getMatchList();

				hsIndexMatchingAtoms.clear();

				for (int[] arrMatch : liMatch) {

					for (int i = 0; i < arrMatch.length; i++) {
						hsIndexMatchingAtoms.add(arrMatch[i]);
					}
				}

				int [] arrIndexMatchingAtomsUnique = ArrayUtilsCalc.toIntArray(hsIndexMatchingAtoms);

				molReduced.deleteAtoms(arrIndexMatchingAtomsUnique);

				molReduced.ensureHelperArrays(Molecule.cHelperRings);

			}
		}


		return molReduced;
	}

	public static StereoMolecule removeSubstructureFromMolecule(StereoMolecule mol, StereoMolecule frag){

		StereoMolecule molReduced = new StereoMolecule(mol);

		molReduced.ensureHelperArrays(Molecule.cHelperRings);

		SSSearcher sss = new SSSearcher();



		sss.setMol(frag, molReduced);

		if(sss.findFragmentInMolecule()>0){

			List<int[]> liMatch = sss.getMatchList();


			if(liMatch.size() > 1){
				throw new RuntimeException("Fragment found more than once!");

			} if(liMatch.size() == 1){

				int [] arrMatchIndex = liMatch.get(0);

				molReduced.deleteAtoms(arrMatchIndex);

				molReduced.ensureHelperArrays(Molecule.cHelperRings);

//				boolean [] arrMatchingAtom = new boolean [molReduced.getAtoms()];
//
//
//
//				for (int i = 0; i < arrMatchIndex.length; i++) {
//					arrMatchingAtom[arrMatchIndex[i]]=true;
//				}
//
//				int indexRGroup = 1;
//
//				for (int i = 0; i < arrMatchIndex.length; i++) {
//
//					int atMatch = arrMatchIndex[i];
//
//					int nConnected = molReduced.getConnAtoms(atMatch);
//
//					for (int j = 0; j < nConnected; j++) {
//
//						int atConnected = molReduced.getConnAtom(atMatch, j);
//
//						if(arrMatchingAtom[atConnected]){
//							continue;
//						}
//
//						int bond = molReduced.getBond(atMatch, nConnected);
//
//						molReduced.deleteBond(bond);
//
//						int atNoRGroup = getAtomicNoRGroup(indexRGroup++);
//
//						int indexAtRGroup = molReduced.addAtom(atNoRGroup);
//
//						molReduced.addBond(atConnected, indexAtRGroup, Molecule.cBondTypeSingle);
//
//						molReduced.ensureHelperArrays(Molecule.cHelperNeighbours);
//
//					}
//
//				}
//


			}



		}



		return molReduced;
	}

	public static StereoMolecule removeWildcards(StereoMolecule mol){

		StereoMolecule molReduced = new StereoMolecule(mol);

		molReduced.ensureHelperArrays(Molecule.cHelperRings);

		for (int i = molReduced.getAtoms()-1; i >= 0; i--) {

			if(molReduced.getAtomicNo(i)==0){
				molReduced.deleteAtom(i);
			}

		}

		molReduced.ensureHelperArrays(Molecule.cHelperRings);

		return molReduced;
	}


	public static void setColorMCS2Molecule(StereoMolecule mol, StereoMolecule molMCS){

		for (int i = 0; i < mol.getAtoms(); i++) {
			mol.setAtomColor(i, Molecule.cAtomColorNone);
		}

		SSSearcher sss = new SSSearcher();

		sss.setMol(molMCS, mol);

		if(sss.findFragmentInMolecule()>0){
			List<int[]> liMatch = sss.getMatchList();

			int [] arrIndexMatch = liMatch.get(0);

			for (int i = 0; i < arrIndexMatch.length; i++) {
				mol.setAtomColor(arrIndexMatch[i], Molecule.cAtomColorBlue);
			}

		}

	}

	public static void setColorMolecule(StereoMolecule mol, int [] arrIndexMatch){

		for (int i = 0; i < mol.getAtoms(); i++) {
			mol.setAtomColor(i, Molecule.cAtomColorNone);
		}

		for (int i = 0; i < arrIndexMatch.length; i++) {
			mol.setAtomColor(arrIndexMatch[i], Molecule.cAtomColorBlue);
		}
	}

	/**
	 *
	 * @param mol
	 * @param arrIndexBonds
	 * @param color i.e. Molecule.cAtomColorBlue
	 */
	public static void setColorMoleculeFromBondIndex(StereoMolecule mol, int [] arrIndexBonds, int color){

		for (int i = 0; i < arrIndexBonds.length; i++) {

			final int indexAtom1 = mol.getBondAtom(0, arrIndexBonds[i]);

			final int indexAtom2 = mol.getBondAtom(1, arrIndexBonds[i]);

			mol.setAtomColor(indexAtom1, color);
			mol.setAtomColor(indexAtom2, color);

		}
	}

	public static void setCoordinatesNull(StereoMolecule mol){
		for (int i = 0; i < mol.getAtoms(); i++) {
			mol.setAtomX(i,0);
			mol.setAtomY(i,0);
			mol.setAtomZ(i,0);
		}
	}

	/**
	 * Returns the colors vector for the substructure in mol.
	 * @param mol
	 * @param frag
	 * @param atomColor Molecule.cAtomColor
	 * @return
	 */
	public static String getColorVectorSubstructure(StereoMolecule mol, StereoMolecule frag, int atomColor){

		SSSearcher sss = new SSSearcher();

		sss.setMol(frag, mol);

		HashSet<Integer> hsAtomIndex = new HashSet<Integer>();

		if(sss.findFragmentInMolecule()>0){
			List<int[]> liMatch = sss.getMatchList();


			int[] arr = liMatch.get(0);
			for (int i = 0; i < arr.length; i++) {
				hsAtomIndex.add(arr[i]);
			}

//			for (int[] arr : liMatch) {
//				for (int i = 0; i < arr.length; i++) {
//					hsAtomIndex.add(arr[i]);
//				}
//			}
		}

		String sAtomColor = ExtendedMoleculeFunctions.getColorRecord(mol, hsAtomIndex, atomColor);

		return sAtomColor;

	}

	public static float getSimilarity(StereoMolecule m1, StereoMolecule m2, DescriptorHandler dh){

		Object d1 = dh.createDescriptor(m1);
		Object d2 = dh.createDescriptor(m2);

		return dh.getSimilarity(d1, d2);
	}

	/**
	 * Taken from Thomas Sander SkeletonSpheres descriptor
	 * @param mol
	 * @param rootAtom atom index to start
	 * @param depth so many spheres are taken
	 * @return Fragment containing the spheres started at rootAtom
	 */
	public static StereoMolecule getSphere(StereoMolecule mol, int rootAtom, int depth){

		mol.ensureHelperArrays(Molecule.cHelperRings);
		StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());

		int[] atomList = new int[mol.getAtoms()];
		boolean[] atomMask = new boolean[mol.getAtoms()];
		if (rootAtom != 0)
			Arrays.fill(atomMask, false);

		int min = 0;
		int max = 0;

		for (int sphere=0; sphere<depth && max<mol.getAtoms(); sphere++) {
			if (max == 0) {
				atomList[0] = rootAtom;
				atomMask[rootAtom] = true;
				max = 1;
			}
			else {
				int newMax = max;
				for (int i=min; i<max; i++) {
					int atom = atomList[i];
					for (int j=0; j<mol.getConnAtoms(atom); j++) {
						int connAtom = mol.getConnAtom(atom, j);
						if (!atomMask[connAtom]) {
							atomMask[connAtom] = true;
							atomList[newMax++] = connAtom;
						}
					}
				}
				min = max;
				max = newMax;
			}

			mol.copyMoleculeByAtoms(fragment, atomMask, true, null);
		}

		return fragment;

	}


	/**
	 * Starts with 1 and goes until 16
	 * @return
	 */
	public static int getAtomicNoRGroup(int r){

		if(r > arrRGroupsAtomicNo.length){
			throw new RuntimeException("");
		}

		return arrRGroupsAtomicNo[r-1];

	}
}
