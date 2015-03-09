/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.chem;

public class AtomFunctionAnalyzer {
	private static int getOxoCount(StereoMolecule mol, int atom) {
		int count = 0;
		for (int i=0; i<mol.getConnAtoms(atom); i++)
			if (mol.getConnBondOrder(atom, i) == 2
			 && (mol.getAtomicNo(mol.getConnAtom(atom, i)) == 7
			  || mol.getAtomicNo(mol.getConnAtom(atom, i)) == 8
			  || mol.getAtomicNo(mol.getConnAtom(atom, i)) == 16))
				count++;

		return count;
		}

	public static int getNegativeNeighbourCount(StereoMolecule mol, int atom) {
		int count = 0;
		for (int i=0; i<mol.getConnAtoms(atom); i++)
			if (mol.isElectronegative(mol.getConnAtom(atom, i)))
				count++;

		return count;
		}

	public static boolean isAlkylAmine(StereoMolecule mol, int atom) {
		if (mol.getAtomicNo(atom) != 7 || mol.isAromaticAtom(atom) || mol.getAtomPi(atom) != 0)
			return false;

		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int conn = mol.getConnAtom(atom, i);
			if (mol.getAtomicNo(conn) != 6
			 || mol.getAtomPi(conn) != 0
			 || mol.isAromaticAtom(conn)
			 || getNegativeNeighbourCount(mol, conn) != 1) {
				return false;
				}
			}
		return true;
		}

	private static boolean isStabilized(StereoMolecule mol, int atom, boolean twice) {
		boolean aleadyFound = false;
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			if (!mol.isAromaticBond(mol.getConnBond(atom, i)) && mol.getConnBondOrder(atom, i) == 1) {
				int conn = mol.getConnAtom(atom, i);
				if (!mol.isAromaticAtom(conn)
				 && ((mol.getAtomicNo(conn) == 6 && getOxoCount(mol, conn) == 1)
				  || (mol.getAtomicNo(conn) == 16 && getOxoCount(mol, conn) == 2))) {
					if (aleadyFound || !twice)
						return true;
					aleadyFound = true;
					}
				}
			}

		return false;
		}

	public static boolean isAmide(StereoMolecule mol, int atom) {
		if (mol.getAtomicNo(atom) != 7 || mol.getAtomPi(atom) != 0)
			return false;

		boolean carbonylEquivalentFound = false;
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int conn = mol.getConnAtom(atom, i);
			if ((mol.getAtomicNo(conn) == 6 && getOxoCount(mol, conn) == 1)
			 || (mol.getAtomicNo(conn) == 16 && getOxoCount(mol, conn) == 2))
				carbonylEquivalentFound = true;
			else if (mol.getAtomicNo(conn) != 6)
				return false;
			}

		return carbonylEquivalentFound;
		}

	public static boolean isAmine(StereoMolecule mol, int atom) {
		if (mol.getAtomicNo(atom) != 7 || mol.isAromaticAtom(atom) || mol.getAtomPi(atom) != 0)
			return false;

		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int conn = mol.getConnAtom(atom, i);
			if (mol.getAtomicNo(conn) != 6
			 || (mol.getAtomPi(conn) != 0 && !mol.isAromaticAtom(conn))
			 || getNegativeNeighbourCount(mol, conn) != 1) {
				return false;
				}
			}

		return true;
		}

	public static boolean isArylAmine(StereoMolecule mol, int atom) {
		if (mol.getAtomicNo(atom) != 7 || mol.isAromaticAtom(atom) || mol.getAtomPi(atom) != 0)
			return false;

		boolean aromaticNeighbourFound = false;
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int conn = mol.getConnAtom(atom, i);
			if (mol.getAtomicNo(conn) != 6
			 || (mol.getAtomPi(conn) != 0 && !mol.isAromaticAtom(conn))
			 || getNegativeNeighbourCount(mol, conn) != 1)
				return false;
			if (mol.isAromaticAtom(conn))
				aromaticNeighbourFound = true;
			}

		return aromaticNeighbourFound;
	}
	
	public static boolean hasUnbalancedAtomCharge(StereoMolecule mol, int atom) {
		
		if (mol.getAtomCharge(atom)==0) {
			return false;
		} 
		
		boolean unbalanced=true;
		
		int chargeCenterAtom = mol.getAtomCharge(atom);
		
		int nConnected = mol.getConnAtoms(atom);
		
		int sumChargeConnectedAtoms = 0;
		
		for (int i = 0; i < nConnected; i++) {
			
			int indexAtom = mol.getConnAtom(atom, i);
			
			sumChargeConnectedAtoms += mol.getAtomCharge(indexAtom);
		}
		
		if(Math.abs(chargeCenterAtom) <= Math.abs(sumChargeConnectedAtoms)) {
			
			if(Math.signum(chargeCenterAtom) != Math.signum(sumChargeConnectedAtoms)) {
				unbalanced=false;
			}
			
		}
			
		return unbalanced;
	}

	public static boolean isAcidicOxygen(StereoMolecule mol, int atom) {
		
		boolean acidic=false;
		
		if (mol.getAtomicNo(atom) != 8)
			return false;
		
		if (mol.getConnAtoms(atom)!=1)
			return false;
		
		int indexConnected = mol.getConnAtom(atom, 0);
		
		// COOH
		if(mol.getAtomicNo(indexConnected)==6){
			int nConnected2C = mol.getConnAtoms(indexConnected);
			for (int i = 0; i < nConnected2C; i++) {
				int indexAtom = mol.getConnAtom(indexConnected, i);
				
				if(indexAtom==atom){
					continue;
				}
				
				if(mol.getAtomicNo(indexAtom) != 8){
					continue;
				}
				
				int indexBond = mol.getBond(indexConnected, indexAtom);
				
				if(mol.getBondType(indexBond)==Molecule.cBondTypeDouble){
					acidic=true;
					break;
				}
			}
		} else if(mol.getAtomicNo(indexConnected)==8){ // CSOOOH
			int nConnected2S = mol.getConnAtoms(indexConnected);
			
			int nDoubleBondedO2S=0;
			for (int i = 0; i < nConnected2S; i++) {
				int indexAtom = mol.getConnAtom(indexConnected, i);
				
				if(indexAtom==atom){
					continue;
				}
				
				if(mol.getAtomicNo(indexAtom) != 8){
					continue;
				}
				
				int indexBond = mol.getBond(indexConnected, indexAtom);
				
				if(mol.getBondType(indexBond)==Molecule.cBondTypeDouble){
					nDoubleBondedO2S++;
				}
			}
			
			if(nDoubleBondedO2S==2){
				acidic=true;
			}
			
		} else if(isAcidicOxygenAtPhosphoricAcid(mol, atom)){ // CP=O(OH)(OH)
			
			acidic=true;
				
		}
		
		return acidic;
	}
	
	public static boolean isAcidicOxygenAtPhosphoricAcid(StereoMolecule mol, int atom) {
		boolean acidic=false;
		
		if (mol.getAtomicNo(atom) != 8)
			return false;
		
		if (mol.getConnAtoms(atom)!=1)
			return false;
		
		int indexConnected = mol.getConnAtom(atom, 0);
		
		if(mol.getAtomicNo(indexConnected)==15){ // CP=O(OH)(OH)
			int nConnected2P = mol.getConnAtoms(indexConnected);
			
			for (int i = 0; i < nConnected2P; i++) {
				int indexAtom = mol.getConnAtom(indexConnected, i);
				
				if(indexAtom==atom){
					continue;
				}
				
				if(mol.getAtomicNo(indexAtom) != 8){
					continue;
				}
				
				int indexBond = mol.getBond(indexConnected, indexAtom);
				
				if(mol.getBondType(indexBond)==Molecule.cBondTypeDouble){
					acidic=true;
					break;
				}
			}
		}
		
		return acidic;
	}
	
	
	public static boolean isMemberOfNitroGroup(StereoMolecule mol, int atom) {
		
		boolean member = false;
		
		if ((mol.getAtomicNo(atom) != 7) && (mol.getAtomicNo(atom) != 8))
			return false;
		
		
		if (mol.getAtomicNo(atom) == 7) {
			
			if(isNitroGroupN(mol, atom)) {
				member = true;
			}
		} else if (mol.getAtomicNo(atom) == 8) {
			
			int nConnAts = mol.getConnAtoms(atom);
			for (int i = 0; i < nConnAts; i++) {
				
				int indexAt = mol.getConnAtom(atom, i);
				
				if (mol.getAtomicNo(indexAt) == 7) {
					if(isNitroGroupN(mol, indexAt)) {
						member = true;
						break;
					}
				}
			}
		}
			
		return member;
	}
	
	
	public static boolean isNitroGroupN(StereoMolecule mol, int atom) {
		
		boolean nitro = false;
		
		if ((mol.getAtomicNo(atom) != 7))
			return false;
		
		int nConnAts = mol.getConnAtoms(atom);
		
		int indexSingleBondedO = -1;
		
		int indexDoubleBondedO = -1;
		
		for (int i = 0; i < nConnAts; i++) {
			int indexAt = mol.getConnAtom(atom, i);
			
			if(mol.getAtomicNo(indexAt)==8){
				int indexBnd = mol.getBond(atom, indexAt);
				
				if(mol.getBondOrder(indexBnd)==1) {
					indexSingleBondedO = indexAt;
				} else if(mol.getBondOrder(indexBnd)==2) {
					indexDoubleBondedO = indexAt;
				}
			}
		}
		
		if((indexSingleBondedO>-1)&&(indexDoubleBondedO>-1)) {
			nitro = true;
		}
			
		return nitro;
	}

	
	
	public static boolean isBasicNitrogen(StereoMolecule mol, int atom) {
		if (mol.getAtomicNo(atom) != 7)
			return false;

		if (mol.getConnAtoms(atom) + mol.getAtomPi(atom) > 3)
			return false;

		if (mol.isAromaticAtom(atom)) {
			if (mol.getAtomPi(atom) != 1)
				return false;	// pyrrol type

			if (mol.getAtomRingCount(atom, 7) != 1)
				return false;

			RingCollection rc = mol.getRingSet();
			for (int r=0; r<rc.getSize(); r++) {
				if (rc.isAtomMember(r, atom)) {
					if (rc.getRingSize(r) == 5 || rc.getRingSize(r) == 6) {
						int[] ring = rc.getRingAtoms(r);

						int nIndex = -1;
						for (int i=0; i<ring.length; i++) {
							if (ring[i] == atom) {
								nIndex = i;
								break;
								}
							}
	
						int enablerCount = 0;

						int[] opi = null;	// ortho,para influences
						int[] mi = null;	// meta influences
						if (ring.length == 5) {
							opi = new int[2];
							opi[0] = ring[(nIndex-1 < 0) ? nIndex+4 : nIndex-1];
							opi[1] = ring[(nIndex-4 < 0) ? nIndex+1 : nIndex-4];
							mi = new int[2];
							mi[0] = ring[(nIndex-2 < 0) ? nIndex+3 : nIndex-2];
							mi[1] = ring[(nIndex-3 < 0) ? nIndex+2 : nIndex-3];
							}

						if (ring.length == 6) {
							opi = new int[3];
							opi[0] = ring[(nIndex-1 < 0) ? nIndex+5 : nIndex-1];
							opi[1] = ring[(nIndex-3 < 0) ? nIndex+3 : nIndex-3];
							opi[2] = ring[(nIndex-5 < 0) ? nIndex+1 : nIndex-5];
							mi = new int[2];
							mi[0] = ring[(nIndex-2 < 0) ? nIndex+4 : nIndex-2];
							mi[1] = ring[(nIndex-4 < 0) ? nIndex+2 : nIndex-4];
							}

						for (int i=0; i<ring.length; i++)
							if (atom != ring[i] && mol.getAtomicNo(ring[i]) == 7 && mol.getAtomPi(ring[i]) == 1)
								enablerCount--;

						for (int i=0; i<opi.length; i++) {
							int exoCyclicAtom = -1;
							int exoCyclicBond = -1;
							for (int j=0; j<mol.getConnAtoms(opi[i]); j++) {
								if (!mol.isAromaticBond(mol.getConnBond(opi[i], j))) {
									exoCyclicAtom = mol.getConnAtom(opi[i], j);
									exoCyclicBond = mol.getConnBond(opi[i], j);
									break;
									}
								}

							if (exoCyclicAtom != -1) {
								if (mol.getAtomicNo(exoCyclicAtom) == 7
								 && mol.getAtomPi(exoCyclicAtom) == 0
								 && (mol.getConnAtoms(exoCyclicAtom) + mol.getAtomPi(exoCyclicAtom) <= 3)
								 && !isStabilized(mol, exoCyclicAtom, false)) {
									enablerCount++;
									continue;
									}
								if (mol.getAtomicNo(exoCyclicAtom) == 8
								 && mol.getConnAtoms(exoCyclicAtom) == 1) {
									enablerCount += 2;
									continue;
									}
								if (mol.isAromaticBond(exoCyclicBond)) {
									for (int s=0; s<rc.getSize(); s++) {
										if (rc.isAromatic(s) && rc.isAtomMember(s, exoCyclicAtom)) {
											int[] ratom = rc.getRingAtoms(s);
											for (int j=0; j<ratom.length; j++) {
												if (mol.getAtomicNo(ratom[j]) == 7 && mol.getAtomPi(ratom[j]) == 1) {
													enablerCount--;
													break;
													}
												}
											break;
											}
										}
									}
								}
							}
						for (int i=0; i<mi.length; i++) {
							int exoCyclicAtom = -1;	// find exocyclic atom
							for (int j=0; j<mol.getConnAtoms(mi[i]); j++)
								if (!mol.isAromaticBond(mol.getConnBond(mi[i], j)))
									exoCyclicAtom = mol.getConnAtom(mi[i], j);

							if (mol.getAtomicNo(mi[i]) == 7) {
								if (mol.getAtomPi(mi[i]) == 0 && (exoCyclicAtom == -1 || getOxoCount(mol, exoCyclicAtom) == 0))	// imidazol type
									enablerCount++;
								continue;
								}
							if (exoCyclicAtom != -1 && getOxoCount(mol, exoCyclicAtom) != 0) {
								enablerCount--;
								continue;
								}
							}
						return enablerCount > 0;
						}
					break;
					}
				}
			return false;
			}

		if (mol.getAtomPi(atom) > 1)
			return false;	// nitrile

		if (mol.getAtomPi(atom) == 1) {
			int imineC = -1;
			int supporterCount = 0;
			for (int i=0; i<mol.getConnAtoms(atom); i++) {
				int conn = mol.getConnAtom(atom, i);
				if (mol.getConnBondOrder(atom, i) == 2) {
					if (mol.getAtomicNo(conn) != 6)
						return false;	// N=O, N=N, N=S, etc.
					imineC = conn;
					continue;
					}
				if (mol.getAtomicNo(conn) == 8)	// C=N-O
					return false;
				if (mol.getAtomicNo(conn) == 7) {
					supporterCount--;
					if (isStabilized(mol, conn, false))	// C=N-N-C=O
						supporterCount--;
					continue;
					}
				if (mol.isAromaticAtom(conn))
					supporterCount--;
				}

			if (imineC == -1)
				return false;

			int aromaticNeighborCount = 0;
			for (int i=0; i<mol.getConnAtoms(imineC); i++) {
				if (mol.getConnBondOrder(imineC, i) == 1) {
					int conn = mol.getConnAtom(imineC, i);
					if (getOxoCount(mol, conn) != 0)
						return false;

					if (mol.isAromaticAtom(conn))
						aromaticNeighborCount++;
					if (mol.getAtomicNo(conn) == 7 && !isStabilized(mol, conn, true))
						supporterCount++;
					if (mol.getAtomicNo(conn) == 8 || mol.getAtomicNo(conn) == 16)
						supporterCount--;	// S-C=N or O-C=N mostly below pKa=7
					}
				}
			if (aromaticNeighborCount == 2)	// two pheny substituents reduce pKa to 6
				supporterCount--;

			return (supporterCount >= 0);
			}

		// non-aromatic nitrogen without pi-bond
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int conn = mol.getConnAtom(atom, i);
			if (mol.isAromaticAtom(conn))
				return false;

			if (mol.getAtomicNo(conn) != 6)
				return false;

			if (getOxoCount(mol, conn) != 0)
				return false;

			if (mol.getAtomPi(conn) != 0 && isVinylogOxo(mol, conn))
				return false;
			}

		return true;
		}

	private static boolean isVinylogOxo(StereoMolecule mol, int atom) {
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			if (mol.getConnBondOrder(atom, i) != 1) {
				int conn = mol.getConnAtom(atom, i);
				for (int j=0; j<mol.getConnAtoms(conn); j++)
					if (mol.getConnBondOrder(conn, j) == 1
					 && getOxoCount(mol, mol.getConnAtom(conn, j)) != 0)
						return true;
				}
			}
		return false;
		}
	}
