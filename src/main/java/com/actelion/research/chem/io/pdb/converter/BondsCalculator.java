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
 * @author Joel Freyss
 */

package com.actelion.research.chem.io.pdb.converter;

import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.VDWRadii;
import com.actelion.research.util.IntQueue;

import java.util.*;

/**
 * BondsCalculator is used to recreate the bonds and / or calculate the bonds orders 
 * based on the 3D coordinates of the atoms
 */
public class BondsCalculator {
	/**
	 * Calculates the bonds of a molecule by checking the distance between
	 * all atoms. The bond order is not set with this function.
	 * 
	 * Complexity O(nAtoms)
	 * Memory O(nAtoms)
	 * 
	 * @param mol
	 * @param lenient
	 * @param atomToGroup
	 * @throws Exception
	 */
	public static void createBonds(StereoMolecule mol, boolean lenient, Map<Integer,String> atomToGroup) throws Exception {
		if(mol.getAllAtoms()==0) return;
		boolean displayWarning = true;

		//1. Create a grid
		MoleculeGrid grid = new MoleculeGrid(mol);
		TreeSet<Integer> atomsToRemove = new TreeSet<Integer>();
	 	List<int[]> potentialBonds = new ArrayList<int[]>();
		int[] neighborCount = new int[mol.getAllAtoms()];
	 	
		//2. For each atom, check the neighbours and
		//   Create a connection if the distance is close to the sum of VDW
		for(int i=0; i<mol.getAllAtoms(); i++) {			
			if(atomsToRemove.contains(i)) continue;
			if(!mol.isOrganicAtom(i)) continue;
			
			//Get the neighbours
			Set<Integer> set = grid.getNeighbours(getCoordinates(mol, i), 3.2, false);
			for(int j:set) {
				if(i>=j || atomsToRemove.contains(j)) continue;
				if(!mol.isOrganicAtom(j)) continue;
				
				double dist = Math.sqrt(mol.getCoordinates(i).distanceSquared(mol.getCoordinates(j)));
				double idealDist = VDWRadii.COVALENT_RADIUS[mol.getAtomicNo(i)] + VDWRadii.COVALENT_RADIUS[mol.getAtomicNo(j)];
				if(atomToGroup!=null) {
					if(!match(atomToGroup.get(i), atomToGroup.get(j))
						|| (mol.getAllAtoms()>200 && ((j-i)>12 && (j-i)>mol.getAllAtoms()/50))) {
						if(dist>idealDist + .45) continue;
						potentialBonds.add(new int[]{i, j});
						continue;
					}
				}				
				else if(dist>idealDist + .45) continue;


				if(neighborCount[i] >= maxNeighborCount(mol, i)
				|| neighborCount[j] >= maxNeighborCount(mol, j)) {
					if(!lenient)
						throw new Exception("Valence exceeded "+i+" or "+j);
					
					if(!displayWarning) {
						System.err.println("Valence exceeded "+i+" or "+j);
						displayWarning = false;
					} 						
					continue;
				}
				try {
					mol.addBond(i, j, 1);
					neighborCount[i]++;
					neighborCount[j]++;
				} catch (Exception e) {
					if(!lenient) throw e;
				}
			}

		}
		//System.out.println(potentialBonds.size()+" potential bonds");
		if(potentialBonds.size()<mol.getAllAtoms()/30) {
			for (int[] pair: potentialBonds) {
				try {
					mol.addBond(pair[0], pair[1], 1);
					neighborCount[pair[0]]++;
					neighborCount[pair[1]]++;
				} catch (Exception e) {
					if(!lenient) throw e;
				}				
			}
		}
		
		if(atomsToRemove.size()>0) {
			if(!lenient && atomsToRemove.size()>4) throw new Exception(atomsToRemove.size()+" atoms in close proximity");
			System.err.println(atomsToRemove.size()+" atoms in too close proximity");

			for (int atom:atomsToRemove)
				mol.markAtomForDeletion(atom);
			mol.deleteMarkedAtomsAndBonds();
		}
		//System.out.println("Create bonds in "+(System.currentTimeMillis()-s)+"ms");
	}

	/**
	 * Calculates for the organic subset of atoms the maximum number of neighbors
	 * @param mol
	 * @param atom
	 * @return
	 */
	private static int maxNeighborCount(StereoMolecule mol, int atom) {
		int atomicNo = mol.getAtomicNo(atom);
		return atomicNo == 5 ? 6    // B
			 : atomicNo <= 7 ? 4    // C
			 : atomicNo == 8 ? 2    // N
			 : atomicNo == 9 ? 1    // O
			 : atomicNo == 14 ? 4   // Si
			 : atomicNo == 15 ? 4   // P
			 : atomicNo == 16 ? 4   // S
			 : atomicNo == 17 ? 1   // Cl
			 : atomicNo == 33 ? 5   // As
			 : atomicNo == 34 ? 6   // Se
			 : atomicNo == 35 ? 6   // Br
			 : atomicNo == 52 ? 6
			 : atomicNo == 53 ? 6 : 8;
	}

	private static Coordinates getCoordinates(StereoMolecule mol, int atom) {
		return new Coordinates(mol.getAtomX(atom), mol.getAtomY(atom), mol.getAtomZ(atom));
	}

	private static boolean match(String g1, String g2) {
		return g1.equals(g2);
		/*
		String s1[] = g1.split(" ");
		String s2[] = g2.split(" ");
		for (String ss1 : s1) {
			for (String ss2 : s2) {
				if(ss1.equals(ss2)) return true;
			}			
		}
		return false;
		*/
	}
	
	/**
	 * Calculate the bond orders of the molecule (without knowing the hydrogens).
	 * The calculation is based on the bond distance between each atoms.
	 * 
	 * 
	 * http://www.ccp14.ac.uk/ccp/web-mirrors/i_d_brown/valence.txt
	 * s = exp((Ro - R)/B)
	 *
	 * @param mol
	 */
	public static void calculateBondOrders(StereoMolecule mol, boolean lenient) throws Exception {
		boolean[] visited = new boolean[mol.getAllBonds()];
		mol.ensureHelperArrays(Molecule.cHelperRings);

		//Hybridization State Determination
		int[] spOrder = new int[mol.getAllAtoms()];
		for(int atom=0; atom<mol.getAllAtoms(); atom++) {
			
			if(mol.getConnAtoms(atom)<=1) {
				spOrder[atom] = 1;
			} else if(mol.getConnAtoms(atom)==2) {
				double angle = GeometryCalculator.getAngle(mol, mol.getConnAtom(atom, 0), atom, mol.getConnAtom(atom, 1));
				if(Math.abs(angle-Math.PI)<Math.PI/6) spOrder[atom] = 1;
				else spOrder[atom] = 2;
			} else if(mol.getConnAtoms(atom)==3) {
				Coordinates c = mol.getCoordinates(atom);
				Coordinates u = c.subC(mol.getCoordinates(mol.getConnAtom(atom, 0)));
				Coordinates v = c.subC(mol.getCoordinates(mol.getConnAtom(atom, 1)));
				Coordinates w = c.subC(mol.getCoordinates(mol.getConnAtom(atom, 2)));
				Coordinates normal = u.cross(v);
				if(normal.distSq()>0) { 
					double proj = normal.unitC().dot(w) / w.dist();
					if(Math.abs(proj)<0.3) spOrder[atom] = 2;
					else spOrder[atom] = 3;
				} else {
					spOrder[atom] = 3;
				}
			}
		}
		
		//////////////////////////////////
		// Functional Group Recognition
		
		for(int atom=0; atom<mol.getAllAtoms(); atom++) {
			if(mol.getAllConnAtoms(atom)==2) {
				//RN=N=N
				
				
			} else if(mol.getAllConnAtoms(atom)==3) {
				int a, b;
				int a1 = mol.getConnAtom(atom, 0);
				int a2 = mol.getConnAtom(atom, 1);
				int a3 = mol.getConnAtom(atom, 2);
				if(mol.getAtomicNo(atom)==6 && spOrder[atom]==2) {
					if(mol.getAtomRingSize(atom)>0) continue;

					//C(R)(O)(=O)
					if( (mol.getAtomicNo(a2)==8 && mol.getAtomicNo(a3)==8 && mol.getAllConnAtoms(a2)==1 && mol.getAllConnAtoms(a3)==1) ||
						(mol.getAtomicNo(a1)==8 && mol.getAtomicNo(a3)==8 && mol.getAllConnAtoms(a1)==1 && mol.getAllConnAtoms(a3)==1) ||
						(mol.getAtomicNo(a1)==8 && mol.getAtomicNo(a2)==8 && mol.getAllConnAtoms(a1)==1 && mol.getAllConnAtoms(a2)==1)) {
							mol.setBondOrder(shortestBond(mol, atom, 8, false), 2); continue;
					}										
					
					//C(R)(OR)(=O)
					a = connectedAtom(mol, atom, 8, 2, 0, 0);
					b = connectedBond(mol, atom, 8, 1);
					if(a>=0 && b>=0) {mol.setBondOrder(b, 2); continue;} 
					
					//C(R)(SR)(=O)
					a = connectedAtom(mol, atom, 16, 2, 0, 0);
					b = connectedBond(mol, atom, 8, 1);
					if(a>=0 && b>=0) { mol.setBondOrder(b, 2); continue;} 
					
					//C(R)(NR)(=O)
					a = connectedAtom(mol, atom, 7, 2, 0, 0);
					b = connectedBond(mol, atom, 8, 1);
					if(a>=0 && b>=0) {mol.setBondOrder(b, 2); continue;} 
					
					//C(R)(SR)(=S)
					a = connectedAtom(mol, atom, 16, 2, 0, 0);
					b = connectedBond(mol, atom, 16, 1);
					if(a>=0 && b>=0) {mol.setBondOrder(b, 2); continue;} 
					
					//C(R)(NR)(=S)
					a = connectedAtom(mol, atom, 7, 2, 0, 0);
					b = connectedBond(mol, atom, 16, 1);
					if(a>=0 && b>=0) {mol.setBondOrder(b, 2); continue;}
					
					
					//C(CR)(N)(=N)
					if((mol.getAtomicNo(a1)==6 && mol.getAtomicNo(a2)==7 && mol.getAtomicNo(a3)==7 && mol.getAllConnAtoms(a2)==1 && mol.getAllConnAtoms(a3)==1) ||
						(mol.getAtomicNo(a1)==7 && mol.getAtomicNo(a2)==6 && mol.getAtomicNo(a3)==7 && mol.getAllConnAtoms(a1)==1 && mol.getAllConnAtoms(a3)==1) ||
						(mol.getAtomicNo(a1)==7 && mol.getAtomicNo(a2)==7 && mol.getAtomicNo(a3)==6 && mol.getAllConnAtoms(a1)==1 && mol.getAllConnAtoms(a2)==1)) {
							mol.setBondOrder(shortestBond(mol, atom, 7, true), 2); continue;
					}				
					//C(NR)(N)(=N) -> Arginin
					if(mol.getAtomicNo(a1)==7 && mol.getAtomicNo(a2)==7 && mol.getAtomicNo(a3)==7) {
						if(mol.getConnAtoms(a1)==2 && mol.getConnAtoms(a2)==1 && mol.getConnAtoms(a3)==1) {
							if(mol.getCoordinates(atom).distSquareTo(mol.getCoordinates(a2))<mol.getCoordinates(atom).distSquareTo(mol.getCoordinates(a3))) {mol.setBondOrder(mol.getConnBond(atom, 1), 2); continue;}
							else {mol.setBondOrder(mol.getConnBond(atom, 2), 2); continue;}
						} else if(mol.getConnAtoms(a1)==1 && mol.getConnAtoms(a2)==2 && mol.getConnAtoms(a3)==1) {
							if(mol.getCoordinates(atom).distSquareTo(mol.getCoordinates(a1))<mol.getCoordinates(atom).distSquareTo(mol.getCoordinates(a3))) {mol.setBondOrder(mol.getConnBond(atom, 0), 2); continue;}
							else {mol.setBondOrder(mol.getConnBond(atom, 2), 2); continue;}
						} else if(mol.getConnAtoms(a1)==1 && mol.getConnAtoms(a2)==1 && mol.getConnAtoms(a3)==2) {
							if(mol.getCoordinates(atom).distSquareTo(mol.getCoordinates(a1))<mol.getCoordinates(atom).distSquareTo(mol.getCoordinates(a2))) {mol.setBondOrder(mol.getConnBond(atom, 0), 2); continue;}
							else {mol.setBondOrder(mol.getConnBond(atom, 1), 2); continue;}
						}
					}								
					/*
					//C(NR)(NR)(=NR)
					if(mol.getAtomicNo(a1)==7 && mol.getAtomicNo(a2)==7 && mol.getAtomicNo(a3)==7 ) {
						
						mol.setBondOrder(shortestBond(mol, i, 7, true), 2); continue;
					}
					*/

				} else if(mol.getAtomicNo(atom)==7) {
					//N(R)(R)C=O -> Amide
					a = connectedAtom(mol, atom, 6, 2, 8, 1);
					b = connectedBond(mol, a, 8, 1);
					if(a>=0 && b>=0) {mol.setBondOrder(b, 2); continue;} 					

					//N(=O)(=O) -> Nitro
					for (int j = 0; j < mol.getAllConnAtoms(atom); j++) {
						if(mol.getAtomicNo(a1)==8 && mol.getAllConnAtoms(a1)==1 && mol.getAtomicNo(a2)==8 && mol.getAllConnAtoms(a2)==1) {
							mol.setBondOrder(mol.getConnBond(atom,0), 2);
							mol.setBondOrder(mol.getConnBond(atom,1), 2);
						} else if(mol.getAtomicNo(a1)==8 && mol.getAllConnAtoms(a1)==1 && mol.getAtomicNo(a3)==8 && mol.getAllConnAtoms(a3)==1) {
							mol.setBondOrder(mol.getConnBond(atom,0), 2);
							mol.setBondOrder(mol.getConnBond(atom,2), 2);
						} else if(mol.getAtomicNo(a2)==8 && mol.getAllConnAtoms(a2)==1 && mol.getAtomicNo(a3)==8 && mol.getAllConnAtoms(a3)==1) {
							mol.setBondOrder(mol.getConnBond(atom,1), 2);
							mol.setBondOrder(mol.getConnBond(atom,2), 2);
						} 
					}
					
					if(a>=0 && b>=0) {mol.setBondOrder(b, 2); continue;} 					
					
				} 
			} else if(mol.getAllConnAtoms(atom)==4) {
				if(mol.getAtomicNo(atom)==16) {
					int count = 0;
					for(int j=0; count<2 && j<mol.getAllConnAtoms(atom); j++) {
						if(mol.getAtomicNo(mol.getConnAtom(atom, j))==8 && mol.getAllConnAtoms(mol.getConnAtom(atom, j))==1) {
							mol.setBondOrder(mol.getConnBond(atom, j), 2);
							count++;
						}
					}
					for(int j=0; count<2 && j<mol.getAllConnAtoms(atom); j++) {
						if(mol.getAtomicNo(mol.getConnAtom(atom, j))==7 && mol.getAllConnAtoms(mol.getConnAtom(atom, j))==1) {
							mol.setBondOrder(mol.getConnBond(atom, j), 2);
							count++;
						}
					}
				} else if(mol.getAtomicNo(atom)==15) {
/*					int b = shortestBond(mol, i, 8, false);
					if( b>=0) {
						if((mol.getBondAtom(0, b)==i && mol.getAllConnAtoms(mol.getBondAtom(1, b))==1) ||
							(mol.getBondAtom(1, b)==i && mol.getAllConnAtoms(mol.getBondAtom(0, b))==1)) { 
								mol.setBondOrder(b, 2);
						}
					}*/
				}
			}
		}
		//Preliminary pass: process obvious bonds outside rings		
		for (int bond = 0; bond < mol.getAllBonds(); bond++) {
			int a1 = mol.getBondAtom(0, bond);
			int a2 = mol.getBondAtom(1, bond);

//			if(atomToRings[a1].size()>0 || atomToRings[a2].size()>0) continue;
			if (mol.isRingBond(bond)) continue;    // instead

			if(mol.getImplicitHydrogens(a1)==0 || mol.getImplicitHydrogens(a2)==0) continue;
			if(!isPlanar(mol, a1, a2)) continue;

			double order = getLikelyOrder(mol, a1, a2);
			if(order>3.0 && spOrder[a1]==1 && spOrder[a2]==1 && mol.getImplicitHydrogens(a1)>=2 && mol.getImplicitHydrogens(a2)>=2) {
				mol.setBondOrder(bond, 3);
				visited[bond] = true;
			} else if(order>2.6 && spOrder[a1]<=2 && spOrder[a2]<=2 ) {
				mol.setBondOrder(bond, 2);
				visited[bond] = true;
			}
			
		}

		/////////////////////////////////////////////////////////
		// Aromatic Ring Perception
		// This procedure calculates a normal to the ring and check that all 
		// atoms in the ring and their neighbours are within the plane defined by the normal
		mol.ensureHelperArrays(Molecule.cHelperRingsSimple);
		RingCollection ringSet = mol.getRingSetSimple();
		ArrayList<Integer>[] atomToRings = getAtomToRings(mol);
		boolean[] aromaticRing = new boolean[ringSet.getSize()];


		//int[] pyroles = new int[allRings.size()];
		//int[] oxo = new int[allRings.size()];
		//int[] toDo = new int[mol.getAllAtoms()];

		for (int size = 5; size <= 6; size++)
		for (int ringNo = 0; ringNo < ringSet.getSize(); ringNo++) {
			int[] ring = ringSet.getRingAtoms(ringNo);
			if(ring.length!=size) continue;  
						
			Coordinates c0 = mol.getCoordinates(ring[0]);
			Coordinates c1 = mol.getCoordinates(ring[1]);
			Coordinates c2 = mol.getCoordinates(ring[2]);
			Coordinates normal = c1.subC(c0).cross(c1.subC(c2));
			if(normal.distSq()==0) continue;
			//int startIndex = 0;
			boolean planar = true;
			
			for(int i=0; i<ring.length && planar; i++) {
				Coordinates c3 = mol.getCoordinates(ring[i]);
				Coordinates w = c1.subC(c3);
				if(Math.abs(normal.unitC().dot(w) / w.dist())>0.3) {planar=false;}
				
				//Make sure that all carbon in the ring are planar
				if(mol.getAtomicNo(ring[i])==6 || mol.getAtomicNo(ring[i])==7) {
					if(spOrder[ring[i]]!=2) {planar=false;}
				} else if(mol.getAtomicNo(ring[i])>16) {
					planar=false;  //continue if the ring has a metal atom
				}
			}			
			if(!planar) continue;

			//
			//Special case 1: Histidine (some obfuscated code in order to avoid a SS search)
			if(ring.length==5) {
				//Central C:			
				int start = -1;
				int[] posN = {-1, -1};
				boolean ok = true;
				for(int i=0; ok && i<ring.length; i++) {
					if(mol.getAtomicNo(ring[i])==6 && mol.getAllConnAtoms(ring[i])==3) {start = i;}
					else if(mol.getAllConnAtoms(ring[i])!=2) ok = false;
					else if(mol.getAtomicNo(ring[i])==7) {
						if(posN[0]<0) posN[0] = i;
						else if(posN[1]<0) posN[1] = i;
						else ok = false;
					}
				}
				if(ok && start>=0 && posN[1]>=0) {
					if((start+2)%5==posN[0] && (start+4)%5==posN[1]) {
						mol.setBondOrder(mol.getBond(ring[start], ring[(start+1)%5]), 2);
						mol.setBondOrder(mol.getBond(ring[(start+3)%5], ring[(start+4)%5]), 2);
						continue;	
					} else if((start+2)%5==posN[1] && (start+4)%5==posN[0]) {
						mol.setBondOrder(mol.getBond(ring[start], ring[(start+1)%5]), 2);
						mol.setBondOrder(mol.getBond(ring[(start+3)%5], ring[(start+4)%5]), 2);
						continue;	
					} else if((start+3)%5==posN[0] && (start+1)%5==posN[1]) {
						mol.setBondOrder(mol.getBond(ring[start], ring[(start+4)%5]), 2);
						mol.setBondOrder(mol.getBond(ring[(start+1)%5], ring[(start+2)%5]), 2);
						continue;	
					} else if((start+3)%5==posN[1] && (start+1)%5==posN[0]) {
						mol.setBondOrder(mol.getBond(ring[start], ring[(start+4)%5]), 2);
						mol.setBondOrder(mol.getBond(ring[(start+1)%5], ring[(start+2)%5]), 2);
						continue;	
					}
					
				}
				
			}

			// 
			//Check Huckel's rule and Find the starting position
			int start = -1;
			int nElectrons = 0;
			int nAmbiguousN = 0;
			int nAmbiguousC = 0;
			for(int i=0; i<ring.length; i++) {
				int a1 = ring[(i)%ring.length];				
				int a2 = ring[(i+1)%ring.length];				
				int a0 = ring[(i-1+ring.length)%ring.length];				
				int bnd1 = mol.getBond(a1, a2);
				int bnd2 = mol.getBond(a1, a0);
				if(mol.getAtomicNo(a1)==6) {
					if(mol.getAllConnAtoms(a1)==3 && (connectedAtom(mol, a1, 8, -1, 0, 0)>=0 || connectedAtom(mol, a1, 16, -1, 0, 0)>=0) ) {
						int valence = mol.getConnBondOrder(a1, 0) + mol.getConnBondOrder(a1, 1) + mol.getConnBondOrder(a1, 2);
						if(valence==4 && (mol.getBondOrder(bnd1)==2 || mol.getBondOrder(bnd2)==2)) nElectrons++;
						else if(valence==4) nElectrons += 0;
						else {nAmbiguousC++; nElectrons++;/*if(start<0) start = i; */} 											
					} else { 
						if(mol.getConnAtoms(a1)==3 && start<0) start=i;
						nElectrons++;
					}
				} else if(mol.getAtomicNo(a1)==7) {
					if(mol.getConnAtoms(a1)==3) {
						nElectrons+=2;
					} else if(mol.getConnAtoms(a1)==2) {
						nAmbiguousN++; nElectrons++; 
					} else {
						nElectrons++;
					}
				} else {
					nElectrons+=2;
				}

				if(mol.getBondOrder(bnd2)>1) start = i;
				else if(mol.getImplicitHydrogens(a1)>0
					 && mol.getImplicitHydrogens(a0)>0
					 &&	(mol.getAtomRingBondCount(a1)==2 || mol.getAtomRingBondCount(a0)==2)) {
					if(mol.getConnAtoms(a1)==3 || mol.getConnAtoms(a0)==3)  start = i;
					else if(start<0) start = i;
				}				
			}
			
			int nPyroles = 0;
			int nOxo = 0;
			int diff = nElectrons%4-2; 
			if(diff<0) {
				nPyroles+=Math.min(-diff, Math.max(0, nAmbiguousN)); nElectrons+=nPyroles;
			} else if(diff>0) {
				nOxo+=Math.min(diff, Math.max(0, nAmbiguousC)); nElectrons-=nOxo;
			}

//			if(ringNo==29) {
			
			if(nElectrons%4!=2) {
				if(ring.length==3) continue; //cyclopropane is of course planar but not aromatic
				boolean ok = false;
				if(diff>0) {
					for (int i = 0; i < ring.length; i++) {
						if(mol.getAtomicNo(ring[i])==7) {							
							//toDo[ring[i]]=2;//Protonated N?
							ok=true;
						}
					}					
				} 
				if(!ok) {
					if(!lenient) throw new Exception("Huckel's rule not verified");
					continue;
				}
			}
			aromaticRing[ringNo] = true;
			/*
			if(start<0) start = 0;
			pyrolles[ringNo] = nPyroles;
			oxo[ringNo] = nOxo;
			/*
			for(int i=0; i<ring.length; i++) {
				int a1 = ring[i];
				if(StructureCalculator.getImplicitHydrogens(mol, a1)==0) continue;
				if(mol.getAtomicNo(a1)==7) {
					if(nPyroles>0) {
						toDo[a1] = 2; //This N may not need a double bond
					} else {
						toDo[a1] = 1; //this N needs a double bond
					}
				} else if(mol.getAtomicNo(a1)==6) {
					double doub = 0; for (int j = 0; j < mol.getAllConnAtoms(a1); j++) if(mol.getConnBondOrder(a1, j)>1) doub++;
					if(doub==0) {
						toDo[a1] = 1;
						if(nOxo>0) { 
							for (int j = 0; j < mol.getAllConnAtoms(a1); j++) {
								int a2 = mol.getConnAtom(a1, j);
								if(mol.getAtomicNo(a2)==8 && mol.getAllConnAtoms(a2)==1) {
									toDo[a2] = 2;
								}
							}
						}
					}
				}				
			}
			*/
			
		}

		//Aromatizer
		//Initialize the visited atoms 
		boolean[] visited2 = new boolean[mol.getAllAtoms()];
		Set<Integer> nonAromaticAtom = new HashSet<Integer>();
		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			if(connected(mol, atom, -1, 2)>=0) visited2[atom] = true; //This atom has been processed above
			if(mol.getAtomicNo(atom)==6) {
				boolean ok = false;
				if (atomToRings[atom] != null)
					for(int r: atomToRings[atom])
						if(aromaticRing[r]) ok = true;
				if(!ok) nonAromaticAtom.add(atom);
			}
		}
		
		for (int i = 0; i < aromaticRing.length; i++) {			
			if(aromaticRing[i]) {
				boolean success = aromatize(mol, atomToRings,ringSet, i, aromaticRing, nonAromaticAtom, visited2, 0, ringSet.getRingSize(i)%2, new ArrayList<Integer>(), true);
				if(!success) success = aromatize(mol, atomToRings, ringSet, i, aromaticRing, nonAromaticAtom, visited2, 0, ringSet.getRingSize(i)%2, new ArrayList<Integer>(), false);
				if(!success) {
					System.out.println("Could not aromatize ring "+i);
					aromaticRing[i] = false;
				}
			}
		}
		boolean[] aromaticAtoms = new boolean[mol.getAllAtoms()];
		for (int i = 0; i < ringSet.getSize(); i++) {
			if(aromaticRing[i]) {
				for(int atm: ringSet.getRingAtoms(i)) {
					aromaticAtoms[atm] = true;
				}
			}
		}

		/////////////////////////////////
		//2nd pass: find obvious double bonds on sp2 carbons outside aromatic rings
		for(int atom=0; atom<mol.getAllAtoms(); atom++) {
			if(spOrder[atom]==2 && !aromaticAtoms[atom] && mol.getAtomicNo(atom)==6 && mol.getAllConnAtoms(atom)==3 && mol.getImplicitHydrogens(atom)>0 && connected(mol, atom, -1, 2)>=0) {
				int a1 = mol.getConnAtom(atom, 0);
				int a2 = mol.getConnAtom(atom, 1);
				int a3 = mol.getConnAtom(atom, 2);
				double order1, order2, order3;

				if(mol.getImplicitHydrogens(a1)==0 && connected(mol, a1, -1, 2)>=0) order1 = 1;
				else order1 = getLikelyOrder(mol, atom, a1);

				if(mol.getImplicitHydrogens(a2)==0 && connected(mol, a2, -1, 2)>=0) order2 = 1;
				else order2 = getLikelyOrder(mol, atom, a2);
				
				if(mol.getImplicitHydrogens(a3)==0 && connected(mol, a3, -1, 2)>=0) order3 = 1;
				else order3 = getLikelyOrder(mol, atom, a3);
				
				
				//the highest is the most likely to have a double bond
				int connBond = -1;
				if(order1>order2 && order1>order3 && order1>1 /*&& ((mol.getAtomicNo(i)!=6 && mol.getAtomicNo(a1)!=6) || isPlanar(mol, i, a1))*/) {
					connBond = 0;  					
				} else if(order2>order1 && order2>order3 && order2>1 /*&& ((mol.getAtomicNo(i)!=6 && mol.getAtomicNo(a2)!=6) || isPlanar(mol, i, a2))*/) {
					connBond = 1;  					
				} else if(order3>order1 && order3>order2 && order3>1 /*&& ((mol.getAtomicNo(i)!=6 && mol.getAtomicNo(a3)!=6) || isPlanar(mol, i, a3))*/) {
					connBond = 2;  					
				}  	
				
				if(connBond>=0) {
					mol.setBondOrder(mol.getConnBond(atom, connBond), 2);
				} 
			}
		}		

		//3rd pass, double bonds inside non-aromatic rings
		IntQueue queue = new IntQueue();
		for(int bond=0; bond<mol.getAllBonds(); bond++) {
			if(!visited[bond]) queue.push(bond);
			while(!queue.isEmpty()) {
				int bnd = queue.pop();
				if(visited[bnd]) continue;
				visited[bnd] = true;
				int atm1 = mol.getBondAtom(0, bnd);
				int atm2 = mol.getBondAtom(1, bnd);
				
				//Push the neighbour bonds into the queue
				for(int j=0; j<mol.getAllConnAtoms(atm1); j++) {
					queue.push(mol.getConnBond(atm1, j));
				}
				for(int j=0; j<mol.getAllConnAtoms(atm2); j++) {
					queue.push(mol.getConnBond(atm2, j));
				}
				
				//Compute the free valence and increase the bond order if needed
				double order = getLikelyOrder(mol, atm1, atm2);

				if(order>2 && !aromaticAtoms[atm1] && !aromaticAtoms[atm2]) {
					if(mol.getAtomPi(atm1) != 0) continue; //no adjacent double bonds
					if(mol.getAtomPi(atm2) != 0) continue; //no adjacent double bonds
					
					//Special case CS
					if(mol.getAtomicNo(atm1)==16 && mol.getAllConnAtoms(atm1)<=2) continue;
					if(mol.getAtomicNo(atm2)==16 && mol.getAllConnAtoms(atm2)<=2) continue;
					
					int freeValence1 = getMaxFreeValence(mol, atm1);
					int freeValence2 = getMaxFreeValence(mol, atm2);
						
					//boolean planar = (spOrder[atm1]<3) && (spOrder[atm2]<3); 
					boolean aligned = spOrder[atm1]==1 && spOrder[atm2]==1;					
					if(order>3.0 && freeValence1>1 && freeValence2>1 && aligned) {
						mol.setBondOrder(bnd, 3);
					} else if(freeValence1>0  && freeValence2>0) {
						if((mol.getAtomicNo(atm1)==6 && mol.getAtomicNo(atm2)==6) && !isPlanar(mol, atm1, atm2)) continue;
						if(mol.getAtomicNo(atm1)==6 && spOrder[atm1]>2) continue;
						if(mol.getAtomicNo(atm2)==6 && spOrder[atm2]>2) continue;						 
						mol.setBondOrder(bnd, 2);
					}
				}
			}
		}
	}

	private static ArrayList<Integer>[] getAtomToRings(StereoMolecule mol) {
		ArrayList<Integer>[] atomToRings = new ArrayList[mol.getAllAtoms()];
		RingCollection ringSet = mol.getRingSetSimple();
		for (int r=0; r<ringSet.getSize(); r++) {
			int[] ringAtom = ringSet.getRingAtoms(r);
			for (int atom:ringAtom) {
				if (atomToRings[atom] == null)
					atomToRings[atom] = new ArrayList<>();
				atomToRings[atom].add(r);
			}
		}
		return atomToRings;
	}
	
	public static boolean aromatize(StereoMolecule mol, Set<Integer> aromaticAtoms, Set<Integer> aromaticBonds) {
		ArrayList<Integer>[] atomToRings = getAtomToRings(mol);
		RingCollection ringSet = mol.getRingSetSimple();
		return aromatize(mol,atomToRings,ringSet,aromaticAtoms,aromaticBonds);
	}

	public static boolean aromatize(StereoMolecule mol, ArrayList<Integer>[] atomToRings, RingCollection ringSet, Set<Integer> aromaticAtoms, Set<Integer> aromaticBonds) {
		//RingCollection ringSet = mol.getRingSetSimple();

		//Flag the aromatic rings
		boolean[] aromaticRings = new boolean[ringSet.getSize()];
		for (int i = 0; i < ringSet.getSize(); i++) {
			//Is is an aromatic ring
			boolean isAromatic = true;
			int ringSize = ringSet.getRingSize(i);
			for (int j = -1; j < ringSize; j++) {
				int[] ringAtom = ringSet.getRingAtoms(i);
				int a1 = j==-1? ringAtom[ringSize-1]: ringAtom[j];
				int a2 = j==ringSize-1? ringAtom[0]: ringAtom[j+1];
				
				int b = mol.getBond(a1, a2);
				if(!aromaticBonds.contains(b)) {
					isAromatic = false;
				}
			}
				
			aromaticRings[i] = isAromatic;
		}
		Set<Integer> nonAromaticAtoms = new HashSet<Integer>();
		for (int i=0;i<mol.getAllAtoms();i++) nonAromaticAtoms.add(i);
		nonAromaticAtoms.removeAll(aromaticAtoms);		
		
		//Launch the aromatizer
		boolean ok = true;
		for (int i = 0; i < aromaticRings.length; i++) {
			if(aromaticRings[i]) {
				boolean success = aromatize(mol, atomToRings, ringSet, i, aromaticRings, nonAromaticAtoms, new boolean[mol.getAllAtoms()], 0, ringSet.getRingSize(i)%2, new ArrayList<Integer>(), true);
				if(!success) success = aromatize(mol, atomToRings, ringSet, i, aromaticRings, nonAromaticAtoms, new boolean[mol.getAllAtoms()], 0, ringSet.getRingSize(i)%2, new ArrayList<Integer>(), false);
				if(!success) {
					System.out.println("Could not aromatize ring "+i);
					aromaticRings[i] = false;
					ok = false;
				}
			}
		}
		return ok; 
	}
	
	private static boolean aromatize(StereoMolecule mol, ArrayList<Integer>[] atomToRings, RingCollection ringSet, int index, boolean[] aromatic, Set<Integer> nonAromaticAtoms, boolean[] visited, int seed, int left, List<Integer> bondsMade, boolean easy) {
		//RingCollection ringSet = mol.getRingSetSimple();

		//Ends if the ring has been fully visited
		int[] ring = (int[]) ringSet.getRingAtoms(index);
		boolean allVisited = true;
		int bestSeed = -1;
		for (int i = 0; i < ring.length; i++) {
			if(!visited[ring[(seed + i)%ring.length]]) {				 
				if(bestSeed<0) bestSeed = seed + i; 
				allVisited = false;				
			}
		}
		if(allVisited) {
			return true;
		} else  {
			seed = bestSeed;
		}
		
		
		int a = ring[seed%ring.length];
		int ap = ring[(seed+1)%ring.length];
		if(visited[a]) { //already treated
			
			return aromatize(mol, atomToRings, ringSet, index, aromatic, nonAromaticAtoms, visited, seed+1, left, bondsMade, easy);
			
		} else {
		
			//Try to create a double bond from the atom a to a connected one
			for (int j = -1; j < mol.getAllConnAtoms(a); j++) {
				int a2 = j==-1? ap: mol.getConnAtom(a, j);
				
				if(visited[a2]) continue;
				if(j>=0 && a2==ap) continue;
				
				if(nonAromaticAtoms.contains(a2)) continue;
				if(nonAromaticAtoms.contains(a)) continue;
				
				if(mol.getAtomicNo(a)==8 || mol.getAtomicNo(a)==16) continue;
				if(mol.getAtomicNo(a2)==8 || mol.getAtomicNo(a2)==16) continue;
				if(easy && mol.getFreeValence(a) <= 0) continue;
				if(easy && mol.getFreeValence(a2) <= 0) continue;
				if(connected(mol, a, -1, 2)>=0) continue;
				if(connected(mol, a2, -1, 2)>=0) continue;
				
				visited[a] = visited[a2] = true;
				int b = mol.getBond(a, a2);
				mol.setBondOrder(b, 2);
				
				//Test whole ring 
				List<Integer> trackBondsMade = new ArrayList<Integer>();
				boolean success = aromatize(mol, atomToRings, ringSet, index, aromatic, nonAromaticAtoms, visited, seed+1, left, trackBondsMade, easy);

				//Test connecting rings
				if(success) {
					List<Integer> rings = atomToRings[a2];
					if(rings.size()>1) {
						for (int r : rings) {
							if(r!=index && r>=0 && r<aromatic.length && aromatic[r]) {
								int newSeed;
								for(newSeed=0; ringSet.getRingAtoms(r)[newSeed]!=a2; newSeed++) {}
							
								//System.out.println("try connected ring "+r+" / "+newSeed);
								success = aromatize(mol, atomToRings, ringSet, r, aromatic, nonAromaticAtoms, visited, newSeed, ringSet.getRingSize(r)%2, trackBondsMade, easy);
								//System.out.println("try connected ring "+r +" -> " +success+" "+trackBondsMade.size());
								if(!success) break;
							}
						}
					}
				}
				
				if(success) {
					//It works!!!
					bondsMade.add(b);
					bondsMade.addAll(trackBondsMade);
					return true;
				} else {
					//Backtrack changes
					visited[a] = visited[a2] = false;
					mol.setBondOrder(b, 1);
					for (int b2 : trackBondsMade) {
						//System.out.println("retrack "+mol.getBondAtom(0, b2)+"-"+mol.getBondAtom(1, b2));
						mol.setBondOrder(b2, 1);
						visited[mol.getBondAtom(0, b2)] = visited[mol.getBondAtom(1, b2)] = false;
					}
				} 
			}
			
			//Try to skip this atom
			if(left>0 && (mol.getAtomicNo(a)!=6)) {
				visited[a] = true;
				List<Integer> trackBondsMade = new ArrayList<Integer>();
				boolean success = aromatize(mol, atomToRings, ringSet, index, aromatic, nonAromaticAtoms, visited, seed+1, left-1, trackBondsMade, easy);
				if(success) {
					bondsMade.addAll(trackBondsMade);
					return true;
				} else {
					visited[a] = false;
					for (int b2 : trackBondsMade) {
						mol.setBondOrder(b2, 1);
						visited[mol.getBondAtom(0, b2)] = visited[mol.getBondAtom(1, b2)] = false;
					}
				}
			}
			
			return false;
		}
	}
	
	private static boolean isPlanar(StereoMolecule mol, int a1, int a2) {
		Coordinates ci = mol.getCoordinates(a1);		
		Coordinates u = null, v =null;
		
		for (int i = 0; v==null && i < mol.getAllConnAtoms(a1); i++) {
			if(u==null) u = mol.getCoordinates(mol.getConnAtom(a1, i)).subC(ci);
			else {v = mol.getCoordinates(mol.getConnAtom(a1, i)).subC(ci);} 
		}
		for (int i = 0; v==null && i < mol.getAllConnAtoms(a2); i++) {
			if(u==null) u = mol.getCoordinates(mol.getConnAtom(a2, i)).subC(ci);
			else {v = mol.getCoordinates(mol.getConnAtom(a2, i)).subC(ci);} 
		}
		
		if(u==null) return false;
		
		Coordinates normal = u.cross(v);
		if(normal.distSq()==0) return false; //what to do?
		normal = normal.unitC();

		Coordinates cj = mol.getCoordinates(a2);					
		for(int k=0; k<mol.getAllConnAtoms(a2); k++) {
			Coordinates ck = mol.getCoordinates(mol.getConnAtom(a2, k));
			if(Math.abs(ck.subC(cj).dot(normal))>0.2) return false;
		}					
		for(int k=0; k<mol.getAllConnAtoms(a1); k++) {
			Coordinates ck = mol.getCoordinates(mol.getConnAtom(a1, k));
			if(Math.abs(ck.subC(cj).dot(normal))>0.2) return false;
		}					
		return true;		
	}
	
	private static double getLikelyOrder(StereoMolecule mol, int atm1, int atm2) {
		int k;				
		for(k=0; k<PARAMS.length; k++) if((PARAMS[k][0]==mol.getAtomicNo(atm1) && PARAMS[k][1]==mol.getAtomicNo(atm2)) || (PARAMS[k][1]==mol.getAtomicNo(atm1) && PARAMS[k][0]==mol.getAtomicNo(atm2))) break;
		if(k>=PARAMS.length) return 1;
						
		//Calculate the order
		double r = mol.getCoordinates(atm1).distance(mol.getCoordinates(atm2));
		return Math.exp((PARAMS[k][2] - r) / PARAMS[k][3]);
	}

	/**
	 * Util function for substructure searches
	 * @param mol
	 * @param a
	 * @param atomicNo
	 * @param valence
	 * @param otherAtomicNo
	 * @param otherValence
	 * @return
	 */
	private final static int connectedAtom(StereoMolecule mol, int a, int atomicNo, int valence, int otherAtomicNo, int otherValence) {
		loop: for(int i=0; i<mol.getAllConnAtoms(a); i++) {
			int atm = mol.getConnAtom(a, i);
			if(atomicNo>0 && mol.getAtomicNo(atm)!=atomicNo) continue;
			if(valence>0 && mol.getAllConnAtoms(atm)!=valence) continue;
			if(otherAtomicNo>0 || otherValence>0) {
				for (int j = 0; j < mol.getAllConnAtoms(atm); j++) {
					int otherAtm = mol.getConnAtom(atm, j);
					if(otherAtm==a) continue loop;
					if(otherAtomicNo>0 && mol.getAtomicNo(otherAtm)!=otherAtomicNo) continue loop;
					if(otherValence>0 && mol.getAllConnAtoms(otherAtm)!=otherValence) continue loop;					
				}	
			}
			
			return atm;
		}
		return -1;
	}
	
	private final static int connectedBond(StereoMolecule mol, int a, int atomicNo, int valence) {
		if(a<0) return -1;
		for(int i=0; i<mol.getAllConnAtoms(a); i++) {
			int atm = mol.getConnAtom(a, i);
			if(atomicNo>0 && mol.getAtomicNo(atm)!=atomicNo) continue;
			if(valence>0 && mol.getAllConnAtoms(atm)!=valence) continue;
			return mol.getConnBond(a, i);
		}
		return -1;
	}
	private final static int shortestBond(StereoMolecule mol, int a, int toAtomicNo, boolean privilegeRing) {
		int bestBond = -1;
		double bestDist = Double.MAX_VALUE;
		for(int i=0; i<mol.getAllConnAtoms(a); i++) {
			int atm = mol.getConnAtom(a, i);					
			if(toAtomicNo>0 && mol.getAtomicNo(atm)!=toAtomicNo) continue;
			if(getMaxFreeValence(mol, atm)==0) continue;
			double dist = mol.getCoordinates(a).distance(mol.getCoordinates(atm));
			if (privilegeRing && mol.isRingBond(mol.getConnBond(a, i)))
				dist -= 2;
			if(dist<bestDist) {
				bestDist = dist;
				bestBond = mol.getConnBond(a, i);
			}			
		}
		return bestBond;
	}

	/**
	 * Does consider implicit or explicit hydrogens as part of the free valence
	 * @return
	 */
	private static int getMaxFreeValence(StereoMolecule mol, int a) {
		return mol.getFreeValence(a) + mol.getImplicitHydrogens(a);
	}

	public static int connected(StereoMolecule mol, int a, int atomicNo, int bondOrder) {
		for(int i=0; i<mol.getAllConnAtoms(a); i++) {
			int atm = mol.getConnAtom(a, i);
			if(atomicNo>=0 && mol.getAtomicNo(atm)!=atomicNo) continue;
			if(bondOrder>0 && mol.getConnBondOrder(a, i)!=bondOrder) continue;
			return atm;
		}
		return -1;
	}

	private final static double[][] PARAMS = new double[][] {
		//AtomicNo, AtomicNo, R0, B
		//s = exp((Ro - R)/B)

		{ 6, 6, 1.523, 0.1855}, //1 > 1.523(16)  2 > 1.3944(25)  3 > 1.212(2)  
		{ 6, 7, 1.470, 0.2458}, //1 > 1.47(28)  2 > 1.2996(10)  3 > 1.158(2)  
		{ 6, 8, 1.410, 0.21}, //1 > 1.435(24)  2 > 1.2884(10)  
		{ 6, 16, 1.815, 0.0690}, //1 > 1.815(12)  2 > 1.7672(20)  
		{ 7, 7, 1.381, 0.1919}, //1 > 1.381(49)  2 > 1.248(4)  3 > 1.23(1)  
		{ 7, 8, 1.310, 0.1500}, //1 > 1.31(42)    
		{ 8, 8, 1.428, 0.0000}, //1 > 1.428(36)  
		{ 8, 15, 1.696, 0.22},   
		{ 8, 16, 1.657, 0.24}, //1 > 1.657  2 > 1.48(8)  
		{ 16, 16, 2.024, 0.36}, //1 > 2.024(9)  2 > 1.771(16)  

	};
	
	/*
	public static void main(String[] args) throws Exception {
		PDBFileParser parser = new PDBFileParser();
		parser.setLoadMainStructureOnly(false);
		Molecule3D m = parser.load("D:\\NoBackup\\PDB\\ENTRIES\\e4\\pdb1e47.ent.gz");
		FFViewer.viewMolecule(m);
	}
	*/
}
