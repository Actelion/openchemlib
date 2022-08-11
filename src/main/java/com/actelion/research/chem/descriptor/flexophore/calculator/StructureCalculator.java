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
 * @author Modest v. Korff
 */

package com.actelion.research.chem.descriptor.flexophore.calculator;

import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.VDWRadii;
import com.actelion.research.chem.io.pdb.converter.MoleculeGrid;
import com.actelion.research.util.ArrayUtils;

import java.util.*;


/**
 * Utility function based on the molecule's connections (groups, rings,...)  
 */
public class StructureCalculator {

	/**
	 * Returns a List of all Connex Components of the graph (List of List of Integer)
	 * <pre>
	 * Example:
	 * 
	 *  The following molecule:
	 *
	 * 			  -- 3             
	 * 	   1 -- 2                 5 -- 6              7
	 * 	          -- 4
	 * will return 
	 *  [[1,2,3,4],[5,6],[7]]
	 * </pre>
	 * Complexity: O(nAtoms)
	 * Memory: O(nAtoms)
	 */
	public static List<List<Integer>> getConnexComponents(Molecule3D mol) {
		int[] groups = getAtomToGroups(mol);
		int nGroups = 0;
		for(int i=0; i<groups.length; i++) {
			if(groups[i]>nGroups) nGroups = groups[i]; 
		}
		List<List<Integer>> r = new ArrayList<List<Integer>>();
		for(int i=0; i<nGroups; i++) r.add(new ArrayList<Integer>()); 
		for(int i=0; i<groups.length; i++) {
			r.get(groups[i]-1).add(i);
		}
		return r;
	}
	
	public static int[] getAtomToGroups(Molecule3D mol) {
		return getAtomToGroups(mol, null);
	}
	/**
	 * For each molecule in <code>mol</code> one group is created. 
	 * The group index starts with 1. The size of the array returned is equal to the number of atoms.
	 * 
	 * @param mol
	 * @param seeds
	 * @return a int[] array such as array[atom] = groupNo >=1
	 */
	public static int[] getAtomToGroups(Molecule3D mol, List<Integer> seeds) {
		int[] res = new int[mol.getAtoms()];
		IntQueue q =  new IntQueue();
		int nGroups = 0;
		for(int i=0; i<res.length; i++) {
			if(res[i]>0) continue;
			if(seeds!=null) seeds.add(i);
			q.push(i);
			nGroups++;
			while(!q.isEmpty()) {
				int a = q.pop();
				res[a] = nGroups;
				for(int j=0; j<mol.getAllConnAtoms(a); j++) {
					int a2 = mol.getConnAtom(a, j);
					if(a2<mol.getAtoms()){
						if(res[a2]==0) {
							q.push(a2);
							res[a2]=-1; //safeguard
						}
					}
				}
			}
			q.clear();
		}
		return res;
	}
	
	private static class Item {
		public Item(int a) {this.a = a;path.add(a);}
		public int a;
		public ArrayList<Integer> path = new ArrayList<Integer> ();
	}
	
	/**
	 * Return the longest molecular chain starting at atm
	 * @param mol
	 * @param atm
	 * @return a List of Integer
	 */
	@SuppressWarnings("unchecked")
	public static List<Integer> getLongestChain(Molecule3D mol, int atm) {
		LinkedList<Item> q = new LinkedList<Item>();
		q.add(new Item(atm));
		boolean[] seen = new boolean[mol.getAllBonds()*2];
		try {
		while(true) {
			Item item = q.removeFirst();
			int a = item.a;
			for(int i=0; i<mol.getAllConnAtoms(a); i++) {
				int b = mol.getConnBond(a, i);
				int a2 = mol.getConnAtom(a, i);
				int bi = b*2 + (a2>a?1:0);
				if(seen[bi]) continue;
				seen[bi] = true;
				if(mol.getAtomicNo(a2)<=1) continue;
				Item ni = new Item(a2);
				ni.a = a2;
				ni.path = (ArrayList<Integer>) item.path.clone();				
				ni.path.add(a2);
				q.add(ni);
			}			
			if(q.isEmpty()) return item.path;
		}
		}catch(Exception e ){
			e.printStackTrace();
			return null; 
		}		
	}
		
	public static boolean[] getBackbones(Molecule3D molecule) {
		return getBackbones(molecule, 70);
	}
	
	/**
	 * Gets the backbone 
	 * @param molecule
	 * @return an array of boolean 
	 */
	public static boolean[] getBackbones(Molecule3D molecule, int minChain) {
		boolean[] inBackbone = new boolean[molecule.getAllAtoms()];
		
		//Find the fragments (connex components of the graph)
		List<List<Integer>> fragments = getConnexComponents(molecule);
		
		
		//For each fragment, find the longest chain
		//List res[] = new List[fragments.size()];
		for(int i=0; i<fragments.size(); i++) {			
			boolean markAll = false;
			
			//Find an extremity
			List<Integer> fragment = fragments.get(i);
			List<Integer> l = null;
			int root = fragment.get(0);			
			
			if(markAll) {
				for(int j=0; j<fragment.size(); j++) {
				   inBackbone[((Integer)fragment.get(j)).intValue()] = true;
				}						
			} else {
				l = StructureCalculator.getLongestChain(molecule, root);
				root = ((Integer)l.get(l.size()-1)).intValue();
				l = StructureCalculator.getLongestChain(molecule, root);				
				if(l.size()<minChain && molecule.getAllAtoms()>80) markAll=true;
				for(int j=0; j<l.size(); j++) {
				   inBackbone[((Integer)l.get(j)).intValue()] = true;
				}		
			}
		}
		
		return inBackbone;
	}
	

	/**
	 * Computes a matrix of distances between all the atoms in the graph.
	 * Complexity: O(m*n*maxBonds) m = number of bonds, n = number of atoms
	 * @param mol
	 * @param maxBonds
	 * @return an array A[i][j] = nBonds if i and j are connected by less than maxBonds
	 * 						      or -1  otherwise  
	 */
	public static int[][] getNumberOfBondsBetweenAtoms(Molecule3D mol, int maxBonds, int[][] dist) {
		//Initialization of the data
		if(dist==null)  dist = new int[mol.getAllAtoms()][mol.getAllAtoms()];
		int N = dist.length;  
		for(int i=0; i<N; i++) {
			dist[i][i] = 0; 
			for(int j=i+1; j<N; j++) {
				dist[i][j] = dist[j][i] = -1;
			}
		}
		
		//Algo: Dynamic Programming 
		for(int j=0; j<maxBonds; j++) { //Maximum of nBonds bonds
			for(int i=0; i<mol.getAllBonds(); i++) {
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
	 * Computes a matrix of distances between all the bonds in the graph.
	 * Complexity: O(m*n*maxBonds) m = number of bonds, n = number of atoms
	 * @param mol
	 * @param maxAtoms
	 * @return an array B[i][j] = nAtoms if i and j are connected by less than maxAtoms
	 * 						      or -1  otherwise  
	 */
	public static int[][] getNumberOfAtomsBetweenBonds(Molecule3D mol, int maxAtoms, int[][] dist) {
		//Initialization of the data
		if(dist==null)  dist = new int[mol.getAllBonds()][mol.getAllBonds()];
		int N = dist.length;  
		for(int i=0; i<N; i++) {
			dist[i][i] = 0; 
			for(int j=i+1; j<N; j++) {
				dist[i][j] = dist[j][i] = -1;
			}
		}
		
		//Algo: Dynamic Programming 
		// Dynamic Programming: dist(b0,b2) = min(dist(b0, b2), dist(b0,b1)+1) [if dist(b0,b1)>0] 
		for(int iter=0; iter<maxAtoms; iter++) { 
			for(int a=0; a<mol.getAllAtoms(); a++) {
				for (int k1 = 0; k1 < mol.getConnAtoms(a); k1++) {
					int b1 = mol.getConnBond(a, k1);
					for (int k2 = k1+1; k2 < mol.getConnAtoms(a); k2++) {
						int b2 = mol.getConnBond(a, k2);
						if(b1>=N || b2>=N) continue;

						for(int b0=0; b0<N; b0++) {
							// Dynamic Programming: dist(a0,a2) = min(dist(a0, a2), dist(a0,a1)+1) [if dist(a0,a1)>0] 
							if(dist[b0][b1]>=0 && (dist[b0][b2]==-1 || dist[b0][b1]+1<dist[b0][b2]) && dist[b0][b1]<maxAtoms) {
								dist[b2][b0] = dist[b0][b2] = dist[b0][b1] + 1;
							}
							if(dist[b0][b2]>=0 && (dist[b0][b1]==-1 || dist[b0][b2]+1<dist[b0][b1]) && dist[b0][b2]<maxAtoms) {
								dist[b1][b0] = dist[b0][b1] = dist[b0][b2] + 1;
							}
						}
					}
				}
			}
		}
		return dist;		
	}
	
	public static void main(String[] args) throws Exception {
		StereoMolecule sm = new StereoMolecule();
		new SmilesParser().parse(sm, "C1CCC1OCCCC");
		Molecule3D m = new Molecule3D(sm);
		int[][] B = getNumberOfAtomsBetweenBonds(m, m.getAllAtoms(), null);
		
		System.out.print("\t");
		for (int j = 0; j < B.length; j++) {
			System.out.print("_"+j+"_\t");
		}
		System.out.println();
		for (int i = 0; i < B.length; i++) {
			System.out.print("_"+i+"_\t");
			for (int j = 0; j < B[i].length; j++) {
				System.out.print((B[i][j]>=0?B[i][j]:"")+"\t");
			}
			System.out.println();
		}
	}
	
	/**
	 * @param atomicNo
	 * @return
	 */
	public static int getMaxValence(int atomicNo) {
		if (atomicNo >= 171 && atomicNo <= 190) return 2;
		switch (atomicNo) {
			case  1: return 1; 	// H
			case  5: return 3; 	// B
			case  6: return 4; 	// C
			case  7: return 3; 	// N
			case  8: return 2; 	// O
			case  9: return 1; 	// F
			case 13: return 3; 	// Al
			case 14: return 4; 	// Si
			case 15: return 3;  // P
			case 16: return 2; 	// S
			case 17: return 1; 	// Cl
			case 33: return 3; 	// As
			case 34: return 2; 	// Se
			case 35: return 1; 	// Br
			case 52: return 2; 	// Te
			case 53: return 1; 	// I
			default: return 0;
		}
	}
	
	/**
	 * Return the number of implicit hydrogens, ie. the number of hydrogens to
	 * add to fullfill the valence requirements (<0 if there are too many)
	 * @param mol
	 * @param atm
	 * @return
	 */
	public static int getImplicitHydrogens(Molecule3D mol, int atm) {
		int atomicNo = mol.getAtomicNo(atm);
		int hydrogensToAdd = getMaxValence(atomicNo);		
		if (atomicNo == 6) {
			hydrogensToAdd -= Math.abs(mol.getAtomCharge(atm));
		} else if (Molecule.isAtomicNoElectronegative(atomicNo)) {
			hydrogensToAdd += mol.getAtomCharge(atm);
		} else {
			hydrogensToAdd -= mol.getAtomCharge(atm);
		}
			
		for(int i=0; i<mol.getAllConnAtoms(atm); i++) {
			if(mol.getAtomicNo(mol.getConnAtom(atm, i))>=1) hydrogensToAdd -= mol.getConnBondOrder(atm, i);
		}

		if(atomicNo==15) return 0;
		else if(atomicNo==16 && hydrogensToAdd<0) hydrogensToAdd+=Math.min(4, -hydrogensToAdd); //S has 2, 3 or 4 conn atoms
		else if(atomicNo==8 && mol.getAllConnAtoms(atm)==1 && mol.getAtomicNo(mol.getConnAtom(atm,0))==15) return 0;
		
		return hydrogensToAdd;
	}
	
	public static int getMaxFreeValence(Molecule3D mol, int atm) {
		switch(mol.getAtomicNo(atm)) {
			case 15:
				return 5 - getValence(mol, atm);
			case 16:
				return 6 - getValence(mol, atm);
			default:
				return getImplicitHydrogens(mol, atm) + getExplicitHydrogens(mol, atm);			
		}
		
	}
	
	/**
	 * Return all connex components with more than 5 atoms
	 * @param mol
	 * @return
	 */
	public static List<Molecule3D> extractFragments(Molecule3D mol) {
//		List<List<Integer>> comps = getConnexComponents(mol);
//		List<Molecule3D> res = new ArrayList<Molecule3D>();
//		int count = 0;
//		for (List<Integer> comp : comps) {
//
//			Molecule3D m = new Molecule3D(mol.getName() + " / F" + (++count));
//
//			extractFragment(mol, m, ArrayUtils.toIntArray(comp));
//			if(m.getAllAtoms()>5) res.add(m);
//		}
//
//		return res;
		return null;
	}
	
	/**
	 * Return all ligands with more than 5 atoms
	 * @param mol
	 * @return
	 */
	public static List<Molecule3D> extractLigands(Molecule3D mol) {
		List<List<Integer>> comps = getConnexComponents(mol);
		List<Molecule3D> res = new ArrayList<Molecule3D>();
		for (List<Integer> comp : comps) {
			
			Molecule3D m = new Molecule3D();
			
			extractFragment(mol, m, ArrayUtils.toIntArray(comp));
			if(m.getAllAtoms()>5 && m.isAtomFlag(0, Molecule3D.LIGAND)) res.add(m);
		}
		
		return res;

	}
	
	public static <T extends Molecule3D> T extractFragment(Molecule3D mol, T res, List<Integer> atoms) {
		return extractFragment(mol, res, ArrayUtils.toIntArray(atoms));
		
	}
	/**
	 * Extract the Ligand from mol and copy it into res 
	 * @param mol
	 * @param res
	 * @return
	 */
	public static <T extends Molecule3D> T extractFragment(Molecule3D mol, T res, int[] atoms) {
		res.clear();
		int[] oldToNew = new int[mol.getAllAtoms()];
		Arrays.fill(oldToNew, -1);
		for(int i=0; i<atoms.length; i++) {
			oldToNew[atoms[i]] = res.addAtom(mol, atoms[i]);
			
			for(int j=0; j<mol.getAllConnAtoms(atoms[i]); j++) {
				int a = mol.getConnAtom(atoms[i], j);
				if(oldToNew[a]<0) continue;
				res.addBond(oldToNew[atoms[i]], oldToNew[a], mol.getConnBondOrder(atoms[i], j));
			}
		}
		return res;
	}	
	

	private final static class Node {
		int atm;
		int bond; //bond going from parent.atm to atm
		Node parent;
		public Node(int value, Node parent, int bondToParent) {
			this.atm = value;
			this.parent = parent;
			this.bond = bondToParent;
		}
		@Override
		public String toString() {
			Node n = this.parent;
			String r = "{"+atm;
			while(n!=null) {
				r += ","+n.atm; n = n.parent;					
			}
			return r+"}";
		}
		
	}

	public static List<Integer>[] getRings(Molecule3D mol, List<int[]> allRings) {
		if(mol.getAllAtoms()<100) return getRingsAccurate(mol, allRings);
		return getRingsFast(mol, allRings);
	}

	/**
	 * Find the smallest covering set of rings for the Molecule. This function considers the smallest subset
	 * of covering rings of any size.
	 * 
	 * The complexity of this function is O(natoms * avgRingSize^nImbricatedRings)
	 * The memory required is O(nAtoms * nBonds)!!!!
	 * 
	 * @param mol
	 * @param allRings - a List ringNo -> atom No in the ring as int[] (output)
	 * @return a List[] atom -> List of rings the atom belongs to
	 */
	@SuppressWarnings("unchecked")
	public static List<Integer>[] getRingsAccurate(Molecule3D mol, List<int[]> allRings) {
		//long s = System.currentTimeMillis();
		//a List[] of Node representing the leaves that can be reached from i 
		List<Node>[] connectables = new ArrayList[mol.getAllAtoms()];
		int maxRingSize = 0;
		for(int atm=0; atm<mol.getAllAtoms(); atm++) {
			if(mol.getAtomicNo(atm)==1) continue;
			maxRingSize++;
			Node node = new Node(atm, null, 0);
			connectables[atm] = new ArrayList<Node>();
			connectables[atm].add(node);
		}
		if(maxRingSize>50) maxRingSize = 50;
		
		int N = mol.getAllAtoms();
		int M = mol.getAllBonds();

		boolean coveredBonds[] = new boolean[M];
		boolean visitedBonds[][] = new boolean[N][M*2];
		int iteration = 0;
		
		while(iteration<maxRingSize) {
			iteration++;
			//Create a connection graph simultaneously for each atom until all bonds
			//have been checked. As a consequence, the smallest rings will be found first
			for(int atm=0; atm<N; atm++) {

				List<Node> lastListOfConnectables = connectables[atm];
				if(lastListOfConnectables==null) continue;
				List<Node> newListOfConnectables = new ArrayList<Node>();
				Iterator<Node> iter = lastListOfConnectables.iterator();
				loopPaths: while(iter.hasNext()) {
					Node node = iter.next(); 
					int atm2 = node.atm;
					
					for(int j=0; j<mol.getAllConnAtoms(atm2); j++) {
						int atm3 = mol.getConnAtom(atm2, j);
						int bond = mol.getConnBond(atm2, j);
						int orientedBond = bond + (atm3<atm2?mol.getAllBonds():0);

						if(mol.getAtomicNo(atm3)==1) continue;
						if(iteration==1 && atm3<atm) continue;
						if(node.parent!=null && node.parent.atm==atm3) continue;
						if(visitedBonds[atm][orientedBond]) continue;
						visitedBonds[atm][orientedBond] = true;						
						
						//Check for rings going from atm to atm3
						if(atm3==atm) {
							//We found a ring
							Node n = node;
							List<Integer> ring = new ArrayList<Integer>();
							boolean isCovering = true;
							while(n!=null) {
								ring.add(0, n.atm);
								if(!coveredBonds[n.bond]) {
									isCovering=false;
									coveredBonds[n.bond] = true;
								}
								n = n.parent;
							}
							if(isCovering) continue;
							if(ring.size()<3) {
								System.err.println("ring of size "+ring.size()+"???");
								continue;
							}

							int[] tmp = ArrayUtils.toIntArray(ring);
							if(!contains(allRings, tmp)) {								
								allRings.add(tmp);
								//System.out.println("ring "+tmp.length);								
								//the maximum ring size is now reduced
								maxRingSize = maxRingSize - (tmp.length-1)/2; 
							}
						} else {
							Node n = node.parent;
							while(n!=null) {
								if(n.atm==atm3) continue loopPaths;
								n = n.parent;
							}
							
							//update the new paths	
							n = new Node(atm3, node, bond);					
							newListOfConnectables.add(n);
						}
					}
				}
				connectables[atm] = newListOfConnectables;
			}		 
		}
		//Assign the rings to each atom
		List<Integer>[] atomToRing = new ArrayList[mol.getAllAtoms()];
		for(int j=0; j<atomToRing.length; j++) atomToRing[j] = new ArrayList<Integer>();
		for(int i=0; i<allRings.size(); i++) {
			int[] atoms = allRings.get(i);
			for(int j=0; j<atoms.length; j++) {
				atomToRing[atoms[j]].add(i);
			}			
		}
		//System.out.println(System.currentTimeMillis()-s+"ms");
		
		return atomToRing;
	}

	/**
	 * Find the list of rings for the Molecule. This function considers all rings
	 * of size  <=7
	 * 
	 * Complexity O(nAtoms) for the common structures O(nAtoms^nRings) at worse 
	 * Memory O(nBonds)
	 * 
	 * @param mol
	 * @param allRings - a List ringNo -> atom No in the ring as int[] 
	 * @return a List[] atom -> List of rings the atom belongs to
	 */
	@SuppressWarnings("unchecked")
	public static List<Integer>[] getRingsFast(Molecule3D mol, List<int[]> allRings) {
		//boolean[] visitedAtoms = new boolean[mol.getAllAtoms()];
		boolean[] visitedBonds = new boolean[mol.getAllBonds()];
		for (int atm = 0; atm < mol.getAllAtoms(); atm++) {
			if(mol.getAtomicNo(atm)<=1) continue;
			List<Integer> previous = new ArrayList<Integer>();
			previous.add(atm);
			
			doRings(mol, allRings, previous, visitedBonds, 7);
		}
		//Remove 6 and 7-rings that are contained in 2 smaller rings
		for (int i = 0; i < allRings.size(); i++) {
			if(allRings.get(i).length<6) continue;
			int found = 0;
			for (int j = 0; j < allRings.size(); j++) {
				if(i==j || allRings.get(j).length>=6) continue;
				
				int atomNotFound=0;
				nextAtom: for(int a: allRings.get(j)) {
					for (int a2: allRings.get(i)) {
						if(a==a2) continue nextAtom;
					}
					atomNotFound++;
				}				
				if(atomNotFound<=1) found++;
			}
			if(found>=2) {
				allRings.remove(i--);
			}
		}
		
		//Sort the rings according to their size
		Collections.sort(allRings, new Comparator<int[]>() {
			@Override
			public int compare(int[] o1, int[] o2) {
				return o1.length - o2.length;
			}
		});
		
		//Assign the rings to each atom
		List<Integer>[] atomToRing = new ArrayList[mol.getAllAtoms()];
		for(int j=0; j<atomToRing.length; j++) atomToRing[j] = new ArrayList<Integer>();
		for(int r=0; r<allRings.size(); r++) {
			for(int a: allRings.get(r)) {
				atomToRing[a].add(r);
			}			
		}		
		return atomToRing;
	}
	
	private static void doRings(Molecule3D mol, List<int[]> allRings, List<Integer> previous, boolean[] visitedBonds, int depth) {
		if(depth<=0) return;
		int firstAtom = previous.get(0);
		int lastAtom = previous.get(previous.size()-1);
		for(int i=0; i<mol.getAllConnAtoms(lastAtom); i++) {
			int b = mol.getConnBond(lastAtom, i); 
			int a = mol.getConnAtom(lastAtom, i);
			if(visitedBonds[b]) continue;
			if(mol.getAtomicNo(a)<=1) continue;
			if(a<firstAtom) continue;
			
			if(a==firstAtom && previous.size()>2) {
				//We found a ring
				if(previous.get(1)>previous.get(previous.size()-1)) continue;
				int[] ring = ArrayUtils.toIntArray(previous);
				allRings.add(ring);
			} else {
				visitedBonds[b]=true;
				previous.add(a);

				doRings(mol, allRings, previous, visitedBonds, depth-1);
				
				previous.remove(previous.size()-1);
				visitedBonds[b]=false;			
			}
		}
	}
	
	/**
	 * Checks if rings (List of int[]) contains ring (int[])
	 * independantly of the order  
	 */	
	private static boolean contains(List<int[]> rings, int[] ring) {
		Iterator<int[]> iter = rings.iterator();
		next: while(iter.hasNext()) {
			int[] r = iter.next();
			if(r.length!=ring.length) continue;
			//compare r and ring
			for(int i=0; i<r.length; i++) {
				boolean found = false;
				for(int j=0; j<ring.length; j++) {
					if(ring[j]==r[i]) found = true;
				}
				if(!found) continue next;
			}
			return true;
		}
		return false;
	}
	
	/**
	 * Return a List of int[] representing all the atom-atom pairs having
	 * an intermolecular interactions (distance close to sum of VDW)  
	 * @param mol
	 * @return
	 */
	public static List<int[]> getInterMolecularInteractions(Molecule3D mol) {
//		//long s = System.currentTimeMillis();
//		int[] atomToGroups = getAtomToGroups(mol);
//		List<int[]> res = new ArrayList<int[]>();
//		MoleculeGrid grid = new MoleculeGrid(mol);
//		for(int a1=0; a1<mol.getAllAtoms(); a1++) {
//			if(mol.getAtomicNo(a1)<=1) continue;
//			Set<Integer> neighbours = grid.getNeighbours(new Coordinates(mol.getAtomX(a1),mol.getAtomY(a1),mol.getAtomZ(a1)), 4);
//			Iterator<Integer> iter = neighbours.iterator();
//			while(iter.hasNext()) {
//				int a2 = ((Integer)iter.next()).intValue();
//				if(a2<a1) continue;
//				if(mol.getAtomicNo(a2)<=1) continue;
//				if(mol.getAllConnAtoms(a1)>0 && mol.getAllConnAtoms(a2)>0 && atomToGroups[a1]==atomToGroups[a2]) continue;
//				double r2 = mol.getCoordinates(a1).distanceSquared(mol.getCoordinates(a2));
//				double vdw = VDWRadii.VDW_RADIUS[mol.getAtomicNo(a1)] + VDWRadii.VDW_RADIUS[mol.getAtomicNo(a2)];
//
//				if( r2>(vdw-.5)*(vdw-.5) && r2 < vdw*vdw) {
//					res.add(new int[]{a1, a2});
//				}
//			}
//		}
//		//System.out.println("inter: "+(System.currentTimeMillis()-s)+"ms");
//		return res;

		return null;
	}
	
	/**
	 * Finds on all atoms going on a path from a1 to a2.
	 * If a1 and a2 are in a small ring (<6), all atoms in the rings will be considered
	 * Complexity: O(m) m=allBonds 
	 * @param mol
	 * @param a1
	 * @param a2
	 * @return
	 */
	public static int[] getAtomsOnPath(Molecule3D mol, int a1, int a2) {
		
		if(a1>=mol.getAllAtoms() || a2>=mol.getAllAtoms()) throw new IllegalArgumentException("Invalid atom number");
		
		boolean[] visited = new boolean[mol.getAllAtoms()];
		SortedSet<Integer> set = new TreeSet<Integer>();
		
		//Push a1
		List<Node> leaves = new ArrayList<Node>();
		Node n = new Node(a1, null, -1);
		leaves.add(n);
				
		int timeout = -1;
		for(int loop=0; loop<mol.getAllBonds() && timeout!=0 && leaves.size()>0; loop++, timeout--) {			
			n = leaves.remove(0);
			int a = n.atm;
			
			if(a==a2) {
				for(Node tmp = n; tmp.parent!=null; tmp = tmp.parent) {
					set.add(tmp.atm);
				}
				timeout = 5; 
				continue;
			}			
			
			if(visited[a]) continue;
			visited[a] = true;
			
			for(int j=0; j<mol.getAllConnAtoms(a); j++) {
				int atm = mol.getConnAtom(a, j);
				if(visited[atm]) continue;				
				Node leaf = new Node(atm, n, mol.getConnBond(a, j));
				leaves.add(leaf);
			}
		}

		if(timeout>=0) {
			set.add(a1);
			set.add(a2);
		}
		int[] res = new int[set.size()];
		int count = 0;
		for(int a: set) {
			res[count++] = a;
		}
				
		return res;
	}
	
	
	public static int getValence(Molecule3D mol, int atm) {
		int res = 0;
		for(int i=0; i<mol.getAllConnAtoms(atm); i++) {
			if(mol.getAtomicNo(mol.getConnAtom(atm, i))>1) res += mol.getConnBondOrder(atm, i);
		}
		return res;
	}


	public static final int getExplicitHydrogens(Molecule3D mol, int atm) {
		int n = 0;
		for(int i=0; i<mol.getAllConnAtoms(atm); i++) {
			if(mol.getAtomicNo(mol.getConnAtom(atm, i))==1) n++;
		}
		return n;		
	}

	public static int getStructureCenter(Molecule3D lig, int[] rotatables, int[][] dists) {
		return getStructureCenter(lig, 0, rotatables, dists);
	}
	/**
	 * The Structure Center Atom is defined as the atom with the biggest no of rotatables bonds and the closest to the center
	 * @param lig
	 * @return
	 */	
	public static int getStructureCenter(Molecule3D lig, int a, int[] rotatables, int[][] dists) {
		List<Integer> backbone = StructureCalculator.getBackbone(lig, a);
		backbone = StructureCalculator.getBackbone(lig, backbone.get(0));
		int bbCenter = backbone.get(backbone.size()/2);
		return bbCenter;		
	}
	
	public static List<Integer> getBackbone(Molecule3D lig, int a) {
		List<Integer> backbone = StructureCalculator.getLongestChain(lig, a);
		backbone = StructureCalculator.getLongestChain(lig, backbone.get(backbone.size()-1));
		return backbone;
	}

	public static int dfs(Molecule3D lig, int start, Set<Integer> seen) {
		return dfs(lig, start, seen, 99, false);
	
	}
	
	public static int dfs(Molecule3D lig, int start, Set<Integer> seen, int depth, boolean takeRingsAsWhole) {
//		if(seen==null) seen = new HashSet<Integer>();
//		List<Integer> stack = new ArrayList<Integer>();
//		List<Integer> newStack = new ArrayList<Integer>();
//		List<Integer> inRings = new ArrayList<Integer>();
//		stack.add(start);
//		int size = 0;
//		while(depth>0 && !stack.isEmpty()) {
//			while(!stack.isEmpty() && depth>0) {
//				int a = stack.remove(0);
//				if(seen.contains(a)) continue;
//				size++;
//				seen.add(a);
//				for (int i = 0; i < lig.getAllConnAtoms(a); i++) {
//					int a2 = lig.getConnAtom(a, i);
//					if(seen.contains(a2)) continue;
//					newStack.add(a2);
//					List<Integer> ring = lig.getAtomToRings()[a2];
//					inRings.addAll(ring);
//				}
//			}
//			depth--;
//			stack = newStack;
//			newStack = new ArrayList<Integer>();
//		}
//
//		if(takeRingsAsWhole) {
//			for (Integer a : inRings) {
//				if(!seen.contains(a)) {
//					size++;
//					seen.add(a);
//				}
//			}
//		}
//
//		return size;

		return -1;
	}
	


	/**
	 * Extract the Ligand from mol and copy it into res 
	 * @param mol
	 * @return
	 */
	public static Molecule3D extractLigand(Molecule3D mol) {
		if(mol==null) throw new IllegalArgumentException("The mol is null");
		if(mol.getNMovables()==mol.getAllAtoms()) return new Molecule3D(mol);
		Molecule3D res = new Molecule3D(mol);
		int[] oldToNew = new int[mol.getAllAtoms()];
		Arrays.fill(oldToNew, -1);
		
		for(int i=0; i<mol.getAllAtoms(); i++) {
			if(mol.isAtomFlag(i, Molecule3D.LIGAND) && oldToNew[i]<0) {
				oldToNew[i] = res.addAtom(mol, i);
				
				for(int j=0; j<mol.getAllConnAtoms(i); j++) {
					int a = mol.getConnAtom(i, j);
					if(oldToNew[a]<0) continue;
					res.addBond(oldToNew[i], oldToNew[a], mol.getConnBondOrder(i, j));
				}
			}
		}
//		res.setName((mol.getName()==null?"": mol.getName()) + (res.getAllAtoms()<mol.getAllAtoms()?" Ligand":""));
		return res;
	}

	public static void insertLigand(Molecule3D mol, Molecule3D lig, Coordinates center) {
		//Add the new ligand
  		int index = mol.fusion(lig);
  		for(int i=index; i<mol.getAllAtoms(); i++) {
	  		mol.setAtomFlag(i, Molecule3D.LIGAND | Molecule3D.PREOPTIMIZED, true);
  		}

  		//Recenter the ligand if needed
  		if(center!=null && GeometryCalculator.getCenterGravity(lig).distance(center)>5 ) {
	  		translateLigand(mol, center.subC(GeometryCalculator.getCenterGravity(lig)));
	  }		
	}		

	public static void replaceLigand(Molecule3D mol, Molecule3D lig) {
		replaceLigand(mol, lig, null);
	}
	
	/**
	 * Replaces the preoptimized ligand in mol
	 * @param mol
	 * @param lig
	 */
	public static void replaceLigand(Molecule3D mol, Molecule3D lig, Coordinates center) {
		
		if(center==null) center = StructureCalculator.getLigandCenter(mol); 
		//Delete the previous ligand
		removeLigand(mol);
		
		if(lig!=null) insertLigand(mol, lig, center);
	}

	
	public static boolean deleteHydrogens(Molecule3D mol) {
		boolean changed = false;
		for (int i = mol.getAllAtoms()-1; i >=0; i--) if(mol.getAtomicNo(i)<=1) {changed = true; mol.deleteAtom(i);}
		return changed;
	}					

	/**
	 * Add the missing hydrogens for the ligand
	 * @param mol
	 * @return
	 */
	public static boolean addHydrogens(Molecule3D mol) {
		return addHydrogens(mol, false);
	}
	/**
	 * Adds the missing hydrogens (no optimization)
	 * @param mol
	 * @return
	 */
	public static boolean addHydrogens(Molecule3D mol, boolean alsoRigidAtoms) {
		boolean changed = false;
		int N = mol.getAllAtoms();
		List<Integer> atomsToBeDeleted = new ArrayList<Integer>();
		
		Coordinates[] bounds = null;
		if(alsoRigidAtoms) {
			bounds = StructureCalculator.getLigandBounds(mol);
			bounds[0].sub(new Coordinates(11,11,11));
			bounds[1].add(new Coordinates(11,11,11));
		}
		 
		
		for(int i=0; i<N; i++) {
			if(!mol.isAtomFlag(i, Molecule3D.RIGID)) {
				//Ok, add H
			} else if(!alsoRigidAtoms) {
				continue; //Don't add H
			} else {
				//Add H if close to ligand
				if(!mol.getCoordinates(i).insideBounds(bounds)) continue;
				
			}
			int n = StructureCalculator.getImplicitHydrogens(mol, i);
			
			if(n<0) {
				for (int j = 0; j < mol.getAllConnAtoms(i) && n<0; j++) {
					int a = mol.getConnAtom(i, j);
					if(mol.getAtomicNo(a)==1) {
						atomsToBeDeleted.add(a);
						changed = true;
						n++;
					}
				}				
			} else if(n>0) {
				for(int j=0; j<n; j++) {
					int a = mol.addAtom(1);
					mol.setAtomFlags(a, mol.getAtomFlags(i) & ~Molecule3D.PREOPTIMIZED);
					mol.addBond(i, a, 1);
					mol.setCoordinates(a, mol.getCoordinates(i));
				}
				changed=true;			
			}
		}
		
		mol.deleteAtoms(atomsToBeDeleted);
		
		return changed;		
	}
	
	public static boolean addHydrogensAroundLigand(Molecule3D mol, double radius) {
//		boolean changed = false;
//		int N = mol.getAllAtoms();
//		MoleculeGrid grid = new MoleculeGrid(mol);
//		List<Integer> atomsToBeDeleted = new ArrayList<Integer>();
//		for(int i=0; i<N; i++) {
//			if(mol.isAtomFlag(i, Molecule3D.RIGID)) {
//				int neighbours = grid.getNeighbours(mol.getCoordinates(i), radius, true, Molecule3D.LIGAND).size();
//				if(neighbours==0) continue;
//			}
//			int n = StructureCalculator.getImplicitHydrogens(mol, i);
//
//			if(n<0) {
//				for (int j = 0; j < mol.getAllConnAtoms(i) && n<0; j++) {
//					int a = mol.getConnAtom(i, j);
//					if(mol.getAtomicNo(a)==1) {
//						atomsToBeDeleted.add(a);
//						changed = true;
//						n++;
//					}
//				}
//
//			} else {
//				for(int j=0; j<n; j++) {
//					int a = mol.addAtom(1);
//					mol.setAtomFlags(a, mol.getAtomFlags(i) & ~Molecule3D.PREOPTIMIZED);
//					mol.addBond(i, a, 1);
//					mol.setCoordinates(a, mol.getCoordinates(i));
//					changed=true;
//				}
//			}
//		}
//
//		mol.deleteAtoms(atomsToBeDeleted);
//
//		return changed;

		return false;
	}
	
	
	
	public static void copyAtoms(Molecule3D mol, Molecule3D res, List<Integer> atomsToBeAdded) {
		
		Map <Integer, Integer> molToRes = new HashMap<Integer, Integer>();		
		//Copy the atoms
		for(int i=0; i<atomsToBeAdded.size(); i++) {
			int a = atomsToBeAdded.get(i);
			if(molToRes.containsKey(a)) continue;
			int n = res.addAtom(mol, a);
			molToRes.put(a,n);
		}
		
		//Add the bonds
		for(int i=0; i<mol.getAllBonds(); i++) {
			int mol1 = mol.getBondAtom(0, i);
			int mol2 = mol.getBondAtom(1, i);
			
			Integer a1 = molToRes.get(mol1);
			Integer a2 = molToRes.get(mol2);
			if(a1!=null && a2!=null) {
				 res.addBond(a1.intValue(), a2.intValue(), mol.getBondOrder(i));				 
			}			
		}
	}	
	

	
	public static void removeWater(Molecule3D mol) {
		removeWater(mol, true);
	}

	public static void removeWater(Molecule3D mol, boolean removeAlsoImportant) {
		if(mol==null) return;	
		for(int i = 0; i<mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i)==0) {
				mol.deleteAtom(i); i--;			
			}			
		}	
		for(int i = 0; i<mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i)==8 && (removeAlsoImportant || !mol.isAtomFlag(i, Molecule3D.IMPORTANT))) {
				if(mol.getAllConnAtoms(i)==0) {
					mol.deleteAtom(i); i--;
				} else if(mol.getAllConnAtoms(i)==2 && mol.getAtomicNo(mol.getConnAtom(i,0))==1 && mol.getAtomicNo(mol.getConnAtom(i,1))==1) {
					int[] atms = new int[]{mol.getAllConnAtoms(i), mol.getConnAtom(i,0), mol.getConnAtom(i,1)};
					Arrays.sort(atms);
					for(int j=2; j>=0; j--) mol.deleteAtom(atms[j]);
					i-=3;				
				}				
			}
		}
	}
		
	
	public static void removeWaterSalts(Molecule3D mol) {
		List<List<Integer>> groups = StructureCalculator.getConnexComponents(mol);
		List<Integer> toRemove = new ArrayList<Integer>();
		for (int i = 0; i < groups.size(); i++) {
			List<Integer> group = groups.get(i);
			if(group.size()>20 || group.size()>mol.getAtoms()/30 || mol.isAtomFlag(group.get(0), Molecule3D.LIGAND)) continue;
			toRemove.addAll(groups.get(i));
		}
		mol.deleteAtoms(toRemove);
	}
	
	public static void removeLigand(Molecule3D mol) {
		for(int i=mol.getAllAtoms()-1; i>=0; i--) {
			if(mol.isAtomFlag(i, Molecule3D.LIGAND)) {
				mol.deleteAtom(i);
			}
		}
	}
	
	public static Molecule3D crop(Molecule3D mol, Coordinates center, double radius) {
		Molecule3D res = new Molecule3D();
		Map<Integer, Integer> molToCrop = new HashMap<Integer, Integer>();
		Set<Integer> nonIndivuals = new TreeSet<Integer>();
		
		
		//Crop the atoms
		for(int i=0; i<mol.getAllAtoms(); i++) {
			Coordinates c = mol.getCoordinates(i);
			if(center.distance(c)<=radius) {
				int n = res.addAtom(mol, i);
				molToCrop.put(i, n);
				
				if(mol.getAllConnAtoms(i)!=0) nonIndivuals.add(n);
			}
		}
		
		//Add the bonds
		for(int i=0; i<mol.getAllBonds(); i++) {
			int mol1 = mol.getBondAtom(0, i);
			int mol2 = mol.getBondAtom(1, i);
			
			Integer crop1 = molToCrop.get(mol1);
			Integer crop2 = molToCrop.get(mol2);
			if(crop1!=null && crop2!=null) {
				 res.addBond(crop1.intValue(), crop2.intValue(), mol.getBondOrder(i));				 
			}			
		}
		
		//Remove the individual atoms
		for(int i=res.getAllAtoms()-1; i>=0; i--) {
			if(res.getAllConnAtoms(i)==0 && nonIndivuals.contains(i)) {
				res.deleteAtom(i);			
			}
		}
		res.setName(mol.getName()+ (res.getAllAtoms()<mol.getAllAtoms()?" Crop":""));		
		res.getAuxiliaryInfos().putAll(mol.getAuxiliaryInfos());
		return res;
	}
	
	public static Molecule3D cropLigand(Molecule3D mol, double radius) {
//		Molecule3D lig = extractLigand(mol);
//		if(lig.getAllAtoms()==0) return new Molecule3D(mol);
//		MoleculeGrid grid = new MoleculeGrid(lig);
//		Molecule3D res = new Molecule3D();
//		Map<Integer, Integer> molToCrop = new HashMap<Integer, Integer>();
//		Set<Integer> nonIndivuals = new TreeSet<Integer>();
//
//
//		//Crop the atoms
//		for(int i=0; i<mol.getAllAtoms(); i++) {
//			Coordinates c = mol.getCoordinates(i);
//			boolean add = false;
//			for (Iterator<Integer> iter = grid.getNeighbours(c, radius).iterator(); !add && iter.hasNext(); ) {
//				Integer a = iter.next();
//				if(lig.getCoordinates(a.intValue()).distanceSquared(c)<radius*radius) {
//					add = true;
//				}
//			}
//			if(add) {
//				int n = res.addAtom(mol, i);
//				molToCrop.put(i, n);
//
//				if(mol.getAllConnAtoms(i)!=0) nonIndivuals.add(n);
//			}
//		}
//
//		//Add the bonds
//		for(int i=0; i<mol.getAllBonds(); i++) {
//			int mol1 = mol.getBondAtom(0, i);
//			int mol2 = mol.getBondAtom(1, i);
//
//			Integer crop1 = molToCrop.get(mol1);
//			Integer crop2 = molToCrop.get(mol2);
//			if(crop1!=null && crop2!=null) {
//				 res.addBond(crop1.intValue(), crop2.intValue(), mol.getBondOrder(i));
//			}
//		}
//
//		//Remove the individual atoms
//		for(int i=res.getAllAtoms()-1; i>=0; i--) {
//			if(res.getAllConnAtoms(i)==0 && nonIndivuals.contains(i)) {
//				res.deleteAtom(i);
//			}
//		}
//		res.setName((mol.getName().indexOf(" Crop")<0? mol.getName(): mol.getName().substring(0, mol.getName().indexOf(" Crop"))) + (res.getAllAtoms()<mol.getAllAtoms()?" Crop"+((int)radius):""));
//		res.getAuxiliaryInfos().putAll(mol.getAuxiliaryInfos());
//
//		return res;
		return null;
	}	
	
	/**
	 * Marks a numbered ligand
	 * @param mol
	 * @param n >=0
	 */
	public static int markLigand(Molecule3D mol, int n) {
//		List<List<Integer>> l = StructureCalculator.getConnexComponents(mol);
//		List<Integer> groups = new ArrayList<Integer>();
//		mol.setAllAtomFlag(Molecule3D.LIGAND, false);
//		mol.setAllAtomFlag(Molecule3D.RIGID, true);
//		for(int i=0; i<l.size(); i++) {
//			List<Integer> group = l.get(i);
//			if(group.size()>70 || group.size()<7) continue;
//
//			for (int a :group) mol.setAtomFlag(a, Molecule3D.LIGAND, true);
//			int S = SurfaceCalculator.getLigandSurface(mol);
//			int SAS = SurfaceCalculator.getLigandSAS(mol);
//			double buried = (1.0*S-SAS)/S;
//			for (int a :group) mol.setAtomFlag(a, Molecule3D.LIGAND, false);
//
//			if(buried<.7) continue;	//Not a potential ligand
//			groups.add(i);
//		}
//		if(n<0 || n>=groups.size()) return 0;
//
//		for(int a: l.get(groups.get(n))) {
//			mol.setAtomFlag(a, Molecule3D.LIGAND, true);
//			mol.setAtomFlag(a, Molecule3D.RIGID, false);
//		}
//		return l.get(groups.get(n)).size();

		return -1;
	}

	
	/**
	 * Find the most likely ligands (in term of size)
	 */
	public static void markLigands(Molecule3D mol) {
		List<List<Integer>> l = StructureCalculator.getConnexComponents(mol);
		
		//double bestScore = -1;
		Set<Integer> potentials = new HashSet<Integer>();
		for(int i=0; i<l.size(); i++) {
			List<Integer> group = l.get(i);
			if(l.size()>1 && (group.size()>90 || group.size()<7)) continue;
			potentials.add(i);

		}

		
		for(int i=0; i<l.size(); i++) {
			List<Integer> group = l.get(i);		
			for(int j=0; j<group.size(); j++) {
				int atom = group.get(j);
				if(potentials.contains(i)) {
					mol.setAtomFlag(atom, Molecule3D.LIGAND, true);
					mol.setAtomFlag(atom, Molecule3D.RIGID, false);
				} else {
					mol.setAtomFlag(atom, Molecule3D.RIGID, true);
					mol.setAtomFlag(atom, Molecule3D.LIGAND, false);
				}
			}
		}
	}

	
	/**
	 * Find the most likely ligand (in term of size, rings, polar atoms)
	 * and mark it
	 * @param mol
	 */
	public static boolean markLigand(Molecule3D mol) {
//		List<List<Integer>> l = StructureCalculator.getConnexComponents(mol);
//		int bestGroup = -1;
//
//		//double bestScore = -1;
//		List<Integer> potentials = new ArrayList<Integer>();
//		if(bestGroup<0) {
//			for (int i = 0; i < mol.getAllAtoms(); i++) {
//				mol.setAtomFlag(i, Molecule3D.LIGAND, false);
//				mol.setAtomFlag(i, Molecule3D.RIGID, true);
//			}
//			for(int i=0; i<l.size(); i++) {
//				List<Integer> group = l.get(i);
//				if(group.size()>80 || group.size()<7) continue;
//				potentials.add(i);
//			}
//		}
//
//		if(bestGroup>=0) {
//			//Ok
//		} else if(potentials.size()==0) {
//			return false;
//		} else if(potentials.size()==1) {
//			bestGroup = potentials.get(0);
//		} else {
//			//what can be considered as the best ligand?
//			MoleculeGrid grid = new MoleculeGrid(mol);
//			double bestDensity = 0;
//			for(int potential : potentials) {
//				Set<Integer> count = new TreeSet<Integer>();
//
//				for(int a : l.get(potential)) {
//					count.addAll(grid.getNeighbours(mol.getCoordinates(a), 5));
//				}
//				double density = (double)count.size()/Math.sqrt(l.get(potential).size());
//				System.out.println("Ligand "+potential+" has " +l.get(potential).size()+" atoms  and a buried density of "+density);
//				if(density>bestDensity) {
//					bestDensity = density;
//					bestGroup = potential;
//				}
//			}
//
//		}
//
//		for(int i=0; i<l.size(); i++) {
//			List<Integer> group = l.get(i);
//			for(int j=0; j<group.size(); j++) {
//				int atom = group.get(j);
//				if(i==bestGroup) {
//					mol.setAtomFlag(atom, Molecule3D.LIGAND, true);
//					mol.setAtomFlag(atom, Molecule3D.RIGID, false);
//				} else {
//					mol.setAtomFlag(atom, Molecule3D.RIGID, true);
//					mol.setAtomFlag(atom, Molecule3D.LIGAND, false);
//				}
//			}
//		}
		return true;
	}
	
	public static Coordinates getCenter(Molecule3D mol) {
		Coordinates sum = new Coordinates();
		int n = 0;
		for(int i = 0; i<mol.getAllAtoms(); i++) {
			sum = sum.addC(mol.getCoordinates(i));
			n++;
		}				
		return n==0? null: sum.scaleC(1.0/n);
	}
	
	public static Coordinates getLigandCenter(Molecule3D mol) {
		Coordinates sum = new Coordinates();
		int n = 0;
		for(int i = 0; i<mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i)>1 && mol.isAtomFlag(i, Molecule3D.LIGAND)) {
				sum = sum.addC(mol.getCoordinates(i));
				n++;
			}
		}				
		return n==0? null: sum.scaleC(1.0/n);
	}
	
	public static Coordinates[] getLigandBounds(Molecule3D mol) {
		return getBounds(mol, Molecule3D.LIGAND);
	}

	public static Coordinates[] getBounds(Molecule3D mol) {
		return getBounds(mol, 0);
	}
	private static Coordinates[] getBounds(Molecule3D mol, int flags) {
		Coordinates min = new Coordinates(Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE);
		Coordinates max = new Coordinates(-Double.MAX_VALUE, -Double.MAX_VALUE, -Double.MAX_VALUE);
		for(int i = 0; i<mol.getAllAtoms(); i++) {
			if(flags==0 || (mol.getAtomFlags(i) & flags)>0) {
				min.x = Math.min(min.x, mol.getAtomX(i));
				min.y = Math.min(min.y, mol.getAtomY(i));
				min.z = Math.min(min.z, mol.getAtomZ(i));

				max.x = Math.max(max.x, mol.getAtomX(i));
				max.y = Math.max(max.y, mol.getAtomY(i));
				max.z = Math.max(max.z, mol.getAtomZ(i));
			}
		}
		if(min.x>max.x) return null;		
		return new Coordinates[]{min, max};
	}

	
	/**
	 * Gets the root mean square deviation between 2 positions.
	 * Only atomicNo>1 are considered
	 * @param mol1
	 * @param mol2
	 * @return
	 */
	public static double getLigandRMSD(Molecule3D mol1, Molecule3D mol2) {
		double sumSq = 0;
		int n = 0;
		for(int i1=0; i1<mol1.getAllAtoms(); i1++) {
			if(!mol1.isAtomFlag(i1, Molecule3D.LIGAND) || mol1.getAtomicNo(i1)<=1) continue;
			long s1 = getAtomHashkey(mol1, i1);
			double bestDistSq = Double.MAX_VALUE;
			int bestAtom = -1;
			for(int i2=0; i2<mol2.getAllAtoms(); i2++) {
				if(!mol2.isAtomFlag(i2, Molecule3D.LIGAND) || mol2.getAtomicNo(i2)!=mol1.getAtomicNo(i1)) continue;
				long s2 = getAtomHashkey(mol2, i2);
				if(s1!=s2) continue;
				double distSq = mol1.getCoordinates(i1).distSquareTo(mol2.getCoordinates(i2));
				if(distSq<bestDistSq) {
					bestDistSq = distSq;
					bestAtom = i2;
				}
			}
			if(bestAtom<0) return -1;
			sumSq += bestDistSq;			
			n++;
		}
		double res = Math.sqrt(sumSq/n);
		return res;
	}
	
	public static void translateLigand(Molecule3D mol, Coordinates v) {
		for(int i = 0; i<mol.getAllAtoms(); i++) {
			if(mol.isAtomFlag(i, Molecule3D.LIGAND)) {
				mol.setCoordinates(i, mol.getCoordinates(i).addC(v));
			}
		}				
	}
		
	public static void rotateLigand(Molecule3D mol, double angle, Coordinates normal, Coordinates center) {
		for(int i = 0; i<mol.getAllAtoms(); i++) {
			if(mol.isAtomFlag(i, Molecule3D.LIGAND)) {
				mol.setCoordinates(i, mol.getCoordinates(i).subC(center).rotate(normal, angle).addC(center));
			}
		}				
	}			
	
	/**
	 * 
	 * @param mol
	 * @param radius maximum moving range in Angstrom.
	 */
	public static void vibrateLigand(Molecule3D mol, double radius) {
		
		for(int i = 0; i<mol.getAllAtoms(); i++) {
			if(mol.isAtomFlag(i, Molecule3D.LIGAND)) {
				Coordinates v = new Coordinates(Math.random()-.5, Math.random()-.5, Math.random()-.5);
				if(v.dist()>0) {
					v = v.unitC().scaleC(Math.random()*radius);
				}
				mol.setCoordinates(i, mol.getCoordinates(i).addC(v));
			}
		}		
	}
		
	/**
	 * Return a long representing the atom and its neighbours
	 * 
	 * @param mol
	 * @param a
	 * @return
	 */
	public static long getAtomHashkey(Molecule3D mol, int a) {
		long res = 0;
		for(int i=0; i<mol.getAllConnAtoms(a); i++) {
			int a1 = mol.getConnAtom(a, i);
			if(mol.getAtomicNo(a1)<=1) continue;						
			res += (mol.getAtomicNo(a1)*10 + mol.getConnBondOrder(a, i) )*100;			
			for(int j=0; j<mol.getAllConnAtoms(a1); j++) {
				int a2 = mol.getConnAtom(a1, j);
				if(mol.getAtomicNo(a2)<=1) continue;						
				res += (mol.getAtomicNo(a2)*10 + mol.getConnBondOrder(a1, j));							
			}
		}
		res = res*100 + mol.getAtomicNo(a);
		return res;
	}
	
//	public static Coordinates[] getLigandCoordinates(Molecule3D mol) {
//		Coordinates[] coords = new Coordinates[mol.getNMovables()];
//		System.arraycopy(mol.getCoordinates(), 0, coords, 0, coords.length);
//		return coords;
//	}
	
//	public static void setLigandCoordinates(Molecule3D mol, Coordinates[] coords) {
//		setLigandCoordinates(mol, coords, true);
//	}
	
//	public static void setLigandCoordinates(Molecule3D mol, Coordinates[] coords, boolean copyReferences) {
//		int l = Math.min(mol.getNMovables(), Math.min(mol.getCoordinates().length, coords.length));
//		if(copyReferences) {
//			System.arraycopy(coords, 0, mol.getCoordinates(), 0, l);
//		} else
//			for (int i = 0; i < l; i++) {
//				mol.getCoordinates(i).x = coords[i].x;
//				mol.getCoordinates(i).y = coords[i].y;
//				mol.getCoordinates(i).z = coords[i].z;
//			}
//	}
	
	public static int makeProteinFlexible(Molecule3D mol, Coordinates center, double radius, boolean keepBackboneRigid) {
//		boolean[] inBackbone = keepBackboneRigid? StructureCalculator.getBackbones(mol, 0): null;
//
//
//		Set<Integer> neighbours = new TreeSet<Integer>();
//
//		if(center==null) {
//			MoleculeGrid grid = new MoleculeGrid(mol);
//			for (int i = 0; i < mol.getAllAtoms(); i++) {
//				if(!mol.isAtomFlag(i, Molecule3D.LIGAND)) continue;
//				neighbours.addAll(grid.getNeighbours(mol.getCoordinates(i), radius));
//			}
//
//		} else {
//			for (int i = 0; i < mol.getAllAtoms(); i++) {
//				if(!mol.isAtomFlag(i, Molecule3D.LIGAND) && mol.getCoordinates(i).distanceSquared(center)<radius*radius) neighbours.add(i);
//			}
//		}
//
//		int count = 0;
//		for(int a:neighbours) {
//			if(inBackbone!=null && inBackbone[a]) continue;
//			mol.setAtomFlag(a, Molecule3D.RIGID, false);
//			count++;
//		}
//		return count;

		return -1;
	}
	
	
	public static void makeProteinRigid(Molecule3D mol) {
		for (int a = mol.getAllAtoms()-1; a >= 0; a--) {
			if(!mol.isAtomFlag(a, Molecule3D.LIGAND)) {
				if(mol.getAtomicNo(a)<=1) mol.deleteAtom(a);
				mol.setAtomFlag(a, Molecule3D.RIGID, true);
			}
		}
	}
	
	/**
	 * Returns a array of 2 elements [number of donors, number of acceptors]
	 * @param m
	 * @return
	 */
	public static int[] getHDonorsAcceptors(ExtendedMolecule m) {
		int[] res = new int[2];
		for (int i = 0; i < m.getAllAtoms(); i++) {
			if(m.getAtomicNo(i)!=8 && m.getAtomicNo(i)!=7) continue;
			if(m.getAllHydrogens(i)>0) res[0]++;
			res[1]++; 
		}		
		return res;
	}
	
	public static boolean isDonor(Molecule3D m, int a) {
		int atomicNo = m.getAtomicNo(a);
		boolean isPolar = (atomicNo==7 || atomicNo==8 || atomicNo==15 || atomicNo==16);
		if(!isPolar) return false;
		return getExplicitHydrogens(m, a)>0 || getImplicitHydrogens(m, a)>0;
	} 
	
	public static boolean isAcceptor(Molecule3D m, int a) {
//		int atomicNo = m.getAtomicNo(a);
//		boolean isPolar = (atomicNo==7 || atomicNo==8 || atomicNo==15);
//		if(!isPolar) return false;
//
//		if(StructureCalculator.connected(m, a, -1, 2)>=0) return true;
//		if(StructureCalculator.connected(m, a, 0, 1)>=0) return true;
//		if(MM2Parameters.getNLonePairs(m, a)>0) return true;
//
//
//
		return false;
	}


	/**
	 * Return a List of int[] representing all the atom-atom pairs having
	 * an intermolecular interactions (distance close to sum of VDW)
	 *   
	 * @param
	 * @return
	 */
//	public static List<int[]> getHBonds(Molecule3D m) {
//		int[] donor = new int[m.getAllAtoms()];
//		int[] acceptor = new int[m.getAllAtoms()];
//		for (int i = 0; i < m.getAllAtoms(); i++) {
//			if(m.getAtomicNo(i)!=8 && m.getAtomicNo(i)!=7 && m.getAtomicNo(i)!=16) continue;
//			donor[i] = StructureCalculator.getImplicitHydrogens(m, i)+StructureCalculator.getExplicitHydrogens(m, i);
//			acceptor[i] = m.getAtomicNo(i)==8?2:m.getAtomicNo(i)==7?1:0; //Approximative
//		}
//
//
//		List<int[]> res = new ArrayList<int[]>();
//		MoleculeGrid grid = new MoleculeGrid(m);
//		for (int i = 0; i < m.getNMovables(); i++) {
//			if(!m.isAtomFlag(i, Molecule3D.LIGAND)) continue;
//			if(donor[i]==0 && acceptor[i]==0) continue;
//
//			Set<Integer> neighbours = grid.getNeighbours(m.getCoordinates(i), 5);
//			for (Iterator<Integer> iter = neighbours.iterator(); iter.hasNext();) {
//				int a =  ((Integer)iter.next()).intValue();
//				if(m.isAtomFlag(a, Molecule3D.LIGAND)) continue;
//				if(!((donor[i]>0 && acceptor[a]>0) || (donor[a]>0 && acceptor[i]>0))) continue;
//
//				double d = Math.sqrt(m.getCoordinates(a).distSquareTo(m.getCoordinates(i)));
//				double vdw = VDWRadii.VDW_RADIUS[m.getAtomicNo(a)]+VDWRadii.VDW_RADIUS[m.getAtomicNo(i)];
//
//				if(d>vdw-.5 && d<vdw+.5) { //H-Bonds
//					boolean hbond = false;
//					for (int j = 0; j < m.getAllConnAtoms(i); j++) {
//						if(m.getAtomicNo(m.getConnAtom(i, j))<=1) continue;
//						for (int k = 0; k < m.getAllConnAtoms(a); k++) {
//							if(m.getAtomicNo(m.getConnAtom(a, k))<=1) continue;
//
//							Coordinates c1 = m.getCoordinates(i).subC(m.getCoordinates(m.getConnAtom(i, j)));
//							Coordinates c2 = m.getCoordinates(a).subC(m.getCoordinates(m.getConnAtom(a, k)));
//							double angle = c1.getAngle(c2);
//							if(Math.abs(2*Math.PI/3-angle)<Math.PI/10) hbond = true;
//							if(Math.abs(Math.PI/3-angle)<Math.PI/10) hbond = true;
//						}
//					}
//					if(hbond) {
//						res.add(new int[]{i, a});
//						/*if(donor[i]>0 && acceptor[a]>0) {
//							donor[i]--; acceptor[a]--;
//						} else {
//							donor[a]--; acceptor[i]--;
//						}*/
//					}
//				}
//			}
//
//		}
//		return res;
//	}
	
//	public static final void createArtificialLigand(Molecule3D protein, Coordinates c, double radius) {
//		MoleculeGrid grid = new MoleculeGrid(protein);
//
//		List<Coordinates> l = new ArrayList<Coordinates>();
//		l.add(c);
//		Set<Coordinates> seen = new TreeSet<Coordinates>();
//		double d = 2;
//		while(!l.isEmpty()) {
//			Coordinates cc = l.remove(0);
//			if(seen.contains(cc)) continue;
//			if(cc.distSquareTo(c)>radius*radius) continue;
//			seen.add(cc);
//			if(!grid.hasNeighbours(cc, 2.1)) {
//				int a = protein.addAtom(6);
//				protein.setCoordinates(a, cc);
//				protein.setAtomFlags(a, Molecule3D.LIGAND);
//				l.add(new Coordinates(cc.x-d, cc.y, cc.z));
//				l.add(new Coordinates(cc.x, cc.y-d, cc.z));
//				l.add(new Coordinates(cc.x, cc.y, cc.z-d));
//				l.add(new Coordinates(cc.x+d, cc.y, cc.z));
//				l.add(new Coordinates(cc.x, cc.y+d, cc.z));
//				l.add(new Coordinates(cc.x, cc.y, cc.z+d));
//			}
//		}
//	}
	

	/*
	public static void normalizeCoordinates(Molecule3D mol) {
		int a1 = StructureCalculator.getCenterAtom(mol);
		int a2 = mol.getConnAtom(a1, 0);
		int a3 = mol.getAllConnAtoms(a1)>1? mol.getConnAtom(a1, 1):
			mol.getAllConnAtoms(a2)>0? mol.getConnAtom(a2, 0): a2;
			System.out.println(a1+" "+a2+" "+a3);
		normalizeCoordinates(mol, a1, a2, a3);
	}
	public static void normalizeCoordinates(Molecule3D mol, int a1, int a2, int a3) {
		//Move a1 gravity to 0
		Coordinates G = mol.getCoordinates(a1);
		FFCalculator.translateLigand(mol, new Coordinates().sub(G));
		
		Coordinates ex = new Coordinates(1,0,0);
		Coordinates ey = new Coordinates(0,1,0);
		Coordinates ez = new Coordinates(0,0,1);

		//move the atom a2 on the X axis (rotation / Z): c0.y = 0
		Coordinates u = new Coordinates(mol.getCoordinates(a2).x, mol.getCoordinates(a2).y, 0);
		double angle = u.getAngle(ex);
		if(u.cross(ex).dot(ez)>0) angle=-angle;
		FFCalculator.rotateLigand(mol, angle, ez, new Coordinates());
		
		//c0.z = 0
		u = new Coordinates(mol.getCoordinates(a2).x, 0, mol.getCoordinates(a2).z);
		angle = u.getAngle(ex);
		if(u.cross(ex).dot(ey)>0) angle=-angle;
		FFCalculator.rotateLigand(mol, angle, ey, new Coordinates());

		
		//move the atom a3 on the Y axis (rotation / ex ): c1.z = 0
		//Coordinates e = mol.getCoordinates(0).sub(G).unit();
		u = new Coordinates(0, mol.getCoordinates(a3).y, mol.getCoordinates(a3).z); 
		angle = u.getAngle(ey);
		if(u.cross(ey).dot(ex)>0) angle=-angle;		
		FFCalculator.rotateLigand(mol, angle, ex, new Coordinates());
	}	
	*/
	
	
	////////////////////////////////////////////////////////////////
	public static void flagBackbone(Molecule3D mol) {
		boolean[] backbone = getBackbones(mol);
		for (int i = 0; i < backbone.length; i++) {
			mol.setAtomFlag(i, Molecule3D.BACKBONE, backbone[i]);
		}
	}

//	public static int[] getRotatableBonds(Molecule3D mol, boolean considerHydrogens) {
//		return getRotatableBonds(mol, -1, considerHydrogens);
//	}
	/**
	 * Returns the rotatable bonds (sorted). A rotatable bond is defined here 
	 * as a simple bond that is not inside a ring (and not connected to only hydrogen).
	 * @param mol
	 * @return
	 */
//	public static int[] getRotatableBonds(Molecule3D mol, int groupSeed, boolean considerHydrogens) {
//		int[] a2g = groupSeed>=0? getAtomToGroups(mol): null;
//
//		List<Integer> bonds = new ArrayList<Integer>();
//
//		List<Integer>[] atomToRings = mol.getAtomToRings();
//
//		for (int i = 0; i < mol.getAllBonds(); i++) {
//			if(mol.getBondOrder(i)!=1) continue;
//			int a1 = mol.getBondAtom(0, i);
//			int a2 = mol.getBondAtom(1, i);
//
//			if(a2g!=null && a2g[groupSeed]!=a2g[a1]) continue;
//
//			if(mol.isAtomFlag(a1, Molecule3D.RIGID) && mol.isAtomFlag(a2, Molecule3D.RIGID)) continue;
//			if(mol.getAtomicNo(a1)<=1 || mol.getAtomicNo(a2)<=1) continue;
//
//			//to be a rotatable bond, each atoms need to be connected to more atoms
//			int nOthers = 0;
//			for(int j=0; j<mol.getAllConnAtoms(a1); j++) {
//				int a = mol.getConnAtom(a1, j);
//				if(a==a2) continue;
//				if(((considerHydrogens && (mol.getAtomicNo(a1)==8 && mol.getAtomicNo(a1)==7)) || mol.getAtomicNo(a)>1)) nOthers++;
//			}
//			if(nOthers==0) continue;
//
//			nOthers = 0;
//			for(int j=0; j<mol.getAllConnAtoms(a2); j++) {
//				int a = mol.getConnAtom(a2, j);
//				if(a==a1) continue;
//				if(((considerHydrogens && (mol.getAtomicNo(a2)==8 && mol.getAtomicNo(a2)==7)) || mol.getAtomicNo(a)>1)) nOthers++;
//			}
//			if(nOthers==0) continue;
//
//			//Amides groups RN(H)C(=O)R are not rotatable
//			if(mol.getAtomicNo(a1)==7 && mol.getConnAtoms(a1)==2 && mol.getAtomicNo(a2)==6 && connected(mol, a2, 8, 2)>=0) continue;
//			if(mol.getAtomicNo(a2)==7 && mol.getConnAtoms(a2)==2 && mol.getAtomicNo(a1)==6 && connected(mol, a1, 8, 2)>=0) continue;
//
//			//Nitrile groups are not N#R rotatable
//			if(connected(mol, a1, 7, 3)>=0) continue;
//			if(connected(mol, a2, 7, 3)>=0) continue;
//
//			//Intersection of rings, if the atoms belongs to the same rings, they cannot be freely rotated
//			Set<Integer> set = new TreeSet<Integer>();
//			set.addAll(atomToRings[a1]);
//			set.retainAll(atomToRings[a2]);
//
//			if(set.size()>0) continue;
//
//			bonds.add(i);
//		}
//		int[] res = new int[bonds.size()];
//		for (int i = 0; i < res.length; i++) res[i] = bonds.get(i);
//
//		return res;
//	}
	
	public static int connected(Molecule3D mol, int a, int atomicNo, int bondOrder) {
		for(int i=0; i<mol.getAllConnAtoms(a); i++) {
			int atm = mol.getConnAtom(a, i);
			if(atomicNo>=0 && mol.getAtomicNo(atm)!=atomicNo) continue;
			if(bondOrder>0 && mol.getConnBondOrder(a, i)!=bondOrder) continue;
			return atm;
		}
		return -1;
	}

//	public static int getNRotatableBonds(Molecule3D mol) {
//		return getRotatableBonds(mol, false).length;
//	}

//	public static boolean fixParities(Molecule3D mol, StereoMolecule ref) {
//		return fixOrCheckParities(mol, ref, true);
//	}
//	public static boolean checkParities(Molecule3D mol, StereoMolecule ref) {
//		return fixOrCheckParities(mol, ref, false)==false;
//	}
	/**
	 * Make sure the parity of mol match the one from ref, otherwise fix it.
	 * 
	 * @param mol
	 * @param
	 * @return a boolean - true if something has been changed
	 */
//	private static boolean fixOrCheckParities(Molecule3D mol, StereoMolecule ref, boolean allowChange) {
//		mol.reorderAtoms();
//		Molecule3D copy = new Molecule3D(mol);
//		StructureCalculator.deleteHydrogens(copy);
//		StereoMolecule sm = copy.toStereoMolecule();
//
//		ref.ensureHelperArrays(ExtendedMolecule.cHelperSymmetrySimple);
//		sm.ensureHelperArrays(ExtendedMolecule.cHelperSymmetrySimple);
//
//
//
//		int N = Math.min(sm.getAtoms(), ref.getAtoms());
//
//		//Make sure the 2 mols matches
//		for (int i = 0; i < N; i++) {
//			if(sm.getAtomicNo(i)!=ref.getAtomicNo(i)) {
//				System.err.println("The 2 molecules don't match (atomic# "+i+")");
//				return false;
//			}
//			if(sm.getConnAtoms(i)!=ref.getConnAtoms(i)) {
//				System.err.println("The 2 molecules don't match (conn# "+i+") "+sm.getConnAtoms(i)+"<>"+ref.getConnAtoms(i));
//				return false;
//			}
//			for (int j = 0; j < sm.getConnAtoms(i); j++) {
//				if(sm.getConnAtom(i, j)!=ref.getConnAtom(i, j)) {
//					System.err.println("The 2 molecules don't match (conn# "+i+","+j+") "+sm.getConnAtom(i, j)+"<>"+ref.getConnAtom(i, j));
//					return false;
//
//				}
//			}
//
//			if(sm.getSymmetryRank(i)!=ref.getSymmetryRank(i)) {
//				System.err.println("Symmetric molecules-> don't fix  parity");
//				return false;
//			}
//		}
//
//		//Compare the parity
//		int[] groupNo = new int[10];
//		for (int i = 0; i < N; i++) {
//			if(ref.getAtomESRGroup(i)>=0 && ref.getAtomESRGroup(i)<10) groupNo[ref.getAtomESRGroup(i)]++;
//		}
////		for (int i = 0; i < 10; i++) {
////			if(groupNo[i]>0) return false; //TODO: Cannot fix ESR group for the moment;
////		}
//
//		boolean shouldChange = false;
//		for (int i = 0; i < N; i++) {
//			if(sm.getConnAtoms(i)<3 || ref.getAtomParity(i)==0 || ref.getAtomParity(i)>2) continue;
//			if(sm.getAtomCIPParity(i)==ref.getAtomCIPParity(i)) continue; //Same parity -> ok
////			if(ref.getAtomESRGroup(i)>=0 && ref.getAtomESRGroup(i)<10 && groupNo[ref.getAtomESRGroup(i)]>0) {
////				//TODO: ESR group, harder to check, not done for the moment
////				continue;
////			}
////			System.out.println("Parity does not match "+sm.getAtomParity(i)+"<>"+ref.getAtomParity(i));
////			System.out.println("CIP does not match "+sm.getAtomCIPParity(i)+"<>"+ref.getAtomCIPParity(i));
//			shouldChange = true;
//			if(allowChange) {
//				//mirror image, keep 3 neighboring atoms fixed and use the mirror image
//				int a1 = sm.getConnAtom(i, 0);
//				int a2 = sm.getConnAtom(i, 1);
//				int a3 = sm.getConnAtom(i, 2);
//				Set<Integer> seen = new HashSet<Integer>();
//				seen.add(a1); seen.add(a2); seen.add(a3);
//				dfs(mol, i, seen);
//				seen.remove(a1); seen.remove(a2); seen.remove(a3);
//
//				for (int j:seen) {
//					mol.setCoordinates(j, Coordinates.getMirror(mol.getCoordinates(j), mol.getCoordinates(a1), mol.getCoordinates(a2), mol.getCoordinates(a3)));
//					if(j<mol.getAllAtoms()) copy.setCoordinates(j, mol.getCoordinates(j));
//				}
//
//				sm = copy.toStereoMolecule();
//				sm.ensureHelperArrays(ExtendedMolecule.cHelperParities);
//			}
//		}
//		//if(changed) StructureCalculator.deleteHydrogens(mol);
//
//		return shouldChange;
//	}
	
	
	public static void rotateBond(Molecule3D mol, int a2, int a3, double angle) {
		Set<Integer> seen = new HashSet<Integer>();
		seen.add(a2);			
		dfs(mol, a3, seen);
		Coordinates normal = mol.getCoordinates(a3).subC(mol.getCoordinates(a2)).unitC();
		for (int j: seen) {
			Coordinates c = mol.getCoordinates(j).subC(mol.getCoordinates(a3));
			mol.setCoordinates(j, c.rotate(normal, angle).addC(mol.getCoordinates(a3)));
		}
	}
	
	public static boolean removeLonePairs(Molecule3D mol) {
		boolean changed = false;
		for(int i=0; i<mol.getAllAtoms(); i++) {
			if(mol.getAtomicNo(i)==0) {mol.deleteAtom(i);i--; changed=true;}
		}
		return changed;
	}
	

	public static final double[][] getDistanceMatrix(Molecule3D m1, Molecule3D m2) {
		double dist[][] = new double[m1.getAllAtoms()][m2.getAllAtoms()];
		for (int a1 = 0; a1 < m1.getAllAtoms(); a1++) {
			for (int a2 = 0; a2 < m2.getAllAtoms(); a2++) {
				dist[a1][a2] = m1.getCoordinates(a1).distance(m2.getCoordinates(a2));					
			}
		}
		return dist;
		
	}
		
//	public static double getSimilarity3D(Molecule3D m1, Molecule3D m2) {
//		return 1-Math.sqrt((1-getSimilarity3DIntransitive(m1, m2, null, null, 0))*(1-getSimilarity3DIntransitive(m2,m1, null, null, 0)));
//	}
//
//	public static double getSimilarity3D(Molecule3D m1, Molecule3D m2, Coordinates center, double radius) {
//		return 1-Math.sqrt((1-getSimilarity3DIntransitive(m1, m2, null, center, radius))*(1-getSimilarity3DIntransitive(m2, m1, null, center, radius)));
//	}
	
	/**
	 * (Very) Simple 3D similarity evaluation
	 * @param m1
	 * @param m2
	 * @return the % of similarities
	 */
//	protected static double getSimilarity3DIntransitive(Molecule3D m1, Molecule3D m2, MoleculeGrid grid2, Coordinates center, double radius) {
//		double sim = 0;
//		if(grid2==null) grid2 = new MoleculeGrid(m2, 3);
//		int count = 0;
//		for (int i = 0; i < m1.getAllAtoms(); i++) {
//
//			if(center!=null && m1.getCoordinates(i).distanceSquared(center)>radius*radius) continue;
//			count++;
//			int neighbour = grid2.getClosestNeighbour(m1.getCoordinates(i), 1.5);
//			if(neighbour<0) continue;
//			double d = m1.getCoordinates(i).distance(m2.getCoordinates(neighbour));
//			if(d>1.5) continue;
//
//			double coeff;
//			if(m1.getAtomicNo(i)==m2.getAtomicNo(neighbour)) coeff=1;
//			else if(m1.getAtomicNo(i)!=6 && m2.getAtomicNo(neighbour)!=6) coeff=.6;
//			else coeff=.3;
//
//			sim+= coeff * (1.0-Math.max(d-.5, 0))/1.0;
//
//		}
//		return sim / count;
//	}
	
	/**
	 * 
	 * @param m1
	 * @param m2
	 * @return an int[][] 0-> number of closer atom matching, 2->number of atom clashes  
	 */
	public static final int[] getOverlap(Molecule3D m1, Molecule3D m2) {
		int N1 = m1.getAllAtoms();
		int N2 = m2.getAllAtoms();
		double[][] dist = StructureCalculator.getDistanceMatrix(m1, m2);
		for (int a1 = 0; a1 < N1; a1++) {
			for (int a2 = 0; a2 < N2; a2++) {
				if(m1.getAtomicNo(a1)!=m2.getAtomicNo(a2)) dist[a1][a2]+=.3; 
			}
		}
		
		boolean seen[] = new boolean[N2];
		int sum1 = 0;
		int sum2 = 0;
		while(true) {
			int[] min = null;
			for (int a2 = 0; a2 < N2; a2++) {
				if(m2.getAtomicNo(a2)<=1) continue;
				if(seen[a2]) continue;
				for (int a1 = 0; a1 < N1; a1++) {
					if(m1.getAtomicNo(a1)<=1) continue;
					if(min==null || dist[a1][a2]<dist[min[0]][min[1]]) min = new int[] {a1, a2};
				}
			}
			if(min==null) break;

			seen[min[1]] = true;
			double d = dist[min[0]][min[1]];
			if(d<.75) {
				sum1++; 
			} else if(d<1.4) {
				sum2++;
			}
			
		}
		return new int[] {sum1, sum2};
	}
	
//	public static String testLipinskiRules(StereoMolecule m) {
//		return testLipinskiRules(m, 600, 5, 10, 10, -1, 5.5, -20, 20);
//	}

//	public static String testLipinskiRules(StereoMolecule m, int weight, int donors, int acceptors, int rot, double clogpmin, double clogpmax, double clogsmin, double clogsmax) {
//		if(m.getMolweight()>weight) return "MolWeight";
//		int[] donorsAcceptors = getHDonorsAcceptors(m);
//		if(donorsAcceptors[0]>donors) return "donors ("+donorsAcceptors[0]+")";
//		if(donorsAcceptors[1]>acceptors) return "acceptors ("+donorsAcceptors[1]+")";
//		int rota = getNRotatableBonds(new Molecule3D(m));
//		if(rota>rot) return "rot ("+rota+")";
//		double cLogP = new CLogPPredictor().assessCLogP(m);
//		if(cLogP<clogpmin) return "cLogP ("+(int)cLogP+")";
//		if(cLogP>clogpmax) return "cLogP ("+(int)cLogP+")";
//		double cLogS = new SolubilityPredictor().assessSolubility(m);
//		if(cLogS<clogsmin) return "cLogS ("+(int)cLogS+")";
//		if(cLogS>clogsmax) return "cLogS ("+(int)cLogS+")";
//		//NastyFunctionDetector nst = new NastyFunctionDetector();
//		//if(nst.getNastyFunctionCount(m, null, null)>0) return nst.getNastyFunctionString();
//		return null;
//	}

//	public static boolean testLipinski(StereoMolecule m) {
//		return testLipinskiRules(m)==null;
//	}
//	public static boolean testLipinskiWide(StereoMolecule m) {
//		return testLipinskiRules(m, 650, 7, 13, 14, -3, 7.5,-20,20)==null;
//	}
	
	/**
	 * Returns the number of paths of the given length (gives an idea of the number of branches) 
	 * @param mol
	 * @param length
	 * @return
	 */
	public static int getMolecularPathCount(Molecule3D mol, int length) {
		int[][] dists = getNumberOfBondsBetweenAtoms(mol, length, null);
		int count = 0;
		for (int i = 0; i < dists.length; i++) {
			for (int j = i; j < dists.length; j++) {
				if(dists[i][j]==length) count++;				
			}
		}
		return count;
	}
	
	
	/**
	 * Find all fragment occurences in mol (only the movable atoms are considered).
	 * Return a List<int[]> of occurences
	 * list.get(n)[fragmentAtom] points to the n occurence of fragments and is equal to the matching atomIndex in mol     
	 * @param fragment
	 * @param mol
	 * @return
	 */
	public static List<int[]> getFragmentMatches(Molecule3D fragment, Molecule3D mol) {
		List<int[]> res = new ArrayList<int[]>();		
		getFragmentMatchesRec(fragment, mol, new int[fragment.getAllAtoms()], 0, res);
		return res;
	}
	
	private static void getFragmentMatchesRec(Molecule3D fragment, Molecule3D mol, int[] match, int fragmentMatchIndex, List<int[]> res) {
		//System.out.println("StructureCalculator.getFragmentMatchesRec() "+matchIndex+" :"+fragment.getAllAtoms()+" "+mol.getNMovables());
		if(fragmentMatchIndex>=fragment.getAllAtoms()) {
			//Terminal Condition, we found a match
			res.add(match.clone());
		} else {
			//Find an atom in mol that matches fragment.getAtom(fragmentMatchIndex)			
			atomLoop: for (int a = 0; a<mol.getNMovables(); a++) {
				//it must have the same atomicNo
				if(fragment.getAtomicNo(fragmentMatchIndex)!=6 && fragment.getAtomicNo(fragmentMatchIndex)!=mol.getAtomicNo(a)) continue atomLoop;
				
				//it must not be already selected
				for (int i = 0; i < fragmentMatchIndex; i++) {
					if(match[i]==a) continue atomLoop;
				}
									
				//it must have the same connection
				for (int j = 0; j < fragment.getConnAtoms(fragmentMatchIndex); j++) {					
					int connAtom = fragment.getConnAtom(fragmentMatchIndex, j);
					if(connAtom>=fragmentMatchIndex) continue;	//Continue to next connection				
					int b = mol.getBond(a, match[connAtom]);
					if(b<0) continue atomLoop;
				}
		

				match[fragmentMatchIndex] = a;
				getFragmentMatchesRec(fragment, mol, match, fragmentMatchIndex+1, res);
				
			}
		}		
	}

	/**
	 * 
	 * @param mols
	 * @param tresh
	 * @return
	 */
//	public static int cluster3D(List<Molecule3D> mols, final double tresh) {
//		return cluster3D(mols, tresh, null, 0);
//	}
	
	/**
	 * Clusters the molectules based solely on the atomtypes and their 3D coordinates. (no superimposition performed)
	 * 
	 * @param sorted
	 * @param tresh
	 * @return
	 */
//	public static int cluster3D(List<Molecule3D> sorted, final double tresh, Coordinates center, double radius) {
//		int cluster = 0;
//		for (int i = 0; i < sorted.size(); i++) {
//			if(sorted.get(i).getAuxiliaryInfo("cluster3D")!=null) continue;
//			sorted.get(i).setAuxiliaryInfo("cluster3D", ++cluster);
//			sorted.get(i).setAuxiliaryInfo("prefered", "true");
//			//MoleculeGrid grid = new MoleculeGrid(sorted.get(i));
//			for (int j = i+1; j < sorted.size(); j++) {
//				if(sorted.get(j).getAuxiliaryInfo("cluster3D")!=null) continue;
//				//double v = getSimilarity3DIntransitive(sorted.get(j), sorted.get(i), grid, center, radius);
//				double v = getSimilarity3D(sorted.get(j), sorted.get(i), center, radius);
//				if(v>tresh) {
//					sorted.get(j).setAuxiliaryInfo("cluster3D", cluster);
//					sorted.get(j).setAuxiliaryInfo("prefered", "false");
//				}
//			}
//		}
//		return cluster;
//	}

	/**
	 * Return a List of int[] representing all the atom-atom pairs having
	 * an intermolecular interactions (distance close to sum of VDW)
	 *
	 * @param m
	 * @return
	 */
	public static List<int[]> getHBonds(Molecule3D m) {
		int[] donor = new int[m.getAllAtoms()];
		int[] acceptor = new int[m.getAllAtoms()];
		for (int i = 0; i < m.getAllAtoms(); i++) {
			if(m.getAtomicNo(i)!=8 && m.getAtomicNo(i)!=7 && m.getAtomicNo(i)!=16) continue;
			donor[i] = StructureCalculator.getImplicitHydrogens(m, i)+StructureCalculator.getExplicitHydrogens(m, i);
			acceptor[i] = m.getAtomicNo(i)==8?2:m.getAtomicNo(i)==7?1:0; //Approximative
		}


		List<int[]> res = new ArrayList<int[]>();
		MoleculeGrid grid = new MoleculeGrid(m);
		for (int i = 0; i < m.getNMovables(); i++) {
			if(!m.isAtomFlag(i, Molecule3D.LIGAND)) continue;
			if(donor[i]==0 && acceptor[i]==0) continue;

			Set<Integer> neighbours = grid.getNeighbours(m.getCoordinates(i), 5);
			for (Iterator<Integer> iter = neighbours.iterator(); iter.hasNext();) {
				int a =  ((Integer)iter.next()).intValue();
				if(m.isAtomFlag(a, Molecule3D.LIGAND)) continue;
				if(!((donor[i]>0 && acceptor[a]>0) || (donor[a]>0 && acceptor[i]>0))) continue;

				double d = Math.sqrt(m.getCoordinates(a).distSquareTo(m.getCoordinates(i)));
				double vdw = VDWRadii.VDW_RADIUS[m.getAtomicNo(a)]+VDWRadii.VDW_RADIUS[m.getAtomicNo(i)];

				if(d>vdw-.5 && d<vdw+.5) { //H-Bonds
					boolean hbond = false;
					for (int j = 0; j < m.getAllConnAtoms(i); j++) {
						if(m.getAtomicNo(m.getConnAtom(i, j))<=1) continue;
						for (int k = 0; k < m.getAllConnAtoms(a); k++) {
							if(m.getAtomicNo(m.getConnAtom(a, k))<=1) continue;

							Coordinates c1 = m.getCoordinates(i).subC(m.getCoordinates(m.getConnAtom(i, j)));
							Coordinates c2 = m.getCoordinates(a).subC(m.getCoordinates(m.getConnAtom(a, k)));
							double angle = c1.getAngle(c2);
							if(Math.abs(2*Math.PI/3-angle)<Math.PI/10) hbond = true;
							if(Math.abs(Math.PI/3-angle)<Math.PI/10) hbond = true;
						}
					}
					if(hbond) {
						res.add(new int[]{i, a});
						/*if(donor[i]>0 && acceptor[a]>0) {
							donor[i]--; acceptor[a]--;
						} else {
							donor[a]--; acceptor[i]--;
						}*/
					}
				}
			}

		}
		return res;
	}


}
