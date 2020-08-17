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
 * @author Joel Freyss
 */
package com.actelion.research.chem.io.pdb.converter;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;

import java.util.Set;
import java.util.TreeSet;

/**
 * Class used to speed up the calculation of neighbours by creating a grid. 
 * Each atom is indexed according to its 3D coordinates.
 * 
 */
public class MoleculeGrid {
	
	protected final StereoMolecule mol;
	protected final double gridWidth;
	protected final Coordinates min;
	protected final Coordinates max;
	protected final int[] gridSize = new int[3];
	protected final Set<Integer>[][][] grid;
	

	
	public MoleculeGrid(StereoMolecule mol) {
		this(mol, 1.1,new Coordinates(0.0,0.0,0.0));
	}

	/**
	 * Creates the Grid: Complexity O(nAtoms)
	 * @param mol
	 */
	@SuppressWarnings("unchecked")
	public MoleculeGrid(StereoMolecule mol, double gridWidth, Coordinates extension) {
		this.mol = mol;
		this.gridWidth = gridWidth;
		//1. Find the Molecule's bounds
		Coordinates[] bounds = GeometryCalculator.getBounds(mol);
		
		min = bounds[0];
		max = bounds[1];
		
		min.x -= extension.x;
		min.y -= extension.y;
		min.z -= extension.z;
		
		max.x += extension.x;
		max.y += extension.y;
		max.z += extension.z;
		
		//2. Creates the grid 
		gridSize[0] = (int)((max.x-min.x)/gridWidth)+1;
		gridSize[1] = (int)((max.y-min.y)/gridWidth)+1;
		gridSize[2] = (int)((max.z-min.z)/gridWidth)+1;
		grid = new Set[Math.max(0, gridSize[0])][Math.max(0, gridSize[1])][Math.max(0, gridSize[2])];		
		
		//3. Put each atom in the grid

		int atoms = mol.getAtoms();
		for (int i = 0; i < atoms; i++) {
			int x = (int)((mol.getAtomX(i)-min.x)/gridWidth);
			int y = (int)((mol.getAtomY(i)-min.y)/gridWidth);
			int z = (int)((mol.getAtomZ(i)-min.z)/gridWidth);			
			if(grid[x][y][z]==null) grid[x][y][z] = new TreeSet<Integer>();
			grid[x][y][z].add(i); 
		}
	}
	
	/**
	 * Gets a Set of all neigbouring atoms. This class ensures that at least
	 * all atoms within maxDist are returned (+more)
	 * @param c
	 * @param maxDist
	 * @return
	 */
	public Set<Integer> getNeighbours(Coordinates c, double maxDist) {
		int radius = (int)(maxDist / gridWidth)  + 1;
		int x = (int)((c.x-min.x)/gridWidth);
		int y = (int)((c.y-min.y)/gridWidth);
		int z = (int)((c.z-min.z)/gridWidth);			

		Set<Integer> res = new TreeSet<Integer>();
		for (int i = Math.max(0, x-radius); i<= Math.min(gridSize[0]-1, x+radius); i++) {
			for (int j = Math.max(0, y-radius); j<= Math.min(gridSize[1]-1, y+radius); j++) {
				for (int k = Math.max(0, z-radius); k<= Math.min(gridSize[2]-1, z+radius); k++) {
				//	int dx = i-x;
				//	int dy = j-y;
				//	int dz = k-z;
				//	if(dx*dx+dy*dy+dz*dz>3+radius*radius) continue;
					if(grid[i][j][k]!=null) res.addAll(grid[i][j][k]);
				}
			}			
		}
		return res;
	}

	
	public Set<Integer> getNeighbours(Coordinates c, double maxDist, boolean enforceDist) {
		return getNeighbours(c, maxDist, enforceDist, -1);
	}
	
	/**
	 * Gets a Set of all neigbouring atoms. This class ensures that at least
	 * all atoms within maxDist are returned (+more)
	 * @param c
	 * @param maxDist
	 * @return
	 */
	public Set<Integer> getNeighbours(Coordinates c, double maxDist, boolean enforceDist, int requiredFlags) {
		int radius = (int)(maxDist / gridWidth)  + 1;
		int x = (int)((c.x-min.x)/gridWidth);
		int y = (int)((c.y-min.y)/gridWidth);
		int z = (int)((c.z-min.z)/gridWidth);			

		Set<Integer> res = new TreeSet<Integer>();
		for (int i = Math.max(0, x-radius); i<= Math.min(gridSize[0]-1, x+radius); i++) {
			for (int j = Math.max(0, y-radius); j<= Math.min(gridSize[1]-1, y+radius); j++) {
				for (int k = Math.max(0, z-radius); k<= Math.min(gridSize[2]-1, z+radius); k++) {
					if(grid[i][j][k]!=null) {
						if(enforceDist) {
							for(int elt: grid[i][j][k]) {

// TODO check whether this was used somewhere
//                  			if(requiredFlags>=0 && !mol.isAtomFlag(elt, requiredFlags)) continue;

								if(mol.getCoordinates(elt).distSquareTo(c)>maxDist*maxDist ) continue;
								res.add(elt);
							}
						} else {
							res.addAll(grid[i][j][k]);
						}
					}
				}
			}			
		}
		return res;
	}	
	
	public boolean hasNeighbours(Coordinates c, double maxDist) {
		int radius = (int)(maxDist / gridWidth)  + 1;
		int x = (int)((c.x-min.x)/gridWidth);
		int y = (int)((c.y-min.y)/gridWidth);
		int z = (int)((c.z-min.z)/gridWidth);			

		for (int i = Math.max(0, x-radius); i<= Math.min(gridSize[0]-1, x+radius); i++) {
			for (int j = Math.max(0, y-radius); j<= Math.min(gridSize[1]-1, y+radius); j++) {
				for (int k = Math.max(0, z-radius); k<= Math.min(gridSize[2]-1, z+radius); k++) {
					if(grid[i][j][k]!=null) {						
						for(int a : grid[i][j][k]) {
							if(mol.getCoordinates(a).distSquareTo(c)<=maxDist*maxDist ) {
								return true;
							}
						}
					}
				}
			}			
		}
		return false;			
	}
	
	
	/**
	 * Gets a Set of all neigbouring atoms. This class ensures that at least
	 * all atoms within maxDist are returned (+more)
	 * @param bounds
	 * @param maxDist
	 * @return
	 */
	public Set<Integer> getNeighbours(Coordinates[] bounds, double maxDist) {
		int radius = (int)(maxDist / gridWidth)  + 1;
		int x1 = (int)((bounds[0].x-min.x)/gridWidth);
		int y1 = (int)((bounds[0].y-min.y)/gridWidth);
		int z1 = (int)((bounds[0].z-min.z)/gridWidth);
					
		int x2 = (int)((bounds[1].x-min.x)/gridWidth);
		int y2 = (int)((bounds[1].y-min.y)/gridWidth);
		int z2 = (int)((bounds[1].z-min.z)/gridWidth);			

		int xMin = Math.min(x1, x2);
		int xMax = Math.max(x1, x2);
		int yMin = Math.min(y1, y2);
		int yMax = Math.max(y1, y2);
		int zMin = Math.min(z1, z2);
		int zMax = Math.max(z1, z2);
			
		Set<Integer> res = new TreeSet<Integer>();
		for (int i = Math.max(0, xMin-radius); i<= Math.min(gridSize[0]-1, xMax+radius); i++) {
			for (int j = Math.max(0, yMin-radius); j<= Math.min(gridSize[1]-1, yMax+radius); j++) {
				for (int k = Math.max(0, zMin-radius); k<= Math.min(gridSize[2]-1, zMax+radius); k++) {
					//int dx = i-x;
					//int dy = j-y;
					//int dz = k-z;
					//if(dx*dx+dy*dy+dz*dz>gridWidth+radius*radius) continue;
					if(grid[i][j][k]!=null) res.addAll(grid[i][j][k]);
				}
			}			
		}
		return res;
	}
	
	/**
	 * Returns the closest neighbour. maxDist has to be given to take advantage of the speed
	 * @param c
	 * @param maxDist
	 * @return
	 */
	public int getClosestNeighbour(Coordinates c, double maxDist) {
		Set<Integer> set = getNeighbours(c, maxDist);
		int closest = -1;
		double bestDist = maxDist;
		for (int a : set) {
			double d = mol.getCoordinates(a).distanceSquared(c);
			if(d<bestDist) {
				closest = a;
				bestDist = d;
			}
		}
		return closest;
	}
	
	public void updateGrid(StereoMolecule mol) {
		for (int i = 0; i < mol.getAllAtoms(); i++) {
			int x = (int)((mol.getAtomX(i)-min.x)/gridWidth);
			int y = (int)((mol.getAtomY(i)-min.y)/gridWidth);
			int z = (int)((mol.getAtomZ(i)-min.z)/gridWidth);			
			if(grid[x][y][z]==null) grid[x][y][z] = new TreeSet<Integer>();
			if(grid[x][y][z].contains(i))
				continue;
			else {
				removeAtom(i);
				grid[x][y][z].add(i);
			}
				
			}
	}
	
	private void removeAtom(int index) {
		int l = Math.max(0, gridSize[0]);
		int m = Math.max(0, gridSize[1]);
		int n = Math.max(0, gridSize[2]);
		for(int i=0;i<l;i++) {
			for(int j=0;j<m;j++) {
				for(int k=0;k<n;n++) {
					if (grid[i][k][k].contains(index))
						grid[i][k][k].remove(index);
				}
			}
		}
		
	}
	
	public     Set<Integer> getNeighbours(Molecule mol, int atom, double maxDist) {
        return getNeighbours(mol, atom, maxDist, false);
    }
	
	 public Set<Integer> getNeighbours(Molecule mol, int atom, double maxDist, boolean enforceDist) {
	        Set<Integer> res = getNeighbours(mol.getCoordinates(atom), maxDist, enforceDist);
	        res.remove(Integer.valueOf(atom));
	        return res;
	    }
	
	public int[] getGridCoordinates(Coordinates c) {
		int[] gridCoords = new int[3];
		gridCoords[0] = (int)((c.x-min.x)/gridWidth);
		gridCoords[1] = (int)((c.y-min.y)/gridWidth);
		gridCoords[2] = (int)((c.z-min.z)/gridWidth);	
		return gridCoords;
	}
	
	public Coordinates getCartCoordinates(int[] gridCoords) {
		int gridX = gridCoords[0];
		int gridY = gridCoords[1];
		int gridZ = gridCoords[2];
		Coordinates cartCoords = new Coordinates();
		cartCoords.x = min.x + gridX*gridWidth;
		cartCoords.y = min.y + gridY*gridWidth;
		cartCoords.z = min.z + gridZ*gridWidth;
		
		return cartCoords;
	}
	
	public int[] getGridSize() {
		return gridSize;
	}

}
