/*
 * Copyright 2013-2020 Thomas Sander, openmolecules.org
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.so;

import com.actelion.research.calc.SingularValueDecomposition;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

public abstract class ConformationRule {
	public static final int RULE_TYPE_DISTANCE = 0;
	public static final int RULE_TYPE_PLANE = 1;
	public static final int RULE_TYPE_LINE = 2;
	public static final int RULE_TYPE_TORSION = 3;
	public static final int RULE_TYPE_STEREO = 4;
	public static final int RULE_TYPE_BINAP = 5;

	public static final String[] RULE_NAME = { "distance", "plane", "line", "torsion", "stereo", "binap" };

	protected int[] mAtom;
	protected boolean mIsEnabled;

	public ConformationRule(int[] atom) {
		mAtom = atom;
		mIsEnabled = true;
		}

	public boolean isEnabled() {
		return mIsEnabled;
		}

	public void setEnabled(boolean b) {
		mIsEnabled = b;
		}

	public abstract boolean apply(Conformer conformer, double cycleFactor);
	public abstract double addStrain(Conformer conformer, double[] atomStrain);
	public abstract int getRuleType();
	public abstract String toString();

	protected int[] getAtomList() {
		return mAtom;
		}

	protected void addAtomList(StringBuilder sb) {
		if (mAtom == null) {
		    sb.append(" atoms:<null>");
			}
		else {
		    sb.append(" atoms:"+mAtom[0]);
	        for (int i=1; i<mAtom.length; i++)
	            sb.append(","+mAtom[i]);
			}
		}

	/**
	 * Move one atom and all non-ring-bond substituents,
	 * provided attachment atom is not part of the notList.
	 * @param conformer
	 * @param atom
	 * @param notList atom list that are not moved
	 * @param dx
	 * @param dy
	 * @param dz
	 */
	protected void moveGroup(Conformer conformer, int atom, int[] notList, double dx, double dy, double dz) {
//System.out.print("move atom: "+atom);
		conformer.getCoordinates(atom).add(dx, dy, dz);
	    StereoMolecule mol = conformer.getMolecule();
        for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
        	int connAtom = mol.getConnAtom(atom, i);
        	boolean found = false;
        	if (notList != null) {
	            for (int j=0; j<notList.length; j++) {
	            	if (connAtom == notList[j]) {
	            		found = true;
	            		break;
	            		}
	            	}
        		}
            if (!found) {
            	if (mol.getAllConnAtoms(connAtom) == 1) {
		            conformer.getCoordinates(connAtom).add(dx, dy, dz);
//System.out.print(" conn"+i+":"+connAtom);
            		}
            	else if (mol.isRingBond(mol.getConnBond(atom, i))) {
		            conformer.getCoordinates(connAtom).add(dx/3, dy/3, dz/3);
//System.out.print(" rconn"+i+":"+connAtom);
            		}
            	else {
            		moveSubstituent(conformer, atom, connAtom, dx, dy, dz);
//System.out.print(" subst"+i+":"+connAtom);
            		}
            	}
        	}
/*
if (notList != null) {
 System.out.print(" not:");
 for (int a:notList)
  System.out.print(" "+a);
 }
System.out.println();
*/
		}

	/**
	 * Moves all atoms of the substituent connected to rootAtom starting with firstAtom.
	 * @param conformer
	 * @param rootAtom
	 * @param firstAtom
	 * @param dx
	 * @param dy
	 * @param dz
	 */
	protected void moveSubstituent(Conformer conformer, int rootAtom, int firstAtom, double dx, double dy, double dz) {
		conformer.getCoordinates(firstAtom).add(dx, dy, dz);

		StereoMolecule mol = conformer.getMolecule();

	    boolean[] isMemberAtom = new boolean[mol.getAllAtoms()];
	    int[] graphAtom = new int[mol.getAllAtoms()];

		graphAtom[0] = firstAtom;
		isMemberAtom[rootAtom] = true;
		isMemberAtom[firstAtom] = true;
		int current = 0;
		int highest = 0;
	 	while (current <= highest) {
			for (int i=0; i<mol.getAllConnAtoms(graphAtom[current]); i++) {
				int candidate = mol.getConnAtom(graphAtom[current], i);
				if (!isMemberAtom[candidate]) {
					isMemberAtom[candidate] = true;
					graphAtom[++highest] = candidate;
					conformer.getCoordinates(candidate).add(dx, dy, dz);
					}
				}
			current++;
			}
		}

	protected void moveAtomWithUnboundedNeighbors(Conformer conformer, int atom, double dx, double dy, double dz) {
		conformer.getCoordinates(atom).add(dx, dy, dz);
		StereoMolecule mol = conformer.getMolecule();
	    for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
	    	int connAtom = mol.getConnAtom(atom, i);
        	if (mol.getAllConnAtoms(connAtom) == 1)
		        conformer.getCoordinates(connAtom).add(dx, dy, dz);
	    	}
	    }

	/**
	 * Rotate one atom and all connected atoms except those in the notList and beyond.
	 * @param conformer
	 * @param atom
	 * @param notList atom list that are not moved
	 * @param p center of rotation
	 * @param n normal vector
	 * @param theta angle
	 */
	protected static void rotateGroup(Conformer conformer, int atom, int[] notList, Coordinates p, Coordinates n, double theta) {
		StereoMolecule mol = conformer.getMolecule();

		boolean[] avoidAtom = new boolean[mol.getAllAtoms()];
		int[] graphAtom = new int[mol.getAllAtoms()];

		graphAtom[0] = atom;
		avoidAtom[atom] = true;
		for (int notAtom:notList)
			avoidAtom[notAtom] = true;

		rotateAtom(conformer, atom, p, n, theta);

		int current = 0;
		int highest = 0;
		while (current <= highest) {
			for (int i=0; i<mol.getAllConnAtoms(graphAtom[current]); i++) {
				int candidate = mol.getConnAtom(graphAtom[current], i);
				if (!avoidAtom[candidate]) {
					avoidAtom[candidate] = true;
					graphAtom[++highest] = candidate;
					rotateAtom(conformer, candidate, p, n, theta);
					}
				}
			current++;
			}
		}

	protected static void rotateAtom(Conformer conformer, int atom, Coordinates p, Coordinates unit, double theta) {
		double x = unit.x;
		double y = unit.y;
		double z = unit.z;
		double c = Math.cos(theta);
		double s = Math.sin(theta);
		double t = 1-c;
		double mx = conformer.getX(atom) - p.x;
		double my = conformer.getY(atom) - p.y;
		double mz = conformer.getZ(atom) - p.z;
		conformer.setX(atom, p.x + (t*x*x+c)*mx + (t*x*y+s*z)*my + (t*x*z-s*y)*mz);
		conformer.setY(atom, p.y + (t*x*y-s*z)*mx + (t*y*y+c)*my + (t*y*z+s*x)*mz);
		conformer.setZ(atom, p.z + (t*x*z+s*y)*mx + (t*z*y-s*x)*my + (t*z*z+c)*mz);
		}

	/**
	 * @param conformer
	 * @param atom
	 * @param cog receives the center of gravity
	 * @param n receives normal vector of plane nearest to all atoms
	 * @param coords receives original atom coordinates minus the center of gravity [atom count][3]
	 */
	public static void calculateNearestPlane(Conformer conformer, int[] atom, Coordinates cog, Coordinates n, double[][] coords) {
		for (int i=0; i<atom.length; i++)
			cog.add(conformer.getCoordinates(atom[i]));
		cog.scale(1.0 / atom.length);

		for (int i=0; i<atom.length; i++) {
			coords[i][0] = conformer.getX(atom[i]) - cog.x;
			coords[i][1] = conformer.getY(atom[i]) - cog.y;
			coords[i][2] = conformer.getZ(atom[i]) - cog.z;
		}

		double[][] squareMatrix = new double[3][3];
		for (int i=0; i<atom.length; i++)
			for (int j=0; j<3; j++)
				for (int k=0; k<3; k++)
					squareMatrix[j][k] += coords[i][j] * coords[i][k];

		SingularValueDecomposition svd = new SingularValueDecomposition(squareMatrix, null, null);
		double[] S = svd.getSingularValues();
		int minIndex = 0;
		for (int i=1; i<3; i++)
			if (S[i] < S[minIndex])
				minIndex = i;

		double[][] U = svd.getU();
		n.x = U[0][minIndex];
		n.y = U[1][minIndex];
		n.z = U[2][minIndex];
		}
	}
