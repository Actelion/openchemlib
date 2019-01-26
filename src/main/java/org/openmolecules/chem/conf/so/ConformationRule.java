/*
 * @(#)ConformationRule.java
 *
 * Copyright 2013 openmolecules.org, Inc. All Rights Reserved.
 * 
 * NOTICE: All information contained herein is, and remains the property
 * of openmolecules.org.  The intellectual and technical concepts contained
 * herein are proprietary to openmolecules.org.
 * Actelion Pharmaceuticals Ltd. is granted a non-exclusive, non-transferable
 * and timely unlimited usage license.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.so;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

public abstract class ConformationRule {
	public static final int RULE_TYPE_DISTANCE = 0;
	public static final int RULE_TYPE_PLANE = 1;
	public static final int RULE_TYPE_LINE = 2;
	public static final int RULE_TYPE_TORSION = 3;
	public static final int RULE_TYPE_STEREO = 4;

	public static final String[] RULE_NAME = { "distance", "plane", "line", "torsion", "stereo" };

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
	}
