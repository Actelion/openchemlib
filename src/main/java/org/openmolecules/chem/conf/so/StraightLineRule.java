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
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

import java.util.ArrayList;

public class StraightLineRule extends ConformationRule {

	/**
	 * Test code...
	 * @param args
	 */
	public static void main(String args[]) {
		final int LINE_ATOMS = 7;
		final int FIRST_LINE_ATOM = 4;
		final int LAST_LINE_ATOM = FIRST_LINE_ATOM+LINE_ATOMS-1;
		int[] atom = new int[LINE_ATOMS];
		for (int i=0; i<LINE_ATOMS; i++)
			atom[i] = FIRST_LINE_ATOM+i;
		new StraightLineRule(atom);
		StereoMolecule mol = new IDCodeParser(true).getCompactMolecule("fkA@@@DjYfYhIbnRtZjjjjjXAPbIbDUD@");
		SelfOrganizedConformer conformer = new SelfOrganizedConformer(mol);
		for (int i=0; i<LINE_ATOMS; i++) {	// straight line
			conformer.setX(FIRST_LINE_ATOM+i, i);
			conformer.setY(FIRST_LINE_ATOM+i, i);
			conformer.setZ(FIRST_LINE_ATOM+i, i);
			}
		conformer.getCoordinates(FIRST_LINE_ATOM+1).y += 1;
		conformer.getCoordinates(FIRST_LINE_ATOM+2).y -= 1;
		conformer.getCoordinates(FIRST_LINE_ATOM+4).y -= 1;
		conformer.getCoordinates(FIRST_LINE_ATOM+5).y += 1;

		System.out.println("------------------- you may copy and paste to DataWarrior ---------------------");
		System.out.println("p\tx\ty\tz\ttime");
		for (int i=FIRST_LINE_ATOM; i<=LAST_LINE_ATOM; i++)
			System.out.println("p"+i+"\t"+conformer.getX(i)+"\t"+conformer.getY(i)+"\t"+conformer.getZ(i)+"\tbefore");
		new StraightLineRule(atom).apply(conformer, 1f);
		for (int i=FIRST_LINE_ATOM; i<=LAST_LINE_ATOM; i++)
			System.out.println("p"+i+"\t"+conformer.getX(i)+"\t"+conformer.getY(i)+"\t"+conformer.getZ(i)+"\tafter");

		for (int i=0; i<LINE_ATOMS; i++) {	// straight line
			conformer.setX(FIRST_LINE_ATOM+i, i);
			conformer.setY(FIRST_LINE_ATOM+i, 0);
			conformer.setZ(FIRST_LINE_ATOM+i, 0);
			}
		conformer.getCoordinates(FIRST_LINE_ATOM+1).y -= 1;
		conformer.getCoordinates(FIRST_LINE_ATOM+2).y += 1;
		conformer.getCoordinates(FIRST_LINE_ATOM+4).y += 1;
		conformer.getCoordinates(FIRST_LINE_ATOM+5).y -= 1;

		conformer.getCoordinates(FIRST_LINE_ATOM+0).z -= 1;
		conformer.getCoordinates(FIRST_LINE_ATOM+2).z += 1;
		conformer.getCoordinates(FIRST_LINE_ATOM+4).z += 1;
		conformer.getCoordinates(FIRST_LINE_ATOM+6).z -= 1;

		System.out.println("------------------- you may copy and paste to DataWarrior ---------------------");
		System.out.println("p\tx\ty\tz\ttime");
		for (int i=FIRST_LINE_ATOM; i<=LAST_LINE_ATOM; i++)
			System.out.println("p"+i+"\t"+conformer.getX(i)+"\t"+conformer.getY(i)+"\t"+conformer.getZ(i)+"\tbefore");
		new StraightLineRule(atom).apply(conformer, 1f);
		for (int i=FIRST_LINE_ATOM; i<=LAST_LINE_ATOM; i++)
			System.out.println("p"+i+"\t"+conformer.getX(i)+"\t"+conformer.getY(i)+"\t"+conformer.getZ(i)+"\tafter");
		}

	public StraightLineRule(int[] atom) {
		super(atom);
		}

    public static void calculateRules(ArrayList<ConformationRule> ruleList, StereoMolecule mol) {
	    boolean[] atomHandled = new boolean[mol.getAllAtoms()];
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (!atomHandled[atom] && mol.getAtomPi(atom) == 2 && mol.getAtomicNo(atom) <= 8) {
				int[] atomList = getLineFragmentAtoms(atom, mol);
				for (int i=0; i<atomList.length; i++)
					atomHandled[atomList[i]] = true;

				if (atomList[0] != atomList[atomList.length-1]) // we don't add cycles
					ruleList.add(new StraightLineRule(atomList));
				}
			}
		}

    /**
     * Compiles a list atom indexes of a straight line atom strand in strand order.
     * It also sets the atomHandled flags for all strand atoms.
     * @param seedAtom
     * @param mol
     * @return
     */
	private static int[] getLineFragmentAtoms(int seedAtom, StereoMolecule mol) {
		// crawl to the end of the chain
		int firstAtom = seedAtom;
		int backAtom = mol.getConnAtom(firstAtom, 0);
		while (mol.getAllConnAtoms(firstAtom) != 1 && mol.getAtomPi(firstAtom) == 2) {	// while we are not at the chain end
			int tempAtom = backAtom;
			backAtom = firstAtom;
			firstAtom = (mol.getConnAtom(firstAtom, 0) == tempAtom) ?
					mol.getConnAtom(firstAtom, 1) : mol.getConnAtom(firstAtom, 0);

			if (firstAtom == seedAtom)  // cycle found
				break;
			}

		// count number of atoms in straight atom strand
		int rearAtom = firstAtom;	// invert direction for counting
		int headAtom = backAtom;
		int count = 2;
		while (mol.getAllConnAtoms(headAtom) != 1 && mol.getAtomPi(headAtom) == 2) {
			int tempAtom = rearAtom;
			rearAtom = headAtom;
			headAtom = (mol.getConnAtom(headAtom, 0) == tempAtom) ?
					mol.getConnAtom(headAtom, 1) : mol.getConnAtom(headAtom, 0);
			count++;

			if (headAtom == firstAtom)  // cycle found
				break;
			}

		int[] atomList = new int[count];
		atomList[0] = firstAtom;
		atomList[1] = backAtom;
		for (int i=2; i<count; i++)
			atomList[i] = (mol.getConnAtom(atomList[i-1], 0) == atomList[i-2]) ?
					mol.getConnAtom(atomList[i-1], 1) : mol.getConnAtom(atomList[i-1], 0);
//		int index = 1;
//		while (mol.getAllConnAtoms(atomList[index]) != 1 && mol.getAtomPi(atomList[index]) == 2) {
//			atomList[index+1] = (mol.getConnAtom(atomList[index], 0) == atomList[index-1]) ?
//					mol.getConnAtom(atomList[index], 1) : mol.getConnAtom(atomList[index], 0);
//			index++;
//			}

		return atomList;
		}

	@Override
	public boolean apply(Conformer conformer, double cycleFactor) {
		double[] cog = new double[3];	// center of gravity
		for (int i=0; i<mAtom.length; i++) {
			cog[0] += conformer.getX(mAtom[i]);
			cog[1] += conformer.getY(mAtom[i]);
			cog[2] += conformer.getZ(mAtom[i]);
			}
		for (int j=0; j<3; j++)
			cog[j] /= mAtom.length;

		double[][] A = new double[mAtom.length][3];
		for (int i=0; i<mAtom.length; i++) {
			A[i][0] = conformer.getX(mAtom[i]) - cog[0];
			A[i][1] = conformer.getY(mAtom[i]) - cog[1];
			A[i][2] = conformer.getZ(mAtom[i]) - cog[2];
			}

		double[][] squareMatrix = new double[3][3];
		for (int i=0; i<mAtom.length; i++)
			for (int j=0; j<3; j++)
				for (int k=0; k<3; k++)
					squareMatrix[j][k] += A[i][j] * A[i][k];

		SingularValueDecomposition svd = new SingularValueDecomposition(squareMatrix, null, null);
		double[] S = svd.getSingularValues();
		int maxIndex = 0;
		for (int i=1; i<3; i++)
			if (S[i] > S[maxIndex])
				maxIndex = i;

		double[][] U = svd.getU();
		double[] n = new double[3];	// normal vector of fitted line
		for (int i=0; i<3; i++)
			n[i] = U[i][maxIndex];

		double[] lambda = new double[mAtom.length];
		for (int i=0; i<mAtom.length; i++)
			lambda[i] = n[0]*A[i][0]+n[1]*A[i][1]+n[2]*A[i][2];

		// in early cycles check and correct order of atoms in strand
		if (cycleFactor == 1.0f) {
			boolean isInconsistent = false;
			boolean isIncreasing = (lambda[0] < lambda[1]);
			for (int i=2; i<mAtom.length; i++) {
				if (isIncreasing != (lambda[i-1] < lambda[i])) {
					isInconsistent = true;
					break;
					}
				}

			if (isInconsistent)
				return false;
/*	re-ordering atoms is problematic if we have multiple straight chains connecting in one point
 * 
 * 			if (isInconsistent) {	// re-arrange lambda values to get atoms back in order
				double s1 = 0;
				double s2 = 0;
				double min = Double.MAX_VALUE;
				double max = Double.MIN_VALUE;
				for (int i=0; i<mAtom.length; i++) {
					s1 += lambda[i] * (mAtom.length - 1 - i);
					s2 += lambda[i] * i;
					min = Math.min(min, lambda[i]);
					max = Math.max(max, lambda[i]);
					}
				isIncreasing = (s2 > s1);
				double delta = isIncreasing ? max-min : min-max;
				lambda[0] = isIncreasing ? min : max;
				for (int i=1; i<mAtom.length; i++)
					lambda[i] = lambda[i-1] + delta;
				}	*/
			}

		// for a point P on the fitted line is: P = COG + lamda * N
		for (int i=0; i<mAtom.length; i++) {
				// calculate lambda that gives the closest point to current atom location
			conformer.getCoordinates(mAtom[i]).add(cycleFactor*(lambda[i]*n[0]-A[i][0]),
												  cycleFactor*(lambda[i]*n[1]-A[i][1]),
												  cycleFactor*(lambda[i]*n[2]-A[i][2]));
			}

		return true;
		}

	@Override
	public double addStrain(Conformer conformer, double[] atomStrain) {
		double[] cog = new double[3];	// center of gravity
		for (int i=0; i<mAtom.length; i++) {
			cog[0] += conformer.getX(mAtom[i]);
			cog[1] += conformer.getY(mAtom[i]);
			cog[2] += conformer.getZ(mAtom[i]);
			}
		for (int j=0; j<3; j++)
			cog[j] /= mAtom.length;

		double[][] A = new double[mAtom.length][3];
		for (int i=0; i<mAtom.length; i++) {
			A[i][0] = conformer.getX(mAtom[i]) - cog[0];
			A[i][1] = conformer.getY(mAtom[i]) - cog[1];
			A[i][2] = conformer.getZ(mAtom[i]) - cog[2];
			}

		double[][] squareMatrix = new double[3][3];
		for (int i=0; i<mAtom.length; i++)
			for (int j=0; j<3; j++)
				for (int k=0; k<3; k++)
					squareMatrix[j][k] += A[i][j] * A[i][k];

		SingularValueDecomposition svd = new SingularValueDecomposition(squareMatrix, null, null);
		double[] S = svd.getSingularValues();
		int maxIndex = 0;
		for (int i=1; i<3; i++)
			if (S[i] > S[maxIndex])
				maxIndex = i;

		double[][] U = svd.getU();
		double[] n = new double[3];	// normal vector of fitted line
		for (int i=0; i<3; i++)
			n[i] = U[i][maxIndex];

		double totalStrain = 0;
		for (int i=0; i<mAtom.length; i++) {
				// calculate lambda that gives the closest point to current atom location
			double lambda = n[0]*A[i][0]+n[1]*A[i][1]+n[2]*A[i][2];
			double dx = lambda*n[0]-A[i][0];
			double dy = lambda*n[1]-A[i][1];
			double dz = lambda*n[2]-A[i][2];
			double strain = 10.0 * dx*dx+dy*dy+dz*dz;
			if (atomStrain != null)
				atomStrain[mAtom[i]] += strain;
			totalStrain += strain;
			}

		return totalStrain;
		}

	@Override
	public int getRuleType() {
		return RULE_TYPE_LINE;
		}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("line rule:");
		super.addAtomList(sb);
		return sb.toString();
		}
	}
