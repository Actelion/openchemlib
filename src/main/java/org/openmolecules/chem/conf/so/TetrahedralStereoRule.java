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

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.MolfileCreator;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;

import java.util.ArrayList;

public class TetrahedralStereoRule extends ConformationRule {
	private int[] mRotatableAtom,mStaticNeighbour,mRotatableNeighbour;

	/**
	 * Expects neighbourAtom[5] as list of 3 or 4 neighbour atoms of the stereo center.
	 * The list must be sorted by atom index. In case of 3 neighbours neighbourAtom[3] is -1.
	 * neighbourAtom[4] contains the stereo center. If desired parity==cAtomParity1, then
	 * neighbourAtom[0] && neighbourAtom[1] have swapped positions.
	 * @param mol
	 * @param neighbourAtom
	 * @param neighbourBond
	 */
	public TetrahedralStereoRule(StereoMolecule mol, int[] neighbourAtom, int[] neighbourBond) {
		super(neighbourAtom);
		calculateRotatableAtoms(mol, neighbourAtom, neighbourBond);
		}

    public static void calculateRules(ArrayList<ConformationRule> ruleList, StereoMolecule mol) {
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mol.getAllConnAtoms(atom) >= 3) {
				int parity = mol.getAtomParity(atom);
				if ((parity == Molecule.cAtomParity1 || parity == Molecule.cAtomParity2)
				 && !mol.isCentralAlleneAtom(atom)) {
					int[] atomList = new int[5];
					int[] bondList = new int[4];
					for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
						int connAtom = mol.getConnAtom(atom, i);
						int connBond = mol.getConnBond(atom, i);
						int index = 0;
						while (index < i && connAtom > atomList[index])
							index++;
						for (int j=i-1; j>=index; j--) {
							atomList[j + 1] = atomList[j];
							bondList[j + 1] = bondList[j];
							}
						atomList[index] = connAtom;
						bondList[index] = connBond;
						}

					if (mol.getAllConnAtoms(atom) == 3) {
						atomList[3] = -1;
						bondList[3] = -1;
						}

					atomList[4] = atom;

					if (parity == Molecule.cAtomParity1) {
						int temp = atomList[1];
						atomList[1] = atomList[0];
						atomList[0] = temp;
						temp = bondList[1];
						bondList[1] = bondList[0];
						bondList[0] = temp;
						}

					ruleList.add(new TetrahedralStereoRule(mol, atomList, bondList));
					}
				}
			}
	    }

	@Override
	public boolean apply(Conformer conformer, double cycleFactor) {
		double[] n = getPlaneVector(conformer);
		// We invert if 3rd substituent is on wrong side (3 substituents)
		// or if 3rd and 4th substituent are on wrong side (4 substituents)
		boolean invert = (mAtom[3] == -1 && !isOnSameSide(conformer, n, mAtom[4], mAtom[2]))
					  || (mAtom[3] != -1 && !isOnSameSide(conformer, n, mAtom[4], mAtom[2])
										 && isOnSameSide(conformer, n, mAtom[4], mAtom[3]));
		if (!invert)
			return false;

		Coordinates rotationAxis = calculateRotationAxis(conformer);
		for (int atom:mRotatableAtom)
			rotateAtom(conformer, atom, conformer.getCoordinates(mAtom[4]), rotationAxis, Math.PI);

		return true;
		}

	@Override
	public double addStrain(Conformer conformer, double[] atomStrain) {
		double totalStrain = 0;
		double[] n = getPlaneVector(conformer);
		// We add a penalty for 3rd and 4th (if existing) substituent that is on the wrong side
		int wrongCount = 0;
		if (mAtom[3] == -1) {
			if (!isOnSameSide(conformer, n, mAtom[4], mAtom[2]))
				wrongCount++;
			}
		else {
			if (!isOnSameSide(conformer, n, mAtom[4], mAtom[2]))
				wrongCount++;
			if (isOnSameSide(conformer, n, mAtom[4], mAtom[3]))
				wrongCount++;
			}
		if (wrongCount != 0) {
			for (int i=0; i<mAtom.length; i++) {
				if (mAtom[i] != -1) {
					double panalty = 0.68 * wrongCount; // arbitrary value: likelyhood factor 10 per two atoms;
					if (atomStrain != null)
						atomStrain[mAtom[i]] += panalty;
					totalStrain += panalty;
					}
				}
			}
		return totalStrain;
		}

	/**
	 * Returns whether atom1 is on the side of the plane defined by atom0 and n
	 * @param conformer
	 * @param n vector of the plane a0,a1,a2
	 * @param atom0	atom on plane
	 * @param atom1 atom not on plane
	 * @return
	 */
	private boolean isOnSameSide(Conformer conformer, double[] n, int atom0, int atom1) {
		// calculate the three vectors leading from atom[0] to the other three atoms
		double[] v = new double[3];
		v[0] = conformer.getX(atom1) - conformer.getX(atom0);
		v[1] = conformer.getY(atom1) - conformer.getY(atom0);
		v[2] = conformer.getZ(atom1) - conformer.getZ(atom0);

		// calculate cos(angle) of a0->a1 to normal vector
		return v[0]*n[0]+v[1]*n[1]+v[2]*n[2] > 0;
		}

	/**
	 * Calculates the vector of plane (a4->a1, a4->a2)
	 * @param conformer
	 * @return
	 */
	private double[] getPlaneVector(Conformer conformer) {
		double[][] coords = new double[2][3];
		for (int i=0; i<2; i++) {
			coords[i][0] = conformer.getX(mAtom[i]) - conformer.getX(mAtom[4]);
			coords[i][1] = conformer.getY(mAtom[i]) - conformer.getY(mAtom[4]);
			coords[i][2] = conformer.getZ(mAtom[i]) - conformer.getZ(mAtom[4]);
			}

		// calculate the vector of the plane (vector product of coords[0] and coords[1])
		double[] v = new double[3];
		v[0] = coords[0][1]*coords[1][2]-coords[0][2]*coords[1][1];
		v[1] = coords[0][2]*coords[1][0]-coords[0][0]*coords[1][2];
		v[2] = coords[0][0]*coords[1][1]-coords[0][1]*coords[1][0];
		return v;
		}

	private void calculateRotatableAtoms(StereoMolecule mol, int[] neighbourAtom, int[] neighbourBond) {
		int stereoCenter = neighbourAtom[4];
		int neighbours = (neighbourAtom[3] == -1) ? 3 : 4;

		int[] fragmentNo = new int[mol.getAllAtoms()];
		boolean[] neglectBond = new boolean[mol.getAllBonds()];
		for (int i=0; i<neighbours; i++)
			neglectBond[neighbourBond[i]] = true;
		int fragmentCount = mol.getFragmentNumbers(fragmentNo, neglectBond, false);

		int[] fragmentSize = new int[fragmentCount];
		for (int f:fragmentNo)
			fragmentSize[f]++;

		int neighbourFragmentCount = 0;
		boolean[] isNeighbourFragment = new boolean[fragmentCount];
		for (int i=0; i<neighbours; i++) {
			int f = fragmentNo[neighbourAtom[i]];
			if (!isNeighbourFragment[f]) {
				isNeighbourFragment[f] = true;
				neighbourFragmentCount++;
				}
			}

		int[] neighbourFragmentNo = new int[neighbourFragmentCount];
		int[] neighbourFragmentSize = new int[neighbourFragmentCount];
		int[] neighbourFragmentBondCount = new int[neighbourFragmentCount];
		int index = 0;
		for (int f=0; f<fragmentCount; f++) {
			if (isNeighbourFragment[f]) {
				neighbourFragmentNo[index] = f;
				neighbourFragmentSize[index] = fragmentSize[f];
				for (int i=0; i<neighbours; i++)
					if (f == fragmentNo[neighbourAtom[i]])
						neighbourFragmentBondCount[index]++;
				index++;
			}
		}

		boolean[] isRotatableAtom = new boolean[mol.getAllAtoms()];
		int rotatableAtomCount = 0;
		int staticNeighbourCount = 2;   // default

		// Now we know about every fragment that touches the stereo center:
		// - number of bonds connecting to stereo center
		// - number of fragment atoms
		// - the fragment number
		if (neighbourFragmentCount == neighbours) { // no rings: we need to rotate neighbours-2 substituents (== fragments)
			int lastSmallestFragmentIndex = -1;
			for (int i=0; i<neighbours-2; i++) {
				int smallestFragmentIndex = -1;
				int smallestFragmentSize = Integer.MAX_VALUE;
				for (int j=0; j<neighbourFragmentCount; j++) {
					if (j != lastSmallestFragmentIndex && neighbourFragmentSize[j] < smallestFragmentSize) {
						smallestFragmentIndex = j;
						smallestFragmentSize = neighbourFragmentSize[j];
					}
				}
				for (int a=0; a<mol.getAllAtoms(); a++) {
					if (fragmentNo[a] == neighbourFragmentNo[smallestFragmentIndex]) {
						isRotatableAtom[a] = true;
						rotatableAtomCount++;
						}
					}
				lastSmallestFragmentIndex = smallestFragmentIndex;
			}
		}
		else {  // If we have one ring (two bonds share the same fragment), we can rotate the ring (if four neighbours) or the rest
			for (int i=0; i<neighbourFragmentNo.length; i++) {
				if (neighbourFragmentBondCount[i] == 1 && neighbours == 3) {
					boolean isSmallerFragment = 2*neighbourFragmentSize[i] < mol.getAllAtoms();
					if (!isSmallerFragment)
						staticNeighbourCount = 1;
					for (int a=0; a<mol.getAllAtoms(); a++) {
						if (a != stereoCenter) {
							if ((fragmentNo[a] == neighbourFragmentNo[i]) == isSmallerFragment) {
								isRotatableAtom[a] = true;
								rotatableAtomCount++;
							}
						}
					}
					break;
				}
				if (neighbourFragmentBondCount[i] == 2 && neighbours == 4) {
					boolean isSmallerFragment = 2*neighbourFragmentSize[i] < mol.getAllAtoms();
					for (int a=0; a<mol.getAllAtoms(); a++) {
						if (a != stereoCenter) {
							if ((fragmentNo[a] == neighbourFragmentNo[i]) == isSmallerFragment) {
								isRotatableAtom[a] = true;
								rotatableAtomCount++;
							}
						}
					}
					break;
				}
			}
		}

		if (rotatableAtomCount == 0) {
			for (int i=0; i<neighbourFragmentNo.length; i++) {
				if (neighbourFragmentBondCount[i] == 1) {
					for (int a=0; a<mol.getAllAtoms(); a++) {
						if (fragmentNo[a] == neighbourFragmentNo[i]) {
							isRotatableAtom[a] = true;
							rotatableAtomCount++;
						}
					}
					break;
				}
			}
			// TODO We have a >=3 neighbour-bond-fragment and to proportionally rotate the closer part of it and not just one atom
			int remainingNeighboursToRotate = neighbours - 2 - (rotatableAtomCount == 0 ? 0 : 1);
			for (int i=0; i<remainingNeighboursToRotate; i++) {
				for (int j=0; j<neighbours; j++) {
					if (!isRotatableAtom[neighbourAtom[j]]) {
						isRotatableAtom[neighbourAtom[j]] = true;
						rotatableAtomCount++;
						break;
					}
				}
			}
		}

		index = 0;
		mRotatableAtom = new int[rotatableAtomCount];
		for (int i=0; i<isRotatableAtom.length; i++)
			if (isRotatableAtom[i])
				mRotatableAtom[index++] = i;

		int staticIndex = 0;
		int rotatableIndex = 0;
		mStaticNeighbour = new int[staticNeighbourCount];
		mRotatableNeighbour = new int[neighbours - staticNeighbourCount];
		for (int i=0; i<neighbours; i++)
			if (isRotatableAtom[neighbourAtom[i]])
{
if (rotatableIndex>=mRotatableNeighbour.length) {
 System.out.println("##### Out of bounds exception:"+new MolfileCreator(mol).getMolfile());
 System.out.print("neighbours:"+neighbours+" ");
 System.out.print("neighbourAtom:"); for (int k:neighbourAtom) System.out.print(" " + k); System.out.println();
 System.out.print("neighbourBond:"); for (int k : neighbourBond) System.out.print(" " + k); System.out.println();
 System.out.print("fragmentNo:"); for (int k : fragmentNo) System.out.print(" " + k); System.out.println();
 System.out.print("isNeighbourFragment:"); for (boolean k : isNeighbourFragment) System.out.print(" " + k); System.out.println();
 System.out.print("isRotatableAtom:"); for (boolean b : isRotatableAtom) System.out.print(" " + b); System.out.println();
}
				mRotatableNeighbour[rotatableIndex++] = neighbourAtom[i];
}
			else
				mStaticNeighbour[staticIndex++] = neighbourAtom[i];
		}

	private Coordinates calculateRotationAxis(Conformer conformer) {
		Coordinates axis = new Coordinates();

		int[] neighbour = (mStaticNeighbour.length == 2) ? mStaticNeighbour : mRotatableNeighbour;
		for (int atom:neighbour)
			axis.add(conformer.getCoordinates(atom));

		axis.scale(0.5);
		axis.sub(conformer.getCoordinates(mAtom[4])).unit();

		return axis;
		}

	@Override
	public int getRuleType() {
		return RULE_TYPE_STEREO;
		}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("stereo rule:");
		super.addAtomList(sb);
		sb.append(" rotatable:");
		for (int ra:mRotatableAtom)
			sb.append(ra+" ");
		return sb.toString();
		}
	}
