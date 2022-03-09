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
 * @author Thomas Sander
 */

package com.actelion.research.chem.coords;

import com.actelion.research.chem.StereoMolecule;

import java.util.ArrayList;
import java.util.Arrays;

public class InventorFragment {
	private static final double cCollisionLimitBondRotation = 0.8;
	private static final double cCollisionLimitAtomMovement = 0.5;

	private final int CIRCULAR_BINS = 36;

	protected int[] mGlobalAtom;
	protected int[] mGlobalBond;
	protected int[] mGlobalToLocalAtom;
	protected int[] mPriority;
	protected double[] mAtomX;
	protected double[] mAtomY;
	protected boolean mKeepMarkedAtoms;

	private StereoMolecule mMol;
	private boolean	mMinMaxAvail;
	private double mMinX;
	private double mMinY;
	private double mMaxX;
	private double mMaxY;
	private double mCollisionPanalty;
	private int[][] mFlipList;
	private int[] mSortedAtom;

	protected InventorFragment(StereoMolecule mol, int atoms, boolean keepMarkedAtoms) {
		mMol = mol;
		mKeepMarkedAtoms = keepMarkedAtoms;
		mGlobalAtom = new int[atoms];
		mPriority = new int[atoms];
		mAtomX = new double[atoms];
		mAtomY = new double[atoms];
	}

	protected InventorFragment(InventorFragment f) {
		mMol = f.mMol;
		mKeepMarkedAtoms = f.mKeepMarkedAtoms;
		mGlobalAtom = new int[f.size()];
		mPriority = new int[f.size()];
		mAtomX = new double[f.size()];
		mAtomY = new double[f.size()];
		for (int i=0; i<f.size(); i++) {
			mGlobalAtom[i]	 = f.mGlobalAtom[i];
			mPriority[i] = f.mPriority[i];
			mAtomX[i]	 = f.mAtomX[i];
			mAtomY[i]	 = f.mAtomY[i];
		}
		if (f.mGlobalBond != null) {
			mGlobalBond = new int[f.mGlobalBond.length];
			for (int i = 0; i<f.mGlobalBond.length; i++)
				mGlobalBond[i] = f.mGlobalBond[i];
		}
		if (f.mGlobalToLocalAtom != null) {
			mGlobalToLocalAtom = new int[f.mGlobalToLocalAtom.length];
			for (int i = 0; i<f.mGlobalToLocalAtom.length; i++)
				mGlobalToLocalAtom[i] = f.mGlobalToLocalAtom[i];
		}
	}

	protected int size() {
		return mGlobalAtom.length;
	}

	protected boolean equals(InventorFragment f) {
		if (f.size() != size())
			return false;

		int[] sorted = getSortedAtoms();
		int[] sortedF = f.getSortedAtoms();
		for (int i=0; i<sorted.length; i++)
			if (sorted[i] != sortedF[i])
				return false;

		return true;
		}

	private int[] getSortedAtoms() {
		if (mSortedAtom == null) {
			mSortedAtom = mGlobalAtom.clone();
			Arrays.sort(mSortedAtom);
			}

		return mSortedAtom;
		}

	protected double getAtomX(int index) {
		return mAtomX[index];
	}

	protected double getAtomY(int index) {
		return mAtomY[index];
	}

	protected double getWidth() {
		calculateMinMax();
		return mMaxX - mMinX + 1.0;	// add half a bond length on every side
	}

	protected double getHeight() {
		calculateMinMax();
		return mMaxY - mMinY + 1.0;	// add half a bond length on every side
	}

	protected boolean isMember(int globalAtom) {
		for (int i = 0; i<mGlobalAtom.length; i++)
			if (globalAtom == mGlobalAtom[i])
				return true;

		return false;
	}

	protected int getGlobalAtom(int localAtom) {
		return mGlobalAtom[localAtom];
	}

	protected int getLocalAtom(int globalAtom) {
		for (int i = 0; i<mGlobalAtom.length; i++)
			if (globalAtom == mGlobalAtom[i])
				return i;

		return -1;
	}

	protected void translate(double dx, double dy) {
		for (int i = 0; i<mGlobalAtom.length; i++) {
			mAtomX[i] += dx;
			mAtomY[i] += dy;
		}
	}

	protected void rotate(double x, double y, double angleDif) {
		for (int i = 0; i< mGlobalAtom.length; i++) {
			double distance = Math.sqrt((mAtomX[i] - x) * (mAtomX[i] - x)
					+ (mAtomY[i] - y) * (mAtomY[i] - y));
			double angle = InventorAngle.getAngle(x, y, mAtomX[i], mAtomY[i]) + angleDif;
			mAtomX[i] = x + distance * Math.sin(angle);
			mAtomY[i] = y + distance * Math.cos(angle);
		}
	}

	protected void flip(double x, double y, double mirrorAngle) {
		for (int i = 0; i< mGlobalAtom.length; i++) {
			double distance = Math.sqrt((mAtomX[i] - x) * (mAtomX[i] - x)
					+ (mAtomY[i] - y) * (mAtomY[i] - y));
			double angle = 2 * mirrorAngle - InventorAngle.getAngle(x, y, mAtomX[i], mAtomY[i]);
			mAtomX[i] = x + distance * Math.sin(angle);
			mAtomY[i] = y + distance * Math.cos(angle);
		}
	}

	protected void flipOneSide(int bond) {
		// The fliplist contains for every bond atoms:
		// [0]->the bond atom that lies on the larger side of the bond
		// [1]->the bond atom on the smaller side of the bond
		// [2...n]->all other atoms on the smaller side of the bond.
		//		  These are the ones getting flipped on the mirror
		//		  line defined by the bond.
		if (mFlipList == null)
			mFlipList = new int[mMol.getAllBonds()][];

		if (mFlipList[bond] == null) {
			int[] graphAtom = new int[mGlobalAtom.length];
			boolean[] isOnSide = new boolean[mMol.getAllAtoms()];
			int atom1 = mMol.getBondAtom(0, bond);
			int atom2 = mMol.getBondAtom(1, bond);
			graphAtom[0] = atom1;
			isOnSide[atom1] = true;
			int current = 0;
			int highest = 0;
			while (current <= highest) {
				for (int i=0; i<mMol.getAllConnAtoms(graphAtom[current]); i++) {
					int candidate = mMol.getConnAtom(graphAtom[current], i);

					if (!isOnSide[candidate] && candidate != atom2) {
						graphAtom[++highest] = candidate;
						isOnSide[candidate] = true;
					}
				}
				if (current == highest)
					break;
				current++;
			}

			// default is to flip the smaller side
			boolean flipOtherSide = (highest+1 > mGlobalAtom.length/2);

			// if we retain core atoms and the smaller side contains core atoms, then flip the larger side
			if (mKeepMarkedAtoms) {
				boolean coreOnSide = false;
				boolean coreOffSide = false;
				for (int i = 0; i< mGlobalAtom.length; i++) {
					int atom = mGlobalAtom[i];
					if (mMol.isMarkedAtom(atom) && atom != atom1 && atom != atom2) {
						if (isOnSide[mGlobalAtom[i]])
							coreOnSide = true;
						else
							coreOffSide = true;
					}
				}
				if (coreOnSide != coreOffSide)
					flipOtherSide = coreOnSide;
			}

			int count = 2;
			mFlipList[bond] = new int[flipOtherSide ? mGlobalAtom.length-highest : highest+2];
			for (int i = 0; i< mGlobalAtom.length; i++) {
				if (mGlobalAtom[i] == atom1)
					mFlipList[bond][flipOtherSide ? 0 : 1] = i;
				else if (mGlobalAtom[i] == atom2)
					mFlipList[bond][flipOtherSide ? 1 : 0] = i;
				else if (flipOtherSide ^ isOnSide[mGlobalAtom[i]])
					mFlipList[bond][count++] = i;
			}
		}

		double x = mAtomX[mFlipList[bond][0]];
		double y = mAtomY[mFlipList[bond][0]];
		double mirrorAngle = InventorAngle.getAngle(x, y, mAtomX[mFlipList[bond][1]],
				mAtomY[mFlipList[bond][1]]);

		for (int i=2; i<mFlipList[bond].length; i++) {
			int index = mFlipList[bond][i];
			double distance = Math.sqrt((mAtomX[index] - x) * (mAtomX[index] - x)
					+ (mAtomY[index] - y) * (mAtomY[index] - y));
			double angle = 2 * mirrorAngle - InventorAngle.getAngle(x, y, mAtomX[index], mAtomY[index]);
			mAtomX[index] = x + distance * Math.sin(angle);
			mAtomY[index] = y + distance * Math.cos(angle);
		}
	}

	protected void arrangeWith(InventorFragment f) {
		double maxGain = 0.0;
		int maxCorner = 0;
		for (int corner=0; corner<4; corner++) {
			double gain = getCornerDistance(corner) + f.getCornerDistance((corner>=2) ? corner-2 : corner+2);
			if (maxGain < gain) {
				maxGain = gain;
				maxCorner = corner;
			}
		}

		double sumHeight = getHeight() + f.getHeight();
		double sumWidth = 0.75 * (getWidth() + f.getWidth());
		double maxHeight = Math.max(getHeight(), f.getHeight());
		double maxWidth = 0.75 * Math.max(getWidth(), f.getWidth());

		double bestCornerSize = Math.sqrt((sumHeight - maxGain) * (sumHeight - maxGain)
				+ (sumWidth - 0.75 * maxGain) * (sumWidth - 0.75 * maxGain));
		double toppedSize = Math.max(maxWidth, sumHeight);
		double besideSize = Math.max(maxHeight, sumWidth);

		if (bestCornerSize < toppedSize && bestCornerSize < besideSize) {
			switch(maxCorner) {
				case 0:
					f.translate(mMaxX - f.mMinX - maxGain + 1.0, mMinY - f.mMaxY + maxGain - 1.0);
					break;
				case 1:
					f.translate(mMaxX - f.mMinX - maxGain + 1.0, mMaxY - f.mMinY - maxGain + 1.0);
					break;
				case 2:
					f.translate(mMinX - f.mMaxX + maxGain - 1.0, mMaxY - f.mMinY - maxGain + 1.0);
					break;
				case 3:
					f.translate(mMinX - f.mMaxX + maxGain - 1.0, mMinY - f.mMaxY + maxGain - 1.0);
					break;
			}
		}
		else if (besideSize < toppedSize) {
			f.translate(mMaxX - f.mMinX + 1.0, (mMaxY + mMinY - f.mMaxY - f.mMinY) / 2);
		}
		else {
			f.translate((mMaxX + mMinX - f.mMaxX - f.mMinX) / 2, mMaxY - f.mMinY + 1.0);
		}
	}

	private void calculateMinMax() {
		if (mMinMaxAvail)
			return;

		mMinX = mAtomX[0];
		mMaxX = mAtomX[0];
		mMinY = mAtomY[0];
		mMaxY = mAtomY[0];
		for (int i = 0; i< mGlobalAtom.length; i++) {
			double surplus = getAtomSurplus(i);

			if (mMinX > mAtomX[i] - surplus)
				mMinX = mAtomX[i] - surplus;
			if (mMaxX < mAtomX[i] + surplus)
				mMaxX = mAtomX[i] + surplus;
			if (mMinY > mAtomY[i] - surplus)
				mMinY = mAtomY[i] - surplus;
			if (mMaxY < mAtomY[i] + surplus)
				mMaxY = mAtomY[i] + surplus;
		}

		mMinMaxAvail = true;
	}

	private double getCornerDistance(int corner) {
		double minDistance = 9999.0;
		for (int atom = 0; atom< mGlobalAtom.length; atom++) {
			double surplus = getAtomSurplus(atom);
			double d = 0.0;
			switch (corner) {
				case 0:	// top right
					d = mMaxX - 0.5 * (mMaxX + mMinY + mAtomX[atom] - mAtomY[atom]);
					break;
				case 1:	// bottom right
					d = mMaxX - 0.5 * (mMaxX - mMaxY + mAtomX[atom] + mAtomY[atom]);
					break;
				case 2:	// bottom left
					d = 0.5 * (mMinX + mMaxY + mAtomX[atom] - mAtomY[atom]) - mMinX;
					break;
				case 3:	// top left
					d = 0.5 * (mMinX - mMinY + mAtomX[atom] + mAtomY[atom]) - mMinX;
					break;
			}

			if (minDistance > d - surplus)
				minDistance = d - surplus;
		}

		return minDistance;
	}

	private double getAtomSurplus(int atom) {
		return (mMol.getAtomQueryFeatures(mGlobalAtom[atom]) != 0) ? 0.6
				: (mMol.getAtomicNo(mGlobalAtom[atom]) != 6) ? 0.25 : 0.0;
	}

	protected ArrayList<int[]> getCollisionList() {
		mCollisionPanalty = 0.0;
		ArrayList<int[]> collisionList = new ArrayList<int[]>();
		for (int i = 1; i< mGlobalAtom.length; i++) {
			for (int j=0; j<i; j++) {
				double xdif = Math.abs(mAtomX[i]-mAtomX[j]);
				double ydif = Math.abs(mAtomY[i]-mAtomY[j]);
				double dist = Math.sqrt(xdif * xdif + ydif * ydif);
				if (dist < cCollisionLimitBondRotation) {
					int[] collidingAtom = new int[2];
					collidingAtom[0] = mGlobalAtom[i];
					collidingAtom[1] = mGlobalAtom[j];
					collisionList.add(collidingAtom);
				}
				double panalty = 1.0 - Math.min(dist, 1.0);
				mCollisionPanalty += panalty * panalty;
			}
		}
		return collisionList;
	}

	protected double getCollisionPanalty() {
		return mCollisionPanalty;
	}

	protected void locateBonds() {
		int fragmentBonds = 0;
		for (int i=0; i<mGlobalAtom.length; i++) {
			int atom = mGlobalAtom[i];
			int connAtoms = mMol.getAllConnAtoms(atom);
			for (int j=0; j<connAtoms; j++)
				if (mMol.getConnAtom(atom, j) > atom)
					fragmentBonds++;
//			for (int j=connAtoms; j<mMol.getAllConnAtomsPlusMetalBonds(atom); j++)
//				if (mMol.getConnAtom(atom, j) > atom && isMember(mMol.getConnAtom(atom, j)))
//					fragmentBonds++;
		}

		mGlobalBond = new int[fragmentBonds];
		mGlobalToLocalAtom = new int[mMol.getAllAtoms()];

		fragmentBonds = 0;
		for (int i=0; i<mGlobalAtom.length; i++) {
			int atom = mGlobalAtom[i];
			int connAtoms = mMol.getAllConnAtoms(atom);
			mGlobalToLocalAtom[atom] = i;
			for (int j = 0; j<connAtoms; j++)
				if (mMol.getConnAtom(atom, j) > atom)
					mGlobalBond[fragmentBonds++] = mMol.getConnBond(atom, j);
//			for (int j=connAtoms; j<mMol.getAllConnAtomsPlusMetalBonds(atom); j++)
//				if (mMol.getConnAtom(atom, j) > atom && isMember(mMol.getConnAtom(atom, j)))
//					mGlobalBond[fragmentBonds++] = mMol.getConnBond(atom, j);
		}
	}

	protected void optimizeAtomCoordinates(int atom) {
		double x = mAtomX[atom];
		double y = mAtomY[atom];

		InventorAngle[] collisionForce = new InventorAngle[4];

		int forces = 0;
		for (int i = 0; i< mGlobalBond.length; i++) {
			if (forces >= 4)
				break;

			if (atom == mGlobalToLocalAtom[mMol.getBondAtom(0, mGlobalBond[i])]
			 || atom == mGlobalToLocalAtom[mMol.getBondAtom(1, mGlobalBond[i])])
				continue;

			double x1 = mAtomX[mGlobalToLocalAtom[mMol.getBondAtom(0, mGlobalBond[i])]];
			double y1 = mAtomY[mGlobalToLocalAtom[mMol.getBondAtom(0, mGlobalBond[i])]];
			double x2 = mAtomX[mGlobalToLocalAtom[mMol.getBondAtom(1, mGlobalBond[i])]];
			double y2 = mAtomY[mGlobalToLocalAtom[mMol.getBondAtom(1, mGlobalBond[i])]];
			double d1 = Math.sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y));
			double d2 = Math.sqrt((x2-x)*(x2-x)+(y2-y)*(y2-y));
			double bondLength = Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

			if (d1<bondLength && d2<bondLength) {
				if (x1 == x2) {
					double d = Math.abs(x-x1);
					if (d<cCollisionLimitAtomMovement)
						collisionForce[forces++] = new InventorAngle(InventorAngle.getAngle(x1,y,x,y),
								(cCollisionLimitAtomMovement-d)/2);
				}
				else if (y1 == y2) {
					double d = Math.abs(y-y1);
					if (d<cCollisionLimitAtomMovement)
						collisionForce[forces++] = new InventorAngle(InventorAngle.getAngle(x,y1,x,y),
								(cCollisionLimitAtomMovement-d)/2);
				}
				else {
					double m1 = (y2-y1)/(x2-x1);
					double m2 = -1/m1;
					double a1 = y1-m1*x1;
					double a2 = y-m2*x;
					double xs = (a2-a1)/(m1-m2);
					double ys = m1*xs+a1;
					double d = Math.sqrt((xs-x)*(xs-x)+(ys-y)*(ys-y));
					if (d<cCollisionLimitAtomMovement)
						collisionForce[forces++] = new InventorAngle(InventorAngle.getAngle(xs,ys,x,y),
								(cCollisionLimitAtomMovement-d)/2);
				}
				continue;
			}

			if (d1<cCollisionLimitAtomMovement) {
				collisionForce[forces++] = new InventorAngle(InventorAngle.getAngle(x1,y1,x,y),
						(cCollisionLimitAtomMovement-d1)/2);
				continue;
			}

			if (d2<cCollisionLimitAtomMovement) {
				collisionForce[forces++] = new InventorAngle(InventorAngle.getAngle(x2,y2,x,y),
						(cCollisionLimitAtomMovement-d2)/2);
				continue;
			}
		}

		if (forces > 0) {
			InventorAngle force = CoordinateInventor.getMeanAngle(collisionForce, forces);
			mAtomX[atom] += force.mLength * Math.sin(force.mAngle);
			mAtomY[atom] += force.mLength * Math.cos(force.mAngle);
		}
	}

	/**
	 * @param x
	 * @param y
	 * @return angle
	 */
	protected double calculatePreferredAttachmentAngle(double x, double y, int neighbourAtomCount, double padding) {
		if (size() == 1)
			return 0;

		final double BIN_ANGLE = 2.0 * Math.PI / CIRCULAR_BINS;

		double neighbourRadius = padding
				+ Math.sqrt(neighbourAtomCount);	// assume a little large, because they neighbour exposes its wide side

		double[] distance = new double[CIRCULAR_BINS];
		for (int i = 0; i< mGlobalAtom.length; i++) {
			double angle = InventorAngle.getAngle(x, y, mAtomX[i], mAtomY[i]);
			int bin = correctBin((int)Math.round(angle * CIRCULAR_BINS / (2.0*Math.PI)));
			double dx = x - mAtomX[i];
			double dy = y - mAtomY[i];
			double sd = dx*dx + dy*dy;
			if (distance[bin] < sd)
				distance[bin] = sd;
			}
		double maxDistance = -1;
		int maxBin = -1;
		for (int i=0; i<CIRCULAR_BINS; i++) {
			distance[i] = Math.sqrt(distance[i]);
			if (maxDistance < distance[i]) {
				maxDistance = distance[i];
				maxBin = i;
				}
			}

		int preferredBin = correctBin(maxBin - CIRCULAR_BINS / 2);
		// give the bin opposite of the farthest atom a tiny bias
		for (int i=0; i<=CIRCULAR_BINS / 2; i++) {
			distance[correctBin(preferredBin+i)] += 0.01 * i;
			distance[correctBin(preferredBin-i)] += 0.01 * i;
			}

		int neighbourCount = CIRCULAR_BINS / 4;
		double[] sin = new double[neighbourCount];
		double[] cos = new double[neighbourCount];
		for (int i=1; i<neighbourCount; i++) {
			sin[i] = Math.sin(i * BIN_ANGLE);
			cos[i] = Math.cos(i * BIN_ANGLE);
			}

		double squareRadius = neighbourRadius * neighbourRadius;

		double minDistance = Double.MAX_VALUE;
		int minBin = -1;
		for (int bin=0; bin<CIRCULAR_BINS; bin++) {
			if (distance[bin] >= minDistance)
				continue;

			double localMinDistance = distance[bin];

			// check, whether localMinDistance is compatible with adjacent bins and adapt, if needed
			for (int i=1; i<neighbourCount; i++) {
				for (int j=-1; j<=1; j+=2) {	// both sides
					int neighbourBin = correctBin(bin + j*i);

					if (distance[neighbourBin] * cos[i] <= localMinDistance)	// this does not conflict
						continue;

					double d = cos[i] * Math.min(distance[neighbourBin], neighbourRadius / sin[i]);
					if (localMinDistance < d) {
						localMinDistance = d;
						if (minDistance <= localMinDistance)
							break;
						}
					}
				if (minDistance <= localMinDistance)
					break;
				}

			if (minDistance > localMinDistance) {
				minDistance = localMinDistance;
				minBin = bin;
				}
			}

		return Math.PI * 2 * minBin / CIRCULAR_BINS;
	}

	private int correctBin(int bin) {
		return bin < 0 ? bin + CIRCULAR_BINS : bin >= CIRCULAR_BINS ? bin - CIRCULAR_BINS : bin;
	}
}
