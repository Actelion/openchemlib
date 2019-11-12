package com.actelion.research.chem.conf;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.phesa.PheSAAlignment;

public class BondRotationHelper {
	
	private StereoMolecule mMol;
	private int[] mRotatableBonds;
	private boolean[] mIsRotatableBond;
	private int[][] mSmallerSideAtomLists;
	private int[][] mTorsionAtoms;
	private int[][] mRearAtoms;
	private int[] mRotationCenters;

	
	
	public BondRotationHelper(StereoMolecule mol) {
		mMol = mol;
		initialize();
		
		
	}
	
	public void initialize() {
		int[] disconnectedFragmentNo = new int[mMol.getAllAtoms()];
		int disconnectedFragmentCount = mMol.getFragmentNumbers(disconnectedFragmentNo, false, true);
		int disconnectedFragmentSize[] = new int[disconnectedFragmentCount];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			disconnectedFragmentSize[disconnectedFragmentNo[atom]]++;
		mIsRotatableBond = new boolean[mMol.getBonds()];
		TorsionDB.findRotatableBonds(mMol,true, mIsRotatableBond);
		List<Integer> rotBonds = new ArrayList<Integer>();
		IntStream.range(0, mIsRotatableBond.length).forEach(e -> {
			if(mIsRotatableBond[e])
				rotBonds.add(e);
		});
		mRotatableBonds = rotBonds.stream().mapToInt(i->i).toArray();
		mTorsionAtoms = new int[mRotatableBonds.length][];
		mRearAtoms = new int[mRotatableBonds.length][];
		mSmallerSideAtomLists = new int[mRotatableBonds.length][];
		mRotationCenters = new int[mRotatableBonds.length];
		
		for(int i=0;i<mRotatableBonds.length;i++) {
			int bond = mRotatableBonds[i];
			int[] torsionAtoms = new int[4];
			int[] rearAtoms = new int[2];
			TorsionDetail detail = new TorsionDetail();
			if (TorsionDB.getTorsionID(mMol, bond, torsionAtoms, detail) != null) {
				rearAtoms[0] = detail.getRearAtom(0);
				rearAtoms[1] = detail.getRearAtom(1);
				}
			else {
				predictAtomSequence(bond, torsionAtoms, rearAtoms);
				}
			mTorsionAtoms[i] = torsionAtoms;
			mRearAtoms[i] = rearAtoms;
			findSmallerSideAtomList(disconnectedFragmentSize[disconnectedFragmentNo[mMol.getBondAtom(0, bond)]],
					disconnectedFragmentNo, i);
			
		}

		
	}
	
	//populates smallerSideAtomList array and returns the rotation center
	private void findSmallerSideAtomList(int disconnectedFragmentSize, int[] disconnectedFragmentNo, int bondIndex) {
		boolean[] isMember = new boolean[mMol.getAllAtoms()];
		int memberCount = mMol.getSubstituent(mRearAtoms[bondIndex][0], mTorsionAtoms[bondIndex][1], isMember, null, null);
		
		int alkyneAtoms = 0;	// if we have an extended linear sp-atom strain
		if (mRearAtoms[bondIndex][0] != mTorsionAtoms[bondIndex][2])
			alkyneAtoms = mMol.getPathLength(mRearAtoms[bondIndex][0], mTorsionAtoms[bondIndex][2]);

		boolean invert = false;
		if (memberCount > disconnectedFragmentSize-alkyneAtoms-memberCount) {
			memberCount = disconnectedFragmentSize-alkyneAtoms-memberCount;
			invert = true;
			}

		// if invert, then flag all linear alkyne atoms to be avoided
		if (invert && alkyneAtoms != 0) {
			int spAtom = mRearAtoms[bondIndex][0];
			int backAtom = mTorsionAtoms[bondIndex][1];
        	while (mMol.getAtomPi(spAtom) == 2
           		&& mMol.getConnAtoms(spAtom) == 2
           		&& mMol.getAtomicNo(spAtom) < 10) {
        		isMember[spAtom] = true;
           		for (int j=0; j<2; j++) {
           			int connAtom = mMol.getConnAtom(spAtom, j);
           			if (connAtom != backAtom) {
           				backAtom = spAtom;
           				spAtom = connAtom;
           				break;
           				}
           			}
           		}
			}

		int memberNo = 0;
		int fragmentNo = disconnectedFragmentNo[mTorsionAtoms[bondIndex][1]];
		int [] smallerSideAtomList = new int[memberCount];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			if (disconnectedFragmentNo[atom] == fragmentNo && (isMember[atom] ^ invert))
				smallerSideAtomList[memberNo++] = atom;


		int rotationCenter = mTorsionAtoms[bondIndex][invert ? 2 : 1];
		mRotationCenters[bondIndex] = rotationCenter;
		mSmallerSideAtomLists[bondIndex] = smallerSideAtomList;
		
		}
	
	
	private void predictAtomSequence(int bond, int[] torsionAtoms, int[] rearAtoms) {
        for (int i=0; i<2; i++) {
    		int centralAtom = mMol.getBondAtom(i, bond);
        	int rearAtom = mMol.getBondAtom(1-i, bond);

        	// walk along sp-chains to first sp2 or sp3 atom
        	while (mMol.getAtomPi(centralAtom) == 2
        		&& mMol.getConnAtoms(centralAtom) == 2
        		&& mMol.getAtomicNo(centralAtom) < 10) {
        		for (int j=0; j<2; j++) {
        			int connAtom = mMol.getConnAtom(centralAtom, j);
        			if (connAtom != rearAtom) {
        				rearAtom = centralAtom;
        				centralAtom = connAtom;
        				break;
        				}
        			}
        		}

        	torsionAtoms[i+1] = centralAtom;
           	rearAtoms[i] = rearAtom;
        	}

    	// A TorsionPrediction does not distinguish hetero atoms from carbons a positions 0 and 3.
        // Therefore we can treat two sp2 neighbors as equivalent when predicting torsions.
        if (mMol.getAtomPi(torsionAtoms[1]) == 0 && mMol.getConnAtoms(torsionAtoms[1]) == 3) {
        	torsionAtoms[0] = -1;
        	}
        else {
			for (int i=0; i<mMol.getConnAtoms(torsionAtoms[1]); i++) {
				int connAtom = mMol.getConnAtom(torsionAtoms[1], i);
				if (connAtom != torsionAtoms[2]) {
					torsionAtoms[0] = connAtom;
					break;
					}
				}
        	}
        if (mMol.getAtomPi(torsionAtoms[2]) == 0 && mMol.getConnAtoms(torsionAtoms[2]) == 3) {
        	torsionAtoms[3] = -1;
        	}
        else {
			for (int i=0; i<mMol.getConnAtoms(torsionAtoms[2]); i++) {
				int connAtom = mMol.getConnAtom(torsionAtoms[2], i);
				if (connAtom != torsionAtoms[1]) {
					torsionAtoms[3] = connAtom;
					break;
					}
				}
        	}
		}
	
	public boolean isRotatableBond(int bond) {
		return mIsRotatableBond[bond];
	}
	
	public void rotateSmallerSide(int bond,double alpha) {
		if(!isRotatableBond(bond))
			return;
		int bondIndex=-1;
		for(int i=0;i<mRotatableBonds.length;i++) {
			if(mRotatableBonds[i]==bond)
				bondIndex = i;
		}
		if(bondIndex==-1)
			return;
		Coordinates t2 = mMol.getCoordinates(mTorsionAtoms[bondIndex][2]);
		Coordinates unit = t2.subC(mMol.getCoordinates(mTorsionAtoms[bondIndex][1])).unit();
		double[][] m = new double[3][3];
		int rotationCenter = mRotationCenters[bondIndex];
		PheSAAlignment.getRotationMatrix((rotationCenter == mTorsionAtoms[bondIndex][1]) ? alpha : -alpha,unit,m);

		for (int atom:mSmallerSideAtomLists[bondIndex]) {
			if (atom != rotationCenter) {
				double x = mMol.getAtomX(atom) - t2.x;
				double y = mMol.getAtomY(atom) - t2.y;
				double z = mMol.getAtomZ(atom) - t2.z;
				mMol.setAtomX(atom, x*m[0][0]+y*m[0][1]+z*m[0][2] + t2.x);
				mMol.setAtomY(atom, x*m[1][0]+y*m[1][1]+z*m[1][2] + t2.y);
				mMol.setAtomZ(atom, x*m[2][0]+y*m[2][1]+z*m[2][2] + t2.z);
				}
			}
		}


}
