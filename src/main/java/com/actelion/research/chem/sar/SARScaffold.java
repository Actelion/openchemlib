package com.actelion.research.chem.sar;

import com.actelion.research.chem.*;
import com.actelion.research.chem.coords.CoordinateInventor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;

public class SARScaffold {
	private StereoMolecule mQuery,mCore,mScaffold;
	private int[] mCoreToQueryAtom,mQueryToCoreAtom;
	private String mIDCodeWithRGroups, mIDCoordsWithRGroups;
	private TreeMap<String,String> mOldToNewMap;
	private int mBridgeAtomRGroupCount;
	private SARScaffoldGroup mScaffoldGroup;
	private ExitVector[] mBridgeAtomExitVector;
	private boolean[] mHasSeenSubstituentOnScaffold;
	private int[] mSeenBondOrdersOnScaffold;

	protected SARScaffold(StereoMolecule query, StereoMolecule core, int[] coreToQueryAtom, int[] queryToCoreAtom, SARScaffoldGroup scaffoldGroup) {
		mQuery = query;
		mCore = core;
		mCoreToQueryAtom = coreToQueryAtom;
		mQueryToCoreAtom = queryToCoreAtom;
		mScaffoldGroup = scaffoldGroup;
		mOldToNewMap = new TreeMap<>();
		mBridgeAtomRGroupCount = -1;
		analyzeAtomBridgeExitVectors(coreToQueryAtom);
		mHasSeenSubstituentOnScaffold = new boolean[getExitVectorCount()];
		mSeenBondOrdersOnScaffold = new int[getExitVectorCount()];
	}

	private void analyzeAtomBridgeExitVectors(int[] coreToQueryAtom) {
		ArrayList<ExitVector> evList = new ArrayList<>();
		for (int atom=0; atom<mCore.getAtoms(); atom++) {
			if (coreToQueryAtom[atom] == -1) {
				int exitVectorCount = mCore.getImplicitHydrogens(atom);
				for (int i=0; i<exitVectorCount; i++)
					evList.add(new ExitVector(atom, false, i, i));
				}
			}
		mBridgeAtomExitVector = evList.toArray(new ExitVector[0]);
		}

	/**
	 * Checks, whether the coreAtom that carries the substituent is part of a bridge bond
	 * or whether it is shared by the entire scaffold group. In the first case it is checked,
	 * whether the substituent was not yet seen on the level of this CoreInfo.
	 * In the second case it is checked, whether the substituent was not yet seen on the level
	 * of the scaffold group, i.e. all CoreInfos belonging to the same scaffold query.
	 * @param substituent
	 * @param exitVectorIndex which includes bridge bond atoms
	 * @param bondOrder
	 */
	protected void checkSubstituent(String substituent, int exitVectorIndex, int bondOrder) {
		if (substituent != null) {
			mHasSeenSubstituentOnScaffold[exitVectorIndex] = true;
			mSeenBondOrdersOnScaffold[exitVectorIndex] |= (1 << bondOrder);
			}

		getExitVector(exitVectorIndex).checkSubstituent(substituent, bondOrder);
		}

	public StereoMolecule getCoreStructure() {
		return mCore;
	}

	public int getCoreAtom(int exitVectorIndex) {
		return getExitVector(exitVectorIndex).getCoreAtom(mQueryToCoreAtom);
	}

	/**
	 * Determines that atom of the molecule that represents the given exitVectorIndex, which, if the
	 * root atom for that exit vector is a stereo center, considers the exit vector's topicity.
	 * While the root atom has a counterpart in the core structure, the exit vector atom has none.
	 * @param mol
	 * @param coreToMolAtom
	 * @param molToCoreAtom
	 * @param exitVectorIndex
	 * @return the exit vector atom index on the given molecule or -1, if there is no substituent at that exit vector
	 */
	public int getExitVectorAtom(StereoMolecule mol, int[] coreToMolAtom, int[] molToCoreAtom, int exitVectorIndex) {
		ExitVector exitVector = getExitVector(exitVectorIndex);
		int rootAtom = coreToMolAtom[exitVector.getCoreAtom(mQueryToCoreAtom)];

		// If we don't have a stereo center at rootAtom and if one of the exiting bonds is a double or triple bond,
		// then we associate the first exit vector at rootAtom with the double/triple bond and a second
		// (single bonded) exit vector with the remaining bond if one exists.
		boolean hasExitPiBond = false;
		for (int i=0; i<mol.getConnAtoms(rootAtom); i++) {
			if (molToCoreAtom[mol.getConnAtom(rootAtom, i)] == -1
			 && mol.getConnBondOrder(rootAtom, i) > 1) {
				hasExitPiBond = true;
				break;
			}
		}
		if (hasExitPiBond
		 && !mol.isAtomStereoCenter(rootAtom)) {
			for (int i=0; i<mol.getConnAtoms(rootAtom); i++) {
				int connAtom = mol.getConnAtom(rootAtom, i);
				if (molToCoreAtom[connAtom] == -1
				 && ((exitVector.getIndex() == 0 && mol.getConnBondOrder(rootAtom, i) > 1)
				  || (exitVector.getIndex() == 1 && mol.getConnBondOrder(rootAtom, i) == 1)))
					return connAtom;
			}
		}

		int count = 0;
		for (int i=0; i<mol.getConnAtoms(rootAtom); i++) {
			int connAtom = mol.getConnAtom(rootAtom, i);
			if (molToCoreAtom[connAtom] == -1) {
				int topicity = (exitVector.getTopicity() == -1) ? -1 : calculateTopicity(mol, rootAtom, connAtom, molToCoreAtom);
				if (topicity != -1) {
					if (topicity == exitVector.getTopicity())
						return connAtom;
				}
				else if (count == exitVector.getIndex())
					return connAtom;

				count++;
			}
		}
		return -1;
	}

	/**
	 * @param exitVectorIndex which includes bridge bond atoms
	 * @return
	 */
	public ExitVector getExitVector(int exitVectorIndex) {
		return exitVectorIndex < mScaffoldGroup.getExitVectorCount() ?
				  mScaffoldGroup.getExitVector(exitVectorIndex)
				: mBridgeAtomExitVector[exitVectorIndex-mScaffoldGroup.getExitVectorCount()];
	}

	/**
	 * @return total exit vector count including exit vectors on matching bridge bond atoms
	 */
	public int getExitVectorCount() {
		return mScaffoldGroup.getExitVectorCount() + mBridgeAtomExitVector.length;
	}

	/**
	 * Returns the exit vector index from given coreAtom and connIndex or topicity.
	 * @param coreAtom the atom on the core structure that carries one or more exit vectors
	 * @param connIndex index to be used in case of multiple homotopic exit vectors at the coreAtom
	 * @param topicity stereo descriptor that distinguishes two exit vectors in case of a stereo center
	 * @return total index into combined list of exit vectors including those on matching bridge bond atoms
	 */
	public int getExitVectorIndex(int coreAtom, int connIndex, int topicity) {
		if (mCoreToQueryAtom[coreAtom] != -1)
			return mScaffoldGroup.getExitVectorIndex(mCoreToQueryAtom[coreAtom], connIndex, topicity);

		for (int i=0; i<mBridgeAtomExitVector.length; i++)
			if (mBridgeAtomExitVector[i].getCoreAtom(null) == coreAtom
			 && (((topicity == -1 || mBridgeAtomExitVector[i].getTopicity() == -1) && mBridgeAtomExitVector[i].getIndex() == connIndex)
			  || (topicity != -1 && mBridgeAtomExitVector[i].getTopicity() == topicity)))
				return mScaffoldGroup.getExitVectorCount() + i;

		return -1;
	}

	public StereoMolecule getScaffoldWithRGroups() {
		return mScaffold;
	}

	public String getIDCodeWithRGroups() {
		if (mIDCodeWithRGroups == null)
			buildIDCodeAndCoords();

		return mIDCodeWithRGroups;
	}

	public String getIDCoordsWithRGroups() {
		if (mIDCodeWithRGroups == null)
			buildIDCodeAndCoords();

		return mIDCoordsWithRGroups;
	}

	protected TreeMap<String,String> getOldToNewMap() {
		return mOldToNewMap;
	}

	protected int assignRGroupsToBridgeAtoms(int firstBridgeAtomRGroup) {
		if (mBridgeAtomRGroupCount == -1) {
			mBridgeAtomRGroupCount = 0;
			for (ExitVector exitVector:mBridgeAtomExitVector)
				if (exitVector.substituentVaries())
					exitVector.setRGroupNo(++mBridgeAtomRGroupCount + firstBridgeAtomRGroup - 1);
		}
		return mBridgeAtomRGroupCount;
	}

	protected void addRGroupsToCoreStructure() {
		mScaffold = new StereoMolecule(mCore);
		mScaffold.ensureHelperArrays(Molecule.cHelperNeighbours);

		double coreAVBL = mScaffold.getAverageBondLength();

		int exitVectorCount = getExitVectorCount();
		boolean[] closureCovered = new boolean[exitVectorCount];
		for (int exitVectorIndex=0; exitVectorIndex<exitVectorCount; exitVectorIndex++) {
			//	if substituent varies within the scaffold group => attach an R group
			ExitVector exitVector = getExitVector(exitVectorIndex);
			if (exitVector.substituentVaries()) {
				// But don't attach an R-group to one scaffold of a scaffold group,
				// if that particular scaffold has never substituents at that position.
				if (mHasSeenSubstituentOnScaffold[exitVectorIndex]) {
					int rGroupNo = exitVector.getRGroupNo();
					int newAtom = mScaffold.addAtom((rGroupNo<=3) ? 141 + rGroupNo : 125 + rGroupNo);
					int bondType = calculateExitVectorCoordsAndBondType(exitVectorIndex, mScaffold.getAtomCoordinates(newAtom));
					mScaffold.addBond(exitVector.getCoreAtom(mQueryToCoreAtom), newAtom, bondType);
				}
			}
			else {	//	else => attach the non-varying substituent (if it is not null = 'unsubstituted')
				if (!closureCovered[exitVectorIndex] && exitVector.getConstantSubstituent() != null) {
					StereoMolecule substituent = new IDCodeParser(true).getCompactMolecule(exitVector.getConstantSubstituent());

					// Substitutions, which connect back to the core fragment are decorated with atomicNo=0 atoms that
					// carry a label with the respective exit vector index. Here we just copy those connection atoms,
					// but mark the exit vector indexes as already attached (closureCovered) to avoid processing from the
					// other end again. When all substituents are attached, we convert those labelled atoms into
					// proper closure connections.
					for (int atom=0; atom<substituent.getAllAtoms(); atom++) {
						String label = substituent.getAtomCustomLabel(atom);
						if (label != null)
							closureCovered[Integer.parseInt(label)] = true;
					}

					// If we have a stereo bond to connect the substituent
					Coordinates coords = new Coordinates();
					int bondType = calculateExitVectorCoordsAndBondType(exitVectorIndex, coords);

					int rootAtom = exitVector.getCoreAtom(mQueryToCoreAtom);
					double wantedAngle = Molecule.getAngle(mScaffold.getAtomX(rootAtom), mScaffold.getAtomY(rootAtom), coords.x, coords.y);

					mScaffold.addSubstituent(substituent, rootAtom, wantedAngle, bondType);
				}
			}
		}

		// Now we convert all pseudo atoms carrying closure labels into proper ring closures.
		boolean labelsFound = false;
		mScaffold.ensureHelperArrays(Molecule.cHelperNeighbours);
		for (int atom=mCore.getAllAtoms(); atom<mScaffold.getAllAtoms(); atom++) {
			String label = mScaffold.getAtomCustomLabel(atom);
			if (label != null) {
				labelsFound = true;
				int exitVectorIndex = Integer.parseInt(label);
				int firstAtom = mScaffold.getConnAtom(atom, 0);
				int bondType = calculateExitVectorCoordsAndBondType(exitVectorIndex, mScaffold.getAtomCoordinates(firstAtom));

				// Correct parity of first substituent atom at ring closure if we have an odd number of neighbours
				// with an atom index between old and new attachment index.
				int parity = mScaffold.getAtomParity(firstAtom);
				if (parity == Molecule.cAtomParity1 || parity == Molecule.cAtomParity2) {
					boolean invert = false;
					for (int i=0; i<mScaffold.getConnAtoms(firstAtom); i++) {
						int connAtom = mScaffold.getConnAtom(firstAtom, i);
						if (connAtom != atom && connAtom<atom)
							invert = !invert;
					}
					if (invert)
						mScaffold.setAtomParity(firstAtom, parity == Molecule.cAtomParity1 ?
										Molecule.cAtomParity2 : Molecule.cAtomParity1, mScaffold.isAtomParityPseudo(firstAtom));
				}
				mScaffold.addBond(getExitVector(exitVectorIndex).getCoreAtom(mQueryToCoreAtom), firstAtom, bondType);
				mScaffold.markAtomForDeletion(atom);
			}
		}
		if (labelsFound)
			mScaffold.deleteMarkedAtomsAndBonds();

		// To close the ring and to ensure the stereo center correctness at core atom of the closure bond
		// we have relocated the substituent atom of the closure bond potentially inverting the stereo configuration
		// if that atom is a stereo center. We rebuild up/down bonds to reflect new atom coordinates.
		mScaffold.ensureHelperArrays(Molecule.cHelperRings);
		for (int atom=mCore.getAllAtoms(); atom<mScaffold.getAllAtoms(); atom++)
			mScaffold.setStereoBondFromAtomParity(atom);

		mScaffold.ensureHelperArrays(Molecule.cHelperParities);

		for (int atom=0; atom<mCore.getAllAtoms(); atom++)
			mScaffold.setAtomMarker(atom, true);
		new CoordinateInventor(CoordinateInventor.MODE_KEEP_MARKED_ATOM_COORDS).invent(mScaffold);
		for (int atom=0; atom<mCore.getAllAtoms(); atom++)
			mScaffold.setAtomMarker(atom, false);
	}

	public int getRGroupCount() {
		return mScaffoldGroup.getRGroupCount() + mBridgeAtomRGroupCount;
	}

	private void buildIDCodeAndCoords() {
		Canonizer canonizer = new Canonizer(mScaffold);
		mIDCodeWithRGroups = canonizer.getIDCode();
		mIDCoordsWithRGroups = canonizer.getEncodedCoordinates();
	}

	/**
	 * If an R-group is attached to a stereo center (or an E or Z double bond) then this method
	 * calculates a topicity value (0 or 1), which assigns an exit vector to one of the two
	 * possible stereo variants. The calculation uses atom coordinates and up/down bond information
	 * from the molecule mapped to the corresponding atom indexes of the query structure.
	 * If the stereo center in a query uses a bridge bond to one of its neighbours, then
	 * that particular neighbour is using the atom coordinates of the first molecule atom in the
	 * bridge.
	 * @param mol
	 * @param rootAtom
	 * @param exitAtom
	 * @param molToCoreAtom
	 * @return 0 or 1 (-1 if topicity cannot be determined)
	 */
	protected int calculateTopicity(StereoMolecule mol, int rootAtom, int exitAtom, int[] molToCoreAtom) {
		if (!mol.isAtomStereoCenter(rootAtom))
			return mol.getAtomPi(rootAtom) == 0 ? -1 : calculateEZTopicity(mol, rootAtom, exitAtom, molToCoreAtom);

		int[] neighbourRank = new int[4];
		int[] neighbourBond = new int[4];
		double[] neighbourAngle = new double[4];
		int exitBond = -1;
		boolean otherExitAtomFound = false;

		int stereoBond = mol.getStereoBond(rootAtom);
		if (stereoBond == -1)
			return -1;

		// Create array of rootAtom neighbour bond angles sorted by the relevant neighbour atom indexes.
		// Included neighbours are all neighbours that are part of the core structure plus the defined exit atom.
		// A potential second exit atom is not considered here.
		int totalNeighbourCount = 0;
		int coreNeighbourCount = 0;
		for (int i=0; i<mol.getConnAtoms(rootAtom); i++) {
			int connAtom = mol.getConnAtom(rootAtom, i);
			int connBond = mol.getConnBond(rootAtom, i);

			int rank;
			if (connAtom == exitAtom) {
				rank = Integer.MAX_VALUE-1;
				exitBond = connBond;
			}
			else if (molToCoreAtom[connAtom] == -1) {
				if (otherExitAtomFound) // don't allow more than one other exit atom
					return -1;

				rank = Integer.MAX_VALUE;
				otherExitAtomFound = true;
			}
			else {
				rank = getTopicityRelevantAtomIndex(molToCoreAtom[rootAtom], molToCoreAtom[connAtom]);
				coreNeighbourCount++;
			}

			int index = totalNeighbourCount;
			while (index > 0 && neighbourRank[index-1] > rank) {
				neighbourRank[index] = neighbourRank[index-1];
				neighbourBond[index] = neighbourBond[index-1];
				neighbourAngle[index] = neighbourAngle[index-1];
				index--;
			}

			neighbourRank[index] = rank;
			neighbourBond[index] = connBond;
			neighbourAngle[index] = mol.getBondAngle(rootAtom, connAtom);

			totalNeighbourCount++;
		}

		if (totalNeighbourCount < 3 || totalNeighbourCount > 4)
			return -1;

		int stereoType = (mol.getBondType(stereoBond) == Molecule.cBondTypeUp) ? 2 : 1;

		// Here we have one of the following neighbour counts in addition to the defined exitAtom:
		// A: 1 neighbour that exists in core structure; 1 additional exit atom
		// B: 2 neighbours that exists in core structure; no additional exit atom
		// C: 2 neighbours that exists in core structure; 1 additional exit atom
		// D: 3 neighbours that exists in core structure; no additional exit atom

		// Case A and B:
		// 3 non-H neighbours at stereo center.
		// Thus, we don't need to change up/down bond type when shifting stereo bond to other neighbour.

		// Cases C and D:
		// 4 non-H neighbours at stereo center.
		// Thus, we need to invert up/down bond type when shifting stereo bond to direct neighbour bond.
		if (mol.getConnAtoms(rootAtom) == 4) {
			if (otherExitAtomFound && stereoBond == neighbourBond[3]) {
				if (areDirectNeighbours(neighbourAngle, 2, 3))
					stereoType = 3 - stereoType;
			}
			else if (stereoBond != exitBond) {
				// if the stereo bond is one of the core bonds and if the exit bond
				int stereoBondIndex = -1;
				for (int i=0; i<totalNeighbourCount; i++) {
					if (stereoBond == neighbourBond[i]) {
						stereoBondIndex = i;
						break;
					}
				}

				if (areDirectNeighbours(neighbourAngle, stereoBondIndex, 2))
					stereoType = 3 - stereoType;
			}
		}

		return calculateTHTopicity(neighbourAngle, coreNeighbourCount, totalNeighbourCount, stereoType);
	}

	/**
	 * Assuming, we have 4 neighbour angles, this method checks, whether angle[index1]
	 * and angle[index2] are direct neighbours, which means that the other two angles
	 * are direct neighbours as well and lie together on one side between angle[index1]
	 * and angle[index2].
	 * @param angle
	 * @param index1
	 * @param index2
	 * @return
	 */
	private boolean areDirectNeighbours(double[] angle, int index1, int index2) {
		double angle1 = Math.min(angle[index1], angle[index2]);
		double angle2 = Math.max(angle[index1], angle[index2]);
		int count = 0;
		for (int i=0; i<4; i++)
			if (i != index1 && i != index2 && (angle[i] > angle1) && (angle[i] < angle2))
				count++;
		return count != 1;
	}

	/**
	 * If a substituent can be attached to a core structure in two distinguishable ways regarding
	 * stereo configuration, then this method calculates in a reproducible way a topicity (0 or 1)
	 * reflecting a given up/down bond stereo type and the bond angles from the stereo center to
	 * all neighbour atoms. The bond angle array is expected to contain sorted core neighbours first,
	 * followed by one or two exit vector angles. The last exit vector is considered to carry the
	 * stereo bond. If in reality the stereo bond is a different one, then the calling method is
	 * responsible to compensate for a stereo bond shift and/or to compensate for the second exit
	 * vector by potentially inverting stereoType.
	 * Topicity is defined as follows:<br>
	 * - 1 neighbour atoms in core structure and 2 exit atoms:<br>
	 *   If walking from first atom (lowest relevant index) via root atom to first exit atom
	 *   making a left turn, then an up-bond connecting second exit atom gives topicity=1<br>
	 * - 2 neighbour atoms in core structure and 1 or 2 exit atoms:<br>
	 *   If walking from first atom (lowest relevant index) via root atom to second atom
	 *   making a left turn, then an up-bond connecting first/only exist atom gives topicity=1<br>
	 * - 3 neighbour atoms in core structure and one exit atom:<br>
	 *   If neighbours 1,2,3 are in counter-clockwise order and forth neighbour is connected
	 *   with an up-bond, then topicity=1<br>
	 * @param angle bond angles at stereo center sorted by relevant atom index
	 * @param coreNeighbourCount number of bonds at stereo center that have a counterpart in the core structure
	 * @param totalNeighbourCount number of bonds at stereo center including a second exit atom
	 * @param stereoType 1 (down) or 2 (up); corrected if original stereo bond is not the exit bond or/and second exit atom exists
	 * @return topicity 0 or 1
	 */
	private int calculateTHTopicity(double[] angle, int coreNeighbourCount, int totalNeighbourCount, int stereoType) {
		for (int i=1; i<totalNeighbourCount; i++)
			if (angle[i] < angle[0])
				angle[i] += Math.PI*2;

		if (coreNeighbourCount <= 2) {
			boolean leftTurn = (angle[1] - angle[0] > Math.PI);
			return leftTurn ^ (stereoType == 2) ? 1 : 0;
		}

		boolean clockwise = (angle[1] > angle[2]);
		return clockwise ^ (stereoType == 2) ? 1 : 0;
	}

	private int calculateEZTopicity(StereoMolecule mol, int rootAtom, int exitAtom, int[] molToCoreAtom) {
		if (mol.getAtomPi(rootAtom) != 1)
			return -1;

		if (mol.getBondOrder(mol.getBond(rootAtom, exitAtom)) != 1)
			return -1;

		int doubleBond = -1;
		int rearDBAtom = -1;
		for (int i=0; i<mol.getConnAtoms(rootAtom); i++) {
			if (mol.getConnBondOrder(rootAtom, i) == 2
			 && molToCoreAtom[mol.getConnAtom(rootAtom, i)] != -1) {
				doubleBond = mol.getConnBond(rootAtom, i);
				rearDBAtom = mol.getConnAtom(rootAtom, i);
				break;
			}
		}

		if (doubleBond == -1
		 || mol.getBondParity(doubleBond) == Molecule.cBondParityNone
		 || mol.getBondParity(doubleBond) == Molecule.cBondParityUnknown)
			return -1;

		double exitAngle = mol.getBondAngle(rootAtom, exitAtom);
		double dbAngle = mol.getBondAngle(rootAtom, rearDBAtom);
		double angleDif1 = Molecule.getAngleDif(exitAngle, dbAngle);

		int oppositeAtom = Integer.MAX_VALUE;
		double oppositeAngle = 0;
		for (int i=0; i<mol.getConnAtoms(rearDBAtom); i++) {
			int connAtom = mol.getConnAtom(rearDBAtom, i);
			if (connAtom != rootAtom) {
				int candidate = getTopicityRelevantAtomIndex(molToCoreAtom[rearDBAtom], molToCoreAtom[connAtom]);
				if (oppositeAtom > candidate) {
					oppositeAtom = candidate;
					oppositeAngle = mol.getBondAngle(connAtom, rearDBAtom);
				}
			}
		}

		if (oppositeAtom == Integer.MAX_VALUE)
			return -1;

		double angleDif2 = Molecule.getAngleDif(oppositeAngle, dbAngle);

		return (angleDif1 < 0) ^ (angleDif2 < 0) ? 0 : 1;   // E:0  Z:1
	}

	/**
	 * If coreRootAtom matches a stereo center in the molecule and if coreConnAtom is one of coreRootAtom's
	 * neighbours in the core structure, then this method returns for this neighbour that atom index, which
	 * is used to determine the topicity for any exit vector (neighbours in mol, which are not part of the
	 * core structure). Typically, we use query structure atom indexes for this, i.e. the relevant atom index
	 * of a coreConnAtom is the index of its correscponding atom in the query structure. If, however,
	 * coreConnAtom is an atom of a matching bridge bond, then there is no corresponding query structure atom.
	 * In that case we walk along the bridge bond atoms in the core structure until we hit an atom that exists
	 * in the query, which is the remote bridge bond atom, whose index is then returned.
	 * @param coreRootAtom
	 * @param coreConnAtom
	 * @return
	 */
	private int getTopicityRelevantAtomIndex(int coreRootAtom, int coreConnAtom) {
		int queryRoot = mCoreToQueryAtom[coreRootAtom];

		// If the stereo center itself is not part of the query, then it is within a bridge bond
		// and does not exist in all scaffolds of the scaffold group.
		// In this case we use atom index of the core rather than the query as reference.
		if (queryRoot == -1)
			return coreConnAtom;

		int queryAtom = mCoreToQueryAtom[coreConnAtom];
		if (queryAtom != -1)
			return queryAtom;

		// If the stereo center neighbour in the core does not exist in the query, then it is part of
		// a bridge bond. In this case we have to find that stereo center neighbour in the query
		// that is connected with that bridge bond in the core, which copntains coreConnAtom.
		// For that we build a graph from the core root atom adding only atoms that don't exist in
		// the query until we hit a query core neighbour, which we return.

		int[] bridgeNeighbour = new int[mQuery.getConnAtoms(queryRoot)];
		int bridgeNeighbourCount = 0;
		for (int i=0; i<mQuery.getConnAtoms(queryRoot); i++)
			if (mQuery.isBondBridge(mQuery.getConnBond(queryRoot, i))
			 && mQueryToCoreAtom[mQuery.getConnAtom(queryRoot, i)] != -1)
				bridgeNeighbour[bridgeNeighbourCount++] = mQueryToCoreAtom[mQuery.getConnAtom(queryRoot, i)];

		int[] graphAtom = new int[mCore.getAtoms()];
		boolean[] atomUsed = new boolean[mCore.getAtoms()];

		graphAtom[0] = coreConnAtom;
		atomUsed[coreRootAtom] = true;
		atomUsed[coreConnAtom] = true;

		int current = 0;
		int highest = 0;
		while (current <= highest) {
			int parent = graphAtom[current];
			for (int i=0; i<mCore.getConnAtoms(parent); i++) {
				int candidate = mCore.getConnAtom(parent, i);
				for (int j=0; j<bridgeNeighbourCount; j++)
					if (candidate == bridgeNeighbour[j])
						return mCoreToQueryAtom[bridgeNeighbour[j]];

				if (mCoreToQueryAtom[candidate] == -1
				 && !atomUsed[candidate]) {
					graphAtom[++highest] = candidate;
					atomUsed[candidate] = true;
				}
			}
			current++;
		}

		return -1;  // should not happen
	}

	/**
	 * When attaching an R-group or another substituent to the unsubstituted core structure,
	 * we need to use a stereo bond (up or down) if we create a stereo center.
	 * From the exit vector's topicity information we can deduce whether to use up or down:<br>
	 * We define: If we have increasing atom indexes of query bonds in clockwise order,
	 * then topicity=0 is associated with an UP-bond and topicity=1 is associated with a DOWN-bond.
	 * @param exitVectorIndex
	 * @param coords receives suggested coordinates for first exit atom
	 * @return
	 */
	private int calculateExitVectorCoordsAndBondType(int exitVectorIndex, Coordinates coords) {
		if ((mSeenBondOrdersOnScaffold[exitVectorIndex] & 2) == 0)
			return ((mSeenBondOrdersOnScaffold[exitVectorIndex] & 4) == 0) ? 3 : 2;

		ExitVector exitVector = getExitVector(exitVectorIndex);
		int rootAtom = exitVector.getCoreAtom(mQueryToCoreAtom);

		int[] neighbour = new int[3];
		double[] angle = new double[3];

		int coreNeighbourCount = 0;
		int piBondSum = 0;
		for (int i=0; i<mScaffold.getConnAtoms(rootAtom); i++) {
			int connAtom = mScaffold.getConnAtom(rootAtom, i);
			if (connAtom < mCoreToQueryAtom.length) {
				int neighbourAtom = getTopicityRelevantAtomIndex(rootAtom, connAtom);

				int index = coreNeighbourCount;
				while (index > 0 && neighbour[index-1] > neighbourAtom) {
					neighbour[index] = neighbour[index-1];
					angle[index] = angle[index-1];
					index--;
				}

				neighbour[index] = neighbourAtom;
				angle[index] = mScaffold.getBondAngle(rootAtom, connAtom);
				piBondSum += mScaffold.getConnBondOrder(rootAtom, i) - 1;
				coreNeighbourCount++;
			}
		}

		if (piBondSum != 0) {
			if (coreNeighbourCount == 1)
				calculateSP2ExitVectorCoords(rootAtom, piBondSum, exitVector.getTopicity(), angle, coords);
			else
				calculateSP3ExitVectorCoords(rootAtom, Arrays.copyOf(angle, coreNeighbourCount), coords);

			// we assume that we don't have a stereo center with double bonds, e.g. at S or P
			return Molecule.cBondTypeSingle;
		}

		calculateSP3ExitVectorCoords(rootAtom, Arrays.copyOf(angle, coreNeighbourCount), coords);

		if (exitVector.getTopicity() == -1)
			return Molecule.cBondTypeSingle;

		angle[coreNeighbourCount] = Molecule.getAngle(mScaffold.getAtomX(rootAtom), mScaffold.getAtomY(rootAtom), coords.x, coords.y);

		int totalNeighbourCount = coreNeighbourCount + 1;

		if (coreNeighbourCount == 1) {
			angle[coreNeighbourCount+1] = angle[coreNeighbourCount] + Math.PI * 2 / 3;
			totalNeighbourCount++;
		}

		int topicity = calculateTHTopicity(angle, coreNeighbourCount, totalNeighbourCount, 1);
		return (topicity == -1) ? Molecule.cBondTypeSingle : (topicity == exitVector.getTopicity()) ? Molecule.cBondTypeDown : Molecule.cBondTypeUp;
	}

	/**
	 * If there are at least two existing neighbours (angle.length >= 2), this method
	 * places the new neighbour at a position furthest away from any existing neighbour
	 * using the scaffolds average bond length. If there is only one neighbour, the new
	 * neighbour will be placed with a bond angle 120 degrees larger.
	 * @param atom
	 * @param angle
	 * @param coords
	 */
	private void calculateSP3ExitVectorCoords(int atom, double[] angle, Coordinates coords) {
		double exitAngle = Math.PI * 2 / 3;

		if (angle.length >= 2) {
			Arrays.sort(angle);
			double largestDiff = -1.0;
			for (int i=0; i<angle.length; i++) {
				double angleDiff = (i == 0) ? angle[0] + 2*Math.PI - angle[angle.length-1] : angle[i] - angle[i-1];
				if (largestDiff<angleDiff) {
					largestDiff = angleDiff;
					exitAngle = angle[i] - 0.5 * angleDiff;
				}
			}
		}
		else if (angle.length == 1) {
			exitAngle = angle[0] + Math.PI * 2 / 3;
		}

		double avbl = mScaffold.getAverageBondLength();
		coords.x = mScaffold.getAtomX(atom) + avbl * Math.sin(exitAngle);
		coords.y = mScaffold.getAtomY(atom) + avbl * Math.cos(exitAngle);
	}

	/**
	 * Assuming exactly one existing neighbour connected with a pi bond,
	 * and using the scaffolds average bond length, this method
	 * places the new neighbour at a suitable position considering hybridisation and topicity.
	 * @param atom
	 * @param piBondSum
	 * @param angle
	 * @param coords
	 */
	private void calculateSP2ExitVectorCoords(int atom, int piBondSum, int topicity, double[] angle, Coordinates coords) {
		double exitAngle = 0.0;

		if (piBondSum == 2) { // triple bond
			exitAngle = angle[0] + Math.PI;
		}
		else if (topicity == -1) {
			exitAngle = angle[0] + 0.6667 * Math.PI;
		}
		else {
			int rearDBAtom = mScaffold.getConnAtom(atom, 0);
			double dbAngle = mScaffold.getBondAngle(rearDBAtom, atom);
			int oppositeAtom = Integer.MAX_VALUE;
			double oppositeAngle = 0;
			for (int i=0; i<mScaffold.getConnAtoms(rearDBAtom); i++) {
				int connAtom = mScaffold.getConnAtom(rearDBAtom, i);
				if (connAtom != atom) {
					int candidate = getTopicityRelevantAtomIndex(rearDBAtom, connAtom);
					if (oppositeAtom > candidate) {
						oppositeAtom = candidate;
						oppositeAngle = mScaffold.getBondAngle(rearDBAtom, connAtom);
					}
				}
			}

			double angleDif = Molecule.getAngleDif(oppositeAngle, dbAngle);
			exitAngle = angle[0] + ((angleDif < 0) ^ (topicity == 1) ? 0.6667 : 1.3333) * Math.PI;
		}

		double avbl = mScaffold.getAverageBondLength();
		coords.x = mScaffold.getAtomX(atom) + avbl * Math.sin(exitAngle);
		coords.y = mScaffold.getAtomY(atom) + avbl * Math.cos(exitAngle);
	}
}
