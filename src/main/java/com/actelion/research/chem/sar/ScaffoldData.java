package com.actelion.research.chem.sar;

import com.actelion.research.chem.*;
import com.actelion.research.chem.coords.CoordinateInventor;

import java.util.ArrayList;
import java.util.TreeMap;

import static com.actelion.research.chem.coords.CoordinateInventor.MODE_PREFER_MARKED_ATOM_COORDS;

public class ScaffoldData {
	private StereoMolecule mCore,mScaffold;
	private int[] mCoreToQueryAtom,mQueryToCoreAtom;
	private String mIDCodeWithRGroups, mIDCoordsWithRGroups;
	private TreeMap<String,String> mOldToNewMap;
	private int mBridgeAtomRGroupCount;
	private ScaffoldGroup mScaffoldGroup;
	private ExitVector[] mBridgeAtomExitVector;
	private boolean[] mHasSeenSubstituentOnScaffold;

	protected ScaffoldData(StereoMolecule core, int[] coreToQueryAtom, int[] queryToCoreAtom, ScaffoldGroup scaffoldGroup) {
		mCore = core;
		mCoreToQueryAtom = coreToQueryAtom;
		mQueryToCoreAtom = queryToCoreAtom;
		mScaffoldGroup = scaffoldGroup;
		mOldToNewMap = new TreeMap<>();
		mBridgeAtomRGroupCount = -1;
		analyzeAtomBridgeExitVectors(coreToQueryAtom);
		mHasSeenSubstituentOnScaffold = new boolean[getExitVectorCount()];
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
	 */
	protected void checkSubstituent(String substituent, int exitVectorIndex) {
		if (substituent != null)
			mHasSeenSubstituentOnScaffold[exitVectorIndex] = true;

		getExitVector(exitVectorIndex).checkSubstituent(substituent);
		}

	public StereoMolecule getCoreStructure() {
		return mCore;
	}

	public int getCoreAtom(int exitVectorIndex) {
		return getExitVector(exitVectorIndex).getCoreAtom(mQueryToCoreAtom);
	}

	public int getExitVectorAtom(StereoMolecule mol, int[] coreToMolAtom, int[] molToCoreAtom, int exitVectorIndex) {
		ExitVector exitVector = getExitVector(exitVectorIndex);
		int coreAtom = exitVector.getCoreAtom(mQueryToCoreAtom);
		int rootAtom = coreToMolAtom[coreAtom];
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
	 * @param coreAtom
	 * @param topicity
	 * @return total index into combined list of exit vectors including those on matching bridge bond atoms
	 */
	public int getExitVectorIndex(int coreAtom, int connIndex, int topicity) {
		if (mCoreToQueryAtom[coreAtom] != -1)
			return mScaffoldGroup.getExitVectorIndex(mCoreToQueryAtom[coreAtom], connIndex, topicity);

		for (int i=0; i<mBridgeAtomExitVector.length; i++)
			if (mBridgeAtomExitVector[i].getCoreAtom(null) == coreAtom
			 && ((topicity == -1 && mBridgeAtomExitVector[i].getIndex() == connIndex)
			  || (topicity != -1 && mBridgeAtomExitVector[i].getTopicity() == topicity)))
				return i;

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
		// We use the same StereoMolecule for both, first the core structure, then for the scaffold,
		// which is the decorated core.
		// By using two references and setting the other one to null, we make sure to generate errors,
		// if we try accessing the wrong one at a wrong time.
		mScaffold = mCore;
		mCore = null;

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
					int bondType = calculateExitVectorCoordsAndBondType(exitVector, mScaffold.getAtomCoordinates()[newAtom]);
					mScaffold.addBond(exitVector.getCoreAtom(mQueryToCoreAtom), newAtom, bondType);

Canonizer can = new Canonizer(new StereoMolecule(mScaffold));
System.out.println(can.getIDCode()+" "+can.getEncodedCoordinates());
				}
			}
			else {	//	else => attach the non-varying substituent (if it is not null = 'unsubstituted')
				if (!closureCovered[exitVectorIndex] && exitVector.getConstantSubstituent() != null) {
					StereoMolecule substituent = new IDCodeParser(true).getCompactMolecule(exitVector.getConstantSubstituent());

					// If we have a stereo bond to connect the substituent
					Coordinates coords = new Coordinates();
					int bondType = calculateExitVectorCoordsAndBondType(exitVector, coords);

					// Translate substituent to correct attachment position
					for (int a=0; a<substituent.getAllAtoms(); a++) {
						if (substituent.getAtomicNo(a) == 0 && substituent.getAtomCustomLabel(a) == null)
							substituent.translateCoords(coords.x - substituent.getAtomX(a), coords.y - substituent.getAtomY(a));
						break;
					}

					// Substitutions, which connect back to the core fragment are encoded with labels on the connecting atoms: "core atom index".
					// Now we translate labels back to atomMapNos, which are used by addSubstituent() to create back connections.
					for (int a=0; a<substituent.getAllAtoms(); a++) {
						String label = substituent.getAtomCustomLabel(a);
						if (label != null) {
							int atom = Integer.parseInt(label);
							substituent.setAtomCustomLabel(a, (String)null);
							substituent.setAtomicNo(a, 0);
							substituent.setAtomMapNo(a, atom+1, false);
							closureCovered[atom] = true;
						}
					}

// TODO ringClosure in mapNo is outdated. Do better...
					int coreAtom = exitVector.getCoreAtom(mQueryToCoreAtom);
					mScaffold.addSubstituent(substituent, coreAtom, true);

					if ((bondType & Molecule.cBondTypeMaskStereo) != 0) {
						for (int bond=mScaffold.getAllBonds()-substituent.getAllBonds(); bond<mScaffold.getAllBonds(); bond++) {
							if (mScaffold.getBondAtom(0, bond) == coreAtom
							 || mScaffold.getBondAtom(1, bond) == coreAtom) {
								mScaffold.setBondType(bond, bondType);
								if (mScaffold.getBondAtom(1, bond) == coreAtom) {
									mScaffold.setBondAtom(1, bond, mScaffold.getBondAtom(0, bond));
									mScaffold.setBondAtom(0, bond, coreAtom);
									break;
								}
							}
						}
					}
				}
			}
		}

		mScaffold.ensureHelperArrays(Molecule.cHelperParities);

		// TODO we may need to remove overspecifying up/down bonds
	}

	public int getRGroupCount() {
		return mScaffoldGroup.getRGroupCount() + mBridgeAtomRGroupCount;
	}

	private void buildIDCodeAndCoords() {
		mScaffold.ensureHelperArrays(Molecule.cHelperParities); // provide parities for CoordinateInventor
		new CoordinateInventor(MODE_PREFER_MARKED_ATOM_COORDS).invent(mScaffold);

		Canonizer canonizer = new Canonizer(mScaffold);
		mIDCodeWithRGroups = canonizer.getIDCode();
		mIDCoordsWithRGroups = canonizer.getEncodedCoordinates();
	}

	public int calculateTopicity(StereoMolecule mol, int rootAtom, int exitAtom, int[] molToCoreAtom) {
		if (!mol.isAtomStereoCenter(rootAtom))
			return -1;

		int stereoBond = mol.getStereoBond(rootAtom);
if (stereoBond == -1) System.out.println("ERROR: No stereobond found"); // TODO remove

		int[] neighbour = new int[3];
		double[] angle = new double[3];

		int count = 0;
		for (int i=0; i<mol.getConnAtoms(rootAtom); i++) {
			neighbour[count] = mol.getConnAtom(rootAtom, i);
			if (molToCoreAtom[neighbour[count]] != -1) {
				angle[count] = mol.getBondAngle(rootAtom, neighbour[count]);
				count++;
			}
		}

		int bondType = mol.getBondType(stereoBond);

		if (mol.getConnAtoms(rootAtom) == 4) {
			for (int i=0; i<mol.getConnAtoms(rootAtom); i++) {
				int connAtom = mol.getConnAtom(rootAtom, i);
				int connBond = mol.getConnBond(rootAtom, i);

				if (molToCoreAtom[connAtom] == -1
				 && connAtom != exitAtom
				 && stereoBond == connBond) {   // stereoBond is other exit atom
					bondType = bondType == Molecule.cBondTypeDown ? Molecule.cBondTypeUp
							 : bondType == Molecule.cBondTypeUp ? Molecule.cBondTypeDown : bondType;
					break;
					}

				// If we have 2 neighbours in the core and two exit bonds, we assume that both exit bonds roughly point
				// into the same direction and no correction is needed!?

				// TODO if we have 3 neighbours in the core and if stereoBond is not the exit bond, then we also need to correct
				}
			}

		boolean isClockWise = getAngleParity(angle, count) == getOrderParity(neighbour, count);
		return isClockWise ^ bondType == Molecule.cBondTypeDown ? 0 : 1;
	}

	/**
	 * When attaching an R-group or another substituent to the unsubstituted core structure,
	 * we need to use a stereo bond (up or down) if we create a stereo center.
	 * From the exit vector's topicity information we can deduce whether to use up or down:<br>
	 * We define: If we have increasing atom indexes of query bonds in clockwise order,
	 * then topicity=0 is associated with a UP-bond and topicity=1 is associated with a DOWN-bond.
	 * @param exitVector
	 * @param coords receives suggested coordinates for first exit atom
	 * @return
	 */
	private int calculateExitVectorCoordsAndBondType(ExitVector exitVector, Coordinates coords) {
		// TODO what about double and triple

		if (exitVector.getTopicity() == -1)
			return Molecule.cBondTypeSingle;

		int coreAtom = exitVector.getCoreAtom(mQueryToCoreAtom);

		int[] neighbour = new int[3];
		double[] angle = new double[3];

		int count = 0;
		double dx = 0;
		double dy = 0;
		for (int i=0; i<mScaffold.getConnAtoms(coreAtom); i++) {
			neighbour[count] = mScaffold.getConnAtom(coreAtom, i);
			if (neighbour[count] < mCoreToQueryAtom.length) {
				dx += mScaffold.getAtomX(neighbour[count]) - mScaffold.getAtomX(coreAtom);
				dy += mScaffold.getAtomY(neighbour[count]) - mScaffold.getAtomY(coreAtom);
				angle[count] = mScaffold.getBondAngle(coreAtom, neighbour[count]);
				count++;
			}
		}

		// We don't care to have a perfect distance, only the direction is important for stereo centers
		coords.x = mScaffold.getAtomX(coreAtom) - dx;
		coords.y = mScaffold.getAtomY(coreAtom) - dy;
System.out.println("coreAtom:"+coreAtom+" dx:"+dx+" dy:"+dy); // TODO remove

		boolean isClockWise = getAngleParity(angle, count) == getOrderParity(neighbour, count);
		return isClockWise ^ exitVector.getTopicity() == 1 ? Molecule.cBondTypeUp : Molecule.cBondTypeDown;
	}

	/**
	 * @param angle
	 * @param count 2 or 3
	 * @return false if angles are in increasing order
	 */
	private boolean getAngleParity(double[] angle, int count) {
		if (count == 2)
			return Molecule.getAngleDif(angle[0], angle[1]) > 0;

		return (Molecule.getAngleDif(angle[0], angle[1]) > 0 && Molecule.getAngleDif(angle[2], angle[0]) > 0)
			|| (Molecule.getAngleDif(angle[2], angle[0]) > 0 && Molecule.getAngleDif(angle[1], angle[2]) > 0)
			|| (Molecule.getAngleDif(angle[1], angle[2]) > 0 && Molecule.getAngleDif(angle[0], angle[1]) > 0);
	}

	/**
	 * @param atomIndex
	 * @param count 2 or 3
	 * @return false if atom indexes are in natural order or if we have an even number of exchanges to get to natural order
	 */
	private boolean getOrderParity(int[] atomIndex, int count) {
		if (count == 2)
			return atomIndex[0] > atomIndex[1];

		return (atomIndex[0] > atomIndex[1] && atomIndex[2] > atomIndex[0])
			|| (atomIndex[2] > atomIndex[0] && atomIndex[1] > atomIndex[2])
			|| (atomIndex[1] > atomIndex[2] && atomIndex[0] > atomIndex[1]);
	}
}
