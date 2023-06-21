package com.actelion.research.chem.sar;

import com.actelion.research.chem.*;
import com.actelion.research.chem.coords.CoordinateInventor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;

public class CoreBasedSARAnalyzer {
	public static final boolean DISTINGUISH_STEREO_CENTERS = true;
	private static final int MAX_R_GROUPS = 16;

	private StereoMolecule mQuery,mFragment;
	private SSSearcher mSearcher;
	private SSSearcherWithIndex mSearcherWithIndex;
	private ArrayList<MoleculeData> mMoleculeDataList;  // contains info of analyzed molecules, e.g. substituents and corresponding ScaffoldData
	private TreeMap<String,ScaffoldData> mScaffoldMap;  // map of idccodes of core structures to corresponding ScaffoldData
	private ScaffoldGroup mScaffoldGroup;

	/**
	 * This class runs a complete structure-activity-relationship (SAR) analysis from molecules that share
	 * a common similar scaffold. For one given query substructure this class analyses many molecules, whether
	 * the query substructure is found and which substituents are connected at which positions.
	 * Then, R-group numbers are assigned to those query structure positions, which carry varying substituents
	 * throughout matches of different molecule. Every molecule is broken down into a scaffold structure with
	 * R-group pseudo atoms (R1, R2, ...) and the substituent structures, which are present in the particular
	 * molecule.<br>
	 * For a given molecule a scaffold structure is constructed the following way:<br>
	 * - A substructure search locates all matches of the query structure and selects a preferred match based
	 *   on the substitution pattern. If no match is found than this molecule is skipped from the analysis<br>
	 * - All matching atoms and bonds are taken as a 'core' structure. If the query contains bridge bonds, then
	 *   those molecule atoms that match on the bridge bonds also belong to the core structure.<br>
	 * - When the core structure is determined for a molecule, then all remaining atoms belong to substituents.
	 *   For every exit vector of the core structure, i.e. a bond that connects to a substituent atom, the
	 *   substituent is determined. A substituent atom may connect back to another exit vector of the core
	 *   structure, causing a ring closure.<br>
	 * - To all exit vectors that have at least two different substituents (e.g. -H and -Me) an R-group is
	 *   attached to the core structure. If all molecules have the same substituent at one position, then that
	 *   substituent is attached to the core structure. R-groups are numbered accordingly.<br>
	 * - The core structure with attached R-groups and attached constant substituents constitutes the scaffold
	 *   structure for a particular molecule.<br>
	 * - All molecules that match to the same query structure don't necessarily share the same scaffold structure,
	 *   e.g. if a query bridge bond matches on a different chain length or if a wildcard atom matches on a
	 *   different atom type, then the constructed scaffold structure differs in these aspects. Also the
	 *   number of attached R-groups may be different for different scaffolds, if all molecules belonging
	 *   to one scaffold have no substituent at a position that is substituted on a related scaffold's molecules.
	 *   However, it is assured, that R-group numbering is always the same among a scaffold group.<br>
	 *   <b>Summary:</b> If molecules match to the same query structure, they belong to the same scaffold group,
	 *   but nonetheless, their assigned scaffold structures may differ concerning atom types, ring sizes, and
	 *   count of attached R-groups. The R-group numbering (R1, R2, ...), however, is compatible, i.e. R-groups
	 *   at equivalent positions have the same number.<br>
	 * @param query
	 */
	public CoreBasedSARAnalyzer(StereoMolecule query) {
		mQuery = query;
		mQuery.ensureHelperArrays(Molecule.cHelperNeighbours);

		mMoleculeDataList = new ArrayList<>();
		mScaffoldMap = new TreeMap<>();
		mScaffoldGroup = new ScaffoldGroup(query);

		mFragment = new StereoMolecule();   // used as molecule buffer
	}

	public int addMolecule(StereoMolecule mol) {
		if (mSearcher == null) {
			mSearcher = new SSSearcher();
			mSearcher.setFragment(mQuery);
		}

		mSearcher.setMolecule(mol);
		int matchCount = mSearcher.findFragmentInMolecule(SSSearcher.cCountModeRigorous, SSSearcher.cDefaultMatchMode);
		if (matchCount == 0) {
			mMoleculeDataList.add(null);
			return 0;
			}

		addMolecule(mol, mSearcher);

		return matchCount;
	}

	public int addMolecule(StereoMolecule mol, long[] ffp) {
		if (mSearcherWithIndex == null) {
			mSearcherWithIndex = new SSSearcherWithIndex();
			mSearcherWithIndex.setFragment(mQuery, (long[])null);
		}

		mSearcherWithIndex.setMolecule(mol, ffp);
		int matchCount = mSearcherWithIndex.findFragmentInMolecule(SSSearcher.cCountModeRigorous, SSSearcher.cDefaultMatchMode);
		if (matchCount == 0) {
			mMoleculeDataList.add(null);
			return 0;
			}

		addMolecule(mol, mSearcherWithIndex.getGraphMatcher());

		return matchCount;
	}

	public int addMolecule(byte[] idcode, byte[] coords, long[] ffp) {
		if (mSearcherWithIndex == null) {
			mSearcherWithIndex = new SSSearcherWithIndex();
			mSearcherWithIndex.setFragment(mQuery, (long[])null);
		}

		mSearcherWithIndex.setMolecule(idcode, ffp);
		int matchCount = mSearcherWithIndex.findFragmentInMolecule(SSSearcher.cCountModeRigorous, SSSearcher.cDefaultMatchMode);
		if (matchCount == 0) {
			mMoleculeDataList.add(null);
			return 0;
			}

		addMolecule(new IDCodeParser(true).getCompactMolecule(idcode, coords), mSearcherWithIndex.getGraphMatcher());

		return matchCount;
	}

	private void addMolecule(StereoMolecule mol, SSSearcher searcher) {
		int match = findPreferredMatch(mol, searcher.getMatchList());

		int[] queryToMolAtom = searcher.getMatchList().get(match);

		// Mark all atoms belonging to core fragment
		boolean[] isCoreAtom = new boolean[mol.getAtoms()];
		for (int i=0; i<queryToMolAtom.length; i++)
			if (queryToMolAtom[i] != -1)
				isCoreAtom[queryToMolAtom[i]] = true;

		boolean[] isBridgeAtom = searcher.getMatchingBridgeBondAtoms(match);
		if (isBridgeAtom != null)
			for (int i=0; i<isBridgeAtom.length; i++)
				if (isBridgeAtom[i])
					isCoreAtom[i] = true;

		// This (core) is basically the query match, cut out of the real molecules.
		// In case of bridge bonds, it has more atoms than the query itself!
		// In case of exclude groups in the query, it has fewer atoms than the query itself!
		StereoMolecule core = new StereoMolecule();
		int[] molToCoreAtom = new int[isCoreAtom.length];
		mol.copyMoleculeByAtoms(core, isCoreAtom, true, molToCoreAtom);
		Canonizer coreCanonizer = new Canonizer(core);
		int[] coreGraphIndex = coreCanonizer.getGraphIndexes();
		for (int i=0; i<molToCoreAtom.length; i++)
			if (molToCoreAtom[i] != -1)
				molToCoreAtom[i] = coreGraphIndex[molToCoreAtom[i]];
		core = coreCanonizer.getCanMolecule(false);

/*  TODO copy ESR features to scaffold stereo centers with Rn-groups

		core.ensureHelperArrays(Molecule.cHelperNeighbours);
		String extendedCoreIDCode = null;
		int[] coreAtomParity = null;

		if (DISTINGUISH_STEREO_CENTERS) {
			boolean[] isExtendedCoreAtom = new boolean[mol.getAtoms()];	// core plus direct neighbours
			for (int coreAtom=0; coreAtom<core.getAtoms(); coreAtom++) {
				int molAtom = coreToMolAtom[coreAtom];
				isExtendedCoreAtom[molAtom] = true;
				for (int j=0; j<mol.getConnAtoms(molAtom); j++)
					isExtendedCoreAtom[mol.getConnAtom(molAtom, j)] = true;
				}

			StereoMolecule extendedCore = new StereoMolecule();	// core plus direct neighbours
			int[] molToExtendedCoreAtom = new int[mol.getAtoms()];
			mol.copyMoleculeByAtoms(extendedCore, isExtendedCoreAtom, true, molToExtendedCoreAtom);

			// Mark atomicNo of non-core atoms with atomicNo=0 to make sure that the atom order of the
			// core atoms within canonical extendedCore matches the one in canonical core:
			// This way we can directly copy parities from extendedCore to the new molecule constructed from core.
			for (int atom=0; atom<mol.getAtoms(); atom++)
				if (isExtendedCoreAtom[atom] && !isCoreAtom[atom])
					extendedCore.setAtomicNo(molToExtendedCoreAtom[atom], 0);	// '?'

			Canonizer extendedCoreCanonizer = new Canonizer(extendedCore);
			extendedCoreIDCode = extendedCoreCanonizer.getIDCode();
			int[] extendedCoreGraphIndex = extendedCoreCanonizer.getGraphIndexes();
			for (int i=0; i<molToExtendedCoreAtom.length; i++)
				if (molToExtendedCoreAtom[i] != -1)
					molToExtendedCoreAtom[i] = extendedCoreGraphIndex[molToExtendedCoreAtom[i]];
			extendedCore = extendedCoreCanonizer.getCanMolecule(false);

			extendedCore.ensureHelperArrays(Molecule.cHelperParities);

			boolean stereoCenterFound = false;
			coreAtomParity = new int[core.getAtoms()];
			for (int coreAtom=0; coreAtom<core.getAtoms(); coreAtom++) {
				int molAtom = coreToMolAtom[coreAtom];
				int ecAtom = molToExtendedCoreAtom[molAtom];
				if (extendedCore.isAtomStereoCenter(ecAtom)) {
					int atomParity = extendedCore.getAtomParity(ecAtom);
					if (atomParity != Molecule.cAtomParityNone)
						stereoCenterFound = true;
					coreAtomParity[coreAtom] = atomParity;
					if (atomParity == Molecule.cAtomParity1
							|| atomParity == Molecule.cAtomParity2) {
						int esrType = extendedCore.getAtomESRType(ecAtom);
						if (esrType != Molecule.cESRTypeAbs) {
							int esrEncoding = (extendedCore.getAtomESRGroup(ecAtom) << 4)
									+ ((esrType == Molecule.cESRTypeAnd) ? 4 : 8);
							coreAtomParity[coreAtom] += esrEncoding;
							}
						}
					}
				}
			if (!stereoCenterFound)
				coreAtomParity = null;
			else
				extendedCoreIDCode = new Canonizer(extendedCore).getIDCode();
			}*/

		core.setFragment(false);

		// Get existing or create new ScaffoldData object for this particular core structure.
		ScaffoldData scaffoldData = getScaffoldData(core, queryToMolAtom, molToCoreAtom, isBridgeAtom != null);

		mMoleculeDataList.add(new MoleculeData(mol, scaffoldData, molToCoreAtom, mFragment));
		}

	private ScaffoldData getScaffoldData(StereoMolecule core, int[] queryToMolAtom, int[] molToCoreAtom, boolean hasBridgeAtoms) {
		// Create one CoreInfo for every distinct core structure and assign all rows
		// having the same core structure to the respective CoreInfo.
		String coreIDCode = new Canonizer(core).getIDCode();
		ScaffoldData scaffoldData = mScaffoldMap.get(coreIDCode);

		if (scaffoldData == null) {
			int[] queryToCoreAtom = new int[mQuery.getAtoms()];
			int[] coreToQueryAtom = new int[core.getAtoms()];
			Arrays.fill(queryToCoreAtom, -1);	// account for exclude atoms
			Arrays.fill(coreToQueryAtom, -1);	// account for bridge atoms
			for (int queryAtom=0; queryAtom<mQuery.getAtoms(); queryAtom++) {
				int molAtom = queryToMolAtom[queryAtom];
				if (molAtom != -1) {
					int coreAtom = molToCoreAtom[molAtom];
					queryToCoreAtom[queryAtom] = coreAtom;
					coreToQueryAtom[coreAtom] = queryAtom;
					}
				}

// TODO make sure, we don't destroy stereo centers when changin coordinates
			adaptCoreAtomCoordsFromQuery(mQuery, core, queryToCoreAtom, hasBridgeAtoms);

			System.out.println(coreIDCode);
			scaffoldData = new ScaffoldData(core, coreToQueryAtom, queryToCoreAtom, mScaffoldGroup);
			mScaffoldMap.put(coreIDCode, scaffoldData);
			mScaffoldGroup.add(scaffoldData);
			}

		return scaffoldData;
		}

	/**
	 * @param firstRGroup
	 * @return true if the maximum R-group count was exceeded in some cases
	 */
	public boolean analyze(int firstRGroup) {
		boolean rGroupCountExceeded = false;

		// check for varying substituents to require a new column
		for (MoleculeData moleculeData:mMoleculeDataList)
			if (moleculeData != null)
				moleculeData.checkSubstituents();

		// Remove substituents from atoms, which didn't see a varying substitution
		for (MoleculeData moleculeData:mMoleculeDataList)
			if (moleculeData != null)
				moleculeData.removeUnchangingSubstituents();

		int scaffoldRGroupCount = mScaffoldGroup.assignRGroups(firstRGroup);

		for (ScaffoldData scaffoldData:mScaffoldGroup) {
			int rGroupsOnBridges = scaffoldData.assignRGroupsToBridgeAtoms(firstRGroup + scaffoldRGroupCount);

			if (scaffoldRGroupCount + rGroupsOnBridges > MAX_R_GROUPS) {
				for (MoleculeData moleculeData:mMoleculeDataList)
					if (moleculeData != null
					 && moleculeData.getScaffoldData().getRGroupCount() > MAX_R_GROUPS)
						moleculeData.clear();

				rGroupCountExceeded = true;
				}

			scaffoldData.addRGroupsToCoreStructure();
			}

		for (MoleculeData moleculeData:mMoleculeDataList)
			if (moleculeData != null)
				moleculeData.correctSubstituentRingClosureLabels();

		return rGroupCountExceeded;
		}

	/**
	 * Uses a simple strategy to determine the preferred match:
	 * It preferrers matches that bears substituent at low atom indexes.
	 * @param mol
	 * @param matchList
	 * @return
	 */
	private int findPreferredMatch(StereoMolecule mol, ArrayList<int[]> matchList) {
		if (matchList.size() == 1)
			return 0;

		int bestMatch = -1;
		int bestScore = Integer.MAX_VALUE;

		mol.ensureHelperArrays(Molecule.cHelperNeighbours);

		for (int i=0; i<matchList.size(); i++) {
			int[] match = matchList.get(i);
			int score = 0;
			for (int atom:match)
				if (atom != -1)
					score += atom * mol.getConnAtoms(atom);
			if (bestScore > score) {
				bestScore = score;
				bestMatch = i;
			}
		}

		return bestMatch;
	}

	private void adaptCoreAtomCoordsFromQuery(StereoMolecule query, StereoMolecule core, int[] queryToCoreAtom, boolean hasBridgeAtoms) {
		if (!hasBridgeAtoms) {
			// just copy query atom coordinates and mark them to be untouched for later coordinate invention
			for (int queryAtom = 0; queryAtom<queryToCoreAtom.length; queryAtom++) {
				if (queryToCoreAtom[queryAtom] != -1) {
					int coreAtom = queryToCoreAtom[queryAtom];
					core.setAtomX(coreAtom, query.getAtomX(queryAtom));
					core.setAtomY(coreAtom, query.getAtomY(queryAtom));
					core.setAtomMarker(coreAtom, true);  // to later keep the original query coordinates
				}
			}
		} else {
			// Generate new core coordinates and flip and rotate to closely match query orientation
			core.ensureHelperArrays(Molecule.cHelperParities);
			new CoordinateInventor().invent(core);
			double[] cogQuery = new double[2];
			double[] cogCore = new double[2];
			int sharedAtomCount = 0;
			for (int queryAtom = 0; queryAtom<queryToCoreAtom.length; queryAtom++) {
				if (queryToCoreAtom[queryAtom] != -1) {
					int coreAtom = queryToCoreAtom[queryAtom];
					cogCore[0] += core.getAtomX(coreAtom);
					cogCore[1] += core.getAtomY(coreAtom);
					cogQuery[0] += query.getAtomX(queryAtom);
					cogQuery[1] += query.getAtomY(queryAtom);
					sharedAtomCount++;
				}
			}
			cogCore[0] /= sharedAtomCount;
			cogCore[1] /= sharedAtomCount;
			cogQuery[0] /= sharedAtomCount;
			cogQuery[1] /= sharedAtomCount;

			double[] weight = new double[sharedAtomCount];
			double[] rotation = new double[sharedAtomCount];
			double[] flippedRotation = new double[sharedAtomCount];
			int index = 0;
			for (int queryAtom = 0; queryAtom<queryToCoreAtom.length; queryAtom++) {
				if (queryToCoreAtom[queryAtom] != -1) {
					int coreAtom = queryToCoreAtom[queryAtom];

					double cx = core.getAtomX(coreAtom) - cogCore[0];
					double cy = core.getAtomY(coreAtom) - cogCore[1];
					double squareDistCore = cx * cx + cy * cy;

					double qx = query.getAtomX(queryAtom) - cogQuery[0];
					double qy = query.getAtomY(queryAtom) - cogQuery[1];
					double squareDistQuery = qx * qx + qy * qy;

					weight[index] = Math.sqrt(squareDistCore * squareDistQuery);

					double angleQuery = Molecule.getAngle(cogQuery[0], cogQuery[1], query.getAtomX(queryAtom), query.getAtomY(queryAtom));
					double angleCore = Molecule.getAngle(cogCore[0], cogCore[1], core.getAtomX(coreAtom), core.getAtomY(coreAtom));
					rotation[index] = Molecule.getAngleDif(angleCore, angleQuery);
					flippedRotation[index] = Molecule.getAngleDif(-angleCore, angleQuery);

					index++;
				}
			}

			double meanRotation = 0.0;
			double meanFlippedRotation = 0.0;
			double weightSum = 0.0;
			for (int i = 0; i<index; i++) {
				meanRotation += weight[i] * rotation[i];
				meanFlippedRotation += weight[i] * flippedRotation[i];
				weightSum += weight[i];
			}
			meanRotation /= weightSum;
			meanFlippedRotation /= weightSum;

			double penalty = 0.0;
			double flippedPanalty = 0.0;
			for (int i = 0; i<index; i++) {
				penalty += weight[i] * Math.abs(Molecule.getAngleDif(rotation[i], meanRotation));
				flippedPanalty += weight[i] * Math.abs(Molecule.getAngleDif(flippedRotation[i], meanFlippedRotation));
			}

			if (penalty < flippedPanalty) {
				core.zoomAndRotateInit(cogCore[0], cogCore[1]);
				core.zoomAndRotate(1.0, meanRotation, false);
			}
			else {
				for (int coreAtom=0; coreAtom<core.getAllAtoms(); coreAtom++)
					core.setAtomX(coreAtom, 2.0 * cogCore[0] - core.getAtomX(coreAtom));
				for (int coreBond=0; coreBond<core.getAllBonds(); coreBond++)
					if (core.isStereoBond(coreBond))
						core.setBondType(coreBond, core.getBondType(coreBond) == Molecule.cBondTypeUp ? Molecule.cBondTypeDown : Molecule.cBondTypeUp);
				core.zoomAndRotateInit(cogCore[0], cogCore[1]);
				core.zoomAndRotate(1.0, meanFlippedRotation, false);
			}

			for (int coreAtom=0; coreAtom<core.getAllAtoms(); coreAtom++)
				core.setAtomMarker(coreAtom, true);  // to later keep the original query coordinates
		}
	}

	public int getRGroupCount() {
		int count = 0;
		for (ScaffoldData scaffoldData:mScaffoldGroup)
			count = Math.max(count, scaffoldData.getRGroupCount());

		return count;
	}

	public ScaffoldGroup getScaffolds() {
		return mScaffoldGroup;
	}

	public ArrayList<MoleculeData> getMoleculeData() {
		return mMoleculeDataList;
	}
}
