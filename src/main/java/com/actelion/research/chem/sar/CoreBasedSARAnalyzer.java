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
	private SARMolecule[] mSARMolecule;  // contains info of analyzed molecules, e.g. substituents and corresponding SARScaffold
	private TreeMap<String, SARScaffold> mScaffoldMap;  // map of core structure idcodes to corresponding SARScaffold
	private SARScaffoldGroup mScaffoldGroup;
	private int[] mPreferredQueryAtomRGroupMatch;

	/**
	 * This class runs a complete structure-activity-relationship (SAR) analysis from molecules that share
	 * one or multiple common similar scaffold(s). For one given query substructure this class analyses many
	 * molecules, whether the query substructure is found and which substituents are connected at which positions.
	 * If the query substructure is found in (matches) multiple molecules, then the matching substructures
	 * may be the same in all cases, or they may differ between some molecules. This happens, for instance,
	 * if the <i>query</i> contains wildcard atoms that match to multiple atom types, or if a bridge bond
	 * matches to atom chains of different lengths.<br>
	 * The matching substructure of a molecule is called <i>core structure</i>. All molecules that share the same
	 * <i>core structure</i> are analyzed regarding which core structure atom carries which substituents within
	 * these molecules. If the substitution is varying within these molecules for a core structrure atom, then
	 * an R-group is assigned to that position. An atom may get multiple R-groups if attachment positions are
	 * diastereotop or if it sees multiple substituents in some of its molecules. The core structure with numbered
	 * attached R-groups at all exit vectors with changing substitution is called a <i>scaffold</i>.
	 * Thus, a specific <i>core structure</i> always gives rise to a specific <i>scaffold</i> structure.
	 * Thus, all N <i>scaffolds</i> derived from all N <i>core structures</i> caused by one <i>query structure</i>
	 * form a <i>scaffold group</i>. R-group numbering between all <i>scaffolds</i> within the same <i>scaffold group</i>
	 * is compatible. This means that R-groups at equivalent exit vectors of two different scaffolds within
	 * the same group have the same number. Two exit vectors are considered equivalent, if they are connected
	 * to atoms, which match to the same query structure atom, and if they are diastereotop, e.g. both connected
	 * with an up-stereo bond when super-positioning their atom coordinates.<br>
	 * For any given molecule a scaffold structure is constructed the following way:<br>
	 * - A substructure search locates all matches of the query structure and selects a preferred match based
	 *   on the substitution pattern. If no match is found, then this molecule is skipped from the analysis.<br>
	 * - All matching atoms and bonds are taken as a 'core' structure. If the query contains bridge bonds, then
	 *   those molecule atoms that match on the bridge bonds also belong to the core structure.<br>
	 * - When the core structure is determined for a molecule, then all remaining atoms of the molecule that
	 *   don't belong to the core structure, are part of substituents.<br>
	 * - For every exit vector of the core structure, i.e. a bond that connects to a substituent atom, the
	 *   substituent structure is determined. A substituent atom may connect back to another exit vector of
	 *   the core structure, causing a ring closure.<br>
	 * - After the analysis of all molecules that share the same <i>core structure</i>, a copy of the
	 *   <i>core structure</i> is created. Then all exit vectors that have at least two different substituents
	 *   throughout all molecules (e.g. -H and -Me) an R-group is attached to the <i>core structure</i>.
	 *   If all molecules have the same substituent at one core structure position, then no R-group is attached.
	 *   Instead, that substituent itself is attached to the core structure.<br>
	 * - The core structure with attached R-groups and attached constant substituents constitutes the scaffold
	 *   structure for a particular molecule.<br>
	 * - All molecules that match to the same query structure don't necessarily share the same scaffold structure,
	 *   e.g. if a query bridge bond matches on a different chain length or if a wildcard atom matches on a
	 *   different atom type, then the constructed scaffold structure differs in these aspects. Also, the
	 *   number of attached R-groups may be different for different scaffolds, if all molecules belonging
	 *   to one scaffold have no substituent at a position that is substituted on a related scaffold's molecules.
	 *   However, it is assured, that R-group numbering is always the same among the entire scaffold group.<br>
	 *   <b>Summary:</b> If molecules match to the same query structure, they belong to the same scaffold group,
	 *   but nonetheless, their assigned scaffold structures may differ concerning atom types, ring sizes, and
	 *   count of attached R-groups. The R-group numbering (R1, R2, ...), however, is compatible, i.e. R-groups
	 *   at equivalent positions have the same number.<br>
	 * @param query substructure with valid atom coordinates that defines one or multiple scaffolds (e.g. via atom lists or bond bridges)
	 */
	public CoreBasedSARAnalyzer(StereoMolecule query, int moleculeCount) {
		mQuery = query;
		mQuery.ensureHelperArrays(Molecule.cHelperNeighbours);

		mSARMolecule = new SARMolecule[moleculeCount];
		mScaffoldMap = new TreeMap<>();
		mScaffoldGroup = new SARScaffoldGroup(query);

		mFragment = new StereoMolecule();   // used as molecule buffer

		mPreferredQueryAtomRGroupMatch = new int[MAX_R_GROUPS];
		Arrays.fill(mPreferredQueryAtomRGroupMatch, -1);
	}

	/**
	 * Adds a molecule to the SAR-analyzer:<br>
	 * - determines core structure from the preferred query match<br>
	 * - if this core structure was not seen yet, creates a new scaffold object for this core structure<br>
	 * - creates new SAR-molecule data object with scaffold and substituent information<br>
	 * Use this version of addMolecule() if you don't have pre-calculated fragment fingerprints
	 * for your molecules available.
	 * @param mol
	 * @param index
	 * @return
	 */
	public int setMolecule(StereoMolecule mol, int index) {
		if (mSearcher == null) {
			mSearcher = new SSSearcher();
			mSearcher.setFragment(mQuery);
		}

		mSearcher.setMolecule(mol);
		int matchCount = mSearcher.findFragmentInMolecule(SSSearcher.cCountModeRigorous, SSSearcher.cDefaultMatchMode);
		if (matchCount == 0)
			return 0;

		setMolecule(mol, mSearcher, index);

		return matchCount;
	}

	/**
	 * Adds a molecule to the SAR-analyzer:<br>
	 * - determines core structure from the preferred query match<br>
	 * - if this core structure was not seen yet, creates a new scaffold object for this core structure<br>
	 * - creates new SAR-molecule data object with scaffold and substituent information<br>
	 * Use this version of addMolecule() if you have in memory molecules and pre-calculated fragment fingerprints.
	 * @param mol
	 * @param ffp
	 * @param index
	 * @return
	 */
	public int setMolecule(StereoMolecule mol, long[] ffp, int index) {
		if (mSearcherWithIndex == null) {
			mSearcherWithIndex = new SSSearcherWithIndex();
			mSearcherWithIndex.setFragment(mQuery, (long[])null);
		}

		mSearcherWithIndex.setMolecule(mol, ffp);
		int matchCount = mSearcherWithIndex.findFragmentInMolecule(SSSearcher.cCountModeRigorous, SSSearcher.cDefaultMatchMode);
		if (matchCount == 0)
			return 0;

		setMolecule(mol, mSearcherWithIndex.getGraphMatcher(), index);

		return matchCount;
	}

	/**
	 * Adds a molecule to the SAR-analyzer:<br>
	 * - determines core structure from the preferred query match<br>
	 * - if this core structure was not seen yet, creates a new scaffold object for this core structure<br>
	 * - creates new SAR-molecule data object with scaffold and substituent information<br>
	 * Use this version of addMolecule() if you have idcodes, coords, and pre-calculated fragment fingerprints
	 * of your molecules.
	 * @param idcode
	 * @param coords
	 * @param ffp
	 * @param index
	 * @return
	 */
	public int setMolecule(byte[] idcode, byte[] coords, long[] ffp, int index) {
		if (mSearcherWithIndex == null) {
			mSearcherWithIndex = new SSSearcherWithIndex();
			mSearcherWithIndex.setFragment(mQuery, (long[])null);
		}

		mSearcherWithIndex.setMolecule(idcode, ffp);
		int matchCount = mSearcherWithIndex.findFragmentInMolecule(SSSearcher.cCountModeRigorous, SSSearcher.cDefaultMatchMode);
		if (matchCount == 0)
			return 0;

		setMolecule(new IDCodeParser(true).getCompactMolecule(idcode, coords), mSearcherWithIndex.getGraphMatcher(), index);

		return matchCount;
	}

	private void setMolecule(StereoMolecule mol, SSSearcher searcher, int index) {
		int match = findPreferredMatch(mol, searcher.getMatchList());

		int[] queryToMolAtom = searcher.getMatchList().get(match);

		// Mark all atoms belonging to core fragment
		boolean[] isCoreAtom = new boolean[mol.getAtoms()];
		for (int i=0; i<queryToMolAtom.length; i++)
			if (queryToMolAtom[i] != -1)
				isCoreAtom[queryToMolAtom[i]] = true;

		// Add all R-group atom directly attached to core (in case of successive SAR steps)
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (isCoreAtom[atom]) {
				for (int i=0; i<mol.getConnAtoms(atom); i++) {
					int connAtom = mol.getConnAtom(atom, i);
					if (!isCoreAtom[connAtom]) {
						int atomicNo = mol.getAtomicNo(connAtom);
						if (atomicNo >= 129 && atomicNo <= 144)
							isCoreAtom[connAtom] = true;
					}
				}
			}
		}

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

		core.setFragment(false);

		// In order to build always the same scaffold molecule from a core structure, we need a canonical core as basis!
		Canonizer coreCanonizer = new Canonizer(core);
		int[] coreGraphIndex = coreCanonizer.getGraphIndexes();
		for (int i=0; i<molToCoreAtom.length; i++)
			if (molToCoreAtom[i] != -1)
				molToCoreAtom[i] = coreGraphIndex[molToCoreAtom[i]];
		StereoMolecule canonicalCore = coreCanonizer.getCanMolecule(false);

		// As key for equal core structures resulting in equal scaffolds, we use the idcode of the core structure.
		// However, we want to have different scaffolds, if substituents are connected with different bond orders.
		// Therefore, we use radical info to mark connection atoms with substituents connected by double/triple bonds.
		for (int atom=0; atom<molToCoreAtom.length; atom++) {
			if (molToCoreAtom[atom] != -1) {
				int exoPiCount = 0;
				for (int i=0; i<mol.getConnAtoms(atom); i++) {
					int connAtom = mol.getConnAtom(atom, i);
					if (molToCoreAtom[connAtom] == -1)
						exoPiCount += mol.getConnBondOrder(atom, i) - 1;
				}

				// encode exo-core bonds with pi-electrons as radicals
				if (exoPiCount != 0)
					core.setAtomRadical(molToCoreAtom[atom], exoPiCount == 1 ?
							Molecule.cAtomRadicalStateD : Molecule.cAtomRadicalStateT);
			}
		}


/* TODO copy ESR features to scaffold stereo centers with Rn-groups???
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

		// Get existing or create new SARScaffold object for this particular core structure.
		String key = new Canonizer(core).getIDCode();
		SARScaffold scaffold = getScaffold(key, canonicalCore, queryToMolAtom, molToCoreAtom, isBridgeAtom != null);

		mSARMolecule[index] = new SARMolecule(mol, scaffold, molToCoreAtom, mFragment);
		}

	/**
	 * Create one SARScaffold for every distinct core structure and assign all rows
	 * having the same core structure to the respective SARScaffold.
	 * @param canonicalCore
	 * @param queryToMolAtom
	 * @param molToCoreAtom
	 * @param hasBridgeAtoms
	 * @return
	 */
	private SARScaffold getScaffold(String key, StereoMolecule canonicalCore, int[] queryToMolAtom, int[] molToCoreAtom, boolean hasBridgeAtoms) {
		SARScaffold scaffold = mScaffoldMap.get(key);

		if (scaffold == null) {
			canonicalCore.ensureHelperArrays(Molecule.cHelperNeighbours);
			int[] queryToCoreAtom = new int[mQuery.getAtoms()];
			int[] coreToQueryAtom = new int[canonicalCore.getAtoms()];
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

			adaptCoreAtomCoordsFromQuery(mQuery, canonicalCore, queryToCoreAtom, hasBridgeAtoms);

			scaffold = new SARScaffold(mQuery, canonicalCore, coreToQueryAtom, queryToCoreAtom, mScaffoldGroup);
			mScaffoldMap.put(key, scaffold);
			mScaffoldGroup.addScaffold(scaffold);
			}

		return scaffold;
		}

	/**
	 * Call this after all molecules have been added to the CoreBasedSARAnalyzer.
	 * This method determines the real scaffold structure(s) and those exit vectors
	 * that carry varying substituents in the molecule set. Then R-group numbers
	 * are assigned to those exit vectors, which carry varying substituents.
	 * If the highest assigned R-group number exceeds the maximum (16) for one or
	 * more scaffolds, then molecules belonging to those scaffolds are not processed.
	 * @param firstRGroupNumber 1 or a higher number, e.g. if the molecules contain R-groups themselves
	 * @return true if all molecules could be processed, i.e. R-group numbers didn't exceed the maximum
	 */
	public boolean analyze(int firstRGroupNumber) {
		boolean rGroupCountExceeded = false;

		// check for varying substituents to require a new column
		for (SARMolecule molecule: mSARMolecule)
			if (molecule != null)
				molecule.checkSubstituents();

		// Remove substituents from atoms, which didn't see a varying substitution
		for (SARMolecule molecule: mSARMolecule)
			if (molecule != null)
				molecule.removeUnchangingSubstituents();

		int scaffoldRGroupCount = mScaffoldGroup.assignRGroups(firstRGroupNumber);

		for (SARScaffold scaffold :mScaffoldGroup.getScaffoldList()) {
			int rGroupsOnBridges = scaffold.assignRGroupsToBridgeAtoms(firstRGroupNumber + scaffoldRGroupCount);

			if (scaffoldRGroupCount + rGroupsOnBridges > MAX_R_GROUPS) {
				for (SARMolecule molecule: mSARMolecule)
					if (molecule != null
					 && molecule.getScaffold().getRGroupCount() > MAX_R_GROUPS)
						molecule.clear();

				rGroupCountExceeded = true;
				}

			scaffold.addRGroupsToCoreStructure();
			}

		for (SARMolecule molecule: mSARMolecule)
			if (molecule != null)
				molecule.correctSubstituentRingClosureLabels();

		return !rGroupCountExceeded;
		}

	/**
	 * Uses a simple strategy to determine the preferred match:
	 * It preferrers matches that carry substituents at low atom indexes.
	 * @param mol
	 * @param matchList
	 * @return
	 */
	private int findPreferredMatch(StereoMolecule mol, ArrayList<int[]> matchList) {
		if (matchList.size() == 1)
			return 0;

		int bestMatch = -1;
		int bestScore = Integer.MIN_VALUE;
		int[] bestQueryAtomRGroupMatch = null;

		mol.ensureHelperArrays(Molecule.cHelperNeighbours);

		for (int i=0; i<matchList.size(); i++) {
			int[] match = matchList.get(i);

			// flag used atoms
			boolean[] isUsedAtom = new boolean[mol.getAtoms()];
			for (int atom:match)
				if (atom != -1)
					isUsedAtom[atom] = true;

			int rGroupMatchScore = mol.getAtoms() * match.length;   // larger than all potential differences from substitution penalties

			int score = 0;
			for (int m=0; m<match.length; m++) {
				int atom = match[m];
				if (atom != -1) {
					// Add penalty for external neighbours that increases with core atom index
					// This prefers atoms with low atom indexes to carry neighbours.
					// We could do better by also looking at substituent sizes, or even optimize
					// substituent structure - atom index combination counts.
					int addedValence = mol.getConnAtoms(atom) - mQuery.getConnAtoms(m);
					if (addedValence>0)
						score -= atom * addedValence;
				}
			}

			// In case of 2-step SAR deconvolutions, where we may have R-groups as substituents,
			// we try to choose those matches, which have the same R-groups at
			// the same positions.
			int[] queryAtomRGroupMatch = getExistingQueryAtomRGroupMatch(match, isUsedAtom, mol);
			if (queryAtomRGroupMatch != null) {
				int matchingRGroupCount = 0;
				for (int k=0; k<MAX_R_GROUPS; k++) {
					if (mPreferredQueryAtomRGroupMatch[k] != -1 && queryAtomRGroupMatch[k] != -1) {
						if (mPreferredQueryAtomRGroupMatch[k] != queryAtomRGroupMatch[k]) {
							matchingRGroupCount = -1;
							break;
						}
						else {
							matchingRGroupCount++;
						}
					}
				}

				score += matchingRGroupCount * rGroupMatchScore;
			}

			if (bestScore < score) {
				bestScore = score;
				bestMatch = i;
				bestQueryAtomRGroupMatch = queryAtomRGroupMatch;
			}
		}

		if (bestQueryAtomRGroupMatch != null)
			for (int i=0; i<MAX_R_GROUPS; i++)
				if (bestQueryAtomRGroupMatch[i] != -1)
					mPreferredQueryAtomRGroupMatch[i] = bestQueryAtomRGroupMatch[i];

		return bestMatch;
	}

	private int[] getExistingQueryAtomRGroupMatch(int[] match, boolean[] isUsedAtom, StereoMolecule mol) {
		int[] rGroupToQueryAtom = null;
		for (int i=0; i<match.length; i++) {
			if (match[i] != -1) {
				for (int j=0; j<mol.getConnAtoms(match[i]); j++) {
					int connAtom = mol.getConnAtom(match[i], j);
					if (!isUsedAtom[connAtom]) {
						// In case of 2-step SAR deconvolutions, where we may have R-groups as substituents,
						// we try to choose those matches, which have the same R-groups at
						// the same positions:
						int atomicNo = mol.getAtomicNo(connAtom);
						if (atomicNo >= 129 && atomicNo <= 144) {
							if (rGroupToQueryAtom == null) {
								rGroupToQueryAtom = new int[MAX_R_GROUPS];
								Arrays.fill(rGroupToQueryAtom, -1);
							}
							int rGroupNo = (atomicNo >= 142) ? atomicNo - 142 : atomicNo - 126; // 0-based
							rGroupToQueryAtom[rGroupNo] = i;
						}
					}
				}
			}
		}
		return rGroupToQueryAtom;
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
		for (SARScaffold scaffold :mScaffoldGroup.getScaffoldList())
			count = Math.max(count, scaffold.getRGroupCount());

		return count;
	}

	public SARScaffoldGroup getScaffoldGroup() {
		return mScaffoldGroup;
	}

	public SARMolecule getMoleculeData(int index) {
		return mSARMolecule[index];
	}
}
