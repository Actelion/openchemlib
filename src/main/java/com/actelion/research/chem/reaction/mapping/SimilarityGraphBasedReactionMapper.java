package com.actelion.research.chem.reaction.mapping;

import com.actelion.research.chem.*;
import com.actelion.research.chem.reaction.Reaction;
import com.actelion.research.util.ByteArrayComparator;

import java.util.Arrays;

public class SimilarityGraphBasedReactionMapper {
	// When building the mapping graph, next neighbors are added based on directed environment similarity.
	// To determine this similarity, directed radial substructures are precalculated for every atom from
	// the perspective of every of one its neighbour atoms (fromAtom). The smallest fragment just contains
	// the neighbor and the atom itself. The next fragment is build from the first by adding all direct
	// neighbours to the atom (except the fromAtom). of neighbor is considered to belong to the graph
	// already
	private static final int SIMILARITY_SHIFT = 8;  // bit shift to allow for lower significant features to distinguish equal environment similarities
	private static final int SKELETON_PENALTY = 2;
	private static final int MAX_ENVIRONMENT_RADIUS = 8;
	private static final int MAX_ENVIRONMENT_SIMILARITY = MAX_ENVIRONMENT_RADIUS << SIMILARITY_SHIFT;
	private static final int MAX_SKELETON_SIMILARITY = (MAX_ENVIRONMENT_RADIUS - SKELETON_PENALTY) << SIMILARITY_SHIFT;

	// Fine-tune environment similarity score modifiers
	private static final int PI_AND_HETERO_PLUS = 1;    // only for skeleton similarity
	private static final int STEREO_MATCH_PLUS = 64;
	private static final int ENDO_RING_PLUS = 128;

	private StereoMolecule mReactant,mProduct;
	private int mMapNo,mGraphMapNoCount,mMappableAtomCount;
	private int[] mReactantMapNo,mProductMapNo,mReactantRingMembership,mProductRingMembership;
	private float mScore;
	private ByteArrayComparator mSimilarityComparator;
	private byte[][][][] mReactantConnAtomEnv,mProductConnAtomEnv;   // indexes: atom,connIndex,radius,idcode bytes
	private byte[][][][] mReactantSkelAtomEnv,mProductSkelAtomEnv;   // indexes: atom,connIndex,radius,idcode bytes

	/**
	 * Entirely maps the given reaction by setting all reactant's and product's mapping numbers.
	 * If the reaction already contains manually mapped atoms, then these are untouched and used
	 * as seed atoms when building the mapping graph.
	 * @param rxn reaction with proper atom coordinates
	 */
	public void map(Reaction rxn) {
		mergeReactantsAndProducts(rxn);
		map(mReactant, mProduct, new int[mReactant.getAtoms()], new int[mProduct.getAtoms()]);
		copyMapNosToReaction(rxn, mReactantMapNo, mProductMapNo, mGraphMapNoCount);
		System.out.println("mapped:"+mMapNo+"("+mGraphMapNoCount+") of "+mMappableAtomCount+" score:"+calculateScore());
		}

	/**
	 * Entirely maps the given reaction by setting all reactant's and product's mapping numbers.
	 * @param reactant all reactants merged into one molecule; may contain seed mapping
	 * @param product all products merged into one molecule; may contain seed mapping
	 * @param reactantMapNo array to receive the created mapping numbers
	 * @param productMapNo array to receive the created mapping numbers
	 */
	public void map(StereoMolecule reactant, StereoMolecule product, int[] reactantMapNo, int[] productMapNo) {
		mReactant = reactant;
		mProduct = product;

		mReactantMapNo = new int[reactantMapNo.length];
		mProductMapNo = new int[productMapNo.length];

		mReactantConnAtomEnv = classifyNeighbourAtomEnvironment(mReactant, false);
		mReactantSkelAtomEnv = classifyNeighbourAtomEnvironment(mReactant, true);
		mProductConnAtomEnv = classifyNeighbourAtomEnvironment(mProduct, false);
		mProductSkelAtomEnv = classifyNeighbourAtomEnvironment(mProduct, true);

		initializeRingMembership();

		mSimilarityComparator = new ByteArrayComparator();
		mScore = -1e10f;

		RootAtomPairSource rootAtomPairSource = new RootAtomPairSource(reactant, product, mReactantMapNo, mProductMapNo);

		while (rootAtomPairSource.hasNextPairSequence()) {
			mMapNo = rootAtomPairSource.getManualMapCount();
			mMappableAtomCount = rootAtomPairSource.getMappableAtomCount();

			RootAtomPair pair = rootAtomPairSource.nextPair();
			while (pair != null) {
				mapFromRootAtoms(pair);
				pair = rootAtomPairSource.nextPair();
				}

			mGraphMapNoCount = mMapNo;

			float score = 0;
			if (mMapNo < mMappableAtomCount) {
				ReactionCenterMapper centerMapper = new ReactionCenterMapper(mReactant, mProduct, mReactantMapNo, mProductMapNo, mMapNo);
				score = centerMapper.completeAndScoreMapping();
				mMapNo += centerMapper.getMappedAtomCount();
				}
			else {
				MappingScorer scorer = new MappingScorer(mReactant, mProduct);
				score = scorer.scoreMapping(scorer.createReactantToProductAtomMap(mReactantMapNo, mProductMapNo));
				}
// TODO confirm that MappingScorer can replace internal scoring
			if (mScore < score) {
				mScore = score;
				for (int i=0; i<reactantMapNo.length; i++)
					reactantMapNo[i] = mReactantMapNo[i];
				for (int i=0; i<productMapNo.length; i++)
					productMapNo[i] = mProductMapNo[i];
				}
			}
		}

	/**
	 * @return number of atoms that could be mapped by graph similarity matching (no reaction center mapping)
	 */
	public int getGraphMapNoCount() {
		return mGraphMapNoCount;
		}

	/**
	 * Calculates and returns a score <= 0 for the current mapping. Higher value (closer to 0) are better.
	 * The score adds -2 for every broken or added bond, and -1 for every changed bond order. It also
	 * adds -1 for every inverted stereo center.
 	 * @return
	 */
	public float getScore() {
		return mScore;
		}

	private int calculateScore() {
		int score = 0;

		int[] mapNoToProduct = new int[mMapNo+1];
		Arrays.fill(mapNoToProduct, -1);
		for (int productAtom=0; productAtom<mProductMapNo.length; productAtom++)
			if (mProductMapNo[productAtom] != 0)
				mapNoToProduct[mProductMapNo[productAtom]] = productAtom;

		boolean[] productBondHandled = new boolean[mProduct.getBonds()];

		for (int reactantBond=0; reactantBond<mReactant.getBonds(); reactantBond++) {
			boolean found = false;
			int productAtom1 = mapNoToProduct[mReactantMapNo[mReactant.getBondAtom(0, reactantBond)]];
			int productAtom2 = mapNoToProduct[mReactantMapNo[mReactant.getBondAtom(1, reactantBond)]];
			if (productAtom1 != -1 && productAtom2 != -1) {
				int productBond = mProduct.getBond(productAtom1, productAtom2);
				if (productBond != -1) {
					found = true;
					productBondHandled[productBond] = true;

					// we give one point for a changed bond type, which includes delocalized as separate bond type
					if (getBondType(mReactant, reactantBond) != getBondType(mProduct, productBond))
						score--;
					}
				}

			// we give two points for removed or added bonds, no matter what type they are
			if (!found)
				score -= 2;
			}

		for (boolean pbh:productBondHandled)
			if (!pbh)
				score -= 2;

		// if we change a stereo center's parity, we must have broken or formed bonds
		for (int reactantAtom=0; reactantAtom<mReactant.getAtoms(); reactantAtom++) {
			int productAtom = mapNoToProduct[mReactantMapNo[reactantAtom]];
			if (productAtom != -1) {
				int reactantParity = mReactant.getAtomParity(reactantAtom);
				if (reactantParity != Molecule.cAtomParityNone) {
					int productParity = mProduct.getAtomParity(productAtom);
					if (reactantParity == Molecule.cAtomParityUnknown) {
						if (productParity == Molecule.cAtomParity1
						 || productParity == Molecule.cAtomParity2)
							score -= 5; // one broken and one formed bond plus additional panelty!
						// must be more expensive than Mitsunobu, which itself must be more expensive than simple esterification (one broken and one formed bond)
						}
					else {
						if (productParity == Molecule.cAtomParityUnknown
						 || isTHParityInversion(reactantAtom, mapNoToProduct) == (reactantParity == productParity))
							score -= 5; // one broken and one formed bond plus additional panelty!
						}
					}
				}
			}

		return score;
		}

	private boolean isTHParityInversion(int reactantAtom, int[] mapNoToProduct) {
		boolean inversion = false;
		if (mReactant.getAtomPi(reactantAtom) == 0) {
			for (int i=1; i<mReactant.getConnAtoms(reactantAtom); i++) {
				for (int j=0; j<i; j++) {
					int connAtom1 = mReactant.getConnAtom(reactantAtom,i);
					int connAtom2 = mReactant.getConnAtom(reactantAtom,j);
					int connMapNo1 = mReactantMapNo[connAtom1];
					int connMapNo2 = mReactantMapNo[connAtom2];
					if (connMapNo1 != -1 && connMapNo2 != -1
					 && ((mapNoToProduct[connMapNo1] > mapNoToProduct[connMapNo2]) ^ (connAtom1 > connAtom2)))
						inversion = !inversion;
					}
				}
			}
/*		else {  // no allene parities for now
			for (int i=0; i<mReactant.getConnAtoms(reactantAtom); i++) {
				int connAtom = mReactant.getConnAtom(reactantAtom, i);
				int neighbours = 0;
				int[] neighbour = new int[3];
				int[] neighbourMapNo = new int[3];
				for (int j=0; j<mReactant.getConnAtoms(connAtom); j++) {
					neighbour[neighbours] = mReactant.getConnAtom(connAtom, j);
					neighbourMapNo[neighbours] = mReactantMapNo[neighbour[neighbours]];
					if (neighbour[neighbours] != reactantAtom)
						neighbours++;
					}
				if (neighbours == 2
				 && neighbourMapNo[0] != -1 && neighbourMapNo[1] != -1
				 && ((mapNoToProduct[neighbourMapNo[0]] > mapNoToProduct[neighbourMapNo[1]])
					^(neighbour[0] > neighbour[1])))
					inversion = !inversion;
				}
			}*/
		return inversion;
		}

	private void initializeRingMembership() {
		mReactantRingMembership = initializeRingMembership(mReactant);
		mProductRingMembership = initializeRingMembership(mProduct);
		}

	private int[] initializeRingMembership(StereoMolecule mol) {
		int[] ringMembership = new int[mol.getAtoms()];
		RingCollection ringSet = mol.getRingSet();
		for (int ring=0; ring<Math.min(32, ringSet.getSize()); ring++) {
			int[] ringAtom = ringSet.getRingAtoms(ring);
			for (int atom:ringAtom)
				ringMembership[atom] |= (1 << ring);
			}
		return ringMembership;
		}

	/**
	 * @param atom1
	 * @param atom2
	 * @param atom3
	 * @return true, if atom3 is not member of any ring that is shared by atom1 and atom2
	 */
	private boolean leavesRing(int atom1, int atom2, int atom3, int[] ringMembership) {
		return (ringMembership[atom1] & ringMembership[atom2] & ~ringMembership[atom3]) != 0;
		}

	/**
	 * Compares reactantAtom to productAtom starting from and extending away from the respective root atoms.
	 * Similarity is increased, if reactantRoot and productRoot are stereo centers and reactantAtom matches
	 * productAtom in regard to the stereo configuration.
	 * @param reactantRoot
	 * @param reactantAtom
	 * @param productRoot
	 * @param productAtom
	 * @return radius<<SIMILARITY_SHIFT; radius 0: mismatch, 1:atomicNo-match, 2:directNeighbourMatch, 3:TwoShellMatch, etc.
	 */
	private int getAtomSimilarity(int reactantRoot, int reactantAtom, int productRoot, int productAtom) {
		int reactantConnIndex = -1;
		for (int i=0; i<mReactant.getConnAtoms(reactantAtom); i++) {
			if (mReactant.getConnAtom(reactantAtom, i) == reactantRoot) {
				reactantConnIndex = i;
				break;
				}
			}

		int productConnIndex = -1;
		for (int i=0; i<mProduct.getConnAtoms(productAtom); i++) {
			if (mProduct.getConnAtom(productAtom, i) == productRoot) {
				productConnIndex = i;
				break;
				}
			}

		for (int radius=0; radius<MAX_ENVIRONMENT_RADIUS; radius++)
			if (mSimilarityComparator.compare(mReactantConnAtomEnv[reactantAtom][reactantConnIndex][radius], mProductConnAtomEnv[productAtom][productConnIndex][radius]) != 0)
				return radius << SIMILARITY_SHIFT;

		return MAX_ENVIRONMENT_SIMILARITY;
		}

	/**
	 * Compares reactantAtom to productAtom starting from and extending away from the respective root atoms.
	 * Similarity is increased, if reactantRoot and productRoot are stereo centers and reactantAtom matches
	 * productAtom in regard to the stereo configuration.
	 * @param reactantRoot
	 * @param reactantAtom
	 * @param productRoot
	 * @param productAtom
	 * @return radius<<SIMILARITY_SHIFT; radius 0: mismatch, 1:atomicNo-match, 2:directNeighbourMatch, 3:TwoShellMatch, etc.
	 */
	private int getSkeletonSimilarity(int reactantRoot, int reactantAtom, int productRoot, int productAtom) {
		int reactantConnIndex = -1;
		for (int i=0; i<mReactant.getConnAtoms(reactantAtom); i++) {
			if (mReactant.getConnAtom(reactantAtom, i) == reactantRoot) {
				reactantConnIndex = i;
				break;
				}
			}

		int productConnIndex = -1;
		for (int i=0; i<mProduct.getConnAtoms(productAtom); i++) {
			if (mProduct.getConnAtom(productAtom, i) == productRoot) {
				productConnIndex = i;
				break;
				}
			}

		for (int radius=SKELETON_PENALTY; radius<MAX_ENVIRONMENT_RADIUS; radius++) {
			if (mReactantSkelAtomEnv[reactantAtom][reactantConnIndex][radius] == null
			 || mSimilarityComparator.compare(mReactantSkelAtomEnv[reactantAtom][reactantConnIndex][radius], mProductSkelAtomEnv[productAtom][productConnIndex][radius]) != 0)
				return (radius - SKELETON_PENALTY) << SIMILARITY_SHIFT;
			}

		return MAX_SKELETON_SIMILARITY;
		}

	private int getPiAndHeteroBondCount(StereoMolecule mol, int atom) {
		int count = mol.getAtomPi(atom);
		for (int i=0; i<mol.getConnAtoms(atom); i++)
			if (mol.isElectronegative(mol.getConnAtom(atom, i)))
				count++;

		return count;
		}

	private byte[][][][] classifyNeighbourAtomEnvironment(StereoMolecule mol, boolean skeletonOnly) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
		StereoMolecule fragment = new StereoMolecule(mol.getAtoms(), mol.getBonds());

		byte[][][][] environment = new byte[mol.getAtoms()][MAX_ENVIRONMENT_RADIUS][][];

		int[] atomList = new int[mol.getAtoms()];
		int[] atomMap = new int[mol.getAtoms()];
		boolean[] atomMask = new boolean[mol.getAtoms()];
		for (int rootAtom=0; rootAtom<mol.getAtoms(); rootAtom++) {
			environment[rootAtom] = new byte[mol.getConnAtoms(rootAtom)][MAX_ENVIRONMENT_RADIUS][];

			if (skeletonOnly && mol.getAtomicNo(rootAtom) != 6)
				continue;

			for (int connIndex=0; connIndex<mol.getConnAtoms(rootAtom); connIndex++) {
				int fromAtom = mol.getConnAtom(rootAtom, connIndex);

				if (atomMask == null)
					atomMask = new boolean[mol.getAtoms()];
				else
					Arrays.fill(atomMask, false);

				int min = 1;
				int max = 2;

				// we need to mark the root atom, because otherwise close-by root atoms may end up with the same fragment
				mol.setAtomSelection(fromAtom, true);

				for (int sphere=0; sphere<MAX_ENVIRONMENT_RADIUS && max<mol.getAtoms(); sphere++) {
					if (sphere == 0) {
						atomList[0] = fromAtom;
						atomMask[fromAtom] = true;
						atomList[1] = rootAtom;
						atomMask[rootAtom] = true;
						}
					else {
						int newMax = max;
						for (int i=min; i<max; i++) {
							int atom = atomList[i];
							for (int j=0; j<mol.getConnAtoms(atom); j++) {
								int connAtom = mol.getConnAtom(atom, j);
								if (!atomMask[connAtom]
								 && (!skeletonOnly || mol.getAtomicNo(connAtom) == 6)) {
									atomMask[connAtom] = true;
									atomList[newMax++] = connAtom;
									}
								}
							}

						if (newMax == max) {
							// We consider small, but exactly matching groups to be of highest possible similarity
							if (!skeletonOnly)
								for (int i=sphere; i<MAX_ENVIRONMENT_RADIUS; i++)
									environment[rootAtom][connIndex][i] = environment[rootAtom][connIndex][i-1];
							break;
							}

						min = max;
						max = newMax;
						}

					if (sphere == 0) {
						environment[rootAtom][connIndex][sphere] = new byte[2];
						environment[rootAtom][connIndex][sphere][0] = (byte)mol.getAtomicNo(rootAtom);
						environment[rootAtom][connIndex][sphere][1] = (byte)mol.getAtomMass(rootAtom);
						}
					else {
						mol.copyMoleculeByAtoms(fragment, atomMask, true, atomMap);
						fragment.setAtomCharge(atomMap[fromAtom], 0);
						fragment.setAtomRadical(atomMap[fromAtom], 0);
						for (int atom=0; atom<mol.getAtoms(); atom++)
							if (atomMap[atom] != -1
							 && mol.getConnAtoms(atom) > fragment.getConnAtoms(atomMap[atom]))
								fragment.setAtomQueryFeature(atomMap[atom], Molecule.cAtomQFMoreNeighbours, true);
						if (skeletonOnly)
							for (int bond=0; bond<fragment.getBonds(); bond++)
								fragment.setBondType(bond, Molecule.cBondTypeSingle);
						environment[rootAtom][connIndex][sphere] = new Canonizer(fragment, Canonizer.ENCODE_ATOM_SELECTION).getIDCode().getBytes();
						}
					}

				mol.setAtomSelection(fromAtom, false);
				}
			}
		return environment;
		}

	/**
	 * Copies all reactant/product molecules into one StereoMolecule for easier processing
	 * Use getReactant() and getProduct() from outside to access these.
	 * @param rxn
	 */
	public void mergeReactantsAndProducts(Reaction rxn) {
		mReactant = new StereoMolecule();
		for (int i=0; i<rxn.getReactants(); i++) {
			StereoMolecule reactant = rxn.getReactant(i);
			reactant.ensureHelperArrays(Molecule.cHelperNeighbours);
			mReactant.addMolecule(reactant, reactant.getAtoms(), reactant.getBonds());
			}

		mProduct = new StereoMolecule();
		for (int i=0; i<rxn.getProducts(); i++) {
			StereoMolecule product = rxn.getProduct(i);
			product.ensureHelperArrays(Molecule.cHelperNeighbours);
			mProduct.addMolecule(product, product.getAtoms(), product.getBonds());
			}

		mReactant.normalizeAmbiguousBonds();
		mProduct.normalizeAmbiguousBonds();
		}

	public StereoMolecule getReactant() {
		return mReactant;
		}

	public StereoMolecule getProduct() {
		return mProduct;
		}

	/**
	 * Copies the generated mapping numbers from the merged reactant/product molecules into the original reaction.
	 * @param rxn
	 */
	public void copyMapNosToReaction(Reaction rxn, int[] reactantMapNo, int[] productMapNo, int graphMapNoCount) {
		int reactantIndex = 0;
		int reactantAtom = 0;
		for (int atom=0; atom<mReactant.getAtoms(); atom++) {
			StereoMolecule reactant = rxn.getReactant(reactantIndex);
			if (reactantAtom == reactant.getAtoms()) {
				reactantAtom = 0;
				reactant = rxn.getReactant(++reactantIndex);
				}
			reactant.setAtomMapNo(reactantAtom, reactantMapNo[atom], reactantMapNo[atom] <= graphMapNoCount);
			reactantAtom++;
			}

		int productIndex = 0;
		int productAtom = 0;
		for (int atom=0; atom<mProduct.getAtoms(); atom++) {
			StereoMolecule product = rxn.getProduct(productIndex);
			if (productAtom == product.getAtoms()) {
				productAtom = 0;
				product = rxn.getProduct(++productIndex);
				}
			product.setAtomMapNo(productAtom, productMapNo[atom], productMapNo[atom] <= graphMapNoCount);
			productAtom++;
			}
		}

	/**
	 * Maps a many atoms a possible starting from the given root atom pair.
	 * Mapping is done by incrementally adding the most similar neighbours to
	 * the currently mapped area provided the the respective bond orders are
	 * exactly matching.
	 * @param root
	 * @return number mapped atoms
	 */
	private int mapFromRootAtoms(RootAtomPair root) {
		int[] graphAtom = new int[mReactant.getAtoms()];
		int[] productAtom = new int[mProduct.getAtoms()];

		int[] graphParent = new int[mReactant.getAtoms()];
		int[] productParent = new int[mProduct.getAtoms()];

		graphAtom[0] = root.reactantAtom;
		graphParent[root.reactantAtom] = root.reactantAtom; // no -1 to simplify leavesRing()
		productAtom[0] = root.productAtom;
		productParent[root.productAtom] = root.productAtom; // no -1 to simplify leavesRing()

		if (mReactantMapNo[root.reactantAtom] == 0) {
			mMapNo++;
			mReactantMapNo[root.reactantAtom] = mMapNo;
			mProductMapNo[root.productAtom] = mMapNo;
			}

		int current = 0;
		int highest = 0;
		while (current <= highest) {
			// we take the best match of reaction & product neighbours as next graph members
			int reactantRoot = graphAtom[current];
			int productRoot = productAtom[current];
			int[][] match = new int[mReactant.getConnAtoms(reactantRoot)][mProduct.getConnAtoms(productRoot)];

			for (int i=0; i<mReactant.getConnAtoms(reactantRoot); i++) {
				int reactantCandidate = mReactant.getConnAtom(reactantRoot, i);
				if (mReactantMapNo[reactantCandidate] != 0)
					continue;

				int bondType = getBondType(mReactant, mReactant.getConnBond(reactantRoot, i));

				for (int j=0; j<mProduct.getConnAtoms(productRoot); j++) {
					int productCandidate = mProduct.getConnAtom(productRoot, j);
					if (mProductMapNo[productCandidate] == 0
					 && mReactant.getAtomicNo(reactantCandidate) == mProduct.getAtomicNo(productCandidate)) {
						int candidateBond = mProduct.getConnBond(productRoot, j);
						int skelSimilarity = getSkeletonSimilarity(reactantRoot, reactantCandidate, productRoot, productCandidate);
						if (bondType == getBondType(mProduct, candidateBond)
						 || skelSimilarity != 0) {
							if (passesBasicRules(reactantRoot, reactantCandidate, productRoot, productCandidate)) {
								int envSimilarity = getAtomSimilarity(reactantRoot, reactantCandidate, productRoot, productCandidate);
								int similarity = Math.max(skelSimilarity, envSimilarity);
								if (similarity != 0) {
									boolean isStereoMatch = matchesStereo(graphParent[reactantRoot], reactantRoot, reactantCandidate, productParent[productRoot], productRoot, productCandidate);
									if (passesSimilarityDependentRules(reactantRoot, reactantCandidate, productRoot, productCandidate, similarity, isStereoMatch)) {
										match[i][j] = similarity;
										// we prefer that atom, which restores the correct TH parity
										if (isStereoMatch)
											match[i][j] += STEREO_MATCH_PLUS;
										// we preferably endo/endo or exo/exo matches
										if (leavesRing(graphParent[reactantRoot], reactantRoot, reactantCandidate, mReactantRingMembership)
										 == leavesRing(productParent[productRoot], productRoot, productCandidate, mProductRingMembership))
											match[i][j] += ENDO_RING_PLUS;
										// we add a bonus if the pi-electron-plus-hetero-neighbour-count value matches at the first atom
										if (getPiAndHeteroBondCount(mReactant, reactantCandidate) == getPiAndHeteroBondCount(mProduct, productCandidate))
											match[i][j] += PI_AND_HETERO_PLUS;
/*
if (reactantRoot == 4) {
 System.out.print("graph:"); for (int k=0; k<current; k++) System.out.print(" "+graphAtom[k]+"("+productAtom[k]+")");
 System.out.println(" match:"+match[i][j]+ " r:"+reactantCandidate+" p:"+productCandidate);
}*/
										}
									}
								}
							}
						// if we have a changing bond order, but the rest of the substituent has maximal similarity
						else if (getAtomSimilarity(reactantRoot, reactantCandidate, productRoot, productCandidate) == MAX_ENVIRONMENT_SIMILARITY) {
//							match[i][j] = MAX_ENVIRONMENT_SIMILARITY - 1;
							match[i][j] = 3 << SIMILARITY_SHIFT;
							}
						}   // TODO may track multiple options and try all recursively; then end this method with best option mapped
					}
				}

			while (true) {
				int bestMatch = 0;
				int bestReactantAtom = -1;
				int bestProductAtom = -1;
				for (int i=0; i<match.length; i++) {
					int reactantCandidate = mReactant.getConnAtom(reactantRoot, i);
					if (mReactantMapNo[reactantCandidate] == 0) {
						for (int j=0; j<match[i].length; j++) {
							int productCandidate = mProduct.getConnAtom(productRoot, j);
							if (mProductMapNo[productCandidate] == 0) {
								if (bestMatch < match[i][j]) {
									bestMatch = match[i][j];
									bestReactantAtom = reactantCandidate;
									bestProductAtom = productCandidate;
									}
								}
							}
						}
					}

				if (bestMatch == 0)
					break;

				highest++;
				graphAtom[highest] = bestReactantAtom;
				graphParent[bestReactantAtom] = graphAtom[current];
				productAtom[highest] = bestProductAtom;
				productParent[bestProductAtom] = productAtom[current];
				mMapNo++;
				mReactantMapNo[bestReactantAtom] = mMapNo;
				mProductMapNo[bestProductAtom] = mMapNo;
				}

			current++;
			}

		return highest+1;
		}

	private int getBondType(StereoMolecule mol, int bond) {
		if (mol.isDelocalizedBond(bond))
			return 0;
		return mol.getBondTypeSimple(bond);
		}

	// Typically we accept the most similar neighbour pair as next atom pair in the mapping graph,
	// provided that bond order matches. In some cases, however, we apply some chemical knowledge,
	// to apply a veto, e.g. to prevent the match of the wrong oxygen atom in carboxylic acids
	// or esters.
	private boolean passesBasicRules(int reactantAtom, int reactantConn, int productAtom, int productConn) {
		// if the next carbon atom looses/gains more than one neighbour, it may be a different carbon, e.g. from rearrangement
		if (mReactant.getAtomicNo(reactantConn) == 6
		 && Math.abs(mReactant.getConnAtoms(reactantConn) - mProduct.getConnAtoms(productConn)) > 1)
			return false;

		// reductive ester cleavage
		if (mReactant.getAtomicNo(reactantConn) == 8
		 && mReactant.getConnAtoms(reactantConn) == 2
		 && connAtomsOfAtomicNo(mReactant, reactantAtom, 8) > connAtomsOfAtomicNo(mProduct, productAtom, 8))
			return false;
		if (mProduct.getAtomicNo(productConn) == 8
		 && mProduct.getConnAtoms(productConn) == 2
		 && connAtomsOfAtomicNo(mReactant, reactantAtom, 8) < connAtomsOfAtomicNo(mProduct, productAtom, 8))
			return false;

		// If the number of neighbours at an atom behind an oxygen changes, then the oxygen itself may in question
		if (mReactant.getAtomicNo(reactantConn) == 8
		 && mReactant.getConnAtoms(reactantConn) == 2
		 && mProduct.getConnAtoms(productConn) == 2) {
//			int reactantConnConn = mReactant.getConnAtom(reactantConn, mReactant.getConnAtom(reactantConn, 0) == reactantAtom ? 1 : 0);
//			int productConnConn = mProduct.getConnAtom(productConn, mProduct.getConnAtom(productConn, 0) == productAtom ? 1 : 0);
//			if (mReactant.getConnAtoms(reactantConnConn) != mProduct.getConnAtoms(productConnConn))
//				return false;

			int reactantConnIndex = mReactant.getConnAtom(reactantConn, 0) == reactantAtom ? 0 : 1;
			int productConnIndex = mProduct.getConnAtom(productConn, 0) == productAtom ? 0 : 1;
			if (mSimilarityComparator.compare(mReactantSkelAtomEnv[reactantConn][reactantConnIndex][3],
											  mProductSkelAtomEnv[productConn][productConnIndex][3]) != 0)
				return false;
			}

/* strangely this causes lots of false mappings
		// attached hydroxies may not be the same
		if (mReactant.getAtomPi(reactantAtom) == 0
		 && mReactant.getAtomicNo(reactantConn) == 8
		 && mReactant.getAtomPi(reactantConn) == 0
		 && mReactant.getConnAtoms(reactantConn) == 1
		 && mProduct.getConnAtoms(productConn) > 1)
			return false;
		if (mProduct.getAtomPi(productAtom) == 0
		 && mProduct.getAtomicNo(productConn) == 8
		 && mProduct.getAtomPi(productConn) == 0
		 && mProduct.getConnAtoms(productConn) == 1
		 && mReactant.getConnAtoms(reactantConn) > 1)
			return false;
*/
		// imine formation
		if (mReactant.getAtomicNo(reactantConn) == 7
		 && connAtomsOfAtomicNo(mReactant, reactantAtom, 7) < connAtomsOfAtomicNo(mProduct, productAtom, 7))
			return false;
		if (mProduct.getAtomicNo(productConn) == 7
		 && connAtomsOfAtomicNo(mReactant, reactantAtom, 7) > connAtomsOfAtomicNo(mProduct, productAtom, 7))
			return false;

		// if a 3-membered ring is formed or broken
		if (mReactant.getBondRingSize(mReactant.getBond(reactantAtom, reactantConn)) == 3
		  ^ mProduct.getBondRingSize(mProduct.getBond(productAtom, productConn)) == 3)
			return false;

		return true;
		}

	private boolean passesSimilarityDependentRules(int reactantAtom, int reactantConn, int productAtom, int productConn, int similarity, boolean isStereoMatch) {
		if (mReactant.getAtomicNo(reactantConn) == 8
		 && mReactant.getAtomPi(reactantConn) == 0
		 && hasOxo(mReactant, reactantAtom)
		 && hasOxo(mProduct, productAtom)
		 && similarity != MAX_ENVIRONMENT_SIMILARITY)
			return false;

		// Since carboxylates,sulfonates, etc are potential leaving groups, we dont follow a bond to oxygen,
		// if on the other side is a keto on the reactant side and similarity is not high
		if (mReactant.getAtomicNo(reactantConn) == 8
		 && mReactant.getConnAtoms(reactantConn) == 2
		 && hasOxo(mReactant, getNextNeighbour(mReactant, reactantAtom, reactantConn))
		 && similarity < (3 << SIMILARITY_SHIFT))
			return false;

		// A bond from oxygen to keto (C=O, S=O, C=N, etc) won't be followed,
		// unless the first target fragment shells are identical
		if (mReactant.getAtomicNo(reactantAtom) == 8
		 && (hasOxo(mReactant, reactantConn)
		  || hasOxo(mProduct, productConn))
		 && similarity < (2 << SIMILARITY_SHIFT))
			return false;

		if (!isStereoMatch
		 && (mReactant.getAtomicNo(reactantConn) != 6
		  || !hasNonCarbonNeighbour(mReactant, reactantAtom)))
			return false;

		// potential replacement of -O-R (R:H,C,any) at non-carbon atoms
		if (mReactant.getAtomicNo(reactantConn) == 8
		 && mReactant.getAtomicNo(reactantAtom) != 6
		 && similarity != MAX_ENVIRONMENT_SIMILARITY)
			return false;

		if (mReactant.getAtomicNo(reactantAtom) == 5
		 && mReactant.getAtomicNo(reactantConn) == 6
		 && similarity < (3 << SIMILARITY_SHIFT))
			return false;

		return true;
		}

	/**
	 * Returns true, if we have no stereo situation or if the stereo center (or EZ-bond) is retained
	 * by continuing the graph with the candidate.
	 * @param reactantParent
	 * @param reactantAtom
	 * @param reactantConn
	 * @param productParent
	 * @param productAtom
	 * @param productConn
	 * @return
	 */
	private boolean matchesStereo(int reactantParent, int reactantAtom, int reactantConn, int productParent, int productAtom, int productConn) {
		if (mReactant.getConnAtoms(reactantAtom) == 3
		 && (mReactant.getAtomParity(reactantAtom) == Molecule.cAtomParity1 || mReactant.getAtomParity(reactantAtom) == Molecule.cAtomParity2)
		 && mProduct.getConnAtoms(productAtom) == 3
		 && (mProduct.getAtomParity(productAtom) == Molecule.cAtomParity1 || mProduct.getAtomParity(productAtom) == Molecule.cAtomParity2)) {
			boolean reactantInversion = (reactantParent > reactantConn);
			int lastNeighbour = -1;
			for (int i=0; i<mReactant.getConnAtoms(reactantAtom); i++) {
				lastNeighbour = mReactant.getConnAtom(reactantAtom, i);
				if (lastNeighbour != reactantParent && lastNeighbour != reactantConn) {
					if ((lastNeighbour > reactantConn && lastNeighbour < reactantParent)
					 || (lastNeighbour < reactantConn && lastNeighbour > reactantParent))
						reactantInversion = !reactantInversion;
					break;
					}
				}
			boolean productInversion = (productParent > productConn);
			for (int i=0; i<mProduct.getConnAtoms(productAtom); i++) {
				lastNeighbour = mProduct.getConnAtom(productAtom, i);
				if (lastNeighbour != productParent && lastNeighbour != productConn) {
					if ((lastNeighbour > productConn && lastNeighbour < productParent)
					 || (lastNeighbour < productConn && lastNeighbour > productParent))
						productInversion = !productInversion;
					break;
					}
				}

			return (reactantInversion == productInversion)
				== (mReactant.getAtomParity(reactantAtom) == mProduct.getAtomParity(productAtom));
			}

		return true;
		}

	public static boolean hasOxo(StereoMolecule mol, int atom) {
		for (int i=0; i<mol.getConnAtoms(atom); i++)
			if (mol.getConnBondOrder(atom, i) == 2
			 && mol.getAtomicNo(mol.getConnAtom(atom, i)) > 6)
				return true;
		return false;
		}

	private boolean hasNonCarbonNeighbour(StereoMolecule mol, int atom) {
		for (int i=0; i<mol.getConnAtoms(atom); i++)
			if (mol.getAtomicNo(mol.getConnAtom(atom, i)) != 6)
				return true;
		return false;
		}

	private boolean isOrganicAcid(StereoMolecule mol, int atom) {
		boolean oxoLikeFound = false;
		boolean hydroxyFound = false;
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int connAtom = mol.getConnAtom(atom, i);
			if (mol.getAtomicNo(connAtom) == 8
			 && mol.getConnAtoms(connAtom) == 1
			 && mol.getConnBondOrder(atom, i) == 1) {
				hydroxyFound = true;
				}
			if (mol.isElectronegative(connAtom)
			 && mol.getConnBondOrder(atom, i) == 2) {
				oxoLikeFound = true;
				}
			}
		return oxoLikeFound && hydroxyFound;
		}

	private int getNextNeighbour(StereoMolecule mol, int rootAtom, int atom) {
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int connAtom = mol.getConnAtom(atom, i);
			if (connAtom != rootAtom)
				return connAtom;
			}
		return -1;
		}

	private int connAtomsOfAtomicNo(StereoMolecule mol, int atom, int atomicNo) {
		int count = 0;
		for (int i=0; i<mol.getConnAtoms(atom); i++)
			if (mol.getAtomicNo(mol.getConnAtom(atom, i)) == atomicNo)
				count++;

		return count;
		}
	}
