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

import com.actelion.research.chem.*;

import java.util.*;

public class CoordinateInventor {
	public static final int MODE_SKIP_DEFAULT_TEMPLATES = 1;
	public static final int MODE_REMOVE_HYDROGEN = 2;
	public static final int MODE_KEEP_MARKED_ATOM_COORDS = 4;
	public static final int MODE_PREFER_MARKED_ATOM_COORDS = 8;
	private static final int MODE_CONSIDER_MARKED_ATOMS = MODE_KEEP_MARKED_ATOM_COORDS | MODE_PREFER_MARKED_ATOM_COORDS;
	public static final int MODE_DEFAULT = MODE_REMOVE_HYDROGEN;

	private static final byte FLIP_AS_LAST_RESORT = 1;
	private static final byte FLIP_POSSIBLE = 2;
	private static final byte FLIP_PREFERRED = 3;
	private static final int  PREFERRED_FLIPS = 32;
	private static final int  POSSIBLE_FLIPS = 64;
	private static final int  LAST_RESORT_FLIPS = 128;
	private static final int  TOTAL_FLIPS = PREFERRED_FLIPS + POSSIBLE_FLIPS + LAST_RESORT_FLIPS;
	private static final float JOIN_DISTANCE_METAL_BONDS = 1.2f;
	private static final float JOIN_DISTANCE_CHARGED_ATOMS = 1.4f;
	public static final float JOIN_DISTANCE_UNCHARGED_FRAGMENTS = 1.6f;

	private static volatile List<InventorTemplate> sDefaultTemplateList;

	private StereoMolecule mMol;
	private long[]		mFFP;
	private Random		mRandom;
	private boolean[]	mAtomHandled;
	private boolean[]	mBondHandled;
	private boolean[]	mAtomIsPartOfCustomTemplate;
	private boolean     mAbsoluteOrientationTemplateFound;  // we just use the first matching template, which is set to define absolute orientation
	private int[]		mUnPairedCharge;
	private final int   mMode;
	private List<InventorFragment> mFragmentList;
	private List<InventorTemplate> mCustomTemplateList;

	private static synchronized void buildDefaultTemplateList() {
		if (sDefaultTemplateList == null)
			sDefaultTemplateList = new InventorDefaultTemplateList();
	}


	/**
	 * Creates an CoordinateInventor, which removes unneeded hydrogen atoms
	 * and creates new atom coordinates for all(!) atoms.
	 * This constructor creates a CoordinateInventor that uses templates
	 * of the InventorDefaultTemplateList to create 3D-projection derived coordinates for
	 * polycyclic structures from these templates (adamantanes, cubane, etc.).
	 */
	public CoordinateInventor () {
		this(MODE_DEFAULT);
		}


	/**
	 * Creates an CoordinateInventor, which removes unneeded hydrogens, if mode flags include
	 * MODE_REMOVE_HYDROGEN. If mode includes MODE_KEEP_MARKED_ATOM_COORDS, then marked atoms
	 * keep their coordinates. If mode includes MODE_PREFER_MARKED_ATOM_COORDS, then coordinates
	 * of marked atoms are changed only, if perfect coordinates are not possible without.
	 * Unless mode includes MODE_SKIP_DEFAULT_TEMPLATES, the CoordinateInventor uses templates
	 * of the InventorDefaultTemplateList to create 3D-projection derived coordinates for
	 * polycyclic structures from these templates (adamantanes, cubane, etc.).
	 * @param mode
	 */
	public CoordinateInventor (int mode) {
		mMode = mode;
		if ((mode & MODE_SKIP_DEFAULT_TEMPLATES) == 0 && sDefaultTemplateList == null)
			buildDefaultTemplateList();
		}


	public void setRandomSeed(long seed) {
		mRandom = new Random(seed);
		}


	/**
	 * A custom template list contains substructures with predefined atom coordinates.
	 * When such a list is provided, and if a molecules contains one of the list's substructures,
	 * then the matching atoms will receive the relative atom coordinates of the provided template,
	 * unless the substructure shares two or more atoms with another earlier found template or
	 * if mode is MODE_????_MARKED_ATOM_COORDS and the substructure match contains two or more
	 * non-marked (and therefore untouchable) atoms. If a template substructure contains an E- or Z-
	 * double bound, then the query feature 'match EZ-parity' should be set.
	 * @param templateList
	 */
	public void setCustomTemplateList(List<InventorTemplate> templateList) {
		mCustomTemplateList = templateList;
		for (InventorTemplate template:templateList)
			template.normalizeCoordinates();
		}

	/**
	 * Creates new atom 2D-coordinates for a molecule or a part of a molecule.
	 * Typically, the molecule has defined TH- and EZ-parities (even if unknown or none), which were not
	 * calculated, but taken from a SMILES or from an IDCode. In these cases setParitiesValid() should have
	 * been called to indicate that a parity calculation is not needed and even would destroy given parities.
	 * New coordinates will correctly reflect E/Z double bond parities, unless the double bond is in a small ring.
	 * If atom parities are available, this call is typically followed by calling mol.setStereoBondsFromParity();
	 * Unneeded explicit hydrogens are removed, if mode includes MODE_REMOVE_HYDROGEN.
	 * The relative orientation of all marked atoms is retained, if mode includes MODE_KEEP_MARKED_ATOM_COORDS.
	 * The relative orientation of all marked atoms is changed as last resort only, if mode includes MODE_PREFER_MARKED_ATOM_COORDS.
	 * If setTemplateList() was called, then any substructures matching any of the templates will be shown
	 * with the relative atom orientation provided with the template. If many molecule's coordinates are invented
	 * with templates, then you should also provide the molecules' fragment fingerprint to speed up template search
	 * using invent(mol, ffp).
	 * @param mol the molecule that gets new 2D coordinates in place
	 */
	public void invent(StereoMolecule mol) {
		invent(mol, null);
		}


	/**
	 * Creates new atom 2D-coordinates for a molecule or a part of a molecule.
	 * Typically, the molecule has defined TH- and EZ-parities (even if unknown or none), which were not
	 * calculated, but taken from a SMILES or from an IDCode. In these cases setParitiesValid() should have
	 * been called to indicate that a parity calculation is not needed and even would destroy given parities.
	 * New coordinates will correctly reflect E/Z double bond parities, unless the double bond is in a small ring.
	 * If atom parities are available, this call is typically followed by calling mol.setStereoBondsFromParity();
	 * Unneeded explicit hydrogens are removed, if mode includes MODE_REMOVE_HYDROGEN.
	 * The relative orientation of all marked atoms is retained, if mode includes MODE_KEEP_MARKED_ATOM_COORDS.
	 * The relative orientation of all marked atoms is changed as last resort only, if mode includes MODE_PREFER_MARKED_ATOM_COORDS.
	 * If setTemplateList() was called, then any substructures matching any of the templates will be shown
	 * with the relative atom orientation provided with the template. If many molecule's coordinates are invented
	 * with templates, then you should also provide the molecules' fragment fingerprint to speed up template search.
	 * @param mol the molecule that gets new 2D coordinates in place
	 * @parem ffp null or fragment fingerprint of the molecule, which is used (if available) for faster template location
	 */
	public void invent(StereoMolecule mol, long[] ffp) {
		boolean paritiesPresent = (mol.getHelperArrayStatus() & Molecule.cHelperParities) != 0;
		int parityState = mol.getHelperArrayStatus() & Molecule.cHelperBitsStereo;

		if (mRandom == null)
			mRandom = new Random();

		if ((mMode & MODE_REMOVE_HYDROGEN) != 0)
			mol.removeExplicitHydrogens(false, false);

		mMol = mol;
		mMol.ensureHelperArrays(Molecule.cHelperRings);

		mFFP = ffp;

		mFragmentList = new ArrayList<InventorFragment>() {
			@Override
			public boolean add(InventorFragment f) {
				for (InventorFragment ff:this)
					if (ff.equals(f))
						return false;
				return super.add(f);
				}
			};
		mAtomHandled = new boolean[mMol.getAllAtoms()];
		mBondHandled = new boolean[mMol.getAllBonds()];

		mUnPairedCharge = new int[mMol.getAllAtoms()];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			mUnPairedCharge[atom] = mMol.getAtomCharge(atom);

		if ((mMode & MODE_CONSIDER_MARKED_ATOMS) != 0) {
			locateMarkedFragments();
			}

		if (mCustomTemplateList != null)
			mAtomIsPartOfCustomTemplate = locateTemplateFragments(mCustomTemplateList, 512);

		if ((mMode & MODE_SKIP_DEFAULT_TEMPLATES) == 0 && sDefaultTemplateList != null)
			locateTemplateFragments(sDefaultTemplateList, 256);

		locateInitialFragments();
		joinOverlappingFragments();

		locateChainFragments();
		joinOverlappingFragments();

		for (InventorFragment f:mFragmentList)
			f.locateBonds();	// mGlobalBond & mGlobalToLocalAtom needed in correctChainEZParities() & optimizeFragments()

		correctChainEZParities();
		optimizeFragments();

		locateSingleAtoms();

		joinMetalBondedFragments();
		joinChargedFragments();

		// using one-by-one rotate and approximate strategy
		joinRemainingFragments();

		for (InventorFragment f : mFragmentList) {
			for (int j = 0; j<f.size(); j++) {
				mMol.setAtomX(f.mGlobalAtom[j], f.mAtomX[j]);
				mMol.setAtomY(f.mGlobalAtom[j], f.mAtomY[j]);
				mMol.setAtomZ(f.mGlobalAtom[j], 0.0);
			}
		}

		if (paritiesPresent) {
			mMol.setParitiesValid(parityState);
			mMol.setStereoBondsFromParity();
			}

		if (mAbsoluteOrientationTemplateFound)
			mMol.removeAtomMarkers();
		}


	public boolean[] getCustomTemplateAtomMask() {
		return mAtomIsPartOfCustomTemplate;
		}


	private boolean[] locateTemplateFragments(List<InventorTemplate> templateList, int priority) {
		boolean useFFP = (mFFP != null && !templateList.isEmpty() && templateList.get(0).getFFP() != null);

		SSSearcher searcher = null;
		SSSearcherWithIndex searcherWithIndex = null;
		if (useFFP) {
			searcherWithIndex = new SSSearcherWithIndex();
			searcherWithIndex.setMolecule(mMol, mFFP);
			}
		else {
			searcher = new SSSearcher();
			searcher.setMolecule(mMol);
			}

		boolean[] atomIsPartOfTemplate = new boolean[mMol.getAtoms()];

		for (InventorTemplate template: templateList) {
			ArrayList<int[]> matchList = null;
			StereoMolecule templateMol = template.getFragment();
			if (useFFP) {
				searcherWithIndex.setFragment(templateMol, template.getFFP());
				if (searcherWithIndex.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cDefaultMatchMode) != 0)
					matchList = searcherWithIndex.getGraphMatcher().getMatchList();
				}
			else {
				searcher.setFragment(templateMol);
				if (searcher.findFragmentInMolecule(SSSearcher.cCountModeOverlapping, SSSearcher.cDefaultMatchMode) != 0)
					matchList = searcher.getMatchList();
				}

			if (matchList != null) {
				for (int[] match:matchList) {
					int templateAtomCount = 0;
					for (int atom:match)
						if (atomIsPartOfTemplate[atom])
							templateAtomCount++;
					if (templateAtomCount <= 1) {
						// we just use the first matching template that is supposed to keep its absolute orientation
						boolean definesAbsoluteOrientation = template.keepAbsoluteOrientation();
						if (mAbsoluteOrientationTemplateFound)
							definesAbsoluteOrientation = false;
						else
							mAbsoluteOrientationTemplateFound = true;

						InventorFragment fragment = new InventorFragment(mMol, match.length, definesAbsoluteOrientation);

						for (int i=0; i<match.length; i++) {
							int atom = match[i];

							// translate templates, which need to retain absolute coordinates, into atom markers
							if (definesAbsoluteOrientation)
								mMol.setAtomMarker(atom, true);

							fragment.mPriority[i] = priority;
							fragment.mGlobalAtom[i] = atom;
							fragment.mAtomX[i] = template.getNormalizedAtomX(i);
							fragment.mAtomY[i] = template.getNormalizedAtomY(i);
							atomIsPartOfTemplate[atom] = true;
							mAtomHandled[atom] = true;
							}

						for (int b=0; b<templateMol.getBonds(); b++)
							mBondHandled[mMol.getBond(match[templateMol.getBondAtom(0, b)], match[templateMol.getBondAtom(1, b)])] = true;

						mFragmentList.add(fragment);
						}
					}
				}
			}
		return atomIsPartOfTemplate;
		}


	private void locateMarkedFragments() {
		int atomCount = 0;
		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			if (mMol.isMarkedAtom(atom))
				atomCount++;

		// there are no relative coordinates with less than 2 points
		if (atomCount < 2)
			return;

		int bondCount = 0;
		double avbl = 0;
		for (int bond=0; bond<mMol.getAllBonds(); bond++) {
			int atom1 = mMol.getBondAtom(0, bond);
			int atom2 = mMol.getBondAtom(1, bond);
			if (mMol.isMarkedAtom(atom1)
			 && mMol.isMarkedAtom(atom2)) {
				mBondHandled[bond] = true;
				mAtomHandled[atom1] = true;
				mAtomHandled[atom2] = true;
				avbl += mMol.getBondLength(bond);
				bondCount++;
				}
			}

		if (bondCount != 0 && avbl != 0.0)
			avbl /= bondCount;
		else {	// if we don't have an avbl from the marked bonds, we take it from all bonds
			avbl = mMol.getAverageBondLength();
			}

		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			if (mMol.isMarkedAtom(atom) && !mAtomHandled[atom])
				atomCount--;

		if (atomCount < 2)
			return;

		int[] fragmentNo = new int[mMol.getAllAtoms()];
		int coreFragmentCount = mMol.getFragmentNumbers(fragmentNo, true, true);

		int[] fragmentAtomCount = new int[coreFragmentCount];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			if (fragmentNo[atom] != -1)
				fragmentAtomCount[fragmentNo[atom]]++;

		InventorFragment[] fragment = new InventorFragment[coreFragmentCount];
		for (int f=0; f<coreFragmentCount; f++)
			fragment[f] = new InventorFragment(mMol, fragmentAtomCount[f], true);

		int[] atomIndex = new int[coreFragmentCount];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			int f = fragmentNo[atom];
			if (f != -1) {
				fragment[f].mPriority[atomIndex[f]] = 1024;
				fragment[f].mGlobalAtom[atomIndex[f]] = atom;
				fragment[f].mAtomX[atomIndex[f]] = mMol.getAtomX(atom) / avbl;
				fragment[f].mAtomY[atomIndex[f]] = mMol.getAtomY(atom) / avbl;
				atomIndex[f]++;
				}
			}

			// Find the largest core fragment and retain its orientation
			// by adding it first to the fragment list
		int maxFragment = -1;
		int maxFragmentAtoms = 0;
		for (int f=0; f<coreFragmentCount; f++) {
			if (maxFragmentAtoms < fragmentAtomCount[f]) {
				maxFragmentAtoms = fragmentAtomCount[f];
				maxFragment = f;
				}
			}

		mFragmentList.add(fragment[maxFragment]);
		for (int f=0; f<coreFragmentCount; f++)
			if (f != maxFragment)
				mFragmentList.add(fragment[f]);
		}


	private void locateInitialFragments() {
		// take every atom with more than 4 neighbours including first neighbour shell
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mMol.getAllConnAtoms(atom) > 4) {
				InventorFragment f = new InventorFragment(mMol, 1+mMol.getAllConnAtoms(atom), false);

				f.mAtomX[mMol.getAllConnAtoms(atom)] = 0.0;
				f.mAtomY[mMol.getAllConnAtoms(atom)] = 0.0;
				f.mPriority[mMol.getAllConnAtoms(atom)] = 32;
				f.mGlobalAtom[mMol.getAllConnAtoms(atom)] = atom;
				mAtomHandled[atom] = true;

				for (int i=0; i<mMol.getAllConnAtoms(atom); i++) {
					int connAtom = mMol.getConnAtom(atom, i);
					f.mAtomX[i] = Math.sin(Math.PI/3*i-Math.PI/3*2);
					f.mAtomY[i] = Math.cos(Math.PI/3*i-Math.PI/3*2);
					f.mPriority[i] = 32;
					f.mGlobalAtom[i] = connAtom;
					mAtomHandled[connAtom] = true;
					mBondHandled[mMol.getConnBond(atom, i)] = true;
					}

				mFragmentList.add(f);
				}
			}


			// take every small ring whose atoms are not a superset of another small ring
		RingCollection ringSet = mMol.getRingSet();
		for (int ringNo=0; ringNo<ringSet.getSize(); ringNo++) {
			int ringSize = ringSet.getRingSize(ringNo);
			int[] ringAtom = ringSet.getRingAtoms(ringNo);

				// skip rings that are entirely in the core fragment, if retainCore is true
			boolean skipRing = false;
			if ((mMode & MODE_CONSIDER_MARKED_ATOMS) != 0) {
				skipRing = true;
				for (int i=0; i<ringSize; i++) {
					if (!mMol.isMarkedAtom(ringAtom[i])) {
						skipRing = false;
						break;
						}
					}
				}

			if (!skipRing) {
				boolean isElementaryRing = false;
				for (int i=0; i<ringSize; i++) {
					if (mMol.getAtomRingSize(ringAtom[i]) == ringSize) {
						isElementaryRing = true;
						break;
						}
					}
				if (isElementaryRing) {
					int[] ringBond = ringSet.getRingBonds(ringNo);

					addRingFragment(ringAtom, ringBond);

					for (int i=0; i<ringSize; i++) {
						mAtomHandled[ringAtom[i]] = true;
						mBondHandled[ringBond[i]] = true;
						}
					}
				}
			}

			// take every large ring that has ring bonds that are not member of a fragment added already
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.isRingBond(bond) && !mBondHandled[bond]) {
				InventorChain theRing = getSmallestRingFromBond(bond);
				int[] ringAtom = theRing.getRingAtoms();
				int[] ringBond = theRing.getRingBonds();
				addRingFragment(ringAtom, ringBond);

				for (int i=0; i<theRing.getChainLength(); i++) {
					mAtomHandled[ringAtom[i]] = true;
					mBondHandled[ringBond[i]] = true;
					}
				}
			}

			// take every triple bond including first level attached atoms
		for (int bond=0; bond<mMol.getAllBonds(); bond++) {
			if (!mBondHandled[bond] && mMol.getBondOrder(bond) == 3) {
				int atom1 = mMol.getBondAtom(0, bond);
				int atom2 = mMol.getBondAtom(1, bond);
				int members = mMol.getAllConnAtoms(atom1) + mMol.getAllConnAtoms(atom2);
				if (members > 2) {
					InventorFragment f = new InventorFragment(mMol, members, false);
					int count = 0;
					for (int i=0; i<mMol.getAllConnAtoms(atom1); i++) {
						int connAtom = mMol.getConnAtom(atom1, i);
						if (connAtom != atom2) {
							f.mGlobalAtom[count++] = connAtom;
							mAtomHandled[connAtom] = true;
							mBondHandled[mMol.getConnBond(atom1, i)] = true;
							}
						}
					f.mGlobalAtom[count++] = atom1;
					f.mGlobalAtom[count++] = atom2;
					for (int i=0; i<mMol.getAllConnAtoms(atom2); i++) {
						int connAtom = mMol.getConnAtom(atom2, i);
						if (connAtom != atom1) {
							f.mGlobalAtom[count++] = connAtom;
							mAtomHandled[connAtom] = true;
							mBondHandled[mMol.getConnBond(atom2, i)] = true;
							}
						}
					for (int i=0; i<members; i++) {
						f.mAtomX[i] = i;
						f.mAtomY[i] = 0.0;
						f.mPriority[i] = 1;
						}
					mAtomHandled[atom1] = true;
					mAtomHandled[atom2] = true;
					mBondHandled[bond] = true;
					mFragmentList.add(f);
					}
				}
			}

			// take cumulated double bonds including first level single bonded atoms
		for (int bond=0; bond<mMol.getAllBonds(); bond++) {
			if (!mBondHandled[bond] && mMol.getBondOrder(bond) == 2) {
				int[] alleneAtom = new int[mMol.getAllAtoms()];
				for (int i=0; i<2; i++) {
					alleneAtom[0] = mMol.getBondAtom(i, bond);
					alleneAtom[1] = mMol.getBondAtom(1-i, bond);
					if (mMol.getAtomPi(alleneAtom[0]) == 1
					 && mMol.getAtomPi(alleneAtom[1]) == 2
					 && mMol.getAllConnAtoms(alleneAtom[1]) == 2) { // found start of cumulated double bonds
						mAtomHandled[alleneAtom[0]] = true;
						mAtomHandled[alleneAtom[1]] = true;
						mBondHandled[bond] = true;
						int last = 1;
						do {
							int nextIndex = (mMol.getConnAtom(alleneAtom[last], 0)
											 == alleneAtom[last-1]) ? 1 : 0;
							alleneAtom[last+1] = mMol.getConnAtom(alleneAtom[last], nextIndex);

								// stop at centers like C=Cr(Rn)=N
							if (mMol.getAtomPi(alleneAtom[last+1]) == 2
							 && mMol.getAllConnAtoms(alleneAtom[last+1]) > 2)
								break;

							mAtomHandled[alleneAtom[last+1]] = true;
							mBondHandled[mMol.getConnBond(alleneAtom[last], nextIndex)] = true;
							last++;
							} while (mMol.getAtomPi(alleneAtom[last]) == 2
								  && mMol.getAllConnAtoms(alleneAtom[last]) == 2);

						int members = mMol.getAllConnAtoms(alleneAtom[0])
									+ mMol.getAllConnAtoms(alleneAtom[last])
									+ last - 1;
						InventorFragment f = new InventorFragment(mMol, members, false);
						for (int j=0; j<=last; j++) {
							f.mAtomX[j] = j;
							f.mAtomY[j] = 0.0;
							f.mPriority[j] = 64;
							f.mGlobalAtom[j] = alleneAtom[j];
							}

						int current = last+1;
						boolean found = false;
						for (int j=0; j<mMol.getAllConnAtoms(alleneAtom[0]); j++) {
							int connAtom = mMol.getConnAtom(alleneAtom[0], j);
							if (connAtom != alleneAtom[1]) {
								f.mAtomX[current] = -0.5;
								f.mAtomY[current] = (found) ? Math.sin(Math.PI/3) : -Math.sin(Math.PI/3);
								f.mPriority[current] = 64;
								f.mGlobalAtom[current] = connAtom;
								current++;
								found = true;
								}
							}

						found = false;
						for (int j=0; j<mMol.getAllConnAtoms(alleneAtom[last]); j++) {
							int connAtom = mMol.getConnAtom(alleneAtom[last], j);
							if (connAtom != alleneAtom[last-1]) {
								f.mAtomX[current] = (double)last + 0.5;
								f.mAtomY[current] = (found) ? -Math.sin(Math.PI/3) : Math.sin(Math.PI/3);
								f.mPriority[current] = 64;
								f.mGlobalAtom[current] = connAtom;
								current++;
								found = true;
								}
							}

						mFragmentList.add(f);
						}
					}
				}
			}

			// predefine quarternary centers with exactly 2 not further subtituted substituents
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			if (mMol.getAllConnAtoms(atom) == 4) {
				int[] primaryConnAtom = new int[4];
				int[] primaryConnBond = new int[4];
				int primaryConns = 0;
				for (int i=0; i<4; i++) {
					primaryConnAtom[primaryConns] = mMol.getConnAtom(atom, i);
					primaryConnBond[primaryConns] = mMol.getConnBond(atom, i);
					if (mMol.getAllConnAtoms(primaryConnAtom[primaryConns]) == 1
					 && !mBondHandled[primaryConnBond[primaryConns]])
						primaryConns++;
					}

				if (primaryConns == 2) {
//					mAtomHandled[atom] = true;	don't break zig-zag of chains that are handled later
					InventorFragment f = new InventorFragment(mMol, 3, false);
					for (int i=0; i<2; i++) {
						mAtomHandled[primaryConnAtom[i]] = true;
						mBondHandled[primaryConnBond[i]] = true;
						f.mGlobalAtom[i] = primaryConnAtom[i];
						f.mPriority[i] = 32;
						}

					f.mAtomX[0] = -0.5;
					f.mAtomY[0] = 0.866;
					f.mAtomX[1] = 0.5;
					f.mAtomY[1] = 0.866;

					f.mAtomX[2] = 0.0;
					f.mAtomY[2] = 0.0;
					f.mPriority[2] = 32;
					f.mGlobalAtom[2] = atom;

					mFragmentList.add(f);
					}
				if (primaryConns == 3) {
					// if there is a single bond make sure that primaryConnBond[2] is one
					for (int i=0; i<2; i++) {
						if (mMol.getBondOrder(primaryConnBond[i]) == 1) {
							int temp = primaryConnAtom[i];
							primaryConnAtom[i] = primaryConnAtom[2];
							primaryConnAtom[2] = temp;
							temp = primaryConnBond[i];
							primaryConnBond[i] = primaryConnBond[2];
							primaryConnBond[2] = temp;
							}
						}

//					mAtomHandled[atom] = true;	don't break zig-zag of chains that are handled later
					InventorFragment f = new InventorFragment(mMol, 4, false);
					for (int i=0; i<3; i++) {
						mAtomHandled[primaryConnAtom[i]] = true;
						mBondHandled[primaryConnBond[i]] = true;
						f.mGlobalAtom[i] = primaryConnAtom[i];
						f.mPriority[i] = 32;
						}

					f.mAtomX[0] = -1.0;
					f.mAtomY[0] = 0.0;
					f.mAtomX[1] = 1.0;
					f.mAtomY[1] = 0.0;
					f.mAtomX[2] = 0.0;
					f.mAtomY[2] = 1.0;

					f.mAtomX[3] = 0.0;
					f.mAtomY[3] = 0.0;
					f.mPriority[3] = 32;
					f.mGlobalAtom[3] = atom;

					mFragmentList.add(f);
					}
				}
			}

/*	The current implementation does not (!!!) retain E/Z geometry of double bonds in a ring...
	Use the following to create E/Z fragments of ring double bonds with E/Z geometry reflecting
	coordinates. Retaining Reliably E/Z geometries, however, will need in addition:
	- consider these fragments coordinates when joining even if all fragments atom are already
	  part of the other fragment
	- a more capable joining algorithm that prevents E/Z inversions of double bonds by relocating
	  atoms that are part of the joint fragment

			// take stero defined ring-double-bonds including first level attached atoms
		for (int bond=0; bond<mMol.getAllBonds(); bond++) {
			int bondParity = mMol.getBondParity(bond);
			if (mMol.isRingBond(bond)
			 && (bondParity == Molecule.cBondParityE
			  || bondParity == Molecule.cBondParityZ)) {
			 	int[] bondAtom = new int[2];
				bondAtom[0] = mMol.getBondAtom(0, bond);
				bondAtom[1] = mMol.getBondAtom(1, bond);
				int members = mMol.getAllConnAtoms()[bondAtom[0]] + mMol.getAllConnAtoms()[bondAtom[1]];

				InventorFragment f = new InventorFragment(mMol, members);
				int count = 0;
				boolean[] secondAtomCounts = new boolean[2];
				for (int i=0; i<2; i++) {
					mAtomHandled[bondAtom[i]] = true;
					f.mAtomX[count] = (double)i-0.5;
					f.mAtomY[count] = 0.0;
					f.mPriority[count] = 128;
					f.mAtom[count++] = bondAtom[i];
					int neighbours = 0;
					for (int j=0; j<mMol.getAllConnAtoms()[bondAtom[i]]; j++) {
						int connAtom = mMol.getConnAtom(bondAtom[i], j);
						if (connAtom != bondAtom[1-i]) {
							if (neighbours == 1 && f.mAtom[count-1] > connAtom)
								secondAtomCounts[i] = true;
							f.mAtomX[count] = (i==0) ? -1.0 : 1.0;
							f.mAtomY[count] = (neighbours==0) ? -Math.sin(Math.PI/3) : Math.sin(Math.PI/3);
							f.mAtom[count++] = connAtom;
							mAtomHandled[connAtom] = true;
							mBondHandled[mMol.getConnBond(bondAtom[i], j)] = true;
							neighbours++;
							}
						}
					}

				if ((bondParity == Molecule.cBondParityE) ^ (secondAtomCounts[0] ^ secondAtomCounts[1]))
					for (int i=1; i<mMol.getAllConnAtoms()[mMol.getBondAtom(0, bond)]; i++)
						f.mAtomY[i] *= -1.0;

				mBondHandled[bond] = true;
				mFragmentList.addElement(f);
				}
			}	*/
		}


	private void locateChainFragments() {
		while (true) {
			InventorChain longestChain = null;

			for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
				int unhandledBonds = 0;
				for (int i=0; i<mMol.getAllConnAtoms(atom); i++)
					if (!mBondHandled[mMol.getConnBond(atom, i)])
						unhandledBonds++;

				if (unhandledBonds == 1) {
					InventorChain theChain = getLongestUnhandledChain(atom);
					if (longestChain == null
					 || theChain.getChainLength() > longestChain.getChainLength())
						longestChain = theChain;
					}
				}

			if (longestChain == null)
				break;

			InventorFragment f = new InventorFragment(mMol, longestChain.getChainLength(), false);
			for (int i=0; i<longestChain.getChainLength(); i++) {
				mAtomHandled[longestChain.mAtom[i]] = true;
				if (i < longestChain.getChainLength() - 1)
					mBondHandled[longestChain.mBond[i]] = true;
				f.mGlobalAtom[i] = longestChain.mAtom[i];
				f.mAtomX[i] = Math.cos(Math.PI / 6) * i;
				f.mAtomY[i] = ((i & 1) == 1) ? 0.0 : 0.5;
				f.mPriority[i] = 128 + longestChain.getChainLength();
				}
			mFragmentList.add(f);
			}
		}


	private void locateSingleAtoms() {
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			if (!mAtomHandled[atom] && mMol.getAllConnAtoms(atom) == 0) {
				InventorFragment f = new InventorFragment(mMol, 1, false);
				mAtomHandled[atom] = true;
				f.mGlobalAtom[0] = atom;
				f.mAtomX[0] = 0.0;
				f.mAtomY[0] = 0.0;
				f.mPriority[0] = 0;
				mFragmentList.add(f);
				}
			}
		}


	private void addRingFragment(int[] ringAtom, int[] ringBond) {
		int ringSize = ringAtom.length;
		InventorFragment f = new InventorFragment(mMol, ringSize, false);
		f.mAtomX[0] = 0.0;
		f.mAtomY[0] = 0.0;
		for (int i=0; i<ringSize; i++) {
			f.mPriority[i] = 128 - ringSize;
			f.mGlobalAtom[i] = ringAtom[i];
			}

		if (ringSize < 8)
			createRegularRingFragment(f);
		else  // create a large ring considering E-double bonds and other constraints
			createLargeRingFragment(f, ringAtom, ringBond);

		mFragmentList.add(f);
		}


	private void createRegularRingFragment(InventorFragment f) {
		double angleChange = Math.PI - (Math.PI * (f.size()-2))/f.size();
		for (int i=1; i<f.size(); i++) {
			f.mAtomX[i] = f.mAtomX[i-1] + Math.sin(angleChange*(i-1));
			f.mAtomY[i] = f.mAtomY[i-1] + Math.cos(angleChange*(i-1));
			}
		}


	private void createRegularRingFragment(InventorFragment f, int bondEConstraint, int bondZConstraint) {
		if (bondEConstraint == 0    // no E-bonds
		 || (bondEConstraint & bondZConstraint) != 0) { // contradictory constraints
			createRegularRingFragment(f);
			return;
			}

		int startIndex = -1;
		int startPriority = 0;
		int bitMinus2 = 1 << (f.size()-2);
		int bitMinus1 = 1 << (f.size()-1);
		int currentBit = 1;
		int bitPlus1 = 2;
		for (int i=0; i<f.size(); i++) {
			// if none of previous bond nor current bond are Z and at least one of them is E and the bond  before is not E
			if ((bondZConstraint & (bitMinus1 | currentBit)) == 0
			 && (bondEConstraint & (bitMinus1 | currentBit)) != 0
			 && (bondEConstraint & bitMinus2) == 0) {
				int priority = 0;
				if ((bondZConstraint & bitMinus2) != 0)
					priority += 4;
				if ((bondEConstraint & bitMinus1) != 0)
					priority += 2;
				if ((bondEConstraint & currentBit) != 0)
					priority += 1;
				if (startPriority < priority) {
					startPriority = priority;
					startIndex = i;
					}
				}

			bitMinus2 = bitMinus1;
			bitMinus1 = currentBit;
			currentBit = bitPlus1;
			bitPlus1 = 1 << (i+2 < f.size() ? i+2 : i+2-f.size());
			}

		if (startIndex == -1) {
			createRegularRingFragment(f);
			return;
			}

		int moveToCenter = 0;   // bit mask marking atoms to be moved towards ring center to satisfy E/Z constraints
		moveToCenter |= (1 << startIndex);
		int offset=2;
		while (offset<f.size()-1) {
			int index = (startIndex + offset < f.size()) ? startIndex+offset : startIndex+offset-f.size();
			bitMinus1 = 1 << (index==0 ? f.size()-1 : index-1);

			if ((bondZConstraint & bitMinus1) != 0) {
				offset++;
				continue;
				}

			currentBit = 1 << index;
			if ((bondEConstraint & bitMinus1) != 0) {
				if ((bondZConstraint & currentBit) != 0) {
					createRegularRingFragment(f);
					return;
					}
				moveToCenter |= currentBit;
				offset += 2;
				continue;
				}

			bitPlus1 = 1 << (index+1 < f.size() ? index+1 : index+1-f.size());
			if ((bondEConstraint & currentBit) != 0
			 && (bondZConstraint & bitPlus1) != 0) {
				moveToCenter |= currentBit;
				offset += 3;
				continue;
				}
			offset++;
			}

		if (moveToCenter == 0) {
			createRegularRingFragment(f);
			return;
			}

		double angleChange = Math.PI - (Math.PI * (f.size()-2))/f.size();
		for (int i=1; i<f.size(); i++) {
			f.mAtomX[i] = f.mAtomX[i-1] + Math.sin(angleChange*(i-1));
			f.mAtomY[i] = f.mAtomY[i-1] + Math.cos(angleChange*(i-1));
			}

		currentBit = 1;
		double shift = 2*Math.sin(angleChange/2);
		for (int i=0; i<f.size(); i++) {
			if ((moveToCenter & currentBit) != 0) {
				f.mAtomX[i] += shift * Math.cos(angleChange * ((double) i - 0.5));
				f.mAtomY[i] -= shift * Math.sin(angleChange * ((double) i - 0.5));
				}
			currentBit <<= 1;
			}
		}


	private void createLargeRingFragment(InventorFragment f, int[] ringAtom, int[] ringBond) {
		final int FIRST_RING_SIZE = 9;
		final int LAST_RING_SIZE = 25;
		final double[][] cAngleCorrection = { {20}, null, null, {0,10}, null, null, {-4,12},
				{0,0,-7.5}, null, null, null, null, {60.0/7, -60.0/7}, null, null, null, {-2.4} };
		final int[][] cBondZList = { // sequence of E/Z parities in rings (E=0, Z=1)
				{   // 9-membered ring
					0x00000092,  // 010010010 sym
				},
				{   // 10-membered ring
					0x00000273,  // 1001110011 sym
				},
				null,
				{   // 12-membered rings
					0x00000999,  // 100110011001 sym
					0x00000492   // 010010010010 sym
				},
				null,
				{   // 14-membered rings
					0x00000993,  // 00100110010011 sym
					0x000021C3,  // 10000111000011 sym
					0x000009D7   // 00100111010111 sym
				},
				{   // 15-membered rings
					0x00002492,  // 010010010010010 sym
					0x000039CE,  // 011100111001110 sym
				},
				{   // 16-membered rings
					0x00008649,  // 1000011001001001 sym
					0x80008759,  // 1000011101011001 asy
					0x00006666,  // 0110011001100110 sym
				},
				null,
				{   // 18-membered rings
					0x00009249,  // 001001001001001001 sym
					0x00021861,  // 100001100001100001 sym
					0x000175D7,  // 010111010111010111 sym
					0x00008643,  // 001000011001000011 sym
					0x000093B7,  // 001001001110110111 sym
					0x0000D66B,  // 001101011001101011 sym
					0x00020703,  // 100000011100000011 sym
					0x8002A753,  // 101010011101010011 asy
					0x0000D649,  // 001101011001001001 sym
					0x0000D759,  // 001101011101011001 sym
					0x80008753,  // 001000011101010011 asy
					0x80008717   // 001000011100010111 asy
				},
				null,
				{   // 20-membered rings
					0x00081909,  // 10000001100100001001 sym
					0x00081D6B,  // 10000001110101101011 sym
					0x000DB861,  // 11011011100001100001 sym
					0x00021849,  // 00100001100001001001 sym
					0x000A9959,  // 10101001100101011001 sym
					0x80081D49,  // 10000001110101001001 asy
					0x800819A3,  // 10000001100110100011 asy
					0x80084ED9,  // 10000100111011011001 asy
					0x80087475,  // 10000111010001110101 asy
					0x80087464,  // 10000111010001100100 asy
					0x800D19A9,  // 11010001100110101001 asy
					0x80086BA9,  // 10000110101110101001 asy
					0x800849A9,  // 10000100100110101001 asy
					0x80086B21   // 10000110101100100001 asy
				},
				{   // 21-membered rings
					0x000F5EBD,  // 011110101111010111101 sym
					0x00095263   // 010010101001010100101 sym
				},
				{   // 22-membered rings
					0x00084909,  // 0010000100100100001001 sym
					0x00021843,  // 0000100001100001000011 sym
					0x00206121,  // 1000000110000100100001 sym
					0x00081903,  // 0010000001100100000011 sym
					0x0021AC35,  // 1000011010110000110101 sym
					0x802A4D49,  // 1010100100110101001001 asy
					0x00035849,  // 0000110101100001001001 sym
					0x002B5909,  // 1010110101100100001001 sym
					0x00021953,  // 0000100001100101010011 sym
					0x80095909,  // 0010010101100100001001 asy
					0x80035959,  // 0000110101100101011001 asy
					0x00095D49,  // 0010010101110101001001 sym
					0x80206561,  // 1000000110010101100001 asy
					0x800D1909,  // 0011010001100100001001 asy
					0x000A9953,  // 0010101001100101010011 sym
					0x00257535,  // 1001010111010100110101 sym
					0x80207461,  // 1000000111010001100001 asy
					0x80021D13,  // 0000100001110100010011 asy
					0x800876C9,  // 0010000111011011001001 asy
					0x80086BA3,  // 0010000110101110100011 asy
					0x802B5D49,  // 1010110101110101001001 asy
					0x80081D43,  // 0010000001110101000011 asy
					0x800D192B,  // 0011010001100100101011 asy
					0x800D1D49,  // 0011010001110101001001 asy
					0x002B5D6B,  // 1010110101110101101011 sym
					0x001066D9,  // 1000000110011011011001 sym
					0x800D19A3,  // 0011010001100110100011 asy
					0x002AB953,  // 1010101011100101010011 sym
					0x802A1D43,  // 1010100001110101000011 asy
					0x00021D57,  // 0000100001110101010111 sym
					0x000D1C59,  // 0011010001110001011001 sym
					0x8021DB35,  // 1000011101101100110101 asy
					0x80229903,  // 1000101001100100000011 asy
					0x800D1D6B,  // 0011010001110101101011 asy
					0x802A76C9,  // 1010100111011011001001 asy
					0x800876EB,  // 0010000111011011101011 asy
					0x80369909,  // 1101101001100100001001 asy
					0x80347535,  // 1101000111010100110101 asy
					0x800A9917,  // 0010101001100100010111 asy
					0x0022EBA3,  // 1000101110101110100011 sym
					0x00084E97,  // 0010000100111010010111 sym
					0x00201C03,  // 1000000001110000000011 sym
					0x8008B917,  // 0010001011100100010111 asy
					0x802DD753,  // 1011011101011101010011 asy
					0x00377249,  // 1101110111001001001001 sym
					0x80095CB7,  // 0010010101110010110111 asy
					0x80081C17   // 0010000001110000010111 asy
				},
				null,
				{   // 24-membered rings
					0x00818181,  // 100000011000000110000001 sym
					0x002126D9,  // 001000010010011011011001 sym
					0x00204C03,  // 001000000100110000000011 sym
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //
					0x00000000,  //

					0x0086BB75   // 100001101011101101110101 sym
				},
				{   // 25-membered rings
					0x00D6B5AD   // 0110101101011010110101101 sym
				},
			};

//		if (f.size() >= FIRST_RING_SIZE)
//			printRingStats(f, ringAtom, ringBond);

		int maxBit = (1 << f.size());
		int bondEConstraint = 0;
		int bondZConstraint = 0;
		if (f.size() > 7) {
			for (int i=0; i<f.size(); i++) {
				int bondParity = getLargeRingBondParity(ringAtom, ringBond, i);
				if (bondParity == Molecule.cBondParityEor1)
					bondEConstraint += maxBit;
				else if (bondParity == Molecule.cBondParityZor2)
					bondZConstraint += maxBit;
				bondEConstraint >>>= 1;
				bondZConstraint >>>= 1;
				}
			}

		int ringIndex = f.size()-FIRST_RING_SIZE;
		if (f.size() >= FIRST_RING_SIZE && f.size() <= LAST_RING_SIZE && cBondZList[ringIndex] != null) {
			for (int zList=0; zList<cBondZList[ringIndex].length; zList++) {
				boolean isSymmetrical = ((0x80000000 & cBondZList[ringIndex][zList]) == 0);
				int bondZList = (0x7FFFFFFF & cBondZList[ringIndex][zList]);
				for (boolean inverted=false; !inverted; inverted = !inverted) {
					if (inverted) {
						if (isSymmetrical)
							break;

						int newBondZList = 0;
						for (int bit=1; bit!=maxBit; bit<<=1) {
							newBondZList <<= 1;
							if ((bondZList & bit) != 0)
								newBondZList |= 1;
							}
						bondZList = newBondZList;
						}
					for (int rotation=0; rotation<f.size(); rotation++) {
						if ((bondZList & bondEConstraint) == 0
						 && (~bondZList & bondZConstraint) == 0) {
							// constraints are satisfied with current E/Z sequence
							double bondAngle = 0.0;
							double correction = Math.PI / 180 * (cAngleCorrection[ringIndex] == null ?
									0 : cAngleCorrection[ringIndex][zList]);

							// we determine whether we have more right or left turns if we start with a right turn
							int rightTurns = 0;
							int tempZList = bondZList;
							boolean isRightTurn = true;
							for (int i=0; i<f.size(); i++) {
								if (isRightTurn)
									rightTurns++;
								if ((tempZList & 1) == 0) // is E-bond
									isRightTurn = !isRightTurn;
								tempZList >>>= 1;
								}

							// We choose the starting turn direction such that we end up
							// with more right turns and, thus, close the ring in clockwise direction.
							boolean wasRightTurn = (rightTurns > f.size() / 2); // create a ring that closes with right turns

							for (int i=1; i<f.size(); i++) {
								f.mAtomX[i] = f.mAtomX[i-1] + Math.sin(bondAngle);
								f.mAtomY[i] = f.mAtomY[i-1] + Math.cos(bondAngle);
								if ((bondZList & 1) == 0) // is E-bond
									wasRightTurn = !wasRightTurn;
								bondAngle += correction + (wasRightTurn ? Math.PI/3.0 : -Math.PI/3.0);
								bondZList >>>= 1;
								}
							return;
							}
						if ((bondZList & 1) != 0)
							bondZList |= maxBit;
						bondZList >>>= 1;
						}
					}
				}
			}

			// if not successful so far
		createRegularRingFragment(f, bondEConstraint, bondZConstraint);
		}


	private int getLargeRingBondParity(int[] ringAtom, int[] ringBond, int index) {
		int higherIndex = (index == ringAtom.length-1) ? 0 : index+1;    // second ringAtom index
		int lowerIndex = (index == 0) ? ringAtom.length-1 : index-1;
		int highestIndex = (higherIndex == ringAtom.length-1) ? 0 : higherIndex+1;
		if (mMol.getBondOrder(ringBond[index]) == 2) {
			int bondParity = mMol.getBondParity(ringBond[index]);
			if (bondParity == Molecule.cBondParityEor1
			 || bondParity == Molecule.cBondParityZor2) {
				// translate bond configuration from lowest atom index to ring members
				if (isLowestIndexNeighbour(ringAtom[lowerIndex], ringAtom[index], ringAtom[higherIndex])
				  ^ isLowestIndexNeighbour(ringAtom[highestIndex], ringAtom[higherIndex], ringAtom[index]))
						bondParity = (bondParity == Molecule.cBondParityEor1) ?
								Molecule.cBondParityZor2 : Molecule.cBondParityEor1;
				return bondParity;
				}
			}

		// If the bond is also a member of a small ring
		if (mMol.isSmallRingBond(ringBond[index])) {
			int sharedRing1 = mMol.getRingSet().getSharedRing(ringBond[lowerIndex], ringBond[index]);
			int sharedRing2 = mMol.getRingSet().getSharedRing(ringBond[higherIndex], ringBond[index]);

			// If more than one bond is shared by the large and the small ring,
			// then the first and the last shared bond are E and the ones between Z.
			if (sharedRing1 != -1 || sharedRing2 != -1)
				return sharedRing1 == sharedRing2 ? Molecule.cBondParityZor2 : Molecule.cBondParityEor1;

			// If only this one bond is shared, then we have a Z-configuration
			return Molecule.cBondParityZor2;
			}

		return Molecule.cBondParityNone;
		}


	/**
	 * @param atom
	 * @param rootAtom
	 * @param excludeAtom
	 * @return whether atom has the lowest atom index of rootAtom's neighbours not considering excludeAtom
	 */
	private boolean isLowestIndexNeighbour(int atom, int rootAtom, int excludeAtom) {
		for (int i=0; i<mMol.getConnAtoms(rootAtom); i++) {
			int connAtom = mMol.getConnAtom(rootAtom, i);
			if (connAtom != excludeAtom && connAtom < atom)
				return false;
			}
		return true;
		}


/*	private void printRingStats(InventorFragment f, int[] ringAtom, int[] ringBond) {
		int maxBit = (1 << f.size());
		int bondEConstraint = 0;
		int bondZConstraint = 0;
		for (int i = 0; i < f.size(); i++) {
			int bondParity = getLargeRingBondParity(ringAtom, ringBond, i);
			if (bondParity == Molecule.cBondParityEor1)
				bondEConstraint += maxBit;
			else if (bondParity == Molecule.cBondParityZor2)
				bondZConstraint += maxBit;
			bondEConstraint >>>= 1;
			bondZConstraint >>>= 1;
			}
		System.out.println("ringSize:"+f.size()
				+" E-list:"+Integer.toBinaryString(bondEConstraint)
				+" Z-list:"+Integer.toBinaryString(bondZConstraint));
		}*/


	/**
	 * Find the smallest ring of given bond
	 * @param bond
	 * @return
	 */
	private InventorChain getSmallestRingFromBond(int bond) {
		int atom1 = mMol.getBondAtom(0, bond);
		int atom2 = mMol.getBondAtom(1, bond);
		int[] graphAtom = new int[mMol.getAllAtoms()];
		int[] graphBond = new int[mMol.getAllAtoms()];
		int[] graphLevel = new int[mMol.getAllAtoms()];
		int[] graphParent = new int[mMol.getAllAtoms()];
		graphAtom[0] = atom1;
		graphAtom[1] = atom2;
		graphBond[1] = bond;
		graphLevel[atom1] = 1;
		graphLevel[atom2] = 2;
		graphParent[0] = -1;
		graphParent[1] = 0;
		int current = 1;
		int highest = 1;
		while (current <= highest) {
			for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
				int candidate = mMol.getConnAtom(graphAtom[current], i);
				if ((current > 1) && candidate == atom1) {
					InventorChain theRing = new InventorChain(graphLevel[graphAtom[current]]);
					graphBond[0] = mMol.getConnBond(graphAtom[current], i);
					int index = current;
					for (int j=0; j<theRing.getChainLength(); j++) {
						theRing.mAtom[j] = graphAtom[index];
						theRing.mBond[j] = graphBond[index];
						index = graphParent[index];
						}
					return theRing;
					}
				if (graphLevel[candidate] == 0 && mMol.isRingAtom(candidate)) {
					graphAtom[++highest] = candidate;
					graphBond[highest] = mMol.getConnBond(graphAtom[current], i);
					graphLevel[candidate] = graphLevel[graphAtom[current]] + 1;
					graphParent[highest] = current;
					}
				}
			current++;
			}
		return null;
		}


	private int getSmallestRingSize(int atom1, int atom2, int atom3) {
			// return size of smallest ring where atom1, atom2 and atom3 occurr in this order
		int[] graphAtom = new int[mMol.getAllAtoms()];
		int[] graphLevel = new int[mMol.getAllAtoms()];
		graphAtom[0] = atom2;
		graphAtom[1] = atom1;
		graphLevel[atom2] = 1;
		graphLevel[atom1] = 2;
		int current = 1;
		int highest = 1;
		while (current <= highest) {
			for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
				int candidate = mMol.getConnAtom(graphAtom[current], i);
				if (candidate == atom3)
					return 1 + graphLevel[graphAtom[current]];

				if (graphLevel[candidate] == 0 && mMol.isRingAtom(candidate)) {
					graphAtom[++highest] = candidate;
					graphLevel[candidate] = graphLevel[graphAtom[current]] + 1;
					}
				}
			current++;
			}
		return 0;
		}


	private void joinOverlappingFragments() {
		while (true) {
			int maxJoinPriority = 0;
			int maxCommonAtoms = 0;
			InventorFragment maxFragment1 = null;
			InventorFragment maxFragment2 = null;
			for (int i=1; i<mFragmentList.size(); i++) {
				InventorFragment f1 = mFragmentList.get(i);
				for (int j=0; j<i; j++) {
					InventorFragment f2 = mFragmentList.get(j);
					int commonAtom = 0;
					int commonAtoms = 0;
					int maxF1Priority = 0;
					int maxF2Priority = 0;
					for (int k=0; k<f1.size(); k++) {
						for (int l=0; l<f2.size(); l++) {
							if (f1.mGlobalAtom[k] == f2.mGlobalAtom[l]) {
								commonAtoms++;
								commonAtom = f1.mGlobalAtom[k];
								if (maxF1Priority < f1.mPriority[k])
									maxF1Priority = f1.mPriority[k];
								if (maxF2Priority < f2.mPriority[l])
									maxF2Priority = f2.mPriority[l];
								}
							}
						}

					if (commonAtoms > 0) {
						int handlePreferred = (commonAtoms == 1
											&& getConnAtoms(f1, commonAtom) == 1
											&& getConnAtoms(f2, commonAtom) == 1) ? 0 : 1;

						int joinPriority;
						if (maxF1Priority > maxF2Priority)
							joinPriority = (handlePreferred << 24)
										 + (maxF1Priority << 16)
										 + (maxF2Priority << 8)
										 + commonAtoms;
						else
							joinPriority = (handlePreferred << 24)
										 + (maxF2Priority << 16)
										 + (maxF1Priority << 8)
										 + commonAtoms;

						if (maxJoinPriority < joinPriority) {
							maxJoinPriority = joinPriority;
							maxCommonAtoms = commonAtoms;

							// retain coordinates of fragment with highest priority atom
							maxF1Priority = 0;
							maxF2Priority = 0;
							for (int k=0; k<f1.size(); k++)
								if (maxF1Priority < f1.mPriority[k])
									maxF1Priority = f1.mPriority[k];
							for (int k=0; k<f2.size(); k++)
								if (maxF2Priority < f2.mPriority[k])
									maxF2Priority = f2.mPriority[k];

							if (maxF1Priority > maxF2Priority) {
								maxFragment1 = f1;
								maxFragment2 = f2;
								}
							else {
								maxFragment1 = f2;
								maxFragment2 = f1;
								}
							}
						}
					}
				}

			if (maxJoinPriority == 0)
				break;

			if (maxCommonAtoms == maxFragment1.size())
				mFragmentList.remove(maxFragment1);
			else if (maxCommonAtoms == maxFragment2.size())
				mFragmentList.remove(maxFragment2);
			else
				joinFragments(maxFragment1, maxFragment2, maxCommonAtoms);
			}
		}


	private void joinFragments(InventorFragment f1, InventorFragment f2, int commonAtoms) {
		int[] commonAtom = new int[commonAtoms];
		int count = 0;
		for (int i = 0; i<f1.mGlobalAtom.length; i++)
			for (int j = 0; j<f2.mGlobalAtom.length; j++)
				if (f1.mGlobalAtom[i] == f2.mGlobalAtom[j])
					commonAtom[count++] = f1.mGlobalAtom[i];

		InventorFragment joinedFragment = (commonAtoms == 1) ?
				  getFusedFragment(f1, f2, commonAtom[0])
				: getFusedFragment(f1, f2, commonAtom, commonAtoms);

		updateFragmentList(f1, f2, joinedFragment);
		}


	private InventorFragment getFusedFragment(InventorFragment f1, InventorFragment f2, int commonAtom) {
		int index1 = f1.getLocalAtom(commonAtom);
		int index2 = f2.getLocalAtom(commonAtom);

		f2.translate(f1.mAtomX[index1]-f2.mAtomX[index2],
					 f1.mAtomY[index1]-f2.mAtomY[index2]);

		double angle1 = suggestNewBondAngle(f1, commonAtom);
		double angle2 = suggestNewBondAngle(f2, commonAtom);

		double angleInc = 0.0;
		if (getConnAtoms(f1, commonAtom) == 1
		 && getConnAtoms(f2, commonAtom) == 1)
			angleInc = Math.PI / 3;

		f2.rotate(f2.mAtomX[index2], f2.mAtomY[index2], angle1 - angle2 + angleInc + Math.PI);

		return getMergedFragment(f1, f2, 1);
		}


	private InventorFragment getFusedFragment(InventorFragment f1, InventorFragment f2,
											  int[] commonAtom, int commonAtoms) {
		int[] index1 = new int[commonAtoms];
		int[] index2 = new int[commonAtoms];
		for (int i=0; i<commonAtoms; i++) {
			index1[i] = f1.getLocalAtom(commonAtom[i]);
			index2[i] = f2.getLocalAtom(commonAtom[i]);
			}

		double meanX1 = 0.0;
		double meanY1 = 0.0;
		double meanX2 = 0.0;
		double meanY2 = 0.0;

		for (int i=0; i<commonAtoms; i++) {
			meanX1 += f1.mAtomX[index1[i]];
			meanY1 += f1.mAtomY[index1[i]];
			meanX2 += f2.mAtomX[index2[i]];
			meanY2 += f2.mAtomY[index2[i]];
			}
		meanX1 /= commonAtoms;
		meanY1 /= commonAtoms;
		meanX2 /= commonAtoms;
		meanY2 /= commonAtoms;
		f2.translate(meanX1 - meanX2, meanY1 - meanY2);

		InventorAngle[] f1Angle = new InventorAngle[commonAtoms];
		InventorAngle[] f2Angle = new InventorAngle[commonAtoms];
		InventorAngle[] angleDif = new InventorAngle[commonAtoms];
		InventorAngle[] angleDifFlip = new InventorAngle[commonAtoms];
		for (int i=0; i<commonAtoms; i++) {
			f1Angle[i] = new InventorAngle(meanX1, meanY1, f1.mAtomX[index1[i]], f1.mAtomY[index1[i]]);
			f2Angle[i] = new InventorAngle(meanX1, meanY1, f2.mAtomX[index2[i]], f2.mAtomY[index2[i]]);
			angleDif[i] = new InventorAngle(f1Angle[i].mAngle - f2Angle[i].mAngle,
											f1Angle[i].mLength * f2Angle[i].mLength);
			angleDifFlip[i] = new InventorAngle(f1Angle[i].mAngle + f2Angle[i].mAngle,
												f1Angle[i].mLength * f2Angle[i].mLength);
			}
		InventorAngle meanAngleDif = getMeanAngle(angleDif, commonAtoms);
		InventorAngle meanAngleDifFlip = getMeanAngle(angleDifFlip, commonAtoms);

		int neighbourCountF1 = 0;
		int neighbourCountF2 = 0;
		for (int i=0; i<commonAtoms; i++) {
			for (int j=0; j<mMol.getAllConnAtoms(commonAtom[i]); j++) {
				int connAtom = mMol.getConnAtom(commonAtom[i], j);

				if (f1.isMember(connAtom) && !f2.isMember(connAtom))
					neighbourCountF1++;

				if (!f1.isMember(connAtom) && f2.isMember(connAtom))
					neighbourCountF2++;
				}
			}

		InventorAngle[] f1NeighbourAngle = new InventorAngle[neighbourCountF1];
		InventorAngle[] f2NeighbourAngle = new InventorAngle[neighbourCountF2];
		InventorAngle[] f2NeighbourAngleFlip = new InventorAngle[neighbourCountF2];
		neighbourCountF1 = 0;
		neighbourCountF2 = 0;
		for (int i=0; i<commonAtoms; i++) {
			for (int j=0; j<mMol.getAllConnAtoms(commonAtom[i]); j++) {
				int connAtom = mMol.getConnAtom(commonAtom[i], j);

				if (f1.isMember(connAtom) && !f2.isMember(connAtom)) {
					int connIndex = f1.getLocalAtom(connAtom);
					f1NeighbourAngle[neighbourCountF1] = new InventorAngle(f1.mAtomX[index1[i]],
																		   f1.mAtomY[index1[i]],
																		   f1.mAtomX[connIndex],
																		   f1.mAtomY[connIndex]);
					neighbourCountF1++;
					}

				if (!f1.isMember(connAtom) && f2.isMember(connAtom)) {
					int connIndex = f2.getLocalAtom(connAtom);
					InventorAngle neighbourAngle = new InventorAngle(f2.mAtomX[index2[i]],
																	 f2.mAtomY[index2[i]],
																	 f2.mAtomX[connIndex],
																	 f2.mAtomY[connIndex]);
					f2NeighbourAngle[neighbourCountF2] = new InventorAngle(meanAngleDif.mAngle
																		   + neighbourAngle.mAngle,
																		   neighbourAngle.mLength);
					f2NeighbourAngleFlip[neighbourCountF2] = new InventorAngle(meanAngleDifFlip.mAngle
																			   - neighbourAngle.mAngle,
																			   neighbourAngle.mLength);
					neighbourCountF2++;
					}
				}
			}

		InventorAngle meanNeighbourAngleF1 = getMeanAngle(f1NeighbourAngle, neighbourCountF1);
		InventorAngle meanNeighbourAngleF2 = getMeanAngle(f2NeighbourAngle, neighbourCountF2);
		InventorAngle meanNeighbourAngleF2Flip = getMeanAngle(f2NeighbourAngleFlip, neighbourCountF2);

		if (Math.abs(getAngleDif(meanNeighbourAngleF1.mAngle, meanNeighbourAngleF2.mAngle))
		  > Math.abs(getAngleDif(meanNeighbourAngleF1.mAngle, meanNeighbourAngleF2Flip.mAngle))) {
			f2.rotate(meanX1, meanY1, meanAngleDif.mAngle);
			}
		else {
			f2.flip(meanX1, meanY1, 0.0);
			f2.rotate(meanX1, meanY1, meanAngleDifFlip.mAngle);
			}

		return getMergedFragment(f1, f2, commonAtoms);
		}


	private InventorFragment getMergedFragment(InventorFragment f1, InventorFragment f2, int commonAtoms) {
			// merges all atoms of two fragments into a new one retaining original coordinates
		InventorFragment f = new InventorFragment(mMol, f1.mGlobalAtom.length + f2.mGlobalAtom.length - commonAtoms, f1.mKeepMarkedAtoms | f2.mKeepMarkedAtoms);
		int count = 0;
		for (int i = 0; i<f1.mGlobalAtom.length; i++) {
			f.mGlobalAtom[count] = f1.mGlobalAtom[i];
			f.mPriority[count] = f1.mPriority[i];
			f.mAtomX[count] = f1.mAtomX[i];
			f.mAtomY[count++] = f1.mAtomY[i];
			}
		for (int i = 0; i<f2.mGlobalAtom.length; i++) {
			int index = f1.getLocalAtom(f2.mGlobalAtom[i]);
			if (index == -1) {
				f.mGlobalAtom[count] = f2.mGlobalAtom[i];
				f.mPriority[count] = f2.mPriority[i];
				f.mAtomX[count] = f2.mAtomX[i];
				f.mAtomY[count++] = f2.mAtomY[i];
				}
			else {
				if (f.mPriority[index] < f2.mPriority[i]) {
					f.mPriority[index] = f2.mPriority[i];
					f.mAtomX[index] = f2.mAtomX[i];
					f.mAtomY[index] = f2.mAtomY[i];
					}
				}
			}
		return f;
		}


	private InventorChain getLongestUnhandledChain(int atom) {
		int[] graphAtom = new int[mMol.getAllAtoms()];
		int[] graphBond = new int[mMol.getAllAtoms()];
		int[] graphLevel = new int[mMol.getAllAtoms()];
		int[] graphParent = new int[mMol.getAllAtoms()];
		graphAtom[0] = atom;
		graphLevel[atom] = 1;
		graphParent[0] = -1;
		int current = 0;
		int highest = 0;
		while (current <= highest) {
			if (current == 0 || !mAtomHandled[graphAtom[current]]) {
				for (int i=0; i<mMol.getAllConnAtoms(graphAtom[current]); i++) {
					int candidate = mMol.getConnAtom(graphAtom[current], i);
					int theBond = mMol.getConnBond(graphAtom[current], i);
					if (graphLevel[candidate] == 0 && !mBondHandled[theBond]) {
						graphAtom[++highest] = candidate;
						graphBond[highest] = theBond;
						graphLevel[candidate] = graphLevel[graphAtom[current]] + 1;
						graphParent[highest] = current;
						}
					}
				}
			if (current == highest) {
				InventorChain theChain = new InventorChain(graphLevel[graphAtom[current]]);
				int index = current;
				for (int j=0; j<theChain.getChainLength(); j++) {
					theChain.mAtom[j] = graphAtom[index];
					theChain.mBond[j] = graphBond[index];
					index = graphParent[index];
					}
				return theChain;
				}
			current++;
			}
		return null;
		}


	private double suggestNewBondAngle(InventorFragment f, int atom) {
		double[] connAngle = new double[mMol.getAllConnAtoms(atom)+1];
		int[] connAtom = new int[mMol.getAllConnAtoms(atom)+1];
		int[] connBond = new int[mMol.getAllConnAtoms(atom)+1];
		int rootIndex = f.getLocalAtom(atom);
		int connAngles = 0;
		for (int i=0; i<mMol.getAllConnAtoms(atom); i++) {
			connAtom[connAngles] = mMol.getConnAtom(atom, i);
			connBond[connAngles] = mMol.getConnBond(atom, i);
			int index = f.getLocalAtom(connAtom[connAngles]);
			if (index != -1)
				connAngle[connAngles++] = InventorAngle.getAngle(f.mAtomX[rootIndex],
																 f.mAtomY[rootIndex],
																 f.mAtomX[index],
																 f.mAtomY[index]);
			}

		if (connAngles == 1)
			return connAngle[0] + Math.PI;

		for (int i=connAngles-1; i>0; i--) {	// bubble sort
			for (int j=0; j<i; j++) {
				if (connAngle[j] > connAngle[j+1]) {
					double tempAngle = connAngle[j];
					connAngle[j] = connAngle[j+1];
					connAngle[j+1] = tempAngle;
					int tempAtom = connAtom[j];
					connAtom[j] = connAtom[j+1];
					connAtom[j+1] = tempAtom;
					int tempBond = connBond[j];
					connBond[j] = connBond[j+1];
					connBond[j+1] = tempBond;
					}
				}
			}

		connAngle[connAngles] = connAngle[0] + 2 * Math.PI;
		connAtom[connAngles] = connAtom[0];
		connBond[connAngles] = connBond[0];

		double maxAngleDif = -100.0;
		int maxIndex = 0;
		for (int i=0; i<connAngles; i++) {
			double angleDif = connAngle[i+1] - connAngle[i];
			if (connAngles > 2
			 && mMol.isRingBond(connBond[i])
			 && mMol.isRingBond(connBond[i+1])) {
				int ringSize = getSmallestRingSize(connAtom[i], atom, connAtom[i+1]);
				if (ringSize != 0)
					angleDif -= 100.0 - ringSize;
				}

			if (maxAngleDif < angleDif) {
				maxAngleDif = angleDif;
				maxIndex = i;
				}
			}

		return (connAngle[maxIndex] + connAngle[maxIndex+1]) / 2;
		}


	private double getAngleDif(double angle1, double angle2) {
		double angleDif = angle1 - angle2;
		while (angleDif < -Math.PI)
			angleDif += 2 * Math.PI;
		while (angleDif > Math.PI)
			angleDif -= 2 * Math.PI;
		return angleDif;
		}


	protected static InventorAngle getMeanAngle(InventorAngle[] angle, int noOfAngles) {
		// adds noOfAngles vectors of length=1 with angles angle[i]
		// and returns angle of sum-vector
		// length of sum-vector is criteria for deviation
		double sinSum = 0;
		double cosSum = 0;
		for (int i=0; i<noOfAngles; i++) {
			sinSum += angle[i].mLength * Math.sin(angle[i].mAngle);
			cosSum += angle[i].mLength * Math.cos(angle[i].mAngle);
			}

		double meanAngle;
		if (cosSum == 0)
			meanAngle = (sinSum > 0) ? Math.PI/2 : -Math.PI/2;
		else {
			meanAngle = Math.atan(sinSum/cosSum);
			if (cosSum < 0)
				meanAngle += Math.PI;
			}

		double length = Math.sqrt(sinSum * sinSum + cosSum * cosSum) / noOfAngles;

		return new InventorAngle(meanAngle, length);
		}


	private int getConnAtoms(InventorFragment f, int atom) {
		int connAtoms = 0;
		for (int i=0; i<mMol.getAllConnAtoms(atom); i++) {
			if (f.isMember(mMol.getConnAtom(atom, i)))
				connAtoms++;
			}
		return connAtoms;
		}


	private void correctChainEZParities() {
		for (InventorFragment f : mFragmentList) {
			for (int i = 0; i<f.mGlobalBond.length; i++) {
				int bond = f.mGlobalBond[i];

				if (mMol.getBondOrder(bond) == 2) {
					if (!mMol.isSmallRingBond(bond)
							&& (mMol.getBondParity(bond) == Molecule.cBondParityUnknown
							|| mMol.getBondParity(bond) == Molecule.cBondParityNone))
						mMol.setBondParityUnknownOrNone(bond);

					if (!mMol.isRingBond(bond)
							&& (mMol.getConnAtoms(mMol.getBondAtom(0, bond))>1)
							&& (mMol.getConnAtoms(mMol.getBondAtom(1, bond))>1)
							&& (mMol.getBondParity(bond) == Molecule.cBondParityEor1
							|| mMol.getBondParity(bond) == Molecule.cBondParityZor2)) {
						int[] minConnAtom = new int[2];
						int[] bondAtom = new int[2];
						for (int j = 0; j<2; j++) {
							minConnAtom[j] = mMol.getMaxAtoms();
							bondAtom[j] = mMol.getBondAtom(j, bond);
							for (int k = 0; k<mMol.getAllConnAtoms(bondAtom[j]); k++) {
								int connAtom = mMol.getConnAtom(bondAtom[j], k);
								if (connAtom != mMol.getBondAtom(1 - j, bond)
										&& minConnAtom[j]>connAtom)
									minConnAtom[j] = connAtom;
								}
							}

						double dbAngle = InventorAngle.getAngle(f.mAtomX[f.mGlobalToLocalAtom[bondAtom[0]]],
								f.mAtomY[f.mGlobalToLocalAtom[bondAtom[0]]],
								f.mAtomX[f.mGlobalToLocalAtom[bondAtom[1]]],
								f.mAtomY[f.mGlobalToLocalAtom[bondAtom[1]]]);
						double angle1 = InventorAngle.getAngle(f.mAtomX[f.mGlobalToLocalAtom[minConnAtom[0]]],
								f.mAtomY[f.mGlobalToLocalAtom[minConnAtom[0]]],
								f.mAtomX[f.mGlobalToLocalAtom[bondAtom[0]]],
								f.mAtomY[f.mGlobalToLocalAtom[bondAtom[0]]]);
						double angle2 = InventorAngle.getAngle(f.mAtomX[f.mGlobalToLocalAtom[bondAtom[1]]],
								f.mAtomY[f.mGlobalToLocalAtom[bondAtom[1]]],
								f.mAtomX[f.mGlobalToLocalAtom[minConnAtom[1]]],
								f.mAtomY[f.mGlobalToLocalAtom[minConnAtom[1]]]);

						if (((getAngleDif(dbAngle, angle1)<0)
								^ (getAngleDif(dbAngle, angle2)<0))
								^ (mMol.getBondParity(bond) == Molecule.cBondParityZor2)) {
							f.flipOneSide(bond);
							}
						}
					}
				}
			}
		}


	private void optimizeFragments() {
		int[] atomSymRank = calculateAtomSymmetries();

		byte[] bondFlipPriority = new byte[mMol.getAllBonds()];
		locateFlipBonds(bondFlipPriority, atomSymRank);
		for (int bond=0; bond<mMol.getAllBonds(); bond++)
			if (bondFlipPriority[bond] == FLIP_POSSIBLE
			 && (mMol.isRingAtom(mMol.getBondAtom(0, bond))
			  || mMol.isRingAtom(mMol.getBondAtom(1, bond))))
				bondFlipPriority[bond] = FLIP_PREFERRED;

		for (int fragmentNo=0; fragmentNo<mFragmentList.size(); fragmentNo++) {
			InventorFragment f = mFragmentList.get(fragmentNo);
			ArrayList<int[]> collisionList = f.getCollisionList();
			double minCollisionPanalty = f.getCollisionPanalty();
			InventorFragment minCollisionFragment = new InventorFragment(f);

			int lastBond = -1;
			for (int flip = 0; flip<TOTAL_FLIPS && !collisionList.isEmpty(); flip++) {
				int collisionNo = mRandom.nextInt(collisionList.size());
				int[] collidingAtom = collisionList.get(collisionNo);
				int[] bondSequence = getShortestConnection(collidingAtom[0], collidingAtom[1]);
				int[] availableBond = new int[bondSequence.length];
				int availableBonds = 0;

				if (flip < PREFERRED_FLIPS) {
					for (int i=1; i<bondSequence.length-1; i++)
						if (bondFlipPriority[bondSequence[i]] == FLIP_PREFERRED)
							availableBond[availableBonds++] = bondSequence[i];
					}
				else if (flip < PREFERRED_FLIPS + POSSIBLE_FLIPS) {
					for (int i=1; i<bondSequence.length-1; i++)
						if (bondFlipPriority[bondSequence[i]] >= FLIP_POSSIBLE)
							availableBond[availableBonds++] = bondSequence[i];
					}
				else {
					for (int i=1; i<bondSequence.length-1; i++)
						if (bondFlipPriority[bondSequence[i]] >= FLIP_AS_LAST_RESORT)
							availableBond[availableBonds++] = bondSequence[i];
					}

				if (availableBonds != 0) {
					int theBond = availableBond[0];
					if (availableBonds > 1) {
						do {	// don't rotate 2 times around same bond
							theBond = availableBond[mRandom.nextInt(availableBonds)];
							} while (theBond == lastBond);
						}

					if (theBond != lastBond) {
						lastBond = theBond;
		
						f.flipOneSide(theBond);

							// getCollisionList() is necessary to update collision panalty
						collisionList = f.getCollisionList();
						if (minCollisionPanalty > f.getCollisionPanalty()) {
							minCollisionPanalty = f.getCollisionPanalty();
							minCollisionFragment = new InventorFragment(f);
							}
						}
					}
				}

			mFragmentList.set(fragmentNo, minCollisionFragment);
				// finished optimization by rotating around single bonds

				// starting optimization by moving individual atoms
/*
double avbl = mMol.getAverageBondLength();
for (int i=0; i<f.mAtom.length; i++) {
f.mAtomX[i] = mMol.getAtomX(f.mAtom[i]) / avbl;
f.mAtomY[i] = mMol.getAtomY(f.mAtom[i]) / avbl;
}*/
			f = minCollisionFragment;
			int currentRank = 1;
			int nextAvailableRank;
			do {
				nextAvailableRank = 9999;
				for (int i=0; i<f.size(); i++) {
					int theRank = atomSymRank[f.mGlobalAtom[i]];
					if (theRank == currentRank)
						f.optimizeAtomCoordinates(i);
					else if (theRank > currentRank && theRank < nextAvailableRank)
						nextAvailableRank = theRank;
					}
				currentRank = nextAvailableRank;
				} while (nextAvailableRank != 9999);
			}
		}


	private int[] getShortestConnection(int atom1, int atom2) {
		int[] graphAtom = new int[mMol.getAllAtoms()];
		int[] graphBond = new int[mMol.getAllAtoms()];
		int[] graphLevel = new int[mMol.getAllAtoms()];
		int[] graphParent = new int[mMol.getAllAtoms()];
		graphAtom[0] = atom2;
		graphLevel[atom2] = 1;
		graphParent[0] = -1;
		int current = 0;
		int highest = 0;
		while (current <= highest) {
			for (int i=0; i<mMol.getAllConnAtoms(graphAtom[current]); i++) {
				int candidate = mMol.getConnAtom(graphAtom[current], i);
				int theBond = mMol.getConnBond(graphAtom[current], i);

				if (candidate == atom1) {
					int chainLength = graphLevel[graphAtom[current]];
					int[] bondSequence = new int[chainLength];
					bondSequence[0] = theBond;
					for (int j=1; j<chainLength; j++) {
						bondSequence[j] = graphBond[current];
						current = graphParent[current];
						}

					return bondSequence;
					}

				if (graphLevel[candidate] == 0) {
					graphAtom[++highest] = candidate;
					graphBond[highest] = theBond;
					graphLevel[candidate] = graphLevel[graphAtom[current]] + 1;
					graphParent[highest] = current;
					}
				}

			if (current == highest)
				return null;

			current++;
			}
		return null;
		}


	private void locateFlipBonds(byte[] bondFlipPriority, int[] atomSymRank) {
		for (int bond=0; bond<mMol.getAllBonds(); bond++) {
			int atom1 = mMol.getBondAtom(0, bond);
			int atom2 = mMol.getBondAtom(1, bond);

			if (mMol.isRingBond(bond)
			 || mMol.getBondOrder(bond) != 1
			 || mMol.getAllConnAtoms(atom1) == 1
			 || mMol.getAllConnAtoms(atom2) == 1)
				continue;

			if ((mMode & MODE_KEEP_MARKED_ATOM_COORDS) != 0
			 && mMol.isMarkedAtom(atom1)
			 && mMol.isMarkedAtom(atom2))
			 	continue;

			boolean oneBondEndIsSymmetric = false;
			for (int i=0; i<2; i++) {
				int bondAtom = mMol.getBondAtom(i, bond);
				if (mMol.getAllConnAtoms(bondAtom) > 2) {
					boolean symmetricEndFound = true;
					int connSymRank = -1;
					for (int j=0; j<mMol.getAllConnAtoms(bondAtom); j++) {
						int connAtom = mMol.getConnAtom(bondAtom, j);
						if (connAtom != mMol.getBondAtom(1-i, bond)) {
							if (connSymRank == -1)
								connSymRank = atomSymRank[connAtom];
							else if (connSymRank != atomSymRank[connAtom])
								symmetricEndFound = false;
							}
						}
					if (symmetricEndFound) {
						oneBondEndIsSymmetric = true;
						break;
						}
					}
				}
			if (!oneBondEndIsSymmetric) {
				if ((mMode & MODE_PREFER_MARKED_ATOM_COORDS) != 0
				 && mMol.isMarkedAtom(atom1)
				 && mMol.isMarkedAtom(atom2))
					bondFlipPriority[bond] = FLIP_AS_LAST_RESORT;
				else
					bondFlipPriority[bond] = FLIP_POSSIBLE;
				}
			}
		}


	private int[] calculateAtomSymmetries() {
		int atomBits = Canonizer.getNeededBits(mMol.getAtoms());
		int maxConnAtoms = 2;
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			maxConnAtoms = Math.max(maxConnAtoms, mMol.getAllConnAtoms(atom));
		int baseValueSize = (62 + 2 * atomBits + maxConnAtoms * (atomBits+1)) / 63;

		CanonizerBaseValue[] baseValue = new CanonizerBaseValue[mMol.getAllAtoms()];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			baseValue[atom] = new CanonizerBaseValue(baseValueSize);
			baseValue[atom].init(atom);
			}

		int[] symRank = new int[mMol.getAllAtoms()];

			// For calculating atom symmetries all atoms are considered the same.
			// Only the different connectivity makes an atom different from another.
			// However, double bonds of different parities must be taken into account.
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			int bondParity = mMol.getBondParity(bond);
			if (bondParity == Molecule.cBondParityEor1
			 || bondParity == Molecule.cBondParityZor2) {
				baseValue[mMol.getBondAtom(0, bond)].add(bondParity);
				baseValue[mMol.getBondAtom(1, bond)].add(bondParity);
				}
			}

		int oldNoOfRanks;
		int newNoOfRanks = consolidateRanks(baseValue, symRank);
		do {
			oldNoOfRanks = newNoOfRanks;
			calcNextBaseValues(baseValue, symRank, atomBits, maxConnAtoms);
			newNoOfRanks = consolidateRanks(baseValue, symRank);
			} while (oldNoOfRanks != newNoOfRanks);

		return symRank;
		}


	private void calcNextBaseValues(CanonizerBaseValue[] baseValue, int[] symRank, int atomBits, int maxConnAtoms) {
		int[] connRank = new int[maxConnAtoms];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
								// generate sorted list of ranks of neighbours
			for (int i=0; i<mMol.getAllConnAtoms(atom); i++) {
				int rank = symRank[mMol.getConnAtom(atom,i)];
				int j;
				for (j=0; j<i; j++)
					if (rank < connRank[j])
						break;
				for (int k=i; k>j; k--)
					connRank[k] = connRank[k-1];
				connRank[j] = rank;
				}

			int neighbours = mMol.getAllConnAtoms(atom);
			baseValue[atom].init(atom);
			baseValue[atom].add(atomBits, symRank[atom]);
			baseValue[atom].add((maxConnAtoms - neighbours)*(atomBits + 1), 0);
			for (int i=0; i<neighbours; i++)
				baseValue[atom].add(atomBits + 1, connRank[i]);
			}
		}


	private int consolidateRanks(CanonizerBaseValue[] baseValue, int[] symRank) {
		int rank = 0;
		Arrays.sort(baseValue);
		for (int i=0; i<baseValue.length; i++) {
			if (i == 0 || baseValue[i].compareTo(baseValue[i-1]) != 0)
				rank++;
			symRank[baseValue[i].getAtom()] = rank;
			}

		return rank;
		}


	private void joinMetalBondedFragments() {
		ArrayList<FragmentAssociation> associationList = createMetalBondAssociations();
		while (associationList != null) {
			FragmentAssociation association = getMaxPriorityAssociation(associationList);
			joinAssociatedFragments(association, JOIN_DISTANCE_METAL_BONDS);
			associationList = createMetalBondAssociations();
			}
		}

	private void joinChargedFragments() {
		FragmentAssociation association = createChargeAssociation();
		while (association != null) {
			joinAssociatedFragments(association, JOIN_DISTANCE_CHARGED_ATOMS);
			association = createChargeAssociation();
			}
		}

	private void joinRemainingFragments() {
		FragmentAssociation association = createDisconnectedAssociation();
		while (association != null) {
			joinAssociatedFragments(association, JOIN_DISTANCE_UNCHARGED_FRAGMENTS);
			association = createDisconnectedAssociation();
			}
		}

	private void joinAssociatedFragments(FragmentAssociation association, double minDistance) {
		association.arrange(minDistance, (mMode & MODE_CONSIDER_MARKED_ATOMS) != 0);
		InventorFragment mergedFragment = getMergedFragment(association.getFragment(0), association.getFragment(1), 0);
		updateFragmentList(association.getFragment(0), association.getFragment(1), mergedFragment);
		}

	private FragmentAssociation getMaxPriorityAssociation(ArrayList<FragmentAssociation> associationList) {
		int maxPriority = 0;
		FragmentAssociation maxFA = null;
		for (FragmentAssociation association:associationList) {
			if (maxPriority < association.getPriority()) {
				maxPriority = association.getPriority();
				maxFA = association;
				}
			}
		return maxFA;
		}

	private ArrayList<FragmentAssociation> createMetalBondAssociations() {
		ArrayList<FragmentAssociation> associationList = null;
		FragmentAssociation[][] fa = null;
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondType(bond) == Molecule.cBondTypeMetalLigand) {
				int atom1 = mMol.getBondAtom(0, bond);
				int atomIndex1 = -1;
				int f1 = 0;
				for (; f1<mFragmentList.size(); f1++) {
					atomIndex1 = mFragmentList.get(f1).getLocalAtom(atom1);
					if (atomIndex1 != -1)
						break;
					}
				int atom2 = mMol.getBondAtom(1, bond);
				int atomIndex2 = -1;
				int f2 = 0;
				for (; f2<mFragmentList.size(); f2++) {
					atomIndex2 = mFragmentList.get(f2).getLocalAtom(atom2);
					if (atomIndex2 != -1)
						break;
					}

				if (f1 != f2) {
					if (f1 > f2) {	// f1 must be the smaller one (and potentially the core fragment with fixed coords)
						int temp = f1;
						f1 = f2;
						f2 = temp;
						temp = atomIndex1;
						atomIndex1 = atomIndex2;
						atomIndex2 = temp;
						}
					if (fa == null)
						fa = new FragmentAssociation[mFragmentList.size()][];
					if (fa[f2] == null)
						fa[f2] = new FragmentAssociation[f2];
					if (fa[f2][f1] != null)
						fa[f2][f1].add(atomIndex1, atomIndex2);
					else {
						fa[f2][f1] = new FragmentAssociation(mFragmentList.get(f1), mFragmentList.get(f2), atomIndex1, atomIndex2);
						if (associationList == null)
							associationList = new ArrayList<FragmentAssociation>();
						associationList.add(fa[f2][f1]);
						}
					}
				}
			}

		return associationList;
		}

	private FragmentAssociation createChargeAssociation() {
		ArrayList<InventorCharge> negChargeList = new ArrayList<InventorCharge>();
		ArrayList<InventorCharge> posChargeList = new ArrayList<InventorCharge>();
		ArrayList<InventorCharge> chargeList = new ArrayList<InventorCharge>();

		for (InventorFragment f:mFragmentList) {
			int fragmentCharge = 0;
			chargeList.clear();
			for (int i=0; i<f.size(); i++) {
				int atom = f.getGlobalAtom(i);
				int charge = mUnPairedCharge[atom];
				if (charge != 0) {
					chargeList.add(new InventorCharge(f, i, charge));
					fragmentCharge += charge;
					}
				}
			if (fragmentCharge != 0) {
				//@Override Annotation incompatible with 1.5
				chargeList.sort((o1, o2) -> {
					int c1 = Math.abs(o1.charge);
					int c2 = Math.abs(o2.charge);
					return Integer.compare(c1, c2);
				});
				for (InventorCharge ic:chargeList) {
					if (fragmentCharge * ic.charge > 0) {	// charges have same sign
						int charge = (Math.abs(fragmentCharge) >= Math.abs(ic.charge)) ? ic.charge : fragmentCharge;
						fragmentCharge -= charge;
						(charge < 0 ? negChargeList : posChargeList).add(new InventorCharge(f, ic.atom, charge));
						if (fragmentCharge == 0)
							break;
						}
					}
				}
			}

		if (negChargeList.isEmpty() || posChargeList.isEmpty())
			return null;

		// with positive charges we have the large fragments first
		posChargeList.sort(new Comparator<InventorCharge>() {
			//@Override Annotation incompatible with 1.5
			public int compare(InventorCharge o1, InventorCharge o2) {
				int c1 = o1.fragment.size();
				int c2 = o1.fragment.size();
				return Integer.compare(c2, c1);
			}
		});

		// with negative charges we have the small o1.fragments first
		negChargeList.sort(new Comparator<InventorCharge>() {
			//@Override Annotation incompatible with 1.5
			public int compare(InventorCharge o1, InventorCharge o2) {
				int c1 = o1.fragment.size();
				int c2 = o1.fragment.size();
				return Integer.compare(c1, c2);
			}
		});

		// we combine preferably small with large
		for (InventorCharge pc:posChargeList) {
			for (InventorCharge nc : negChargeList) {
				if (pc.charge == -nc.charge) {
					mUnPairedCharge[pc.fragment.getGlobalAtom(pc.atom)] -= pc.charge;
					mUnPairedCharge[nc.fragment.getGlobalAtom(nc.atom)] -= nc.charge;
					return new FragmentAssociation(pc.fragment, nc.fragment, pc.atom, nc.atom);
					}
				}
			}

		for (InventorCharge pc:posChargeList) {
			for (InventorCharge nc : negChargeList) {
				if (pc.charge > -nc.charge) {
					mUnPairedCharge[pc.fragment.getGlobalAtom(pc.atom)] += nc.charge;
					mUnPairedCharge[nc.fragment.getGlobalAtom(nc.atom)] -= nc.charge;
					return new FragmentAssociation(pc.fragment, nc.fragment, pc.atom, nc.atom);
					}
				}
			}

		for (InventorCharge pc:posChargeList) {
			for (InventorCharge nc : negChargeList) {
				if (pc.charge < -nc.charge) {
					mUnPairedCharge[pc.fragment.getGlobalAtom(pc.atom)] -= pc.charge;
					mUnPairedCharge[nc.fragment.getGlobalAtom(nc.atom)] += pc.charge;
					return new FragmentAssociation(pc.fragment, nc.fragment, pc.atom, nc.atom);
					}
				}
			}

		return null;
		}

	private FragmentAssociation createDisconnectedAssociation() {
		if (mFragmentList.size() < 2)
			return null;

		return new FragmentAssociation(mFragmentList.get(0), mFragmentList.get(1));
		}

	/* ingenious solution to whether we have a cyclic relationship between fragments and metal bond associations
	private boolean isCyclic(FragmentAssociation[][] fa) {
		for (int i=1; i<fa.length; i++) {
			for (int j=0; j<i; j++) {
				if (fa[i][j] != null) {
					for (int ii=i+1; ii<fa.length; ii++) {
						if (fa[ii][j] != null) {
							if (fa[ii][i] != null)
								return true;
							fa[ii][i] = fa[i][j];	// indicate the higher order connection
							}
						}
					for (int jj=j+1; jj<i; jj++) {
						if (fa[i][jj] != null) {
							if (fa[jj][j] != null)
								return true;
							fa[jj][j] = fa[i][j];	// indicate the higher order connection
							}
						}
					}
				}
			}
		return false;
		}*/

	private void updateFragmentList(InventorFragment fOld1, InventorFragment fOld2, InventorFragment fJoined) {
		int index = Math.min(mFragmentList.indexOf(fOld1), mFragmentList.indexOf(fOld2));
		mFragmentList.add(index, fJoined);
		mFragmentList.remove(fOld1);
		mFragmentList.remove(fOld2);
		}
	}
