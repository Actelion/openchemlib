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
*/

package com.actelion.research.chem.reaction;

import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;

import com.actelion.research.chem.DrawingObjectList;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

public class Reaction implements java.io.Serializable {
	static final long serialVersionUID = 0x2006CAFE;

	private ArrayList<StereoMolecule> mReactant;
	private ArrayList<StereoMolecule> mProduct;
	private ArrayList<StereoMolecule> mCatalyst;
	private DrawingObjectList mDrawingObjectList;
	private String mName;
	private int mMaxMapNo;
	private boolean mIsFragment;	// if there are molecules, then there fragment status takes precedence over this flag

	public Reaction() {
		mReactant = new ArrayList<StereoMolecule>();
		mProduct = new ArrayList<StereoMolecule>();
		mCatalyst = new ArrayList<StereoMolecule>();
		mMaxMapNo = -1;
		mIsFragment = false;
		}

	public Reaction(String name) {
		this();
		mName = name;
		}

	public void clear() {
		mReactant.clear();
		mProduct.clear();
		mCatalyst.clear();
		mDrawingObjectList = null;
		mMaxMapNo = -1;
		}

	public void removeCatalysts() {
		mCatalyst.clear();
		}

	public void removeDrawingObjects() {
		mDrawingObjectList = null;
		}

	/**
	 * @return true, if neither of reactants, products, or catalysts contain any atoms
	 */
	public boolean isEmpty() {
		for (StereoMolecule mol:mReactant)
			if (mol.getAllAtoms() != 0)
				return false;

		for (StereoMolecule mol:mProduct)
			if (mol.getAllAtoms() != 0)
				return false;

		for (StereoMolecule mol:mCatalyst)
			if (mol.getAllAtoms() != 0)
				return false;

		return true;
		}

	/**
	 * Sets all reactants and products of this reaction to the given fragment state, i.e. whether they are considered
	 * molecules with all valences satisfied with hydrogens or whether they are substructure fragments that can have query features.
	 * @param f
	 */
	public void setFragment(boolean f) {
		mIsFragment = f;

		for (StereoMolecule mol:mReactant)
			mol.setFragment(f);

		for (StereoMolecule mol:mProduct)
			mol.setFragment(f);
		}

	/**
	 * The naming of this method is in analogy to the corresponding method of the Molecule class.
	 * Returns whether at least one of the reactants or products is marked as substructure fragment.
	 * Only if there are no molecules, then the reaction's explicit fragment status is returned.
	 * @return fragment status of molecules or reaction
	 */
	public boolean isFragment() {
		return getMolecules() == 0 ? mIsFragment : determineFragment();
		}

	/**
	 * @return whether at least one of the reactants or products is marked as substructure fragment.
	 */
	private boolean determineFragment() {
		for (StereoMolecule mol:mReactant)
			if (mol.isFragment())
				return true;
		for (StereoMolecule mol:mProduct)
			if (mol.isFragment())
				return true;

		return false;
		}

	public Reaction(Reaction rxn) {
		this();
		int r = (rxn == null) ? 0 : (rxn.mReactant == null ? 0 : rxn.mReactant.size());
		int p = (rxn == null) ? 0 : (rxn.mProduct == null ? 0 : rxn.mProduct.size());
		int c = (rxn == null) ? 0 : (rxn.mCatalyst == null ? 0 : rxn.mCatalyst.size());
		for (int i = 0; i < r; i++)
			mReactant.add(new StereoMolecule(rxn.getReactant(i)));
		for (int i = 0; i < p; i++)
			mProduct.add(new StereoMolecule(rxn.getProduct(i)));
		for (int i = 0; i < c; i++)
			mCatalyst.add(new StereoMolecule(rxn.getCatalyst(i)));
		mDrawingObjectList = new DrawingObjectList(rxn.getDrawingObjects());
		if (rxn.mName != null)
			mName = rxn.mName;
		mIsFragment = rxn.mIsFragment;
		}

	public Reaction(StereoMolecule[] mol, int reactantCount) {
		this();
		if (mol != null) {
			for (int i = 0; i < reactantCount; i++)
				mReactant.add(mol[i]);
			for (int i = reactantCount; i < mol.length; i++)
				mProduct.add(mol[i]);
			}
		mIsFragment = determineFragment();
		}

	public StereoMolecule getReactant(int no) {
		return mReactant.get(no);
		}

	public int getReactants() {
		return mReactant.size();
		}

	public StereoMolecule getProduct(int no) {
		return mProduct.get(no);
		}

	public int getProducts() {
		return mProduct.size();
		}

	public StereoMolecule getCatalyst(int no) {
		return mCatalyst.get(no);
	}

	public int getCatalysts() {
		return mCatalyst.size();
	}

	/**
	 * @return count of reactants and products
	 */
	public int getMolecules() {
		return mReactant.size() + mProduct.size();
		}

	public StereoMolecule getMolecule(int no) {
		return (no < mReactant.size()) ?
			mReactant.get(no)
			: mProduct.get(no - mReactant.size());
		}

	public void addReactant(StereoMolecule reactant) {
		mReactant.add(reactant);
		mMaxMapNo = -1;
		}

	public void addReactant(StereoMolecule reactant, int position) {
		mReactant.add(position, reactant);
		mMaxMapNo = -1;
		}

	public void addProduct(StereoMolecule product) {
		mProduct.add(product);
		mMaxMapNo = -1;
		}

	public void addProduct(StereoMolecule product, int position) {
		mProduct.add(position, product);
		mMaxMapNo = -1;
		}

	public void addCatalyst(StereoMolecule catalyst) {
		mCatalyst.add(catalyst);
	}

	public void addCatalyst(StereoMolecule catalyst, int position) {
		mCatalyst.add(position, catalyst);
	}

	public String getName() {
		return (mName == null) ? "Unknown Reaction" : mName;
		}

	public void setName(String name) {
		mName = name;
		}

	public DrawingObjectList getDrawingObjects() {
		return mDrawingObjectList;
		}

	public void setDrawingObjects(DrawingObjectList l) {
		mDrawingObjectList = l;
		}

	public double getAverageBondLength() {
		int bondCount = 0;
		double avbl = 0.0;
		for (int i=0; i<getMolecules(); i++) {
			StereoMolecule mol = getMolecule(i);
			if (mol.getAllBonds() != 0) {
				bondCount += mol.getAllBonds();
				avbl += mol.getAverageBondLength() * mol.getAllBonds();
				}
			}

		return bondCount == 0 ? Molecule.getDefaultAverageBondLength() : avbl / bondCount;
		}

	/**
	 * @return whether the molecule's atom coordinate bounds touch or overlap
	 */
	public boolean isReactionLayoutRequired() {
		if (getMolecules() <= 1)
			return false;

		double avbl = getAverageBondLength();

		Rectangle2D.Double[] r = new Rectangle2D.Double[getMolecules()];

		for (int i=0; i<getMolecules(); i++) {
			r[i] = getMolecule(i).getBounds(null);
			if (r[i] != null) {
				for (int j=0; j<i; j++) {
					if (r[j] != null) {
						// make new layout, if molecule bounds overlap
						if (r[i].x + r[i].width >= r[j].x && r[i].x <= r[j].x + r[j].width)
							return true;
						if (r[i].y + r[i].height >= r[j].y && r[i].y <= r[j].y + r[j].height)
							return true;
						}
					}

				// make new layout, molecule bounds are unreasonably far from each other
				if (i != 0 && r[i-1] != null) {
					if (r[i].x - r[i-1].x - r[i].width > 5 * avbl)
						return true;
					if (r[i].y - r[i-1].y - r[i].height > 5 * avbl)
						return true;
					}
				}
			}

		return false;
		}

	/**
	 * Checks, whether all non-hydrogen atoms are mapped and whether every reactant atom has exactly one assigned product atom.
	 * @return
	 */
	public boolean isPerfectlyMapped() {
		int atoms = 0;
		for (StereoMolecule reactant:mReactant) {
			reactant.ensureHelperArrays(Molecule.cHelperNeighbours);
			atoms += reactant.getAtoms();
			}
		for (StereoMolecule product:mProduct) {
			product.ensureHelperArrays(Molecule.cHelperNeighbours);
			atoms -= product.getAtoms();
			}
		if (atoms != 0)
			return false;	// reactant atom count is different from product atom count

		int maxMapNo = getHighestMapNo();

		boolean[] isUsed = new boolean[maxMapNo+1];

		for (StereoMolecule reactant:mReactant) {
			for (int atom=0; atom<reactant.getAtoms(); atom++) {
				int mapNo = reactant.getAtomMapNo(atom);
				if (isUsed[mapNo])
					return false;
				isUsed[mapNo] = true;
				}
			}

		for (StereoMolecule product:mProduct) {
			product.ensureHelperArrays(Molecule.cHelperNeighbours);
			for (int atom=0; atom<product.getAtoms(); atom++) {
				int mapNo = product.getAtomMapNo(atom);
				if (mapNo >= maxMapNo || !isUsed[mapNo])
					return false;
				isUsed[mapNo] = false;
				}
			}

		return true;
		}

	public int getHighestMapNo() {
		if (mMaxMapNo != -1)
			return mMaxMapNo;

		mMaxMapNo = 0;
		for (int i=0; i<getMolecules(); i++) {
			StereoMolecule mol = getMolecule(i);
			for (int atom=0; atom<mol.getAllAtoms(); atom++) {
				if (mMaxMapNo < mol.getAtomMapNo(atom))
					mMaxMapNo = mol.getAtomMapNo(atom);
				}
			}

		return mMaxMapNo;
		}

	/**
	 * Removes mapping numbers that are only used on one side of the reaction.
	 * Throws an exception if duplicate mapping numbers occur in reactants or products.
	 * @throws Exception
	 */
	public void validateMapping() throws Exception {
		StereoMolecule reactant, product;

		for (int i = 0; i < mReactant.size(); i++) {
			reactant = mReactant.get(i);
			for (int j = 0; j < reactant.getAllAtoms(); j++) {
				int mapNo = reactant.getAtomMapNo(j);
				if (mapNo != 0) {
					int found = 0;
					for (int k = 0; k < mProduct.size(); k++) {
						product = mProduct.get(k);
						for (int l = 0; l < product.getAllAtoms(); l++)
							if (product.getAtomMapNo(l) == mapNo)
								found++;
						}

					if (found == 0)
						reactant.setAtomMapNo(j, 0, false);
					else if (found > 1)
						throw new Exception("Duplicate mapping no in products");
					}
				}
			}

		for (int i = 0; i < mProduct.size(); i++) {
			product = mProduct.get(i);
			for (int j = 0; j < product.getAllAtoms(); j++) {
				int mapNo = product.getAtomMapNo(j);
				if (mapNo != 0) {
					int found = 0;
					for (int k = 0; k < mReactant.size(); k++) {
						reactant = mReactant.get(k);
						for (int l = 0; l < reactant.getAllAtoms(); l++)
							if (reactant.getAtomMapNo(l) == mapNo)
								found++;
						}

					if (found == 0)
						product.setAtomMapNo(j, 0, false);
					else if (found > 1)
						throw new Exception("Duplicate mapping no in reactants");
					}
				}
			}
		}

	/**
	 * This method determines the largest mapping number in use (maxMapNo), creates a boolean array[maxMapNo+1],
	 * and within this array flags every mapping number that refers to atoms, which change bonds in the course
	 * of the reaction. Mapped atoms that are connected to unpammed atoms are also considered being part of the
	 * reaction center. If the reaction is unmapped or has no reactants or products, then null is returned.
	 * @return null or boolean array indicating, which mapping numbers are referring to the reaction center
	 */
	public boolean[] getReactionCenterMapNos() {
		if (getReactants() == 0 || getProducts() == 0)
			return null;

		int maxMapNo = getHighestMapNo();
		if (maxMapNo == 0)
			return null;

		// build mappings from mapNo to atom index in all products
		int[][] mapNo2Atom = new int[getProducts()][];
		for (int i=0; i<getProducts(); i++) {
			StereoMolecule product = getProduct(i);
			product.ensureHelperArrays(Molecule.cHelperParities);
			mapNo2Atom[i] = new int[maxMapNo+1];
			Arrays.fill(mapNo2Atom[i], -1);
			for (int atom=0; atom<product.getAllAtoms(); atom++) {
				int mapNo = product.getAtomMapNo(atom);
				if (mapNo != 0 && mapNo2Atom[i][mapNo] != -1)
					return null;	// same mapNo used twice in same product
				mapNo2Atom[i][mapNo] = atom;
			}
		}

		// find reaction centers as those mapped atoms that change bonding or are connected to unmapped atoms
		boolean[] isReactionCenter = new boolean[maxMapNo+1];
		for (int i=0; i<getReactants(); i++) {
			StereoMolecule reactant = getReactant(i);
			reactant.ensureHelperArrays(Molecule.cHelperParities);
			for (int rAtom=0; rAtom<reactant.getAllAtoms(); rAtom++) {
				int mapNo = reactant.getAtomMapNo(rAtom);
				if (mapNo != 0 && !isReactionCenter[mapNo]) {
					for (int j=0; j<getProducts(); j++) {
						int pAtom = mapNo2Atom[j][mapNo];
						if (pAtom != -1) {
							StereoMolecule product = getProduct(j);
							if (reactant.getConnAtoms(rAtom) != product.getConnAtoms(pAtom)) {
								isReactionCenter[mapNo] = true;
								break;
								}
							if (reactant.getAtomParity(rAtom) != product.getAtomParity(pAtom)) {	// TODO match neighbor positions for proper parity check
								isReactionCenter[mapNo] = true;
								break;
								}
							for (int k=0; k<reactant.getConnAtoms(rAtom); k++) {
								int connMapNo = reactant.getAtomMapNo(reactant.getConnAtom(rAtom, k));
								if (connMapNo == 0) {
									isReactionCenter[mapNo] = true;
									}
								else {
									int rBond = reactant.getConnBond(rAtom, k);
									boolean connMapNoFound = false;
									for (int l=0; l<product.getConnAtoms(pAtom); l++) {
										int productConnMapNo = product.getAtomMapNo(product.getConnAtom(pAtom, l));
										if (productConnMapNo == 0) {
											isReactionCenter[mapNo] = true;
											break;
											}
										if (productConnMapNo == connMapNo) {
											connMapNoFound = true;
											int pBond = product.getConnBond(pAtom, l);
											if ((reactant.isDelocalizedBond(rBond) ^ product.isDelocalizedBond(pBond))
													|| (!reactant.isDelocalizedBond(rBond)
													&& (reactant.getBondOrder(rBond) != product.getBondOrder(pBond)
													|| reactant.getBondParity(rBond) != product.getBondParity(pBond)))) {	// TODO match neighbor positions for proper parity check
												isReactionCenter[mapNo] = true;
												isReactionCenter[connMapNo] = true;
												break;
												}
											break;
											}
										}
									if (!connMapNoFound) {
										isReactionCenter[mapNo] = true;
										}
									}
								}
							}
						}
					}
				}
			}
		return isReactionCenter;
		}

	/**
	 * Fills an array mapping on the atoms for the given molecule of this reaction. Array bits are set if the respective atom
	 * has a mapping number flagged to be part of the reaction center, or if the atom has no mapping number but is connected
	 * to a mapped atom. The isReactionCenterAtom must be initialized with false before calling this method.
	 * @param moleculeNo reaction molecule index
	 * @param isReactionCenterMapNo flagged list of all reaction center mapping numbers, typically from getReactionCenterMapNos()
	 * @param isReactionCenterAtom null or empty array not smaller than the addressed molecule's total atom count
	 * @param reactionCenterAtom null or array not smaller than the addressed molecule's total atom count
	 * @return number of discovered reaction center atoms
	 */
	public int getReactionCenterAtoms(int moleculeNo, boolean[] isReactionCenterMapNo, boolean[] isReactionCenterAtom, int[] reactionCenterAtom) {
		StereoMolecule mol = getMolecule(moleculeNo);

		if (isReactionCenterAtom == null)
			isReactionCenterAtom = new boolean[mol.getAllAtoms()];

		int atomCount = 0;

		// mark atoms with known mapping numbers to be reaction centers
		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			if (isReactionCenterMapNo[mol.getAtomMapNo(atom)]) {
				isReactionCenterAtom[atom] = true;
				if (reactionCenterAtom != null)
					reactionCenterAtom[atomCount] = atom;
				atomCount++;
				}
			}

		// mapped atoms that connect to non mapped atoms are aleady covered, the non mapped ones, however, not
		for (int bond=0; bond<mol.getAllBonds(); bond++) {
			int atom1 = mol.getBondAtom(0, bond);
			int atom2 = mol.getBondAtom(1, bond);
			if (mol.getAtomMapNo(atom1) == 0 ^ mol.getAtomMapNo(atom2) == 0) {
				if (!isReactionCenterAtom[atom1]) {
					isReactionCenterAtom[atom1] = true;
					if (reactionCenterAtom != null)
						reactionCenterAtom[atomCount] = atom1;
					atomCount++;
					}
				if (!isReactionCenterAtom[atom2]) {
					isReactionCenterAtom[atom2] = true;
					if (reactionCenterAtom != null)
						reactionCenterAtom[atomCount] = atom2;
					atomCount++;
					}
				}
			}

		return atomCount;
		}

	/**
	 * Merges all reactants into one molecule and all rpoducts into another ad creates a new Reaction object from those.
	 * @return new reaction from merged reactants and merged products
	 */
	public Reaction getMergedCopy() {
		Reaction mergedReaction = new Reaction();
		if (mReactant.size() != 0) {
			StereoMolecule reactant = new StereoMolecule(mReactant.get(0));
			for (int i=1; i<mReactant.size(); i++)
				mReactant.get(i).copyMolecule(reactant);
			mergedReaction.addReactant(reactant);
			}
		if (mProduct.size() != 0) {
			StereoMolecule product = new StereoMolecule(mProduct.get(0));
			for (int i=1; i<mProduct.size(); i++)
				mProduct.get(i).copyMolecule(product);
			mergedReaction.addProduct(product);
			}
		return mergedReaction;
		}

	/*	public void removeEmptyMolecules() {
		int size = mReactant.size();
		for (int i = size-1; i >= 0; i--) {
			StereoMolecule mol = mReactant.get(i);
			if (mol.getAllAtoms() == 0) {
				mReactant.remove(i);
				}
			}
		size = mProduct.size();
		for (int i = size-1; i >= 0; i--) {
			StereoMolecule mol = mProduct.get(i);
			if (mol.getAllAtoms() == 0) {
				mProduct.remove(i);
				}
			}
		}*/
	}
