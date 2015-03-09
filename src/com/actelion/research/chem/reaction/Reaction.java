/*
 * Copyright 2014 Actelion Pharmaceuticals Ltd., Gewerbestrasse 16, CH-4123 Allschwil, Switzerland
 *
 * This file is part of DataWarrior.
 * 
 * DataWarrior is free software: you can redistribute it and/or modify it under the terms of the
 * GNU General Public License as published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 * 
 * DataWarrior is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with DataWarrior.
 * If not, see http://www.gnu.org/licenses/.
 *
 * @author Thomas Sander
 */

package com.actelion.research.chem.reaction;

import java.util.ArrayList;

import com.actelion.research.chem.DrawingObjectList;
import com.actelion.research.chem.StereoMolecule;

public class Reaction implements java.io.Serializable {
	static final long serialVersionUID = 0x2006CAFE;

	private ArrayList<StereoMolecule> mReactant;
	private ArrayList<StereoMolecule> mProduct;
	private DrawingObjectList mDrawingObjectList;
	private String mName;
	private boolean mReactionLayoutRequired;

	public Reaction() {
		mReactant = new ArrayList<StereoMolecule>();
		mProduct = new ArrayList<StereoMolecule>();
		}

	public Reaction(String name) {
		this();
		mName = name;
		}

	public Reaction(Reaction rxn) {
		this();
		int r = (rxn == null) ? 0 : (rxn.mReactant == null ? 0 : rxn.mReactant.size());
		int p = (rxn == null) ? 0 : (rxn.mProduct == null ? 0 : rxn.mProduct.size());
		for (int i = 0; i < r; i++)
			mReactant.add(new StereoMolecule(rxn.getReactant(i)));
		for (int i = 0; i < p; i++)
			mProduct.add(new StereoMolecule(rxn.getProduct(i)));
		mDrawingObjectList = new DrawingObjectList(rxn.getDrawingObjects());
		}

	public Reaction(StereoMolecule[] mol, int reactantCount) {
		this();
		if (mol != null) {
			for (int i = 0; i < reactantCount; i++)
				mReactant.add(mol[i]);
			for (int i = reactantCount; i < mol.length; i++)
				mProduct.add(mol[i]);
			}
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
		}


	public void addReactant(StereoMolecule reactant, int position) {
		mReactant.add(position, reactant);
		}


	public void addProduct(StereoMolecule product) {
		mProduct.add(product);
		}


	public void addProduct(StereoMolecule product, int position) {
		mProduct.add(position, product);
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


	public boolean isReactionLayoutRequired() {
		return mReactionLayoutRequired;
		}


	public void setReactionLayoutRequired(boolean b) {
		mReactionLayoutRequired = b;
		}

	/**
	 * DISABLED because assumption is not correct. 1st atom's coordinates are only 0,0
	 * if they come from parsing an idcode with(!) (relative) coordinates. If coordinates
	 * are created on the fly by the IDCodeParser, then this is different. Now we use a
	 * dedicated flag.
	 * 
	 * Atom coordinates may be relative or absolute.<br> If coordinates are relative,
	 * then the first atom of every involved molecule is located at x=0.0, y=0.0.
	 * For depicting this reaction, molecules need to be individually translated
	 * and scaled to layout the reaction for display.<br>
	 * If coordinates are absolute, then the relative orientation of molecules
	 * (and drawing objects) in the reaction context is correct, such that entire
	 * reaction can be scaled and translated for display.
	 * @return
	 *
	public boolean hasAbsoluteCoordinates() {
		for (StereoMolecule mol:mReactant)
			if (mol.getAllAtoms() != 0 && (mol.getAtomX(0) != 0.0 || mol.getAtomY(0) != 0.0))
				return true;
		for (StereoMolecule mol:mProduct)
			if (mol.getAllAtoms() != 0 && (mol.getAtomX(0) != 0.0 || mol.getAtomY(0) != 0.0))
				return true;
		return false;
		}	*/


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

	public void normalize() {
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
		}
	}
