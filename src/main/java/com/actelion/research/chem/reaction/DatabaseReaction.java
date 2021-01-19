/*
 * Copyright 2017 Idorsia Pharmaceuticals Ltd., Hegenheimermattweg 91, CH-4123 Allschwil, Switzerland
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


public class DatabaseReaction extends Reaction {
	private static final long serialVersionUID = 1L;

	private int		mReactionRegNo;		// reaction's internal registry number
	private int 	mMoleculeRegNo[];	// molecules internal registry nums
	private int		mYield;				// yield (-1 = not available)

//	private int initialize() {}

//	CReactionText	*getReactionTextP();
//	CDatabase		*getDatabaseP();
//	void setDatabaseP		( CDatabase *databaseP );
//	void reset				();

	public void setReactionRegNo(int regNo) {
		mReactionRegNo = regNo;
		}

	/**
	 * @param yield (-1 = not available)
	 */
	public void setReactionYield(int yield) {
		mYield = yield;
		}

	public void setMoleculeRegNo(int mol, int regNo) {
		if (mMoleculeRegNo == null)
			mMoleculeRegNo = new int[getMolecules()];

		mMoleculeRegNo[mol] = regNo;
		}

	public int getReactionRegNo() {
		return mReactionRegNo;
		}

	public int getReactionYield() {
		return mYield;
		}

	public int getMoleculeRegNo(int mol) {
		return (mMoleculeRegNo == null || mMoleculeRegNo.length<= mol) ? -1 : mMoleculeRegNo[mol];
		}
	}
