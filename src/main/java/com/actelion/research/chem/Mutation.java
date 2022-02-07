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

package com.actelion.research.chem;

import com.actelion.research.util.DoubleFormat;

public class Mutation {
    public static final int MUTATION_ADD_ATOM = 1;
    public static final int MUTATION_INSERT_ATOM = 2;
    public static final int MUTATION_CHANGE_ATOM = 4;
    public static final int MUTATION_CUTOUT_ATOM = 8;
    public static final int MUTATION_DELETE_ATOM = 16;
    public static final int MUTATION_CLOSE_RING = 32;
    public static final int MUTATION_CHANGE_BOND = 64;
    public static final int MUTATION_DELETE_BOND = 128;
    public static final int MUTATION_CHANGE_RING = 256;
	public static final int MUTATION_CLOSE_RING_AND_AROMATIZE = 512;
	public static final int MUTATION_TOGGLE_AMID_SULFONAMID = 1024;
    public static final int MUTATION_MIGRATE = 2048;
    public static final int MUTATION_SWAP_SUBSTITUENT = 4096;
    public static final int MUTATION_DELETE_SUBSTITUENT = 8192;
    public static final int MUTATION_CUTOUT_SFRAGMENT = 16384;
	public static final int MUTATION_INVERT_PARITY = 32768;

	public static final int[][] cAllowedAtomicNo = 
		{ { 5, 6, 7, 8, 9, 15, 16, 17, 35, 53 },
		  { 6, 7, 8, 15, 16 },
		  { 6, 7 } };

	public int mMutationType;
	public int mWhere1;
    public int mWhere2;
	public int mSpecifier1;
	public int mSpecifier2;
	public int[] mAtomList;
	public double mProbability;

	public Mutation() {
		}

	public Mutation(int mutationType,
	                int where1, int where2,
	                int specifier1, int specifier2,
	                double probability) {
		mMutationType	= mutationType;
		mWhere1			= where1;
        mWhere2         = where2;
		mSpecifier1		= specifier1;
		mSpecifier2		= specifier2;
		mProbability = probability;
		}

	public Mutation(int mutationType, int where1, int where2, int specifier1, int[] atomList,
					double probability) {
		mMutationType	= mutationType;
		mWhere1			= where1;
		mWhere2         = where2;
		mSpecifier1		= specifier1;
		mAtomList		= atomList;
		mProbability = probability;
	}

	public String toString() {
        switch (mMutationType) {
        case MUTATION_ADD_ATOM:
            return "Atom Addition; AtAtom:"+mWhere1
                  +" AtomicNo:"+mSpecifier1
                  +" BondType:"+mSpecifier2
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_INSERT_ATOM:
            return "Atom Insertion; AtBond:"+mWhere1
                  +" AtomicNo:"+mSpecifier1
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_CHANGE_ATOM:
            return "Atom Change; Atom:"+mWhere1
                  +" AtomicNo:"+mSpecifier1
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_DELETE_ATOM:
            return "Atom Deletion; Atom:"+mWhere1
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_CUTOUT_ATOM:
            return "Atom CutOut; Atom:"+mWhere1
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_CLOSE_RING:
            return "Bond Addition; FromAtom:"+mWhere1
                  +" ToAtom:"+mWhere2
                  +" BondType:"+mSpecifier1
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_CHANGE_BOND:
            return "Bond Change; Bond:"+mWhere1
                  +" BondType:"+mSpecifier1
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_DELETE_BOND:
            return "Bond Deletion; Bond:"+mWhere1
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_CHANGE_RING:
            return "Ring Change; Ring:"+mWhere1
				  +" Probability:"+DoubleFormat.toString(mProbability);
		case MUTATION_CLOSE_RING_AND_AROMATIZE:
			return "Close Ring And Aromatize; FromAtom:"+mWhere1
				  +" ToAtom:"+mWhere2
				  +" Probability:"+DoubleFormat.toString(mProbability);
		case MUTATION_TOGGLE_AMID_SULFONAMID:
			return "Toggle Amid-Sulfonamid; Atom:"+mWhere1
				  +" Oxygen:"+mWhere1
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_MIGRATE:
            return "Group Migration; Bond:"+mWhere1
                  +" originalAtom:"+mSpecifier1
                  +" newAtom:"+mSpecifier2
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_SWAP_SUBSTITUENT:
            return "Swap Substituent; coreAtom1:"+mWhere1
                  +" firstAtom1:"+mSpecifier1
                  +" coreAtom2:"+mWhere2
                  +" firstAtom2:"+mSpecifier2
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_DELETE_SUBSTITUENT:
            return "Delete Substituent; coreAtom:"+mWhere1
                  +" firstAtom:"+mSpecifier1
				  +" Probability:"+DoubleFormat.toString(mProbability);
        case MUTATION_CUTOUT_SFRAGMENT:
            return "CutOut Fragment; rootAtom:"+mWhere1
                  +" new bond from atom1:"+mSpecifier1+" to atom2:"+mSpecifier2
				  +" Probability:"+DoubleFormat.toString(mProbability);
	        case MUTATION_INVERT_PARITY:
	        return "Invert Parity; "
			      +(mWhere1 != -1 ? "atom:"+mWhere1 : "bond:"+mWhere2)
		          +" Probability:"+DoubleFormat.toString(mProbability);
            }
        return "Unknown Mutation";
        }
    }