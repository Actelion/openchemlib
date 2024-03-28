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

package com.actelion.research.chem.reaction;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;

public class ReactionClassifier {
	private static final boolean DEBUG = true;

//	private static final int MAXAROMS = 19;
	private static final int MAXFATMS = 64;		// maximal no of atoms per fragment
	private static final int MAXATMCNO = 190;	// max allowed atomic number
	private static final int MAXCONNS = 6;

//	private static final int MAXMOLS = 16;
	private static final int MAXATOMS = 256;
//	private static final int MAXBONDS = 256;
	private static final int MAXMPPNGNO = 512;

	private static final int MAXINDICES = Classification.MAXINDICES;
//	private static final int INDEXLEN = Classification.INDEXLEN;
	private static final int MAXFLASH = Classification.MAXFLASH;
	private static final int MAXUNITS = Classification.MAXUNITS;
//	private static final int MAXENDUR = Classification.MAXENDUR;

//	private static final int  RCLASSES	= 13;
	private static final int  S_CONST	=  0;
	private static final int HD_D_CONST	=  1;
	private static final int D_D_CONST	=  2;
	private static final int HD_REAR	=  3;
	private static final int D_REAR		=  4;
	private static final int S_FRAG		=  5;
	private static final int C1_REFU	=  6;
	private static final int LONG_REFU	=  7;
	private static final int HETERO_REFU = 8;
	private static final int HD_D_FRAG	=  9;
	private static final int D_D_FRAG	= 10;
	private static final int E_RING_C	= 11;
	private static final int E_RING_O	= 12;

	private static final int BONDBITS	=  5;		// bits in indexfile: C-C bond specific data
	private static final int HRXNBITS	= 25;		//					half rxn specific data
	private static final int DRXNBITS	= 11;		//			two end half rxn specific data
	private static final int CRXNBITS	=  4;		//	bits for electrocyclic reaction class
	private static final int REARBITS	= 13;		//			half rearrangement specific data
	private static final int PROPBITS	= 11;		// nr of bits of carbon describing properties

//	private static final int MAXPRUNOPTIONS = 58;

	public static final int cErrorNoError = 0;
	public static final int cErrorNoChangeNorEFG = 1;
	public static final int cErrorEductMapNoOverused = 2;
	public static final int cErrorProdMapNoOverused1 = 3;
	public static final int cErrorProdMapNoOverused2 = 4;
	public static final int cErrorMapNoNotInProduct = 5;
	public static final int cErrorMapNoNotInEduct = 6;
	public static final int cErrorDupProdMapNoDifEd = 7;
//	public static final int cErrorMaxNoOfMolsReached = 8;
	public static final int cErrorProdRemapFailed = 9;
	public static final int cErrorEduFragPartMapped = 10;
	public static final int cErrorFragmentAtomLimit = 11;
	public static final int cErrorProdFragPartMapped = 12;
	public static final int cErrorUnexpected = 13;
	public static final int cErrorUnMappedCInConOrFr = 14;
	public static final int cErrorCCBondCleavageLimit = 15;
	public static final int cErrorCCBondCreationLimit = 16;
	public static final int cErrorComplexReaction = 17;
	public static final int cErrorNoChangingAtoms = 18;
	public static final int cErrorForkedOrLongStrand = 19;
	public static final int cError2AlphasInSameStrand = 20;
	public static final int cErrorHRClassifyError = 21;
	public static final int cErrorCRClassifyError = 22;
	public static final int cErrorDRClassifyError = 23;
	public static final int cErrorRAClassifyError = 24;
	public static final int cErrorREClassifyError = 25;
	public static final int cErrorIncoOrLeavMissing = 26;
	public static final int cErrorUnitReactionLimit = 27;
	public static final int cErrorNoDatabaseReaction = 28;

	public static final int cIndexNone = 0;
	public static final int cIndexOnePermToFile = 1;
	public static final int cIndexFullPermutation = 2;

	private static final int cAtomPiChange			= 0x0100;	// TODO use real value
	private static final int cAtomZChange			= 0x0200;	// TODO use real value
	private static final int cAtomSigmaChange		= 0x0400;	// TODO use real value
	private static final int cAtomChanges			= 0x1000;	// TODO use real value
	private static final int cAtomNoCounterAtom		= 0x2000;	// TODO use real value
	private static final int cAtomDBondToHetero		= 0x4000;	// TODO use real value
	private static final int cAtomNotClassifiedYet	= 0x8000;	// TODO use real value

	private static final int CARBON = 1;
//	private static final int HYDROGEN	= 3;

	private static int data,mask;
	private static int indexnum,indexpoin,availbits;

	private ClassificationData mClassificationData;
	private BufferedWriter gErrout;
	private Reaction mRxn;
	private Classification mResult;
	private int mIndexToCreate,mUnitRxn;
	private int[][] mAtomType,mSigma,mPi,mH,mZ,mAtomFlags,mCorProd,mCorAtom;
	private int[][][] mConnMpNo,mConnCMNo,mConnType;

	public ReactionClassifier() {
		mClassificationData = ClassificationData.getInstance();
		}


	public int classify(Reaction theReaction) {
		return classify(cIndexNone, theReaction);
		}


	public int classify(int indexToCreate, /*CBatchClassifier theBatchClassifierP,*/ Reaction theReaction) {
//TODO		CEFGClassifier	efgClassifier;
		
		mIndexToCreate = indexToCreate;
		mRxn = theReaction;
//		mBatchClassifierP = theBatchClassifierP;

		if (indexToCreate != cIndexNone && !(theReaction instanceof DatabaseReaction))
			return cErrorNoDatabaseReaction;

		for (int mol=0; mol<mRxn.getMolecules(); mol++)
			mRxn.getMolecule(mol).ensureHelperArrays(Molecule.cHelperParities);

		int err = ensureCleanMapping();
		if (err != 0)
			return err;
		
		mResult = new Classification();
		
		if (DEBUG)
			try {
				gErrout = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("errout.txt"), StandardCharsets.UTF_8));
				}
			catch (IOException ioe) {}

		err = classrxn();

		if (DEBUG)
			try { gErrout.close(); } catch (IOException ioe) {}
		
		if (err != 0) {
			mResult = null;
			return err;
			}

//TODO		efgClassifier = new CEFGClassifier();
//TODO		err = efgClassifier.classifyEFG( mRxn, this, mResult );
		
		if (mUnitRxn == 0 && mResult.mEnduringFGs == 0) {
			mResult = null;
			return cErrorNoChangeNorEFG;
			}

		return err;
		}


	public Classification getClassificationResult() {
		return mResult;
		}


	/**
	 * Checks if any product molecule contains some mapping numbers twice.
	 * If none of the products contains duplicate mapping numbers 0 is returned.
	 * If more than one product contain duplicate mapping numbers an error is returned.
	 * If exactly one product contains duplicate mapping numbers the educts containing
	 * these mapping numbers are located. If more than one educt contain some of these
	 * mapping numbers an error is returned. Otherwise the located educt is duplicated
	 * and new mapping numbers are created to eliminate duplicate numbers.
	 * @return
	 */
	private int ensureCleanMapping() {
		int[] eduMapNoCount = new int[MAXMPPNGNO];
		int[] proMapNoCount = new int[MAXMPPNGNO];
		int[] eduMol = new int[MAXMPPNGNO];
		int[] eduAtm = new int[MAXMPPNGNO];
		int[] proMol = new int[MAXMPPNGNO];

		  				// set fragments to unmapped if fragment members have
		  				// mapping numbers which are already used in the same educt
		int[] doubleMapped = new int[MAXATOMS];
		for (int e=0; e<mRxn.getReactants(); e++) {
			StereoMolecule mol = mRxn.getReactant(e);
			boolean doubleFound = false;
			for (int i=1; i<MAXMPPNGNO; i++)
				eduMapNoCount[i] = 0;
		
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				int mapNo = mol.getAtomMapNo(atm);
				if (mapNo != 0) {
					eduMapNoCount[mapNo]++;
					if (eduMapNoCount[mapNo] == 2)
						doubleFound = true;
					else if (eduMapNoCount[mapNo] > 2)
						return cErrorEductMapNoOverused;
					}
				}
		
			if (!doubleFound)
				continue;
		
			for (int atm=0; atm<mol.getAtoms(); atm++) {	// flag atoms with ambiguous mppngNo
				doubleMapped[atm] = 0;
				int mapNo = mol.getAtomMapNo(atm);
				if (mapNo != 0)
					if (eduMapNoCount[mapNo] == 2)
						doubleMapped[atm] = 1;
				}
		
			for (int atm=0; atm<mol.getAtoms(); atm++) // unmap one of each ambiguous fragment
				if (doubleMapped[atm] != 0)
					unmapFragment(mol, atm, doubleMapped);
			}
		
		
		for (int i=1; i<MAXMPPNGNO; i++)
			eduMapNoCount[i] = proMapNoCount[i] = 0;
		
		for (int e=0; e<mRxn.getReactants(); e++) {
			StereoMolecule mol = mRxn.getReactant(e);
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				int mapNo = mol.getAtomMapNo(atm);
				if (mapNo != 0) {
					if (eduMapNoCount[mapNo] != 0) {
						if (eduMol[mapNo] != e) {	// if mapping no was already used in another educt
							mol.setAtomMapNo(atm, 0, false);
							continue;
							}							// same mapping no used twice within
						return cErrorEductMapNoOverused;	// one educt; shouldn't happen 
						}								// if fragment stuff above works
		
					eduMapNoCount[mapNo]++;
					eduMol[mapNo] = e;
					eduAtm[mapNo] = atm;
					}
				}
			}

		int dirtyProduct = -1;
		for (int p=0; p<mRxn.getProducts(); p++) {
			StereoMolecule mol = mRxn.getProduct(p);
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				int mapNo = mol.getAtomMapNo(atm);
				if (mapNo != 0) {
					if (proMapNoCount[mapNo] != 0) {
						if (proMol[mapNo] != p) {
							mol.setAtomMapNo(atm, 0, false);
							continue;
							}
						if (proMapNoCount[mapNo] > 1)
							return cErrorProdMapNoOverused1;
						if (dirtyProduct != -1 && dirtyProduct != p)
							return cErrorProdMapNoOverused2;
						dirtyProduct = p;
						}
		
					proMapNoCount[mapNo]++;
					proMol[mapNo] = p;
					}
				}
			}
		
		if (dirtyProduct == -1)  // all products are free of duplicate mapping numbers
			return cErrorNoError;
		
		int highestMappingNo = 0;
		for (int i=1; i<MAXMPPNGNO; i++) {
			if (eduMapNoCount[i] > 0 && proMapNoCount[i] > 0) {
				highestMappingNo = i;
				continue;
				}
			if (eduMapNoCount[i] != 0)
				return cErrorMapNoNotInProduct;
			if (proMapNoCount[i] != 0)
				return cErrorMapNoNotInEduct;
			}
		
		int dirtyEduct = -1;
		for (int i=1; i<=highestMappingNo; i++) {
			if (proMapNoCount[i] > 1) {
				if (dirtyEduct != -1 && dirtyEduct != eduMol[i])
					return cErrorDupProdMapNoDifEd;
				dirtyEduct = eduMol[i];
				}
			}

		StereoMolecule duplicateEduct = new StereoMolecule(mRxn.getReactant(dirtyEduct));
		duplicateEduct.ensureHelperArrays(Molecule.cHelperParities);
		mRxn.addReactant(duplicateEduct, dirtyEduct);

		StereoMolecule mol = mRxn.getProduct(dirtyProduct);
		for (int atm=0; atm<mol.getAtoms(); atm++) { // flag atoms with ambiguous mppngNo
			doubleMapped[atm] = 0;
			int mapNo = mol.getAtomMapNo(atm);
			if (mapNo != 0) {
				if (proMapNoCount[mapNo] == 2)
					doubleMapped[atm] = 1;
				}
			}

		for (int atm=0; atm<mol.getAtoms(); atm++) // remap one of every ambiguous fragment
			if (doubleMapped[atm] != 0)
				if (tryRemapFragment( mRxn.getReactant(dirtyEduct), duplicateEduct,
						mRxn.getProduct(dirtyProduct), atm, highestMappingNo, eduAtm, doubleMapped ))
					return cErrorNoError;

		return cErrorProdRemapFailed;
		}


	private void unmapFragment(StereoMolecule mol, int atm, int[] atmFlags) {
		for (int i=0; i<mol.getAtoms(); i++) {
			if (i == atm)
				continue;
			if (atmFlags[i] != 0) {
				if (mol.getAtomMapNo(atm) == mol.getAtomMapNo(i)) {
					atmFlags[i] = 0;
					break;
					}
				}
			}
		
		atmFlags[atm] = 0;
		mol.setAtomMapNo(atm, 0, false);
		
		for (int i=0; i<mol.getConnAtoms(atm); i++) {
			int connAtm = mol.getConnAtom(atm, i);
			if (atmFlags[connAtm] != 0)
				unmapFragment(mol, connAtm, atmFlags);
			}
		}


	private boolean tryRemapFragment(StereoMolecule educt1, StereoMolecule educt2, StereoMolecule product,
									 int proAtm, int highestMappingNo, int[] eduAtm, int[] atmFlags) {
		int atm,eductAtom,numbersNeeded;
		
		atmFlags[proAtm] |= 2;		// mark product atom to be member of located fragment
		checkAdjacent(educt1, product, proAtm, eduAtm, atmFlags );
									// this recursive call should mark all

		numbersNeeded = 0;
		for (atm=0; atm<product.getAtoms(); atm++)
			if ((atmFlags[atm] & 2) != 0)
				numbersNeeded++;
		if ((highestMappingNo + numbersNeeded) > MAXMPPNGNO)
			return false;
		
		boolean successful = true;  // default
		for (atm=0; atm<product.getAtoms(); atm++) {
			if ((atmFlags[atm] & 2) != 0) {
				eductAtom = eduAtm[product.getAtomMapNo(atm)];
				if ((atmFlags[eductAtom] & 4) != 0) {  // educt atom was already marked
					successful = false;
					break;
					}
				atmFlags[eductAtom] |= 4;
				}
			}

		if (successful) { // check if any mapped educt carbon is not covered by the located product fragment
			for (atm=0; atm<educt1.getAtoms(); atm++) {
				if (educt1.getAtomMapNo(atm) != 0) {
					if ((atmFlags[atm] & 4) == 0) {
						successful = false;
						break;
						}
					}
				}
			}
		
		if (successful) { // do the actual remapping finally
			for (atm=0; atm<product.getAtoms(); atm++) {
				if ((atmFlags[atm] & 2) != 0) {
					highestMappingNo++;
					eductAtom = eduAtm[product.getAtomMapNo(atm)];
					product.setAtomMapNo(atm, highestMappingNo, false);
					educt1.setAtomMapNo(eductAtom, highestMappingNo, false);
					if ((atmFlags[atm] & 1) == 0)
						educt2.setAtomMapNo(eductAtom, 0, false);
					}
				}
			}
		
		for (atm=0; atm<MAXATOMS; atm++)
			atmFlags[atm] &= 0xF9;		// reset educt and product fragment member flags
		
		return successful;
		}


	private void checkAdjacent(StereoMolecule educt, StereoMolecule product, int atm, int[] eduAtm, int[] atmFlags) {
		int i,j,eductAtm,connAtm,connMapNo,otherConn;
		
		eductAtm = eduAtm[product.getAtomMapNo(atm)];
		for (i=0; i<product.getConnAtoms(atm); i++) {
			connAtm = product.getConnAtom(atm, i);
			connMapNo = product.getAtomMapNo(connAtm);
		
			if (connMapNo == 0)		// exclude unmapped atoms from fragments
				continue;
		
			if ((atmFlags[connAtm] & 2) != 0)	// atom is already flagged as member
				continue;
		
			if (connMapNo == product.getAtomMapNo(atm))
				continue;
		
			boolean secondFound = false;
			// check if another atom with the same mapping number is also connected to atm
			for (j=0; j<product.getConnAtoms(atm) && j!=i; j++) {
				otherConn = product.getConnAtom(atm, j);
				if (connMapNo == product.getAtomMapNo(otherConn))
					secondFound = true;
				}
			if (secondFound)
				continue;
		
			boolean found = false;		// check educt if corresponding atom is connected to eductAtom
			for (j=0; j<educt.getConnAtoms(eductAtm); j++)
				if (connMapNo == educt.getAtomMapNo( educt.getConnAtom(eductAtm, j) ))
					found = true;
			if (!found)
				continue;
		
			atmFlags[connAtm] |= 2;			// mark atom to be member of located fragment
			checkAdjacent(educt, product, connAtm, eduAtm, atmFlags);
			}
		}


	/**
     *	This routine relies on the following information:
     *	reaction specific:
     *		mRxn.mEducts, mRxn.mMols, mRxn.mCatalysts,
     *	molecule specific:
     *		mRxn.mAtoms[], mRxn.mBonds[], mRxn.mAllAtoms[], mRxn.mAllBonds[],
     *	atom specific:
     *		mRxn.mAtomicNo[][], mRxn.mAtomMapNo[][], mRxn.mAtomCharge[][],
     *		mRxn.mAtomX[][], mRxn.mAtomY[][], mRxn.mAtomFlags[][],
     *		mConnAtoms[][], mConnAtom[][][], mConnBond[][][], mConnBondOrder[][],
     *	bond specific:
     *		mRxn.mBondAtom[][][], mRxn.mBondOrder[][]
	 * @return
	 */
	private int classrxn() {
		final int[] types = {  0,
		0x03,0x00,0x03,0x03,0x03,0x01,0x82,0x02,0x06,0x00,0x03,0x03,0x03,0x03,0x92,0x42,
		0x2E,0x00,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,
		0xB2,0x4A,0x4E,0x00,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,
		0x03,0x03,0x03,0x6A,0x6E,0x00,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,
		0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,
		0x03,0x03,0x03,0x00,0x00,0x00,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x03,
		0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x01,0x01,0x01,
		0x01,0x01,0x01,0x01,0x00,0x00,0x03,0x03,0xEE,0x01,0x00,0x03,0x00,0x00,0x01,0x00,
		0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x01,0x01,0x01,0x01,0x01,0x01,
		0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01 };

		int[] ccBFmol = new int[10];
		int[] ccBFatm = new int[10];
		int[] ccBCmol = new int[10];
		int[] ccBCatm = new int[10];
		int[] fragAtms = new int[MAXFATMS];

		for (int i=0; i<MAXUNITS; i++) {	// reset stereo information and flash atom info
			mResult.mStereoInfo[i] = 0;
			for (int j=0; j<MAXFLASH; j++) {
				mResult.mFlashMol[i][j] = -1;
				mResult.mFlashAtom[i][j] = -1;
				}
			}
		
		mUnitRxn = 0;

		
		mSigma = new int[mRxn.getMolecules()][];
		mPi = new int[mRxn.getMolecules()][];
		mH = new int[mRxn.getMolecules()][];
		mZ = new int[mRxn.getMolecules()][];
		mAtomType = new int[mRxn.getMolecules()][];
		mAtomFlags = new int[mRxn.getMolecules()][];
		mCorProd = new int[mRxn.getMolecules()][];
		mCorAtom = new int[mRxn.getMolecules()][];
		mConnMpNo = new int[mRxn.getMolecules()][][];
		mConnCMNo = new int[mRxn.getMolecules()][][];
		mConnType = new int[mRxn.getMolecules()][][];
		for (int m=0; m<mRxn.getMolecules(); m++) {	// initialize arrays
			StereoMolecule mol = mRxn.getMolecule(m);
			mSigma[m] = new int[mol.getAtoms()];
			mPi[m] = new int[mol.getAtoms()];
			mH[m] = new int[mol.getAtoms()];
			mZ[m] = new int[mol.getAtoms()];
			mAtomType[m] = new int[mol.getAtoms()];
			mAtomFlags[m] = new int[mol.getAtoms()];
			mCorProd[m] = new int[mol.getAtoms()];
			mCorAtom[m] = new int[mol.getAtoms()];
			mConnMpNo[m] = new int[mol.getAtoms()][MAXCONNS];
			mConnCMNo[m] = new int[mol.getAtoms()][MAXCONNS];
			mConnType[m] = new int[mol.getAtoms()][MAXCONNS];
			}

		for (int m=0; m<mRxn.getMolecules(); m++) {	// convert atomic number to atomtype
			StereoMolecule mol = mRxn.getMolecule(m);
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				int atomicNo = mol.getAtomicNo(atm);
				if ((atomicNo >= 1) && (atomicNo <= MAXATMCNO))
					mAtomType[m][atm] = types[atomicNo];
				else if (gErrout != null) {
					try {
						gErrout.write("Rxn: "+mRxn.getName()+", unknown atomic #: "+atomicNo);
						gErrout.newLine();
						}
					catch (IOException ioe) {}
					}
				}
			}
		
		for (int m=0; m<mRxn.getMolecules(); m++) {	// store no's and by value sorted
			StereoMolecule mol = mRxn.getMolecule(m);
			for (int bnd=0; bnd<mol.getBonds(); bnd++) {
				for (int i=0; i<2; i++) {
					int actlAtm = mol.getBondAtom(i, bnd);
					int connAtm = mol.getBondAtom(1-i, bnd);
		
					int remotTyp = mAtomType[m][connAtm];		// sorted atomtypes
					for (int k=0; k<mol.getBondOrder(bnd); k++) {
						for (i=0; i<MAXCONNS; i++) {
							if (remotTyp > mConnType[m][actlAtm][i]) {
								if (mConnType[m][actlAtm][i] != 0) {
									for (int j=MAXCONNS-1; j>i; j--)
										mConnType[m][actlAtm][j] = mConnType[m][actlAtm][j-1];
									}
								mConnType[m][actlAtm][i] = remotTyp;
								break;
								}
							}
						}
					}
				}

			if (gErrout != null) {
				for (int atm=0; atm<mol.getAtoms(); atm++) {
					try {
						gErrout.write("Atom:"+atm+"; mConnAtom:");
						for (int i=0; i<mol.getConnAtoms(atm); i++)
							gErrout.write(" "+mol.getConnAtom(atm, i));
						gErrout.newLine();
						}
					catch (IOException ioe) {}
					}
				}
			}
		
		for (int e=0; e<mRxn.getReactants(); e++) {	// convert mapping information
			StereoMolecule mol = mRxn.getReactant(e);
			int[] foundMsk = new int[MAXATOMS];

			for (int atm=0; atm<mol.getAtoms(); atm++) {
				int corNo = mol.getAtomMapNo(atm);
				if (corNo == 0) {
					if (mAtomType[e][atm] == CARBON) {	// unmapped educt carbons
						if (foundMsk[atm] == 0) {
							int retval = findFrag(e,atm,fragAtms,foundMsk);
							if (retval != 0)
								return retval;

							for (int i=1; i<=fragAtms[0]; i++)
								if (mol.getAtomMapNo(fragAtms[i]) != 0)
									return cErrorEduFragPartMapped;
							}
						}
					mCorProd[e][atm] = 0;
					mCorAtom[e][atm] = 255;
					continue;
					}
				boolean found = false;
				for (int p=mRxn.getReactants(); p<mRxn.getMolecules(); p++) {
					StereoMolecule product = mRxn.getMolecule(p);
					for (int i=0; i<product.getAtoms(); i++) {
						if (corNo == product.getAtomMapNo(i)) {
							mCorProd[e][atm] = p;
							mCorAtom[e][atm] = i;
							mH[p][i] = 1;				// indicate assignment
							found = true;				// temporarily in h[][]
							break;
							}
						}
					if (found) break;
					}
				if (!found) return cErrorMapNoNotInProduct;	// no mapping partner in product
				}
			}
		
											// pro atms without corresponding edu atms
		for (int p=mRxn.getReactants(); p<mRxn.getMolecules(); p++) {
			StereoMolecule product = mRxn.getMolecule(p);
			int[] foundMsk = new int[MAXATOMS];
		
			for (int atm=0; atm<product.getAtoms(); atm++) {
				if (product.getAtomMapNo(atm) == 0) {
					if (mAtomType[p][atm] != CARBON) {
						mAtomFlags[p][atm] |= cAtomChanges + cAtomNoCounterAtom;
						continue;
						}
					else {
						if (foundMsk[atm] == 0) {
							int retval = findFrag(p,atm,fragAtms,foundMsk);
							if (retval != 0) return retval;
							for (int i=1; i<=fragAtms[0]; i++) {
								mAtomFlags[p][fragAtms[i]] |= cAtomChanges + cAtomNoCounterAtom;
								if (product.getAtomMapNo(fragAtms[i]) != 0) return cErrorProdFragPartMapped;
								}							// unmapped product carbon
							}
						}
					}
				else if (mH[p][atm] != 1) {
					boolean found = false;
					for (int m=mRxn.getReactants(); m<mRxn.getMolecules(); m++) {
						StereoMolecule mol = mRxn.getMolecule(m);
						for (int j=0; j<mol.getAtoms(); j++) {
							if ((m == p) && (j == atm)) continue;
							if (mol.getAtomMapNo(j) == product.getAtomMapNo(atm)) {
								found = true;
								for (int k=0; k<MAXCONNS; k++)
									if (mConnType[p][atm][k] != mConnType[m][j][k])
										found = false;
								if (found)
									break;	// number in product, allowed, when same connType
								}
							}
						if (found) break;
						}
					if (!found) return cErrorMapNoNotInEduct;	// no mapping partner in educt
					}
				}
			}
		
		for (int m=0; m<mRxn.getMolecules(); m++) {
			StereoMolecule mol = mRxn.getMolecule(m);
			for (int bnd=0; bnd<mol.getBonds(); bnd++) {
				for (int atm=0; atm<2; atm++) {
					int actlAtm = mol.getBondAtom(atm, bnd);
					int connAtm = mol.getBondAtom(1-atm, bnd);
		
										// store sorted mapping no's of connected atoms
										// one entry for every (!) sigma or pi bond
					for (int k=0; k<mol.getBondOrder(bnd); k++) {
						for (int i=0; i<MAXCONNS; i++) {
							if (mol.getAtomMapNo(connAtm) > mConnMpNo[m][actlAtm][i]) {
								if (mConnMpNo[m][actlAtm][i] != 0) {
									for (int j=MAXCONNS-1; j>i; j--)
										mConnMpNo[m][actlAtm][j] =
											mConnMpNo[m][actlAtm][j-1];
									}
								mConnMpNo[m][actlAtm][i] = mol.getAtomMapNo(connAtm);
								break;
								}
							}
						}
		
									// store sorted mapping no's of connected carbons
									// independent of bond order only one entry
					if (mAtomType[m][connAtm] == CARBON) {
						for (int i=0; i<MAXCONNS; i++) {
							if (mol.getAtomMapNo(connAtm) > mConnCMNo[m][actlAtm][i]) {
								if (mConnCMNo[m][actlAtm][i] != 0) {
									for (int j=MAXCONNS-1; j>i; j--)
										mConnCMNo[m][actlAtm][j] =
											mConnCMNo[m][actlAtm][j-1];
									}
								mConnCMNo[m][actlAtm][i] = mol.getAtomMapNo(connAtm);
								break;
								}
							}
						}
		
					switch (3 & mAtomType[m][connAtm]) {  // determine z,pi,sigma,h
						case 1:
							mSigma[m][actlAtm] += 1;
							mPi[m][actlAtm] += mol.getBondOrder(bnd) - 1;
							break;
						case 2:
							mZ[m][actlAtm] += mol.getBondOrder(bnd);
							break;
						}
					}
				}
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				if (mAtomType[m][atm] == CARBON) {
					if (mol.getAtomCharge(atm) == 1)  // carbenium
						mZ[m][atm]++;
					mH[m][atm] = 4-mSigma[m][atm]-mPi[m][atm]-mZ[m][atm];
					if (mH[m][atm] < 0) return cErrorUnexpected;
					}
				}
			}

/*		for (int m=0; m<mRxn.getMolecules(); m++) {			// set aromaticity flags
			StereoMolecule mol = mRxn.getMolecule(m);
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				if (cAtomAromatic & mAtomFlags[m][atm]) continue;
				int retval = chkAroma(mol,atm,aromAtms);
				if (retval != 0) return retval;
				if (aromAtms[0] == 0)
					continue;
				for (int i=1; i<=(31 & aromAtms[0]); i++)
					mAtomFlags[m][aromAtms[i]] |= cAtomAromatic;
				flagAromaticBonds( mol, aromAtms );
				}
			}
		
		
		for (int m=0; m<mRxn.getMolecules(); m++) {			// set allylic flags
			StereoMolecule mol = mRxn.getMolecule(m);
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				for (int i=0; i<getConnAtoms(atm); i++) {
					int connAtm = mol.getConnAtom(atm, i);
					if (mol.getAtomicNo(connAtm) != 6) continue;
					if (mPi[m][connAtm] == 0) continue;
					for (int j=0; j<mol.getConnAtoms(connAtm); j++) {
						int actlAtm = mol.getConnAtom(connAtm, j);
						if (actlAtm == atm) continue;
						if (mol.getAtomicNo(actlAtm) != 6) continue;
						if (mol.getConnBondOrder(connAtm, j) == 1) continue;
						mAtomFlags[m][atm] |= cAtomAllylic;
						goto nextAllylAtm;
						}
					}
		nextAllylAtm:;
				}
			}*/
		
		for (int m=0; m<mRxn.getMolecules(); m++) {  // set double_bond_to_hetero flags
			StereoMolecule mol = mRxn.getMolecule(m);
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				for (int i=0; i<mol.getConnAtoms(atm); i++)	{
					if (mol.getConnBondOrder(atm, i) == 1) continue;
					int connAtm = mol.getConnAtom(atm, i);
					if ((mAtomType[m][connAtm] & 3) != 2) continue;
					mAtomFlags[m][atm] |= cAtomDBondToHetero;
					break;
					}
				}
			}
		
		
/*		for (int m=0; m<mRxn.getMolecules(); m++) {			// set stabilized flags
			StereoMolecule mol = mRxn.getMolecule(m);
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				for (int i=0; i<mol.getConnAtoms(atm); i++) {
					int connAtm = mol.getConnAtom(atm,i);
					if (mol.isAromaticAtom(atm)) continue;
					for (int j=0; j<mol.getConnAtoms(connAtm); j++) {
						int actlAtm = mol.getConnAtom(connAtm, j);
						if (actlAtm == atm) continue;
						if ((mAtomType[m][actlAtm] & 3) != 2) continue;
						if (mol.getConnBondOrder(connAtm, j) == 1) continue;
						mAtomFlags[m][atm] |= cAtomStabilized;
						goto nextStabAtm;
						}
					}
		nextStabAtm:;
				}
			}

		for (int m=0; m<mRxn.getMolecules(); m++) {		// set 3-ring member flags
			StereoMolecule mol = mRxn.getMolecule(m);
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				if (mAtomFlags[m][atm] & cAtomThreeRingMember) continue;
				int ringSize = ringMember( mol, atm, 12, member );
				if (ringSize != 0) {
					if (ringSize == 3)
						for (int i=0; i<3; i++)
							mAtomFlags[m][member[i]] |= cAtomThreeRingMember;
					for (int i=0; i<ringSize; i++)
						mAtomFlags[m][member[i]] |= cAtomRingMember;
					}
				}
			}*/
		
		for (int m=0; m<mRxn.getReactants(); m++) {   // mark all mapped changing atoms
			StereoMolecule mol = mRxn.getMolecule(m);
			int chCount = 0;
			for (int atm=0; atm<mol.getAtoms(); atm++) {
//				if ((mol >= mols) && (mRxn.getAtomMappingNo( mol, atm ) == 0)) continue;
				if (mCorAtom[m][atm] == 255) {
					mAtomFlags[m][atm] |= cAtomChanges + cAtomNoCounterAtom;
					continue;							// no mapping info available
					}

				int proMol = mCorProd[m][atm];
				int proAtm = mCorAtom[m][atm];
		
				if (gErrout != null) {
					try {
						gErrout.write("EduAtom: "+atm+" conntype:");
						for (int i=0; i<4; i++)
							gErrout.write(mConnType[m][atm][i]);
						gErrout.newLine();
						gErrout.write("ProAtom: "+atm+" conntype:");
						for (int i=0; i<4; i++)
							gErrout.write(mConnType[proMol][proAtm][i]);
						gErrout.newLine();
						gErrout.write("EduAtom: "+atm+" connmpno:");
						for (int i=0; i<4; i++)
							gErrout.write(mConnMpNo[m][atm][i]);
						gErrout.newLine();
						gErrout.write("ProAtom: "+atm+" connmpno:");
						for (int i=0; i<4; i++)
							gErrout.write(mConnMpNo[proMol][proAtm][i]);
						gErrout.newLine();
						}
					catch (IOException ioe) {}
					}
		
				boolean found = false;
				for (int i=0; i<MAXCONNS; i++)
						if (mConnType[m][atm][i] != mConnType[proMol][proAtm][i])
							found = true;
				for (int i=0; i<MAXCONNS; i++)
						if (mConnMpNo[m][atm][i] != mConnMpNo[proMol][proAtm][i])
							found = true;
								if (found)
				if (found || mZ[m][atm] != mZ[proMol][proAtm]) {
					int flags = cAtomChanges + cAtomNotClassifiedYet
					+ ( (mPi[m][atm] != mPi[proMol][proAtm]) ? cAtomPiChange : 0)
					+ ( (mZ[m][atm] != mZ[proMol][proAtm]) ? cAtomZChange : 0)
					+ ( (mSigma[m][atm] != mSigma[proMol][proAtm]) ? cAtomSigmaChange : 0);
					mAtomFlags[m][atm] |= flags;
					mAtomFlags[proMol][proAtm] |= flags;
					chCount++;
					}
				}
		
			if (chCount > 5) {				// reset changMsk for nonchanging aromatics
				for (int atm=0; atm<mol.getAtoms(); atm++) {
					if (!mol.isAromaticAtom(atm)) continue;
		
					if ((mAtomFlags[m][atm] & cAtomChanges) == 0) continue;
		
					if ((mAtomFlags[m][atm] & (cAtomSigmaChange + cAtomNoCounterAtom)) != 0)
						continue;	// sigma change or no counter atom
		
					int proM = mCorProd[m][atm];
					int proAtm = mCorAtom[m][atm];
					StereoMolecule proMol = mRxn.getMolecule(proM);
		
					if (!proMol.isAromaticAtom(proAtm)) continue;
		
					boolean found = false;
					for (int i=0; i<MAXCONNS; i++)
						if (mConnCMNo[m][atm][i] != mConnCMNo[proM][proAtm][i])
							found = true;
					if (found)
						continue;
		
					if (mPi[m][atm]+mZ[m][atm]
					 != mPi[proM][proAtm]+mZ[proM][proAtm]) continue;
		
							// if at least one bond to an atom, which is non-aromatic
							// either in educt or in product or in both, is changing
							// its order then keep atm marked as changing
					found = false;
					for (int i=0; i<mol.getConnAtoms(atm); i++) {
						int connAtm = mol.getConnAtom(atm, i);
						int proConnAtm = mCorAtom[m][connAtm];
						if (proConnAtm == 255) continue;
						if (mol.isAromaticAtom(connAtm))
							if (mRxn.getMolecule(mCorProd[m][connAtm]).isAromaticAtom(proConnAtm))
								continue;
						for (int j=0; j<proMol.getConnAtoms(proAtm); j++)
							if (proConnAtm == proMol.getConnAtom(proAtm, j))
								if (mol.getBondOrder(mol.getConnBond(atm,i))
								 != proMol.getBondOrder(proMol.getConnBond(proAtm, j)))
									{ found = true; break; }
						if (found) break;
						}
					if (found) continue;
		
					int flags = cAtomChanges + cAtomNotClassifiedYet
						  + cAtomPiChange + cAtomZChange + cAtomSigmaChange;
					mAtomFlags[m][atm] &= ~flags;
					mAtomFlags[proM][proAtm] &= ~flags;
					}
				}
			}
		
		
		if (gErrout != null) {
			for (int m=0; m<mRxn.getMolecules(); m++) {
				try {
					gErrout.newLine();
					gErrout.write("***** mol: "+m+" *****");
					gErrout.newLine();
					gErrout.write("!!!! atm  a# pro  m# cha sig  pi   h   z cAt cPr typ !!!!");
					gErrout.newLine();
					StereoMolecule mol = mRxn.getMolecule(m);
					for (int atm=0; atm<mol.getAtoms(); atm++) {
						gErrout.write(String.format("	%4hd%4hd%4hd%4hd%4hd%4hd%4hd%4hd%4hd%4hd%4hd%4hd",
							atm,mol.getAtomicNo(atm),mol.getAtomCharge(atm),mol.getAtomMapNo(atm),
							mAtomFlags[m][atm],mSigma[m][atm],mPi[m][atm],mH[m][atm],mZ[m][atm],
							mCorAtom[m][atm],mCorProd[m][atm],mAtomType[m][atm]));
						gErrout.newLine();
						}
					gErrout.write("!!!! atm  a# +/-  m# cha sig  pi   h   z cAt cPr typ !!!!");
					gErrout.newLine();
					}
			catch (IOException ioe) {}
				}
			}
		
		int ccBndFrm = 0;
		int ccBndClv = 0;				// determine gross reaction devision
		for (int m=0; m<mRxn.getReactants(); m++) {
			StereoMolecule mol = mRxn.getMolecule(m);
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				if ((mAtomFlags[m][atm] & cAtomChanges) == 0) continue;
				if ((mAtomFlags[m][atm] & cAtomNoCounterAtom) != 0) continue;
				if (mAtomType[m][atm] != CARBON) continue;
				int proM = mCorProd[m][atm];
				int proAtm = mCorAtom[m][atm];

				int deltaS = 0;
				for (int i=0; i<MAXCONNS; i++) if (mConnCMNo[proM][proAtm][i] != 0) deltaS++;
				for (int i=0; i<MAXCONNS; i++) if (mConnCMNo[m][atm][i] != 0) deltaS--;
				if (deltaS != mSigma[proM][proAtm]-mSigma[m][atm]) return cErrorUnMappedCInConOrFr;
		
				for (int i=0; i<MAXCONNS; i++) {							// C-C Bond cleavage
					if (mConnCMNo[m][atm][i] != 0) {
						boolean found = false;
						for (int j=0; j<MAXCONNS; j++) {
							if (mConnCMNo[m][atm][i]
								== mConnCMNo[proM][proAtm][j]) {
								found = true;
								break;
								}
							}
						if (!found) {
							int oppositeAtom = -1;
							for (int l=0; l<mol.getAtoms(); l++) {
								if (mol.getAtomMapNo(l)	== mConnCMNo[m][atm][i]) {
									oppositeAtom = l;
									break;
									}
								}
		
							if (oppositeAtom > atm) {
								if (ccBndClv == 4) return cErrorCCBondCleavageLimit;
								ccBCmol[ccBndClv] = m;	// more than 2 CC-bonds cleaved
								ccBCatm[ccBndClv] = atm;
								ccBndClv++;
								ccBCmol[ccBndClv] = m;
								ccBCatm[ccBndClv] = oppositeAtom;
								ccBndClv++;
								}
							}
						}
					}
		
				for (int i=0; i<MAXCONNS; i++) {						// C-C Bond formation
					if (mConnCMNo[proM][proAtm][i] != 0) {
						boolean found = false;
						for (int j=0; j<MAXCONNS; j++) {
							if (mConnCMNo[proM][proAtm][i] == mConnCMNo[m][atm][j]) {
								found = true;
								break;
								}
							}
						if (!found) {
							int oppositeMol = -1;
							int oppositeAtom = -1;
							for (int k=m; k<mRxn.getMolecules(); k++) {
								if ((k >= mRxn.getReactants())) continue;
								StereoMolecule kMol = mRxn.getMolecule(k);
								for (int l=0; l<kMol.getAtoms(); l++) {
									if (kMol.getAtomMapNo(l) == mConnCMNo[proM][proAtm][i]) {
										oppositeMol = k;
										oppositeAtom = l;
										break;
										}
									}
								}
							if ((oppositeMol > m)
							|| ((oppositeMol == m) && (oppositeAtom > atm))) {
								if (ccBndFrm == 4) return cErrorCCBondCreationLimit;
								ccBFmol[ccBndFrm] = m;	// more than 2 new C-C bonds
								ccBFatm[ccBndFrm] = atm;
								ccBndFrm++;
								ccBFmol[ccBndFrm] = oppositeMol;
								ccBFatm[ccBndFrm] = oppositeAtom;
								ccBndFrm++;
								}
							}
						}
					}
				}
			}
		
		if (gErrout != null) {
			try {
				gErrout.write("Form:"+ccBndFrm+", Clv:"+ccBndClv+"; ccBFmols:"+ccBFmol[0]+","+ccBFmol[1]+"; ccBFatms:"+ccBFatm[0]+","+ccBFatm[0]);
				gErrout.newLine();
				}
			catch (IOException ioe) {}
			}
		
		if (ccBndFrm == 2 && ccBndClv == 0)
			{
		//	gSC_no++;
			int retval = snglCnst(ccBFmol,ccBFatm);
			if (retval != 0) return retval;
		//	gCSC_no++;
			}
		else if (ccBndFrm == 4 && ccBndClv == 0)
			{
		//	gDC_no++;
			int retval = dblCnst(ccBFmol,ccBFatm);
			if (retval != 0) return retval;
		//	gCDC_no++;
			}
		else if (ccBndFrm == 0 && ccBndClv == 2)
			{
		//	gSF_no++;
			int retval = snglFrgm(ccBCmol,ccBCatm);
			if (retval != 0) return retval;
		//	gCSF_no++;
			}
		else if (ccBndFrm == 0 && ccBndClv == 4)
			{
		//	gDF_no++;
			int retval = dblFrgm(ccBCmol,ccBCatm);
			if (retval != 0) return retval;
		//	gCDF_no++;
			}
		else if (ccBndFrm == 2 && ccBndClv == 2)
			{
		//	gRA_no++;
			int retval = rearrang(ccBFmol,ccBFatm,ccBCmol,ccBCatm);
			if (retval != 0) return retval;
		//	gCRA_no++;
			}
		else if (ccBndFrm != 0 || ccBndClv != 0)
			{
			if (gErrout != null) {
				try {
					gErrout.write("Form:"+ccBndFrm+", Clv:"+ccBndClv);
					gErrout.newLine();
					}
				catch (IOException ioe) {}
				}
			mResult.mClassName = "not classified";
			return cErrorComplexReaction;
			}
		
		//if (ccBndFrm == 0 && ccBndClv == 0) gRE_no++;
		int retval = refuncs();
		if (retval != 0) return retval;

		if (mUnitRxn == 0)	// no skeletal nor refunc unit reactions found
			if (mIndexToCreate != cIndexFullPermutation)
				return cErrorNoChangingAtoms;
		
		//if (ccBndFrm == 0 && ccBndClv == 0) gCRE_no++;
		
		mResult.mUnitRxns = mUnitRxn;
		findClassNam();
		markReactionCenters();
		
		return cErrorNoError;
		}


	private int snglCnst(int[] ccBFmol, int[] ccBFatm) {
		int[] strLngth = new int[2];
		int[][] strand = new int[2][8];		// strand[atm][0]: molecule
		int atm;				// strand[atm][1-8]: alpha - delta
		
		for (atm=0; atm<2; atm++)	// get atom no's of changing strand
			{
			strand[atm][0] = ccBFmol[atm];
			strand[atm][1] = ccBFatm[atm];
			strLngth[atm] = findStrand(strand[atm]);

			if (gErrout != null)
				try { gErrout.newLine(); } catch (IOException ioe) {}

			if (strLngth[atm] < 1) return cErrorForkedOrLongStrand;	// branched strand
			}
		return oneConst(strand[0],strand[1],strLngth[0],strLngth[1]);
		}


	/**
     * manage one Construction Reaction with known strands
	 * @return
	 */
	private int oneConst( int[] strand1, int[] strand2, int strLen1, int strLen2 ) {
		int[] indexdata = new int[3];
														// electrocyclic ring closure
		if (strand1[0] == strand2[0] && strand1[1] == strand2[strLen2])
			return ringClosure(strand1,strand2,strLen1,strLen2);
		
		int retval = bndSpcDat(strand1,strand2,2,indexdata,0);
		if (retval != 0) return retval;

		// z-pi-Lists of entire strand
		int deltaLst = getDelta(strand1,strLen1,7);	// classify first half construction
		retval = clssfyHR(indexdata,1,deltaLst,strand1,strLen1,true,6);
		if (retval != 0) return retval;						// classification error
		
		deltaLst = getDelta(strand2,strLen2,7);   // classify second half construction
		retval = clssfyHR(indexdata,2,deltaLst,strand2,strLen2,true,7);
		if (retval != 0) return retval;						// classification error
		
		mResult.mMainClass[mUnitRxn] = S_CONST;

		if (mIndexToCreate != cIndexNone) {
			DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

			putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
			putindexbits(' ',BONDBITS,indexdata[0]);			// write rxn index entry
			putindexbits(' ',HRXNBITS,indexdata[1]);
			putindexbits('l',HRXNBITS,indexdata[2]);
			
	//		if (mBatchClassifierP)
	//			mBatchClassifierP->incIndexEntries( S_CONST );
			
			if (mIndexToCreate == cIndexFullPermutation) {		// supply index for query
				putindexbits('c',32,dbrxn.getReactionRegNo());
				putindexbits(' ',BONDBITS,indexdata[0]);
				putindexbits(' ',HRXNBITS,indexdata[2]);
				putindexbits('l',HRXNBITS,indexdata[1]);
				}
			}

		mUnitRxn++;
		return cErrorNoError;
		}


	/**
	 * manage an electrocyclic ring closure reaction
	 * @return
	 */
	private int ringClosure( int[] strand1, int[] strand2, int strLen1, int strLen2 ) {
		int[] indexdata = new int[4];
		
		int retval = bndSpcDat(strand1,strand2,2,indexdata,0);
		if (retval != 0) return retval;
		
		int deltaLst = getDelta(strand1,strLen1,7);
		int tempDeltaLst = getDelta(strand2,strLen2,7);
		if (tempDeltaLst > deltaLst)
			deltaLst = tempDeltaLst;

		retval = clssfyCR(indexdata,1,deltaLst,strand1,strLen1,true);
		if (retval != 0) return retval;						// classification error

		indexdata[2] = properties(strand1[0],strand1[1]);
		indexdata[3] = properties(strand2[0],strand2[1]);
		
		mResult.mFlashMol[mUnitRxn][6] = strand1[0];
		mResult.mFlashAtom[mUnitRxn][6] = strand1[1];
		mResult.mFlashMol[mUnitRxn][7] = strand2[0];
		mResult.mFlashAtom[mUnitRxn][7] = strand2[1];
		mResult.mMainClass[mUnitRxn] = E_RING_C;
		
		if (mIndexToCreate != cIndexNone) {
			DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

			putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
			putindexbits(' ',BONDBITS,indexdata[0]);			// write rxn index entry
			putindexbits(' ',CRXNBITS,indexdata[1]);
			putindexbits(' ',PROPBITS,indexdata[2]);
			putindexbits('l',PROPBITS,indexdata[3]);
			
	//		if (mBatchClassifierP)
	//			mBatchClassifierP->incIndexEntries( E_RING_C );
			
			if (mIndexToCreate == cIndexFullPermutation)
				{
				putindexbits('c',32,dbrxn.getReactionRegNo());
				putindexbits(' ',BONDBITS,indexdata[0]);
				putindexbits(' ',CRXNBITS,indexdata[1]);
				putindexbits(' ',PROPBITS,indexdata[3]);
				putindexbits('l',PROPBITS,indexdata[2]);
				}
			}

		mUnitRxn++;
		
		return cErrorNoError;
		}


	/**
	 * calculate z-pi-list from strand no's
	 */
	private int getDelta(int[] strand, int strLngth, int maxLngth) {
		int eduZ,proZ;
		int i,mol;
		int eduLst,proLst;
		
		eduLst = proLst = 0;
		mol = strand[0];

		for (i=1; i<=maxLngth; i++)
			{
			eduLst <<= 4;
			proLst <<= 4;
			if (i>strLngth) continue;
			eduZ = mZ[mol][strand[i]];
			proZ = mZ[mCorProd[mol][strand[i]]][mCorAtom[mol][strand[i]]];
			if ((eduZ == 4) || (proZ == 4))
				{
				eduZ -= 1;
				proZ -= 1;
				}
			eduLst += (eduZ<<2) + mPi[mol][strand[i]];
			proLst += (proZ<<2) + mPi[mCorProd[mol][strand[i]]][mCorAtom[mol][strand[i]]];
			}

		if (gErrout != null) {
			try {
				gErrout.write(String.format("eduLst:%8x, proLst:%8x, deltaLst:%8x\n", eduLst,proLst,eduLst-proLst));
				gErrout.newLine();
				}
			catch (IOException ioe) {}
			}

		return eduLst-proLst;
		}


	/**
	 * double Construction Reaction
	 * @return
	 */
	private int dblCnst(int[] ccBFmol, int[] ccBFatm) {
		int[] indexdata = new int[6];
		int[] strLngth = new int[4];
		int strnd;
		int hlfRxn1,hlfRxn2,depRxn;
		int[] dependnc = new int[2];	// indicates for each educt atom of new formed
										// bond, if this an independent half reaction (0)
										// or if the appropriate strand ends with an atom
										// from that the second bond is formed (atm no)
		int[][] strand = new int[4][8];	// strand[atm][0]: molecule
										// strand[atm][1-7]: alpha - delta
		
		for (int atm=0; atm<4; atm++) {		// get atom no's of changing strands
			strand[atm][0] = ccBFmol[atm];
			strand[atm][1] = ccBFatm[atm];
			strLngth[atm] = findStrand(strand[atm]);

			if (gErrout != null)
				try { gErrout.newLine(); } catch (IOException ioe) {}

			if (strLngth[atm] < 1) return cErrorForkedOrLongStrand;	// branched strand
			}
		
		for (int i=0; i<2; i++)
			dependnc[i] =
		 	+ (strand[i][0] == strand[2][0] && strand[i][1] == strand[2][strLngth[2]] ? 2 : 0)
		 	+ (strand[i][0] == strand[3][0] && strand[i][1] == strand[3][strLngth[3]] ? 3 : 0);

		if ((dependnc[0] == 0) && (dependnc[1] == 0)) {
			int retval = oneConst(strand[0],strand[1],strLngth[0],strLngth[1]);
			if (retval != 0) return retval;
			return oneConst(strand[2],strand[3],strLngth[2],strLngth[3]);
			}
		
		if ((dependnc[0] == 0) || (dependnc[1] == 0)) {
			if (dependnc[0] == 0) {
				hlfRxn1 = 0;
				hlfRxn2 = (dependnc[1] == 2) ? 3 : 2;
				depRxn = 1;
				}
			else
				{
				hlfRxn1 = 1;
				hlfRxn2 = (dependnc[0] == 2) ? 3 : 2;
				depRxn = 0;
				}
		
			if (strand[hlfRxn1][0] == strand[depRxn][0])
				if (strand[hlfRxn1][1] == strand[depRxn][strLngth[depRxn]])
					return cError2AlphasInSameStrand;		// both alpha carbons share one strand
			if (strand[hlfRxn2][0] == strand[5-hlfRxn2][0])
				if (strand[hlfRxn2][1] == strand[5-hlfRxn2][strLngth[5-hlfRxn2]])
					return cError2AlphasInSameStrand;		// both alpha carbons share one strand
		
			int retval = bndSpcDat(strand[hlfRxn1],strand[depRxn],0,indexdata,0);
			if (retval != 0) return retval;
			retval = bndSpcDat(strand[hlfRxn2],strand[5-hlfRxn2],2,indexdata,1);
			if (retval != 0) return retval;
		
			int deltaLst = getDelta(strand[depRxn],strLngth[depRxn],strLngth[depRxn]);
			int tempDeltaLst = getDelta(strand[5-hlfRxn2],strLngth[5-hlfRxn2],strLngth[5-hlfRxn2]);
		
			if (tempDeltaLst > deltaLst)	// use strand with higher delta z-pi-list
				{
				deltaLst = tempDeltaLst;
				depRxn = 5-hlfRxn2;
				}
		
												// two end hlf rxn
			retval = clssfyDR(indexdata,2,deltaLst,true,strand[depRxn],strLngth[depRxn],4);
			if (retval != 0) return retval;
		
			deltaLst = getDelta(strand[hlfRxn1],strLngth[hlfRxn1],7);			// standard hlf rxn
			retval = clssfyHR(indexdata,4,deltaLst,strand[hlfRxn1],strLngth[hlfRxn1],true,6);
			if (retval != 0) return retval;						// classification error
		
			deltaLst = getDelta(strand[hlfRxn2],strLngth[hlfRxn2],7);		// standard hlf rxn
			retval = clssfyHR(indexdata,5,deltaLst,strand[hlfRxn2],strLngth[hlfRxn2],true,7);
			if (retval != 0) return retval;						// classification error
		
			mResult.mMainClass[mUnitRxn] = HD_D_CONST;

			if (mIndexToCreate != cIndexNone) {
				DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

				putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
				putindexbits(' ',BONDBITS,indexdata[0]);		// write rxn index entry
				putindexbits(' ',BONDBITS,indexdata[1]);
				putindexbits(' ',DRXNBITS,indexdata[2]);
				putindexbits(' ',HRXNBITS,indexdata[4]);
				putindexbits('l',HRXNBITS,indexdata[5]);
			
	//			if (mBatchClassifierP)
	//				mBatchClassifierP->incIndexEntries( HD_D_CONST );
			
				if (mIndexToCreate == cIndexFullPermutation) {
					putindexbits('c',32,dbrxn.getReactionRegNo());
					putindexbits(' ',BONDBITS,indexdata[1]);
					putindexbits(' ',BONDBITS,indexdata[0]);
					putindexbits(' ',DRXNBITS,indexdata[3]);
					putindexbits(' ',HRXNBITS,indexdata[5]);
					putindexbits('l',HRXNBITS,indexdata[4]);
					}
				}
		
			mUnitRxn++;
			for (int i=1; i<=strLngth[depRxn]; i++)
				mAtomFlags[strand[depRxn][0]][strand[depRxn][i]] &= ~cAtomNotClassifiedYet;
			}
		else
			{
			if (strand[0][0] == strand[1][0])
				if (strand[0][1] == strand[1][strLngth[1]])
					return cError2AlphasInSameStrand;		// both alpha carbons share one strand
			if (strand[2][0] == strand[3][0])
				if (strand[2][1] == strand[3][strLngth[3]])
					return cError2AlphasInSameStrand;		// both alpha carbons share one strand
		
			int retval = bndSpcDat(strand[0],strand[1],0,indexdata,0);
			if (retval != 0) return retval;
			retval = bndSpcDat(strand[2],strand[3],2,indexdata,1);
			if (retval != 0) return retval;
		
			for (strnd=0; strnd<2; strnd++)
				{
				int deltaLst = getDelta(strand[strnd],strLngth[strnd],strLngth[strnd]);
				int tempDeltaLst = getDelta(strand[dependnc[strnd]],strLngth[dependnc[strnd]],
						strLngth[dependnc[strnd]]);
		
				if (tempDeltaLst > deltaLst)
					{
					deltaLst = tempDeltaLst;
					depRxn = dependnc[strnd];
					}
				else depRxn = strnd;
		
				retval = clssfyDR(indexdata,2+2*strnd,deltaLst,true,strand[depRxn],strLngth[depRxn],strnd*2+4);
				if (retval != 0) return retval;
		
				for (int i=1; i<=strLngth[strnd]; i++)
					mAtomFlags[strand[strnd][0]][strand[strnd][i]] &= ~cAtomNotClassifiedYet;
				}
		
			mResult.mMainClass[mUnitRxn] = D_D_CONST;

			if (mIndexToCreate != cIndexNone) {
				DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

				putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
				putindexbits(' ',BONDBITS,indexdata[0]);		// write rxn index entry
				putindexbits(' ',BONDBITS,indexdata[1]);
				putindexbits(' ',DRXNBITS,indexdata[2]);
				putindexbits('l',DRXNBITS,indexdata[4]);
			
	//			if (mBatchClassifierP)
	//				mBatchClassifierP->incIndexEntries( D_D_CONST );
			
				if (mIndexToCreate == cIndexFullPermutation) {
					putindexbits('c',32,dbrxn.getReactionRegNo());
					putindexbits(' ',BONDBITS,indexdata[0]);
					putindexbits(' ',BONDBITS,indexdata[1]);
					putindexbits(' ',DRXNBITS,indexdata[4]);
					putindexbits('l',DRXNBITS,indexdata[2]);
					putindexbits('c',32,dbrxn.getReactionRegNo());
					putindexbits(' ',BONDBITS,indexdata[1]);
					putindexbits(' ',BONDBITS,indexdata[0]);
					putindexbits(' ',DRXNBITS,indexdata[3]);
					putindexbits('l',DRXNBITS,indexdata[5]);
					putindexbits('c',32,dbrxn.getReactionRegNo());
					putindexbits(' ',BONDBITS,indexdata[1]);
					putindexbits(' ',BONDBITS,indexdata[0]);
					putindexbits(' ',DRXNBITS,indexdata[5]);
					putindexbits('l',DRXNBITS,indexdata[3]);
					}
				}

			mUnitRxn++;
			}
		return cErrorNoError;
		}


	/**
	 * single fragmentation Reaction
	 * @return
	 */
	private int snglFrgm(int[] ccBCmol, int[] ccBCatm) {
		int[][] strand = new int[2][8];
		int[] strLngth = new int[2];
		
		for (int atm=0; atm<2; atm++) {			// get atom no's of changing strand
			strand[atm][0] = ccBCmol[atm];
			strand[atm][1] = ccBCatm[atm];
			strLngth[atm] = findStrand(strand[atm]);

			if (gErrout != null)
				try { gErrout.newLine(); } catch (IOException ioe) {}

			if (strLngth[atm] < 1) return cErrorForkedOrLongStrand;	// branched strand
			}
		return oneFrgm(strand[0],strand[1],strLngth[0],strLngth[1]);
		}


	/**
	 * handle one fragmentation Reaction
	 * @return
	 */
	private int oneFrgm(int[] strand1, int[] strand2, int strLen1, int strLen2) {
		int[] indexdata = new int[3];
		
														// electrocyclic ring opening
		if (strand1[0] == strand2[0] && strand1[1] == strand2[strLen2])
			return ringOpening(strand1,strand2,strLen1,strLen2);
		
		int retval = bndSpcDat(strand1,strand2,2,indexdata,0);
		if (retval != 0) return retval;
		
		int deltaLst = getDelta(strand1,strLen1,7);   // classify first half fragmentation
		retval = clssfyHR(indexdata,1,-deltaLst,strand1,strLen1,false,6);
		if (retval != 0) return retval;						// classification error
		
		deltaLst = getDelta(strand2,strLen2,7);  // classify second half fragmentation
		retval = clssfyHR(indexdata,2,-deltaLst,strand2,strLen2,false,7);
		if (retval != 0) return retval;						// classification error
		
		mResult.mMainClass[mUnitRxn] = S_FRAG;

		if (mIndexToCreate != cIndexNone) {
			DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

			putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
			putindexbits(' ',BONDBITS,indexdata[0]);			// write rxn index entry
			putindexbits(' ',HRXNBITS,indexdata[1]);
			putindexbits('l',HRXNBITS,indexdata[2]);
			
	//		if (mBatchClassifierP)
	//			mBatchClassifierP->incIndexEntries( S_FRAG );
			
			if (mIndexToCreate == cIndexFullPermutation) {
				putindexbits('c',32,dbrxn.getReactionRegNo());
				putindexbits(' ',BONDBITS,indexdata[0]);
				putindexbits(' ',HRXNBITS,indexdata[2]);
				putindexbits('l',HRXNBITS,indexdata[1]);
				}
			}

		mUnitRxn++;
		return cErrorNoError;
		}


	/**
	 * manage an electrocyclic ring opening reaction
	 * @return
	 */
	private int ringOpening(int[] strand1, int[] strand2, int strLen1, int strLen2) {
		int[] indexdata = new int[4];
		
		int retval = bndSpcDat(strand1,strand2,2,indexdata,0);
		if (retval != 0) return retval;
		
		int deltaLst = getDelta(strand1,strLen1,7);
		int tempDeltaLst = getDelta(strand2,strLen2,7);
		if (tempDeltaLst < deltaLst)
			deltaLst = tempDeltaLst;
		
		retval = clssfyCR(indexdata,1,-deltaLst,strand1,strLen1,false);
		if (retval != 0) return retval;						// classification error
		
		indexdata[2] = properties(strand1[0],strand1[1]);
		indexdata[3] = properties(strand2[0],strand2[1]);
		
		mResult.mFlashMol[mUnitRxn][6] = strand1[0];
		mResult.mFlashAtom[mUnitRxn][6] = strand1[1];
		mResult.mFlashMol[mUnitRxn][7] = strand2[0];
		mResult.mFlashAtom[mUnitRxn][7] = strand2[1];
		
		mResult.mMainClass[mUnitRxn] = E_RING_O;

		if (mIndexToCreate != cIndexNone) {
			DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

			putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
			putindexbits(' ',BONDBITS,indexdata[0]);			// write rxn index entry
			putindexbits(' ',CRXNBITS,indexdata[1]);
			putindexbits(' ',PROPBITS,indexdata[2]);
			putindexbits('l',PROPBITS,indexdata[3]);
	
	//		if (mBatchClassifierP)
	//			mBatchClassifierP->incIndexEntries( E_RING_O );
			
			if (mIndexToCreate == cIndexFullPermutation)
				{
				putindexbits('c',32,dbrxn.getReactionRegNo());
				putindexbits(' ',BONDBITS,indexdata[0]);
				putindexbits(' ',CRXNBITS,indexdata[1]);
				putindexbits(' ',PROPBITS,indexdata[3]);
				putindexbits('l',PROPBITS,indexdata[2]);
				}
			}

		mUnitRxn++;
		
		return cErrorNoError;
		}


	/**
	 * double Fragmentation Reaction
	 * @return
	 */
	private int dblFrgm(int[] ccBCmol, int[] ccBCatm) {
		int[] indexdata = new int[6];
		int[] strLngth = new int[4];
		int hlfRxn1,hlfRxn2,depRxn;
		int[] dependnc = new int[2];	// indicates for each educt atom of new formed
										// bond, if this an independent half reaction (0)
										// or if the appropriate strand ends with an atom
										// from that the second bond is formed (atm no)

		int[][] strand = new int[4][8];	// strand[atm][0]: molecule
										// strand[atm][1-8]: alpha - delta

		for (int atm=0; atm<4; atm++) {		// get atom no's of changing strands
			strand[atm][0] = ccBCmol[atm];
			strand[atm][1] = ccBCatm[atm];
			strLngth[atm] = findStrand(strand[atm]);

			if (gErrout != null)
				try { gErrout.newLine(); } catch (IOException ioe) {}

			if (strLngth[atm] < 1) return cErrorForkedOrLongStrand;	// branched strand
			}
		
		for (int i=0; i<2; i++)
			dependnc[i] =
			+ (strand[i][0] == strand[2][0] && strand[i][1] == strand[2][strLngth[2]] ? 2 : 0)
			+ (strand[i][0] == strand[3][0] && strand[i][1] == strand[3][strLngth[3]] ? 3 : 0);
		
		if ((dependnc[0] == 0) && (dependnc[1] == 0)) {
			int retval = oneFrgm(strand[0],strand[1],strLngth[0],strLngth[1]);
			if (retval != 0) return retval;
			return oneFrgm(strand[2],strand[3],strLngth[2],strLngth[3]);
			}
		
		if ((dependnc[0] == 0) || (dependnc[1] == 0))
			{
			if (dependnc[0] == 0)
				{
				hlfRxn1 = 0;
				hlfRxn2 = (dependnc[1] == 2) ? 3 : 2;
				depRxn = 1;
				}
			else
				{
				hlfRxn1 = 1;
				hlfRxn2 = (dependnc[0] == 2) ? 3 : 2;
				depRxn = 0;
				}
		
			if (strand[hlfRxn1][0] == strand[depRxn][0])
				if (strand[hlfRxn1][1] == strand[depRxn][strLngth[depRxn]])
					return cError2AlphasInSameStrand;		// both alpha carbons share one strand
			if (strand[hlfRxn2][0] == strand[5-hlfRxn2][0])
				if (strand[hlfRxn2][1] == strand[5-hlfRxn2][strLngth[5-hlfRxn2]])
					return cError2AlphasInSameStrand;		// both alpha carbons share one strand
		
			int retval = bndSpcDat(strand[hlfRxn1],strand[depRxn],0,indexdata,0);
			if (retval != 0) return retval;
			retval = bndSpcDat(strand[hlfRxn2],strand[5-hlfRxn2],2,indexdata,1);
			if (retval != 0) return retval;
		
			int deltaLst = getDelta(strand[depRxn],strLngth[depRxn],strLngth[depRxn]);
			int tempDeltaLst = getDelta(strand[5-hlfRxn2],strLngth[5-hlfRxn2],strLngth[5-hlfRxn2]);
		
			if (tempDeltaLst < deltaLst) {  // lower delta z-pi-list here gives higher after inversion
				deltaLst = tempDeltaLst;
				depRxn = 5-hlfRxn2;
				}
		
															// two end hlf rxn
			retval = clssfyDR(indexdata,2,-deltaLst,false,strand[depRxn],strLngth[depRxn],4);
			if (retval != 0) return retval;
		
			deltaLst = getDelta(strand[hlfRxn1],strLngth[hlfRxn1],7);		// standard hlf rxn
			retval = clssfyHR(indexdata,4,-deltaLst,strand[hlfRxn1],strLngth[hlfRxn1],false,6);
			if (retval != 0) return retval;						// classification error
		
			deltaLst = getDelta(strand[hlfRxn2],strLngth[hlfRxn2],7);		// standard hlf rxn
			retval = clssfyHR(indexdata,5,-deltaLst,strand[hlfRxn2],strLngth[hlfRxn2],false,7);
			if (retval != 0) return retval;						// classification error
		
			mResult.mMainClass[mUnitRxn] = HD_D_FRAG;

			if (mIndexToCreate != cIndexNone) {
				DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

				putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
				putindexbits(' ',BONDBITS,indexdata[0]);		// write rxn index entry
				putindexbits(' ',BONDBITS,indexdata[1]);
				putindexbits(' ',DRXNBITS,indexdata[2]);
				putindexbits(' ',HRXNBITS,indexdata[4]);
				putindexbits('l',HRXNBITS,indexdata[5]);
			
	//			if (mBatchClassifierP)
	//				mBatchClassifierP->incIndexEntries( HD_D_FRAG );
			
				if (mIndexToCreate == cIndexFullPermutation)
					{
					putindexbits('c',32,dbrxn.getReactionRegNo());
					putindexbits(' ',BONDBITS,indexdata[1]);
					putindexbits(' ',BONDBITS,indexdata[0]);
					putindexbits(' ',DRXNBITS,indexdata[3]);
					putindexbits(' ',HRXNBITS,indexdata[5]);
					putindexbits('l',HRXNBITS,indexdata[4]);
					}
				}
		
			mUnitRxn++;
			for (int i=1; i<=strLngth[depRxn]; i++)
				mAtomFlags[strand[depRxn][0]][strand[depRxn][i]] &= ~cAtomNotClassifiedYet;
			}
		else
			{
			if (strand[0][0] == strand[1][0])
				if (strand[0][1] == strand[1][strLngth[1]])
					return cError2AlphasInSameStrand;		// both alpha carbons share one strand
			if (strand[2][0] == strand[3][0])
				if (strand[2][1] == strand[3][strLngth[3]])
					return cError2AlphasInSameStrand;		// both alpha carbons share one strand
		
			int retval = bndSpcDat(strand[0],strand[1],0,indexdata,0);
			if (retval != 0) return retval;
			retval = bndSpcDat(strand[2],strand[3],2,indexdata,1);
			if (retval != 0) return retval;
		
			for (int strnd=0; strnd<2; strnd++) {
				int deltaLst = getDelta(strand[strnd],strLngth[strnd],strLngth[strnd]);
				int tempDeltaLst = getDelta(strand[dependnc[strnd]],strLngth[dependnc[strnd]],strLngth[dependnc[strnd]]);
		
				if (tempDeltaLst < deltaLst) {
					deltaLst = tempDeltaLst;
					depRxn = dependnc[strnd];
					}
				else depRxn = strnd;
		
				retval = clssfyDR(indexdata,2+2*strnd,-deltaLst,false,strand[depRxn],strLngth[depRxn],strnd*2+4);
				if (retval != 0) return retval;
		
				for (int i=1; i<=strLngth[strnd]; i++)
					mAtomFlags[strand[strnd][0]][strand[strnd][i]] &= ~cAtomNotClassifiedYet;
				}
		
			mResult.mMainClass[mUnitRxn] = D_D_FRAG;

			if (mIndexToCreate != cIndexNone) {
				DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

				putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
				putindexbits(' ',BONDBITS,indexdata[0]);		// write rxn index entry
				putindexbits(' ',BONDBITS,indexdata[1]);
				putindexbits(' ',DRXNBITS,indexdata[2]);
				putindexbits('l',DRXNBITS,indexdata[4]);
			
	//			if (mBatchClassifierP)
	//				mBatchClassifierP->incIndexEntries( D_D_FRAG );
			
				if (mIndexToCreate == cIndexFullPermutation)
					{
					putindexbits('c',32,dbrxn.getReactionRegNo());
					putindexbits(' ',BONDBITS,indexdata[0]);
					putindexbits(' ',BONDBITS,indexdata[1]);
					putindexbits(' ',DRXNBITS,indexdata[4]);
					putindexbits('l',DRXNBITS,indexdata[2]);
					putindexbits('c',32,dbrxn.getReactionRegNo());
					putindexbits(' ',BONDBITS,indexdata[1]);
					putindexbits(' ',BONDBITS,indexdata[0]);
					putindexbits(' ',DRXNBITS,indexdata[3]);
					putindexbits('l',DRXNBITS,indexdata[5]);
					putindexbits('c',32,dbrxn.getReactionRegNo());
					putindexbits(' ',BONDBITS,indexdata[1]);
					putindexbits(' ',BONDBITS,indexdata[0]);
					putindexbits(' ',DRXNBITS,indexdata[5]);
					putindexbits('l',DRXNBITS,indexdata[3]);
					}
				}

			mUnitRxn++;
			}
		return cErrorNoError;
		}


	/**
	 * rearrangement
	 * @return
	 */
	private int rearrang(int[] ccBFmol, int[] ccBFatm, int[] ccBCmol, int[] ccBCatm) {
		int[] strLngth = new int[4];
		int[] indexdata = new int[5];
		int[] dependnc = new int[2];// indicates for each educt atom of new formed
									// bond, if this an independent half reaction (0)
									// or if the appropriate strand ends with an atom
									// from that the second bond is formed (atm no)

		int[][] strand = new int[4][8];				// strand[atm][0]: molecule
													// strand[atm][1-8]: alpha - delta
		for (int i=0; i<2; i++) {
			strand[i][0] = ccBFmol[i];
			strand[i][1] = ccBFatm[i];
			strLngth[i] = findStrand(strand[i]);
			if (strLngth[i] < 1) return cErrorForkedOrLongStrand;		// branched strand
		
			int strEnd = strand[i][strLngth[i]];
			for (int j=0; j<2; j++)
				if (strEnd == ccBCatm[j] && ccBFmol[i] == ccBCmol[j])
					dependnc[i] = j+2;
			}
		
		for (int i=0; i<2; i++) {
			strand[i+2][0] = ccBCmol[i];
			strand[i+2][1] = ccBCatm[i];
			strLngth[i+2] = findStrand(strand[i+2]);
			if (strLngth[i+2] < 1) return cErrorForkedOrLongStrand;	// branched strand
			}
		
		if ((dependnc[0] == 0) && (dependnc[1] == 0)) {
			int retval = oneConst(strand[0],strand[1],strLngth[0],strLngth[1]);
			if (retval != 0) return retval;
			return oneFrgm(strand[2],strand[3],strLngth[2],strLngth[3]);
			}
		
		if ((dependnc[0] == 0) || (dependnc[1] == 0)) {
			int hlfCon,hlfFrg,depRxn;
			if (dependnc[0] == 0) {
				hlfCon = 0;
				hlfFrg = (dependnc[1] == 2) ? 3 : 2;
				depRxn = 1;
				}
			else
				{
				hlfCon = 1;
				hlfFrg = (dependnc[0] == 2) ? 3 : 2;
				depRxn = 0;
				}
		
			mResult.mRearStrandLen[0] = strLngth[depRxn];

			int retval = bndSpcDat(strand[hlfCon],strand[depRxn],0,indexdata,0);
			if (retval != 0) return retval;
			retval = bndSpcDat(strand[hlfFrg],strand[5-hlfFrg],2,indexdata,1);
			if (retval != 0) return retval;
		
			int deltaLst = getDelta(strand[depRxn],strLngth[depRxn],7);		// two end hlf rearr
			retval = clssfyRA(indexdata,2,deltaLst,
								strand[depRxn],strLngth[depRxn],4);
			if (retval != 0) return retval;
		
			deltaLst = getDelta(strand[hlfCon],strLngth[hlfCon],7);			// contruction
			retval = clssfyHR(indexdata,3,deltaLst,strand[hlfCon],strLngth[hlfCon],true,6);
			if (retval != 0) return retval;						// classification error
		
			deltaLst = getDelta(strand[hlfFrg],strLngth[hlfFrg],7);			// fragmentation
			retval = clssfyHR(indexdata,4,-deltaLst,strand[hlfFrg],strLngth[hlfFrg],false,7);
			if (retval != 0) return retval;						// classification error
		
			mResult.mMainClass[mUnitRxn] = HD_REAR;

			if (mIndexToCreate != cIndexNone) {
				DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

				putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
				putindexbits(' ',BONDBITS,indexdata[0]);		// write rxn index entry
				putindexbits(' ',BONDBITS,indexdata[1]);
				putindexbits(' ',REARBITS,indexdata[2]);
				putindexbits(' ',HRXNBITS,indexdata[3]);
				putindexbits('l',HRXNBITS,indexdata[4]);
			
	//			if (mBatchClassifierP)
	//				mBatchClassifierP->incIndexEntries( HD_REAR );
				}
		
			mUnitRxn++;
			for (int i=1; i<=strLngth[depRxn]; i++)
				mAtomFlags[strand[depRxn][0]][strand[depRxn][i]] &= ~cAtomNotClassifiedYet;
			}
		else {
			for (int strnd=0; strnd<2; strnd++) {
				mResult.mRearStrandLen[strnd] = strLngth[strnd];
		
				int retval = bndSpcDat( strand[strnd*2],strand[strnd*2+1],
									strnd*2,indexdata,strnd);
				if (retval != 0) return retval;
		
				int deltaLst = getDelta(strand[strnd],strLngth[strnd],7);
				retval = clssfyRA(indexdata,2+strnd,deltaLst,
									strand[strnd],strLngth[strnd],strnd*2+4);
				if (retval != 0) return retval;
		
				for (int i=1; i<=strLngth[strnd]; i++)
					mAtomFlags[strand[strnd][0]][strand[strnd][i]] &= ~cAtomNotClassifiedYet;
				}

			mResult.mMainClass[mUnitRxn] = D_REAR;

			if (mIndexToCreate != cIndexNone) {
				DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

				putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
				putindexbits(' ',BONDBITS,indexdata[0]);		// write rxn index entry
				putindexbits(' ',BONDBITS,indexdata[1]);
				putindexbits(' ',REARBITS,indexdata[2]);
				putindexbits('l',REARBITS,indexdata[3]);
			
	//			if (mBatchClassifierP)
	//				mBatchClassifierP->incIndexEntries( D_REAR );
			
				if (mIndexToCreate == cIndexFullPermutation) {
					putindexbits('c',32,dbrxn.getReactionRegNo());
					putindexbits(' ',BONDBITS,indexdata[0]);
					putindexbits(' ',BONDBITS,indexdata[1]);
					putindexbits(' ',REARBITS,indexdata[3]);
					putindexbits('l',REARBITS,indexdata[2]);
					}
				}

			mUnitRxn++;
			}
		
		return cErrorNoError;
		}


	/**
	 * classify simple half reaction
	 * @return
	 */
	private int clssfyHR(int[] HRindex, int hr, int deltazp,
						 int[] strand, int strLngth, boolean construction, int flshBase) {
		int[] fgroup = new int[4];
		int[] fAtm = new int[4];
		
		mResult.mFlashMol[mUnitRxn][flshBase] = strand[0];
		mResult.mFlashAtom[mUnitRxn][flshBase] = strand[1];
		
		for (int i=0; i<mClassificationData.getHEntries(); i++) {
			if (deltazp != mClassificationData.getHRxnDelta(i)) continue;

			String unitName = mResult.mUnitName[mUnitRxn];
			if (construction) {
				if (unitName == null)
					unitName = mClassificationData.getHConstructionName(i);
				else
					unitName += " " + mClassificationData.getHConstructionName(i);
				}
			else {
				if (unitName == null)
					unitName = mClassificationData.getHFragmentationName(i);
				else
					unitName += " " + mClassificationData.getHFragmentationName(i);
				}
		
			HRindex[hr] = i<<PROPBITS;

			for (int j=1; j<=strLngth; j++)
				mAtomFlags[strand[0]][strand[j]] &= ~cAtomNotClassifiedYet;
		
			int mol = strand[0];
			int atm = strand[1];
		
			HRindex[hr] |= properties(mol,atm);
		
			int proMol = mCorProd[mol][atm];
			int proAtm = mCorAtom[mol][atm];
											// functional group bearing C
			int fatm = strand[mClassificationData.getHRxnGroupDef(i) & 15];
			int fproMol = mCorProd[mol][fatm];
			int fproAtm = mCorAtom[mol][fatm];
		
			HRindex[hr] <<= 8;

			int flashMol = -1;
			int maxgroup = -1;
			int maxatm = -1;
			int nrofGrps = 0;
			int[] dummy = new int[1];
											// if (contruction ^ fgroup is incoming)
			if ((mSigma[mol][atm] < mSigma[proMol][proAtm])
			^ ((mClassificationData.getHRxnGroupDef(i) & 16) == 16))
				{										// get type of leaving atom
				nrofGrps = leaving(mol,fatm,fproMol,fproAtm,fgroup,fAtm,dummy);
				flashMol = mol;
				maxatm = fatm;				// default setting in case nrofGrps = 0
				maxgroup = 2;

				for (int j=0; j<nrofGrps; j++) {  // mark leaving groups as already handled
					mAtomFlags[mol][fAtm[j]] &= ~cAtomNotClassifiedYet;
					int mapNo = mRxn.getMolecule(mol).getAtomMapNo(fAtm[j]);
					if (mapNo != 0) {
						boolean found = false;
						for (int rmol=mRxn.getReactants(); !found && rmol<mRxn.getMolecules(); rmol++) {
							StereoMolecule product = mRxn.getMolecule(rmol);
							for (int ratm=0; ratm<product.getAtoms(); ratm++) {
								if (product.getAtomMapNo(ratm) == mapNo) {
									mAtomFlags[rmol][ratm] &= ~cAtomNotClassifiedYet;
									found = true;
									break;
									}
								}
							}
						}
					}
				}
			else {										// get type of incoming atom
				nrofGrps = leaving(fproMol,fproAtm,mol,fatm,fgroup,fAtm,dummy);
				flashMol = fproMol;
				maxatm = fproAtm;			// default setting in case nrofGrps = 0
				maxgroup = 2;
		
				for (int j=0; j<nrofGrps; j++) { // mark incoming groups as already handled
					mAtomFlags[fproMol][fAtm[j]] &= ~cAtomNotClassifiedYet;
					int mapNo = mRxn.getMolecule(fproMol).getAtomMapNo(fAtm[j]);
					if (mapNo != 0) {
						boolean found = false;
						for (int rmol=0; !found && rmol<mRxn.getReactants(); rmol++) {
							StereoMolecule educt = mRxn.getMolecule(rmol);
							for (int ratm=0; ratm<educt.getAtoms(); ratm++) {
								if (educt.getAtomMapNo(ratm) == mapNo) {
									mAtomFlags[rmol][ratm] &= ~cAtomNotClassifiedYet;
									found = true;
									break;
									}
								}
							}
						}
					}
				}
		
			if (nrofGrps != 0) {
				maxgroup = fgroup[0];
				maxatm = fAtm[0];
				for (int j=1; j<nrofGrps; j++)
					{
					if (fgroup[j] > maxgroup)
						{
						maxgroup = fgroup[j];
						maxatm = fAtm[j];
						}
					}
				}
			HRindex[hr] |= maxgroup;
			mResult.mChngGrps[mUnitRxn][flshBase-6] = maxgroup;
			mResult.mFlashMol[mUnitRxn][flshBase+2] = flashMol;
			mResult.mFlashAtom[mUnitRxn][flshBase+2] = maxatm;
			return cErrorNoError;
			}

//		#if WRITE_DELTA
		if (gErrout != null) {
			try {
				gErrout.write(String.format("HR-deltaZP: %x",deltazp));
				gErrout.newLine();
				}
			catch (IOException ioe) {}
			}
		
		return cErrorHRClassifyError;
		}


	/**
	 * classify ring closure reaction
	 * @return
	 */
	private int clssfyCR(int[] CRindex, int cr, int deltazp, int[] strand, int strLen, boolean ringClosure) {
		for (int i=0; i<mClassificationData.getCEntries(); i++) {
			if (deltazp == mClassificationData.getCRxnDelta(i)) {
				if (ringClosure)
					mResult.mUnitName[mUnitRxn] = mClassificationData.getCRingClosureName(i);
				else
					mResult.mUnitName[mUnitRxn] = mClassificationData.getCRingOpeningName(i);
		
				CRindex[cr] = i;
		
				for (int j=1; j<=strLen; j++)
					mAtomFlags[strand[0]][strand[j]] &= ~cAtomNotClassifiedYet;
		
				return cErrorNoError;
				}
			}
		
//		#if WRITE_DELTA
		if (gErrout != null) {
			try {
				gErrout.write("CR-deltaZP: "+deltazp);
				gErrout.newLine();
				}
			catch (IOException ioe) {}
			}
		
		return cErrorCRClassifyError;
		}


	/**
	 * number of bits of properties is kept in PROPBITS in file 'classes.h'
	 */
	private int properties(int mol, int atm) {
		int prop =  mSigma[mol][atm] << 8;						// sigma
		prop |= mZ[mol][atm] << 5;								//   z
		prop |= mPi[mol][atm] << 3;								//   pi
		if (mRxn.getMolecule(mol).isAromaticAtom(atm))   prop |= 4L;
		if (mRxn.getMolecule(mol).isStabilizedAtom(atm)) prop |= 2L;
		if (mRxn.getMolecule(mol).isAllylicAtom(atm))	 prop |= 1L;
		return prop;
		}


	/**
	 * classify two end half reaction
	 * @return
	 */
	private int clssfyDR(int[] DRindex, int dr, int deltazp,
						 boolean construction, int[] strand, int strLngth, int flshBase) {
		
		mResult.mFlashMol[mUnitRxn][flshBase] = strand[0];
		mResult.mFlashAtom[mUnitRxn][flshBase] = strand[1];
		mResult.mFlashMol[mUnitRxn][flshBase+1] = strand[0];
		mResult.mFlashAtom[mUnitRxn][flshBase+1] = strand[strLngth];
		
		for (int i=0; i<mClassificationData.getDEntries(); i++) {
			if (deltazp == mClassificationData.getDRxnDelta(i)) {
				if (construction) {
					if (mResult.mUnitName[mUnitRxn] == null)
						mResult.mUnitName[mUnitRxn] = mClassificationData.getDConstructionName(i);
					else
						mResult.mUnitName[mUnitRxn] += " " + mClassificationData.getDConstructionName(i);
					}
				else
					{
					if (mResult.mUnitName[mUnitRxn] == null)
						mResult.mUnitName[mUnitRxn] = mClassificationData.getDFragmentationName(i);
					else
						mResult.mUnitName[mUnitRxn] += " " + mClassificationData.getDFragmentationName(i);
					}
		
				int m = strand[0];
				int atm1 = strand[1];
				int atm2 = strand[strLngth];
		
				DRindex[dr] = i;
		
				int sigma1 = mSigma[m][atm1];
				int sigma2 = mSigma[m][atm2];
				if (!construction) {
					sigma1--;
					sigma2--;
					}
				DRindex[dr] <<= 2;
				DRindex[dr] |= sigma1;
				DRindex[dr] <<= 2;
				DRindex[dr] |= sigma2;

				StereoMolecule mol = mRxn.getMolecule(m);

				int donating = 0;
				int wdrawing = 0;
				for (int j=1; j<=strLngth; j++) {
					for (int k=0; k<mol.getConnAtoms(strand[j]); k++) {
						int atm = mol.getConnAtom(strand[j], k);
						if ((mAtomType[m][atm] & 3) == 2) donating = 2;
						for (int l=0; l<mol.getConnAtoms(atm); l++) {
							if ((mAtomType[m][mol.getConnAtom(atm, l)] & 3) == 2)
								if (mol.getConnBondOrder(atm, l) > 1) wdrawing = 1;
							}
						}
					}
				DRindex[dr] <<= 2;
				DRindex[dr] |= donating + wdrawing;
		
				DRindex[dr+1] = (0xFFFFFFC3 & DRindex[dr]) | (16*sigma2 + 4*sigma1);
		
				return cErrorNoError;
				}
			}
		
//		#if WRITE_DELTA
		if (gErrout != null) {
			try {
				gErrout.write(String.format("DR-deltaZP: %x",deltazp));
				gErrout.newLine();
				}
			catch (IOException ioe) {}
			}
		
		return cErrorDRClassifyError;
		}


	/**
	 * classify two end half rearrangement
	 * @return
	 */
	private int clssfyRA(int[] RAindex, int ra, int deltazp, int[] strand, int strLngth, int flshBase) {
		int i,mol,formAtm,clvgAtm;

		mResult.mFlashMol[mUnitRxn][flshBase] = strand[0];
		mResult.mFlashAtom[mUnitRxn][flshBase] = strand[1];
		mResult.mFlashMol[mUnitRxn][flshBase+1] = strand[0];
		mResult.mFlashAtom[mUnitRxn][flshBase+1] = strand[strLngth];
		
		for (i=0; i<mClassificationData.getEEntries(); i++) {
			if (deltazp == mClassificationData.getERxnDelta(i)) {
				if (mResult.mUnitName[mUnitRxn] == null)
					mResult.mUnitName[mUnitRxn] = mClassificationData.getERearrangementName(i);
				else
					mResult.mUnitName[mUnitRxn] += " " + mClassificationData.getERearrangementName(i);
		
				mol = strand[0];
				formAtm = strand[1];
				clvgAtm = strand[strLngth];
		
				RAindex[ra] = i;
				RAindex[ra] <<= 2;
				RAindex[ra] |= mSigma[mol][formAtm];
				RAindex[ra] <<= 3;
				RAindex[ra] |= mSigma[mol][clvgAtm];
				RAindex[ra] <<= 2;
				RAindex[ra] |= mZ[mol][formAtm];
				RAindex[ra] <<= 2;
				RAindex[ra] |= mZ[mol][clvgAtm];
		
				return cErrorNoError;
				}
			}
		
//		#if WRITE_DELTA
		if (gErrout != null) {
			try {
				gErrout.write(String.format("RA-deltaZP: %x",deltazp));
				gErrout.newLine();
				}
			catch (IOException ioe) {}
			}
		
		return cErrorRAClassifyError;
		}


	/**
	 * classify refunctionalisation reaction
	 * @param deltazp
	 * @return
	 */
	private int clssfyRE(int deltazp) {
		int poin = 127;
		int step = 64;

		int entries = mClassificationData.getREntries();
		for (int i=0; i<8; i++) {
			if ((poin >= entries) || (deltazp > mClassificationData.getRRxnDelta(poin))) {
				poin -= step;
				step >>= 1;
				}
			else if (deltazp < mClassificationData.getRRxnDelta(poin)) {
				poin += step;
				step >>= 1;
				}
			else {
				mResult.mUnitName[mUnitRxn] = mClassificationData.getRRefuncName( poin );
				mResult.mClassResult[mUnitRxn] = poin;
				return cErrorNoError;
				}
			}
		
//		#if WRITE_DELTA
		if (gErrout != null) {
			try {
				gErrout.write(String.format("RE-deltaZP: %x",deltazp));
				gErrout.newLine();
				}
			catch (IOException ioe) {}
			}
		
		return cErrorREClassifyError;
		}


	/**
	 * store C-C bond specific data
	 * @return
	 */
	private int bndSpcDat(int[] strand1, int[] strand2, int flshBase, int[] bndIndex, int bi) {
		bndIndex[bi] = 0;
		
		int m1 = strand1[0];
		int m2 = strand2[0];
		int atm1 = strand1[1];
		int atm2 = strand2[1];
		
		mResult.mFlashAtom[mUnitRxn][flshBase] = m1;
		mResult.mFlashAtom[mUnitRxn][flshBase] = atm1;
		mResult.mFlashAtom[mUnitRxn][flshBase+1] = m2;
		mResult.mFlashAtom[mUnitRxn][flshBase+1] = atm2;
		
		if (m1 != m2) return cErrorNoError;
								// intermolecular because of different educt molecules
		
		if (conLngth(m1,atm1,atm2,false) != 0) bndIndex[bi] = 16;		// bit 4: intramolecular

		StereoMolecule mol1 = mRxn.getMolecule(m1);
		for (int i=0; i<mol1.getConnAtoms(atm1); i++) {
			int bnd = mol1.getConnBond(atm1, i);
			if (mol1.getBondAtom(0, bnd) == atm2
			 || mol1.getBondAtom(1, bnd) == atm2)
				{ bndIndex[bi] = 16; break; }
			}

		if (bndIndex[bi] == 0) return cErrorNoError;
						// intermolecular since there is no connection between atoms

		int connLength = conLngth(m1,atm1,atm2,true);
		if (connLength > 0) if (connLength > 15) connLength = 15;
		bndIndex[bi] |= connLength;							// bit 0-3: new ring's size
		return cErrorNoError;
		}


	/**
	 * classify refunctionalisation reactions
	 * @return
	 */
	private int refuncs() {
		int[][] leavType = new int[8][MAXCONNS];	// all leaving groups within one refunc
		int[][] incoType = new int[8][MAXCONNS];	// all incoming groups within one refunc
		int[] nrofLeav = new int[8];
		int[] nrofInco = new int[8];				// appropriate numbers of inc/leav grps
		int[][] leavAtm = new int[8][MAXCONNS];		// corresponding atom numbers			
		int[][] incoAtm = new int[8][MAXCONNS];
		int[][] stereo = new int[8][1];
		int[][] strand = new int[2][8];
		int[][][] fGrpCombs = new int[MAXUNITS][MAXINDICES][4];
		int[][] REindex = new int[MAXUNITS][2];
				// contains index information for all single refunc rxns within one
				// REACCS entry (in case mIndexToCreate == cIndexOnePermToFile) or
				// all possible indices of one refunc rxn due to different combinations
				// of functional groups (mIndexToCreate == cIndexFullPermutation)
		int[] delLst = new int[2];
		int[] dummy = new int[1];

		for (int mol=0; mol<mRxn.getReactants(); mol++) {			// set changMsk of carbon atoms
			for (int atm=0; atm<mRxn.getMolecule(mol).getAtoms(); atm++) {
				if ((mAtomFlags[mol][atm] & cAtomChanges) != 0)
					continue;										// for all mapped and
				if (mRxn.getMolecule(mol).getAtomMapNo(atm) == 0) continue;	// until now not as
				int proMol = mCorProd[mol][atm];						// changing marked atoms
				int proAtm = mCorAtom[mol][atm];
				if (mZ[mol][atm] != mZ[proMol][proAtm]) continue;
						// this change ought to be detected by now. If changMsk is 0,
						// it was reset due to mesomeric bond changes (e.g. pyridine)
				if (leaving(mol,atm,proMol,proAtm,leavType[0],leavAtm[0],dummy) != 0) {
					mAtomFlags[mol][atm] |= (cAtomChanges + cAtomNotClassifiedYet);
					mAtomFlags[proMol][proAtm] |= (cAtomChanges + cAtomNotClassifiedYet);
					continue;
					}
				if (leaving(proMol,proAtm,mol,atm,incoType[0],incoAtm[0],dummy) != 0) {
					mAtomFlags[mol][atm] |= (cAtomChanges + cAtomNotClassifiedYet);
					mAtomFlags[proMol][proAtm] |= (cAtomChanges + cAtomNotClassifiedYet);
					continue;
					}
				}
			}
		
		int cChange = 0;
		int firstRefuUnit = mUnitRxn;
		for (int mol=0; mol<mRxn.getReactants(); mol++) {
			strand[0][0] = strand[1][0] = mol;
			for (int atm=0; atm<mRxn.getMolecule(mol).getAtoms(); atm++) {
				if (mAtomType[mol][atm] != CARBON) continue;
				if ((mAtomFlags[mol][atm] & cAtomNotClassifiedYet) == 0) continue;
				strand[0][1] = atm;
				int strLngth = findStrand(strand[0]);
				if (strLngth == 0) continue;				// if not endcarbon of strand
				if (strLngth == -1) return cErrorForkedOrLongStrand;	// strand exceeds 7 carbons
				delLst[0] = getDelta(strand[0],strLngth,strLngth);
		
				for (int i=1; i<=strLngth; i++)
					strand[1][i] = strand[0][strLngth-i+1];
				delLst[1] = getDelta(strand[1],strLngth,strLngth);
		
				int higher = (delLst[1] > delLst[0]) ? 1 : 0;

				int retval = clssfyRE(delLst[higher]);
				if (retval != 0) return retval;
		
				for (int i=1; i<=strLngth; i++)
					mAtomFlags[mol][strand[0][i]] &= ~cAtomNotClassifiedYet;
		
				cChange++;
		
				if (strLngth == 1) mResult.mMainClass[mUnitRxn] = C1_REFU;
				else			   mResult.mMainClass[mUnitRxn] = LONG_REFU;
		
		        for (int i=1; i<=strLngth; i++)  // initialize arrays for changing groups
					for (int j=0; j<4; j++) {
			            incoType[i][j] = 0;
		    	        leavType[i][j] = 0;
			            }
		
				for (int i=1; i<=strLngth; i++) {   // get class numbers of changing groups
					int actlAtm = strand[higher][i];
					int proMol = mCorProd[mol][actlAtm];
					int proAtm = mCorAtom[mol][actlAtm];
																	// incoming atoms
					nrofInco[i] = leaving(proMol,proAtm,mol,actlAtm,incoType[i],incoAtm[i],stereo[i]);
		
					for (int j=0; j<nrofInco[i]; j++) {		// reset change flag of incoming atoms
						mAtomFlags[proMol][incoAtm[i][j]] &= ~cAtomNotClassifiedYet;
						int mapNo = mRxn.getMolecule(proMol).getAtomMapNo(incoAtm[i][j]);
						if (mapNo != 0) {
							boolean found = false;
							for (int rmol=0; !found && rmol<mRxn.getReactants(); rmol++) {
								StereoMolecule educt = mRxn.getMolecule(rmol);
								for (int ratm=0; ratm<educt.getAtoms(); ratm++) {
									if (educt.getAtomMapNo(ratm) == mapNo) {
										mAtomFlags[rmol][ratm] &= ~cAtomNotClassifiedYet;
										found = true;
										break;
										}
									}
								}
							}
						}

					int expected = 0;
					for (int j=0; j<mClassificationData.getRRxnFGroups(mResult.mClassResult[mUnitRxn] ); j++)
						if ((mClassificationData.getRRxnDef( mResult.mClassResult[mUnitRxn], 0, j ) & 64) != 0)
							if ((mClassificationData.getRRxnDef( mResult.mClassResult[mUnitRxn], 0, j ) & 15) == i) 
								expected++;
					int surplus = nrofInco[i]-expected;
					if (surplus < 0) {
						if (delLst[higher] == 0) surplus = 0;		// simple substitution
						else				return cErrorIncoOrLeavMissing;
						}
					while (surplus != 0) {
						boolean found = false;
						for (int j=0; !found && j<nrofInco[i]; j++) {
							for (int k=j+1; k<nrofInco[i]; k++) {
								if (incoType[i][j] == incoType[i][k]) {
									for (int l=k+1; l<nrofInco[i]; l++) {
										incoType[i][l-1] = incoType[i][l];
										incoAtm[i][l-1] = incoAtm[i][l];
										}
									nrofInco[i]--;
									surplus--;
									found = true;
									break;
									}
								}
							}
						if (!found) {
							surplus = 0;				// cut, if not enough double found
							nrofInco[i] = expected;
							}
						}
																	// leaving atoms
					nrofLeav[i] = leaving(mol,actlAtm,proMol,proAtm,leavType[i],leavAtm[i],stereo[i]);
		
					for (int j=0; j<nrofLeav[i]; j++) {		// reset change flag of leaving atoms
						mAtomFlags[mol][leavAtm[i][j]] &= ~cAtomNotClassifiedYet;
						int mapNo = mRxn.getMolecule(mol).getAtomMapNo(leavAtm[i][j]);
						if (mapNo != 0) {
							boolean found = false;
							for (int rmol=mRxn.getReactants(); !found && rmol<mRxn.getMolecules(); rmol++) {
								StereoMolecule product = mRxn.getMolecule(rmol);
								for (int ratm=0; ratm<product.getAtoms(); ratm++) {
									if (product.getAtomMapNo(ratm) == mapNo) {
										mAtomFlags[rmol][ratm] &= ~cAtomNotClassifiedYet;
										found = true;
										break;
										}
									}
								}
							}
						}
		
					expected = 0;
					for (int j=0; j<mClassificationData.getRRxnFGroups(mResult.mClassResult[mUnitRxn]); j++)
						if ((mClassificationData.getRRxnDef(mResult.mClassResult[mUnitRxn], 0, j ) & 64) == 0)
							if ((mClassificationData.getRRxnDef(mResult.mClassResult[mUnitRxn], 0, j ) & 15) == i) 
								expected++;
					surplus = nrofLeav[i]-expected;
					if (surplus < 0) {
						if (delLst[higher] == 0) surplus = 0;		// simple substitution
						else				return cErrorIncoOrLeavMissing;
						}
					while (surplus != 0) {
						boolean found = false;
						for (int j=0; !found && j<nrofLeav[i]; j++) {
							for (int k=j+1; k<nrofLeav[i]; k++) {
								if (leavType[i][j] == leavType[i][k]) {
									for (int l=k+1; l<nrofLeav[i]; l++) {
										leavType[i][l-1] = leavType[i][l];
										leavAtm[i][l-1] = leavAtm[i][l];
										}
									nrofLeav[i]--;
									surplus--;
									found = true;
									break;
									}
								}
							}
						if (!found) {
							surplus = 0;				// cut, if not enough double found
							nrofLeav[i] = expected;
							}
						}
					}

				if (gErrout != null) {
					try {
						gErrout.write("leaving/incoming info of one refu strand; UnitRxn:"+mUnitRxn);
						gErrout.newLine();
						for (int i=1; i<=strLngth; i++) {
							int actlAtm = strand[higher][i];
							gErrout.write("atm: "+actlAtm+"; leaving:");
							for (int j=0; j<nrofLeav[i]; j++)
								gErrout.write(" "+leavType[i][j]);
							gErrout.write(" incoming: ");
							for (int j=0; j<nrofInco[i]; j++)
								gErrout.write(" "+incoType[i][j]);
							gErrout.newLine();
							}
						}
					catch (IOException ioe) {}
					}
		
				int stereoInfo = 0;
				if (delLst[0] == 0) {		// determine stereo info in simple substitutions
					stereoInfo = stereo[1][0];
					mResult.mStereoInfo[mUnitRxn] = stereo[1][0];
					}
		
				if (mIndexToCreate == cIndexNone) { mUnitRxn++; continue; }
		
				if (strLngth == 1)					// one atom refunctionalisation
					{
					REindex[mUnitRxn][0] = mResult.mClassResult[mUnitRxn]<<PROPBITS;
					REindex[mUnitRxn][0] |= properties(mol,(int)strand[higher][1]);
					REindex[mUnitRxn][0] <<= 2;
					REindex[mUnitRxn][0] |= stereoInfo;			// 2 bit stereo info
					}
				else			// refunctionalisation strand contains more than 1 atom
					{
					REindex[mUnitRxn][0] = mResult.mClassResult[mUnitRxn]<<PROPBITS;
					int properties1 = properties(mol,(int)strand[higher][1]);
					REindex[mUnitRxn][0] |= properties1;
					REindex[mUnitRxn][0] <<= PROPBITS;
					int properties2 = properties(mol,(int)strand[1-higher][1]);
					REindex[mUnitRxn][0] |= properties2;
		
					if (mClassificationData.getRRxnSymmetric(mResult.mClassResult[mUnitRxn]))
						{
						REindex[mUnitRxn][1] = mResult.mClassResult[mUnitRxn]<<PROPBITS;
						REindex[mUnitRxn][1] |= properties2;
						REindex[mUnitRxn][1] <<= PROPBITS;
						REindex[mUnitRxn][1] |= properties1;
						}

					mResult.mFlashMol[mUnitRxn][1] = strand[1-higher][0];
					mResult.mFlashAtom[mUnitRxn][1] = strand[1-higher][1];
					}
		
				mResult.mFlashMol[mUnitRxn][0] = strand[higher][0];
				mResult.mFlashAtom[mUnitRxn][0] = strand[higher][1];
		
				for (int i=0; i<mClassificationData.getRRxnMasks( mResult.mClassResult[mUnitRxn] ); i++) {
					for (int j=0; j<4; j++)
						fGrpCombs[mUnitRxn][i][j] = 0;
					for (int j=0; j<mClassificationData.getRRxnFGroups( mResult.mClassResult[mUnitRxn] ); j++) {
						int strndPos = mClassificationData.getRRxnDef( mResult.mClassResult[mUnitRxn], i, j ) & 15;
						int num = (mClassificationData.getRRxnDef( mResult.mClassResult[mUnitRxn], i, j ) & 48) >> 4;
						if ((mClassificationData.getRRxnDef( mResult.mClassResult[mUnitRxn], i, j ) & 64) != 0)
							fGrpCombs[mUnitRxn][i][j] = incoType[strndPos][num];
						else
							fGrpCombs[mUnitRxn][i][j] = leavType[strndPos][num];
						}
					}
		
				boolean schonda = false;
				for (int j=firstRefuUnit; j<mUnitRxn; j++) {	// forgo this unit if equal unit already detected
					if (mClassificationData.getRRxnSymmetric(mResult.mClassResult[mUnitRxn])) {
						if (REindex[mUnitRxn][0] != REindex[j][0]
						&& REindex[mUnitRxn][0] != REindex[j][1]) continue;
		
						if (REindex[mUnitRxn][0] == REindex[j][0]) {
							for (int k=0; k<(mClassificationData.getRRxnMasks( mResult.mClassResult[mUnitRxn] ) >>1); k++) {
								if (fGrpCombs[mUnitRxn][k][0] == fGrpCombs[j][0][0]
								 && fGrpCombs[mUnitRxn][k][1] == fGrpCombs[j][0][1]
								 && fGrpCombs[mUnitRxn][k][2] == fGrpCombs[j][0][2]
								 && fGrpCombs[mUnitRxn][k][3] == fGrpCombs[j][0][3]) {
									schonda = true;
									break;
									}
								}
							}
						if (REindex[mUnitRxn][0] == REindex[j][1]) {
							for (int k=(mClassificationData.getRRxnMasks( mResult.mClassResult[mUnitRxn] ) >>1);
									k<mClassificationData.getRRxnMasks( mResult.mClassResult[mUnitRxn] ); k++) {
								if (fGrpCombs[mUnitRxn][k][0] == fGrpCombs[j][0][0]
								 && fGrpCombs[mUnitRxn][k][1] == fGrpCombs[j][0][1]
								 && fGrpCombs[mUnitRxn][k][2] == fGrpCombs[j][0][2]
								 && fGrpCombs[mUnitRxn][k][3] == fGrpCombs[j][0][3]) {
									schonda = true;
									break;
									}
								}
							}
						}
					else {
						if (REindex[mUnitRxn][0] != REindex[j][0]) continue;
						for (int k=0; k<mClassificationData.getRRxnMasks( mResult.mClassResult[mUnitRxn] ); k++) {
							if (fGrpCombs[mUnitRxn][k][0] == fGrpCombs[j][0][0]
							 && fGrpCombs[mUnitRxn][k][1] == fGrpCombs[j][0][1]
							 && fGrpCombs[mUnitRxn][k][2] == fGrpCombs[j][0][2]
							 && fGrpCombs[mUnitRxn][k][3] == fGrpCombs[j][0][3]) {
								schonda = true;
								break;
								}
							}
						}
					if (schonda) break;
					}
		
				if (schonda) {
					mResult.mUnitName[mUnitRxn] = null;
					continue;
					}
		
				for (int i=0; i<4; i++)
					mResult.mChngGrps[mUnitRxn][i] = fGrpCombs[mUnitRxn][0][i];
		
				for (int j=0; j<mClassificationData.getRRxnFGroups( mResult.mClassResult[mUnitRxn] ); j++) {
					int strndPos = mClassificationData.getRRxnDef( mResult.mClassResult[mUnitRxn], 0, j ) & 15;
					int num = (mClassificationData.getRRxnDef( mResult.mClassResult[mUnitRxn], 0, j ) & 48) >> 4;
					if ((mClassificationData.getRRxnDef( mResult.mClassResult[mUnitRxn], 0, j ) & 64) != 0) {
						mResult.mFlashMol[mUnitRxn][j+2] = mCorProd[mol][strand[higher][1]];
						mResult.mFlashAtom[mUnitRxn][j+2] = incoAtm[strndPos][num];
						}
					else {
						mResult.mFlashMol[mUnitRxn][j+2] = mol;
						mResult.mFlashAtom[mUnitRxn][j+2] = leavAtm[strndPos][num];
						}
					}
		
				mUnitRxn++;
				if (mUnitRxn == MAXUNITS) return cErrorUnitReactionLimit;
				}
			if (mUnitRxn == MAXUNITS) return cErrorUnitReactionLimit;
			}
		
		if (cChange != 0) {
			if (mIndexToCreate != cIndexNone) {
				DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

				if (mIndexToCreate == cIndexOnePermToFile) {	// write one entry per separated unit rxn to file
					int unitRxns = mUnitRxn;
					for (mUnitRxn=firstRefuUnit; mUnitRxn<unitRxns; mUnitRxn++) {
						int midBits = (mResult.mMainClass[mUnitRxn] == C1_REFU) ? 8+PROPBITS+2 : 8+2*PROPBITS;
						putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
						putindexbits(' ',midBits,REindex[mUnitRxn][0]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][0]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][1]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][2]);
						putindexbits('l',8,fGrpCombs[mUnitRxn][0][3]);
			
	//					mBatchClassifierP->incIndexEntries( mResult.mMainClass[mUnitRxn] );
						}
					}
				else if (mIndexToCreate == cIndexFullPermutation) {		// store all combinations of all unit rxns into array
					int unitRxns = mUnitRxn;
					for (mUnitRxn=firstRefuUnit; mUnitRxn<unitRxns; mUnitRxn++) {
						int midBits = (mResult.mMainClass[mUnitRxn] == C1_REFU) ? 8+PROPBITS+2 : 8+2*PROPBITS;
						putindexbits('n',32,dbrxn.getReactionRegNo());
						putindexbits(' ',midBits,REindex[mUnitRxn][0]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][0]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][1]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][2]);
						putindexbits('l',8,fGrpCombs[mUnitRxn][0][3]);
	//????				mResult.mClassResult[mUnitRxn] = mResult.mClassResult[mUnitRxn];
						for (int j=1; j<mClassificationData.getRRxnMasks( mResult.mClassResult[mUnitRxn] ); j++) {
							if ((mClassificationData.getRRxnSymmetric(mResult.mClassResult[mUnitRxn]))
							&& (j >= (mClassificationData.getRRxnMasks( mResult.mClassResult[mUnitRxn] )/2))) {
								putindexbits('c',32,dbrxn.getReactionRegNo());
								putindexbits(' ',midBits,REindex[mUnitRxn][1]);
								putindexbits(' ',8,fGrpCombs[mUnitRxn][j][0]);
								putindexbits(' ',8,fGrpCombs[mUnitRxn][j][1]);
								putindexbits(' ',8,fGrpCombs[mUnitRxn][j][2]);
								putindexbits('l',8,fGrpCombs[mUnitRxn][j][3]);
								}
							else {
								putindexbits('c',32,dbrxn.getReactionRegNo());
								putindexbits(' ',midBits,REindex[mUnitRxn][0]);
								putindexbits(' ',8,fGrpCombs[mUnitRxn][j][0]);
								putindexbits(' ',8,fGrpCombs[mUnitRxn][j][1]);
								putindexbits(' ',8,fGrpCombs[mUnitRxn][j][2]);
								putindexbits('l',8,fGrpCombs[mUnitRxn][j][3]);
								}
							}
						}
					}
				}
			}
		
		int hChange = 0;						// refunctionalisation at non-carbon atom
		firstRefuUnit = mUnitRxn;
		for (int m=0; m<mRxn.getReactants(); m++) {
			StereoMolecule mol = mRxn.getMolecule(m);
			for (int atm=0; atm<mol.getAtoms(); atm++) {
				if (mAtomType[m][atm] == CARBON) continue;
				if ((mAtomFlags[m][atm] & cAtomNotClassifiedYet) == 0) continue;
				if (mol.getAtomMapNo(atm) == 0) continue;
				int proMol = mCorProd[m][atm];
				int proAtm = mCorAtom[m][atm];
				int deltaZ = mZ[m][atm]-mZ[proMol][proAtm];
				if (deltaZ < 0)
					mResult.mUnitName[mUnitRxn] = "X"+(-deltaZ)+Molecule.cAtomLabel[mol.getAtomicNo(atm)];
				else if (deltaZ > 0)
					mResult.mUnitName[mUnitRxn] = "R"+deltaZ+Molecule.cAtomLabel[mol.getAtomicNo(atm)];
				else
					mResult.mUnitName[mUnitRxn] = "S,"+Molecule.cAtomLabel[mol.getAtomicNo(atm)];
		
				mAtomFlags[m][atm] &= ~cAtomNotClassifiedYet;
		
				mResult.mMainClass[mUnitRxn] = HETERO_REFU;
		
				hChange++;
		
				if (mIndexToCreate == cIndexNone) { mUnitRxn++; continue; }
		
				incoType[0][0] = incoType[0][1] = 0;
				nrofInco[0] = leaving(mCorProd[m][atm],mCorAtom[m][atm],m,atm,incoType[0],incoAtm[0],dummy);
		
				while (nrofInco[0] > 2) {
					boolean found = false;
					for (int j=0; !found && j<nrofInco[0]; j++) {
						for (int k=j+1; k<nrofInco[0]; k++) {
							if (incoType[0][j] == incoType[0][k]) {
								for (int l=k+1; l<nrofInco[0]; l++) {
									incoType[0][l-1] = incoType[0][l];
									incoAtm[0][l-1] = incoAtm[0][l];
									nrofInco[0]--;
									}
								found = true;
								break;
								}
							}
						}
					if (!found)
						nrofInco[0] = 2;
					}
		
				leavType[0][0] = leavType[0][1] = 0;
				nrofLeav[0] = leaving(m,atm,mCorProd[m][atm],mCorAtom[m][atm],leavType[0],leavAtm[0],dummy);
		
				while (nrofLeav[0] > 2) {
					boolean found = false;
					for (int j=0; !found && j<nrofLeav[0]; j++) {
						for (int k=j+1; k<nrofLeav[0]; k++) {
							if (leavType[0][j] == leavType[0][k]) {
								for (int l=k+1; l<nrofLeav[0]; l++) {
									leavType[0][l-1] = leavType[0][l];
									leavAtm[0][l-1] = leavAtm[0][l];
									nrofLeav[0]--;
									}
								found = true;
								break;
								}
							}
						}
					if (!found)
						nrofLeav[0] = 2;
					}
		
				REindex[mUnitRxn][0] = (deltaZ+8) << 8;
				REindex[mUnitRxn][0] |= mol.getAtomicNo(atm);
				REindex[mUnitRxn][0] <<= PROPBITS;
				REindex[mUnitRxn][0] |= properties(m,atm);
				fGrpCombs[mUnitRxn][0][0] = incoType[0][0];
				fGrpCombs[mUnitRxn][0][1] = incoType[0][1];
				fGrpCombs[mUnitRxn][0][2] = leavType[0][0];
				fGrpCombs[mUnitRxn][0][3] = leavType[0][1];
				fGrpCombs[mUnitRxn][1][0] = incoType[0][0];
				fGrpCombs[mUnitRxn][1][1] = incoType[0][1];
				fGrpCombs[mUnitRxn][1][2] = leavType[0][1];
				fGrpCombs[mUnitRxn][1][3] = leavType[0][0];
				fGrpCombs[mUnitRxn][2][0] = incoType[0][1];
				fGrpCombs[mUnitRxn][2][1] = incoType[0][0];
				fGrpCombs[mUnitRxn][2][2] = leavType[0][0];
				fGrpCombs[mUnitRxn][2][3] = leavType[0][1];
				fGrpCombs[mUnitRxn][3][0] = incoType[0][1];
				fGrpCombs[mUnitRxn][3][1] = incoType[0][0];
				fGrpCombs[mUnitRxn][3][2] = leavType[0][1];
				fGrpCombs[mUnitRxn][3][3] = leavType[0][0];
		
				mResult.mChngGrps[mUnitRxn][0] = incoType[0][0];
				mResult.mChngGrps[mUnitRxn][1] = incoType[0][1];
				mResult.mChngGrps[mUnitRxn][2] = leavType[0][0];
				mResult.mChngGrps[mUnitRxn][3] = leavType[0][1];

				mResult.mFlashMol[mUnitRxn][1] = m;
				mResult.mFlashAtom[mUnitRxn][1] = atm;
				mResult.mFlashMol[mUnitRxn][2] = mCorProd[m][atm];
				mResult.mFlashAtom[mUnitRxn][2] = incoAtm[0][0];
				mResult.mFlashMol[mUnitRxn][3] = mCorProd[m][atm];
				mResult.mFlashAtom[mUnitRxn][3] = incoAtm[0][1];
				mResult.mFlashMol[mUnitRxn][4] = m;
				mResult.mFlashAtom[mUnitRxn][4] = leavAtm[0][0];
				mResult.mFlashMol[mUnitRxn][5] = m;
				mResult.mFlashAtom[mUnitRxn][5] = leavAtm[0][1];

				boolean schonda = false;
				for (int i=firstRefuUnit; i<mUnitRxn; i++) {
					if (REindex[i][0] != REindex[mUnitRxn][0]) continue;
					for (int j=0; j<4; j++) {
						if (fGrpCombs[mUnitRxn][0][0] == fGrpCombs[i][j][0]
						 && fGrpCombs[mUnitRxn][0][1] == fGrpCombs[i][j][1]
						 && fGrpCombs[mUnitRxn][0][2] == fGrpCombs[i][j][2]
						 && fGrpCombs[mUnitRxn][0][3] == fGrpCombs[i][j][3]) {
							schonda = true;
							break;
							}
						}
					}
				if (schonda) {
					mResult.mUnitName[mUnitRxn] = null;
					continue;
					}
		
				mUnitRxn++;
				if (mUnitRxn == MAXUNITS) return cErrorUnitReactionLimit;
				}
			if (mUnitRxn == MAXUNITS) return cErrorUnitReactionLimit;
			}
		
		if (hChange != 0) {
			if (mIndexToCreate != cIndexNone) {
				DatabaseReaction dbrxn = (DatabaseReaction)mRxn;

				if (mIndexToCreate == cIndexOnePermToFile) {	// write one entry per separated refunc rxn to file
					int unitRxns = mUnitRxn;
					for (mUnitRxn=firstRefuUnit; mUnitRxn<unitRxns; mUnitRxn++) {
						putindexbits('n',32,(dbrxn.getReactionYield() << 24) + dbrxn.getReactionRegNo());
						putindexbits(' ',12+PROPBITS,REindex[mUnitRxn][0]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][0]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][1]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][2]);
						putindexbits('l',8,fGrpCombs[mUnitRxn][0][3]);
			
	//					mBatchClassifierP->incIndexEntries( HETERO_REFU );
						}
					}
				else if (mIndexToCreate == cIndexFullPermutation) {		// store all combinations of one(!) refunc rxn in array
					int unitRxns = mUnitRxn;
					for (mUnitRxn=firstRefuUnit; mUnitRxn<unitRxns; mUnitRxn++)
						{
						putindexbits('n',32,dbrxn.getReactionRegNo());
						putindexbits(' ',12+PROPBITS,REindex[mUnitRxn][0]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][0]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][1]);
						putindexbits(' ',8,fGrpCombs[mUnitRxn][0][2]);
						putindexbits('l',8,fGrpCombs[mUnitRxn][0][3]);
						for (int j=1; j<4; j++) {
							putindexbits('c',32,dbrxn.getReactionRegNo());
							putindexbits(' ',12+PROPBITS,REindex[mUnitRxn][0]);
							putindexbits(' ',8,fGrpCombs[mUnitRxn][j][0]);
							putindexbits(' ',8,fGrpCombs[mUnitRxn][j][1]);
							putindexbits(' ',8,fGrpCombs[mUnitRxn][j][2]);
							putindexbits('l',8,fGrpCombs[mUnitRxn][j][3]);
							}
						}
					}
				}
			}
		
		//if (cChange == 0 && hChange == 0) return kErrNoChangingAtoms;
		// in OOP Cognos it is not considered an error if no changing groups are found
		return cErrorNoError;
		}


	private void findClassNam() {
		if (mUnitRxn == 0) {
			mResult.mClassName = "no changing atoms";
			return;
			}
		
		int commonClass = mResult.mMainClass[0];
		
		for (int i=1; i<mUnitRxn; i++) {
			if (commonClass != mResult.mMainClass[i])
				{
				commonClass = 127;
				break;
				}
			}
		
		switch (commonClass) {
			case S_CONST:
				mResult.mClassName = "single construction";
				break;
			case E_RING_C:
				mResult.mClassName = "electrocyclic ring closure";
				break;
			case HD_D_CONST:
				mResult.mClassName = "asymmetric double construction";
				break;
			case D_D_CONST:
				mResult.mClassName = "symmetric double construction";
				break;
			case S_FRAG:
				mResult.mClassName = "fragmentation";
				break;
			case E_RING_O:
				mResult.mClassName = "electrocyclic ring opening";
				break;
			case HD_D_FRAG:
				mResult.mClassName = "asymmetric double fragmentation";
				break;
			case D_D_FRAG:
				mResult.mClassName = "symmetric double fragmentation";
				break;
			case HD_REAR:
				mResult.mClassName = "asymmetric rearrangement";
				break;
			case D_REAR:
				mResult.mClassName = "symmetric rearrangement";
				break;
			case C1_REFU:
				mResult.mClassName = "refunctionalisation";
				break;
			case LONG_REFU:
				mResult.mClassName = "refunctionalisation";
				break;
			case HETERO_REFU:
				mResult.mClassName = "refunctionalisation";
				break;
			default:
				mResult.mClassName = "miscellaneous classes";
				break;
			}
		}


	private void putindexbits(int mode, int bits, int datum) {
		if (mIndexToCreate == cIndexNone) return;
		
		if (mode == 'n') {	// start new entry
			indexnum = 0;
			indexpoin = 0;
			data = 0;
			availbits = 32;
			}
		
		if (mode == 'c') {
			indexnum++;
			indexpoin = 0;
			data = 0;
			availbits = 32;
			}
		
		if ((mIndexToCreate == cIndexOnePermToFile) && (indexnum > 0)) return;
		
		if (bits != 0) mask = 1 << (bits-1);
		for (int i=0; i<bits; i++) {
			if (availbits == 0) {
				if (mIndexToCreate == cIndexFullPermutation)
					mResult.mIndex[mUnitRxn][indexnum][indexpoin++] = data;
//				if (mIndexToCreate == cIndexOnePermToFile)
//					mBatchClassifierP->longToIndexFile( &data,
//									 mClassResultP->getMainClass( mUnitRxn ) );
				data = 0;
				availbits = 32;
				}
			data <<= 1;
			if ((mask & datum) != 0) data |= 1;
			mask >>= 1;
			availbits--;
			}
		
		if (mode == 'l') {	// write last bits of entry
			data <<= availbits;
			if (mIndexToCreate == cIndexFullPermutation)
				mResult.mIndex[mUnitRxn][indexnum][indexpoin] = data;
//			if (mIndexToCreate == cIndexOnePermToFile)
//				mBatchClassifierP->longToIndexFile( &data,
//								 mClassResultP->getMainClass( mUnitRxn ) );
			}
		}


	/**
	 * if (checkProduct == 0) then return length of shortest chain between two
	 * 							atoms except direct bond (any mol allowed )
	 * if (checkProduct == 1) then return length of shortest chain between two
	 *							educt atoms, that exists in the product as well
	 *							(only educt atoms allowed)
	 * @param mol
	 * @param atm1
	 * @param atm2
	 * @param checkProduct
	 * @return
	 */
	private int conLngth( int mol, int atm1, int atm2, boolean checkProduct ) {
		int[][] mask = new int[2][MAXATOMS];
		int[] lastAtm = new int[MAXATOMS];
		int[] proAtm = new int[MAXATOMS];

		StereoMolecule mol1 = mRxn.getMolecule(mol);

		for (int atm=0; atm<mol1.getAtoms(); atm++) {
			mask[0][atm] = 0;
			mask[1][atm] = 0;
			lastAtm[atm] = 0;
			}
		
		mask[0][atm1] = 2;
		lastAtm[atm1] = atm2;

		int proM = -1;
		if (checkProduct) {
			proM = mCorProd[mol][atm1];
			proAtm[atm1] = mCorAtom[mol][atm1];
			}

		int msk = 0;
		for (int depth=2; depth<30; depth++) {
			for (int atm=0; atm<mol1.getAtoms(); atm++) {
				if (mask[msk][atm] != depth) continue;
				for (int conn=0; conn<mol1.getConnAtoms(atm); conn++) {
					int connNo = mol1.getConnAtom(atm, conn);
					if (connNo == atm1 || connNo == lastAtm[atm]) continue;
					if (checkProduct) {
						boolean proFound = false;
						StereoMolecule product = mRxn.getMolecule(proM);
						for (int proConn=0;	proConn<product.getConnAtoms(proAtm[atm]); proConn++) {
							int proConnNo = product.getConnAtom(proAtm[atm], proConn);
							if (mol1.getAtomMapNo(connNo) != 0) {
								if (mol1.getAtomMapNo(connNo) == product.getAtomMapNo(proConnNo)) {
									proAtm[connNo] = proConnNo;
									proFound = true;
									break;
									}
								}
							else {
								if (mol1.getAtomicNo(connNo) == product.getAtomicNo(proConnNo)
								&& product.getAtomMapNo(proConnNo ) == 0) {
									proAtm[connNo] = proConnNo;
									proFound = true;
									break;
									}
								}
							}
						if (!proFound) continue;
						}
					if (connNo == atm2) return depth;
					mask[1-msk][connNo] = depth+1;
					lastAtm[connNo] = atm;
					}
				}
			msk = 1-msk;
			}
		
		return 0;
		}


	/**
	 * Searchs for largest connected carbon fragment to carbon atm in molecule mol
	 * and returns list of fragment atom numbers. fragAtms[0] = no of frag atoms.
	 * @param m
	 * @param startAtm
	 * @return
	 */
	private int findFrag(int m, int startAtm, int[] fragAtms, int[] foundMsk ) {
		int[] histAtms = new int[MAXATOMS];
		int[] histConn = new int[MAXATOMS];
		int histpoin=1;
		int frgAtmNo=2;

		StereoMolecule mol = mRxn.getMolecule(m);
		histAtms[0] = MAXATOMS;
		histAtms[1] = startAtm;
		histConn[1] = 0;
		fragAtms[1] = startAtm;
		foundMsk[startAtm] = 1;
		
		if (gErrout != null) {
			try {
				gErrout.newLine();
				gErrout.write("findfrag->");
				}
			catch (IOException ioe) {}
			}

		while (true) {
			boolean overflow = false;
			int atm = -1;
			do {
				overflow = (histConn[histpoin] == mol.getConnAtoms(histAtms[histpoin]));
				if (overflow) break;
				atm = mol.getConnAtom(histAtms[histpoin], histConn[histpoin]++);
				} while ((atm == histAtms[histpoin-1])
						|| (foundMsk[atm] != 0)
						|| (mAtomType[m][atm] != CARBON));
		
			if (!overflow) {
				if (frgAtmNo > MAXFATMS)
					return cErrorFragmentAtomLimit;	// fragment carbons exceed limit
				fragAtms[frgAtmNo] = atm;
				foundMsk[atm] = frgAtmNo++;
				histAtms[++histpoin] = atm;
				histConn[histpoin] = 0;
		
				if (gErrout != null)
					try { gErrout.write(" "+atm); } catch (IOException ioe) {}
				}
			else {
				if (--histpoin == 0) break;
				}
			}
		fragAtms[0] = frgAtmNo-1;
		
		if (gErrout != null)
			try { gErrout.newLine(); } catch (IOException ioe) {}
		
		return cErrorNoError;
		}


/*/////////////////////////////////////////////////////////////////// private //
void flagAromaticBonds( int mol, int aromAtms[] )
////////////////////////////////////////////////////////////////////////////////
{
int atm,aroms,i,j,k;

if (!aromAtms[0])
	return;

aroms = 31 & (int)aromAtms[0];
for (i=2; i<=aroms; i++)
	{
	atm = aromAtms[i];
	for (j=1; j<i; j++)
		{
		for (k=0; k<mConnAtoms[mol][atm]; k++)
			{
			if (mConnAtom[mol][atm][k] == aromAtms[j])
				{
				mBondFlags[mol][mConnBond[mol][atm][k]] |= cBondAromatic;
				if ((aroms & 0x03) == 2)
					mBondFlags[mol][mConnBond[mol][atm][k]] |= cBondAromatic4N2;
				break;
				}
			}
		}
	}
}*/


	/**
	 * find reactive strand from alpha carbon of educt!!
	 * @return
	 */
	private int findStrand(int[] strand) {
		int strLngth = 1;
		int m = strand[0];
		StereoMolecule mol = mRxn.getMolecule(m);

		if (gErrout != null)
			try { gErrout.write("strandmol:"+mol+"; Atms:"+strand[1]+","); } catch (IOException ioe) {}
		
		for (int carbon=1; carbon<=7; carbon++) {
			int actlAtm = strand[carbon];
			boolean found = false;
			for (int connectd=0; connectd<mol.getConnAtoms(actlAtm); connectd++) {
				int tstAtm = mol.getConnAtom(actlAtm, connectd);
		
				if (mAtomType[m][tstAtm] != CARBON) continue; // reject nonCarbon
		
				if ((mAtomFlags[m][tstAtm] & cAtomChanges) == 0)	// reject if not changing
					continue;
		
		// added to terminate strands at bonds which are aromatic in educt and product.
		// TLS 24.6.95
				int tstBnd = mol.getConnBond(actlAtm, connectd);
				int proM = mCorProd[m][actlAtm];
				int proAtm = mCorAtom[m][actlAtm];
				StereoMolecule proMol = mRxn.getMolecule(proM);
				int proTstAtm = -1;
				int proTstBnd = -1;
				for (int i=0; i<proMol.getConnAtoms(proAtm); i++) {
					proTstAtm = proMol.getConnAtom(proAtm, i);
					if (proMol.getAtomMapNo(proTstAtm) != mol.getAtomMapNo(tstAtm))
						continue;
					proTstBnd = proMol.getConnBond(proAtm, i);
					break;
					}
				if (proTstBnd == -1)
					continue;					// reject when no bond in product
												// important for strandrecognition
												// in fragmentation reactions
		
				if (mol.getBondOrder(tstBnd) == proMol.getBondOrder(proTstBnd))
					continue;					// reject when order of bond to
												// adjacent atom not changing
		
				if (mol.isDelocalizedBond(tstBnd)
				 && proMol.isDelocalizedBond(proTstBnd)) {
					boolean found2 = false;
					for (int i=0; i<mol.getConnAtoms(actlAtm); i++) {
						int connAtm = mol.getConnAtom(actlAtm, i);
						int proConnM = mCorProd[m][connAtm];
						int proConnAtm = mCorAtom[m][connAtm];
						StereoMolecule proConnMol = mRxn.getMolecule(proConnM);
						if (proConnAtm == 255) continue;
						int j = 0;
						for (; j<proMol.getConnAtoms(proAtm); j++)
							if (proMol.getConnAtom(proAtm, j) == proConnAtm)
								break;
						if (j == proMol.getConnAtoms(proAtm)) continue;
		
						if (mol.getConnBondOrder(actlAtm, i) == 2
						 || proMol.getConnBondOrder(proAtm, j) == 2) {
							if (mol.isAromaticAtom(connAtm)
							 && proConnMol.isAromaticAtom(proConnAtm))
								continue;
							found2 = true;
							}
						}
		
					for (int i=0; i<mol.getConnAtoms(tstAtm); i++)
						{
						int connAtm = mol.getConnAtom(tstAtm, i);
						int proConnM = mCorProd[m][connAtm];
						int proConnAtm = mCorAtom[m][connAtm];
						StereoMolecule proConnMol = mRxn.getMolecule(proConnM);
						if (proConnAtm == 255) continue;
						int j = 0;
						for (; j<proMol.getConnAtoms(proTstAtm); j++)
							if (proMol.getConnAtom(proTstAtm, j) == proConnAtm)
								break;
						if (j == proMol.getConnAtoms(proTstAtm)) continue;
		
						if (mol.getConnBondOrder(tstAtm, i) == 2
						 || proMol.getConnBondOrder(proTstAtm, j) == 2) {
							if (mol.isAromaticAtom(connAtm)
							 && proConnMol.isAromaticAtom(proConnAtm))
								continue;
							found2 = true;
							}
						}
		
					if (!found2) continue;
					}
						/* reject if bond is aromatic in both educt and product     */
						/* however accept if at least one of both atoms of bond in  */
						/* question has in educt or product a double bond to an     */
						/* atom which is not aromatic in both educt and product     */
		// end of addition
		
		
		/*		EBndOrdr = PBndOrdr = 0;			/* reject when order of bond to **
				for (i=0; i<MAXCONNS; i++)					/* adjacent atom not changing   **
					if (mConnMpNo[mol][actlAtm][i] == mRxn.getAtomMappingNo( mol, tstAtm ))
						EBndOrdr++;
				for (i=0; i<MAXCONNS; i++)
					if (mConnMpNo[mCorProd[mol][actlAtm]]
								[mCorAtm[mol][actlAtm]][i]
						== mRxn.getAtomMappingNo( mol, tstAtm )) PBndOrdr++;
				if (EBndOrdr == PBndOrdr) continue;
		
				if (PBndOrdr == 0) continue;	/* reject when no bond in product  **
												/* important for strandrecognition **
												/* in fragmentation reactions	   **
		this was the old part which was replaced and enhanced by the addition above */
		
				boolean instrand = false;
				for (int i=1; i<carbon; i++) {	// reject if already in strand
					if (tstAtm == strand[i]) {
						instrand = true;
						break;
						}
					}
				if (instrand) continue;
		
				if (found)
					return 0;		// strand is branched
				else {
					if (strLngth == 7)
						return -1; 	// strand length exceeds 7 carbons
					found = true;
					strand[carbon+1] = tstAtm;

					if (gErrout != null)
						try { gErrout.write(","+strand[carbon+1]); } catch (IOException ioe) {}

					strLngth++;
					}
				}
			if (!found) return strLngth;
			}
		return 0; 	// should never reach this
		}


	/**
	 * returns number and types of noncarbon-nonhydrogen substituents at carbon A
	 * that happen not to appear at carbon B (determines incoming or leaving groups)
	 * @param m
	 * @param atm
	 * @param proM
	 * @param proAtm
	 * @return
	 */
	private int leaving(int m, int atm, int proM, int proAtm, int[] group, int[] flshAtm, int[] stereo) {
		boolean[] isUnchangedEType = new boolean[8];
		boolean[] matched = new boolean[8];
		int[] associatedPIndex = new int[8];
		int[][] eSecCrit = new int[8][2];
		int[][] pSecCrit = new int[8][2];
		int[] eAtm = new int[8];
		int[] pAtm = new int[8];
		int[] eType = new int[8];
		int[] pType = new int[8];
		int[] eMapNo = new int[8];
		int[] pMapNo = new int[8];

		int nrofETypes = getFGrpTypes( m, atm, eType, eAtm, eMapNo, eSecCrit );
		int nrofPTypes = getFGrpTypes( proM, proAtm, pType, pAtm, pMapNo, pSecCrit );
		
		stereo[0] = 0;
		int counter = 0;
		
		for (int i=0; i<nrofETypes; i++) {
			int j = 0;
			for (; j<nrofPTypes; j++) {
				if (matched[j]) continue;
				if (eMapNo[i] != pMapNo[j]) continue;
				if (eSecCrit[i][0] == 0)
					{ if (eType[i] != pType[j]) continue; }
				else
					{ if (eSecCrit[i][0] != pSecCrit[j][0] || eSecCrit[i][1] != pSecCrit[j][1]) continue; }
		
				matched[j] = true;
				isUnchangedEType[i] = true;
				associatedPIndex[i] = j;
				break;
				}
		
			if (j == nrofPTypes) {
				flshAtm[counter] = eAtm[i];
				group[counter++] = eType[i];
				if (counter == 4) return 4;	// neglect stereo center check here
				}
			}

		StereoMolecule mol = mRxn.getMolecule(m);
		StereoMolecule proMol = mRxn.getMolecule(proM);

		if (mol.isAtomStereoCenter(atm)		// stereo center in pseudo educt or pseudo product
		 || proMol.isAtomStereoCenter(proAtm)) {
			if (mol.getConnAtoms(atm) == proMol.getConnAtoms(proAtm)) {
				stereo[0] = askInversion(m,atm,proM,proAtm);
				if (stereo[0] != 3 && counter == 0)	{  // if not retention and no prior leaving group detected
					boolean found = false;
					for (int i=0; i<8; i++)
						if (matched[i])
							found = true;
					if (found) {	 // at least one group of same type found
						for (int i=0; i<nrofETypes; i++) {
							if (isUnchangedEType[i]) {
								mAtomFlags[m][eAtm[i]] |= cAtomChanges;
								mAtomFlags[m][pAtm[associatedPIndex[i]]] |= cAtomChanges;
								flshAtm[counter] = eAtm[i];
								group[counter++] = eType[i];
								}
							}
						}
					}
				}
			}
		
		int eCharge = mol.getAtomCharge(atm);
		int pCharge = proMol.getAtomCharge(proAtm);
		if (eCharge < 0) eCharge = 0;	// don't consider negative charges here
		if (pCharge < 0) pCharge = 0;
		if (eCharge > 0)
			if (eCharge > pCharge)
				for (int i=0; i<eCharge-pCharge; i++)
					{
					flshAtm[counter] = atm;
					group[counter++] = 208;	/* carbenium */
					if (counter == 4) return 4;
					}
		
		int leavingH = pCharge - eCharge;
		for (int i=0; i<proMol.getConnAtoms(proAtm); i++)	// calculate leaving hydrogens
			leavingH += proMol.getConnBondOrder(proAtm, i);
		for (int i=0; i<mol.getConnAtoms(atm); i++)
			leavingH -= mol.getConnBondOrder(atm, i);

		for (int i=0; i<leavingH; i++) { // leaving Hydrogen
			flshAtm[counter] = atm;
			group[counter++] = 2;
			if (counter == 4) return 4;
			}
		
		return counter;
		}


	/**
     * (if atm1 belongs to product and atm2 to educt: retval 1 = racemisation)
     * This routine compares a stereocenter in the starting material with it's
     * counterpart in the products. Relying on the supplied parity information
     * it is calculating if retention or conversion occures. In order to do so
     * it sorts first the numbers of atoms connected to atm1 in ascending order.
     * Then it determines the corresponding atoms in mol2 and lists them in the
     * same order. Herefore corresponding[] may provide an atom number
     * corresponding to mConnAtom[mol1][atm1][?]. Otherwise this value must
     * equal 255 which causes this routine to find the corresponding atom by
     * comparing mRxn.mAtomMapNo[][] of both molecules. If a corresponding
     * atom can not be found zero is returned.
     * Then the list of corresponding atoms is converted in such a way that the
     * lowest atom number becomes 0, the next higher one 1 and so on. This list
     * (each entry 2 bit) serves as an index for lookup tables that define, if
     * the same parity value in educt and product means inversion or retention.
	 * @param m1
	 * @param atm1
	 * @param m2
	 * @param atm2
	 * @return 0: no info available; 1: creation of stereochemistry; 2:inversion; 3: retention
	 */
	private int askInversion( int m1, int atm1, int m2, int atm2 ) {
		final int[] parityChange3 =
		 { 2,2,2,2, 2,2,0,2, 2,1,2,2, 2,2,2,2, 2,2,1,2, 2,2,2,2, 0,2,2,2, 2,2,2,2,
		   2,0,2,2, 1,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2 };
		final int[] parityChange4 =
		 { 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,0, 2,2,1,2,
		   2,2,2,2, 2,2,2,1, 2,2,2,2, 2,0,2,2, 2,2,2,2, 2,2,0,2, 2,1,2,2, 2,2,2,2,
		   2,2,2,2, 2,2,2,2, 2,2,2,1, 2,2,0,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2,
		   2,2,2,0, 2,2,2,2, 2,2,2,2, 1,2,2,2, 2,2,1,2, 2,2,2,2, 0,2,2,2, 2,2,2,2,
		   2,2,2,2, 2,2,2,0, 2,2,2,2, 2,1,2,2, 2,2,2,1, 2,2,2,2, 2,2,2,2, 0,2,2,2,
		   2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,0,2,2, 1,2,2,2, 2,2,2,2, 2,2,2,2,
		   2,2,2,2, 2,2,1,2, 2,0,2,2, 2,2,2,2, 2,2,0,2, 2,2,2,2, 1,2,2,2, 2,2,2,2,
		   2,1,2,2, 0,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2 };
		int[] proCorres = new int[MAXCONNS];

		StereoMolecule mol1 = mRxn.getMolecule(m1);
		StereoMolecule mol2 = mRxn.getMolecule(m2);

		if (mol1.getConnAtoms(atm1) != mol2.getConnAtoms(atm2)) return 0;
		if (mol1.getConnAtoms(atm1) < 3 || mol1.getConnAtoms(atm1) > 4) return 0;
		
		if (!mol2.isAtomStereoCenter(atm2)) return 0;
		if (!mol1.isAtomStereoCenter(atm1)) return 1;

		// for (i=0; i<MAXCONNS; i++) proCorres[i] = corresponding[i]; changed to:
		for (int i=0; i<MAXCONNS; i++) proCorres[i] = 255;

		int alreadyUsed = 255;
		for (int i=0; i<mol1.getConnAtoms(atm1); i++) {	//fullfill list of corresponding atoms in mol2
			if (proCorres[i] == 255) {
				if (mol1.getAtomMapNo(mol1.getConnAtom(atm1, i)) != 0) {
					for (int j=0; j<mol2.getConnAtoms(atm2); j++) {
						if (mol1.getAtomMapNo(mol1.getConnAtom(atm1, i))
						 == mol2.getAtomMapNo(mol2.getConnAtom(atm2, j))) {
							proCorres[i] = mol2.getConnAtom(atm2, j);
							break;
							}
						}
					if (proCorres[i] == 255) {
						if (alreadyUsed != 255) return 0;
						alreadyUsed = i;
						}
					}
				else {	// tolerate one unmapped substituent in educt and product
					if (alreadyUsed != 255) return 0;
					alreadyUsed = i;
					}
				}
			}
		
		if (alreadyUsed != 255) {
			for (int i=0; i<mol2.getConnAtoms(atm2); i++) {
				boolean found = false;
				for (int j=0; j<mol1.getConnAtoms(atm1); j++) {
					if (proCorres[j] == mol2.getConnAtom(atm2, i)) {
						found = true;
						break;
						}
					}
				if (!found) {
					proCorres[alreadyUsed] = mol2.getConnAtom(atm2, i);
					break;
					}
				}
			}
		
		int index = 0;								// prepare index for lookup tables
		for (int i=0; i<mol1.getConnAtoms(atm1); i++) {
			index <<= 2;
			int position = 0;
			for (int j=0; j<mol1.getConnAtoms(atm1); j++) {
				if (proCorres[j] < proCorres[i]) position++;
				}
			index |= position;
			}
		
		int lookUpResult = (mol1.getConnAtoms(atm1) == 3) ?
				parityChange3[index] : parityChange4[index];

		if (lookUpResult == 2) return 0;		// should not happen
		
		if (mol1.getAtomParity(atm1) == mol2.getAtomParity(atm2))
			return 3 - lookUpResult;
		else
			return 2 + lookUpResult;
		}


	private int getFGrpTypes(int m, int atm, int[] type, int[] grpAtm, int[] grpMapNo, int[][] secCrit) {

		StereoMolecule mol = mRxn.getMolecule(m);

		int nrofTypes = 0;
		for (int i=0; i<mol.getConnAtoms(atm); i++) {
			int connAtm = mol.getConnAtom(atm, i);
			if (mAtomType[m][connAtm] == CARBON) continue;
			int fGrpBndOrdr = mol.getConnBondOrder(atm, i);
			if (nrofTypes + fGrpBndOrdr > 8) continue;  // not likely to happen
			grpAtm[nrofTypes] = connAtm;
			grpMapNo[nrofTypes] = mol.getAtomMapNo(connAtm);
			type[nrofTypes] = gettyp( m, connAtm, atm, fGrpBndOrdr, secCrit[nrofTypes] );
			nrofTypes++;
			for (int j=1; j<fGrpBndOrdr; j++) {
				type[nrofTypes] = type[nrofTypes - 1];
				grpAtm[nrofTypes] = grpAtm[nrofTypes - 1];
				grpMapNo[nrofTypes] = grpMapNo[nrofTypes - 1];
				secCrit[nrofTypes][0] = secCrit[nrofTypes - 1][0];
				secCrit[nrofTypes][1] = secCrit[nrofTypes - 1][1];
				nrofTypes++;
				}
			}
		
		return nrofTypes;
		}


	/**
     * gettyp() returns basically the 8 bit type of a leaving group.
     *
     * However, since some tasks that inquire leaving groups require more information,
     * additional information - specific for the type of atom Y is returned in secCrit.
     * If the secCrit value of an educts' functional group equals the secCrit of a
     * functional group at the appropriate product atom (and is not 0!), then shall be
     * assumed that both groups are equivalent, that means that different gettyp()-
     * return values are due to remote changes and shall not be considered as a local
     * substitution.
     * In the case Y=oxygen secCrit contains the mapping number of atom Z.
     * In the case Y=nitrogen secCrit contains the mapping number of atom Y as well as
     * the bondorder between A and Y.
     * By this information a routine may investigate, if some slight changes of
     * attached groups happen due to local substitution reactions or are simply remote
     * red/ox or substitution reactions.
     *
     * Example:  C-C-O-CO-Me  -->  C-C-O-CH2-Me
     * comparing the mapping numbers of CO and CH2 answers if a real substitution of
     * acetate by ethanolate or a reduction takes place !
     * @param atm atom Y
     * @param xAtm atom A
     * @param bndOrdr bond order A-Y
	 * @return type of leaving group connected to atom A in A-Y[-Z]
	 */
	private int gettyp(int m, int atm, int xAtm, int bndOrdr, int[] secCrit ) {
		int connAtm,i;
		final int[] fGrpClass = { 0,
		  2,																  0,
		 40, 48,										 56,  0,  0,  0,216,  0,
		 36, 42,										 50, 58,  0,  0,220,  0,
		 37, 44, 64, 68, 70, 80, 84, 85, 86, 87, 72, 76, 52, 59,232,  0,222,  0,
		 38, 45, 65, 69, 71, 82,  0, 89, 90, 91, 73, 77, 54, 60, 62,  0,223,  0,
		 39, 46, 66,	 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
					 69, 71, 83, 92, 93, 94, 95, 74, 78, 55, 61, 63,  0,  0,  0,
		  0, 47, 67,	114,115,116,117,118,119,120,121,122,123,124,125,126,127,
		  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
		  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
		  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  6,222,  0,  0,  0,  0,
		  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
		  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
		  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 };
		
		secCrit[0] = 0;	// default for uncomplicated leaving groups
		secCrit[1] = 0;

		StereoMolecule mol = mRxn.getMolecule(m);
		if (mol.getAtomicNo(atm) == 7) {	// nitrogen functionalities
			if (mol.getAtomMapNo(atm) != 0) {
				secCrit[0] = mol.getAtomMapNo(atm);
				secCrit[1] = bndOrdr;
				}
		
			if (bndOrdr == 3) return 184;					// nitrile
		
			if (bndOrdr == 2) {
				if (mol.getAtomCharge(atm) == 1) {	// positive charge at N
					switch (mZ[m][atm]) {
						case 0: return 176;						// immonium
						case 1: return 178;						// nitrones etc.
						default: return 179;					// iso nitro etc.
						}
					}
				else
					{
					switch (mZ[m][atm]) {
						case 0: return 180;						// imine
						default: return 182;					// oximes etc.
						}
					}		
				}
			else
				{
				if (mol.getAtomCharge(atm) == 1) {	// positive charge at N
					if (mZ[m][atm] == 1) return 172;		// amine oxide
					if (mZ[m][atm] > 1) return 174;		// nitro etc.
					switch (mPi[m][atm]) {
						case 0: return 168;						// ammonium
						case 1: return 170;						// immonium
						default: return 171;					// isonitrile
						}
					}
		
				if (mZ[m][atm] == 0) {					// no hetero atoms at N
					if (mPi[m][atm] == 0) {				// no double bonds to C
						switch (mSigma[m][atm]) {
							case 1: return 160;					// prim amine
							case 2: return 161;					//  sec amine
							default: return 162;				// tert amine
							}
						}
					else return 163;						// imine -N=C
					}
				else if (mZ[m][atm] == 1) {				// z=1 at N
					return 164;
					}
				else {										// z=2 at N
					return 166;
					}
				}
			}
		else if (mol.getAtomicNo(atm) == 8) {	// oxygen functionalities
			if (bndOrdr == 2) return 156;			// carbonyl

			connAtm = -1;
			for (i=0; i<mol.getConnAtoms(atm); i++)
				if (mol.getConnAtom(atm, i) != xAtm)
					connAtm = mol.getConnAtom(atm, i);
		
			if (connAtm == -1) return 152;			// alcohol
		
			if (mol.getAtomMapNo(connAtm) != 0) {
				secCrit[0] = mol.getAtomMapNo(connAtm);
				secCrit[1] = 8;
				}
		
			if (mAtomType[m][connAtm] == CARBON) {	// esters, ethers
				if (mPi[m][connAtm] == 0) {			// saturated alpha atom
					switch (mZ[m][connAtm]) {
					case 1:							// ethers
						for (i=0; i<mol.getConnAtoms(connAtm); i++)
							if (mPi[m][mol.getConnAtom(connAtm, i)] != 0)
								return 132;   		// allyloxy / benzyloxy
		
						switch(mSigma[m][connAtm]) {
							case 0: return 128;				// methoxy
							case 1: return 129;				// prim alkoxy
							case 2: return 130;				// sec alkoxy
							case 3: return 131;				// tert alkoxy
							}
					case 2:							// acetals
						switch(mSigma[m][connAtm]) {
							case 0: return 136;				// -O(CH2)X
							case 1: return 137;				// -O(CHR)X
							case 2: return 138;				// -O(CR2)X
							}
					case 3:							// ester etc.
						switch(mSigma[m][connAtm]) {
							case 0: return 140;				// -O(CH)X2
							case 1: return 141;				// -O(CR)X2
							}
					case 4:							// carbonates etc.
						return 142;
					default:
						return 0;	// shouldn't happen
						}
					}
				else {				// enol and phenol ethers, keten acetals
					if (mol.isAromaticAtom(connAtm))
						return 134;					// phenoxy
					else
						return 135;					// alkenoxy
					}
				}
			else {
				switch (mol.getAtomicNo(connAtm)) {
					case 7 : return 145;					// -O-N
					case 8 : return 144;					// -O-O
					case 14: return 154;					// -O-Si
					case 15: return 147;					// -O-P
					case 16: return 146;					// -O-S
					default: return 148;
					}
				}
			}
		else if (mol.getAtomicNo(atm) == 15) {	// phosphor functionalities
			if (mol.getAtomMapNo(atm) != 0) {
				secCrit[0] = mol.getAtomMapNo(atm);
				secCrit[1] = 15;
				}
		
			if (mZ[m][atm] == 0) {
				if (mol.getAtomCharge(atm) != 0)	return 226; // phosphonium
				else								return 224; // phosphide
				}
			if (mZ[m][atm] <= 2) return 228;
			return 230;
			}
		else if (mol.getAtomicNo(atm) == 16) {	// sulfur functionalities
			if (mol.getAtomMapNo(atm) != 0) {
				secCrit[0] = mol.getAtomMapNo(atm);
				secCrit[1] = 16;
				}
		
			if (bndOrdr == 2) return 249;
			switch (mZ[m][atm]) {
				case 0: return 248;
				case 1: return 242;
				case 2: return 244;
				case 4: return 246;
				default: return 240;
				}
			}
		else if (mol.getAtomicNo(atm) == 34) {	// selenium functionalities
			if (mol.getAtomMapNo(atm) != 0) {
				secCrit[0] = mol.getAtomMapNo(atm);
				secCrit[1] = 34;
				}
		
			switch (mZ[m][atm]) {
				case 0: return 250;
				case 1: return 243;
				case 2: return 245;
				case 4: return 247;
				default: return 241;
				}
			}
		else if (mol.getAtomicNo(atm) == 52) {	// tellurum functionalities
			if (mol.getAtomMapNo(atm) != 0) {
				secCrit[0] = mol.getAtomMapNo(atm);
				secCrit[1] = 52;
				}
		
			if (mZ[m][atm] == 0) return 252;
			else return 254;
			}
		else
			{
			return fGrpClass[mol.getAtomicNo(atm)];
			}
		}

	private void markReactionCenters() {
		for (int m=0; m<mRxn.getMolecules(); m++) {
			StereoMolecule mol = mRxn.getMolecule(m);
			for (int atm=0; atm<mol.getAtoms(); atm++)
				mol.setAtomMarker(atm, (mAtomFlags[m][atm] & cAtomChanges) != 0);
			}
		}
	}
