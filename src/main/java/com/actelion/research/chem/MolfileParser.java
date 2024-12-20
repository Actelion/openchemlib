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

/*
 *       Date            User        Remark
 *       ==========      =========   ===========================================
 *       02/08/2002      CXR         Handle the chiral flag now
 *       12/12/2002      TLS         assumes non-stereo bond in case of missing stereo bond info
 *       02/18/2003      TLS         atom and bond query features added
 *       04/27/2006      TLS         added support for molfile version 3.0
 *       02/22/2007      CXR         Handle Atoms lists in V3 Molfiles
 *       02/07/2011      TLS         added assignment of stereochemical group to bonds as Actelion specific extension to MDL V3 format
 *
 */
package com.actelion.research.chem;

import com.actelion.research.io.BOMSkipper;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.TreeMap;

public class MolfileParser
{
	public static final int MODE_KEEP_HYDROGEN_MAP = 1;
	public static final int ALLOWED_ATOM_LABELS = Molecule.cPseudoAtomsHydrogenIsotops
												| Molecule.cPseudoAtomsRGroups
												| Molecule.cPseudoAtomsAminoAcids;
	public static final int ALLOWED_ATOM_LABELS_IN_LIST = Molecule.cPseudoAtomsHydrogenIsotops;

	public static boolean debug = false;
	private StereoMolecule mMol;
	private TreeMap<Integer,Integer> mAtomIndexMap,mBondIndexMap;
	private boolean mTreatAnyAsMetalBond,mDeduceMissingCharges,mChiralFlag,mIsV3000,mAssumeChiralTrue;
	private int mMode;
	private int[] mHydrogenMap;

	/**
	 * Constructor of a MolFileParser, which will mirror Y,Z coordinates
	 */
	public MolfileParser() {
		mMode = 0;
	}


	public MolfileParser(int mode) {
		mMode = mode;
	}


	/**
	 * If this MoflileParser was instantiated with MODE_KEEP_HYDROGEN_MAP
	 * @return
	 */
	public int[] getHandleHydrogenMap() {
		return mHydrogenMap == null ? mMol.getHandleHydrogenMap() : mHydrogenMap;
	}


	private boolean readMoleculeFromBuffer(BufferedReader reader) {
		int[] valence = null;   // Some toolkit (RDKit) set vvv (valence) for charge atoms with normal valence
								// Thus, we must check, whether we just have a charged normal valence atom
								// before setting an abnormal valence for the atom.

		try{
			String line;
			int natoms,nbonds,nlists;

			mHydrogenMap = null;

			if(mMol != null){
				mMol.clear();
				mMol.setFragment(false);
			}

			/*** Name line ***/
			String name = (line = reader.readLine());
			if(null == name){
				TRACE("readMoleculeFromBuffer: No Header Line\n");
				return false;
			}
			/*** User, program ***/
			if(null == (line = reader.readLine())){
				TRACE("Error [readMoleculeFromBuffer]: No Program Line\n");
				return false;
			}
			/*** Comment ***/
			if(null == (line = reader.readLine())){
				TRACE("Error [readMoleculeFromBuffer]: No Comment Line\n");
				return false;
			}

			mTreatAnyAsMetalBond = line.contains("From CSD data. Using bond type 'Any'");
			mDeduceMissingCharges = line.contains("From CSD data.");

			/*** Counts line ***/
			if(null == (line = reader.readLine())){
				TRACE("Error [readMoleculeFromBuffer]: No Counts Line\n");
				return false;
			}

			mIsV3000 = false;
			mChiralFlag = mAssumeChiralTrue;
			try{
				natoms = Integer.parseInt(line.substring(0,3).trim());
				nbonds = Integer.parseInt(line.substring(3,6).trim());
				nlists = parseIntOrSpaces(line.substring(6,9).trim());
				mChiralFlag |= (1 == parseIntOrSpaces(line.substring(12,15).trim()));
				mIsV3000 = (line.length() >= 39 && line.startsWith("V3000", 34));
			} catch(Exception e){
				TRACE("Warning [readMoleculeFromBuffer]: Unable to interpret counts line\n");
				return false;
			}

			if(mIsV3000){
				boolean res = readMoleculeV3FromBuffer(reader);
				mMol.setName(name);
				return res;
			}

			if(mMol == null){
				mMol = new StereoMolecule(natoms,nbonds);
			}

			mMol.setName(name);

			if(!mChiralFlag){
				mMol.setToRacemate();
			}

			/*** Handle special case of natoms = 0 ***/
			if(0 == natoms){
				while(line != null && (!(line.equals("M  END") || line.equals("$$$$") || line.substring(1).equals("$")))){
					line = reader.readLine();
				}
				return true;
			}

			for(int i = 0;i < natoms;i++){
				if(null == (line = reader.readLine())){
					TRACE("Error [readMoleculeFromBuffer]: No Atom Line\n");
					return false;
				}

				float x = Float.parseFloat(line.substring(0,10).trim());
				float y = Float.parseFloat(line.substring(10,20).trim());
				float z = Float.parseFloat(line.substring(20,30).trim());

				int atom = mMol.addAtom(x, -y, -z);

				String label = line.substring(31,34).trim();
				if(label.equals("A") || label.equals("*")){
					mMol.setAtomQueryFeature(atom,Molecule.cAtomQFAny,true);
				} else if(label.equals("Q")) {    // 'Q' is defined as 'unspecified' for V2000; we use V3000 behaviour (any but not C,H)
					int[] list = new int[1];
					list[0] = 6;
					mMol.setAtomList(atom, list, true);
				} else {
					int atomicNo = Molecule.getAtomicNoFromLabel(label, ALLOWED_ATOM_LABELS);
					mMol.setAtomicNo(atom,atomicNo);
				}

				int massDif = parseIntOrSpaces(line.substring(34,36).trim());
				if(massDif != 0){
					mMol.setAtomMass(atom,Molecule.cRoundedMass[mMol.getAtomicNo(atom)] + massDif);
				}

				int chargeDif = parseIntOrSpaces(line.substring(36,39).trim());
				if(chargeDif != 0){
					if (chargeDif == 4)
						mMol.setAtomRadical(atom, Molecule.cAtomRadicalStateD);
					else
						mMol.setAtomCharge(atom,4 - chargeDif);
				}

				int mapNo = (line.length() < 63) ? 0 : parseIntOrSpaces(line.substring(60,63).trim());
				mMol.setAtomMapNo(atom,mapNo,false);

				//parity = parseIntOrSpaces(line.substring(39, 42).trim());

				int hCount = (line.length() < 45) ? 0 : parseIntOrSpaces(line.substring(42,45).trim());
				switch(hCount){
					case 0:
						break;
					case 1: // no hydrogen
						mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot1Hydrogen
						                             | Molecule.cAtomQFNot2Hydrogen, true);
						break;
					case 2: // at least 1 hydrogen
						mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot0Hydrogen, true);
						break;
					case 3: // at least 2 hydrogens
                        mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot0Hydrogen
                                                     | Molecule.cAtomQFNot1Hydrogen, true);
                        break;
					default: // at least 3,4 hydrogens
						mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot0Hydrogen
						                             | Molecule.cAtomQFNot1Hydrogen
												     | Molecule.cAtomQFNot2Hydrogen, true);
						break;
				}

				if(line.length() >= 48 && line.charAt(47) == '1') {
					mMol.setAtomQueryFeature(atom,Molecule.cAtomQFMatchStereo,true);
				}

                int v = (line.length() < 51) ? 0 : parseIntOrSpaces(line.substring(48,51).trim());
				if (v != 0) {
					if (valence == null)
						valence = new int[natoms];
					valence[atom] = v;
				}
			}

			// Loop all the bonds , read the bond record and generate
			// the internal representation
			for(int i = 0;i < nbonds;i++){
				if(null == (line = reader.readLine())){
					TRACE("Error [readMoleculeFromBuffer]:No Bond Line\n");
					return false;
				}

				int atom1 = Integer.parseInt(line.substring(0,3).trim()) - 1;
				int atom2 = Integer.parseInt(line.substring(3,6).trim()) - 1;
				int bondType = Integer.parseInt(line.substring(6,9).trim());
				int stereo = (line.length() < 12) ? 0 : parseIntOrSpaces(line.substring(9,12).trim());
				int topology = (line.length() < 18) ? 0 : parseIntOrSpaces(line.substring(15,18).trim());

				if (bondType == 8
				 && (mTreatAnyAsMetalBond
				  || mMol.isMetalAtom(atom1)
				  || mMol.isMetalAtom(atom2)))
					bondType = 9;      // metal ligand bond doesn't exist in molfile version 2

				buildBond(atom1,atom2,bondType,stereo,topology);
			}

			// skip atom list block
			for(int i = 0;i < nlists;i++){
				if(null == (line = reader.readLine())){
					TRACE("Error [readMoleculeFromBuffer]: No List Line\n");
					return false;
				}
			}

			/********************************************************************
			 ***  Check for "M  CHG" charge record or "M  ISO" isomer record.
			 ***  --> Must have "M  END" or "$$$$" at end of molecule !
			 ********************************************************************/
			if(null == (line = reader.readLine())){
				TRACE("Error ReadMoleculeFromBuffer Missing M END or $$$$\n");

				if ((mMode & MODE_KEEP_HYDROGEN_MAP) != 0)
					mHydrogenMap = mMol.getHandleHydrogenMap();

				handleValences(valence);

				// to run the racemization scheduled with mMol.setToRacemate()
				if(!mChiralFlag)
					mMol.ensureHelperArrays(Molecule.cHelperParities);

				return true;
			}

			while(line != null && (!(line.equals("M  END") || line.equals("$$$$")))){
				if(line.startsWith("M  CHG")){
					int aaa,vvv;
					int j = Integer.parseInt(line.substring(6,9).trim());
					if(j > 0){
						aaa = 10;
						vvv = 14;
						for(int k = 1;k <= j;k++,aaa += 8,vvv += 8){
							int atom = Integer.parseInt(line.substring(aaa,aaa + 3).trim()) - 1;
							int charge = Integer.parseInt(line.substring(vvv,vvv + 3).trim());
							mMol.setAtomCharge(atom,charge);
						}
					}
				}

				if(line.startsWith("M  ISO")){
					int aaa,vvv;
					int j = Integer.parseInt(line.substring(6,9).trim());
					if(j > 0){
						aaa = 10;
						vvv = 14;
						for(int k = 1;k <= j;k++,aaa += 8,vvv += 8){
							int atom = Integer.parseInt(line.substring(aaa,aaa + 3).trim()) - 1;
							int mass = Integer.parseInt(line.substring(vvv,vvv + 3).trim());
							mMol.setAtomMass(atom,mass);
						}
					}
				}

				if(line.startsWith("M  RAD")){
					int aaa,vvv;
					int j = Integer.parseInt(line.substring(6,9).trim());
					if(j > 0){
						aaa = 10;
						vvv = 14;
						for(int k = 1;k <= j;k++,aaa += 8,vvv += 8){
							int atom = Integer.parseInt(line.substring(aaa,aaa + 3).trim()) - 1;
							int radical = Integer.parseInt(line.substring(vvv,vvv + 3).trim());
							switch(radical){
								case 1:
									mMol.setAtomRadical(atom,Molecule.cAtomRadicalStateS);
									break;
								case 2:
									mMol.setAtomRadical(atom,Molecule.cAtomRadicalStateD);
									break;
								case 3:
									mMol.setAtomRadical(atom,Molecule.cAtomRadicalStateT);
									break;
							}
						}
					}
				}

				if(line.startsWith("M  RBC") || line.startsWith("M  RBD")){
					int j = Integer.parseInt(line.substring(6,9).trim());
					if(j > 0){
						int aaa = 10;
						int vvv = 14;
						for(int k = 1;k <= j;k++,aaa += 8,vvv += 8){
							int atom = Integer.parseInt(line.substring(aaa,aaa + 3).trim()) - 1;
							int ringState = Integer.parseInt(line.substring(vvv,vvv + 3).trim());
							switch(ringState){
								case -1:
									mMol.setAtomQueryFeature(atom,
										Molecule.cAtomQFNot2RingBonds
										| Molecule.cAtomQFNot3RingBonds
										| Molecule.cAtomQFNot4RingBonds,
										true);
									break;
								case 1:
									mMol.setAtomQueryFeature(atom,
										Molecule.cAtomQFNotChain,
										true);
									break;
								case 2:
									mMol.setAtomQueryFeature(atom,
										Molecule.cAtomQFNotChain
										| Molecule.cAtomQFNot3RingBonds
										| Molecule.cAtomQFNot4RingBonds,
										true);
									break;
								case 3:
									mMol.setAtomQueryFeature(atom,
										Molecule.cAtomQFNot2RingBonds
										| Molecule.cAtomQFNot3RingBonds
										| Molecule.cAtomQFNot4RingBonds,
										true);
									break;
								case 4:
									mMol.setAtomQueryFeature(atom,
										Molecule.cAtomQFNotChain
										| Molecule.cAtomQFNot2RingBonds
										| Molecule.cAtomQFNot3RingBonds,
										true);
									break;
							}
						}
					}
				}

				// The Atom list is implemented as an int[] of atomic numbers.
				// NOT Lists are implemented as a sorted vector as negative Integers
				if(line.startsWith("M  ALS")){
					int atom = Integer.parseInt(line.substring(7,10).trim()) - 1;
					if(atom >= 0){
						int no = Integer.parseInt(line.substring(10,13).trim());
						boolean bNotList = (line.charAt(14) == 'T');
						int[] v = new int[no];
						int aaa = 16;
						for(int k = 0;k < no;k++,aaa += 4){
							String sym = line.substring(aaa,aaa + 4).trim();
							v[k] = Molecule.getAtomicNoFromLabel(sym, ALLOWED_ATOM_LABELS_IN_LIST);
						}
						mMol.setAtomicNo(atom, 6);
						mMol.setAtomList(atom,v,bNotList);
					}
				}

				if(line.startsWith("M  SUB")){
					int aaa,vvv;
					int j = Integer.parseInt(line.substring(6,9).trim());
					if(j > 0){
						aaa = 10;
						vvv = 14;
						for(int k = 1;k <= j;k++,aaa += 8,vvv += 8){
							int atom = Integer.parseInt(line.substring(aaa,aaa + 3).trim()) - 1;
							int substitution = Integer.parseInt(line.substring(vvv,vvv + 3).trim());
							if(substitution == -2){
								mMol.setAtomQueryFeature(atom,Molecule.cAtomQFNoMoreNeighbours,true);
							} else if(substitution > 0){
								int substitutionCount = 0;
								for(int bond = 0;bond < mMol.getAllBonds();bond++){
									if(mMol.getBondAtom(0,bond) == atom
									   || mMol.getBondAtom(1,bond) == atom){
										substitutionCount++;
									}
								}
								if(substitution > substitutionCount){
									mMol.setAtomQueryFeature(atom,Molecule.cAtomQFMoreNeighbours,true);
								}
							}
						}
					}
				}

				if(line.startsWith("M  RGP")){
					int aaa,vvv;
					int j = Integer.parseInt(line.substring(6,9).trim());
					if(j > 0){
						aaa = 10;
						vvv = 14;
						for(int k = 1;k <= j;k++,aaa += 8,vvv += 8){
							int atom = Integer.parseInt(line.substring(aaa,aaa + 3).trim()) - 1;
							int rno = Integer.parseInt(line.substring(vvv,vvv + 3).trim());
							if(rno >= 1 && rno <= 20){
								mMol.setAtomicNo(atom, Molecule.getAtomicNoFromLabel("R"+rno, Molecule.cPseudoAtomsRGroups));
							}
						}
					}
				}

				line = reader.readLine();
			}
		} catch(Exception e){
			e.printStackTrace();
			System.err.println("error reading molfile " + e);
			return false;
		}

		if (mDeduceMissingCharges) {
			introduceObviousMetalBonds();
			deduceMissingCharges();
		}

		// needs to be done for molfiles with chiral=0 that have stereo
		// centers which will be assigned to one ESR-AND group
		if ((mMode & MODE_KEEP_HYDROGEN_MAP) != 0)
			mHydrogenMap = mMol.getHandleHydrogenMap();

		handleValences(valence);

		mMol.ensureHelperArrays(Molecule.cHelperParities);

		return true;
	}

	/**
	 * Some software exports mol/sd-files with an unset chiral flag despite
	 * the original molecule is a pure enantiomer. This method allows the
	 * MolfileParser to override the molfile's chiral flag for such cases.
	 * @param b
	 */
	public void setAssumeChiralTrue(boolean b) {
		mAssumeChiralTrue = b;
	}

	/**
	 * @return whether the previous molfile's chiral flag was set
	 */
	public boolean isChiralFlagSet() {
		return mChiralFlag;
	}

	/**
	 * @return whether the previous molfile was a V3000
	 */
	public boolean isV3000() {
		return mIsV3000;
	}

	private void handleValences(int[] valence) {
		if (valence != null) {
			mMol.ensureHelperArrays(Molecule.cHelperNeighbours);
			for (int atom=0; atom<mMol.getAtoms(); atom++) {
				if (valence[atom] != 0) {
					int chargeCorrection = mMol.getElectronValenceCorrection(atom, mMol.getOccupiedValence(atom));
					if (valence[atom] == 15 ) {
						if (chargeCorrection >= 0)
							mMol.setAtomAbnormalValence(atom, 0);
					}
					else {
						if (valence[atom] != mMol.getMaxValence(atom))
							mMol.setAtomAbnormalValence(atom, valence[atom] - chargeCorrection);
					}
				}
			}
		}
	}

	private boolean readMoleculeV3FromBuffer(BufferedReader reader) throws IOException
	{
		final int MODE_CTAB = 1;
		final int MODE_CTAB_ATOM = 2;
		final int MODE_CTAB_BOND = 3;
		final int MODE_CTAB_COLLECTION = 4;

		if (mAtomIndexMap != null)
			mAtomIndexMap.clear();
		if (mBondIndexMap != null)
			mBondIndexMap.clear();

		int mode = 0;
		String line = reader.readLine();
		while(line != null && line.startsWith("M  V30 ")){
			line = line.substring(7).trim();
			while(line.endsWith("-")){
				String cont = reader.readLine();
				if(!cont.startsWith("M  V30 ")){
					return false;
				}
				line = line.substring(0,line.length() - 1).concat(cont.substring(7)).trim();
			}

			if(line.startsWith("BEGIN")){
				String modeString = line.substring(6).trim();
				if(modeString.startsWith("CTAB")){
					mode = MODE_CTAB;
				} else if(modeString.startsWith("ATOM")){
					mode = MODE_CTAB_ATOM;
				} else if(modeString.startsWith("BOND")){
					mode = MODE_CTAB_BOND;
				} else if(modeString.startsWith("COLLECTION")){
					mode = MODE_CTAB_COLLECTION;
				} else{
					TRACE("Error MolfileParser: Unsupported version 3 block\n");
					return false;
				}
			} else if(line.startsWith("END")){
				mode = 0;
			} else if(mode == MODE_CTAB){
				interpretV3CountLine(line);
			} else if(mode == MODE_CTAB_ATOM){
				interpretV3AtomLine(line);
			} else if(mode == MODE_CTAB_BOND){
				interpretV3BondLine(line);
			} else if(mode == MODE_CTAB_COLLECTION){
				interpretV3CollectionLine(line);
			} else{
				TRACE("Error MolfileParser: Unexpected version 3 line\n");
				return false;
			}

			line = reader.readLine();
		}

		while(line != null && (!(line.startsWith("M  END") || line.equals("$$$$")))){
			line = reader.readLine();
		}

		return true;
	}

	private void interpretV3CountLine(String line)
	{
		if(mMol == null){
			if(line.startsWith("COUNTS")){
				int index1 = 7;
				int index2 = indexOfNextItem(line,indexOfWhiteSpace(line,7));
				int natoms = Integer.parseInt(line.substring(index1,indexOfWhiteSpace(line,index1)));
				int nbonds = Integer.parseInt(line.substring(index2,indexOfWhiteSpace(line,index2)));
				mMol = new StereoMolecule(natoms,nbonds);
			}
		}
	}

	private void interpretV3AtomLine(String line) throws IOException
	{
		int index1 = 0;
		int index2 = endOfItem(line,index1);
		int atomIndex = Integer.parseInt(line.substring(index1,index2));

		index1 = indexOfNextItem(line,index2);
		index2 = endOfItem(line,index1);
		String label = line.substring(index1,index2);
//		System.out.println("Atom Index Line is " + line);
		int v[] = null;
		boolean bNotList = false;
		int l = isV3AtomList(line);
		if(l != 0) {
			v = interpretV3AtomList(line);
			if (l < 0)
				bNotList = true;				
			index2 = Math.abs(l);
		} 
		index1 = indexOfNextItem(line,index2);
		index2 = endOfItem(line,index1);
		float x = Float.parseFloat(line.substring(index1,index2));

		index1 = indexOfNextItem(line,index2);
		index2 = endOfItem(line,index1);
		float y = Float.parseFloat(line.substring(index1,index2));

		index1 = indexOfNextItem(line,index2);
		index2 = endOfItem(line,index1);
		float z = Float.parseFloat(line.substring(index1,index2));

		index1 = indexOfNextItem(line,index2);
		index2 = endOfItem(line,index1);
		int mapNo = Integer.parseInt(line.substring(index1,index2));

		int atom = mMol.addAtom(x, -y, -z);
		if(atom + 1 != atomIndex)
			mapAtomIndex(atomIndex, atom);

		if (v != null) {
			mMol.setAtomicNo(atom, 6);
			mMol.setAtomList(atom, v, bNotList);
		}

		if(mapNo != 0){
			mMol.setAtomMapNo(atom,mapNo,false);
		}

		if(label.equals("A") || label.equals("*")){
			mMol.setAtomQueryFeature(atom,Molecule.cAtomQFAny,true);
		} else if(label.equals("Q")){
			int[] list = new int[1];
			list[0] = 6;
			mMol.setAtomList(atom,list,true);
		} else{
			mMol.setAtomicNo(atom,Molecule.getAtomicNoFromLabel(label, ALLOWED_ATOM_LABELS));
		}

		while((index1 = indexOfNextItem(line,index2)) != -1){
			index2 = endOfItem(line,index1);
			String specifier = line.substring(index1,index2);
			int index = specifier.indexOf('=');
			String field = specifier.substring(0,index);
			int value = Integer.parseInt(specifier.substring(index + 1));
			if(field.equals("CHG")){
				mMol.setAtomCharge(atom,value);
			} else if(field.equals("RAD")){
				switch(value){
					case 1:
						mMol.setAtomRadical(atom,Molecule.cAtomRadicalStateS);
						break;
					case 2:
						mMol.setAtomRadical(atom,Molecule.cAtomRadicalStateD);
						break;
					case 3:
						mMol.setAtomRadical(atom,Molecule.cAtomRadicalStateT);
						break;
				}
			} else if(field.equals("CFG")){
				//  don't read parities from molfile, they are calculated from up/down bonds
				//  mMol.setAtomParity(atom, value, false);
			} else if(field.equals("MASS")){
				mMol.setAtomMass(atom,value);
            } else if(field.equals("VAL")){
                mMol.setAtomAbnormalValence(atom, (value==-1) ? 0 : (value==0) ? -1 : value);
			} else if(field.equals("HCOUNT")){
				switch(value){
					case 0:
						break;
					case -1: // no hydrogen
						mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot1Hydrogen
												     | Molecule.cAtomQFNot2Hydrogen
												     | Molecule.cAtomQFNot3Hydrogen, true);
						break;
					case 1: // at least 1 hydrogen
						mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot0Hydrogen, true);
						break;
                    case 2: // at least 2 hydrogen
                        mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot0Hydrogen
                                                     | Molecule.cAtomQFNot1Hydrogen, true);
                        break;
					default: // at least 3,4 hydrogens
						mMol.setAtomQueryFeature(atom, Molecule.cAtomQFNot0Hydrogen
												     | Molecule.cAtomQFNot1Hydrogen
												     | Molecule.cAtomQFNot2Hydrogen, true);
						break;
				}
			} else if(field.equals("SUBST")){
				if(value == -1){
					mMol.setAtomQueryFeature(atom,Molecule.cAtomQFNoMoreNeighbours,true);
				} else if(value > 0){
					int substitutionCount = 0;
					for(int bond = 0;bond < mMol.getAllBonds();bond++){
						if(mMol.getBondAtom(0,bond) == atom
						   || mMol.getBondAtom(1,bond) == atom){
							substitutionCount++;
						}
					}
					if(value > substitutionCount){
						mMol.setAtomQueryFeature(atom,Molecule.cAtomQFMoreNeighbours,true);
					}
				}
			} else if(field.equals("RBCNT")){
				switch(value){
					case -1:
						mMol.setAtomQueryFeature(atom,
												 Molecule.cAtomQFNot2RingBonds
												 | Molecule.cAtomQFNot3RingBonds
												 | Molecule.cAtomQFNot4RingBonds,
												 true);
						break;
					case 1:
						mMol.setAtomQueryFeature(atom,
												 Molecule.cAtomQFNotChain,
												 true);
						break;
					case 2:
						mMol.setAtomQueryFeature(atom,
												 Molecule.cAtomQFNotChain
												 | Molecule.cAtomQFNot3RingBonds
												 | Molecule.cAtomQFNot4RingBonds,
												 true);
						break;
					case 3:
						mMol.setAtomQueryFeature(atom,
												 Molecule.cAtomQFNot2RingBonds
												 | Molecule.cAtomQFNot3RingBonds
												 | Molecule.cAtomQFNot4RingBonds,
												 true);
						break;
					case 4:
						mMol.setAtomQueryFeature(atom,
												 Molecule.cAtomQFNotChain
												 | Molecule.cAtomQFNot2RingBonds
												 | Molecule.cAtomQFNot3RingBonds,
												 true);
						break;
				}
			} else{
				TRACE("Warning MolfileParser: Unused version 3 atom specifier:" + field + "\n");
			}
		}
	}

	private void interpretV3BondLine(String line) throws IOException
	{
		int index1 = 0;
		int index2 = endOfItem(line,index1);
		int bondIndex = Integer.parseInt(line.substring(index1,index2));

		index1 = indexOfNextItem(line,index2);
		index2 = endOfItem(line,index1);
		int bondType = Integer.parseInt(line.substring(index1,index2));

		index1 = indexOfNextItem(line,index2);
		index2 = endOfItem(line,index1);
		int atom1 = getUsedAtomIndex(Integer.parseInt(line.substring(index1,index2)));

		index1 = indexOfNextItem(line,index2);
		index2 = endOfItem(line,index1);
		int atom2 = getUsedAtomIndex(Integer.parseInt(line.substring(index1,index2)));

		int stereo = 0;
		int topology = 0;
		while((index1 = indexOfNextItem(line,index2)) != -1){
			index2 = endOfItem(line,index1);
			String specifier = line.substring(index1,index2);
			int index = specifier.indexOf('=');
			String field = specifier.substring(0,index);
			int value = Integer.parseInt(specifier.substring(index + 1));
			if(field.equals("CFG")){
				switch(value){
					case 1: // up (version3)
						stereo = 1; // up (version2)
						break;
					case 2: // either (version3)
						stereo = (bondType == 2) ? 3 : 4; // 3:cross; 4:either (version2)
						break;
					case 3: // down (version3)
						stereo = 6; // down (version2)
						break;
				}
			} else if(field.equals("TOPO")){
				topology = value;
			} else{
				TRACE("Warning MolfileParser: Unused version 3 bond specifier:" + field + "\n");
			}
		}

		int bond = buildBond(atom1,atom2,bondType,stereo,topology);
		if(bond + 1 != bondIndex)
			mapBondIndex(bondIndex, bond);
	}

	private void interpretV3CollectionLine(String line)
	{
		String objectType = interpretObjectType(line);
		if (objectType != null) {
			int[] list = interpretV3List(line,objectType);
			if(line.startsWith("MDLV30/STEABS")){
				if (objectType.equals("ATOMS"))
					for(int i = 0;i < list.length;i++)
						mMol.setAtomESR(getUsedAtomIndex(list[i]),Molecule.cESRTypeAbs, -1);
				else
					for(int i = 0;i < list.length;i++)
						mMol.setBondESR(getUsedBondIndex(list[i]),Molecule.cESRTypeAbs, -1);
			} else if(line.startsWith("MDLV30/STERAC")){
				int group = Integer.parseInt(line.substring(13,indexOfWhiteSpace(line,13)));
				if (objectType.equals("ATOMS"))
					for(int i = 0;i < list.length;i++)
						mMol.setAtomESR(getUsedAtomIndex(list[i]),Molecule.cESRTypeAnd,group - 1);
				else
					for(int i = 0;i < list.length;i++)
						mMol.setBondESR(getUsedBondIndex(list[i]),Molecule.cESRTypeAnd,group - 1);
			} else if(line.startsWith("MDLV30/STEREL")){
				int group = Integer.parseInt(line.substring(13,indexOfWhiteSpace(line,13)));
				if (objectType.equals("ATOMS"))
					for(int i = 0;i < list.length;i++)
						mMol.setAtomESR(getUsedAtomIndex(list[i]),Molecule.cESRTypeOr,group - 1);
				else
					for(int i = 0;i < list.length;i++)
						mMol.setBondESR(getUsedBondIndex(list[i]),Molecule.cESRTypeOr,group - 1);
			} else if(line.startsWith("MDLV30/HILITE")){
				if (objectType.equals("ATOMS")){
					for(int i = 0;i < list.length;i++)
						mMol.setAtomColor(getUsedAtomIndex(list[i]),Molecule.cAtomColorDarkRed);
				} else {
					for(int i = 0;i < list.length;i++){
						int bond = getUsedBondIndex(list[i]);
						mMol.setAtomColor(mMol.getBondAtom(0, bond),Molecule.cAtomColorDarkRed);
						mMol.setAtomColor(mMol.getBondAtom(1, bond),Molecule.cAtomColorDarkRed);
					}
				}
			} else{
				TRACE("Error [readMoleculeFromBuffer]: Unknown version 3 collection type\n");
			}
		}
	}

	/**
	 * Interprets the object type of a collection and returns it as String
	 * @return object type or null if unsupported type
	 */
	private String interpretObjectType(String line) {
		if (line.contains("ATOMS=("))
			return "ATOMS";
		if (line.contains("BONDS=("))
			return "BONDS";

		TRACE("Error [readMoleculeFromBuffer]: Unknown or missing collection object type\n");
		return null;
	}

	/**
	 * Interprets the atom description line and returns the atom list for this atom
	 * @param line String Atom description line
	 * @return int[] Array contaning the atomic numbers for the list or null if no atom list could be interpreted
	 */
	private int[] interpretV3AtomList(String line)
	{
		int res[] = null;
//		System.out.println("Atom list |" + line + "|");
//		if(line.indexOf("NOT[") >= 0){
//			System.out.println("This is a 'NOT' list");
//		}
		int i1 = line.indexOf("[");
		int i2 = line.indexOf("]",i1);
		if(i1 >= 0 && i2 > 0){
			int atoms[] = new int[16];
			String s = line.substring(i1 + 1,i2);
			int index = 0;
			boolean ok = true;
			while(ok && index < 16){
				i1 = s.indexOf(",");
				String l = null;
				if(i1 == -1){
					l = s;
					ok = false;
				} else{
					l = s.substring(0,i1);
					s = s.substring(i1+1);
				}
				atoms[index++] = Molecule.getAtomicNoFromLabel(l, ALLOWED_ATOM_LABELS_IN_LIST);
			}
			res = new int[index];
			System.arraycopy(atoms,0,res,0,index);
		}
		return res;
	}

	/**
	 * Checks whether or not the atom description contains an atom list
	 * @param line String Atom description line
	 * @return int negative if an exclusion (NOT) list is present, positive if an atom list is present, 0 if no atom list. 
	 * The values for negative and positive results represent the index to the closing ']' bracket
	 */
	private int isV3AtomList(String line)
	{
		
		// simple check for atom list
		if (line.indexOf("[") >= 0) {
			// Detail check for non-quoted version
			int i1 = line.indexOf(" NOT[");
			int i2 = line.indexOf("]",i1);
			if(i1 >= 0 && i2 > 0){
				return -(i2+1); // point after the ]'
			} else{
				i1 = line.indexOf(" [");
				i2 = line.indexOf("]",i1);
				if(i1 >= 0 && i2 > 0){
					return i2+1; // point after the ]'
				}
			} 

			// Detail check for quoted version
			i1 = line.indexOf(" 'NOT[");
			i2 = line.indexOf("]'",i1);
			if(i1 >= 0 && i2 > 0){
				return -(i2+2); // point after the ]'
			} else{
				i1 = line.indexOf(" '[");
				i2 = line.indexOf("]'",i1);
				if(i1 >= 0 && i2 > 0){
					return i2+2; // point after the ]'
				}
			} 
			System.err.println("Warning invalid atom list in line: " + line);
		}
		return 0;
	}

	private int[] interpretV3List(String line,final String type)
	{
		int index1 = line.indexOf(type + "=(") + type.length() + 2;
		int index2 = line.indexOf(')',index1);
		int index = indexOfWhiteSpace(line,index1);
		int count = Integer.parseInt(line.substring(index1,index));
		int[] list = new int[count];
		for(int i = 0;i < count;i++){
			index1 = indexOfNextItem(line,index);
			index = indexOfWhiteSpace(line,index1);
			if(index == -1 || index > index2){
				index = index2;
			}
			list[i] = Integer.parseInt(line.substring(index1,index));
		}
		return list;
	}

	// with a given File, fill a Molecule
	public boolean parse(StereoMolecule mol, File file)
	{
		mMol = mol;
		try{
			BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file), StandardCharsets.UTF_8));
			BOMSkipper.skip(reader);
			return readMoleculeFromBuffer(reader);
		} catch(IOException e){
			System.err.println("Error reading file " + e);
		}
		return false;
	}

	// with a given String, fill a Molecule
	public boolean parse(StereoMolecule mol, String molFile)
	{
		return parse(mol,new BufferedReader(new StringReader(molFile)));
	}

	// with a given StringBuffer, fill a Molecule
	public boolean parse(StereoMolecule mol, StringBuffer molFile)
	{
		return parse(mol,molFile.toString());
	}

	public boolean parse(StereoMolecule m, BufferedReader rd)
	{
		mMol = m;
		return readMoleculeFromBuffer(rd);
	}

	// with a given String, create a compact sized Molecule
	public StereoMolecule getCompactMolecule(String molFile)
	{
		mMol = null;
		return readMoleculeFromBuffer(new BufferedReader(new StringReader(molFile))) ? mMol : null;
	}

	public StereoMolecule getCompactMolecule(BufferedReader reader)
	{
		mMol = null;
		return (readMoleculeFromBuffer(reader)) ? mMol : null;
	}

	public StereoMolecule getCompactMolecule(File file)
	{
		mMol = null;
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			boolean success = readMoleculeFromBuffer(reader);
			try { reader.close(); } catch (IOException ioe) {}
			return success ? mMol : null;
		} catch (FileNotFoundException fnfe) {
			return null;
		}
	}

	private int buildBond(int atom1,int atom2,int bondType,
						  int stereo,int topology)
	{
		int realBondType = Molecule.cBondTypeSingle;
		boolean isAtomESRAnd = false;

		switch(stereo){
			case 1:
				realBondType = Molecule.cBondTypeUp;
				break;
			case 3:
				realBondType = Molecule.cBondTypeCross;
				break;
			case 4: // we interpret 'either' as being racemic
				realBondType = Molecule.cBondTypeUp;
				isAtomESRAnd = true;
				break;
			case 6:
				realBondType = Molecule.cBondTypeDown;
				break;
			default:
				switch(bondType){
					case 1:
						realBondType = Molecule.cBondTypeSingle;
						break;
					case 2:
						realBondType = Molecule.cBondTypeDouble;
						break;
					case 3:
						realBondType = Molecule.cBondTypeTriple;
						break;
					case 4:
						realBondType = Molecule.cBondTypeDelocalized;
						break;
					case 9: // exists in version 3 only
						realBondType = Molecule.cBondTypeMetalLigand;
						break;
				}
				break;
		}

		int bond = mMol.addBond(atom1,atom2,realBondType);
		int queryFeatures = 0;

		if(isAtomESRAnd){
			mMol.setAtomESR(atom1,Molecule.cESRTypeAnd, -1);
		}

		if(bondType > 4){
			switch(bondType){
				case 5:
					queryFeatures |= Molecule.cBondTypeSingle | Molecule.cBondTypeDouble;
					break;
				case 6:
					queryFeatures |= Molecule.cBondTypeSingle | Molecule.cBondTypeDelocalized;
					break;
				case 7:
					queryFeatures |= Molecule.cBondTypeDouble | Molecule.cBondTypeDelocalized;
					break;
				case 8:
					if (realBondType != Molecule.cBondTypeMetalLigand)
						queryFeatures |= Molecule.cBondQFBondTypes;
					break;
			}
		}

		if(topology == 1){
			queryFeatures |= Molecule.cBondQFRing;
		}
		if(topology == 2){
			queryFeatures |= Molecule.cBondQFNotRing;
		}

		if(queryFeatures != 0){
			mMol.setBondQueryFeature(bond,queryFeatures,true);
		}

		return bond;
	}

	private void mapAtomIndex(int sourceAtomIndex, int usedAtomIndex) {
		if (mAtomIndexMap == null)
			mAtomIndexMap = new TreeMap<Integer,Integer>();

		mAtomIndexMap.put(new Integer(sourceAtomIndex), new Integer(usedAtomIndex));
	}

	private void mapBondIndex(int sourceBondIndex, int usedBondIndex) {
		if (mBondIndexMap == null)
			mBondIndexMap = new TreeMap<Integer,Integer>();

		mBondIndexMap.put(new Integer(sourceBondIndex), new Integer(usedBondIndex));
	}

	private int getUsedAtomIndex(int sourceAtomIndex) {
		Integer ui = (mAtomIndexMap == null) ? null : mAtomIndexMap.get(new Integer(sourceAtomIndex));
		return (ui == null) ? sourceAtomIndex-1 : ui.intValue();
	}

	private int getUsedBondIndex(int sourceBondIndex) {
		Integer ui = (mBondIndexMap == null) ? null : mBondIndexMap.get(new Integer(sourceBondIndex));
		return (ui == null) ? sourceBondIndex-1 : ui.intValue();
	}

	private int parseIntOrSpaces(String s) throws NumberFormatException
	{
		return(s.length() == 0) ? 0 : Integer.parseInt(s);
	}

	private int endOfItem(String line,int start)
	{
		int end = indexOfWhiteSpace(line,start + 1);
		return(end == -1) ? line.length() : end;
	}

	private int indexOfWhiteSpace(String line,int fromIndex)
	{
		for(int i = fromIndex;i < line.length();i++){
			if(line.charAt(i) == ' ' || line.charAt(i) == '\t'){
				return i;
			}
		}
		return -1;
	}

	private int indexOfNextItem(String line,int afterPreviousItem)
	{
		if(afterPreviousItem == -1){
			return -1;
		}
		for(int i = afterPreviousItem + 1;i < line.length();i++){
			if(line.charAt(i) != ' ' && line.charAt(i) != '\t'){
				return i;
			}
		}
		return -1;
	}

	void TRACE(String s)
	{
		if(debug){
			System.out.println(s);
		}
	}

	/**
	 * If we have single atoms from a metal to an electronegative atom
	 * that therefore exceeds its max valence, then reduce the bond to a
	 * metal ligand bond.
	 */
	private void introduceObviousMetalBonds() {
		int[] occupiedValence = new int[mMol.getAllAtoms()];

		// initialize with 1 for all delocalized atoms
		for (int bond=0; bond<mMol.getAllBonds(); bond++)
			if (mMol.getBondType(bond) == Molecule.cBondTypeDelocalized)
				for (int i=0; i<2; i++)
					occupiedValence[mMol.getBondAtom(i, bond)] = 1;

		// all bond orders
		for (int bond=0; bond<mMol.getAllBonds(); bond++) {
			int order = mMol.getBondOrder(bond);
			for (int i=0; i<2; i++)
				occupiedValence[mMol.getBondAtom(i, bond)] += order;
		}

		for (int bond=0; bond<mMol.getAllBonds(); bond++) {
			if (mMol.getBondOrder(bond) == 1) {
				for (int i=0; i<2; i++) {
					int metalAtom = mMol.getBondAtom(1-i, bond);
					if (mMol.isMetalAtom(metalAtom)) {
						int atom = mMol.getBondAtom(i, bond);
						if (mMol.isElectronegative(atom)
						 && occupiedValence[atom] > mMol.getMaxValence(atom)) {
							mMol.setBondType(bond, Molecule.cBondTypeMetalLigand);
							continue;
						}
					}
				}
			}
		}
	}

	/**
	 * SD-Files exported from the CSD database contain aromatic bonds rather than single/double bonds.
	 * Charges of aromatic systems are usually not given (e.g. in cyclopentadienyl(-) or pyridinium(+))
	 * and counter ions carry reduced charges to compensate (e.g. Fe in ferrocene wrongly has no charge assigned).
	 * To prevent valence problems and wrong idcode encoding we need to repair.
	 */
	private void deduceMissingCharges() {
		int[] chargeChange = new int[mMol.getAllAtoms()];
		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			chargeChange[atom] = -mMol.getAtomCharge(atom);

		new AromaticityResolver(mMol).locateDelocalizedDoubleBonds(null, true, false);

		for (int atom=0; atom<mMol.getAllAtoms(); atom++)
			chargeChange[atom] += mMol.getAtomCharge(atom);

		for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
			if (chargeChange[atom] != 0) {
				int chargeToDistribute = -chargeChange[atom];

				for (int bond=0; bond<mMol.getAllBonds(); bond++) {
					for (int i=0; i<2; i++) {
						if (chargeToDistribute > 0
						 && mMol.getBondType(bond) == Molecule.cBondTypeMetalLigand
						 && mMol.getBondAtom(1-i, bond) == atom) {
							int metal = mMol.getBondAtom(i, bond);
							if (mMol.isMetalAtom(metal)) {
								int maxCharge = getMaxOxidationState(metal);
								int charge = mMol.getAtomCharge(metal);
								if (charge < maxCharge) {
									int dif = Math.min(chargeToDistribute, maxCharge - charge);
									mMol.setAtomCharge(metal, charge + dif);
									chargeToDistribute -= dif;
								}
							}
						}
					}
				}

			}
		}
	}

	private int getMaxOxidationState(int metal) {
		int atomicNo = mMol.getAtomicNo(metal);
		byte[] os = (atomicNo < Molecule.cCommonOxidationState.length) ?
				Molecule.cCommonOxidationState[atomicNo] : null;
		return (os == null) ? 0 : os[os.length-1];
	}
}
