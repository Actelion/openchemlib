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

package com.actelion.research.chem.conf;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class BondLengthSet {
	private static final String cBondDataFile = "bondLengthData.txt";
	private static String sCustomBondDataFile = null;

	private static boolean isInitialized = false;
	private static int[] BOND_ID,BOND_COUNT;
	private static float[] BOND_LENGTH,BOND_STDDEV;

	public static final float DEFAULT_BOND_LENGTH = 2.0005f;
	public static final float DEFAULT_BOND_STDDEV = 1.0000f;

	private static final boolean CONSIDER_PI[] = { false,
				   false,  false,  false,  false,   true,   true,   //  H,  He, Li, Be, B,  C,
				   true,   true,  false,  false,  false,  false,    //  N,  O,  F,  Ne, Na, Mg,
				   false,  false,   true,   true };	                //  Al, Si, P,  S

	private final float[] mBondLength,mBondStdDev;

	/**
	 * Calculates and caches a list of bond length estimates from molecule.
	 * @param mol
	 */
	public BondLengthSet(final StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
	
		mBondLength = new float[mol.getAllBonds()];
		mBondStdDev = new float[mol.getAllBonds()];
		for (int bond=0; bond<mol.getAllBonds(); bond++) {
			int index = getBondIndex(mol, bond);
			if (index == -1) {
				mBondLength[bond] = getBondLengthFromCovalentRadii(mol, bond);
				mBondStdDev[bond] = getStdDevFromCovalentRadii(mol, bond);
				}
			else {
				mBondLength[bond] = getBondLength(index);
				mBondStdDev[bond] = getBondStdDev(index);
				}
			}
		}

	/**
	 * @param filePathAndName null (for default) or valid path & name to a custom bond length data file
	 */
	public static void setCustomDataFile(String filePathAndName) {
		sCustomBondDataFile = filePathAndName;
		isInitialized = false;
		}

	private static void initialize() {
		if (!isInitialized) {
			synchronized (BondLengthSet.class) {
				try {
					BufferedReader bdr = (sCustomBondDataFile == null) ? TorsionDB.openReader(cBondDataFile)
										: new BufferedReader(new FileReader(sCustomBondDataFile));

					String countString = bdr.readLine();
					int count = (countString == null) ? 0 : Integer.parseInt(countString);

					BOND_ID = new int[count];
					BOND_LENGTH = new float[count];
					BOND_STDDEV = new float[count];
					BOND_COUNT = new int[count];

					for (int i=0; i<count; i++) {
						String line = bdr.readLine();
						if (line != null) {
							String[] item = line.split("\\t");
							if (item.length == 4) {
								try {
									BOND_ID[i] = Integer.parseInt(item[0]);
									BOND_LENGTH[i] = Float.parseFloat(item[1]);
									BOND_STDDEV[i] = Float.parseFloat(item[2]);
									BOND_COUNT[i] = Integer.parseInt(item[3]);
									}
								catch (NumberFormatException nfe) {
									break;
									}
								}
							}
						}

					bdr.close();
					isInitialized = true;
					}
				catch (IOException e) {}
				}
			}
		}

	public float getLength(int bond) {
		return mBondLength[bond];
		}

	public float getStdDev(int bond) {
		return mBondStdDev[bond];
		}

	public static String getBondIDString(int index) {
		if (index == -1)
			return "unknown";
		int id = BOND_ID[index];
		int order = (id & 0xFF000000) >> 24;
		int pi1 = (id & 0x00F00000) >> 20;
		int pi2 = (id & 0x00000F00) >> 8;
		int atomicNo1 = (id & 0x000FF000) >> 12;
		int atomicNo2 = (id & 0x000000FF);
		String piString1 = (isPiConsidered(atomicNo1)) ? Integer.toString(pi1) : "";
		String piString2 = (isPiConsidered(atomicNo2)) ? Integer.toString(pi2) : "";
		String s = (order == 0) ? "d" : ((order > 3) ? "a" : "") + (order & 3);
		return s+Molecule.cAtomLabel[atomicNo1]+piString1+Molecule.cAtomLabel[atomicNo2]+piString2;
		}

	/**
	 * Constructs a bond classification ID from individual parameters and returns the ID's
	 * index in the sorted list of bond length information.
	 * The index can be used to get typical bond length and standard deviation.
	 * @param bondOrder
	 * @param isAromatic
	 * @param isDelocalized
	 * @param atomicNo1
	 * @param atomicNo2
	 * @param atomPi1
	 * @param atomPi2
	 * @return
	 */
	public static int getBondIndex(int bondOrder, boolean isAromatic, boolean isDelocalized, int atomicNo1, int atomicNo2, int atomPi1, int atomPi2) {
		return lookupBondIndex(getBondID(bondOrder, isAromatic, isDelocalized, atomicNo1, atomicNo2, atomPi1, atomPi2));
		}

	/**
	 * Constructs a bond classification index from individual parameters.
	 * @param bondOrder
	 * @param isAromatic
	 * @param isDelocalized
	 * @param atomicNo1
	 * @param atomicNo2
	 * @param atomPi1
	 * @param atomPi2
	 * @return
	 */
	private static int getBondID(int bondOrder, boolean isAromatic, boolean isDelocalized, int atomicNo1, int atomicNo2, int atomPi1, int atomPi2) {
		int pi1 = (atomicNo1 < CONSIDER_PI.length && CONSIDER_PI[atomicNo1]) ? atomPi1 << 8 : 0;
		int pi2 = (atomicNo2 < CONSIDER_PI.length && CONSIDER_PI[atomicNo2]) ? atomPi2 << 8 : 0;
		int atomType1 = pi1 + atomicNo1;
		int atomType2 = pi2 + atomicNo2;
		int bondType = isDelocalized ? 0 : isAromatic ? 4+bondOrder : bondOrder;
		return (bondType<<24)+((atomType1<atomType2)?(atomType1<<12)+atomType2:(atomType2<<12)+atomType1);
		}

	/**
	 * Constructs a bond classification ID from a specific bond in a molecule and returns the ID's
	 * index in the sorted list of bond length information.
	 * The index can be used to get typical bond length and standard deviation.
	 * @param mol
	 * @param bond
	 * @return
	 */
	public static int getBondIndex(StereoMolecule mol, int bond) {
		int atom1 = mol.getBondAtom(0, bond);
		int atom2 = mol.getBondAtom(1, bond);
		int atomicNo1 = mol.getAtomicNo(atom1);
		int atomicNo2 = mol.getAtomicNo(atom2);
		return getBondIndex(mol.getBondOrder(bond), mol.isAromaticBond(bond), mol.isDelocalizedBond(bond), atomicNo1, atomicNo2, getAtomPi(mol, atom1), getAtomPi(mol, atom2));
		}

	/**
	 * Constructs a bond classification index from a specific bond in a molecule.
	 * The index can be used to get typical bond length and standard deviation.
	 * @param mol
	 * @param bond
	 * @return
	 */
	public static int getBondID(StereoMolecule mol, int bond) {
		int atom1 = mol.getBondAtom(0, bond);
		int atom2 = mol.getBondAtom(1, bond);
		int atomicNo1 = mol.getAtomicNo(atom1);
		int atomicNo2 = mol.getAtomicNo(atom2);
		return getBondID(mol.getBondOrder(bond), mol.isAromaticBond(bond), mol.isDelocalizedBond(bond), atomicNo1, atomicNo2, getAtomPi(mol, atom1), getAtomPi(mol, atom2));
		}

	/**
	 * @param atomicNo
	 * @return whether this atomicNo uses the atom's pi count to further distinguish bond type cases
	 */
	public static final boolean isPiConsidered(int atomicNo) {
		return (atomicNo < CONSIDER_PI.length) && CONSIDER_PI[atomicNo];
		}

	/**
	 * Returns an atom's pi electron count for the purpose of classifying a connected bond
	 * to determine its length. The returned value may differ from the formal pi-electron count
	 * if the formal count does not represent the real mesomeric situation well.
	 * @param mol
	 * @param atom
	 * @return pi electron count
	 */
	public static int getAtomPi(StereoMolecule mol, int atom) {
		if (atom >= mol.getAtoms())
			return 0;

		if (mol.isAromaticAtom(atom) && mol.getAtomicNo(atom) == 6 && mol.getAtomCharge(atom) != 0)
			return 1;

		return mol.getAtomPi(atom);
		}

	/**
	 * Returns an estimate of the bond length based on atom and bond characteristics.
	 * The bond is classified based on its characteristics the returned value is
	 * the median of equally classified bonds within the COD database. Statistics of
	 * purely organic bonds (no metal atoms) are taken from the organic subset of the COD.
	 * Requires cHelperRings level of helper arrays.
	 * @param mol
	 * @param bond
	 * @return
	 */
	public static float lookupBondLength(StereoMolecule mol, int bond) {
		int index = getBondIndex(mol, bond);
		return (index == -1) ? getBondLengthFromCovalentRadii(mol, bond) : getBondLength(index);
		}

	/**
	 * Returns an estimate of the bond length based on atom and bond characteristics.
	 * The bond is classified based on its characteristics the returned value is
	 * the median of equally classified bonds within the COD or CSD database. Statistics of
	 * purely organic bonds (no metal atoms) are taken from the organic subset of the COD/CSD.
	 * If there are no similar bond in the crystallographic data sets, the the bond length
	 * is estimated from the covalent radii, which if not available may be estimated from
	 * the van-der-waals radii.
	 * @param mol
	 * @param bond
	 * @return
	 */
	public static float getBondLengthFromCovalentRadii(StereoMolecule mol, int bond) {
		int atomicNo1 = mol.getAtomicNo(mol.getBondAtom(0, bond));
		int atomicNo2 = mol.getAtomicNo(mol.getBondAtom(1, bond));
		return getCovalentRadius(atomicNo1) + getCovalentRadius(atomicNo2);
		}

	/**
	 * Returns an the standard deviation of bond lengths from bonds with similar
	 * characteristics from crystallographic data. If the estimate is based on
	 * covalent radii or van-der-waals radii, the standard deviation reflects that.
	 * @param mol
	 * @param bond
	 * @return
	 */
	public static float getStdDevFromCovalentRadii(StereoMolecule mol, int bond) {
		int atomicNo1 = mol.getAtomicNo(mol.getBondAtom(0, bond));
		int atomicNo2 = mol.getAtomicNo(mol.getBondAtom(1, bond));
		return (atomicNo1 < VDWRadii.COVALENT_RADIUS.length ? 0.05f : 0.125f)
			 + (atomicNo2 < VDWRadii.COVALENT_RADIUS.length ? 0.05f : 0.125f);
		}

	private static float getCovalentRadius(int atomicNo) {
		return (atomicNo < VDWRadii.COVALENT_RADIUS.length) ? VDWRadii.COVALENT_RADIUS[atomicNo]
			 : (atomicNo < VDWRadii.VDW_RADIUS.length) ? 0.6f * VDWRadii.VDW_RADIUS[atomicNo] : 1.8f;
	}

	/**
	 * Returns an estimate of the bond length based on atom and bond characteristics.
	 * The bond is classified based on its characteristics the returned value is
	 * the median of equally classified bonds within the COD or CSD database. Statistics of
	 * purely organic bonds (no metal atoms) are taken from the organic subset of the COD/CSD.
	 * @param index bond id index obtained with getBondIndex()
	 * @return
	 */
	public static float getBondLength(int index) {
		return (index == -1) ? DEFAULT_BOND_LENGTH : BOND_LENGTH[index];
		}

	public static float getBondCount(int index) {
		return (index == -1) ? 0 : BOND_COUNT[index];
	}

	/**
	 * Returns an the standard deviation of bond lengths from bonds with similar
	 * characteristics from crystallographic data.
	 * @param index bond id index obtained with getBondIndex()
	 * @return
	 */
	public static float getBondStdDev(int index) {
		return (index == -1) ? DEFAULT_BOND_STDDEV : BOND_STDDEV[index];
		}


	/**
	 * Returns an estimate of the bond length based on atom and bond characteristics.
	 * Requires cHelperRings level of helper arrays.
	 * If no statistics information is available, then it return DEFAULT_BOND_LENGTH,
	 * which is garanteed to be different from any value in the statistics table.
	 * @param id valid bond classification obtained with one of the getBondType() methods
	 * @return DEFAULT_BOND_LENGTH if no statistics information is available for that bond type
	 *
	private static float lookupBondLength(int id) {
		int index = lookupBondIndex(id);
		return (index == -1) ? DEFAULT_BOND_LENGTH : BOND_LENGTH[index];
		}*/

	private static int lookupBondIndex(int id) {
		if (!isInitialized)
			initialize();

		int index = 2048;
		int increment = 1024;
		for (int i=0; i<12; i++) {
			int comparison = (index >= BOND_ID.length || id < BOND_ID[index]) ? -1
					: (id == BOND_ID[index]) ? 0 : 1;
			if (comparison == 0)
				return index;
			index = (comparison < 0) ? index-increment : index+increment;
			increment /= 2;
			}

		return -1;	// if not found
		}
	}
