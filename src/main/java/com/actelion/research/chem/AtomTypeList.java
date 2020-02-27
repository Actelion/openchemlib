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

package com.actelion.research.chem;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.TreeMap;
import java.util.TreeSet;

import com.actelion.research.chem.io.CompoundFileParser;
import com.actelion.research.chem.io.DWARFileParser;
import com.actelion.research.chem.io.SDFileParser;


public class AtomTypeList {
    private static final String VERSION_STRING = "AtomTypeList v1.1";

    private static final int DEFAULT_MIN_ATOMS = 10;
	private static final int DEFAULT_MAX_ATOMS = 50;

    private TreeMap<Long,Integer>	mCountList;
	private TreeMap<Long,Double>	mProbabilityList;
	private float[]	            	mRingSizeAdjust;
	private int                 	mAtomTypeMode,mMaxAtoms,mMinAtoms;

	/**
	 * Creates an empty AtomTypeList, which must be populated by multiply calling processMolecule()
	 * and finally calling finalizeProcessMolecules() once.
	 * @param mode
	 */
    public AtomTypeList(int mode) {
        mRingSizeAdjust = new float[8];
        mCountList = new TreeMap<Long,Integer>();
		mAtomTypeMode = mode;
		mMinAtoms = DEFAULT_MIN_ATOMS;
		mMaxAtoms = DEFAULT_MAX_ATOMS;
        }


	/**
	 * Creates a new AtomTypeList from a given file using the given mode.
	 * If the the filename references a .typ file, then the mode is checked, whether it matches the file's content.
	 * If the the filename references a compound file, then the molecules are parsed and a new AtomTypeList is created
	 * reflecting the all contained atom types.
	 * @param filename either .typ file or a .dwar or .sdf compound file
	 * @param mode
	 * @throws Exception
	 */
	public AtomTypeList(String filename, int mode) throws Exception {
        this(mode);

        if (filename.endsWith(".typ")) {
	        BufferedReader theReader = new BufferedReader(new InputStreamReader(this.getClass().getResourceAsStream(filename)));
	        String version =theReader.readLine();
	        if (!VERSION_STRING.equals(version)) {
	            throw new Exception("Outdated atom type list file.");
	            }

	        mAtomTypeMode = Integer.parseInt(theReader.readLine());
	        if (mAtomTypeMode != mode) {
	            throw new Exception("Incompatible atom type mode.");
	            }

	        for (int i=0; i<8; i++)
	            mRingSizeAdjust[i] = Float.parseFloat(theReader.readLine());

	        while (true) {
	            String theLine = theReader.readLine();
	            if (theLine == null)
	                break;

	            int tab = theLine.indexOf('\t');
	            mCountList.put(new Long(Long.parseLong(theLine.substring(0, tab))),
	                      new Integer(Integer.parseInt(theLine.substring(tab+1))));
	            }
	        theReader.close();
	        return;
	        }

        CompoundFileParser parser = filename.endsWith(".dwar") ? new DWARFileParser(filename)
								  : filename.endsWith(".sdf")  ? new SDFileParser(filename) : null;

		if (parser != null) {
			TreeSet<String> moleculeCache = new TreeSet<String>();
			try {
				while (parser.next())
					processMolecule(parser.getMolecule(), moleculeCache);
				}
			catch (Exception e) {
				e.printStackTrace();
				}

			finalizeProcessMolecules();
			}
		}


	public void finalizeProcessMolecules() {
		float ringSum = 0.0f;
		for (int i=0; i<8; i++)
			ringSum += mRingSizeAdjust[i];
		if (ringSum != 0)
			for (int i=0; i<8; i++)
				mRingSizeAdjust[i] /= ringSum;
		}

	public synchronized void calculateProbabilities() {
		if (mProbabilityList == null) {
			int averageCount = 0;
			for (Integer count: mCountList.values())
				averageCount += count;
			averageCount /= mCountList.size();

			mProbabilityList = new TreeMap<>();
			for (Long type: mCountList.keySet())
				mProbabilityList.put(type, (double)mCountList.get(type)/averageCount);
			}
		}


	public void writeTypeFile(String filename) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
			writer.write(VERSION_STRING);
			writer.newLine();

			writer.write(""+mAtomTypeMode);
			writer.newLine();

			for (int i=0; i<8; i++) {
				writer.write(""+mRingSizeAdjust[i]);
				writer.newLine();
				}

			for (Long type: mCountList.keySet()) {
				writer.write(type.toString()+"\t"+ mCountList.get(type).toString());
				writer.newLine();
				}
			writer.close();
			}
		catch (Exception e) {
			e.printStackTrace();
			}
		}


	/**
	 * Writes this AtomTypeList into a TAB-delimited text file in human readable form.
	 * Atom type frequency and all properties making up an atom type are written into separate columns.
	 * @param textfilename
	 * @param mode
	 */
	public void writeTextFile(String textfilename, int mode) {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(textfilename));
			writer.write("AtomType\tFrequency\t"+AtomTypeCalculator.getHeaderString(mode));
			writer.newLine();

			for (Long type: mCountList.keySet()) {
				writer.write(Long.toString(type)+"\t"+ mCountList.get(type).toString()+"\t"+AtomTypeCalculator.getTypeString(type, mode));
				writer.newLine();
				}
			writer.close();
			}
		catch (Exception e) {
			e.printStackTrace();
			}
		}


	public TreeMap<Long,Integer> getAtomTypeList() {
	    return mCountList;
	    }


	/**
	 *
	 * @param mol
	 * @param moleculeCache
	 */
	public void processMolecule(StereoMolecule mol, TreeSet<String> moleculeCache) {
		if (mol != null) {
			mol.stripIsotopInfo();
			mol.stripSmallFragments();
			mol.stripStereoInformation();

			mol.ensureHelperArrays(Molecule.cHelperNeighbours);

			boolean containsMetal = false;
			for (int atom=0; atom<mol.getAtoms(); atom++) {
				if (mol.isMetalAtom(atom)) {
					containsMetal = true;
					break;
					}
				}

			if (!containsMetal && mol.getAtoms() >= mMinAtoms && mol.getAtoms() <= mMaxAtoms) {
				String idcode = new SimpleCanonizer(mol).getIDCode();
				if (!moleculeCache.contains(idcode)) {
					moleculeCache.add(idcode);

					for(int atom=0;atom<mol.getAtoms();atom++)
						try {
							addAtomType(AtomTypeCalculator.getAtomType(mol, atom, mAtomTypeMode));
							} catch (Exception e) {}

//for(int atom=0;atom<mol.getAtoms();atom++) { long type = AtomTypeCalculator.getAtomType(mol, atom, mAtomTypeMode); if (type == 2689028 || type == 11077636) System.out.println(idcode); };

					RingCollection ringSet = mol.getRingSet();
					for (int ring=0; ring<ringSet.getSize(); ring++)
						mRingSizeAdjust[ringSet.getRingSize(ring)]++;
					}
				}
			}
		}


	private void addAtomType(Long atomType) {
	    Integer count = mCountList.get(atomType);
	    if (count == null)
	        mCountList.put(atomType, 1);
	    else
            mCountList.put(atomType, count+1);
		}


	public int getCountFromType(long type){
	    return mCountList.get(new Long(type)).intValue();
		}


	public double getProbabilityFromType(long type){
		return mProbabilityList.get(type);
	}


	public float getRingSizeAdjust(int ringSize) {
		return mRingSizeAdjust[ringSize];
		}
	}
