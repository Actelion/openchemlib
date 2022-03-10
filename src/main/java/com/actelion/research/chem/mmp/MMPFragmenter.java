/*
 * Copyright (c) 2017
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
 * 3. Neither the name of the copyright holder nor the
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
 * @author Gregori Gerebtzoff
 */

package com.actelion.research.chem.mmp;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.StereoMolecule;

import java.util.*;
import java.util.Map.Entry;

public class MMPFragmenter {
	private static final int SINGLE_CUT = 1;
	private static final int DOUBLE_CUT = 2;
	private static final int NUMBER_OF_CUTS = DOUBLE_CUT;      // desired number of cuts; currently supported: SINGLE_CUT, DOUBLE_CUT
	public static final Integer KEYS_MIN_ATOMS = 4;            // set to null to bypass
	private static final Integer VALUE_MAX_ATOMS = null;       // set to null to bypass (default:15)
	private static final String r1H = createR1HMoleculeID();   // R1-H molecule IDCode (for Hydrogen replacements)
	public static final String FRAGMENT_DELIMITER = "#";       // Symbol used for fragment delimiter
	public static final int FRAGMENT_ATOMIC_NO = 142;          // Atomic number used for fragment delimiter
	
	private StereoMolecule mol;                                // molecule might get modified upon initialization (Hydrogens are removed)
	private List<String[]> moleculeFragmentsID;                // {R-group, "clean" (without R-groups)} fragmentsID for Hydrogen replacements
	private HashMap<String, StereoMolecule> moleculeFragments; // unique list of R-groups- containing fragmentsID for Hydrogen replacements (IDCode, StereoMol)
	private List<MoleculeIndexID> moleculeIndexesID;           // container of keysID-valueID for single cuts and double cuts
	private List<MoleculeIndexIDByte> moleculeIndexesIDByte;   // container of keysIDBytes-valueIDBytes for single cuts and double cuts
	private Integer nRotBonds = null;
	private ArrayList<Integer> rotBondsIndex;
	
	/**
	 * Generates the idCode for the R-[H] molecule
	 * @return idCode String
	 */
	static public String createR1HMoleculeID() {
		StereoMolecule r1H = new StereoMolecule();
		int atom1Index = r1H.addAtom(FRAGMENT_ATOMIC_NO);
		r1H.setAtomCustomLabel(atom1Index, FRAGMENT_DELIMITER);
		int atom2Index = r1H.addAtom(1);
		r1H.setAtomCustomLabel(atom2Index, "[H]");
		r1H.addBond(atom1Index, atom2Index, StereoMolecule.cBondTypeSingle);
		return getIDCodeWithCustomLabels(r1H);
	}
	
	/**
	 * Helper function to generate idCodes with custom labels
	 * @param mol
	 * @return idCode String
	 */
	static private String getIDCodeWithCustomLabels(StereoMolecule mol) {
		Canonizer canonizer = new Canonizer(mol, Canonizer.ENCODE_ATOM_CUSTOM_LABELS);
		return canonizer.getIDCode();
	}
	
	static public class MoleculeIndexID {
		private String[] keysID;
		private int[] keysIndex;
		private String valueID;
		private int valueIndex;
		private int[] keysIDAtoms;
		private int valueIDAtoms;
		private IDCodeParser idCodeParser = new IDCodeParser();
		private int[] bondIndexes;                                // bond indexes for bonds between two heavy atoms
		private int[] valueAtomIndexes;                           // atom indexes of the key for heavy atom - hydrogen bonds
		private List<Double[]> coordinates;                       // coordinates of the middle of the bond between key and value
		private Integer[] chemicalSpaceSizes;
		
		public MoleculeIndexID(){ }
		
		/**
		 * Adds one keys-value combination
		 * @param keysID Array of one (single cut) or two (double cut) idCodes of the 'key' (constant part of the molecule)
		 * @param keysIndex Array of one (single cut) or two (double cut) indexes of the 'key' (from mmpUniqueFragments)
		 * @param valueID idCode of the 'value' (variable part of the molecule)
		 * @param valueIndex Index of the 'value'
		 * @param keysIDAtoms Number of heavy atoms of the 'key(s)'
		 * @param valueIDAtoms Number of heavy atoms of the 'value'
		 * @param bondIndexes Array of one (single cut) or two (double cut) bond indexes where the cuts occur
		 * @param valueAtomIndexes Atom indexes of the 'key' for heavy atom - hydrogen bonds
		 */
		public MoleculeIndexID(String[] keysID, int[] keysIndex, String valueID, int valueIndex, int[] keysIDAtoms, int valueIDAtoms, int[] bondIndexes, int[] valueAtomIndexes){
			this.keysID = keysID;
			this.keysIndex = keysIndex;
			this.valueID = valueID;
			this.valueIndex = valueIndex;
			this.keysIDAtoms = keysIDAtoms;
			if (keysIDAtoms == null) {
				this.keysIDAtoms = new int[keysID.length];
				for (int i=0; i<keysID.length; i++) {
					this.keysIDAtoms[i] = idCodeParser.getAtomCount(keysID[i]) - 1;
				}
			}
			this.valueIDAtoms = valueIDAtoms;
			this.bondIndexes = bondIndexes;
			this.valueAtomIndexes = valueAtomIndexes;
			this.coordinates = new ArrayList<Double[]>();
		}
		
		/**
		 * Adds one keys-value combination
		 * @param keysID Array of one (single cut) or two (double cut) idCodes of the 'key' (constant part of the molecule)
		 * @param valueID idCode of the 'value' (variable part of the molecule)
		 * @param keysIDAtoms Number of heavy atoms of the 'key(s)'
		 * @param valueIDAtoms Number of heavy atoms of the 'value'
		 * @param bondIndexes Array of one (single cut) or two (double cut) bond indexes where the cuts occur
		 * @param valueAtomIndexes Atom indexes of the 'key' for heavy atom - hydrogen bonds
		 */
		public MoleculeIndexID(String[] keysID, String valueID, int[] keysIDAtoms, int valueIDAtoms, int[] bondIndexes, int[] valueAtomIndexes){
			this.keysID = keysID;
			this.valueID = valueID;
			this.keysIDAtoms = keysIDAtoms;    // number of atoms minus one (because of R1)
			if (keysIDAtoms == null) {
				this.keysIDAtoms = new int[keysID.length];
				for (int i=0; i<keysID.length; i++) {
					this.keysIDAtoms[i] = idCodeParser.getAtomCount(keysID[i]) - 1;
				}
			}
			this.valueIDAtoms = valueIDAtoms;  // number of atoms minus one (because of R1)
			this.bondIndexes = bondIndexes;
			this.valueAtomIndexes = valueAtomIndexes;
			this.coordinates = new ArrayList<Double[]>();
		}
		
		public String[] getKeysID() {
			return keysID;
		}
		
		public String getValueID() {
			return valueID;
		}
		
		public int[] getKeysIDAtoms() {
			return keysIDAtoms;
		}
		
		public int getValueIDAtoms() {
			return valueIDAtoms;
		}
		
		public int[] getKeysIndex() {
			return keysIndex;
		}
		
		public void setKeysIndex(int[] keysIndex) {
			this.keysIndex = keysIndex;
		}
		
		public int getValueIndex() {
			return valueIndex;
		}
		
		public void setValueIndex(int valueIndex) {
			this.valueIndex = valueIndex;
		}
		
		public int[] getBondIndexes() {
			return bondIndexes;
		}
		
		public int[] getValueAtomIndexes() {
			return valueAtomIndexes;
		}
		
		public void setCoordinates(double x, double y) {
			this.coordinates.add(new Double[]{x, y});
		}
		
		public List<Double[]> getCoordinates() {
			return coordinates;
		}
		
		public void setChemicalSpaceSize(Integer[] chemicalSpaceSizes) {
			this.chemicalSpaceSizes = chemicalSpaceSizes; 
		}
		
		public Integer[] getChemicalSpaceSizes() {
			return chemicalSpaceSizes;
		}
	}
	
	public class MoleculeIndexIDByte {
		private byte[][] keysIDByte;
		private byte[] valueIDByte;
		
		public MoleculeIndexIDByte(){ }
		
		public MoleculeIndexIDByte(byte[][] keysIDByte, byte[] valueIDByte){
			this.keysIDByte = keysIDByte;
			this.valueIDByte = valueIDByte;
		}
		
		public byte[][] getKeysIDByte() {
			return keysIDByte;
		}
		
		public byte[] getValueIDByte() {
			return valueIDByte;
		}
	}

	public MMPFragmenter(StereoMolecule mol) {
		this.mol = removeHydrogens(mol);
		this.moleculeFragmentsID = new ArrayList<String[]>();
		this.moleculeFragments = new HashMap<String, StereoMolecule>();
		this.moleculeIndexesID = new ArrayList<MoleculeIndexID>();
		this.moleculeIndexesIDByte = new ArrayList<MoleculeIndexIDByte>();
	}
	
	public List<MoleculeIndexID> getMoleculeIndexesID() {
		return getMoleculeIndexesID(true);
	}

	/**
	 * Returns an ArrayList of MoleculeIndexID containing<br>
	 * IDCodes of keys and values of single and double cuts.
	 * @param generateWholeMoleculeVariations true/false to generate whole molecule variations, used for identification of H-replacements
	 */
	public List<MoleculeIndexID> getMoleculeIndexesID(boolean generateWholeMoleculeVariations) {
		if (nRotBonds == null) {
			fragmentMolecule(generateWholeMoleculeVariations);
		}
		return moleculeIndexesID;
	}
	
	/**
	 * Returns an ArrayList of {R-group, "clean" (without R-groups)}<br>
	 * fragmentsID used for Hydrogen replacements.
	 */
	public List<String[]> getMoleculeFragmentsID() {
		if (nRotBonds == null) {
			fragmentMolecule(false);
		}
		return moleculeFragmentsID;
	}
	
	public List<MoleculeIndexIDByte> getMoleculeIndexesIDByte() {
		return getMoleculeIndexesIDByte(true);
	}
	
	/**
	 * Returns an ArrayList of MoleculeIndexIDByte containing
	 * bytes of IDCodes of keys and values of single and double cuts.
	 * @param generateWholeMoleculeVariations true/false to generate whole molecule variations, used for identification of H-replacements
	 */
	public List<MoleculeIndexIDByte> getMoleculeIndexesIDByte(boolean generateWholeMoleculeVariations) {
		if (nRotBonds == null) {
			fragmentMolecule(generateWholeMoleculeVariations);
		}
		if (nRotBonds > 0 && moleculeIndexesIDByte.size() == 0) {
			for (MoleculeIndexID moleculeIndexID:moleculeIndexesID) {
				String[] keys = moleculeIndexID.keysID;
				String value = moleculeIndexID.valueID;
				byte[][] keysByte = null;
				if (keys.length == SINGLE_CUT) {
					keysByte = new byte[][]{keys[0].getBytes(), keys[1].getBytes()};
				}
				else if (keys.length == DOUBLE_CUT) {
					keysByte = new byte[][]{keys[0].getBytes(), keys[1].getBytes(), keys[2].getBytes()};
				} 
				byte[] valueByte = value.getBytes();
				MoleculeIndexIDByte moleculeIndexByte = new MoleculeIndexIDByte(keysByte, valueByte);
				moleculeIndexesIDByte.add(moleculeIndexByte);
			}
		}
		return moleculeIndexesIDByte;
	}
	
	/**
	 * Returns the indexes of the rotatable bonds
	 * @return Array of integers
	 */
	public ArrayList<Integer> getRotBondsIndex() {
		if (rotBondsIndex == null) {
			rotBondsIndex = new ArrayList<Integer>();
			for (int bond=0; bond<mol.getBonds(); bond++) {
				if (!mol.isRingBond(bond) && mol.getBondOrder(bond) == 1) {
					rotBondsIndex.add(bond);
				}
			}
		}
		return rotBondsIndex;
	}
	
	/**
	 * Counts the number of R-Groups
	 * @param mol
	 * @return number of R-Groups
	 */
	private int countRGroups(StereoMolecule mol) {
		int count = 0;
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mol.getAtomicNo(atom) == 0 || mol.getAtomicNo(atom) >= 142) {
				count++;
			}
		}
		return count;
	}
	
	private StereoMolecule addRGroups(StereoMolecule mol) {
		return addRGroups(mol, false);
	}
	
	/**
	 * Tags R-Groups with correct label (#1, #2)
	 * @param mol
	 * @param inverse true/false if the middle fragment is inverted
	 * @return Modified StereoMolecule
	 */
	private StereoMolecule addRGroups(StereoMolecule mol, boolean inverse) {
		int rGroup = 1; // atomic number of R1
		mol.ensureHelperArrays(StereoMolecule.cHelperNeighbours);
		if (inverse == false) {
			for (int atom=0; atom<mol.getAtoms(); atom++) {
				if (mol.getAtomicNo(atom) == 0 || mol.getAtomicNo(atom) >= 142) {
					mol.setAtomicNo(atom, FRAGMENT_ATOMIC_NO + rGroup - 1);
					mol.setAtomCustomLabel(atom, FRAGMENT_DELIMITER + Integer.toString(rGroup));
					rGroup++;
				}
			}
		}
		else {
			for (int atom=mol.getAtoms()-1; atom>=0; atom--) {
				if (mol.getAtomicNo(atom) == 0 || mol.getAtomicNo(atom) >= 142) {
					mol.setAtomicNo(atom, FRAGMENT_ATOMIC_NO + rGroup - 1);
					mol.setAtomCustomLabel(atom, FRAGMENT_DELIMITER + Integer.toString(rGroup));
					rGroup++;
				}
			}
		}
		return mol;
	}
	
	/**
	 * Creates a new MoleculeIndexID from fragments
	 * @param fragments Array of two (single cut) or three (double cut) fragments
	 * @param bondIndexes Array of one or two bond indexes
	 * @param valueAtomIndexes Array of atom indexes for Hydrogen replacements 
	 * @param cutType SINGLE_CUT or DOUBLE_CUT
	 * @return a MoleculeIndexID
	 */
	private MoleculeIndexID processFragments(StereoMolecule[] fragments, int[] bondIndexes, int[] valueAtomIndexes, int cutType) {
		MoleculeIndexID retVal = null;
		if (cutType == SINGLE_CUT && fragments.length > 1) {
			String idCode1 = getIDCodeWithCustomLabels(fragments[0]);
			String idCode2 = getIDCodeWithCustomLabels(fragments[1]);
			retVal = new MoleculeIndexID(new String[]{idCode1}, idCode2, new int[]{fragments[0].getAtoms()-1}, fragments[1].getAtoms()-1, bondIndexes, valueAtomIndexes);
			if (KEYS_MIN_ATOMS == null || fragments[0].getAtoms() >= KEYS_MIN_ATOMS) {
				moleculeFragments.put(idCode1, fragments[0]);
			}
			if (KEYS_MIN_ATOMS == null || fragments[1].getAtoms() >= KEYS_MIN_ATOMS) {
				moleculeFragments.put(idCode2, fragments[1]);
			}
		}
		else if (cutType == DOUBLE_CUT && fragments.length > 2) {
			for (StereoMolecule mol:fragments) {
				mol.ensureHelperArrays(StereoMolecule.cHelperNeighbours);
			}
			if (countRGroups(fragments[1]) == 2) {
				fragments = new StereoMolecule[]{fragments[1], fragments[0], fragments[2]};
			}
			else if (countRGroups(fragments[2]) == 2) {
				fragments = new StereoMolecule[]{fragments[2], fragments[0], fragments[1]};
			}
			Canonizer canonizer = new Canonizer(fragments[0]);
			int[] graphIndexes = canonizer.getGraphIndexes();
			int[] rGroupIndexes = new int[2];
			int rGroupCounter = 0;
			for (int atom=0; atom<fragments[0].getAtoms(); atom++) {
				if (fragments[0].getAtomicNo(atom) >= 142 || fragments[0].getAtomicNo(atom) == 0) {
					rGroupIndexes[rGroupCounter] = atom;
					rGroupCounter++;
				}
			}
			if (graphIndexes[rGroupIndexes[0]] < graphIndexes[rGroupIndexes[1]]) {
				retVal = new MoleculeIndexID(new String[]{getIDCodeWithCustomLabels(fragments[1]), getIDCodeWithCustomLabels(fragments[2])}, getIDCodeWithCustomLabels(addRGroups(fragments[0])), new int[]{fragments[1].getAtoms()-1, fragments[2].getAtoms()-1}, fragments[0].getAtoms()-1, bondIndexes, valueAtomIndexes);
			}
			else {
				retVal = new MoleculeIndexID( new String[]{getIDCodeWithCustomLabels(fragments[2]), getIDCodeWithCustomLabels(fragments[1])}, getIDCodeWithCustomLabels(addRGroups(fragments[0], true)), new int[]{fragments[2].getAtoms()-1, fragments[1].getAtoms()-1}, fragments[0].getAtoms()-1, new int[]{bondIndexes[1], bondIndexes[0]}, new int[]{valueAtomIndexes[1], valueAtomIndexes[0]});
			}
		}
		return retVal;
	}
	
	/**
	 * Removes hydrogens from a StereoMolecule (in case the source is a SDF file with explicit hydrogens)
	 * @param mol
	 * @return a modified StereoMolecule
	 */
	private StereoMolecule removeHydrogens(StereoMolecule mol) {
		mol.ensureHelperArrays(StereoMolecule.cHelperNeighbours | StereoMolecule.cHelperParities);
		mol.setAllAtoms(mol.getAtoms()); // this way the hydrogens will just be ignored and overwritten by new atoms later
		mol.setAllBonds(mol.getBonds());
		return mol;
	}
	
	public void fragmentMolecule() {
		fragmentMolecule(true);
	}
	
	/**
	 * Fragments one StereoMolecule into keys-value pairs of StereoMolecules;<br>
	 *  - for single cuts, keys is an array containing one fragment<br>
	 *    ('constant' part of the molecule), value the 'variable' part;<br>
	 *  - for double cuts, keys in an array containing two fragments<br> 
	 *    (constant parts of the molecule, i.e. 'left' and 'right' part), value the 'middle' part.<br>
	 * An array containing all 'single cut' fragments is available for Hydrogen replacements.
	 * @param generateWholeMoleculeVariations true/false to generate whole molecule variations, used for identification of H-replacements
	 */
	public void fragmentMolecule(boolean generateWholeMoleculeVariations) {
		int[] rGroupsIndex = new int[NUMBER_OF_CUTS*2];
		int[] newBondsIndex = new int[NUMBER_OF_CUTS];
		Set<String> hFragments = new HashSet<String>();
		StereoMolecule editableMol = new StereoMolecule();
		mol.copyMolecule(editableMol);
		editableMol.ensureHelperArrays(StereoMolecule.cHelperNeighbours);
		// We add R1-R1 and R2-R2 fragments to the molecule;
		// These fragments will be used to decompose the molecule and allow recycling of the molecule
		for (int i=0; i<NUMBER_OF_CUTS; i++) {
			rGroupsIndex[i*2] = editableMol.addAtom(FRAGMENT_ATOMIC_NO); // R1a and R2a
			editableMol.setAtomCustomLabel(rGroupsIndex[i*2], FRAGMENT_DELIMITER);
			rGroupsIndex[i*2+1] = editableMol.addAtom(FRAGMENT_ATOMIC_NO); // R1b and R2b
			editableMol.setAtomCustomLabel(rGroupsIndex[i*2+1], FRAGMENT_DELIMITER);
			newBondsIndex[i] = editableMol.addBond(rGroupsIndex[i*2], rGroupsIndex[i*2+1], StereoMolecule.cBondTypeSingle); // R1a-R1b and R2a-R2b
		}
		rotBondsIndex = getRotBondsIndex();
		nRotBonds = rotBondsIndex.size();
		for (int firstCut=0; firstCut<nRotBonds; firstCut++) { // to do: store the bond type and reuse it
			int atom0 = editableMol.getBondAtom(0, rotBondsIndex.get(firstCut));
			int atom1 = editableMol.getBondAtom(1, rotBondsIndex.get(firstCut));
			editableMol.setBondType(newBondsIndex[0], editableMol.getBondType(rotBondsIndex.get(firstCut)));
			editableMol.setBondAtom(1, rotBondsIndex.get(firstCut), rGroupsIndex[0]); // X-R1a
			editableMol.setBondAtom(0, newBondsIndex[0], atom1); // Y-R1b
			StereoMolecule[] scFragments = editableMol.getFragments();
			MoleculeIndexID singleCutFragments = processFragments(scFragments, new int[]{rotBondsIndex.get(firstCut)}, new int[]{atom1}, SINGLE_CUT);
			if ((KEYS_MIN_ATOMS == null || singleCutFragments.keysIDAtoms[0] >= KEYS_MIN_ATOMS) && (VALUE_MAX_ATOMS == null || singleCutFragments.valueIDAtoms <= VALUE_MAX_ATOMS)) {
				moleculeIndexesID.add(singleCutFragments);
			}
			if ((KEYS_MIN_ATOMS == null || singleCutFragments.valueIDAtoms >= KEYS_MIN_ATOMS) && (VALUE_MAX_ATOMS == null || singleCutFragments.keysIDAtoms[0] <= VALUE_MAX_ATOMS)) {
				MoleculeIndexID invertedSingleCutFragments = new MoleculeIndexID(new String[]{singleCutFragments.valueID}, singleCutFragments.keysID[0], new int[]{singleCutFragments.valueIDAtoms}, singleCutFragments.keysIDAtoms[0], new int[]{rotBondsIndex.get(firstCut)}, new int[]{atom0});
				moleculeIndexesID.add(invertedSingleCutFragments);
			}
			if (KEYS_MIN_ATOMS == null || singleCutFragments.keysIDAtoms[0] >= KEYS_MIN_ATOMS) {
				hFragments.add(singleCutFragments.keysID[0]);
			}
			if (KEYS_MIN_ATOMS == null || singleCutFragments.valueIDAtoms >= KEYS_MIN_ATOMS) {
				hFragments.add(singleCutFragments.valueID);
			}
			if (NUMBER_OF_CUTS >= DOUBLE_CUT) {
				for (int secondCut=firstCut+1; secondCut<nRotBonds; secondCut++) {
					int atom2 = editableMol.getBondAtom(1, rotBondsIndex.get(secondCut));
					editableMol.setBondType(newBondsIndex[1], editableMol.getBondType(rotBondsIndex.get(secondCut)));
					editableMol.setBondAtom(1, rotBondsIndex.get(secondCut), rGroupsIndex[2]);
					editableMol.setBondAtom(0, newBondsIndex[1], atom2);
					StereoMolecule[] dcFragments = editableMol.getFragments();
					MoleculeIndexID doubleCutFragments = processFragments(dcFragments, new int[]{rotBondsIndex.get(firstCut), rotBondsIndex.get(secondCut)}, new int[]{atom1, atom2}, DOUBLE_CUT);
					if ((KEYS_MIN_ATOMS == null || (doubleCutFragments.keysIDAtoms[0] >= KEYS_MIN_ATOMS && doubleCutFragments.keysIDAtoms[1] >= KEYS_MIN_ATOMS)) && (VALUE_MAX_ATOMS == null || doubleCutFragments.valueIDAtoms <= VALUE_MAX_ATOMS)) {
						moleculeIndexesID.add(doubleCutFragments);
					}
					editableMol.setBondAtom(1, rotBondsIndex.get(secondCut), atom2);
					editableMol.setBondAtom(0, newBondsIndex[1], rGroupsIndex[2]);
				}
			}
			editableMol.setBondAtom(1, rotBondsIndex.get(firstCut), atom1); // X-Y
			editableMol.setBondAtom(0, newBondsIndex[0], rGroupsIndex[0]); // R1a-R1b
		}
		for (Entry<String, StereoMolecule> cursor : moleculeFragments.entrySet()) {
			String fragmentID = cursor.getKey();
			StereoMolecule moleculeFragment = cursor.getValue();
			for (int atom=moleculeFragment.getAtoms()-1; atom>=0; atom--) {
				if (moleculeFragment.getAtomicNo(atom) == 0 || moleculeFragment.getAtomicNo(atom) >= 142) {
					moleculeFragment.setAtomicNo(atom, 1); //fragment.deleteAtom(atom); // this doesn't work for Cl-R1: both atoms get deleted!!
					moleculeFragment.setAtomCustomLabel(atom, (String)null); // otherwise the atom still contains the custom label, and the IDCode won't be canonical (not identical to the whole molecule)
					break;
				}
			}
			moleculeFragmentsID.add(new String[]{fragmentID, moleculeFragment.getIDCode()});
		}
		if (generateWholeMoleculeVariations) {
			generateWholeMoleculeVariations();
		}
	}
	
	/**
	 * Generates all variations of Hydrogen replacements on the whole molecule;<br>
	 * fragments are added to the moleculeIndexesID.<br>
	 */
	public void generateWholeMoleculeVariations() {
		HashMap<String, ArrayList<Integer>> wholeMoleculeVariations = new HashMap<String, ArrayList<Integer>>(); // in case of symmetric part of a molecule, two variations might be identical -> we use HashMap to avoid it
		int atomCount = mol.getAtoms();
		StereoMolecule editableMol = new StereoMolecule();
		mol.copyMolecule(editableMol);
		int rGroupIndex = editableMol.addAtom(FRAGMENT_ATOMIC_NO);
		editableMol.setAtomCustomLabel(rGroupIndex, FRAGMENT_DELIMITER);
		editableMol.ensureHelperArrays(StereoMolecule.cHelperNeighbours);
		for (int atom=0; atom<editableMol.getAtoms()-1; atom++) {
			if (mol.getPlainHydrogens(atom) > 0) {
				int newBond = editableMol.addBond(rGroupIndex, atom, StereoMolecule.cBondTypeSingle);
				String idCode = getIDCodeWithCustomLabels(editableMol);
				ArrayList<Integer> atomIndexes = new ArrayList<Integer>();
				if (wholeMoleculeVariations.containsKey(idCode)) {
					atomIndexes = wholeMoleculeVariations.get(idCode);
				}
				atomIndexes.add(atom);
				wholeMoleculeVariations.put(idCode, atomIndexes);
				editableMol.deleteBond(newBond);
			}
		}
		for (Entry<String, ArrayList<Integer>> cursor: wholeMoleculeVariations.entrySet()) {
			int[] atomIndexes = new int[cursor.getValue().size()];
			  for(int i=0; i<atomIndexes.length; i++)
				  atomIndexes[i] = cursor.getValue().get(i);
			MoleculeIndexID hFrag = new MoleculeIndexID(new String[]{cursor.getKey()}, r1H, new int[]{atomCount}, 0, new int[]{-1}, atomIndexes);
			moleculeIndexesID.add(hFrag);
		}
	}
	
	public StereoMolecule getMol() {
		return mol;
	}
}
