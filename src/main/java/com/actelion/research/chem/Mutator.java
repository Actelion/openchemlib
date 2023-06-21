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

import com.actelion.research.chem.coords.CoordinateInventor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;

public class Mutator {
	public static final int MUTATION_GROW = Mutation.MUTATION_ADD_ATOM
										  | Mutation.MUTATION_INSERT_ATOM;
	public static final int MUTATION_SHRINK = Mutation.MUTATION_CUTOUT_ATOM
											| Mutation.MUTATION_DELETE_ATOM
											| Mutation.MUTATION_DELETE_SUBSTITUENT
											| Mutation.MUTATION_CUTOUT_SFRAGMENT;
	public static final int MUTATION_KEEP_SIZE = Mutation.MUTATION_CHANGE_ATOM
											   | Mutation.MUTATION_CLOSE_RING
											   | Mutation.MUTATION_CHANGE_BOND
											   | Mutation.MUTATION_DELETE_BOND
											   | Mutation.MUTATION_CHANGE_RING
											   | Mutation.MUTATION_CLOSE_RING_AND_AROMATIZE
											   | Mutation.MUTATION_TOGGLE_AMID_SULFONAMID
											   | Mutation.MUTATION_MIGRATE
											   | Mutation.MUTATION_SWAP_SUBSTITUENT
											   | Mutation.MUTATION_INVERT_PARITY;

	public static final int MUTATION_ANY = MUTATION_GROW
										 | MUTATION_SHRINK
										 | MUTATION_KEEP_SIZE;

	private static final int cMinRingClosureSize = 3;
	private static final int cMaxRingClosureSize = 7;

	private static final int cDefaultMinAtoms = 4;
	private static final int cDefaultOptAtoms = 9;
	private static final int cDefaultMaxAtoms = 24;
	private static final double cProbabilityFactor = 1.0f;	// larger/smaller values than 1.0 over/under-express natural probabilities

	// Increase back&forth probabilities to better represent natural equilibriums where mutation path to one of two states is unlikely
	private static final double BOOST_CHANGE_RING = 1.0;
	private static final double BOOST_TOGGLE_AMID_SULFONAMID = 10.f;

	private static final double BOOST_CLOSE_RING = 1.0f;

	private static final double cMinEductProbability = 0.00001f;

	private Random mRandom;

	private StereoMolecule	    mMolCopy 	= new StereoMolecule();
	private StereoMolecule	    mMol		= new StereoMolecule();
	private AtomTypeList		mAtomTypeList;
    private Canonizer           mCanonizer;
    private boolean[]           mIsPrimaryAtom;
	private int					mMinAtoms,mOptAtoms,mMaxAtoms;
	private double				mGrowBoost;
	private MutationBiasProvider mBiasProvider;

	/**
	 * Creates a new Mutator that calculates probabilities of individual mutations from
	 * a given atom type file. For every potential mutation the algorithm calculates
	 * a probability considering the frequencies of broken and formed atom types during the conversion.<br>
	 * If no type list is given, all possible mutations are considered equally likely.<br>
	 * New type files can be created from an sdf or dwar file like this:<br>
	 * <i>new AtomTypeList("/somepath/chembl14.dwar", AtomTypeCalculator.cPropertiesForMutator).writeTypeFile("/somepath/chembl14.typ");<i><br>
	 * This would cause a file <i>/somepath/chembl14.typ</i> to be created with the statistics taken from <i>chembl14.dwar</i>.
	 * @param filename null or name of a atomTypeList file ('*.typ')
	 */
	public Mutator(String filename) {
		if (filename != null) {
	        try {
	            mAtomTypeList = new AtomTypeList(filename, AtomTypeCalculator.cPropertiesForMutator);
	            mAtomTypeList.calculateProbabilities();
	            }
	        catch (Exception e) {
	            e.printStackTrace();
	            }
	        }

	    mRandom = new Random();

		mGrowBoost = 1.0;

		mMinAtoms = cDefaultMinAtoms;
		mOptAtoms = cDefaultOptAtoms;
		mMaxAtoms = cDefaultMaxAtoms;
		}

	/**
	 * Creates a new Mutator that calculates probabilities of individual mutations from
	 * a given atom type list. For every potential mutation the algorithm calculates
	 * a probability considering the frequencies of broken and formed atom types during the conversion.<br>
	 * If no type list is given, all possible mutations are considered equally likely.<br>
	 * New type files can be created from an sdf or dwar file like this:<br>
	 * <i>new AtomTypeList().create("/somepath/chembl14.dwar", AtomTypeCalculator.cPropertiesForMutator);<i><br>
	 * This would cause a file <i>/somepath/chembl14.typ</i> to be created with the statistics taken from <i>chembl14.dwar</i>.
	 * Instantiating multiple Mutator objects from the same AtomTypeList is thread-safe.
	 * @param atomTypeList null or name of a atomTypeList file ('*.typ')
	 */
	public Mutator(AtomTypeList atomTypeList) {
		mAtomTypeList = atomTypeList;
		mAtomTypeList.calculateProbabilities();

		mRandom = new Random();

		mGrowBoost = 1.0;

		mMinAtoms = cDefaultMinAtoms;
		mOptAtoms = cDefaultOptAtoms;
		mMaxAtoms = cDefaultMaxAtoms;
		}

	public void setBiasProvider(MutationBiasProvider mbp) {
		mBiasProvider = mbp;
		}

	public void setGrowBoost(double boost) {
		mGrowBoost = boost;
		}

	public void setPreferredSize(int minAtoms, int preferredAtoms, int maxAtoms) {
		mMinAtoms = minAtoms;
		mOptAtoms = preferredAtoms;
		mMaxAtoms = maxAtoms;
		}

	public StereoMolecule[] getMutatedSet(StereoMolecule mol,
                                          int mutationType,
                                          boolean regulateSize,
                                          int count) {
	    ArrayList<Mutation> mutationList = generateMutationList(mol, mutationType, regulateSize);

        if (count > mutationList.size())
	        count = mutationList.size();

        StereoMolecule[] set = new StereoMolecule[count];
	    for (int i=0; i<count; i++) {
	        set[i] = new StereoMolecule(mol);
	        mutate(set[i], mutationList);
	        }
	    return set;
	    }

	/**
	 * Does an in-place mutation of the molecule allowing any kind of mutation
	 * at any of the molecules non-selected atoms aiming for 4-24 non-H atoms with an optimum of 9 atoms.
	 * @param mol
     * @return list of all not-used mutations or null if no mutation was done because no possible mutation was found
	 */
	public ArrayList<Mutation> mutate(StereoMolecule mol) {
		return mutate(mol, MUTATION_ANY, true);
		}

	/**
	 * Does an in-place mutation of the molecule allowing the defined mutation kinds at any of the molecules unselected atoms.
	 * @param mol
	 * @param mutationType MUTATION_ANY, MUTATION_GROW, MUTATION_KEEP_SIZE, or MUTATION_SHRINK or other combination of allowed mutations types
	 * @param regulateSize whether to regulate the molecule size within 4 to 24 non-H atoms or what was defined by setPreferredSize()
	 * @return list of all not-used mutations or null if no mutation was done because no possible mutation was found
	 */
	public ArrayList<Mutation> mutate(StereoMolecule mol,
	                                  int mutationType,
	                                  boolean regulateSize) {
	    ArrayList<Mutation> ml = generateMutationList(mol, mutationType, regulateSize);

	    if (ml.size() == 0) {
	        System.out.println("no possible mutation found. ID-Code:"+new Canonizer(mol, Canonizer.ENCODE_ATOM_SELECTION).getIDCode());
	        return null;
	        }

	    mutate(mol, ml);
	    return ml;
	    }

	/**
     * Selects a likely mutation from the list, performs the mutation and removes it from the list.
	 * If the mutation list is empty, then no mutation is performed.
     * @param mol
     * @param mutationList
     */
    public void mutate(StereoMolecule mol, ArrayList<Mutation> mutationList) {
    	if (!mutationList.isEmpty()) {
    		Mutation mutation = selectLikelyMutation(mutationList);
			if (mutation != null) {
				performMutation(mol, mutation);
				}
			}
        }

	/**
	 * Creates a list of possible mutations and their probabilities
	 * @param mol
	 * @param mutationType MUTATION_ANY, MUTATION_GROW, MUTATION_KEEP_SIZE, or MUTATION_SHRINK or other combination of allowed mutations types
	 * @param regulateSize if true keeps non-H atoms between 4 and 24 with an optimum at 9 or what was defined by setPreferredSize().
	 */
	public ArrayList<Mutation> generateMutationList(StereoMolecule mol,
	                                                int mutationType,
	                                                boolean regulateSize) {
	    mMol = mol;

	    if (mBiasProvider != null)
	    	mBiasProvider.setBiasReference(mol);

	    detectSymmetry();

		ArrayList<Mutation> mutationList = new ArrayList<Mutation>();
	    ArrayList<MutatorSubstituent> substituentList = createSubstituentList();

		if (!regulateSize
		 || (mRandom.nextDouble() > (float)(mMol.getAtoms() - mOptAtoms)
		 						  / (float)(mMaxAtoms - mOptAtoms))) {
			if ((mutationType & Mutation.MUTATION_ADD_ATOM) != 0)
				addProbabilitiesForAddAtom(mutationList);
			if ((mutationType & Mutation.MUTATION_INSERT_ATOM) != 0)
				addProbabilitiesForInsertAtom(mutationList);
			}

		if ((mutationType & Mutation.MUTATION_CHANGE_ATOM) != 0)
			addProbabilitiesForChangeAtom(mutationList);

		if (!regulateSize
		 || (mRandom.nextDouble() < (float)(mMol.getAtoms() - mMinAtoms)
		 						  / (float)(mOptAtoms - mMinAtoms))) {
			if ((mutationType & Mutation.MUTATION_DELETE_ATOM) != 0)
				addProbabilitiesForDeleteAtom(mutationList);
			if ((mutationType & Mutation.MUTATION_CUTOUT_ATOM) != 0)
				addProbabilitiesForCutOutAtom(mutationList);
            if ((mutationType & Mutation.MUTATION_DELETE_SUBSTITUENT) != 0)
                addProbabilitiesForDeleteSubstituent(mutationList, substituentList);
	        if ((mutationType & Mutation.MUTATION_CUTOUT_SFRAGMENT) != 0)
	            addProbabilitiesForCutOutFragment(mutationList);
			}

		if ((mutationType & (Mutation.MUTATION_CLOSE_RING | Mutation.MUTATION_CLOSE_RING_AND_AROMATIZE)) != 0)
			addProbabilitiesForCloseRing(mutationList, mutationType);
		if ((mutationType & (Mutation.MUTATION_TOGGLE_AMID_SULFONAMID)) != 0)
			addProbabilitiesForToggleAmidSulfonamid(mutationList);
		if ((mutationType & Mutation.MUTATION_CHANGE_BOND) != 0)
			addProbabilitiesForChangeBond(mutationList);
		if ((mutationType & Mutation.MUTATION_DELETE_BOND) != 0)
			addProbabilitiesForDeleteBond(mutationList);
		if ((mutationType & Mutation.MUTATION_CHANGE_RING) != 0)
			addProbabilitiesForChangeRing(mutationList);
		if ((mutationType & Mutation.MUTATION_MIGRATE) != 0)
			addProbabilitiesForMigrate(mutationList);
        if ((mutationType & Mutation.MUTATION_SWAP_SUBSTITUENT) != 0)
            addProbabilitiesForSwapSubstituent(mutationList, substituentList);
		if ((mutationType & Mutation.MUTATION_INVERT_PARITY) != 0)
			addProbabilitiesForInvertParity(mutationList);

        if (cProbabilityFactor != 1)
        	tweakProbabilities(mutationList);

		return mutationList;
	    }

    private void tweakProbabilities(ArrayList<Mutation> mutationList) {
		for (Mutation mutation:mutationList)
			mutation.mProbability = Math.pow(mutation.mProbability, cProbabilityFactor);
		}

	private void detectSymmetry() {
        mCanonizer = new Canonizer(mMol, Canonizer.CREATE_SYMMETRY_RANK);
		mIsPrimaryAtom = new boolean[mMol.getAtoms()];

        // here we count unselected atoms per rank
        int[] rankCount = new int[1+mMol.getAtoms()];
        for (int atom=0; atom<mMol.getAtoms(); atom++)
            if (!mMol.isSelectedAtom(atom))
				rankCount[mCanonizer.getSymmetryRank(atom)]++;

        int[] graphAtom = new int[mMol.getAtoms()];

		// then we traverse all non-selected atoms:
		// - if the rank was not used than this is a primary atom to be used as representative of this rank
		// - if a primary atom is discovered, then recursively define all connected unselected atoms with an unused rank also as primary
        for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (!mMol.isSelectedAtom(atom)) {
				int rank = mCanonizer.getSymmetryRank(atom);
				if (rankCount[rank] != 0) {
					mIsPrimaryAtom[atom] = true;
					rankCount[rank] = 0;

					// flag all connected unselected atoms with unassigned ranks
					int current = 0;
					int highest = 0;
					graphAtom[0] = atom;
					while (current <= highest) {
						for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
							int connAtom = mMol.getConnAtom(graphAtom[current], i);
							if (!mMol.isSelectedAtom(connAtom)) {
								int connRank = mCanonizer.getSymmetryRank(connAtom);
								if (rankCount[connRank] != 0) {
									mIsPrimaryAtom[connAtom] = true;
									rankCount[connRank] = 0;
									graphAtom[highest++] = connAtom;
									}
								}
							}
						current++;
						}
					}
				}
			}
	    }

	private void addProbabilitiesForAddAtom(ArrayList<Mutation> mutationList){
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mIsPrimaryAtom[atom]) {
				int freeValence = mMol.getFreeValence(atom);
				int maxBondOrder = (freeValence > 3) ? 3 : freeValence;
				if( maxBondOrder > 0) {
					mMol.copyMolecule(mMolCopy);
					int newAtom = mMolCopy.addAtom(6);
					int newBond = mMolCopy.addBond(atom, newAtom, Molecule.cBondTypeSingle);

					//loop for single float or triple bond to add
					for (int bondOrder=1; bondOrder<=maxBondOrder; bondOrder++) {
						//check if the atom to add has enough bonds
						boolean allowed = false;
						for (int i=0; i<Mutation.cAllowedAtomicNo[bondOrder-1].length; i++) {
							if (Mutation.cAllowedAtomicNo[bondOrder-1][i] == mMol.getAtomicNo(atom)) {
								allowed = true;
								break;
								}
							}
						if (!allowed)
							break;

						for (int i=0; i<Mutation.cAllowedAtomicNo[bondOrder-1].length;i++) {
							double f_Educt = Math.max(cMinEductProbability, getFrequency(mMol, atom));

							int atomicNo	= Mutation.cAllowedAtomicNo[bondOrder-1][i];
							int bondType	= getBondTypeFromOrder(bondOrder);
							mMolCopy.setAtomicNo(newAtom, atomicNo);
							mMolCopy.setBondType(newBond, bondType);
							mMolCopy.ensureHelperArrays(Molecule.cHelperRings);
							double f_Product = getFrequency(mMolCopy, atom);

							double p = mGrowBoost * f_Product / f_Educt * Math.sqrt(getFrequency(mMolCopy, newAtom));
							if (p > 0.0 && isValidStructure(mMolCopy)) {
								if (mBiasProvider != null)
									p *= mBiasProvider.getBiasFactor(mMolCopy);

								mutationList.add(new Mutation(Mutation.MUTATION_ADD_ATOM, atom, -1, atomicNo, bondType, p));
								}
							}
						}
					}
				}
			}
		}

	private void addProbabilitiesForInsertAtom(ArrayList<Mutation> mutationList) {
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.getBondType(bond) == Molecule.cBondTypeSingle
			 && !mMol.isAromaticBond(bond)) {
				for (int i=0; i<2; i++) {
					int atom1 = mMol.getBondAtom(i, bond);
					int atom2 = mMol.getBondAtom(1-i, bond);

					if (mIsPrimaryAtom[atom1]
					 && !mMol.isSelectedAtom(atom2)
					 && mCanonizer.getSymmetryRank(atom1) < mCanonizer.getSymmetryRank(atom2)) {
						double f_Educt1 = Math.max(cMinEductProbability, getFrequency(mMol, atom1));
						double f_Educt2 = Math.max(cMinEductProbability, getFrequency(mMol, atom2));

						mMol.copyMolecule(mMolCopy);
						mMolCopy.deleteBond(bond);
						int newAtom = mMolCopy.addAtom(6);
						mMolCopy.addBond(atom1, newAtom, Molecule.cBondTypeSingle);
						mMolCopy.addBond(atom2, newAtom, Molecule.cBondTypeSingle);

						for (int j=0; j<Mutation.cAllowedAtomicNo[1].length;j++) {
							int atomicNo = Mutation.cAllowedAtomicNo[1][j];
							mMolCopy.setAtomicNo(newAtom, atomicNo);

							double f_Product1 = getFrequency(mMolCopy, atom1);
							double f_Product2 = getFrequency(mMolCopy, atom2);
							double f_NewAtom = getFrequency(mMolCopy, newAtom);

							double p = mGrowBoost * Math.sqrt(f_NewAtom * f_Product1 * f_Product2 / (f_Educt1 * f_Educt2));
							if (p > 0.0 && isValidStructure(mMolCopy)) {
								if (mBiasProvider != null)
									p *= mBiasProvider.getBiasFactor(mMolCopy);

								mutationList.add(new Mutation(Mutation.MUTATION_INSERT_ATOM, bond, -1, atomicNo, -1, p));
								}
							}
						}
					}
				}
			}
		}

	private void addProbabilitiesForChangeAtom(ArrayList<Mutation> mutationList) {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mIsPrimaryAtom[atom]) {
				int maxBondOrder = 1;
				for (int i=0; i<mMol.getConnAtoms(atom); i++)
					if (maxBondOrder < mMol.getConnBondOrder(atom, i))
						maxBondOrder = mMol.getConnBondOrder(atom, i);

				for (int i=0; i<Mutation.cAllowedAtomicNo[maxBondOrder-1].length; i++) {
					int proposedAtomicNo = Mutation.cAllowedAtomicNo[maxBondOrder-1][i];
					int valences = mMol.getConnAtoms(atom) + mMol.getAtomPi(atom);

					if (mMol.getAtomicNo(atom) == proposedAtomicNo)
						continue;

					if (valences > 1
					 && (proposedAtomicNo == 9
					  || proposedAtomicNo == 17
					  || proposedAtomicNo == 35
					  || proposedAtomicNo == 53))
						continue;

					if (valences > 2
					 && proposedAtomicNo == 8)
						continue;

					if (valences > 3
					 && proposedAtomicNo == 5)
						continue;

					if (valences > 4
					 && (proposedAtomicNo == 6
					  || proposedAtomicNo == 7))
						continue;

					if (valences > 5
					 && proposedAtomicNo == 15)
						continue;

					mMol.copyMolecule(mMolCopy);
					mMolCopy.setAtomicNo(atom, proposedAtomicNo);
					mMolCopy.ensureHelperArrays(Molecule.cHelperRings);

					double f_Educt = Math.max(cMinEductProbability, getFrequency(mMol,atom));
					double f_Product = getFrequency(mMolCopy,atom);

					for (int j=0; j<mMol.getConnAtoms(atom); j++) {
						int connAtom = mMol.getConnAtom(atom, j);
						f_Educt *=  Math.max(cMinEductProbability, getFrequency(mMol, connAtom));
						f_Product *= getFrequency(mMolCopy,connAtom);
						}

					double p = Math.pow(f_Product / f_Educt, 1.0 / (1 + mMol.getConnAtoms(atom)));
					if (p > 0.0 && isValidStructure(mMolCopy)) {
						if (mBiasProvider != null)
							p *= mBiasProvider.getBiasFactor(mMolCopy);

						mutationList.add(new Mutation(Mutation.MUTATION_CHANGE_ATOM, atom, -1, proposedAtomicNo, -1, p));
						}
					}
				}
			}
		}

	private void addProbabilitiesForDeleteAtom(ArrayList<Mutation> mutationList) {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mIsPrimaryAtom[atom]) {
				if (mMol.getConnAtoms(atom)==1) {
					int connAtom = mMol.getConnAtom(atom,0);
					if (!mMol.isSelectedAtom(connAtom)) {
						mMol.copyMolecule(mMolCopy);
						mMolCopy.deleteBond(mMol.getConnBond(atom, 0));

						double f_Educt = Math.max(cMinEductProbability, getFrequency(mMol, connAtom));
						double f_Product = getFrequency(mMolCopy, connAtom);
						double f_Deleted = Math.max(cMinEductProbability, getFrequency(mMol, atom));

						double p = f_Product / f_Educt / Math.sqrt(f_Deleted);
						if (p > 0.0 && isValidStructure(mMolCopy)) {
							if (mBiasProvider != null)
								p *= mBiasProvider.getBiasFactor(mMolCopy);

							mutationList.add(new Mutation(Mutation.MUTATION_DELETE_ATOM, atom, -1, -1, -1, p));
							}
						}
					}
				}
			}
		}

	private void addProbabilitiesForCutOutAtom(ArrayList<Mutation> mutationList) {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mIsPrimaryAtom[atom]
			 && !mMol.isAromaticAtom(atom)
			 && mMol.getConnAtoms(atom) == 2
			 && mMol.getAtomPi(atom) == 0
			 && mMol.getAtomRingSize(atom) != 3) {
				int atom1 = mMol.getConnAtom(atom, 0);
				int atom2 = mMol.getConnAtom(atom, 1);

				double f_Educt1 = Math.max(cMinEductProbability, getFrequency(mMol, atom1));
				double f_Educt2 = Math.max(cMinEductProbability, getFrequency(mMol, atom2));
				double f_Deleted = Math.max(cMinEductProbability, getFrequency(mMol, atom));

				mMol.copyMolecule(mMolCopy);
				mMolCopy.addBond(atom1, atom2, Molecule.cBondTypeSingle);
				mMolCopy.deleteAtom(atom);

				int newBond = mMolCopy.getAllBonds() - 1;
				atom1 = mMolCopy.getBondAtom(0, newBond);
				atom2 = mMolCopy.getBondAtom(1, newBond);

				double f_Product1 = getFrequency(mMolCopy, atom1);
				double f_Product2 = getFrequency(mMolCopy, atom2);

				double p = Math.sqrt(f_Product1 * f_Product2 / (f_Educt1 * f_Educt2)) / Math.sqrt(f_Deleted);
				if (p > 0.0 && isValidStructure(mMolCopy)) {
					if (mBiasProvider != null)
						p *= mBiasProvider.getBiasFactor(mMolCopy);

					mutationList.add(new Mutation(Mutation.MUTATION_CUTOUT_ATOM, atom, -1, -1, -1, p));
					}
				}
			}
		}//end_deleting

	private void addProbabilitiesForToggleAmidSulfonamid(ArrayList<Mutation> mutationList) {
		int[] oxygen = new int[4];
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mIsPrimaryAtom[atom]) {
				int atomicNo = mMol.getAtomicNo(atom);
				if (atomicNo != 6 && atomicNo != 16)
					continue;

				int oxoCount = 0;
				for (int i=0; i<mMol.getConnAtoms(atom); i++) {
					oxygen[oxoCount] = mMol.getConnAtom(atom, i);
					if (mMol.getConnBondOrder(atom, i) == 2 && mMol.getAtomicNo(oxygen[oxoCount]) == 8)
						oxoCount++;
					}

				if ((atomicNo == 6 && oxoCount == 1)
				 || (atomicNo == 16 && oxoCount == 2)) {
					mMol.copyMolecule(mMolCopy);
					mMolCopy.setAtomicNo(atom, atomicNo == 6 ? 16 : 6);
					if (atomicNo == 6) {
						int newOxygen = mMolCopy.addAtom(8);
						mMolCopy.addBond(atom, newOxygen, Molecule.cBondTypeDouble);
						}
					else {
						mMolCopy.deleteBond(mMol.getBond(atom, oxygen[1]));
						// don't delete the atom here to not disturb the atom table for scoring
						}

					double f_Educt = Math.max(cMinEductProbability, getFrequency(mMol, atom));
					double f_Product = getFrequency(mMolCopy, atom);
					int count = 1;
					for (int i=0; i<mMol.getConnAtoms(atom); i++) {
						int connAtom = mMol.getConnAtom(atom, i);
						boolean isOxo = false;
						for (int j=0; j<oxoCount; j++) {
							if (connAtom != oxygen[j]) {
								isOxo = true;
								break;
								}
							}
						if (!isOxo) {
							f_Educt *= Math.max(cMinEductProbability, getFrequency(mMol, connAtom));
							f_Product *= getFrequency(mMolCopy, connAtom);
							count++;
							}
						}

					double p = BOOST_TOGGLE_AMID_SULFONAMID * Math.pow(f_Product / f_Educt, 1.0 / count);
					if (p > 0.0 && isValidStructure(mMolCopy)) {
						if (mBiasProvider != null)
							p *= mBiasProvider.getBiasFactor(mMolCopy);

						mutationList.add(new Mutation(Mutation.MUTATION_TOGGLE_AMID_SULFONAMID, atom, oxygen[1], -1, -1, p));
						}
					}
				}
			}
		}


	private void addProbabilitiesForCloseRing(ArrayList<Mutation> mutationList, int mode) {
		for (int atom1=0; atom1<mMol.getAtoms(); atom1++) {
			if (mIsPrimaryAtom[atom1]) {
				int graphAtom[] = new int[mMol.getAtoms()];
				int graphBond[] = new int[mMol.getAtoms()];
				int graphParent[] = new int[mMol.getAtoms()];
				int graphLevel[] = new int[mMol.getAtoms()];
				graphAtom[0] = atom1;
				graphLevel[atom1] = 1;
				int current = 0;
				int highest = 0;
				while (current <= highest) {
					if (graphLevel[graphAtom[current]] >= cMaxRingClosureSize)
						break;

					for (int i=0; i<mMol.getConnAtoms(graphAtom[current]); i++) {
						int candidate = mMol.getConnAtom(graphAtom[current],i);
						if (graphLevel[candidate] == 0) {
							graphParent[candidate] = graphAtom[current];
							graphLevel[candidate] = graphLevel[graphAtom[current]] + 1;
							graphBond[candidate] = mMol.getConnBond(graphAtom[current],i);
							graphAtom[++highest] = candidate;
							}
						}
					current++;
					}

				int candidateIndex = cMinRingClosureSize-1;
				while (candidateIndex <= highest && graphLevel[candidateIndex] < cMinRingClosureSize)
					candidateIndex++;

				while (candidateIndex <= highest) {
					int atom2 = graphAtom[candidateIndex++];

					if (mMol.isSelectedAtom(atom2) || atom2 < atom1)
						continue;

					if (mMol.getBond(atom1, atom2) != -1)
						continue;

					if (!qualifiesForRing(graphBond, graphLevel, graphParent, atom1, atom2))
						continue;

					for (int order=1; order<=2; order++) {
						if (mMol.getFreeValence(atom1) < order
						 || mMol.getFreeValence(atom2) < order)
							break;

						// here we also check, whether we can form a ring and aromatize it
						if ((mode & Mutation.MUTATION_CLOSE_RING_AND_AROMATIZE) != 0 && order == 1) {
							int[] ringAtom = new int[6];
							int ringSize = 1 + mMol.getPath(ringAtom, atom1, atom2, 5, null);
							if (ringSize == 5 || ringSize == 6) {
								boolean selectedAtomFound = false;
								for (int i=1; i<ringSize-1; i++) {
									if (mMol.isSelectedAtom(ringAtom[i])) {
										selectedAtomFound = true;
										break;
										}
									}
								if (!selectedAtomFound) {
									mMol.copyMolecule(mMolCopy);
									mMolCopy.addBond(atom1, atom2, 1);
									if (aromatizeRing(mMolCopy, ringAtom, ringSize)) {
										double f_Educt = 1.0;
										double f_Product = 1.0;
										for (int atom:ringAtom) {
											f_Educt /= Math.max(cMinEductProbability, getFrequency(mMol, atom));
											f_Product *= getFrequency(mMolCopy, atom);
											}

										double p = BOOST_CLOSE_RING * Math.pow(f_Product / f_Educt, 1.0 / ringAtom.length);
										if (p > 0.0 && isValidStructure(mMolCopy)) {
											if (mBiasProvider != null)
												p *= mBiasProvider.getBiasFactor(mMolCopy);

											mutationList.add(new Mutation(Mutation.MUTATION_CLOSE_RING_AND_AROMATIZE, atom1, atom2, ringSize, ringAtom, p));
											}
										}
									}
								}
							}

						if ((mode & Mutation.MUTATION_CLOSE_RING) != 0 && order == 1) {
							mMol.copyMolecule(mMolCopy);
							mMolCopy.addBond(atom1, atom2, getBondTypeFromOrder(order));

							double f_Educt1 = Math.max(cMinEductProbability, getFrequency(mMol, atom1));
							double f_Educt2 = Math.max(cMinEductProbability, getFrequency(mMol, atom2));
							double f_Product1 = getFrequency(mMolCopy, atom1);
							double f_Product2 = getFrequency(mMolCopy, atom2);

							double p = (mAtomTypeList == null) ? 1.0f : mAtomTypeList.getRingSizeAdjust(graphLevel[atom2]);
							p *= BOOST_CLOSE_RING * Math.sqrt(f_Product1 * f_Product2 / (f_Educt1 * f_Educt2));
							if (p > 0.0 && isValidStructure(mMolCopy)) {
								if (mBiasProvider != null)
									p *= mBiasProvider.getBiasFactor(mMolCopy);

								mutationList.add(new Mutation(Mutation.MUTATION_CLOSE_RING, atom1, atom2, getBondTypeFromOrder(order), -1, p));
								}
							}
						}
					}
				}
			}
		}


	private void addProbabilitiesForChangeBond(ArrayList<Mutation> mutationList) {
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			for (int i=0; i<2; i++) {
				int atom1 = mMol.getBondAtom(i, bond);
				if (mIsPrimaryAtom[atom1]) {
					int atom2 = mMol.getBondAtom(1-i, bond);
					if (!mMol.isSelectedAtom(atom2)) {
						if (mIsPrimaryAtom[atom2]) {
							if (atom2 < atom1)
								continue;    // this is handled with other i
							}
						else {	// we must have a symmetrical bond somewhere else
							if (mCanonizer.getSymmetryRank(atom1) < mCanonizer.getSymmetryRank(atom2))
								continue;
							}

						int minFreeValence = mMol.getFreeValence(atom1);
						if (minFreeValence > mMol.getFreeValence(atom2))
							minFreeValence = mMol.getFreeValence(atom2);
						int maxBondOrder = mMol.getBondOrder(bond) + minFreeValence;

						if (maxBondOrder > 3)
							maxBondOrder = 3;

						if (mMol.isAromaticBond(bond)) {
							if (mMol.getBondOrder(bond) == 1)
								continue;
							maxBondOrder = 2;
							}

						if (mMol.isSmallRingBond(bond)) {
							if (mMol.getBondOrder(bond) == 1
							 && mMol.getAtomPi(atom1)+mMol.getAtomPi(atom2) != 0)
								maxBondOrder = 1;
							else if (maxBondOrder > 2)
								maxBondOrder = 2;
							}

						if (maxBondOrder == 2
						 && (mMol.getAtomicNo(atom1) < 5
						  || (mMol.getAtomicNo(atom1) > 8
						   && mMol.getAtomicNo(atom1) != 15
						   && mMol.getAtomicNo(atom1) != 16)
						  || mMol.getAtomicNo(atom2) < 5
						  || (mMol.getAtomicNo(atom2) > 8
						   && mMol.getAtomicNo(atom2) != 15
						   && mMol.getAtomicNo(atom2) != 16)))
							maxBondOrder = 1;

						for (int bondOrder=1; bondOrder<=maxBondOrder; bondOrder++) {
							if (bondOrder == mMol.getBondOrder(bond))
								continue;

							mMol.copyMolecule(mMolCopy);
							mMolCopy.setBondType(bond, getBondTypeFromOrder(bondOrder));

							double f_Educt1 = Math.max(cMinEductProbability, getFrequency(mMol, atom1));
							double f_Educt2 = Math.max(cMinEductProbability, getFrequency(mMol, atom2));
							double f_Product1 = getFrequency(mMolCopy, atom1);
							double f_Product2 = getFrequency(mMolCopy, atom2);

							double p = Math.sqrt(f_Product1 * f_Product2 / (f_Educt1 * f_Educt2));
							if (p > 0.0 && isValidStructure(mMolCopy)) {
								if (mBiasProvider != null)
									p *= mBiasProvider.getBiasFactor(mMolCopy);

								mutationList.add(new Mutation(Mutation.MUTATION_CHANGE_BOND, bond, -1, getBondTypeFromOrder(bondOrder), -1, p));
								}
							}
						}
					}
				}
			}
		}


	private void addProbabilitiesForDeleteBond(ArrayList<Mutation> mutationList) {
		for (int bond=0; bond<mMol.getBonds(); bond++) {
			if (mMol.isRingBond(bond)) {
				for (int i=0; i<2; i++) {
					int atom1 = mMol.getBondAtom(i, bond);
					if (mIsPrimaryAtom[atom1]) {
						int atom2 = mMol.getBondAtom(1-i, bond);
						if (!mMol.isSelectedAtom(atom2)) {
							if (mIsPrimaryAtom[atom2]) {
								if (atom2 < atom1)
									continue;    // this is handled with other i
								}
							else {	// we must have a symmetrical bond somewhere else
								if (mCanonizer.getSymmetryRank(atom1) < mCanonizer.getSymmetryRank(atom2))
									continue;
								}

							mMol.copyMolecule(mMolCopy);
							mMolCopy.deleteBond(bond);

							double f_Educt1 = Math.max(cMinEductProbability, getFrequency(mMol, atom1));
							double f_Educt2 = Math.max(cMinEductProbability, getFrequency(mMol, atom2));
							double f_Product1 = getFrequency(mMolCopy, atom1);
							double f_Product2 = getFrequency(mMolCopy, atom2);

							double p = Math.sqrt(f_Product1 * f_Product2 / (f_Educt1 * f_Educt2));
							if (p > 0.0 && isValidStructure(mMolCopy)) {
								if (mBiasProvider != null)
									p *= mBiasProvider.getBiasFactor(mMolCopy);

								mutationList.add(new Mutation(Mutation.MUTATION_DELETE_BOND, bond, -1, -1, -1, p));
								}
							}
						}
					}
				}
			}
		}


	private void addProbabilitiesForChangeRing(ArrayList<Mutation> mutationList) {
		mMol.ensureHelperArrays(Molecule.cHelperRings);
		RingCollection ringSet = mMol.getRingSet();
		for (int ring=0; ring<ringSet.getSize(); ring++) {
			int ringSize = ringSet.getRingSize(ring);
			if (ringSize == 5 || ringSize == 6) {
				int[] ringAtom = ringSet.getRingAtoms(ring);

				boolean fixedAtomFound = false;
				for (int i=0; i<ringSize; i++)
					if (mMol.isSelectedAtom(ringAtom[i]))
						fixedAtomFound = true;
				if (fixedAtomFound)
					continue;

				int maxRank = 0;
				int maxAtom = -1;
				for (int i=0; i<ringSize; i++) {
					if (maxRank < mCanonizer.getSymmetryRank(ringAtom[i])) {
						maxRank = mCanonizer.getSymmetryRank(ringAtom[i]);
						maxAtom = ringAtom[i];
						}
					}
				if (mIsPrimaryAtom[maxAtom]) {
					int[] ringBond = ringSet.getRingBonds(ring);
					if (hasExocyclicPiBond(ringAtom, ringBond))
						continue;

					for (int heteroPosition=0; heteroPosition<ringSize; heteroPosition++) {
						if (heteroPosition > 0
						 && (ringSize == 6
						  || ringSet.isAromatic(ring)))
							break;

						if (ringSize == 5
						 && mMol.getAtomicNo(ringAtom[heteroPosition]) != 7
						 && mMol.getAtomicNo(ringAtom[heteroPosition]) != 8
						 && mMol.getAtomicNo(ringAtom[heteroPosition]) != 16)
							continue;

						mMol.copyMolecule(mMolCopy);
						mMolCopy.ensureHelperArrays(Molecule.cHelperRings);
						if (!changeAromaticity(mMolCopy, ring, heteroPosition))
							continue;

						double f_Educt = 1.0;
						double f_Product = 1.0;
						for (int atom=0; atom<ringSize; atom++) {
							f_Educt *= Math.max(cMinEductProbability, getFrequency(mMol, atom));
							f_Product *= getFrequency(mMolCopy, atom);
							}

						double p = BOOST_CHANGE_RING * Math.pow(f_Product / f_Educt, 1.0 / ringSize);
						if (p > 0.0 && isValidStructure(mMolCopy)) {
							if (mBiasProvider != null)
								p *= mBiasProvider.getBiasFactor(mMolCopy);

							mutationList.add(new Mutation(Mutation.MUTATION_CHANGE_RING, ring, -1, heteroPosition, -1, p));
							}
						}
					}
				}
			}
		}


	private void addProbabilitiesForMigrate(ArrayList<Mutation> mutationList) {
		for (int atom=0; atom<mMol.getAtoms(); atom++) {
			if (mIsPrimaryAtom[atom] && mMol.getConnAtoms(atom) > 2) {
				for (int i=0; i<mMol.getConnAtoms(atom); i++) {
					int migratingAtom = mMol.getConnAtom(atom, i);
					if (!mMol.isSelectedAtom(migratingAtom)) {
						int migratingBond = mMol.getConnBond(atom, i);
						if (!mMol.isRingBond(migratingBond)
						 && mMol.getBondOrder(migratingBond) == 1) {
							for (int j=0; j<mMol.getConnAtoms(atom); j++) {
								if (i == j)
									continue;

								int destinationAtom = mMol.getConnAtom(atom, j);
								if (mIsPrimaryAtom[destinationAtom]
								 && mMol.getFreeValence(destinationAtom) > 0) {
									mMol.copyMolecule(mMolCopy);
									for (int k=0; k<2; k++)
										if (mMolCopy.getBondAtom(k, migratingBond) == atom)
											mMolCopy.setBondAtom(k, migratingBond, destinationAtom);

									double f_Educt = Math.max(cMinEductProbability, getFrequency(mMol, atom))
												   * Math.max(cMinEductProbability, getFrequency(mMol, migratingAtom))
												   * Math.max(cMinEductProbability, getFrequency(mMol, destinationAtom));
									double f_Product = getFrequency(mMolCopy, atom)
													 * getFrequency(mMolCopy, migratingAtom)
													 * getFrequency(mMolCopy, destinationAtom);

									double p = Math.sqrt(f_Product / f_Educt);	// slight boost
									if (p > 0.0 && isValidStructure(mMolCopy)) {
										if (mBiasProvider != null)
											p *= mBiasProvider.getBiasFactor(mMolCopy);

										mutationList.add(new Mutation(Mutation.MUTATION_MIGRATE, migratingBond, -1, atom, destinationAtom, p));
										}
									}
								}
							}
						}
					}
				}
			}
		}


    private void addProbabilitiesForSwapSubstituent(ArrayList<Mutation> mutationList,
                                                    ArrayList<MutatorSubstituent> substituentList) {
        for (MutatorSubstituent s1:substituentList) {
			if (mIsPrimaryAtom[s1.firstAtom]) {
				for (MutatorSubstituent s2:substituentList) {
					if (mIsPrimaryAtom[s2.firstAtom]
					 && s1.coreAtom != s2.coreAtom
					 && s1.bond != s2.bond) {
						mMol.copyMolecule(mMolCopy);
						for (int k=0; k<2; k++) {
							if (mMolCopy.getBondAtom(k, s1.bond) == s1.firstAtom)
								mMolCopy.setBondAtom(k, s1.bond, s2.firstAtom);
							if (mMolCopy.getBondAtom(k, s2.bond) == s2.firstAtom)
								mMolCopy.setBondAtom(k, s2.bond, s1.firstAtom);
							}

						double f_Educt = Math.max(cMinEductProbability, getFrequency(mMol, s1.coreAtom))
									   * Math.max(cMinEductProbability, getFrequency(mMol, s1.firstAtom))
									   * Math.max(cMinEductProbability, getFrequency(mMol, s2.coreAtom))
									   * Math.max(cMinEductProbability, getFrequency(mMol, s2.firstAtom));
						double f_Product = getFrequency(mMolCopy, s1.coreAtom)
										 * getFrequency(mMolCopy, s1.firstAtom)
										 * getFrequency(mMolCopy, s2.coreAtom)
										 * getFrequency(mMolCopy, s2.firstAtom);

						double p = Math.sqrt(f_Product / f_Educt);	// slight boost
						if (p > 0.0 && isValidStructure(mMolCopy)) {
							if (mBiasProvider != null)
								p *= mBiasProvider.getBiasFactor(mMolCopy);

							mutationList.add(new Mutation(Mutation.MUTATION_SWAP_SUBSTITUENT, s1.coreAtom, s2.coreAtom, s1.firstAtom, s2.firstAtom, p));
							}
						}
					}
                }
            }
        }


    private void addProbabilitiesForDeleteSubstituent(ArrayList<Mutation> mutationList,
                                                      ArrayList<MutatorSubstituent> substituentList) {
        for (MutatorSubstituent s:substituentList) {
            if (mIsPrimaryAtom[s.coreAtom]) {
            	boolean[] isMemberAtom = new boolean[mMol.getAllAtoms()];
            	mMol.getSubstituent(s.coreAtom, s.firstAtom, isMemberAtom, null, null);
            	boolean selectedAtomFound = false;
            	for (int atom=0; atom<mMol.getAllAtoms(); atom++) {
					if (isMemberAtom[atom] && mMol.isSelectedAtom(atom)) {
						selectedAtomFound = true;
						break;
						}
					}

            	if (!selectedAtomFound) {
					mMol.copyMolecule(mMolCopy);
					mMolCopy.deleteBond(s.bond);

					double f_Educt = Math.max(cMinEductProbability, getFrequency(mMol, s.coreAtom));

					double f_Deleted = 1.0;
					for (int atom:s.atoms)
						f_Deleted *= Math.max(cMinEductProbability, getFrequency(mMol, atom));

					double f_Product = getFrequency(mMolCopy, s.coreAtom);

					double p = f_Product / (f_Educt * Math.pow(f_Deleted, 1.0 / s.atoms.length));
					if (p > 0.0 && isValidStructure(mMolCopy)) {
						if (mBiasProvider != null)
							p *= mBiasProvider.getBiasFactor(mMolCopy);

						mutationList.add(new Mutation(Mutation.MUTATION_DELETE_SUBSTITUENT, s.coreAtom, -1, s.firstAtom, -1, p));
						}
					}
				}
            }
        }


    private void addProbabilitiesForCutOutFragment(ArrayList<Mutation> mutationList) {
        for (int atom=0; atom<mMol.getAtoms(); atom++) {
            if (mIsPrimaryAtom[atom]
             && mMol.getConnAtoms(atom) > 2) {
                int ringBondCount = 0;
                for (int i=0; i<mMol.getConnAtoms(atom); i++)
                    if (mMol.isRingBond(mMol.getConnBond(atom, i)))
                        ringBondCount++;
                if (ringBondCount <= 2) {
                    for (int i=1; i<mMol.getConnAtoms(atom); i++) {
                        int atom1 = mMol.getConnAtom(atom, i);
                        int bond1 = mMol.getConnBond(atom, i);
                        if (mMol.getBondOrder(bond1) != 1)
                            continue;
                        for (int j=0; j<i; j++) {
                            int atom2 = mMol.getConnAtom(atom, j);
                            int bond2 = mMol.getConnBond(atom, j);
                            if (mMol.getBondOrder(bond2) != 1)
                                continue;
                            int coveredRingBondCount = (mMol.isRingBond(bond1) ? 1 : 0)
                                                     + (mMol.isRingBond(bond2) ? 1 : 0);
                            if (coveredRingBondCount == ringBondCount) {
                                mMol.copyMolecule(mMolCopy);
                                mMolCopy.setBondAtom(0, bond1, atom1);
                                mMolCopy.setBondAtom(1, bond1, atom2);
                                mMolCopy.deleteBond(bond2);

                                int[] deletedAtoms = mMolCopy.getFragmentAtoms(atom);

                                boolean selectedAtomFound = false;
                                for (int deletedAtom:deletedAtoms) {
                                	if (mMol.isSelectedAtom(deletedAtom)) {
                                		selectedAtomFound = true;
                                		break;
										}
									}
								if (!selectedAtomFound) {
									int[] atomMap = mMolCopy.deleteAtoms(deletedAtoms);

									double f_Educt = Math.max(cMinEductProbability, getFrequency(mMol, atom1))
												   * Math.max(cMinEductProbability, getFrequency(mMol, atom2));
									double f_Deleted = 1.0;
									for (int deletedAtom:deletedAtoms)
										f_Deleted *= Math.max(cMinEductProbability, getFrequency(mMol, deletedAtom));

									double f_Product = getFrequency(mMolCopy, atomMap[atom1])
													 * getFrequency(mMolCopy, atomMap[atom2]);

									double p = Math.sqrt(f_Product / f_Educt) / Math.pow(f_Deleted, 1.0 / deletedAtoms.length);
									if (p > 0.0 && isValidStructure(mMolCopy)) {
										if (mBiasProvider != null)
											p *= mBiasProvider.getBiasFactor(mMolCopy);

										mutationList.add(new Mutation(Mutation.MUTATION_CUTOUT_SFRAGMENT, atom, -1, atom1, atom2, p));
										}
									}
								}
							}
						}
                    }
                }
            }
        }

	private void addProbabilitiesForInvertParity(ArrayList<Mutation> mutationList) {
		for (int atom=0; atom<mMol.getAtoms(); atom++)
			if (mIsPrimaryAtom[atom] && mMol.isAtomStereoCenter(atom))
				mutationList.add(new Mutation(Mutation.MUTATION_INVERT_PARITY, atom, -1, -1, -1, 1.0));
		for (int bond=0; bond<mMol.getBonds(); bond++)
			if ((mIsPrimaryAtom[mMol.getBondAtom(0, bond)] || mIsPrimaryAtom[mMol.getBondAtom(1, bond)])
			 && mMol.getBondParity(bond) != Molecule.cBondParityNone)
				mutationList.add(new Mutation(Mutation.MUTATION_INVERT_PARITY, -1, bond, -1, -1, 1.0));
		}

    private ArrayList<MutatorSubstituent> createSubstituentList() {
		ArrayList<MutatorSubstituent> substituentList = new ArrayList<MutatorSubstituent>();

		boolean[] isCentralAtom = ScaffoldHelper.findMurckoScaffold(mMol);
		if (isCentralAtom != null) {
			for (int bond=0; bond<mMol.getBonds(); bond++) {
				for (int i=0; i<2; i++) {
					int atom1 = mMol.getBondAtom(i, bond);
					int atom2 = mMol.getBondAtom(1-i, bond);
					if (isCentralAtom[atom1]
					 && !isCentralAtom[atom2]
					 && !mMol.isSelectedAtom(atom1)
					 && !mMol.isSelectedAtom(atom2)) {
						substituentList.add(new MutatorSubstituent(atom1, atom2, bond));
						}
					}
				}
			}

        return substituentList;
	    }


	private double getFrequency(StereoMolecule mol, int atom) {
		if (mAtomTypeList != null) {
			try {
				long atomType = AtomTypeCalculator.getAtomType(mol, atom, AtomTypeCalculator.cPropertiesForMutator);
				return mAtomTypeList.getProbabilityFromType(atomType);
				}
			catch (Exception e) {
				return 0.0f;
				}
			}
		return 1.0;
		}


    private synchronized Mutation selectLikelyMutation(ArrayList<Mutation> mutationList) {
		double probabilitySum = 0.0f;
		for (Mutation m:mutationList)
			probabilitySum += m.mProbability;

		double selector = mRandom.nextDouble() * probabilitySum;

		probabilitySum = 0.0f;
        for (Mutation m:mutationList) {
			probabilitySum += m.mProbability;
			if (selector < probabilitySum) {
				mutationList.remove(m);
			    return m;
			    }
			}
		return null;
		}


    /**
     * Performs the given mutation on the molecule, updates atom coordinates,
     * and updates stereo bonds to reflect lost or new stereo centers.
     * @param mol
     * @param mutation
     */
	public void performMutation(StereoMolecule mol, Mutation mutation) {
        mol.ensureHelperArrays(Molecule.cHelperParities);
		switch (mutation.mMutationType) {
		case Mutation.MUTATION_ADD_ATOM:
			int newAtom = mol.addAtom(mutation.mSpecifier1);
			mol.addBond(mutation.mWhere1, newAtom, mutation.mSpecifier2);
			break;
		case Mutation.MUTATION_INSERT_ATOM:
			int atom1 = mol.getBondAtom(0, mutation.mWhere1);
			int atom2 = mol.getBondAtom(1, mutation.mWhere1);
			mol.deleteBond(mutation.mWhere1);
			newAtom = mol.addAtom(mutation.mSpecifier1);
			mol.addBond(atom1, newAtom, Molecule.cBondTypeSingle);
			mol.addBond(atom2, newAtom, Molecule.cBondTypeSingle);
			break;
		case Mutation.MUTATION_CHANGE_ATOM:
			mol.setAtomicNo(mutation.mWhere1, mutation.mSpecifier1);
			break;
		case Mutation.MUTATION_DELETE_ATOM:
			mol.deleteAtom(mutation.mWhere1);	
			break;
		case Mutation.MUTATION_CUTOUT_ATOM:
			atom1 = mol.getConnAtom(mutation.mWhere1, 0);
			atom2 = mol.getConnAtom(mutation.mWhere1, 1);
			mol.addBond(atom1, atom2, Molecule.cBondTypeSingle);
			mol.deleteAtom(mutation.mWhere1);
			break;
		case Mutation.MUTATION_CLOSE_RING:
			mol.addBond(mutation.mWhere1, mutation.mWhere2, mutation.mSpecifier1);
			break;
		case Mutation.MUTATION_CLOSE_RING_AND_AROMATIZE:
			mol.addBond(mutation.mWhere1, mutation.mWhere2, mutation.mSpecifier1);
			aromatizeRing(mol, mutation.mAtomList, mutation.mSpecifier1);
			break;
		case Mutation.MUTATION_TOGGLE_AMID_SULFONAMID:
			if (mol.getAtomicNo(mutation.mWhere1) == 6) {
				mol.setAtomicNo(mutation.mWhere1, 16);
				mol.addBond(mutation.mWhere1, mol.addAtom(8), Molecule.cBondTypeDouble);
				}
			else {
				mol.setAtomicNo(mutation.mWhere1, 6);
				mol.deleteAtom(mutation.mWhere2);
				}
			break;
		case Mutation.MUTATION_CHANGE_BOND:
			mol.setBondType(mutation.mWhere1, mutation.mSpecifier1);
			break;
		case Mutation.MUTATION_DELETE_BOND:
			mol.deleteBond(mutation.mWhere1);
			break;
		case Mutation.MUTATION_CHANGE_RING:
			changeAromaticity(mol, mutation.mWhere1, mutation.mSpecifier1);
			break;
		case Mutation.MUTATION_MIGRATE:
			for (int i=0; i<2; i++)
				if (mol.getBondAtom(i, mutation.mWhere1) == mutation.mSpecifier1)
					mol.setBondAtom(i, mutation.mWhere1, mutation.mSpecifier2);
			break;
		case Mutation.MUTATION_DELETE_SUBSTITUENT:
		    boolean[] atomMask = new boolean[mol.getAtoms()];
		    mol.getSubstituent(mutation.mWhere1, mutation.mSpecifier1, atomMask, null, null);
		    mol.deleteAtoms(atomMask);
		    break;
        case Mutation.MUTATION_CUTOUT_SFRAGMENT:
            int rootAtom = mutation.mWhere1;
            atom1 = mutation.mSpecifier1;
            atom2 = mutation.mSpecifier2;
            int bond1 = -1;
            int bond2 = -1;
            for (int i=0; i<mol.getConnAtoms(rootAtom); i++) {
                if (mol.getConnAtom(rootAtom, i) == atom1)
                    bond1 = mol.getConnBond(rootAtom, i);
                else if (mol.getConnAtom(rootAtom, i) == atom2)
                    bond2 = mol.getConnBond(rootAtom, i);
                }
            if (bond1 != -1 && bond2 != -1) {
                mol.deleteBond(bond1);
                mol.deleteBond(bond2);
                int[] atomMap = mol.deleteAtoms(mol.getFragmentAtoms(rootAtom));
                mol.addBond(atomMap[atom1], atomMap[atom2], Molecule.cBondTypeSingle);
                }
            break;
		case Mutation.MUTATION_INVERT_PARITY:
			if (mutation.mWhere1 != -1)
				mol.setAtomParity(mutation.mWhere1, mol.getAtomParity(mutation.mWhere1) != Molecule.cAtomParity1 ?
						Molecule.cAtomParity1 : Molecule.cAtomParity2, false);
			if (mutation.mWhere2 != -1)
				mol.setBondParity(mutation.mWhere2, mol.getBondParity(mutation.mWhere2) != Molecule.cBondParityEor1 ?
						Molecule.cBondParityEor1 : Molecule.cBondParityZor2, false);
			break;
			}

		repairCharges(mol);

		// Most of the parity flags are still valid; the CoordinateInventor creates up/down bonds from them
		mol.setParitiesValid(0);
		new CoordinateInventor().invent(mol);

		// we need to invalidate to detect all parities correctly now
		mol.invalidateHelperArrays(Molecule.cHelperAll);

		repairStereoChemistry(mol);	// assign random parities to new stereo centers, and change up/down accordingly
	    }


	private int getBondTypeFromOrder(int bondOrder) {
		switch (bondOrder) {
		case 1:
			return Molecule.cBondTypeSingle;
		case 2:
			return Molecule.cBondTypeDouble;
		case 3:
			return Molecule.cBondTypeTriple;
			}
		return 0;
		}


	private boolean qualifiesForRing(int[] graphBond, int[] graphLevel,
									 int[] graphParent, int atom1, int atom2) {
		int ringSize = graphLevel[atom2];
		int firstConsecutiveRingAtom = 0;
		int firstMaxConsecutiveRingAtom = 0;
		int lastMaxConsecutiveRingAtom = 0;
		int consecutiveRingBonds = 0;
		int maxConsecutiveRingBonds = 0;
		int currentAtom = atom2;
		while (currentAtom != atom1) {
			if (mMol.getAtomPi(currentAtom) == 2) {
				if (mMol.getConnBondOrder(currentAtom, 0) == 2) {	// allene
					if (ringSize < 9)
						return false;
					}
				else {	// alkyne
					if (ringSize < 10)
						return false;
					}
				}

			if (mMol.isRingBond(graphBond[currentAtom])) {
				if (consecutiveRingBonds == 0)
					firstConsecutiveRingAtom = currentAtom;
				consecutiveRingBonds++;
				}
			else {
				if (maxConsecutiveRingBonds < consecutiveRingBonds) {
					maxConsecutiveRingBonds = consecutiveRingBonds;
					firstMaxConsecutiveRingAtom = firstConsecutiveRingAtom;
					lastMaxConsecutiveRingAtom = currentAtom;
					}
				consecutiveRingBonds = 0;
				}

			currentAtom = graphParent[currentAtom];
			}

		if (maxConsecutiveRingBonds < consecutiveRingBonds) {
			maxConsecutiveRingBonds = consecutiveRingBonds;
			firstMaxConsecutiveRingAtom = firstConsecutiveRingAtom;
			lastMaxConsecutiveRingAtom = currentAtom;
			}

		int penalty = maxConsecutiveRingBonds;
		if (maxConsecutiveRingBonds > 0
		 && ((mMol.isAromaticAtom(firstMaxConsecutiveRingAtom)
		   || mMol.isAromaticAtom(lastMaxConsecutiveRingAtom))
		  || (mMol.getAtomPi(firstMaxConsecutiveRingAtom) > 0
		   && mMol.getAtomPi(lastMaxConsecutiveRingAtom) > 0)))
			penalty++;

		if (ringSize - penalty < 3)
			return false;
		
		return true;
		}


	/**
	 * Runs a rule based check for bicyclic systems with unacceptable strain.
	 * Checks any bridge between two different atoms of one ring
	 * connected to different atoms of one ring.
	 * @param mol
	 * @return
	 */
	private boolean isValidStructure(StereoMolecule mol) {
	    // The largest bridge length (bonds) between two direct ring neighbours
	    // that would still cause an unacceptable ring strain.
	    // (this is the smallest allowed chain (bond count) from ring atom to ring atom minus 3!!!)
	    // First index is ring size, 2nd index is number of bonds between attachment points
	    // We distinguish 3 types of rings:
	    // - aromatic or both attachment points sp2
	    // - one attachment point sp2
	    // - both attachment points sp3
	    final int[][] MAX_FORBIDDEN_BRIDGE_LENGTH_AROM = { null, null, null,
	            { -1, 2 },
	            { -1, 2, 6 },
	            { -1, 1, 5 },
	            { -1, 1, 4, 6 },
	            { -1, 1, 4, 6 } };
        final int[][] MAX_FORBIDDEN_BRIDGE_LENGTH_PI = { null, null, null,
                { -1, 1 },
                { -1, 1, 5 },
                { -1, 1, 4 },
                { -1, 0, 3, 4 },
                { -1, 0, 2, 3 } };
        final int[][] MAX_FORBIDDEN_BRIDGE_LENGTH_ALIPH = { null, null, null,
                { -1, 0 },
                { -1, -1, 3 },
                { -1, -1, 0 },
                { -1, -1, 0, -1 },
                { -1, -1, 0, -1 } };

	    mol.ensureHelperArrays(Molecule.cHelperRings);

	    RingCollection ringSet = mol.getRingSet();
	    boolean[] neglectAtom = new boolean[mol.getAtoms()];
	    for (int ring=0; ring<ringSet.getSize() && (ringSet.getRingSize(ring)<8); ring++) {
	        int[] ringAtom = ringSet.getRingAtoms(ring);
	        for (int atom:ringAtom)
	            neglectAtom[atom] = true;

	        int ringSize = ringSet.getRingSize(ring);
	        for (int i=1; i<ringSize; i++) {
	            if (mol.getConnAtoms(ringAtom[i]) > 2) {
	                for (int ci=0; ci<mol.getConnAtoms(ringAtom[i]); ci++) {
	                    int atom1 = mol.getConnAtom(ringAtom[i], ci);
	                    if (mol.isRingAtom(atom1) && !ringSet.isAtomMember(ring, atom1)) {
            	            for (int j=0; j<i; j++) {
            	                if (mol.getConnAtoms(ringAtom[j]) > 2) {
            	                    for (int cj=0; cj<mol.getConnAtoms(ringAtom[j]); cj++) {
            	                        int atom2 = mol.getConnAtom(ringAtom[j], cj);
            	                        if (mol.isRingAtom(atom2) && !ringSet.isAtomMember(ring, atom2)) {
            	                            int ringDif = Math.min(i-j, ringSize-i+j);
            	                            int pi1 = mol.getAtomPi(ringAtom[i]);
                                            int pi2 = mol.getAtomPi(ringAtom[j]);
            	                            int maxForbiddenLength =
            	                                (ringSet.isAromatic(ring) || (pi1!=0 && pi2!=0)) ?
            	                                    MAX_FORBIDDEN_BRIDGE_LENGTH_AROM[ringSize][ringDif]
                                              : (pi1!=0 || pi2!=0) ?
                                                    MAX_FORBIDDEN_BRIDGE_LENGTH_PI[ringSize][ringDif]
                                                  : MAX_FORBIDDEN_BRIDGE_LENGTH_ALIPH[ringSize][ringDif];
            	                            if (maxForbiddenLength != -1 && mol.getPathLength(atom1, atom2, maxForbiddenLength, neglectAtom) != -1)
            	                                return false;
            	                            }
            	                        }
            	                    }
            	                }
	                        }
	                    }
	                }
	            }

	        for (int atom:ringAtom)
                neglectAtom[atom] = false;
	        }
	    
	    return true;
	    }


	private boolean hasExocyclicPiBond(int[] ringAtom, int[] ringBond) {
		for (int i=0; i<ringAtom.length; i++) {
			for (int j=0; j<mMol.getConnAtoms(ringAtom[i]); j++) {
				int connBond = mMol.getConnBond(ringAtom[i], j);
				boolean bondIsInRing = false;
				for (int k=0; k<ringBond.length; k++) {
					if (connBond == ringBond[k]) {
						bondIsInRing = true;
						break;
						}
					}
				if (!bondIsInRing
				 && (mMol.isAromaticBond(connBond)
				  || mMol.getBondOrder(connBond) > 1))
					return true;
				}
			}
		return false;
		}


	private boolean changeAromaticity(StereoMolecule mol, int ring, int heteroPosition) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
		RingCollection ringSet = mol.getRingSet();
		int ringSize = ringSet.getRingSize(ring);
		int ringAtom[] = ringSet.getRingAtoms(ring);
		int ringBond[] = ringSet.getRingBonds(ring);
		if (ringSet.isAromatic(ring)) {
			for (int i=0; i<ringSize; i++)
				mol.setBondType(ringBond[i], Molecule.cBondTypeSingle);
			return true;
			}
		else {
			for (int i=0; i<ringSize; i++)
				mol.setBondType(ringBond[i], Molecule.cBondTypeSingle);

			if (ringSize == 5) {
				for (int i=0; i<ringSize; i++) {
					if (heteroPosition != i) {
						if (mol.getFreeValence(ringAtom[i]) < 1)
							return false;
						if (mol.getAtomicNo(ringAtom[i]) != 6
						 && mol.getAtomicNo(ringAtom[i]) != 7)
							return false;
						}
					}
				int doubleBondPosition1 = heteroPosition + 1;
				int doubleBondPosition2 = heteroPosition + 3;
				if (doubleBondPosition1 > 4)
					doubleBondPosition1 -= 5;
				if (doubleBondPosition2 > 4)
					doubleBondPosition2 -= 5;
				mol.setBondType(doubleBondPosition1, Molecule.cBondTypeDouble);
				mol.setBondType(doubleBondPosition2, Molecule.cBondTypeDouble);

				return true;
				}

			if (ringSize == 6) {
				for (int i=0; i<ringSize; i++) {
					if (mol.getFreeValence(ringAtom[i]) < 1)
						return false;
					if (mol.getAtomicNo(ringAtom[i]) != 6
					 && mol.getAtomicNo(ringAtom[i]) != 7)
						return false;
					}

				for (int i=0; i<ringSize; i+=2)
					mol.setBondType(ringBond[i], Molecule.cBondTypeDouble);

				return true;
				}
			}
		return false;
		}

	private void repairCharges(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperRings);
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			// make sure that quarternary nitrogen is charged
			if (mol.getAtomicNo(atom) == 7
			 && (mol.getConnAtoms(atom) + mol.getAtomPi(atom)) == 4) {
				mol.setAtomCharge(atom, 1);
				for (int i=0; i<mol.getConnAtoms(atom); i++) {
					int connAtom = mol.getConnAtom(atom, i);
					if ((mol.getAtomicNo(connAtom) == 7 && (mol.getConnAtoms(connAtom) + mol.getAtomPi(connAtom)) < 3)
					 || (mol.getAtomicNo(connAtom) == 8 && (mol.getConnAtoms(connAtom) + mol.getAtomPi(connAtom)) < 2)) {
						mol.setAtomCharge(connAtom, -1);
						}
					}
				}
			// discharge carbon atoms
			if (mol.getAtomicNo(atom) == 6
			 && mol.getAtomCharge(atom) != 0) {
				mol.setAtomCharge(atom, 0);
				}
			}
		// remove isolated charges where possible
		for (int atom=0; atom<mol.getAtoms(); atom++) {
			if (mol.getAtomCharge(atom) != 0) {
				boolean chargedNeighborFound = false;
				for (int i=0; i<mol.getConnAtoms(atom); i++) {
					if (mol.getAtomCharge(mol.getConnAtom(atom, i)) != 0) {
						chargedNeighborFound = true;
						break;
						}
					}
				if (!chargedNeighborFound) {
					int valence = mol.getConnAtoms(atom) + mol.getAtomPi(atom);
					int maxValence = mol.getMaxValenceUncharged(atom);
					if (valence <= maxValence)
						mol.setAtomCharge(atom, 0);
					}
				}
			}
		}

	private void repairStereoChemistry(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperParities);	// detect over/under-specified stereo information
		for (int bond=0; bond<mol.getAllBonds(); bond++)
			if (mol.isStereoBond(bond))
				mol.setBondType(bond, Molecule.cBondTypeSingle);

        for (int atom=0; atom<mol.getAtoms(); atom++) {
        	int parity = mol.getAtomParity(atom);
        	if (parity == Molecule.cAtomParityUnknown) {
		        parity = (mRandom.nextDouble()<0.5) ? Molecule.cAtomParity1 : Molecule.cAtomParity2;
		        boolean isPseudo = mol.isAtomParityPseudo(atom);
		        mol.setAtomParity(atom, parity, isPseudo);
		        mol.setAtomESR(atom, Molecule.cESRTypeAbs, 0);
	            }
        	if (parity != Molecule.cAtomParityNone)
				mol.setStereoBondFromAtomParity(atom);
            }
        for (int bond=0; bond<mol.getBonds(); bond++) {
        	if (mol.isBINAPChiralityBond(bond)) {
	            switch (mol.getBondParity(bond)) {
	            case Molecule.cBondParityUnknown:
	                int parity = (mRandom.nextDouble() < 0.5) ? Molecule.cBondParityEor1 : Molecule.cBondParityZor2;
	                boolean isPseudo = mol.isBondParityPseudo(bond);
	                mol.setBondParity(bond, parity, isPseudo);
	            case Molecule.cBondParityEor1:
	            case Molecule.cBondParityZor2:
	                mol.setStereoBondFromBondParity(bond);
	                break;
	                }
        		}
            }
        }

    private boolean aromatizeRing(StereoMolecule mol, int[] ringAtom, int ringSize) {
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);

		int leakIndex = -1;
		for (int i=0; i<ringSize; i++) {
			int atomicNo = mol.getAtomicNo(ringAtom[i]);
			if (atomicNo == 8 || atomicNo == 16) {
				if (ringSize == 6)
					return false;
				if (leakIndex != -1)
					return false;
				leakIndex = i;
				}
			else if (atomicNo == 7) {
				if (mol.getFreeValence(ringAtom[i]) == 0 && mol.getAtomPi(ringAtom[i]) == 0) {
					if (leakIndex != -1)
						return false;
					leakIndex = i;
					}
				}
			else if (atomicNo != 6) {
				return false;
				}
			}

		if (leakIndex == -1 && ringSize == 5) {
			for (int i=0; i<ringSize; i++) {
				if (mol.getAtomicNo(ringAtom[i]) == 7) {
					leakIndex = i;
					break;
					}
				}
			}

		int[] ringBond = new int[ringSize];
		for (int i=0; i<ringSize; i++)
			ringBond[i] = mol.getBond(ringAtom[i], ringAtom[i != 0 ? i-1 : ringSize-1]);

		// convert exocyclic double bonds into single bonds
		for (int i=0; i<ringSize; i++) {
			int atom1 = ringAtom[i==0 ? ringSize-1 : i-1];
			int atom2 = ringAtom[i==ringSize-1 ? 0 : i+1];
			for (int j=0; j<mol.getConnAtoms(ringAtom[i]); j++) {
				int connAtom = mol.getConnAtom(ringAtom[i], j);
				if (connAtom != atom1 && connAtom != atom2 && mol.getConnBondOrder(ringAtom[i], j) != 1)
					mol.setBondType(mol.getConnBond(ringAtom[i], j), Molecule.cBondTypeSingle);
				}
			}

		for (int i=0; i<ringSize; i++)
			mol.setBondType(ringBond[i], Molecule.cBondTypeSingle);

		int doubleBondIndex = (leakIndex == -1) ? 0 : leakIndex+2;
		int doubleBondCount = (ringSize == 5) ? 2 : 3;
		for (int i=0; i<doubleBondCount; i++) {
			if (doubleBondIndex >= ringSize)
				doubleBondIndex -= ringSize;
			mol.setBondType(ringBond[doubleBondIndex], Molecule.cBondTypeDouble);
			doubleBondIndex += 2;
			}

		return true;
		}

    public void printMutationList(ArrayList<Mutation> mutationList, boolean sortByPriority) {
		Mutation[] mutation = mutationList.toArray(new Mutation[0]);

		if (sortByPriority)
			Arrays.sort(mutation, new Comparator<Mutation>() {
				@Override public int compare(Mutation m1, Mutation m2) {
					return m1.mProbability > m2.mProbability ? -1 : m1.mProbability < m2.mProbability ? 1 : 0;
				}
			} );

		System.out.println("Mutation list ("+mutationList.size()+" mutations):");
		for (Mutation m:mutation)
			System.out.println(m.toString());
		}

    class MutatorSubstituent {
        public int coreAtom;
        public int firstAtom;
        public int bond;
        public int[] atoms;
//        public int size;

        public MutatorSubstituent(int coreAtom, int firstAtom, int bond) {
            this.coreAtom = coreAtom;
            this.firstAtom = firstAtom;
            this.bond = bond;
            boolean[] isMember = new boolean[mMol.getAllAtoms()];
            int count = mMol.getSubstituent(coreAtom, firstAtom, isMember, null, null);
            atoms = new int[count];
            for (int atom=0, i=0; atom<mMol.getAllAtoms(); atom++)
            	if (isMember[atom])
            		atoms[i++] = atom;

//            this.size = mMol.getSubstituentSize(coreAtom, firstAtom);
// in a first approach size is not considered
            }
        }
    }
