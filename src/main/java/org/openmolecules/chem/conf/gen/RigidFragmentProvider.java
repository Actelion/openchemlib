/*
 * Copyright 2013-2020 Thomas Sander, openmolecules.org
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Thomas Sander
 */

package org.openmolecules.chem.conf.gen;

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Coordinates;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.forcefield.mmff.BadAtomTypeException;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import org.openmolecules.chem.conf.so.ConformationSelfOrganizer;
import org.openmolecules.chem.conf.so.SelfOrganizedConformer;

import java.util.ArrayList;

/**
 * An instance of this class is used by any ConformerGenerator to hand out one or multiple 3D-coordinate sets
 * for rigid fragments within a molecules. RigidFragments are those substructures obtained, when breaking a
 * molecule apart at its non-ring single bonds, which are called 'rotatable bonds'.
 * A RigidFragment consists of core atoms, extended atoms and outer shell atoms.
 * Core atoms are all atoms making up the original substructure. Extended atoms are all atoms of the original
 * molecule, which directly connect to any core atom. Extended atoms are important, because they define the exit
 * vectors of the core atom fragment. Outer shell atoms are another layer of atoms connected the extended atoms.
 * They are needed for the self-organizer to create context depending 3D-coordinates that considers steric and
 * atom type dependent geometries.<br>
 * RigidFragments are fragments rather than molecules. There are no implicit hydrogen atoms and some atoms
 * have unoccupied valences. Nitrogen atoms may carry the flat-nitrogen query feature, which is considered
 * by the self organizer when creating 3D-coordinates.<br>
 * <b>Important:</b> Stereo centers in the molecule may not be a stereo center in the fragment anymore.
 * Nevertheless, the self-organized coordinates must reflect the correct configuration of the originating molecule.
 * Therefore, atom parities are copied from molecule to fragment and defined valid, which is possible, because
 * parities are based on atom indexes and the atom order is kept in-tact, when assembling the fragment.<br>
 * <b>Caching:</b> When generating conformers of many molecules, then many RigidFragments occurr repeatedly.
 * Therefore, caching of fragment's coordinates speeds up conformer generation significantly. However, we need
 * to consider various issues:<br>
 * - As key for the cache we use the fragment's idcode (all atoms). To not loose hydrogen atoms, we convert all
 * plain hydrogen atoms into deuterium.<br>
 * - In the cache we just store 3D-coordinate sets and their likelyhoods based on self-organizer scores or MMFF-energies<br>
 * - When locating the same fragment in a different molecule, then the atom order may be different. Thus, we need to normalize.
 * For this purpose we use the graphIndex of the Canonizer used to create the fragment's idcode and store coordinates in
 * canonical atom order.<br>
 * - For molecule stereo centers, which are gone in the fragment, copying of molecule parities to the fragment ensures
 * proper 3D-coordinates. For the Canonizer graphIndex to reflect the original configuration, we need to make sure,
 * that up/down-bonds are copied (parities won't do), which are now overspecifying the non-stereo center.
 * And we use the Canonizer mode TIE_BREAK_FREE_VALENCE_ATOMS to distinguish symmetrical fragment atoms if they have
 * free valences and, thus, could be differently substituted in molecule matches. This way we locate all potential
 * stereo centers in fragments. If a given molecule does not specify a parity for any of the potential stereo
 * centers, then this fragment is not cached. Otherwise a later hit with a defined stereo center at that position
 * might get coordinates for the wrong stereo configuration.<br>
 * - If a fragment contains stereo centers, then only one of the two possible enantiomers is cached. The other one is
 * constructed by z-coordinate inversion.<br>
 * - If optimizeRigidFragments is true, then the MMFF94s+ forcefield is used to minimize fragments before caching/using
 * them. A different force field may be used by overriding RigidFragmentProvider and all of its forceField...
 * methods and passing an overridden instance to the constructor(s) of the ConformerGenerator to be used.
 */
public class RigidFragmentProvider {
	private static final int MAX_CONFORMERS = 16;
	private static final int MAX_ATOMS_FOR_CACHING = 32;

	// Random seed for initializing the SelfOrganizer.
	private final long mRandomSeed;
	private final boolean mOptimizeFragments;
	private RigidFragmentCache mCache;
	private ThreadMaster mThreadMaster;
	private long mStopMillis;

	public RigidFragmentProvider(long randomSeed, RigidFragmentCache cache, boolean optimizeRigidFragments) {
		mRandomSeed = randomSeed;
		mCache = cache;
		mOptimizeFragments = optimizeRigidFragments;
		if (optimizeRigidFragments)
			forceFieldInitialize();
		}

	/**
	 * If the conformer generation must be stopped from outside, for instance because of user
	 * intervention or because of a defined timeout, then provide a ThreadMaster with this method.
	 * @param tm
	 */
	public void setThreadMaster(ThreadMaster tm) {
		mThreadMaster = tm;
		}

	/**
	 * @param millis time point as system millis after which to gracefully stop self organization even if not successful
	 */
	public void setStopTime(long millis) {
		mStopMillis = millis;
		}

	public void setCache(RigidFragmentCache cache) {
		mCache = cache;
		}

	public RigidFragment createFragment(StereoMolecule mol, int[] fragmentNo, int fragmentIndex) {
		int coreAtomCount = 0;
		int atomCount = 0;
		int extendedAtomCount = 0;

		// mark all atoms with specified fragmentNo and two layers around it
		boolean[] includeAtom = new boolean[mol.getAllAtoms()];
		boolean[] isOuterShellAtom = new boolean[mol.getAllAtoms()];
		boolean[] isCoreFragment = new boolean[mol.getAllAtoms()];

		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			if (fragmentNo[atom] == fragmentIndex) {
				isCoreFragment[atom] = true;
				includeAtom[atom] = true;
				atomCount++;
				coreAtomCount++;
				}
			}
		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			if (isCoreFragment[atom]) {
				for (int i=0; i<mol.getConnAtoms(atom); i++) {
					int connAtom = mol.getConnAtom(atom, i);
					if (!includeAtom[connAtom]) {
						includeAtom[connAtom] = true;
						atomCount++;
						extendedAtomCount++;
						}
					}
				}
			}
		for (int atom=0; atom<mol.getAllAtoms(); atom++) {
			if (includeAtom[atom] && !isCoreFragment[atom] && !isOuterShellAtom[atom]) {
				for (int i=0; i<mol.getAllConnAtoms(atom); i++) {
					int connAtom = mol.getConnAtom(atom, i);
					if (!includeAtom[connAtom]) {
						isOuterShellAtom[connAtom] = true;
						includeAtom[connAtom] = true;
						atomCount++;
						}
					}
				}
			}

		int bondCount = 0;
		for (int bond = 0; bond < mol.getAllBonds(); bond++)
			if (includeAtom[mol.getBondAtom(0, bond)]
			 && includeAtom[mol.getBondAtom(1, bond)])
				bondCount++;

		StereoMolecule fragment = new StereoMolecule(atomCount, bondCount);
		mol.copyMoleculeByAtoms(fragment, includeAtom, false, null);

		fragment.setFragment(true); // Can encode as fragment, because H-atoms are converted to deuterium

		int[] coreToFragmentAtom = new int[coreAtomCount];
		int[] fragmentToOriginalAtom = new int[atomCount];
		int[] extendedToFragmentAtom = new int[coreAtomCount + extendedAtomCount];
		int[] originalToExtendedAtom = new int[mol.getAllAtoms()];

		int coreAtom = 0;
		int fragmentAtom = 0;
		int extendedAtom = 0;
		for (int atom = 0; atom < mol.getAllAtoms(); atom++) {
			if (includeAtom[atom]) {
				if (mol.isFlatNitrogen(atom))
				    fragment.setAtomQueryFeature(fragmentAtom, Molecule.cAtomQFFlatNitrogen, true);

				// The ConformationSelfOrganizer uses all torsions of qualifying rotatable bonds
				// to check, whether a conformer is a new one. Bonds with marked atoms are excluded.
				// We add two shells of atoms to simulate real molecule steric effects, but we don't
				// want to have multiple conformers of rigid fragments that only differ in these
				// unused substituents.
				if (!isCoreFragment[atom])
					fragment.setAtomMarker(fragmentAtom, true);

				if (isCoreFragment[atom] || !isOuterShellAtom[atom]) {
					extendedToFragmentAtom[extendedAtom] = fragmentAtom;
					originalToExtendedAtom[atom] = extendedAtom;
					extendedAtom++;
				}

				if (isCoreFragment[atom]) {
					coreToFragmentAtom[coreAtom] = fragmentAtom;
					coreAtom++;
				}

				fragmentToOriginalAtom[fragmentAtom] = atom;

				// convert all plain hydrogen to deuterium that we don't lose them in the idcode
				if (fragment.getAtomicNo(fragmentAtom) == 1)
					fragment.setAtomMass(fragmentAtom, 2);

				fragmentAtom++;
			}
		}

		ArrayList<SelfOrganizedConformer> conformerList = null;
		double[] likelihood = null;
		Canonizer canonizer = null;
		String key = null;
		boolean invertedEnantiomer = false;

		boolean useCache = (mCache != null && mCache.canAddEntry() && atomCount <= MAX_ATOMS_FOR_CACHING);

		// Generate stereo parities for all potential stereo features in the fragment.
		// If one or more potential stereo features are unknown, then the fragment doesn't qualify to be cached.
		if (useCache) {
			// By distinguishing equal ranking atoms, if they have free valences, we detect all possible stereo features
			canonizer = new Canonizer(fragment, Canonizer.TIE_BREAK_FREE_VALENCE_ATOMS);

			// we don't cache fragments with unspecified stereo configurations
			for (int atom=0; atom<fragment.getAtoms(); atom++) {
				if (canonizer.getTHParity(atom) == Molecule.cAtomParityUnknown) {
					useCache = false;
					break;
					}
				}
			for (int bond=0; bond<fragment.getBonds(); bond++) {
				if (canonizer.getEZParity(bond) == Molecule.cBondParityUnknown) {
					useCache = false;
					break;
					}
				}

			// If the fragment qualifies for caching, then we use Canonizer parities, which are consistent with canonical atom numbering.
			// Otherwise, we keep and use the original parities for coordinate generation.
			if (useCache)
				canonizer.setParities();
			}

		if (mCache != null && !useCache)
			mCache.increaseNonCachableCount();

		// No matter, whether parities were copied from the original molecule, or whether we use freshly calculated parities,
		// we need to define them to be valid for the coordinate self-organization.
		fragment.setParitiesValid(0);

		// Check, whether we have the fragment already in the cache.
		// If yes, then map coordinates from canonical order and mirror coordinates, if needed.
		// Coordinates are store normalized to one enantiomer
		if (useCache) {
			invertedEnantiomer = canonizer.normalizeEnantiomer();
			key = canonizer.getIDCode();

			RigidFragmentCache.CacheEntry cacheEntry = mCache.get(key);
			if (cacheEntry != null) {
				// convert from canonical coordinates back to fragment
				int[] graphIndex = canonizer.getGraphIndexes();
				conformerList = new ArrayList<>();
				for (Coordinates[] coords:cacheEntry.coordinates) {
					for (int j = 0; j<fragment.getAllAtoms(); j++) {
						Coordinates c = coords[graphIndex[j]];
						fragment.setAtomX(j, c.x);
						fragment.setAtomY(j, c.y);
						fragment.setAtomZ(j, invertedEnantiomer ? -c.z : c.z);
						}
					conformerList.add(new SelfOrganizedConformer(fragment));
					}

				likelihood = cacheEntry.likelihood;
				}
			}

		if (conformerList == null) {
			if (mOptimizeFragments && !forceFieldAllowsOpenValences()) {
				fragment.setFragment(false);    // to allow conversion of implicit to explicit hydrogen
				ConformerGenerator.addHydrogenAtoms(fragment);
			}

			ConformationSelfOrganizer selfOrganizer = new ConformationSelfOrganizer(fragment, true);
			selfOrganizer.setThreadMaster(mThreadMaster);
			selfOrganizer.setStopTime(mStopMillis);
			selfOrganizer.initializeConformers(mRandomSeed, MAX_CONFORMERS);

			// Generate multiple low constraint conformers
			conformerList = new ArrayList<>();
			SelfOrganizedConformer bestConformer = selfOrganizer.getNextConformer();
			conformerList.add(bestConformer);
			SelfOrganizedConformer conformer = selfOrganizer.getNextConformer();
			while (conformer != null) {
				conformerList.add(conformer);
				conformer = selfOrganizer.getNextConformer();
				}

			// Calculate fraction values of the population from strain values, which somewhat resemble energies in kcal/mol
			double ENERGY_FOR_FACTOR_10 = 1.36; // The ConformerSelfOrganizer and MMFF use kcal/mol; 1.36 kcal/mol is factor 10

			double minStrain = Double.MAX_VALUE;
			for(Conformer conf:conformerList)
				minStrain = Math.min(minStrain, ((SelfOrganizedConformer)conf).getTotalStrain());

			// Strain values resemble energies in kcal/mol, but not as reliable. Therefore we are less strict and allow factor 1000
			double strainLimit = minStrain + 3.0 * ENERGY_FOR_FACTOR_10;
			for (int i=conformerList.size()-1; i>=0; i--)
				if (conformerList.get(i).getTotalStrain()>strainLimit)
					conformerList.remove(i);

			likelihood = new double[conformerList.size()];
			double likelihoodSum = 0;
			int index = 0;
			for(int i=0; i<conformerList.size(); i++) {
				SelfOrganizedConformer conf = conformerList.get(i);
				likelihood[i] = Math.pow(10, (minStrain - conf.getTotalStrain()) / ENERGY_FOR_FACTOR_10);
				likelihoodSum += likelihood[i];
				}
			for (int i=0; i<conformerList.size(); i++)
				likelihood[i] /= likelihoodSum;

			if(mOptimizeFragments) {
				int validEnergyCount = 0;
				double minEnergy = Double.MAX_VALUE;
				for(Conformer conf:conformerList) {
					double energy = forceFieldMinimize(conf.toMolecule());
					conf.setEnergy(energy);
					if (!Double.isNaN(energy)) {
						minEnergy = Math.min(minEnergy, energy);
						validEnergyCount++;
						}
					conf.copyFrom(fragment);
					}

				double energyLimit = minEnergy + 2.0 * ENERGY_FOR_FACTOR_10;    // population of less than 1% of best conformer
				for(Conformer conf:conformerList) {
					if (!Double.isNaN(conf.getEnergy()) && conf.getEnergy()>energyLimit) {
						conf.setEnergy(Double.NaN);
						validEnergyCount--;
						}
					}

				// If we have no valid MMFF energy values, we keep the likelihoods from the self organizer, otherwise...
				if (validEnergyCount != 0) {
					double[] population = new double[validEnergyCount];
					double populationSum = 0;
					index = 0;
					for(int i=conformerList.size()-1; i>=0; i--) {
						Conformer conf = conformerList.get(i);
						if (Double.isNaN(conf.getEnergy()))
							conformerList.remove(i);
						else {
							population[index] = Math.pow(10, (minEnergy - conf.getEnergy()) / ENERGY_FOR_FACTOR_10);
							populationSum += population[index];
							index++;
							}
						}

					likelihood = new double[validEnergyCount];
					for (int i=0; i<validEnergyCount; i++)
						likelihood[i] = population[i] / populationSum;
					}
				}

			if (useCache) {
				int[] graphIndex = canonizer.getGraphIndexes();
				Coordinates[][] coords = new Coordinates[conformerList.size()][graphIndex.length];
				for (int i=0; i<coords.length; i++) {
					for (int j=0; j<graphIndex.length; j++) {
						Coordinates xyz = conformerList.get(i).getCoordinates(j);
						coords[i][graphIndex[j]] = new Coordinates(xyz.x, xyz.y, invertedEnantiomer ? -xyz.z : xyz.z);
						}
					}

				mCache.put(key, new RigidFragmentCache.CacheEntry(coords, likelihood));
				}
			}

		return new RigidFragment(coreAtomCount, coreToFragmentAtom, fragmentToOriginalAtom,
				extendedToFragmentAtom, originalToExtendedAtom, conformerList.toArray(new SelfOrganizedConformer[0]), likelihood);
	}

	/**
	 * For using a different forcefield you may override all three forceField... methods.
	 */
	public void forceFieldInitialize() {
		ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
	}

	/**
	 * For using a different forcefield you may override all three forceField... methods.
	 */
	public boolean forceFieldAllowsOpenValences() {
		return true;
	}

	/**
	 * For using a different forcefield you may override all three forceField... methods.
	 */
	public double forceFieldMinimize(StereoMolecule mol) {
		try {
			ForceFieldMMFF94 ff = new ForceFieldMMFF94(mol, ForceFieldMMFF94.MMFF94SPLUS);
			ff.minimise();
			return ff.getTotalEnergy();
			}
		catch (BadAtomTypeException bate) {
			return Double.NaN;
		}
		catch(Exception ex) {
			ex.printStackTrace();
			return Double.NaN;
		}
	}

	/**
	 * @param conformer
	 * @param energy_out [0] is the final energy, [1] is the starting energy
	 * @return
	 *
	public static Conformer minimizeConformer( Conformer conformer , double[] energy_out ) {
		StereoMolecule mol = conformer.toMolecule();

		// make fragment to omit missing hydrogens
		mol.setFragment(true);
		int n_atoms = conformer.getSize();//mol.getAtoms();

		try {
			ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
			ForceFieldMMFF94 ff = new ForceFieldMMFF94(mol, ForceFieldMMFF94.MMFF94SPLUS);

			if(energy_out!=null) {
				if (energy_out.length > 1) {
					energy_out[1] = ff.getTotalEnergy();
				}
			}

			ff.minimise();
			if(energy_out!=null) {
				if (energy_out.length > 1) {
					energy_out[0] = ff.getTotalEnergy();
				}
			}

			//MMFFMolecule m_optimized = ff.getMMFFMolecule();
			Conformer conf_optimized = new Conformer(conformer);
			copyFFMolCoordsToConformer(conf_optimized,ff,n_atoms);
			return conf_optimized;
		}
		catch (BadAtomTypeException bate) {
			if(energy_out!=null) {
				if (energy_out.length > 1) {
					energy_out[0]= Double.NaN;
					energy_out[1]= Double.NaN;
				}
			}
			return new Conformer(conformer);
		}
		catch(Exception ex) {
			ex.printStackTrace();

			if(energy_out!=null) {
				if (energy_out.length > 1) {
					energy_out[0]= Double.NaN;
					energy_out[1]= Double.NaN;
				}
			}
			return new Conformer(conformer);
		}
	}

	private static void copyFFMolCoordsToConformer(Conformer mol, ForceFieldMMFF94 ff, int n_atoms) {
		//int n_atoms = mol.getMolecule().getAllAtoms();
		for (int atom=0; atom<n_atoms ; atom++) {
			mol.setX(atom, ff.getCurrentPositions()[atom*3]);
			mol.setY(atom, ff.getCurrentPositions()[atom*3+1]);
			mol.setZ(atom, ff.getCurrentPositions()[atom*3+2]);
		}
	}*/
}
