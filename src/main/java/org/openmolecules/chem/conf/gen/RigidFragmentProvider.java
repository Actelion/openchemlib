package org.openmolecules.chem.conf.gen;

import com.actelion.research.chem.*;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.forcefield.mmff.BadAtomTypeException;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import org.openmolecules.chem.conf.so.ConformationSelfOrganizer;
import org.openmolecules.chem.conf.so.SelfOrganizedConformer;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.TreeSet;

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
 * - When locating the same fragment in a different molecule, then atom order may be different. Thus, we need to normalize.
 * For this purpose we use the graphIndex of the Canonizer used to create the fragment's idcode and store coordinates in
 * canonical atom order.<br>
 * - For molecule stereo centers, which are gone in the fragment, copying of molecule parities to the fragment ensures
 * proper 3D-coordinates. For the Canonizer graphIndex to reflect the original configuration, we need to make sure,
 * that up/down-bonds are copied (parities won't do), which are now overspecifying the non-stereo center.
 * And we use the Canonizer mode CONSIDER_STEREOHETEROTOPICITY to distinguish enantio- and diastereo-topic neighbours,
 * which doesn't change the idcode, but is reflected in the graphindex, because stereoheterotopic atoms are ranked
 * now differently before tie-breaking.<br>
 * - Prochiral fragments in symmetrical environment: If a new fragment with a potential stereo center is found
 * first in a symmetrical molecule, such that the potential stereo center is none, then no parity and no up/down-bond
 * are copied and 3D-coordinates randomly reflect one of the two options. If the same fragment is later retrieved
 * from the cache when found in a chiral situation, then cached coordinates have a 50% change to be wrong. Counter
 * measure: For every pro-chiral atom we introduce an arbitrary parity and according up/down bond before
 * generating coordinates and graphIndex.
 */
public class RigidFragmentProvider {
	public static boolean sPrintParityFragment = false;
	private static int MAX_CONFORMERS = 16;

	private static final boolean DEBUG_INFO_MMFF = false;

	// Random seed for initializing the SelfOrganizer.
	private long mRandomSeed;
	private boolean mOptimizeFragments;
	private RigidFragmentCache mCache;

	private static TreeSet<String> sDebugFragmentSet;
	private static BufferedWriter sDebugWriter;

	public RigidFragmentProvider(long randomSeed, RigidFragmentCache cache, boolean optimizeRigidFragments) {
		mRandomSeed = randomSeed;
		mCache = cache;
		mOptimizeFragments = optimizeRigidFragments;
		}

	/**
	 * @param debugFragmentSet if given, then only listed fragment are generated
	 * @param debugWriter if given, then generated fragments are written as SD-entries into file
	 */
	public static void setDebugMode(TreeSet<String> debugFragmentSet, BufferedWriter debugWriter) {
		sDebugFragmentSet = debugFragmentSet;
		sDebugWriter = debugWriter;
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

		fragment.setFragment(true); // if can encode as fragment, because H-atoms are converted deuterium

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

				if (isOuterShellAtom[atom]) // for the ConformationSelfOrganizer to neglect
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

				// convert all plain hydrogen to deuterium that we don't loose them in the idcode
				if (fragment.getAtomicNo(fragmentAtom) == 1)
					fragment.setAtomMass(fragmentAtom, 2);

				fragmentAtom++;
			}
		}

		Conformer[] conformers = null;
		double[] likelihood = null;
		Canonizer canonizer = null;
		String key = null;
		boolean invertedEnantiomer = false;

		boolean putFragmentIntoCache = false;

		// Generate stereo parities for all potential stereo configurations in fragment.
		// If one or more potential stereo configurations are unknown, then the fragment doesn't qualify to be cached.
		if (mCache != null) {
			// we may use original parities for coordinate generation, if fragment doesn't qualify for caching
			int[] originalAtomParity = new int[atomCount];
			for (int atom=0; atom<atomCount; atom++) {
				originalAtomParity[atom] = fragment.getAtomParity(atom);
				if (fragment.isAtomParityPseudo(atom));
					originalAtomParity[atom] = -originalAtomParity[atom];
				}
			int[] originalBondParity = new int[bondCount];
			for (int bond=0; bond<bondCount; bond++) {
				originalBondParity[bond] = fragment.getBondParity(bond);
				if (fragment.isBondParityPseudo(bond));
					originalBondParity[bond] = -originalBondParity[bond];
				}

			// By distinguishing equal ranking atoms, if they have free valencens, we detect all possible stereo features
			canonizer = new Canonizer(fragment, Canonizer.TIE_BREAK_FREE_VALENCE_ATOMS);
			canonizer.setParities();

			// we don't cache fragments with unspecified stereo configurations
			putFragmentIntoCache = true;

			fragment.ensureHelperArrays(Molecule.cHelperNeighbours);
			for (int atom=0; atom<fragment.getAtoms(); atom++) {
				if (fragment.getAtomParity(atom) == Molecule.cAtomParityUnknown) {
					putFragmentIntoCache = false;
					break;
					}
				}
			for (int bond=0; bond<fragment.getBonds(); bond++) {
				if (fragment.getBondParity(bond) == Molecule.cBondParityUnknown) {
					putFragmentIntoCache = false;
					break;
					}
				}

			if (!putFragmentIntoCache) {    // restore orignal fragment parities
				for (int atom=0; atom<atomCount; atom++)
					fragment.setAtomParity(atom, Math.abs(originalAtomParity[atom]), originalAtomParity[atom] < 0);
				for (int bond=0; bond<bondCount; bond++)
					fragment.setBondParity(bond, Math.abs(originalBondParity[bond]), originalBondParity[bond] < 0);
				}
			}

		// no matter, whether we use the original parities or freshly calculated parities,
		// we need to define them to be valid for the coordinate self-organization
		fragment.setParitiesValid(0);

		// Check, whether we have the fragment already in the cache.
		// If yes, then map coordinates from canonical order and mirror coordinates, if needed.
		// Coordinates are store normalized to one enantiomer
		if (mCache != null) {
			invertedEnantiomer = canonizer.normalizeEnantiomer();
			key = canonizer.getIDCode();

			if (sDebugFragmentSet != null && (!sDebugFragmentSet.contains(key) || mCache.containsKey(key)))
				return null;

			RigidFragmentCache.CacheEntry cacheEntry = mCache.get(key);

			if (cacheEntry != null) {
				// convert from canonical coordinates back to fragment
				int[] graphIndex = canonizer.getGraphIndexes();
				conformers = new Conformer[cacheEntry.coordinates.length];
				for (int i=0; i<conformers.length; i++) {
					for (int j = 0; j<fragment.getAllAtoms(); j++) {
						Coordinates coords = cacheEntry.coordinates[i][graphIndex[j]];
						fragment.setAtomX(j, coords.x);
						fragment.setAtomY(j, coords.y);
						fragment.setAtomZ(j, invertedEnantiomer ? -coords.z : coords.z);
					}
					conformers[i] = new Conformer(fragment);
				}
				likelihood = cacheEntry.likelihood;
			}
		}
		else if (sDebugFragmentSet != null)
			return null;

		if (conformers == null) {
			ConformationSelfOrganizer selfOrganizer = new ConformationSelfOrganizer(fragment, true);
			selfOrganizer.initializeConformers(mRandomSeed, MAX_CONFORMERS);

			// Generate multiple low constrain conformers
			ArrayList<SelfOrganizedConformer> conformerList = new ArrayList<>();
			SelfOrganizedConformer bestConformer = selfOrganizer.getNextConformer();
			conformerList.add(bestConformer);
			SelfOrganizedConformer conformer = selfOrganizer.getNextConformer();
			while (conformer != null) {
				conformerList.add(conformer);
				conformer = selfOrganizer.getNextConformer();
			}

			conformers = conformerList.toArray(new Conformer[0]);
			likelihood = new double[conformers.length];
			double likelyhoodSum = 0.0;
			for (int i = 0; i < conformers.length; i++) {
				likelihood[i] = ((SelfOrganizedConformer) conformers[i]).getLikelyhood();
				likelyhoodSum += likelihood[i];
			}
			if (likelyhoodSum != 0.0) {
				for (int i = 0; i < conformers.length; i++)
					likelihood[i] /= likelyhoodSum;
			}

			if(mOptimizeFragments) {
				// @TODO: update the likelihoods according to the resulting energies..)
				for(int zi=0;zi<conformers.length;zi++) {
					if(DEBUG_INFO_MMFF){
						System.out.println("FFMIN: Minimize Conformer "+zi);
					}
					double[] ff_energies = new double[2];
					conformers[zi] = minimizeConformer(conformers[zi], ff_energies);
					if(DEBUG_INFO_MMFF){
						System.out.println("FFMIN: Result: E_start= "+ff_energies[1]+" , E_end= "+ff_energies[0]);
					}
				}
			}

			if (putFragmentIntoCache) {
				int[] graphIndex = canonizer.getGraphIndexes();
				Coordinates[][] coords = new Coordinates[conformers.length][fragment.getAllAtoms()];
				for (int i=0; i<coords.length; i++) {
					for (int j = 0; j<coords[i].length; j++) {
						Coordinates xyz = conformers[i].getCoordinates(j);
						coords[i][graphIndex[j]] = new Coordinates(xyz.x, xyz.y, invertedEnantiomer ? -xyz.z : xyz.z);
						}
					}

				mCache.put(key, new RigidFragmentCache.CacheEntry(coords, likelihood));

				if (sDebugWriter != null) {
					conformers[0].toMolecule(fragment);
					if (invertedEnantiomer)
						for (int atom=0; atom<fragment.getAllAtoms(); atom++)
							fragment.setAtomZ(atom, -fragment.getAtomZ(atom));
					try {
//						sDebugWriter.write("parent"+new MolfileCreator(mol).getMolfile());
//						sDebugWriter.write("$$$$\n");

						sDebugWriter.write(key+new MolfileCreator(fragment).getMolfile());
						sDebugWriter.write(">  <graph-index>\n");
						for (int gi:graphIndex)
							sDebugWriter.write(" "+gi);
						sDebugWriter.newLine();
						sDebugWriter.newLine();
						sDebugWriter.write("$$$$\n");
					} catch (Exception e) {}
					return null;
				}
			}
		}

		return new RigidFragment(coreAtomCount, coreToFragmentAtom, fragmentToOriginalAtom,
				extendedToFragmentAtom, originalToExtendedAtom,
				conformers, likelihood);
	}

	private boolean isPotentialStereoCenter(StereoMolecule mol, int atom) {
		if (mol.getAtomicNo(atom) != 5
		 && mol.getAtomicNo(atom) != 6
		 && mol.getAtomicNo(atom) != 7
		 && mol.getAtomicNo(atom) != 14
		 && mol.getAtomicNo(atom) != 15
		 && mol.getAtomicNo(atom) != 16)
			return false;

		if (mol.getAtomPi(atom) != 0) {
			if (mol.isCentralAlleneAtom(atom))
				return true;

			if (mol.getAtomicNo(atom) != 15
			 && mol.getAtomicNo(atom) != 16)
				return false;
			}

		if (mol.getConnAtoms(atom) < 3 || mol.getAllConnAtoms(atom) > 4)
			return false;

		// no carbenium
		if (mol.getAtomCharge(atom) > 0 && mol.getAtomicNo(atom) == 6)
			return false;

		// no trivalent boron
		if (mol.getAtomicNo(atom) == 5 && mol.getAllConnAtoms(atom) != 4)
			return false;

		// don't consider tetrahedral nitrogen, unless found to qualify for parity calculation
		if (mol.getAtomicNo(atom) == 7 && mol.getConnAtoms(atom) < 4)
			return false;

		int[] canRank = new int[mol.getConnAtoms(atom)];
		Canonizer canonizer = new Canonizer(mol, Canonizer.CREATE_SYMMETRY_RANK);
		for (int i=0; i<canRank.length; i++)
			canRank[i] = canonizer.getSymmetryRank(mol.getConnAtom(atom, i));

		// create array to remap connAtoms according to canRank order
		boolean[] hasOpenValence = new boolean[mol.getConnAtoms(atom)];
		int[] sortedRank = new int[mol.getConnAtoms(atom)];
		boolean[] neighbourUsed = new boolean[4];
		for (int i=0; i<mol.getConnAtoms(atom); i++) {
			int highestRank = -1;
			int highestConn = 0;
			for (int j=0; j<mol.getConnAtoms(atom); j++) {
				if (!neighbourUsed[j]) {
					if (highestRank < canRank[j]) {
						highestRank = canRank[j];
						highestConn = j;
						}
					}
				}
			neighbourUsed[highestConn] = true;
			sortedRank[i] = highestRank;
			hasOpenValence[i] = isOpenValenceSubstituent(mol, atom, mol.getConnAtom(atom, highestConn));
			}

		for (int i=1; i<canRank.length; i++)
			if (sortedRank[i-1] == sortedRank[i]
			 && !hasOpenValence[i-1]
			 && !hasOpenValence[i])
				return false;

		return true;
		}

	private boolean isOpenValenceSubstituent(StereoMolecule mol, int root, int first) {
		if (mol.getConnAtoms(first) == 1)
			return false;

		boolean[] isUsed = new boolean[mol.getAtoms()];
		isUsed[root] = true;
		isUsed[first] = true;

		int[] graphAtom = new int[mol.getAtoms()];
		graphAtom[0] = first;
		int current = 0;
		int highest = 0;
		while (current <= highest) {
			for (int i=0; i<mol.getConnAtoms(graphAtom[current]); i++) {
				int candidate = mol.getConnAtom(graphAtom[current], i);
				if (!isUsed[candidate]) {
					if (mol.getFreeValence(candidate) != 0)
						return true;
					isUsed[candidate] = true;
					graphAtom[++highest] = candidate;
					}
				}
			current++;
			}
		return false;
		}

	/**
	 * @param conformer
	 * @param energy_out [0] is the final energy, [1] is the starting energy
	 * @return
	 */
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
	}
}
