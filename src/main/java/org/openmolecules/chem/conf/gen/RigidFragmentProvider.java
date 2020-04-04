package org.openmolecules.chem.conf.gen;

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

public class RigidFragmentProvider {
	private static int MAX_CONFORMERS = 16;

	private static final boolean DEBUG_INFO_MMFF = false;

	// Random seed for initializing the SelfOrganizer.
	private long mRandomSeed;
	private boolean mOptimizeFragments;
	private RigidFragmentCache mCache;

	public RigidFragmentProvider(long randomSeed, RigidFragmentCache cache, boolean optimizeRigidFragments) {
		mRandomSeed = randomSeed;
		mCache = cache;
		mOptimizeFragments = optimizeRigidFragments;
	}

	public RigidFragment createFragment(StereoMolecule mol, int[] fragmentNo, int fragmentIndex) {
		int coreAtomCount = 0;
		int atomCount = 0;
		int extendedAtomCount = 0;

		// mark all atoms with specified fragmentNo and two layers around it
		boolean[] includeAtom = new boolean[mol.getAllAtoms()];
		boolean[] isOuterShellAtom = new boolean[mol.getAllAtoms()];
		boolean[] isCoreFragment = new boolean[mol.getAllAtoms()];

		for (int atom = 0; atom < mol.getAllAtoms(); atom++) {
			if (fragmentNo[atom] == fragmentIndex) {
				isCoreFragment[atom] = true;
				includeAtom[atom] = true;
				atomCount++;
				coreAtomCount++;
				for (int i = 0; i < mol.getConnAtoms(atom); i++) {
					int connAtom = mol.getConnAtom(atom, i);
					if (fragmentNo[connAtom] != fragmentIndex) {
						includeAtom[connAtom] = true;
						atomCount++;
						extendedAtomCount++;
						for (int j = 0; j < mol.getAllConnAtoms(connAtom); j++) {
							int connconn = mol.getConnAtom(connAtom, j);
							if (fragmentNo[connconn] != fragmentIndex) {
								isOuterShellAtom[connconn] = true;
								includeAtom[connconn] = true;
								atomCount++;
							}
						}
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
		fragment.setParitiesValid(0);

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

				if (isOuterShellAtom[atom])
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

		if (mCache != null) {
			canonizer = new Canonizer(fragment);
			key = canonizer.getIDCode();

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
						fragment.setAtomZ(j, coords.z);
					}
					conformers[i] = new Conformer(fragment);
				}
				likelihood = cacheEntry.likelihood;
			}
		}

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

			if (mCache != null) {
				int[] graphIndex = canonizer.getGraphIndexes();
				Coordinates[][] coords = new Coordinates[conformers.length][fragment.getAllAtoms()];
				for (int i=0; i<coords.length; i++)
					for (int j=0; j<coords[i].length; j++)
						coords[i][graphIndex[j]] = conformers[i].getCoordinates(j);

				mCache.put(key, new RigidFragmentCache.CacheEntry(coords, likelihood));
			}
		}

		return new RigidFragment(coreAtomCount, coreToFragmentAtom, fragmentToOriginalAtom,
				extendedToFragmentAtom, originalToExtendedAtom,
				conformers, likelihood);
	}

	/**
	 * @param conformer
	 * @param energy_out [0] is the final energy,, [1] is the starting energy
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
