package com.actelion.research.chem.conf;

import com.actelion.research.calc.ThreadMaster;
import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import org.openmolecules.chem.conf.gen.ConformerGenerator;
import org.openmolecules.chem.conf.gen.RigidFragmentCache;
import org.openmolecules.chem.conf.so.ConformationSelfOrganizer;

import java.util.HashMap;
import java.util.Map;

public class ConformerSetGenerator {

	public static final int CONFORMERS=200;

	private int mMaxNrConfs;
	private int mStrategy;
	private boolean mUseFF;
	private long mSeed;
	private static long DEFAULT_SEED = 12345L;
	private ThreadMaster threadMaster;
	
	public ConformerSetGenerator(int maxNrConfs, int strategy, boolean useFF, long seed) {
		mMaxNrConfs = maxNrConfs;
		mStrategy = strategy;
		mUseFF = useFF;
		mSeed = seed;
		RigidFragmentCache fragCache = RigidFragmentCache.getDefaultInstance();
		fragCache.loadDefaultCache();
	}
	
	
	/**
	 * STRATEGY_LIKELY_RANDOM was evaluated to be the best strategy for reproducing bioactive
	 * conformers (J.W. 05/19)
	 * 
	 */
	
	public ConformerSetGenerator(int maxConfs) {
		this(maxConfs,ConformerGenerator.STRATEGY_LIKELY_RANDOM,false,DEFAULT_SEED);
	}
	
	public ConformerSetGenerator() {
		this(CONFORMERS,ConformerGenerator.STRATEGY_LIKELY_RANDOM,false,DEFAULT_SEED);
	}
	
	public ConformerSetGenerator(boolean useFF) {
		this(CONFORMERS,ConformerGenerator.STRATEGY_LIKELY_RANDOM,useFF,DEFAULT_SEED);
	}
	
	public ConformerSetGenerator(boolean useFF, long seed) {
		this(CONFORMERS,ConformerGenerator.STRATEGY_LIKELY_RANDOM,useFF,seed);
	}
	
	public void setThreadMaster(ThreadMaster tm) {
		this.threadMaster = tm;
	}

	/**
	 * Generates a set of distinct conformers of the canonical largest fragment of the passed molecule.
	 * @param mol
	 * @return
	 */
	public ConformerSet generateConformerSet(StereoMolecule mol) {   
		ConformerSet confSet = new ConformerSet();

/* Simplified; TLS 31Jan2023
		StereoMolecule canMol = new StereoMolecule(mol);
		canMol.stripSmallFragments();
		ConformerGenerator.addHydrogenAtoms(canMol);        not necessary, because the ConformerGenerator does that anyway
		Canonizer can = new Canonizer(canMol);
		canMol = can.getCanMolecule(true);
		canMol.ensureHelperArrays(Molecule.cHelperCIP);     not allowed, because the molecule may not have coordinates!
*/
		StereoMolecule largestFragment = mol.getCompactCopy();
		largestFragment.stripSmallFragments();
		StereoMolecule canonicalFragment = new Canonizer(largestFragment).getCanMolecule(true);

		int maxTorsionSets = (int) Math.max(2 * mMaxNrConfs, (1000 * Math.sqrt(mMaxNrConfs)));
		ConformerGenerator cg = new ConformerGenerator(mSeed,false);
		cg.setThreadMaster(this.threadMaster);
		Map<String, Object> ffOptions = null;
		if(mUseFF) {
			ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
			ffOptions = new HashMap<String, Object>();
			ffOptions.put("dielectric constant", 4.0);
		}

		if (cg.initializeConformers(canonicalFragment, mStrategy, maxTorsionSets, false)) {
			for (int i = 0; i < mMaxNrConfs; i++) {
				Conformer conformer = cg.getNextConformer();
				if (conformer == null && i==0) {
					ConformationSelfOrganizer sampler = new ConformationSelfOrganizer(canonicalFragment, true);
					conformer = sampler.generateOneConformer(mSeed);
				}

				if (conformer == null)
					break;

				if(mUseFF) {
					conformer.copyTo(canonicalFragment);
					ForceFieldMMFF94 mmff = new ForceFieldMMFF94(canonicalFragment, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
					mmff.minimise();
					conformer.copyFrom(canonicalFragment);
				}
				confSet.add(conformer);
				if(threadMaster!=null && threadMaster.threadMustDie())
					break;
			}
		}

		return confSet;
	}
}
