package com.actelion.research.chem.conf;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import org.openmolecules.chem.conf.gen.ConformerGenerator;
import org.openmolecules.chem.conf.gen.RigidFragmentCache;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ConformerSetGenerator {
	
	private int mMaxNrConfs;
	private int mStrategy;
	private boolean mUseFF;
	private static  long DEFAULT_SEED = 12345L;
	
	public ConformerSetGenerator(int maxNrConfs, int strategy, boolean useFF, long seed) {
		mMaxNrConfs = maxNrConfs;
		mStrategy = strategy;
		mUseFF = useFF;
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
		this(200,ConformerGenerator.STRATEGY_LIKELY_RANDOM,false,DEFAULT_SEED);
		
	}
	
	public ConformerSetGenerator(boolean useFF) {
		this(200,ConformerGenerator.STRATEGY_LIKELY_RANDOM,useFF,DEFAULT_SEED);
		
	}
	
	public ConformerSetGenerator(boolean useFF, long seed) {
		this(200,ConformerGenerator.STRATEGY_LIKELY_RANDOM,useFF,seed);
		
	}
	
	public ConformerSet generateConformerSet(StereoMolecule mol) {   
		ConformerSet confSet = new ConformerSet();
		StereoMolecule m = new StereoMolecule(mol);
		m.ensureHelperArrays(Molecule.cHelperCIP);
		m.stripSmallFragments();
		ConformerGenerator.addHydrogenAtoms(m);
		Canonizer can = new Canonizer(m);
		m = can.getCanMolecule(true);
		m.ensureHelperArrays(Molecule.cHelperCIP);
		int maxTorsionSets = (int) Math.max(2 * mMaxNrConfs, (1000 * Math.sqrt(mMaxNrConfs)));
		ConformerGenerator cg = new ConformerGenerator(12345L,false);
		cg.initializeConformers(m, mStrategy, maxTorsionSets, false);
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		if(mUseFF)
				ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		for (int i = 0; i < mMaxNrConfs; i++) {
			Conformer conformer = null;
			conformer = cg.getNextConformer();
			if (conformer == null)
				break;
			if(mUseFF) {
				StereoMolecule molecule = conformer.toMolecule(mol);
				ffOptions.put("dielectric constant", 4.0);
				ForceFieldMMFF94 mmff = new ForceFieldMMFF94(molecule, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
				mmff.minimise();
				for(int a=0;a<molecule.getAllAtoms();a++) {
					conformer.setX(a, molecule.getAtomX(a));
					conformer.setY(a, molecule.getAtomY(a));
					conformer.setZ(a, molecule.getAtomZ(a));
				}
				
			}
			confSet.add(conformer);
		}

		return confSet;
	}
}

