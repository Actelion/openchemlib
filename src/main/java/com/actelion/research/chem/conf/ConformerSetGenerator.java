package com.actelion.research.chem.conf;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;
import org.openmolecules.chem.conf.gen.ConformerGenerator;

import java.util.HashMap;
import java.util.Map;

public class ConformerSetGenerator {
	
	private int mMaxNrConfs;
	private int mStrategy;
	private boolean mUseFF;
	private ConformerGenerator mConfGen;
	private static  long DEFAULT_SEED = 12345L;
	
	public ConformerSetGenerator(int maxNrConfs, int strategy, boolean useFF, long seed) {
		mMaxNrConfs = maxNrConfs;
		mStrategy = strategy;
		mUseFF = useFF;
		mConfGen = new ConformerGenerator(seed, useFF);
	}
	
	/**
	 * STRATEGY_LIKELY_RANDOM was evaluated to be the best strategy for reproducing bioactive
	 * conformers (J.W. 05/19)
	 * 
	 */
	
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
		StereoMolecule m = new StereoMolecule(mol);
		m.ensureHelperArrays(Molecule.cHelperCIP);
		m.stripSmallFragments();
		ConformerGenerator.addHydrogenAtoms(m);
		Canonizer can = new Canonizer(m);
		m = can.getCanMolecule(true);
		m.ensureHelperArrays(Molecule.cHelperCIP);
		if(mUseFF) {
			ForceFieldMMFF94.initialize(ForceFieldMMFF94.MMFF94SPLUS);
		}
		Map<String, Object> ffOptions = new HashMap<String, Object>();
		ConformerSet confSet = new ConformerSet();
		mConfGen.initializeConformers(m,mStrategy,100000,false);
		boolean nextConf = true;
		int i=0;
		while(i<mMaxNrConfs && nextConf) {
			m = mConfGen.getNextConformerAsMolecule(m);
			if(m==null) {
				nextConf=false;
			}
			else {
				m.ensureHelperArrays(Molecule.cHelperCIP);
				if(mUseFF) {
					ffOptions.put("dielectric constant", 4.0);
					ForceFieldMMFF94 mmff = new ForceFieldMMFF94(m, ForceFieldMMFF94.MMFF94SPLUS, ffOptions);
					mmff.minimise();
					confSet.add(new Conformer(m));
				}
				else {
					confSet.add(new Conformer(m));
				}
				
			
		}
			i++;
		}
		return confSet;
	}
}

