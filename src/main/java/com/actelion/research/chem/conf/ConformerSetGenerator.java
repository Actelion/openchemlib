package com.actelion.research.chem.conf;

import java.util.HashMap;
import java.util.Map;

import org.openmolecules.chem.conf.gen.ConformerGenerator;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.forcefield.mmff.ForceFieldMMFF94;

public class ConformerSetGenerator {
	
	private int mMaxNrConfs;
	private int mStrategy;
	private boolean mUseFF;
	private ConformerGenerator mConfGen;
	
	public ConformerSetGenerator(int maxNrConfs, int strategy, boolean useFF, long seed) {
		mMaxNrConfs = maxNrConfs;
		mStrategy = strategy;
		mUseFF = useFF;
		mConfGen = new ConformerGenerator(seed);
	}
	
	/**
	 * STRATEGY_LIKELY_RANDOM was evaluated to be the best strategy for reproducing bioactive
	 * conformers (J.W. 05/19)
	 */
	
	public ConformerSetGenerator() {
		this(200,ConformerGenerator.STRATEGY_LIKELY_RANDOM,false,0L);
		
	}
	
	public ConformerSetGenerator(boolean useFF) {
		this(200,ConformerGenerator.STRATEGY_LIKELY_RANDOM,useFF,0L);
		
	}
	
	public ConformerSetGenerator(boolean useFF, long seed) {
		this(200,ConformerGenerator.STRATEGY_LIKELY_RANDOM,useFF,seed);
		
	}
	
	public ConformerSet generateConformerSet(StereoMolecule mol) {   
		StereoMolecule m = new StereoMolecule(mol);
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

