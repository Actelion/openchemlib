package com.actelion.research.chem.combinatorialspace;

import java.util.Map;

import com.actelion.research.chem.StereoMolecule;

/**
 * a Synthon is a reaction-specific representation of a building block, already encoding structural changes
 * that occur upon the reaction (e.g. ring formations). A chemical reaction defined by synthons can easily be conducted
 * by simply forming bonds between the connector atoms
 * reactions: lists the reactions in which this synthon can participate, the 
 * value is the reactant ID the synthon has for the given reaction
 * @author wahljo1
 *
 */

public class Synthon {

	private StereoMolecule origBB;
	private StereoMolecule synthon;
	private Map<String,Integer> reactions; 
	
	public Synthon(StereoMolecule origBB, StereoMolecule synthon) {
		this.origBB = synthon;
		this.synthon = synthon;
	}
	
	public void addReaction(String rxn, int rxnID) {
		reactions.put(rxn, rxnID);
	}

	public StereoMolecule getOrigBB() {
		return origBB;
	}

	public StereoMolecule getSynthon() {
		return synthon;
	}

	public Map<String, Integer> getReactions() {
		return reactions;
	}
	
	public boolean areSynthonsCompatible(Synthon synthon1, Synthon synthon2 ) {
		for(String rxn1 : synthon1.reactions.keySet()) {
			int id1 = synthon1.reactions.get(rxn1);
			if(synthon2.reactions.containsKey(rxn1)) {
				int id2 = synthon2.reactions.get(rxn1);
				if(id1!=id2)
					return true;
			}
		}
		return false;
			
		}
	
}
