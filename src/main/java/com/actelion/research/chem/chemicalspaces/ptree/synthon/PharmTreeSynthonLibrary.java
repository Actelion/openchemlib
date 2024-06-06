package com.actelion.research.chem.chemicalspaces.ptree.synthon;

import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.chemicalspaces.ptree.PharmTreeSynthonReactionHelper;


public class PharmTreeSynthonLibrary {
	
	private List<List<PharmTreeSynthon>> synthons;
	private PharmTreeSynthonReactionHelper rxnHelper;
	private String reactionID;
	
	public PharmTreeSynthonLibrary(List<List<PharmTreeSynthon>> synthons) {
		this.synthons = synthons;
		rxnHelper = new PharmTreeSynthonReactionHelper(getGenericReactants());

	}

	public List<List<PharmTreeSynthon>> getSynthons() {
		return synthons;
	}

	public PharmTreeSynthonReactionHelper getReactionHelper() {
		return rxnHelper;
	}

	public void setReactionHelper(PharmTreeSynthonReactionHelper reactor) {
		this.rxnHelper = reactor;
	}
	
	public List<StereoMolecule> getGenericReactants() {
		List<StereoMolecule> genericReactants = new ArrayList<>();
		for(List<PharmTreeSynthon> synthonList : synthons) {


			StereoMolecule minimalSynthon = synthonList.get(0).createMinimalSynthon();

			try {
				genericReactants.add(minimalSynthon);
			} catch (Exception e) {
				System.err.println(minimalSynthon.getIDCode());
				e.printStackTrace();

			}
		}
		return genericReactants;
	}

	public String getReactionID() {
		return reactionID;
	}

	public void setReactionID(String reactionID) {
		this.reactionID = reactionID;
	}
	
	
	
	

}
