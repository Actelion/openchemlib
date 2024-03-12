package com.actelion.research.chem.chemicalspaces.ptree.synthon;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.chemicalspaces.synthon.SynthonReactor;
import com.actelion.research.chem.descriptor.pharmacophoretree.DescriptorHandlerPTree;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreTree;

public class PharmTreeSynthon {
	
	private String idcode;
	private Map<String,Integer> connectorLabels;
	private String pharmTreeEncoded;
	private DescriptorHandlerPTree dhp;
	private IDCodeParser idcParser;
	private String id;
	
	public PharmTreeSynthon(StereoMolecule mol) {
		connectorLabels = new HashMap<String,Integer>();
		dhp = new DescriptorHandlerPTree();
		idcParser = new IDCodeParser();
		process(mol);
	}
	
	private void process(StereoMolecule mol) {
		mol = new Canonizer(mol).getCanMolecule(true);
		mol.ensureHelperArrays(Molecule.cHelperParities);
		for(int a=0;a<mol.getAtoms();a++) {
			if(mol.getAtomicNo(a)>=SynthonReactor.CONNECTOR_OFFSET) {
				connectorLabels.put(mol.getAtomLabel(a), a);
			}
		}
		this.idcode = mol.getIDCode();

		
	}
	
	public int getLinkerNoFromConnectorLabel(String label) {
		return StereoMolecule.getAtomicNoFromLabel(label)+1-SynthonReactor.CONNECTOR_OFFSET;
		
	}
	
	public Map<String,Integer> getConnectorLabels() {
		return connectorLabels;
	}

	public PharmacophoreTree getPharmTree() {
		return dhp.decode(pharmTreeEncoded);
	}

	public void setPharmacophoreTree(PharmacophoreTree pTree) {
		pharmTreeEncoded = dhp.encode(pTree);
	}


	public StereoMolecule getStructure() {
		return idcParser.getCompactMolecule(idcode);
	}
	
	/**
	 * only take atoms that are part of the shortest paths between connector atoms
	 * @return
	 */
	
	public StereoMolecule createMinimalSynthon() {
		StereoMolecule synthon = getStructure();
		List<Integer> connectorAtoms = new ArrayList<>(connectorLabels.values());
	
		return cutConnectorPaths(synthon, connectorAtoms);
	}
	
	public StereoMolecule cutConnectorPaths(StereoMolecule synthon, List<Integer> connectorAtoms) {

		List<int[]> combis = new ArrayList<>();
		for(int i=0;i<connectorAtoms.size();i++) {
			for(int j=i+1;j<connectorAtoms.size();j++) {
				combis.add(new int[] {connectorAtoms.get(i),connectorAtoms.get(j)});
			}
		}
		Set<Integer> includedAtoms = new HashSet<>();
		int maxPathLength = 50;
		if(combis.size()==0) {
			int connectorAtom = connectorAtoms.get(0);
			includedAtoms.add(connectorAtom);
			includedAtoms.add(synthon.getConnAtom(connectorAtom, 0));
		}
		else {
			for(int[] combi : combis) {
				int[] pathAtoms = new int[50];
				Arrays.fill(pathAtoms, -1);
				synthon.getPath(pathAtoms, combi[0], combi[1], maxPathLength, null);
				includedAtoms.addAll(Arrays.stream(pathAtoms).boxed().filter(e -> e >= 0).collect(Collectors.toList()));
				includedAtoms.add(combi[0]);
				includedAtoms.add(combi[1]);
			}
		}

		boolean[] isAtomIncluded = new boolean[synthon.getAtoms()];
		for(int a=0;a<synthon.getAtoms();a++) {
			if(includedAtoms.contains(a))
				isAtomIncluded[a] = true;
		}
		StereoMolecule minimalSynthon = new StereoMolecule();
		synthon.copyMoleculeByAtoms(minimalSynthon, isAtomIncluded, true,null);
		return minimalSynthon;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	

	
}
