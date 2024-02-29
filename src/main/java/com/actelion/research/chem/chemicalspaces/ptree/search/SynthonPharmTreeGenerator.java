package com.actelion.research.chem.chemicalspaces.ptree.search;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.actelion.research.chem.Canonizer;
import com.actelion.research.chem.IDCodeParser;
import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.SSSearcher;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.chemicalspaces.ptree.PharmTreeSynthonReactionHelper;
import com.actelion.research.chem.chemicalspaces.ptree.synthon.PharmTreeSynthon;
import com.actelion.research.chem.chemicalspaces.ptree.synthon.PharmTreeSynthonLibrary;
import com.actelion.research.chem.chemicalspaces.synthon.SynthonReactor;
import com.actelion.research.chem.descriptor.pharmacophoretree.FeatureCalculator;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreNode;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreTree;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreTreeGenerator;
import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;

/**
 * Generates PharmacophoreTrees (PharmTrees) for building blocks (or fragments) that are part of virtual chemical space.
 * Building blocks are represented as synthon structures, containing reactive linkers (R1,R2,R3,R4). 
 * In order to make the searches in the virtual chemical space using PharmTrees more efficient, atoms that will be part 
 * of the same ring (formed by the reaction), are merged into one node. Thereby, also the number of linkers required
 * to make the ring formation is reduced
 * @author Joel Wahl
 *
 */


public class SynthonPharmTreeGenerator {
	private PharmTreeSynthonLibrary synthonLib;
	private List<StereoMolecule> genericSynthons;
	private PharmTreeSynthonReactionHelper rxnHelper;

	public SynthonPharmTreeGenerator(PharmTreeSynthonLibrary synthonLib) {
		this.synthonLib = synthonLib;
		this.genericSynthons = synthonLib.getGenericReactants();
		this.rxnHelper = synthonLib.getReactionHelper();
		furtherProcessGenericSynthons();
	}
	
	public void processBuildingBlocks() {
		Map<Integer,Set<Integer>> ringReactantIndeces = rxnHelper.getRingReactantIndeces();
		Map<Integer,Set<Integer>> ringWithLinks = rxnHelper.getRingWithLinks();
		StereoMolecule genericProduct = rxnHelper.getGenericProduct();
		Map<Integer,List<Integer>> reactantsWithLinkers = rxnHelper.getReactantsWithLinkers().entrySet().stream()
			    .collect(Collectors.toMap(e -> e.getKey(), e -> new ArrayList<Integer>(e.getValue()))); //deep copy of map, since it will be modified later
		
		Set<Integer> linkersToDelete = new HashSet<Integer>();
		for(int ring : ringReactantIndeces.keySet()) {
			int nLinkers = ringReactantIndeces.get(ring).size()-1; 
			//number of linkers actually needed to describe ring closure reaction is equal to number of reactants that are involved in forming the ring minus 1
			int toDelete = ringWithLinks.get(ring).size()-nLinkers; // number of links that can be removed 
			Iterator<Integer> iterator = ringWithLinks.get(ring).iterator();
			for(int i=0;i<toDelete;i++) {
				int del = iterator.next();
				linkersToDelete.add(del);
				reactantsWithLinkers.forEach((k,v) -> v.remove(Integer.valueOf(del)));		
			}
		}
		for(int i=0;i<synthonLib.getSynthons().size();i++) {
			List<PharmTreeSynthon> synthonList = synthonLib.getSynthons().get(i);
			for(int bb=0;bb<synthonList.size();bb++) {
				StereoMolecule buildingBlock = synthonList.get(bb).getStructure();
				for(int a=0;a<buildingBlock.getAtoms();a++) {
					if(buildingBlock.getAtomLabel(a).equals("U")) {
						int bond = buildingBlock.getConnBond(a, 0);
						buildingBlock.setBondOrder(bond, rxnHelper.getBondOrders()[0]);
					}
					else if(buildingBlock.getAtomLabel(a).equals("Np")) {
						int bond = buildingBlock.getConnBond(a, 0);
						buildingBlock.setBondOrder(bond, rxnHelper.getBondOrders()[1]);
					}
					
					else if(buildingBlock.getAtomLabel(a).equals("Pu")) {
						int bond = buildingBlock.getConnBond(a, 0);
						buildingBlock.setBondOrder(bond, rxnHelper.getBondOrders()[2]);
					}
					
					else if(buildingBlock.getAtomLabel(a).equals("Am")) {
						int bond = buildingBlock.getConnBond(a, 0);
						buildingBlock.setBondOrder(bond, rxnHelper.getBondOrders()[3]);
					}
				}
				PharmacophoreTree bbTree = calculatePTree(buildingBlock,genericProduct,i, ringReactantIndeces, reactantsWithLinkers);
				if(bbTree==null)
					continue;
				List<PharmacophoreNode> nodesToBeDeleted = new ArrayList<PharmacophoreNode>();
				for(PharmacophoreNode node : bbTree.getNodes()) {
					if(node.isLinkNode()) {
						if(linkersToDelete.contains(node.getFunctionalities()[0]))
							nodesToBeDeleted.add(node);
					}
				}
				nodesToBeDeleted.forEach(e -> {
					bbTree.removeNode(e);
				});
				
				synthonList.get(bb).setPharmacophoreTree(bbTree);

			}	
		}
		
		
	}
	
	
	
	public PharmacophoreTree calculatePTree(StereoMolecule bb, StereoMolecule genericProduct, int synthonID,
			Map<Integer,Set<Integer>> ringReactantIndeces, Map<Integer,List<Integer>> reactantsWithLinkers) {
		

		FeatureCalculator calculator = new FeatureCalculator(genericProduct);
		calculator.calculate();
		int[][] productAtomFunctionalities = calculator.getAtomFunctionalities();
		double[] productAtomVolumes = PharmacophoreTreeGenerator.getAtomVolumes(genericProduct);
		
		List<int[]> ringFunctionalities = new ArrayList<int[]>();
		List<List<Double>> ringVolumes = new ArrayList<List<Double>>();
		
		for(Set<Integer> ring : rxnHelper.getFormedRings()) {
			int[] func = new int[PharmacophoreCalculator.LIPO_ID+1];
			List<Double> volumes = new ArrayList<Double>();
			for(int ringAtom : ring) {
				int[] atomFunc = productAtomFunctionalities[ringAtom];
				IntStream.range(0,func.length).forEach(e -> func[e]+=atomFunc[e]);
				volumes.add(productAtomVolumes[ringAtom]);
			}
			ringFunctionalities.add(func);
			ringVolumes.add(volumes);
		}
		//from this info, the node properties of the ring nodes representing rings formed
		//in the reaction are fully defined. Now we need to map them to the correct ring indeces
		// and process the building blocks accordingly
		//one building block involved in the ring formation should bear the full ring node, 
		//the others are dummy nodes

		int[] reactantToBBMap = new int[0];
		Map<Integer,Set<Integer>> rings = new HashMap<Integer,Set<Integer>>(); //formed rings are indiced
		if(rxnHelper.getFormedRings().size()>0) {
			SSSearcher searcher = new SSSearcher();
			searcher.setFragment(genericSynthons.get(synthonID));
			searcher.setMolecule(bb);
			searcher.findFragmentInMolecule();
			ArrayList<int[]> matches = searcher.getMatchList();
			if(matches.size()!=1) {
				System.err.println("could not match generic reactant to building block");
				System.err.println(genericSynthons.get(synthonID).getIDCode());
				System.err.println(bb.getIDCode());
				return null;
			}
			else {
				reactantToBBMap = matches.get(0);

			}
					
			for(int a=0;a<bb.getAtoms();a++) {
				if(bb.getAtomicNo(a)>=SynthonReactor.CONNECTOR_OFFSET)
					continue;
				for(int i=0;i<reactantToBBMap.length;i++) {
					if(reactantToBBMap[i]==a) {
						int reactantIndex = i;
						
						int productIndex = rxnHelper.mapSynthonToGenericProductIndex(synthonID,reactantIndex);
						for(int l=0;l<rxnHelper.getFormedRings().size();l++) {
							if(rxnHelper.getFormedRings().get(l).contains(productIndex)) {
								int ringIndex = l;
								rings.putIfAbsent(ringIndex, new HashSet<Integer>());
								rings.get(ringIndex).add(a);
							}
								
							
						}
					}
				}
			}
							
		}
		
		Map<Integer,List<Integer>> atomToNodes = new HashMap<Integer,List<Integer>>();
		PharmacophoreTree bbTree = PharmacophoreTreeGenerator.generate(bb, atomToNodes,rings.values().stream().collect(Collectors.toList()));

		Map<Integer,Integer> ringMainReactants = new HashMap<Integer,Integer>();
		//for every formed ring in the reaction, one reactant has all the features assigned to it, whereas the remaining ones have zero nodes
		// the main reactant should be the one with the highest number of linkers
		for(int ring : ringReactantIndeces.keySet()) {
			Set<Integer> reactants = ringReactantIndeces.get(ring);
			Iterator<Integer> iterator = reactants.iterator();
			int maxLinkers = 0;
			int mainReactant = -1;
			while(iterator.hasNext()) {
				int reactant = iterator.next();
				if(reactantsWithLinkers.get(reactant).size()>maxLinkers) {
					maxLinkers = reactantsWithLinkers.get(reactant).size();
					mainReactant = reactant;
				}
			}
			ringMainReactants.put(ring,mainReactant);
		}
		
		// for the nodes that represent groups of atoms that will be part of a newly formed ring, do the following:
		// if the building block is considered the main reactant for this ring formation, all ring functionalities/volumes
		// are assigned to it
		for(PharmacophoreNode node : bbTree.getNodes()) {
			if(node.isLinkNode())
				continue;
			else {
				
			boolean nodeProcessed = false;
			for(int index=0;index<node.getAtoms().size() && !nodeProcessed; index++) {
				int a = node.getAtoms().get(index);
				for(int i=0;i<reactantToBBMap.length;i++) {
					if(reactantToBBMap[i]==a) {
						int reactantIndex = i;
						int productIndex = rxnHelper.mapSynthonToGenericProductIndex(synthonID,reactantIndex);
						for(int l=0;l<rxnHelper.getFormedRings().size();l++) {
							Set<Integer> formedRing = rxnHelper.getFormedRings().get(l);
							if(formedRing.contains(productIndex)) {
								nodeProcessed = true;
								int ringIndex = l;
								int mainReactantSynthonID = ringMainReactants.get(ringIndex);
								
								if(synthonID==mainReactantSynthonID) {
									int dummyAtomsToAdd = formedRing.size()-node.getAtoms().size();
									for(int da=0;da<dummyAtomsToAdd;da++) {
										node.getAtoms().add(-1);
										node.getWeights().add(1.0);
									}
									List<Double> weights = new ArrayList<Double>();
									IntStream.range(0, node.getAtoms().size()).forEach(e -> weights.add(1.0));
									node.setFunctionalities(ringFunctionalities.get(l));
									node.setVolumes(ringVolumes.get(l));
								}
								else {
									node.setFunctionalities(new int[PharmacophoreCalculator.LIPO_ID+1]);
									node.setVolumes(new ArrayList<Double>());
									node.setWeights(new ArrayList<Double>());
									node.setAtoms(new ArrayList<Integer>());
									node.setRole(PharmacophoreNode.ZERO_NODE);
								}
								node.calculate();
								
						}
						}
					}
				}
			}
			}
		}

		return bbTree;
		
	}
	
	/**
	 * only keep parts of synthons responsible for ring formation
	 * if one synthon is involved in the formation of more than one ring, it is split into different, disconnected
	 * fragments so that atoms that are not involved in ring formation are omitted
	 * @param genericSynthons
	 * @param rxnHelper
	 */
	private void furtherProcessGenericSynthons() {
		List<StereoMolecule> newGenericSynthons = new ArrayList<>();
		Map<Integer,Set<Integer>> ringWithLinks = rxnHelper.getRingWithLinks();
		Map<Integer,Set<Integer>> reactantsWithRings = new HashMap<Integer,Set<Integer>>();
		Map<Integer,Set<Integer>> ringsWithReactants = rxnHelper.getRingReactantIndeces();
		for(int s=0;s<genericSynthons.size();s++) {
			reactantsWithRings.put(s, new HashSet<Integer>());
			for(Map.Entry<Integer, Set<Integer>> entry : ringsWithReactants.entrySet()) {
				if (entry.getValue().contains(s)) {
					reactantsWithRings.get(s).add(entry.getKey());
				}
			}
			StereoMolecule ringFormationSynthon = new StereoMolecule();
			StereoMolecule mol = genericSynthons.get(s);
			if(reactantsWithRings.get(s).size()==0)
				ringFormationSynthon.addMolecule(mol);
			else {
				PharmTreeSynthon synthon = new PharmTreeSynthon(mol);
				List<Integer> allConnectorAtoms = new ArrayList<>(synthon.getConnectorLabels().values());
				List<Integer> usedConnectorAtoms = new ArrayList<>();
				for(int ring : reactantsWithRings.get(s)) {
					synthon = new PharmTreeSynthon(mol);
					List<Integer> connectorAtomsToKeep = new ArrayList<>();
					
					
					for(Map.Entry<String, Integer> entry : synthon.getConnectorLabels().entrySet()) {
						int linkerNo = synthon.getLinkerNoFromConnectorLabel(entry.getKey());
						Set<Integer> ringLinks = ringWithLinks.get(ring);
							if(ringLinks.contains(linkerNo)) {
								connectorAtomsToKeep.add(entry.getValue());
							}
						
					}
					usedConnectorAtoms.addAll(connectorAtomsToKeep);
					ringFormationSynthon.addMolecule(synthon.cutConnectorPaths(synthon.getStructure(), connectorAtomsToKeep));
				}
				for(int connAtom: allConnectorAtoms) {
					if(!usedConnectorAtoms.contains(connAtom)) {
						List<Integer> atomPath = Arrays.asList(connAtom, mol.getConnAtom(connAtom, 0));
						ringFormationSynthon.addMolecule(synthon.cutConnectorPaths(synthon.getStructure(), atomPath));
					}
				}
			}
			ringFormationSynthon.ensureHelperArrays(Molecule.cHelperParities);
			newGenericSynthons.add(new Canonizer(ringFormationSynthon).getCanMolecule());
		
	
		}

		genericSynthons.clear();
		for(StereoMolecule gs : newGenericSynthons) {
			gs.ensureHelperArrays(Molecule.cHelperParities);
			genericSynthons.add(gs);

		}
		rxnHelper = new PharmTreeSynthonReactionHelper(genericSynthons);


		
	}
	


	

	
	

}
