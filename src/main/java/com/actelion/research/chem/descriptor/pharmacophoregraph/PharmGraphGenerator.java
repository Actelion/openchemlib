package com.actelion.research.chem.descriptor.pharmacophoregraph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.stream.Collectors;

import org.openmolecules.chem.conf.gen.ConformerGenerator;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.descriptor.pharmacophoretree.FeatureCalculator;
import com.actelion.research.chem.descriptor.pharmacophoretree.Graph;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreNode;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreTree;
import com.actelion.research.chem.descriptor.pharmacophoretree.PharmacophoreTreeGenerator;
import com.actelion.research.chem.descriptor.pharmacophoretree.TreeUtils;
import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;

public class PharmGraphGenerator {
	private static final int MAX_RING_SIZE = 8;
	private enum PharmAtomElements {ACCEPTOR(85),DONOR(110),POS_CHARGE(78),NEG_CHARGE(41),AROM(47),LIPO(71),AROM_RING(79),
		RING(44),L1(92),L2(93),L3(94),L4(95), RING3(31), RING4(32), RING5(30), RING6(36), RING7(49), RING8(50), HETEO_GROUP(80);
		int id;
		PharmAtomElements(int id) {
			this.id = id;
		}
		int getID() {
			return this.id;
		}
	}
	
	
	public StereoMolecule generateDescriptor(StereoMolecule mol) {
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		mol.stripSmallFragments();
		mol.ensureHelperArrays(Molecule.cHelperNeighbours);
		ConformerGenerator.addHydrogenAtoms(mol);
		List<PharmacophoreNode> ppNodes = new ArrayList<PharmacophoreNode>();
		List<PharmacophoreTree.BiGramInt> nodeEdges = new ArrayList<PharmacophoreTree.BiGramInt>();
		
		Map<Integer,List<Integer>> atomToNodes = new HashMap<Integer,List<Integer>>();
		
		double[] atomVolumes = PharmacophoreTreeGenerator.getAtomVolumes(mol);
		
		List<PharmacophoreNode> ringNodes = new ArrayList<>();
		List<PharmacophoreTree.BiGramInt> ringNodeEdges = new ArrayList<PharmacophoreTree.BiGramInt>();
		List<Set<Integer>> rings = new ArrayList<>();
		Map<Integer,Boolean> aromaticity = new HashMap<Integer,Boolean>();
		FeatureCalculator calculator = new FeatureCalculator(mol);
		calculator.calculate();
		int[][] functionalities = calculator.getAtomFunctionalities();
		processRings(mol, ringNodes, aromaticity, ringNodeEdges,rings, functionalities, atomVolumes);
	
		for(int i=0;i<ringNodes.size();i++) {
			PharmacophoreNode node = ringNodes.get(i);
			PharmacophoreTreeGenerator.addNode(node,ppNodes,atomToNodes);
			
		}

		for(int i=0;i<ringNodeEdges.size();i++) {
			PharmacophoreTree.BiGramInt newEdge = ringNodeEdges.get(i);
			if(!nodeEdges.contains(newEdge))
				nodeEdges.add(newEdge);	
		}
		
		List<PharmacophoreNode> heteroNodes = new ArrayList<>();
		
		processHeteroGroups(mol, heteroNodes, functionalities, atomVolumes);
		
		for(int i=0;i<heteroNodes.size();i++) {
			PharmacophoreNode node = heteroNodes.get(i);
			PharmacophoreTreeGenerator.addNode(node,ppNodes,atomToNodes);
			
		}

		treeWalk(0,mol,atomToNodes,ppNodes,nodeEdges,functionalities, atomVolumes);


		Map<Integer,Map<Integer,Integer>> adj = TreeUtils.getAdjacencyListWithBondOrders(ppNodes.size(), nodeEdges);
		
		StereoMolecule pharmMol = createPharmMol(ppNodes,adj);
		
		//System.out.println(pharmMol.getIDCode());
		
		return pharmMol;

		
		
		
	}
	
	
	private static StereoMolecule createPharmMol(List<PharmacophoreNode> ppNodes,
			Map<Integer,Map<Integer,Integer>> adjList) {
		StereoMolecule pharmMol = new StereoMolecule();
		Map<Integer,Integer> nodesToGraphAtoms = new HashMap<Integer,Integer>();
		for(int n=0;n<ppNodes.size();n++) {
			int atomId = -1;
			if(!nodesToGraphAtoms.keySet().contains(n)) {
				atomId = addNodeToPharmGraph(pharmMol,ppNodes,n,nodesToGraphAtoms);
			}
			else {//already added --> get id
				atomId = nodesToGraphAtoms.get(n);
			}
			Map<Integer,Integer> neighbours = adjList.get(n);
			for(int nn : neighbours.keySet()) {
				int atomId2 = -1;
				if(nodesToGraphAtoms.keySet().contains(nn))
					atomId2 = nodesToGraphAtoms.get(nn);
				else 
					atomId2 = addNodeToPharmGraph(pharmMol,ppNodes,nn,nodesToGraphAtoms);
				int bondOrder = neighbours.get(nn);
				pharmMol.ensureHelperArrays(Molecule.cHelperNeighbours);
				if(pharmMol.getBond(atomId, atomId2)==-1)	{
					if(bondOrder==2)
						pharmMol.addBond(atomId, atomId2, Molecule.cBondTypeDouble);
					else if(bondOrder==1)
						pharmMol.addBond(atomId, atomId2, Molecule.cBondTypeSingle);

				}
				
			}
			

		}
		return pharmMol;
		
	}
	
	private static int addNodeToPharmGraph(StereoMolecule pharmGraph, List<PharmacophoreNode> nodes, int nodeIndex, Map<Integer,Integer> nodesToGraphAtoms) {
		int atomID = -1;
		PharmacophoreNode node = nodes.get(nodeIndex);
		int acceptors = node.getFunctionalities()[PharmacophoreCalculator.ACCEPTOR_ID];
		int donors = node.getFunctionalities()[PharmacophoreCalculator.DONOR_ID];
		int negCharge = node.getFunctionalities()[PharmacophoreCalculator.CHARGE_NEG_ID];
		int posCharge = node.getFunctionalities()[PharmacophoreCalculator.CHARGE_POS_ID];
		int totPolar = acceptors+donors+negCharge+posCharge;
		if (!node.isRing() && totPolar>0) {
			atomID = pharmGraph.addAtom(PharmAtomElements.HETEO_GROUP.id);
			nodesToGraphAtoms.put(nodeIndex, atomID);
			addPolarInfo(pharmGraph, atomID, acceptors, donors, negCharge, posCharge);
		}
		else if(node.isAromatic()) {
			atomID = pharmGraph.addAtom(PharmAtomElements.AROM_RING.id);
			nodesToGraphAtoms.put(nodeIndex, atomID);
			addRingSizeInfo(pharmGraph,atomID,node.getAtoms().size());
			if(totPolar>0)
				addPolarInfo(pharmGraph, atomID, acceptors, donors, negCharge, posCharge);
		}
		else if(node.isRing()) {
			atomID = pharmGraph.addAtom(PharmAtomElements.RING.id);
			nodesToGraphAtoms.put(nodeIndex, atomID);
			addRingSizeInfo(pharmGraph,atomID,node.getAtoms().size());
			if(totPolar>0)
				addPolarInfo(pharmGraph, atomID, acceptors, donors, negCharge, posCharge);
		}
		else if(node.isLinkNode()) {
			int linkerID = node.getFunctionalities()[0];
			switch(linkerID) {
				case 1:
					atomID = pharmGraph.addAtom(PharmAtomElements.L1.id);
					nodesToGraphAtoms.put(nodeIndex, atomID);
					break;
				case 2:
					atomID = pharmGraph.addAtom(PharmAtomElements.L2.id);
					nodesToGraphAtoms.put(nodeIndex, atomID);
					break;
				case 3:
					atomID = pharmGraph.addAtom(PharmAtomElements.L3.id);
					nodesToGraphAtoms.put(nodeIndex, atomID);
					break;
				case 4:
					atomID = pharmGraph.addAtom(PharmAtomElements.L4.id);
					nodesToGraphAtoms.put(nodeIndex, atomID);
					break;
			}
			
		}

		else {
			atomID = pharmGraph.addAtom(PharmAtomElements.LIPO.id);
			nodesToGraphAtoms.put(nodeIndex, atomID);
		}
		return atomID;
	}
	
	private static void addRingSizeInfo(StereoMolecule pharmGraph, int center, int ringSize) {
		pharmGraph.ensureHelperArrays(Molecule.cHelperNeighbours);
		switch(ringSize) {
			case 3: 
				int id = pharmGraph.addAtom(PharmAtomElements.RING3.id);
				pharmGraph.addBond(center, id, Molecule.cBondTypeMetalLigand);
				break;
			case 4: 
				id = pharmGraph.addAtom(PharmAtomElements.RING4.id);
				pharmGraph.addBond(center, id, Molecule.cBondTypeMetalLigand);
				break;
			case 5: 
				id = pharmGraph.addAtom(PharmAtomElements.RING5.id);
				pharmGraph.addBond(center, id, Molecule.cBondTypeMetalLigand);
				break;
			case 6: 
				id = pharmGraph.addAtom(PharmAtomElements.RING6.id);
				pharmGraph.addBond(center, id, Molecule.cBondTypeMetalLigand);
				break;
			case 7: 
				id = pharmGraph.addAtom(PharmAtomElements.RING7.id);
				pharmGraph.addBond(center, id, Molecule.cBondTypeMetalLigand);
				break;
			case 8: 
				id = pharmGraph.addAtom(PharmAtomElements.RING8.id);
				pharmGraph.addBond(center, id, Molecule.cBondTypeMetalLigand);
				break;
		}
		pharmGraph.ensureHelperArrays(Molecule.cHelperNeighbours);
	}
	
	private static void addPolarInfo(StereoMolecule pharmGraph, int center, int acceptors, int donors, int chargeNeg, 
			int chargePos) {
		pharmGraph.ensureHelperArrays(Molecule.cHelperNeighbours);
		for(int i=0;i<acceptors;i++) {
			int id = pharmGraph.addAtom(PharmAtomElements.ACCEPTOR.id);
			pharmGraph.addBond(center, id, Molecule.cBondTypeMetalLigand);
		}
		for(int i=0;i<donors;i++) {
			int id = pharmGraph.addAtom(PharmAtomElements.DONOR.id);
			pharmGraph.addBond(center, id, Molecule.cBondTypeMetalLigand);
		}
		for(int i=0;i<chargeNeg;i++) {
			int id = pharmGraph.addAtom(PharmAtomElements.NEG_CHARGE.id);
			pharmGraph.addBond(center, id, Molecule.cBondTypeMetalLigand);
		}
		for(int i=0;i<chargePos;i++) {
			int id = pharmGraph.addAtom(PharmAtomElements.POS_CHARGE.id);
			pharmGraph.addBond(center, id, Molecule.cBondTypeMetalLigand);
		}
		pharmGraph.ensureHelperArrays(Molecule.cHelperNeighbours);
	}
	
	private static void processRings(StereoMolecule mol, List<PharmacophoreNode> ringNodes, Map<Integer,Boolean> aromaticity, List<PharmacophoreTree.BiGramInt> ringNodeEdges,
			List<Set<Integer>> rings, int[][] functionalities, double[] atomVolumes) {
		RingCollection rc = new RingCollection(mol,RingCollection.MODE_SMALL_AND_LARGE_RINGS_AND_AROMATICITY, MAX_RING_SIZE);
		for(int atom=0;atom<mol.getAtoms();atom++) {
			//if(!mol.isRingAtom(atom))
			//	continue;
			//else {
			List<Integer> smallestRings = PharmacophoreTreeGenerator.getSmallestRingsOfAtom(rc,atom);
			for(Integer smallRingIndex: smallestRings) {
				Set<Integer> smallRing = Arrays.stream(rc.getRingAtoms(smallRingIndex)).boxed().collect(Collectors.toSet());
				if(!rings.contains(smallRing)) {
					rings.add(smallRing);
					boolean isAromatic = rc.isAromatic(smallRingIndex);
					aromaticity.put(rings.indexOf(smallRing), isAromatic);
				}
					
				}
				//}
		}
		for(int ringIndex=0;ringIndex<rings.size();ringIndex++) {
			boolean isAromatic = aromaticity.get(ringIndex);
			Set<Integer> ring = rings.get(ringIndex);
			PharmacophoreNode node = new PharmacophoreNode(new ArrayList<Integer>(ring),functionalities, atomVolumes, true,isAromatic);
			ringNodes.add(node);
		}
		for(int r1=0;r1<rings.size();r1++) {
			Set<Integer> ring1 = rings.get(r1);
			for(int r2=r1+1;r2<rings.size();r2++) {
				Set<Integer> ring2 = rings.get(r2);
				long commonElements = ring1.stream().filter(e -> ring2.contains(e)).count();
				if(commonElements>0) {
					ringNodeEdges.add(new PharmacophoreTree.BiGramInt(new int[] {r1,r2},2));
				}
			}
		}
		
			
	}
	
	/**
	 * merge atoms with terminal polar substituents into hetero nodes
	 * @param mol
	 * @param heteroNodes
	 */
	private static void processHeteroGroups(StereoMolecule mol, List<PharmacophoreNode> heteroNodes,
			int[][] functionalities, double[] atomVolumes) {
		for(int atom=0;atom<mol.getAtoms();atom++) {
			if(mol.isRingAtom(atom))
				continue;
			Set<Integer> mergedAtoms = new HashSet<Integer>();
			mergedAtoms.add(atom);
			for(int n=0;n<mol.getConnAtoms(atom);n++) {
				int aa = mol.getConnAtom(atom, n);
				if(mol.getAtomicNo(aa)==7 || mol.getAtomicNo(aa)==8)
				if(mol.getConnAtoms(aa)==1)
					mergedAtoms.add(aa);
			}
			if(mergedAtoms.size()>1) {
				PharmacophoreNode node = new PharmacophoreNode(new ArrayList<Integer>(mergedAtoms),functionalities, atomVolumes, false, false);
				heteroNodes.add(node);
			}
		}
		
	}
	

	
	/**
	 * Walk through the molecule and create nodes, edges and update the atom-to-nodes mapping.
	 * @param atom
	 * @param mol
	 * @param atomToNodes
	 * @param ppNodes
	 * @param nodeEdges
	 * @param functionalities
	 * @param atomVolumes
	 */
	private static void treeWalk(int atom, StereoMolecule mol,Map<Integer,List<Integer>> atomToNodes,List<PharmacophoreNode> ppNodes, List<PharmacophoreTree.BiGramInt> nodeEdges,
			int[][] functionalities, double[] atomVolumes) {
		boolean[] visited = new boolean[mol.getAtoms()];
		int[] parents = new int[mol.getAtoms()];
		Queue<Integer> atoms = new LinkedList<Integer> ();
		atoms.add(atom);
		parents[atom] = -1; //root does not have parents
		visited[atom] = true;
		while(!atoms.isEmpty()) {
			int currentAtom  = atoms.poll();
			visited[currentAtom] = true;
			if(atomToNodes.get(currentAtom)==null) {
					String label = mol.getAtomLabel(currentAtom);
					PharmacophoreNode node;
					if(PharmacophoreTreeGenerator.RGROUPS.contains(label)) { //is a link atom
						node = new PharmacophoreNode(Arrays.asList(currentAtom),functionalities,atomVolumes,7,false,false);
						node.getFunctionalities()[0] = Integer.parseInt(label.split("R")[1]);
					}
					else { //not a link, check functional groups;
						
						node = new PharmacophoreNode(Arrays.asList(currentAtom),functionalities,atomVolumes,false,false);
					}
					PharmacophoreTreeGenerator.addNode(node, ppNodes, atomToNodes);
					int index = ppNodes.size()-1;
					if(parents[currentAtom]!=-1) {
						List<Integer> connNodes = atomToNodes.get(parents[currentAtom]);
							for(int connNode : connNodes) {
								PharmacophoreTree.BiGramInt newEdge = new PharmacophoreTree.BiGramInt(new int[] {connNode,index});
							if(!nodeEdges.contains(newEdge))
								nodeEdges.add(newEdge);
							}
					}
			}
			
			else { 
				List<Integer> nodes = atomToNodes.get(currentAtom); //get all nodes this atom is part of
				if(parents[currentAtom]!=-1) { //not root
					List<Integer> connNodes = atomToNodes.get(parents[currentAtom]); //get nodes that the parent atom is part of
					for(int node : nodes) {
						for(int connNode : connNodes) {
							if(connNode!=node) {
								PharmacophoreTree.BiGramInt newEdge = new PharmacophoreTree.BiGramInt(new int[] {connNode,node});
								if(!nodeEdges.contains(newEdge))
									nodeEdges.add(newEdge);
							}
						}
					}
				}
				
			}
			
			for(int a=0;a<mol.getConnAtoms(currentAtom);a++) {
				int nextAtom = mol.getConnAtom(currentAtom, a);
				if(visited[nextAtom])
					continue;
				else {
					atoms.add(nextAtom);
					parents[nextAtom] = currentAtom;
				}
			}
		}
		
		
	}

}
