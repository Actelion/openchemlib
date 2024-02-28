package com.actelion.research.chem.descriptor.pharmacophoretree;
import java.util.ArrayList;
import java.util.Arrays;
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
import com.actelion.research.chem.conf.BondLengthSet;
import com.actelion.research.chem.conf.VDWRadii;
import com.actelion.research.chem.RingCollection;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.chemicalspaces.synthon.SynthonReactor;
import com.actelion.research.chem.phesa.pharmacophore.IonizableGroupDetector;
import com.actelion.research.chem.phesa.pharmacophore.PharmacophoreCalculator;
import com.actelion.research.chem.phesa.pharmacophore.pp.IPharmacophorePoint;


/**
 * Generates a PharmacophoreTree from a StereoMolecule. Calculates chemical features and steric properties of the nodes spanning the tree.
 * Cycles are resolved by merging nodes and introducing zero-nodes. 
 * @author joel
 *
 */


public class PharmacophoreTreeGenerator { 
	
	public static final int MAX_RING_SIZE = 50;
	
	public static final Set<String> RGROUPS = new HashSet<String>(Arrays.asList("U", "Np", "Pu", "Am"));
	

	private PharmacophoreTreeGenerator() {}
	
	public static PharmacophoreTree generate(StereoMolecule mol) {
		return generate(mol,new HashMap<Integer,List<Integer>>(),new ArrayList<Set<Integer>>());
	}
	
	

	/**
	 * the parameter rings can be used for fragments/building blocks to submit atoms that will belong to a ring after they react
	 * @param mol
	 * @param atomToNodes
	 * @param rings
	 * @return
	 */
	
	public static PharmacophoreTree generate(StereoMolecule mol,Map<Integer,List<Integer>> atomToNodes, List<Set<Integer>> rings) {
		mol.ensureHelperArrays(Molecule.cHelperCIP);
		ConformerGenerator.addHydrogenAtoms(mol);
		mol.ensureHelperArrays(Molecule.cHelperCIP);

		double[] atomVolumes = getAtomVolumes(mol);
		
		List<PharmacophoreNode> ppNodes = new ArrayList<PharmacophoreNode>();
		List<PharmacophoreTree.BiGramInt> nodeEdges = new ArrayList<PharmacophoreTree.BiGramInt>();
		
		List<Set<Integer>> ringNodes = new ArrayList<Set<Integer>>();
		List<int[]> ringNodeEdges = new ArrayList<int[]>();
		createRingGraphs(mol, ringNodes, ringNodeEdges,rings);
		FeatureCalculator calculator = new FeatureCalculator(mol);
		calculator.calculate();
		int[][] functionalities = calculator.getAtomFunctionalities();
		
		for(int i=0;i<ringNodes.size();i++) {
			Set<Integer> ring = ringNodes.get(i);
			PharmacophoreNode node = new PharmacophoreNode(new ArrayList<Integer>(ring),functionalities, atomVolumes,false,false);
			addNode(node,ppNodes,atomToNodes);
			
		}

		for(int i=0;i<ringNodeEdges.size();i++) {
			PharmacophoreTree.BiGramInt newEdge = new PharmacophoreTree.BiGramInt(ringNodeEdges.get(i));
			if(!nodeEdges.contains(newEdge))
				nodeEdges.add(newEdge);	
		}

		treeWalk(0,mol,atomToNodes,ppNodes,nodeEdges,functionalities, atomVolumes);


		Map<Integer,List<Integer>> adj = TreeUtils.getAdjacencyList(ppNodes.size(), nodeEdges.stream().map(e -> e.edge).collect(Collectors.toList()));

		Graph g = new Graph(adj);
		//resolve remaining cycles by introducing zero-nodes
		Set<PharmacophoreTree.BiGramInt> edgesToDelete = new HashSet<PharmacophoreTree.BiGramInt>();
		List<List<int[]>> bccs = g.bcc();
		for(List<int[]> bcc : bccs) {
			if(bcc.size()>1) {
				Set<Integer> bccNodes = new HashSet<Integer>();
				for(int[] c : bcc ) {
					bccNodes.add(c[0]);
					bccNodes.add(c[1]);
					edgesToDelete.add(new PharmacophoreTree.BiGramInt(c));
				}
				PharmacophoreNode node = new PharmacophoreNode(new ArrayList<Integer>(), functionalities,
					atomVolumes,1,false,false);
				addNode(node,ppNodes,atomToNodes);
				int index = ppNodes.size()-1;
				for(Integer connNode : bccNodes) {
					nodeEdges.add(new PharmacophoreTree.BiGramInt(new int[] {connNode,index}));
				}
		}
		}
		nodeEdges.removeAll(edgesToDelete);
		for(PharmacophoreNode ppNode : ppNodes)
			ppNode.updateWeights(atomToNodes);
		return new PharmacophoreTree(ppNodes,nodeEdges.stream().map(e -> e.edge).collect(Collectors.toList()));
		
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
		parents[atom] = -1;
		visited[atom] = true;
		while(!atoms.isEmpty()) {
			int currentAtom  = atoms.poll();
			visited[currentAtom] = true;
			if(!mol.isRingAtom(currentAtom) && atomToNodes.get(currentAtom)==null) {
					String label = mol.getAtomLabel(currentAtom);
					PharmacophoreNode node;
					if(RGROUPS.contains(label)) { //is a link atom
						node = new PharmacophoreNode(Arrays.asList(currentAtom),functionalities,atomVolumes,7,false,false);
						node.getFunctionalities()[0] = StereoMolecule.getAtomicNoFromLabel(label)+1-SynthonReactor.CONNECTOR_OFFSET;
					}
					else 
						node = new PharmacophoreNode(Arrays.asList(currentAtom),functionalities,atomVolumes,false,false);
					addNode(node, ppNodes, atomToNodes);
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
			else { //is ringAtom 
				List<Integer> nodes = atomToNodes.get(currentAtom);
				if(parents[currentAtom]!=-1) {
					List<Integer> connNodes = atomToNodes.get(parents[currentAtom]);
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
	
	public static void addNode(PharmacophoreNode node, List<PharmacophoreNode> nodes,Map<Integer,List<Integer>> atomToNodes) {
		int index = nodes.size();
		nodes.add(node);
		for(int atom:node.getAtoms()) {
			atomToNodes.putIfAbsent(atom, new ArrayList<Integer>());
			atomToNodes.get(atom).add(index);
		}
	}
	
	
	/**
	 * ring graphs are created as follows: for every atom in a ring, the size s_min of the smallest ring containing this
	 * atom is determined. All rings with size s_min containing this atom are then evaluated and added to the set of rings.
	 * Rings that share atoms are subsequently connected through an edge. In order to convert the ring graph effectively to a tree,
	 * ring systems that form a cyclic graph of nodes are merged into a single node. Cyclic systems are detected using a bicyclic-component algorithm.
	 * @param mol
	 * @param nodes
	 * @param edges
	 */
	
	private static void createRingGraphs(StereoMolecule mol, List<Set<Integer>> ringNodes, List<int[]> ringNodeEdges,
			List<Set<Integer>> rings) {
		RingCollection rc = new RingCollection(mol,RingCollection.MODE_SMALL_AND_LARGE_RINGS,MAX_RING_SIZE);
		List<int[]> ringEdges = new ArrayList<int[]>();
		for(int atom=0;atom<mol.getAtoms();atom++) {
			//if(!mol.isRingAtom(atom))
			//	continue;
			//else {
				List<Integer> smallestRings = getSmallestRingsOfAtom(rc,atom);
				for(Integer smallRingIndex: smallestRings) {
					Set<Integer> smallRing = Arrays.stream(rc.getRingAtoms(smallRingIndex)).boxed().collect(Collectors.toSet());
					if(!rings.contains(smallRing))
						rings.add(smallRing);
						
					}
				//}
		}
		for(int r1=0;r1<rings.size();r1++) {
			Set<Integer> ring1 = rings.get(r1);
			for(int r2=r1+1;r2<rings.size();r2++) {
				Set<Integer> ring2 = rings.get(r2);
				long commonElements = ring1.stream().filter(e -> ring2.contains(e)).count();
				if(commonElements>0) {
					ringEdges.add(new int[] {r1,r2});
				}
			}
		}


		Map<Integer,List<Integer>> adj = TreeUtils.getAdjacencyList(rings.size(), ringEdges);

		Graph g = new Graph(adj);
		
		List<List<int[]>> bccs = g.bcc();
		
		// merge ring systems that form biconnected components
		List<Set<Integer>> toDelete = new ArrayList<Set<Integer>>();
		for(List<int[]> bcc : bccs) {
			if (bcc.size()>1) {
				Set<Integer> ring = new HashSet<Integer>();
				for(int[] edge : bcc) {
					ring.addAll(rings.get(edge[0]));
					ring.addAll(rings.get(edge[1]));
					toDelete.add(rings.get(edge[0]));
					toDelete.add(rings.get(edge[1]));
				}
				ringNodes.add(ring);				
			}
		}
		rings.removeAll(toDelete);
		for(Set<Integer> ring : rings) {
			ringNodes.add(ring);
		}


		for(int r1=0;r1<ringNodes.size();r1++) {
			Set<Integer> ring1 = ringNodes.get(r1);
			for(int r2=r1+1;r2<ringNodes.size();r2++) {
				Set<Integer> ring2 = ringNodes.get(r2);
				long commonElements = ring1.stream().filter(e -> ring2.contains(e)).count();
				if(commonElements>0) {
					ringNodeEdges.add(new int[] {r1,r2});
				}
			}
		}

		
		
		
		 
		
			
	}
	
	public static List<Integer> getSmallestRingsOfAtom(RingCollection rc, int atom){
		List<Integer> rings = new ArrayList<Integer>();
		int smallestRingSize = rc.getAtomRingSize(atom);
		for(int r=0;r<rc.getSize();r++) {
			Set<Integer> ringAtoms = Arrays.stream(rc.getRingAtoms(r)).boxed().collect(Collectors.toSet());
			if(ringAtoms.size()>smallestRingSize)
				continue;
			if(ringAtoms.contains(atom)) 
				rings.add(r);
				
			
		}
				
				
		return rings;
		}
	

	
	public static double[] getAtomVolumes(StereoMolecule mol) {
		double[] atomVolumes = new double[mol.getAtoms()];
		for(int a=0;a<mol.getAtoms();a++) {
			atomVolumes[a] = getAtomicVdWVolume(mol.getAtomicNo(a));
		}
		for(int b=0;b<mol.getBonds();b++) {
			int a1 = mol.getBondAtom(0, b);
			int a2 = mol.getBondAtom(1, b);
			double[] capCorrection = calculateCapVolumesIntersectingAtomSpheres(mol,a1,a2);
			atomVolumes[a1]-=capCorrection[0];
			atomVolumes[a2]-=capCorrection[1];
		}
		return atomVolumes;
		
		
	}
	
	private static double getAtomicVdWVolume(int atomicNo) {
		double r=0;
		if(atomicNo< VDWRadii.VDW_RADIUS.length)
			r = VDWRadii.VDW_RADIUS[atomicNo];
		return (4.0/3)*Math.PI*Math.pow(r, 3);
	}
	
	/**
	 * sphere-sphere intersection: http://mathworld.wolfram.com/Sphere-SphereIntersection.html
	 * @param mol
	 * @param a
	 * @param b
	 * @return
	 */
	private static double[] calculateCapVolumesIntersectingAtomSpheres(StereoMolecule mol, int a, int b) {
		double r1 = 0.0;
		double r2 = 0.0;
		double v1 = 0.0;
		double v2 = 0.0;
		if(mol.getAtomicNo(a)<VDWRadii.VDW_RADIUS.length && mol.getAtomicNo(b)<VDWRadii.VDW_RADIUS.length) {
			r1 = VDWRadii.VDW_RADIUS[mol.getAtomicNo(a)];
			r2 = VDWRadii.VDW_RADIUS[mol.getAtomicNo(b)];
			int bond = mol.getBond(a, b);
			double d = BondLengthSet.lookupBondLength(mol, bond);
			//x coordinate of intersection:
			double x = (d*d-r1*r1+r2*r2)/(2*d);
			// distance of sphere centers to cap bases of intersection:
			double d1 = x;
			double d2 = d-x;
			//height of caps:
			double h1 = r1-d1;
			double h2 = r2-d2;
			// calculate volumes of the caps:
			v1 = (1.0/3)*Math.PI*h1*h1+(3*r1-h1);
			v2 = (1.0/3)*Math.PI*h2*h2+(3*r2-h2);
		}
		
		return new double[] {v1,v2};
		
	}

	


}
