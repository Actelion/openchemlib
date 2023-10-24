package com.actelion.research.chem.descriptor.pharmacophoretree;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.stream.Collectors;

import com.actelion.research.chem.StereoMolecule;


/**
 * A PharmacophoreTree closely follows the concepts of FeatureTrees as described by Rarey and Dixon
 * (DOI:10.1023/a:1008068904628). A PharmacophoreTree consists of connected nodes and it's data structure 
 * is an undirected tree. No cycles are allowed which distinguishes it from a molecular graph. 
 * PharmacophoreTrees facilitates the comparison of chemical structures with respect to their features in 2D. 
 * @author joel
 *
 */



public class PharmacophoreTree {
	
	public static final int CUT_NONE = 0;
	public static final int CUT_RIGHT = 1;
	public static final int CUT_LEFT = -1;
	
	public static final int MAX_EXTENSION_CUTS = 500;

	private List<int[]> edges;
	private List<PharmacophoreNode> nodes;
	private Map<Integer,List<Integer>> adjacencyList;
	private int linkerNodes;

	
	public PharmacophoreTree(List<PharmacophoreNode> nodes, List<int[]> edges) {
		this.nodes = nodes;
		this.edges = edges;
		linkerNodes = 0;
		update();
	}
	
	private void update() {
		linkerNodes = 0;
		for(PharmacophoreNode node : nodes) {
			if(node.isLinkNode())
				linkerNodes++;
		}
		adjacencyList = new HashMap<Integer,List<Integer>>();
		for(int i=0;i<nodes.size();i++) {
			adjacencyList.putIfAbsent(i, new ArrayList<Integer>());
			for(int j=0;j<edges.size();j++) {
				int[] edge = edges.get(j);
				if(edge[0]==i || edge[1] == i)
					adjacencyList.get(i).add(j);
			}
		}
		
	}
	


	/**
	 * 
	 * a cut divides a edge into a target node and a source node: 
	 * an edge is defined as a tuple of two nodes (a,b), a left cut results in node a
	 * being the source node and b being the target node: a-->--b 
	 * a right cut results in the right node being the source node:
	 *  a--<--b
	 *  returns int array with source node as 0th element and target node as 1st element
	 *
	 * @param cut
	 * @param edge
	 * @param sourceTreeEdges
	 * @param sourceTreeEdgeParents
	 * @param targetTreeEdges
	 * @param targetTreeEdgeParents
	 * @return
	 */
	public int[] initialCut(int cut, int edge, List<Integer> sourceTreeEdges,List<Integer> sourceTreeEdgeParents, 
			List<Integer> targetTreeEdges, List<Integer> targetTreeEdgeParents) {
		int sourceNode = -1;
		int targetNode = -1;
		if(cut==CUT_LEFT) {
			sourceNode = edges.get(edge)[0];
			targetNode = edges.get(edge)[1];
		}
		else if (cut==CUT_RIGHT) {
			sourceNode = edges.get(edge)[1];
			targetNode = edges.get(edge)[0];
		}
		else {
			throw new IllegalArgumentException();
		}

		treeWalkBFS(sourceNode,edge,sourceTreeEdges,sourceTreeEdgeParents);
		treeWalkBFS(targetNode,edge,targetTreeEdges,targetTreeEdgeParents);
		
		return new int[] {sourceNode, targetNode};
	}
	
	
	
	/**
	 * Walks a subtree of the PharmacophoreTree in breadth-first manner starting from a deleted edge and a designated
	 * head node. Returns a list of edge indeces and the parents of the edges. 
	 * @param headNode
	 * @param deletedEdgeIndex
	 * @param treeEdgesBFS
	 * @param edgeParentsBFS
	 */
	public void treeWalkBFS(int headNode,int deletedEdgeIndex,List<Integer> treeEdgesBFS, List<Integer> edgeParentsBFS) {

		Set<Integer> visitedNodes = new HashSet<Integer>();
		Set<Integer> visitedEdges = new HashSet<Integer>();
		Queue<Integer> pqNodes = new LinkedList<Integer> ();
		Map<Integer,Integer> parentEdges = new HashMap<Integer,Integer>();
		pqNodes.add(headNode);
		parentEdges.put(headNode,deletedEdgeIndex);
		visitedEdges.add(deletedEdgeIndex);
		treeEdgesBFS.add(deletedEdgeIndex);
		edgeParentsBFS.add(-1);
		while(!pqNodes.isEmpty()) {
			int node = pqNodes.poll();
			visitedNodes.add(node);
			int parentEdge = parentEdges.get(node);
			for(int edgeNo : adjacencyList.get(node)) {
				int[] edge = edges.get(edgeNo);

				if(visitedEdges.contains(edgeNo))
					continue;
				visitedEdges.add(edgeNo);
				int nextNode = -1;
				if(edge[0] == node) 
					nextNode = edge[1];
				else if(edge[1] == node) 
					nextNode = edge[0];
				pqNodes.add(nextNode);
				parentEdges.put(nextNode,edgeNo);
				treeEdgesBFS.add(edgeNo);
				edgeParentsBFS.add(parentEdge);
			
			}
		}
			
	}
	
	/**
	 * 	
	 * retrieve nodes that are part of extension match ("extension nodes") as well as the nodes that
	 * are part of a cut subtree (source nodes) from a cut-string
	 *
	 * @param head
	 * @param cut
	 * @param subTreeEdgeIndeces
	 * @param extensionNodes
	 * @param sourceNodes
	 */
	public void enumerateExtensionCutFast(int head,int[] cut, List<Integer> subTreeEdgeIndeces,
			Set<Integer> extensionNodes, Set<Integer> sourceNodes) {
		for(int i=0;i<cut.length;i++) {
			int[] edge = edges.get(subTreeEdgeIndeces.get(i));
			if(cut[i]==0) {
				extensionNodes.add(edge[0]);
				extensionNodes.add(edge[1]);
			}
			else if(cut[i]==1) {
				if(extensionNodes.contains(edge[0]))
					sourceNodes.add(edge[1]);
				else 
					sourceNodes.add(edge[0]);
			}
			else if(cut[i]==-1) {
				sourceNodes.add(edge[0]);
				sourceNodes.add(edge[1]);
			}
				
			}
	}
	
/**
 * retrieves nodes that are part of extension match ("extension nodes") as well as the different subtrees
 * resulting from the cuts with a list of their edges (in bfs order) as well as the parents of the edges
 * @param head
 * @param cut
 * @param subTreeEdgeIndeces
 * @param subTreeParentEdgeIndeces
 * @param sourceTreeEdgeIndeces
 * @param sourceTreeEdgeParentIndeces
 * @param sourceTreeHeadNodes
 * @param extensionNodes
 * @param cutEdges
 * @param cutDirections
 */
	
	public void enumerateExtensionCutFull(int head,int[] cut, List<Integer> subTreeEdgeIndeces,
			List<Integer> subTreeParentEdgeIndeces, List<List<Integer>> sourceTreeEdgeIndeces,List<List<Integer>> sourceTreeEdgeParentIndeces, 
			List<Integer> sourceTreeHeadNodes,Set<Integer> extensionNodes, List<Integer> cutEdges, List<Integer> cutDirections) {
		extensionNodes.add(head);
		for(int i=1;i<cut.length;i++) { //first edge is the cut edge, don't assign full edge to extension cut
			int edgeIndex = subTreeEdgeIndeces.get(i);
			int[] edge = edges.get(edgeIndex);
			if(cut[i]==0) {
				extensionNodes.add(edge[0]);
				extensionNodes.add(edge[1]);
			}
			else if(cut[i]==1) { // a cut node serves as a parent node for a subtree
				cutEdges.add(edgeIndex);
				if(extensionNodes.contains(edge[0])) { // right-cut
					sourceTreeHeadNodes.add(edge[1]);
					cutDirections.add(CUT_RIGHT);
					
				}
				else {
					sourceTreeHeadNodes.add(edge[0]);
					cutDirections.add(CUT_LEFT);
				}
				List<Integer> edgeIndeces = new ArrayList<Integer>();
				edgeIndeces.add(edgeIndex);
				sourceTreeEdgeIndeces.add(edgeIndeces);
				List<Integer> parentEdgeIndeces = new ArrayList<Integer>();
				parentEdgeIndeces.add(-1);
				sourceTreeEdgeParentIndeces.add(parentEdgeIndeces);
			}
			else if(cut[i]==-1) {
				int parentEdgeIndex = subTreeParentEdgeIndeces.get(i);
				int subTreeIndex = cutEdges.indexOf(parentEdgeIndex);
				if(subTreeIndex<0) {
					for(int j=0;j<sourceTreeEdgeIndeces.size();j++) {
						if(sourceTreeEdgeIndeces.get(j).contains(parentEdgeIndex)) {
							subTreeIndex = j;
							break;
						}
					}
					
				}
				sourceTreeEdgeIndeces.get(subTreeIndex).add(edgeIndex);
				sourceTreeEdgeParentIndeces.get(subTreeIndex).add(parentEdgeIndex);
				
			}	
			}
	}
	
	/**
	 * get a list of int[] arrays that define an extension cut:
	 * Given a subtree with a designated head-node, an extension cut separates an extension match from
	 * the remaining subtrees. The extension cut defines for every edge of the input subtrees their 
	 * status: 0 -> edge contained in extension match   1 -> edge is cut   -1 -> edge is part of a cut subtree
	 * @param subtreeEdgeIndeces
	 * @param subtreeEdgeParentIndeces
	 * @return
	 */
	
	public List<int[]> getExtensionCuts(List<Integer> subtreeEdgeIndeces, 
			List<Integer> subtreeEdgeParentIndeces) {
		List<int[]> cuts = new ArrayList<int[]> ();
		int lowerBound = 0;
		int[] previousCut = new int[subtreeEdgeIndeces.size()];
		Arrays.fill(previousCut, -1);
		previousCut[0] = 1;
		while(lowerBound>=0) {
			if(cuts.size()>MAX_EXTENSION_CUTS) // prevent combinatorial explosion
				return cuts;
			int[] nextCut = previousCut.clone();
			lowerBound = getNextCut(previousCut,nextCut,subtreeEdgeIndeces,subtreeEdgeParentIndeces);
			previousCut = nextCut;
			boolean alreadyPresent = false;
			for(int i=0;i<cuts.size();i++) {
				if(Arrays.equals(previousCut, cuts.get(i)))
					alreadyPresent = true;
			}
			if(!alreadyPresent)
				cuts.add(previousCut);
			
			
		}

		return cuts;
		
		
	}
	
	public Set<Integer> getNodesFromEdges(List<Integer> edgeIndeces) {
		Set<Integer> nodes = new HashSet<Integer>();
		for(int i=1;i<edgeIndeces.size();i++) { // skip first node, since it is the cut node (its head node is added later)
			int edgeIndex =edgeIndeces.get(i);
			int[] edge = edges.get(edgeIndex);
			nodes.add(edge[0]);
			nodes.add(edge[1]);
		}
		return nodes;
	}
	
	
	
	/**
	 * get an array representing the next extension-cut from a previous one
	 * adapt from: 
	 * DOI: 10.1023/a:1008068904628 "Fig. 25 INCREASE_CUTSTRING"
	 * @param previousCut
	 * @param nextCut
	 * @param subtreeEdgeIndeces
	 * @param subtreeEdgeParentIndeces
	 * @return
	 */
	
	private int getNextCut(int[] previousCut, int[] nextCut, List<Integer> subtreeEdgeIndeces, 
			List<Integer> subtreeEdgeParentIndeces) {

		int lowerBound = -1;
		int head = Integer.MAX_VALUE;
		int tail = Integer.MIN_VALUE;
		for(int i=0;i<previousCut.length;i++) {
			if(previousCut[i]==1) {
				if(i<head) 
					head=i;
				if(i>tail)
					tail = i;
			}
		}
		if(head<Integer.MAX_VALUE && tail>Integer.MIN_VALUE) {
			Set<Integer> childNodes = new HashSet<Integer>();
			for(int i=0;i<head;i++) {
				childNodes.add(edges.get(subtreeEdgeIndeces.get(i))[0]);
				childNodes.add(edges.get(subtreeEdgeIndeces.get(i))[1]);
			}
			lowerBound = childNodes.size();
			if(lowerBound>TreeMatcher.MATCH_NODE_NR_LIMIT)  //extension match too large
				lowerBound = -1;
			else {
				nextCut[tail] = 0;
				for(int i=tail+1;i<nextCut.length;i++) {
					int parent = subtreeEdgeParentIndeces.get(i);
					int parentIndexInCutArray = subtreeEdgeIndeces.indexOf(parent);
					if(parentIndexInCutArray==tail)
						nextCut[i] = 1;
					else if(nextCut[parentIndexInCutArray]==0)
						nextCut[i] = 1;
					else
						nextCut[i] = -1;
				}
			}
			
		}

		return lowerBound;
		
	}
	
	public List<PharmacophoreNode> getNodes(Collection<Integer> indeces) {
		List<PharmacophoreNode> nodeList = new ArrayList<PharmacophoreNode>();
		for(int i : indeces) {
			 nodeList.add(nodes.get(i));
		}
		return nodeList;
	}
	
	public List<PharmacophoreNode> getNodes() {
		return nodes;
	}
	
	public List<int[]> getEdges() {
		return edges;
	}
	
	public int getLinkNodes() {
		return linkerNodes;
	}
	
	public void removeNode(PharmacophoreNode node) {

		int nodeIndex = nodes.indexOf(node);
		List<Integer> edgesToBeDeleted = new ArrayList<Integer>();
		/*
		 * remove edge that connects the two merged nodes, redirect edges that point to the node that is 
		 * to be merged into the target node to the newly merged node
		 */
		for(int e=edges.size()-1;e>=0;e--) {
			if(edges.get(e)[0]==nodeIndex || edges.get(e)[1] == nodeIndex)
				edgesToBeDeleted.add(e);
			
		}
		
		nodes.remove(nodeIndex);
		
		// update edges so that they point to the new,updated node indeces
		for(int e=0;e<edges.size();e++) { //
			int[] edge = edges.get(e);
			if(edge[0]>nodeIndex)
				edge[0]+=-1;
			if(edge[1]>nodeIndex)
				edge[1]+=-1;
		}

		for(int edgeToDelete : edgesToBeDeleted) {
			edges.remove(edgeToDelete);
		}
		update();

	}
	
	public double getSize() {
		double size = nodes.stream().mapToDouble(i -> i.getSize()).
				reduce((a,b) -> a+b).getAsDouble();
		
		return size;
	}
	
	public double getSubTreeSize(List<Integer> edges, int headNode) {
		Collection<Integer> indeces = getNodesFromEdges(edges);
		Collection<PharmacophoreNode> n = getNodes(indeces);
		n.add(nodes.get(headNode));
		double size = n.stream().mapToDouble(i -> i.getSize()).
				reduce((a,b) -> a+b).getAsDouble();

		
		return size;
	}
	
	public double getDirectSim(PharmacophoreTree pTree2) {
		return PharmacophoreNode.getSimilarity(nodes, pTree2.getNodes());
	}
	
	
	public List<Set<Integer>> getAllSubtrees() {
		Set<Set<Integer>> allSubTrees = new HashSet<Set<Integer>>();
		for(int e=0;e<edges.size();e++) {
			for(int i=0;i<2;i++) {
				int[] edge = edges.get(e);
				List<Integer> treeEdgeIndeces = new ArrayList<Integer>();
				List<Integer> treeEdgeParentIndeces = new ArrayList<Integer>();
				treeWalkBFS(edge[i],e,treeEdgeIndeces,treeEdgeParentIndeces);
				List<int[]> cuts = getExtensionCuts(treeEdgeIndeces,treeEdgeParentIndeces);
				for(int[] cut : cuts) {
					List<List<Integer>> sourceTreeEdgeIndeces = new ArrayList<List<Integer>>();
					List<List<Integer>> sourceTreeEdgeParentIndeces = new ArrayList<List<Integer>>();
					List<Integer> sourceTreeHeadNodes = new ArrayList<Integer>();
					Set<Integer> extensionNodes = new HashSet<Integer>();
					List<Integer> cutEdges = new ArrayList<Integer>();
					List<Integer> cutDirs = new ArrayList<Integer>();
					enumerateExtensionCutFull(edge[i],cut, treeEdgeIndeces,
						treeEdgeParentIndeces, sourceTreeEdgeIndeces,sourceTreeEdgeParentIndeces, 
						sourceTreeHeadNodes,extensionNodes, cutEdges, cutDirs);
					allSubTrees.add(extensionNodes);
				}
			}
		}
		//System.out.println(allSubTrees);
		return new ArrayList<Set<Integer>>(allSubTrees);
	}
	
	public void getPathsFromHeadNode(int headNode, List<Set<Integer>> paths, Set<Integer> visitedEdges) {
		Set<Integer> lastPath = paths.get(paths.size()-1);
		Set<Integer> nextPath = new HashSet<Integer> (lastPath);
		nextPath.add(headNode);
		paths.add(nextPath);
		List<Integer> adjacentEdges = adjacencyList.get(headNode);
		for(Integer edge : adjacentEdges) {
			if(visitedEdges.contains(edge))
				continue;
			else {
				int[] e = edges.get(edge);
				int u = e[0];
				int v = e[1];
				if(u==headNode)
					getPathsFromHeadNode(v,paths,visitedEdges);
				else
					getPathsFromHeadNode(u,paths,visitedEdges);
			}
		}
		
	}
	

	
	/**
	 * container of an undirected edge
	 * @author joel
	 *
	 */
	
	public static class BiGramInt {
		public int[] edge;
		public int order;
		int u,v;
		
		public BiGramInt(int[] edge) {
			this(edge,1);
		}
		
		public BiGramInt(int[] edge, int order) {
			this.edge = edge;
			this.order = order;
			if(edge[0]<edge[1]) {
				u = edge[0];
				v = edge[1];
			}
				
			else {
				u = edge[1];
				v = edge[0];
			}

				
		}
		public BiGramInt(int i, int j) {
			this(new int[] {i,j},1);
		}
		@Override
		public boolean equals(Object obj) {

		    if (obj == null) return false;

		    if (!(obj instanceof BiGramInt))

		        return false;

		    if (obj == this)

		        return true;

		    return ((this.u == ((BiGramInt) obj).u) && (this.v == ((BiGramInt) obj).v));

		}
		
		@Override
		public int hashCode() {
		    return u * 31 + v;
		}
		
		}
	
	
		
	
	
}
