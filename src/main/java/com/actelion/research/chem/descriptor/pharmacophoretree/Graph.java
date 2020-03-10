package com.actelion.research.chem.descriptor.pharmacophoretree;

// A Java program to find biconnected components in a given 
//undirected graph 
//taken and slightly modified from: https://www.geeksforgeeks.org/biconnected-components/


import java.io.*; 
import java.util.*; 

//This class represents a directed graph using adjacency 
//list representation 
public class Graph { 

	private Map<Integer,List<Integer>> adj; // Adjacency List 
	private int V; // number of knodes
	private List<List<int[]>> bccs;
	private List<int[]> bcc;

	// Count is number of biconnected components. time is 
	// used to find discovery times 
	private static int count = 0, time = 0; 

	public class Edge { 
		int u; 
		int v; 
     
		Edge(int u, int v) { 
			this.u = u; 
			this.v = v; 
     } 
	}; 

	// Constructor 
	public Graph(Map<Integer,List<Integer>> adj) {
		this.adj = adj;
		V = adj.size();
	} 



 // A recursive function that finds and prints strongly connected 
 // components using DFS traversal 
 // u --> The vertex to be visited next 
 // disc[] --> Stores discovery times of visited vertices 
 // low[] -- >> earliest visited vertex (the vertex with minimum 
 // discovery time) that can be reached from subtree 
 // rooted with current vertex 
 // *st -- >> To store visited edges 
	 public void bccutil(int u, int disc[], int low[], LinkedList<Edge> st, 
              int parent[]) 
 { 

     // Initialize discovery time and low value 
     disc[u] = low[u] = ++time; 
     int children = 0; 

     // Go through all vertices adjacent to this 
     Iterator<Integer> it = adj.get(u).iterator(); 
     while (it.hasNext()) { 
         int v = it.next(); // v is current adjacent of 'u' 

         // If v is not visited yet, then recur for it 
         if (disc[v] == -1) { 
             children++; 
             parent[v] = u; 

             // store the edge in stack 
             st.add(new Edge(u, v)); 
             bccutil(v, disc, low, st, parent); 

             // Check if the subtree rooted with 'v' has a 
             // connection to one of the ancestors of 'u' 
             // Case 1 -- per Strongly Connected Components Article 
             if (low[u] > low[v]) 
                 low[u] = low[v]; 

             // If u is an articulation point, 
             // pop all edges from stack till u -- v 
             if ((disc[u] == 1 && children > 1) || (disc[u] > 1 && low[v] >= disc[u])) { 
                 while (st.getLast().u != u || st.getLast().v != v) { 
                	 bcc.add(new int[] {st.getLast().u,st.getLast().v});
                     st.removeLast(); 
                 } 
                 

                 bcc.add(new int[] {st.getLast().u,st.getLast().v});
                 bccs.add(bcc);
                 bcc = new ArrayList<int[]>();
                 st.removeLast(); 

                 count++; 
             } 
         } 

         // Update low value of 'u' only if 'v' is still in stack 
         // (i.e. it's a back edge, not cross edge). 
         // Case 2 -- per Strongly Connected Components Article 
         else if (v != parent[u] && disc[v] < disc[u] ) { 
             if (low[u] > disc[v]) 
                 low[u] = disc[v]; 

             st.add(new Edge(u, v)); 
         } 
     } 
 } 

 // The function to do DFS traversal. It uses BCCUtil() 
 	public List<List<int[]>> bcc() {
 		 bccs = new ArrayList<List<int[]>>();
 		 bcc = new ArrayList<int[]>();
	     int disc[] = new int[V]; 
	     int low[] = new int[V]; 
	     int parent[] = new int[V]; 
	     LinkedList<Edge> st = new LinkedList<Edge>(); 
	
	     // Initialize disc and low, and parent arrays 
	     for (int i = 0; i < V; i++) { 
	         disc[i] = -1; 
	         low[i] = -1; 
	         parent[i] = -1; 
	     } 
	     for (int i = 0; i < V; i++) { 
	         if (disc[i] == -1) 
	             bccutil(i, disc, low, st, parent); 
	
	         int j = 0; 
	
	         // If stack is not empty, pop all edges from stack 
	         while (st.size() > 0) { 
	             j = 1; 
	             bcc.add(new int[] {st.getLast().u,st.getLast().v});
	             st.removeLast(); 
	         } 
	         if (j == 1) { 
	        	 bccs.add(bcc);
	        	 bcc = new ArrayList<int[]>();
	             count++; 
	         } 
	     } 
	     return bccs;
	 } 
}
